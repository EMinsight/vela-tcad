#!/usr/bin/env python3
"""Raw-input semantic verifier for the Minimal6 inverse-audit package.

This module does not import the Task 8 report, plot, or diagnose modules.  It
reconstructs expected artifacts from the sealed Task 2 input bundle and uses
Task 4-7 diagnostic primitives only for independently specified calculations.
"""

from __future__ import annotations

import argparse
from dataclasses import asdict
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import statistics
import sys


REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from PIL import Image

from scripts.pn2d_minimal6_diagnostics.inverse_contracts import (
    AcceptanceThresholds, Identifiability, Observation, SampleStatus, SupportKind,
)
from scripts.pn2d_minimal6_diagnostics.inverse_inputs import (
    InputBundle, field_inventory, load_input_bundle,
)

from scripts.pn2d_minimal6_diagnostics.independent_science import (
    recompute_science,
)

REPORT_EXCLUSIONS = {"report_manifest.json", "verification.json", "package_manifest.json"}
PACKAGE_EXCLUSIONS = {"package_manifest.json"}
FIGURES = ("potential_field", "qf_gradient", "current_density",
           "alpha_generation", "replacement_matrix")
OBSERVATION_COLUMNS = (
    "solver", "topology", "bias_V", "support_kind", "support_id", "quantity",
    "component", "raw_value", "raw_unit", "value_si", "unit_si",
    "coordinate_frame", "orientation", "conversion", "status", "source_path",
    "source_sha256",
)
CANDIDATE_COLUMNS = (
    "candidate", "quantity", "carrier", "split", "topology", "bias_V",
    "support_kind", "valid_count", "median_abs_error", "p95_abs_error",
    "median_angle_deg", "classification",
)
REPLACEMENT_COLUMNS = (
    "sequence", "step", "factor", "value", "incremental_dex", "closure_abs_dex",
)

_MIRROR_NODE_MAP = {
    "0": "4", "1": "5", "2": "3",
    "3": "2", "4": "0", "5": "1",
}
_MIRROR_VECTOR_QUANTITIES = {
    "ElectricField", "eCurrentDensity", "hCurrentDensity",
}
INVERSE_DEPENDENCIES = (
    "gradient_recovery", "mobility", "current_semantics",
    "impact_driving_field", "alpha_law", "geometric_integration",
    "source_to_node_mapping",
)
_ELEMENTARY_CHARGE_C = 1.602176634e-19
_BOLTZMANN_EV_PER_K = 8.617333262145e-5
_VELA_HEADER_RELATIVE = "source/tracked/include/vela/physics/ImpactIonizationModel.h"
_VELA_MESH_RELATIVE = "source/topologies/{topology}/mesh.json"
_VELA_DECK_PREFIX = "source/decks/"
_IMPACT_STRING_MEMBERS = (
    ("model", "model"), ("parameterSet", "parameter_set"),
    ("drivingForce", "driving_force"), ("generation", "generation"),
    ("currentApproximation", "current_approximation"),
    ("currentMagnitudeMode", "current_magnitude_mode"),
    ("cellReconstructedMidpointDensity", "cell_reconstructed_midpoint_density"),
    ("drivingForceInterpolation", "driving_force_interpolation"),
    ("quasiFermiGradientDiscretization", "quasi_fermi_gradient_discretization"),
    ("sourceVolumePolicy", "source_volume_policy"),
    ("sourceMappingMode", "source_mapping_mode"),
    ("edgeSourcePartition", "edge_source_partition"),
)
_IMPACT_NUMBER_MEMBERS = (
    ("electronDrivingForceRefDensity", "electron_driving_force_ref_density_m3"),
    ("holeDrivingForceRefDensity", "hole_driving_force_ref_density_m3"),
    ("sourceGeometryScale", "source_geometry_scale"),
    ("sourceVolumeFactor", "source_volume_factor"),
    ("quasiFermiCarrierTruncation", "quasi_fermi_carrier_truncation"),
    ("minimumField", "minimum_field_V_m"),
    ("aScale", "A_scale"), ("bScale", "B_scale"),
    ("electronA", "electron_A_m_inv"), ("electronB", "electron_B_V_m"),
    ("holeA", "hole_A_m_inv"), ("holeB", "hole_B_V_m"),
    ("carrierVelocity", "carrier_velocity_m_s"),
    ("electronALow", "electron_a_low_m_inv"),
    ("electronAHigh", "electron_a_high_m_inv"),
    ("electronBLow", "electron_b_low_V_m"),
    ("electronBHigh", "electron_b_high_V_m"),
    ("holeALow", "hole_a_low_m_inv"),
    ("holeAHigh", "hole_a_high_m_inv"),
    ("holeBLow", "hole_b_low_V_m"),
    ("holeBHigh", "hole_b_high_V_m"),
    ("switchField", "switch_field_V_m"),
    ("phononEnergy", "phonon_energy_eV"),
    ("referenceTemperature_K", "reference_temperature_K"),
    ("temperature_K", "temperature_K"),
)
_IMPACT_BOOL_MEMBERS = (
    ("debugRawVanOverstraeten", "debug_raw_vanoverstraeten"),
)
_NESTED_REF_KEYS = {
    "electron_ref_density_m3": "electron_driving_force_ref_density_m3",
    "hole_ref_density_m3": "hole_driving_force_ref_density_m3",
}
_DIRECT_STRING_KEYS = frozenset(
    key for _, key in _IMPACT_STRING_MEMBERS
    if key != "driving_force_interpolation"
)
_DIRECT_NUMBER_KEYS = frozenset(
    key for _, key in _IMPACT_NUMBER_MEMBERS
    if key not in set(_NESTED_REF_KEYS.values())
)
_DIRECT_BOOL_KEYS = frozenset(key for _, key in _IMPACT_BOOL_MEMBERS)
_IMPACT_DECK_KEYS = (
    _DIRECT_STRING_KEYS | _DIRECT_NUMBER_KEYS | _DIRECT_BOOL_KEYS
    | {"driving_force_interpolation", "quasi_fermi_carrier_trucation"}
)
_POSITIVE_CONFIG_KEYS = frozenset({
    "source_geometry_scale", "A_scale", "B_scale", "electron_A_m_inv",
    "electron_B_V_m", "hole_A_m_inv", "hole_B_V_m", "carrier_velocity_m_s",
    "electron_a_low_m_inv", "electron_a_high_m_inv", "electron_b_low_V_m",
    "electron_b_high_V_m", "hole_a_low_m_inv", "hole_a_high_m_inv",
    "hole_b_low_V_m", "hole_b_high_V_m", "switch_field_V_m",
    "phonon_energy_eV", "reference_temperature_K", "temperature_K",
})
_ALLOWED_STRINGS = {
    "model": {"van_overstraeten"},
    "parameter_set": {"default"},
    "driving_force": {
        "electric_field", "quasi_fermi_gradient", "grad_potential_parallel_j",
        "effective_field_parallel_j",
    },
    "generation": {"carrier_density", "current_density"},
    "current_approximation": {
        "mobility_density_gradient", "density_gradient", "grad_qf",
        "cell_reconstructed", "psi_gradient_proxy",
        "cell_current_reconstructed", "cell_vector_current_reconstructed",
        "conserved_total_current",
    },
    "current_magnitude_mode": {"edge_scalar_abs", "dual_face_vector_mag"},
    "cell_reconstructed_midpoint_density": {
        "bernoulli", "arithmetic", "gss_logistic",
    },
    "driving_force_interpolation": {"none", "quasi_fermi_to_electric_field"},
    "quasi_fermi_gradient_discretization": {"edge_difference", "cell_gradient"},
    "source_volume_policy": {"genius_truncated", "edge_half_box", "edge_box"},
    "source_mapping_mode": {
        "node_F_node_alpha_node_G", "edge_F_edge_alpha_edge_G_to_node",
        "cell_F_cell_alpha_cell_G_to_node", "triangle_gss_gradqf_truncated",
    },
    "edge_source_partition": {"symmetric", "qf_gradient"},
}

_NUMBER = r"([0-9.eE+-]+)"


def _independent_triangle_gradient(points, values):
    (x0, y0), (x1, y1), (x2, y2) = (
        (float(point[0]), float(point[1])) for point in points
    )
    f0, f1, f2 = (float(value) for value in values)
    determinant = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)
    if not math.isfinite(determinant) or abs(determinant) <= 1.0e-300:
        raise ValueError("degenerate triangle")
    return (
        ((f1 - f0) * (y2 - y0) - (f2 - f0) * (y1 - y0)) / determinant,
        ((x1 - x0) * (f2 - f0) - (x2 - x0) * (f1 - f0)) / determinant,
    )


def _independent_current_inverted_gradient(carrier, density, mobility, current):
    sign = -1.0 if carrier == "electron" else 1.0
    denominator = sign * _ELEMENTARY_CHARGE_C * float(mobility) * float(density)
    if denominator == 0.0 or not math.isfinite(denominator):
        raise ValueError("invalid independent current inversion denominator")
    return float(current[0]) / denominator, float(current[1]) / denominator


def _independent_impact_generation(alpha_n, current_n, alpha_p, current_p):
    result = (
        float(alpha_n) * math.hypot(float(current_n[0]), float(current_n[1]))
        + float(alpha_p) * math.hypot(float(current_p[0]), float(current_p[1]))
    ) / _ELEMENTARY_CHARGE_C
    if not math.isfinite(result):
        raise ValueError("independent impact generation is nonfinite")
    return result


def _independent_van_overstraeten_alpha(
    alpha, *, prefactor, critical_field, gamma, branch=None,
    switch_field=None,
):
    alpha = float(alpha)
    prefactor = float(prefactor)
    critical_field = float(critical_field)
    gamma = float(gamma)
    if min(alpha, prefactor, critical_field, gamma) <= 0.0:
        return None, SampleStatus.BELOW_FLOOR
    if alpha >= gamma * prefactor:
        return None, SampleStatus.BRANCH_AMBIGUOUS
    field = -gamma * critical_field / math.log(alpha / (gamma * prefactor))
    if branch == "low" and switch_field is not None and field >= switch_field:
        return None, SampleStatus.BRANCH_AMBIGUOUS
    if branch == "high" and switch_field is not None and field < switch_field:
        return None, SampleStatus.BRANCH_AMBIGUOUS
    return field, SampleStatus.VALID


def _independent_forward_van_overstraeten(field, parameters, carrier):
    field = abs(float(field))
    if field <= 0.0:
        return 0.0, "low"
    branch = "low" if field < parameters["switch_field_V_m"] else "high"
    prefactor, critical = parameters[carrier][branch]
    gamma = parameters["gamma"]
    return (
        gamma * prefactor * math.exp(-gamma * critical / field),
        branch,
    )

_MIRROR_ABSOLUTE_TOLERANCE_BY_UNIT_SI = {"V/m": 1.0e-8}
def _config_number(value, label, *, positive=False):
    if isinstance(value, bool):
        raise ValueError(f"{label} must be numeric")
    try:
        number = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{label} must be numeric") from error
    if not math.isfinite(number) or (positive and number <= 0.0):
        qualifier = " and positive" if positive else ""
        raise ValueError(f"{label} must be finite{qualifier}")
    return number


def _verify_header_defaults(text):
    config = {}
    for member, key in _IMPACT_STRING_MEMBERS:
        matches = re.findall(
            rf'\bstd::string\s+{re.escape(member)}\s*=\s*"([^"]*)"\s*;', text
        )
        if len(matches) != 1:
            raise ValueError(
                f"sealed Vela header lacks unique string {member}"
            )
        config[key] = matches[0]
    for member, key in _IMPACT_NUMBER_MEMBERS:
        matches = re.findall(
            rf"\bReal\s+{re.escape(member)}\s*=\s*{_NUMBER}\s*;", text
        )
        if len(matches) != 1:
            raise ValueError(
                f"sealed Vela header lacks unique numeric {member}"
            )
        config[key] = _config_number(
            matches[0], f"sealed Vela header {member}"
        )
    for member, key in _IMPACT_BOOL_MEMBERS:
        matches = re.findall(
            rf"\bbool\s+{re.escape(member)}\s*=\s*(true|false)\s*;", text
        )
        if len(matches) != 1:
            raise ValueError(
                f"sealed Vela header lacks unique bool {member}"
            )
        config[key] = matches[0] == "true"
    if len(config) != 38:
        raise ValueError("sealed Vela header default config is incomplete")
    return config


def _verify_cross_field_config(config: dict[str, object], label: str) -> None:
    debug = bool(config["debug_raw_vanoverstraeten"])
    driving = config["driving_force"]
    current_aligned = (
        not debug
        and driving in {"grad_potential_parallel_j", "effective_field_parallel_j"}
    )
    if current_aligned and (
        config["generation"] != "current_density"
        or config["current_approximation"] not in {"density_gradient", "grad_qf"}
    ):
        raise ValueError(
            f"{label} cross-field config invalid: current-aligned driving force"
        )
    if (
        config["driving_force_interpolation"] != "none"
        and driving != "quasi_fermi_gradient"
        and not debug
    ):
        raise ValueError(
            f"{label} cross-field config invalid: driving-force interpolation"
        )
    if (
        config["quasi_fermi_gradient_discretization"] == "cell_gradient"
        and driving != "quasi_fermi_gradient"
        and not debug
    ):
        raise ValueError(
            f"{label} cross-field config invalid: quasi-Fermi discretization"
        )
    triangle_gss = config["source_mapping_mode"] == "triangle_gss_gradqf_truncated"
    if triangle_gss:
        required = (
            config["generation"] == "current_density"
            and driving == "quasi_fermi_gradient"
            and config["current_approximation"] == "cell_reconstructed"
            and config["current_magnitude_mode"] == "edge_scalar_abs"
            and config["cell_reconstructed_midpoint_density"] == "gss_logistic"
            and config["quasi_fermi_gradient_discretization"] == "cell_gradient"
            and config["source_volume_policy"] == "genius_truncated"
            and config["source_volume_factor"] == 0.0
            and config["source_geometry_scale"] == 1.0
            and config["edge_source_partition"] == "symmetric"
            and config["driving_force_interpolation"] == "none"
            and config["quasi_fermi_carrier_truncation"] == 0.0
            and config["minimum_field_V_m"] == 0.0
            and config["electron_driving_force_ref_density_m3"] == 0.0
            and config["hole_driving_force_ref_density_m3"] == 0.0
        )
        if not required:
            raise ValueError(
                f"{label} cross-field config invalid: canonical GSS mapping"
            )
    if (
        config["cell_reconstructed_midpoint_density"] == "gss_logistic"
        and not triangle_gss
    ):
        raise ValueError(
            f"{label} cross-field config invalid: GSS midpoint mapping"
        )


def _verify_effective_config(config, label):
    expected = (
        {key for _, key in _IMPACT_STRING_MEMBERS}
        | {key for _, key in _IMPACT_NUMBER_MEMBERS}
        | {key for _, key in _IMPACT_BOOL_MEMBERS}
    )
    if set(config) != expected or len(config) != 38:
        raise ValueError(
            f"{label} complete canonical effective config mismatch"
        )
    for key, allowed in _ALLOWED_STRINGS.items():
        if not isinstance(config[key], str) or config[key] not in allowed:
            raise ValueError(f"{label} unsupported {key}")
    for _, key in _IMPACT_NUMBER_MEMBERS:
        number = _config_number(
            config[key],
            f"{label} {key}",
            positive=key in _POSITIVE_CONFIG_KEYS,
        )
        if key not in _POSITIVE_CONFIG_KEYS and number < 0.0:
            raise ValueError(f"{label} {key} must be non-negative")
        config[key] = number
    volume = config["source_volume_factor"]
    if volume != 0.0 and not 0.5 <= volume <= 1.0:
        raise ValueError(
            f"{label} source_volume_factor must be 0 or within [0.5, 1.0]"
        )
    if type(config["debug_raw_vanoverstraeten"]) is not bool:
        raise ValueError(f"{label} invalid debug_raw_vanoverstraeten")
    _verify_cross_field_config(config, label)



def _verify_impact_entries(value):
    entries = []
    if isinstance(value, dict):
        for key, child in value.items():
            if key == "impact_ionization":
                entries.append(child)
            else:
                entries.extend(_verify_impact_entries(child))
    elif isinstance(value, list):
        for child in value:
            entries.extend(_verify_impact_entries(child))
    return entries


def _verify_deck_effective(path, relative, defaults):
    entries = _verify_impact_entries(_load_json(path))
    if len(entries) != 1:
        raise ValueError(
            f"sealed Vela deck {relative} impact config mismatch"
        )
    raw = entries[0]
    if isinstance(raw, str):
        overrides = {"model": raw}
    elif isinstance(raw, dict):
        unknown = set(raw) - _IMPACT_DECK_KEYS
        if unknown:
            raise ValueError(
                f"sealed Vela deck {relative} unknown impact override "
                f"{sorted(unknown)[0]}"
            )
        overrides = {}
        for key in _DIRECT_STRING_KEYS:
            if key in raw:
                if not isinstance(raw[key], str):
                    raise ValueError(
                        f"sealed Vela deck {relative} invalid {key}"
                    )
                overrides[key] = raw[key]
        for key in _DIRECT_NUMBER_KEYS:
            if key in raw:
                overrides[key] = _config_number(
                    raw[key], f"sealed Vela deck {relative} {key}"
                )
        for key in _DIRECT_BOOL_KEYS:
            if key in raw:
                if type(raw[key]) is not bool:
                    raise ValueError(
                        f"sealed Vela deck {relative} invalid {key}"
                    )
                overrides[key] = raw[key]
        if "quasi_fermi_carrier_trucation" in raw:
            if "quasi_fermi_carrier_truncation" in raw:
                raise ValueError(
                    f"sealed Vela deck {relative} declares both truncation spellings"
                )
            overrides["quasi_fermi_carrier_truncation"] = _config_number(
                raw["quasi_fermi_carrier_trucation"],
                f"sealed Vela deck {relative} quasi_fermi_carrier_trucation",
            )
        if "driving_force_interpolation" in raw:
            interpolation = raw["driving_force_interpolation"]
            if isinstance(interpolation, str):
                overrides["driving_force_interpolation"] = interpolation
            elif isinstance(interpolation, dict):
                unknown_nested = set(interpolation) - (
                    {"mode"} | set(_NESTED_REF_KEYS)
                )
                if unknown_nested:
                    raise ValueError(
                        f"sealed Vela deck {relative} unknown interpolation "
                        f"override {sorted(unknown_nested)[0]}"
                    )
                if "mode" in interpolation:
                    if not isinstance(interpolation["mode"], str):
                        raise ValueError(
                            f"sealed Vela deck {relative} invalid interpolation mode"
                        )
                    overrides["driving_force_interpolation"] = interpolation["mode"]
                for source, target in _NESTED_REF_KEYS.items():
                    if source in interpolation:
                        overrides[target] = _config_number(
                            interpolation[source],
                            f"sealed Vela deck {relative} {source}",
                        )
            else:
                raise ValueError(
                    f"sealed Vela deck {relative} invalid "
                    "driving_force_interpolation"
                )
    else:
        raise ValueError(
            f"sealed Vela deck {relative} impact config mismatch"
        )
    effective = dict(defaults)
    effective.update(overrides)
    _verify_effective_config(effective, f"sealed Vela deck {relative}")
    return effective


def _verify_parameters_from_effective(config):
    reference = float(config["reference_temperature_K"])
    temperature = float(config["temperature_K"])
    phonon = float(config["phonon_energy_eV"])
    gamma = (
        math.tanh(phonon / (2.0 * _BOLTZMANN_EV_PER_K * reference))
        / math.tanh(phonon / (2.0 * _BOLTZMANN_EV_PER_K * temperature))
    )
    a_scale = float(config["A_scale"])
    b_scale = float(config["B_scale"])
    parameters = {
        "gamma": gamma,
        "switch_field_V_m": float(config["switch_field_V_m"]),
        "electron": {
            "low": [
                a_scale * float(config["electron_a_low_m_inv"]),
                b_scale * float(config["electron_b_low_V_m"]),
            ],
            "high": [
                a_scale * float(config["electron_a_high_m_inv"]),
                b_scale * float(config["electron_b_high_V_m"]),
            ],
        },
        "hole": {
            "low": [
                a_scale * float(config["hole_a_low_m_inv"]),
                b_scale * float(config["hole_b_low_V_m"]),
            ],
            "high": [
                a_scale * float(config["hole_a_high_m_inv"]),
                b_scale * float(config["hole_b_high_V_m"]),
            ],
        },
        "phonon_energy_eV": phonon,
        "reference_temperature_K": reference,
        "temperature_K": temperature,
    }
    return parameters, _BOLTZMANN_EV_PER_K * temperature


def _independent_vela_context(roots):
    root = Path(roots["vela_root"]).resolve()
    manifest = _load_json(root / "manifest.json")
    member_hashes = manifest.get("member_sha256", {})
    if not isinstance(member_hashes, dict):
        raise ValueError("sealed Vela manifest member hashes mismatch")

    def bound_path(relative):
        path = root / relative
        if relative not in member_hashes or _sha256(path) != member_hashes[relative]:
            raise ValueError(f"sealed Vela member hash mismatch: {relative}")
        return path

    header_path = bound_path(_VELA_HEADER_RELATIVE)
    defaults = _verify_header_defaults(
        header_path.read_text(encoding="utf-8")
    )
    mesh, mesh_sources = {}, {}
    for topology in ("sketch", "mirror"):
        relative = _VELA_MESH_RELATIVE.format(topology=topology)
        path = bound_path(relative)
        payload = _load_json(path)
        if not isinstance(payload, dict) or payload.get("coordinate_unit") != "um":
            raise ValueError("sealed Vela mesh coordinate contract mismatch")
        coordinates = {
            str(item["id"]): (float(item["x"]), float(item["y"]))
            for item in payload.get("nodes", [])
        }
        triangles = {}
        for item in payload.get("triangles", []):
            nodes = tuple(str(node) for node in item.get("node_ids", []))
            if len(nodes) != 3 or len(set(nodes)) != 3:
                raise ValueError("sealed Vela mesh triangle contract mismatch")
            points = tuple(coordinates[node] for node in nodes)
            twice_area = (
                (points[1][0] - points[0][0])
                * (points[2][1] - points[0][1])
                - (points[2][0] - points[0][0])
                * (points[1][1] - points[0][1])
            )
            if not math.isfinite(twice_area) or abs(twice_area) <= 1.0e-300:
                raise ValueError("sealed Vela mesh has degenerate triangle")
            triangles[str(item["id"])] = nodes
        if not triangles:
            raise ValueError("sealed Vela mesh has no triangles")
        mesh[topology] = {"triangles": triangles}
        mesh_sources[topology] = {
            "logical_id": f"vela:{relative}",
            "sha256": member_hashes[relative],
        }

    deck_relatives = sorted(
        relative for relative in member_hashes
        if relative.startswith(_VELA_DECK_PREFIX) and relative.endswith(".json")
    )
    if len(deck_relatives) != 40:
        raise ValueError(
            "sealed Vela context requires exactly 40 hash-bound simulation decks"
        )
    deck_sources = []
    effective_configs = []
    for relative in deck_relatives:
        path = bound_path(relative)
        effective = _verify_deck_effective(path, relative, defaults)
        effective_configs.append(effective)
        deck_sources.append({
            "logical_id": f"vela:{relative}",
            "sha256": member_hashes[relative],
            "effective_config": dict(effective),
        })
    effective_config = effective_configs[0]
    if any(config != effective_config for config in effective_configs[1:]):
        raise ValueError("sealed Vela deck effective config mismatch")
    parameters, thermal_voltage = _verify_parameters_from_effective(
        effective_config
    )
    provenance = {
        "parameter_identity": "sealed_vela_production_defaults",
        "parameter_source": f"vela:{_VELA_HEADER_RELATIVE}",
        "parameter_sha256": member_hashes[_VELA_HEADER_RELATIVE],
        "mesh_sources": mesh_sources,
        "deck_sources": deck_sources,
        "effective_impact_ionization_config": dict(effective_config),
    }
    return {
        "mesh_by_topology": mesh,
        "parameters": parameters,
        "thermal_voltage_V": thermal_voltage,
        "provenance": provenance,
    }


def _expected_mirror_invariance(rows: tuple[Observation, ...]) -> dict[str, object]:
    """Independently reconstruct the labelled vertical-reflection gate."""
    node_rows = tuple(
        row for row in rows
        if row.support_kind is SupportKind.NODE
        and row.topology in {"sketch", "mirror"}
    )
    index = {}
    duplicates = 0
    for row in node_rows:
        key = (
            row.solver, row.topology, float(row.bias_V), str(row.support_id),
            row.quantity, row.component,
        )
        if key in index:
            duplicates += 1
        else:
            index[key] = row

    valid_pairs = 0
    matching_nonvalid_pairs = 0
    mismatches = duplicates
    unpaired = 0
    consumed_mirror_keys = set()
    for sketch in sorted(
        (row for row in node_rows if row.topology == "sketch"),
        key=lambda row: (
            row.solver, float(row.bias_V), str(row.support_id),
            row.quantity, row.component,
        ),
    ):
        reflected_node = _MIRROR_NODE_MAP.get(str(sketch.support_id))
        reflected_key = (
            sketch.solver, "mirror", float(sketch.bias_V), reflected_node,
            sketch.quantity, sketch.component,
        )
        mirror = index.get(reflected_key) if reflected_node is not None else None
        if mirror is None:
            unpaired += 1
            continue
        consumed_mirror_keys.add(reflected_key)
        if sketch.unit_si != mirror.unit_si or sketch.status is not mirror.status:
            mismatches += 1
            continue
        if sketch.status is not SampleStatus.VALID:
            matching_nonvalid_pairs += 1
            continue
        if sketch.value_si is None or mirror.value_si is None:
            mismatches += 1
            continue
        sketch_value = float(sketch.value_si)
        mirror_value = float(mirror.value_si)
        if not math.isfinite(sketch_value) or not math.isfinite(mirror_value):
            mismatches += 1
            continue
        if sketch.quantity == "coordinate" and sketch.component == "y":
            reflected_value = 0.5e-6 - sketch_value
        elif (sketch.quantity in _MIRROR_VECTOR_QUANTITIES
              and sketch.component == "component1"):
            reflected_value = -sketch_value
        else:
            reflected_value = sketch_value
        valid_pairs += 1
        if not math.isclose(
            mirror_value,
            reflected_value,
            rel_tol=1.0e-9,
            abs_tol=_MIRROR_ABSOLUTE_TOLERANCE_BY_UNIT_SI.get(sketch.unit_si, 0.0),
        ):
            mismatches += 1

    mirror_keys = {key for key in index if key[1] == "mirror"}
    unpaired += len(mirror_keys - consumed_mirror_keys)
    passed = valid_pairs > 0 and mismatches == 0 and unpaired == 0
    return {
        "status": "pass" if passed else "fail",
        "reflection": "(x,y)->(x,0.5um-y)",
        "node_map": dict(_MIRROR_NODE_MAP),
        "vector_transform": "(vx,vy)->(vx,-vy)",
        "relative_tolerance": 1.0e-9,
        "absolute_tolerance_by_unit_si": dict(_MIRROR_ABSOLUTE_TOLERANCE_BY_UNIT_SI),
        "valid_pair_count": valid_pairs,
        "matching_nonvalid_pair_count": matching_nonvalid_pairs,
        "mismatch_count": mismatches,
        "unpaired_count": unpaired,
    }


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _pixel_sha256(path: Path) -> str:
    with Image.open(path) as image:
        return hashlib.sha256(image.convert("RGBA").tobytes()).hexdigest()


def _load_json(path: Path) -> object:
    try:
        return json.loads(path.read_text(encoding="utf-8"),
                          parse_constant=lambda value: (_ for _ in ()).throw(ValueError(value)))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"invalid JSON: {path.name}") from error


def _write_json_bytes(payload: object) -> bytes:
    return (json.dumps(payload, sort_keys=True, indent=2, allow_nan=False) + "\n").encode("utf-8")


def _format_csv(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError("expected CSV value is non-finite")
        return format(value, ".17g")
    if hasattr(value, "value"):
        return str(value.value)
    return str(value)


def _read_csv(path: Path, columns: tuple[str, ...]) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if tuple(reader.fieldnames or ()) != columns:
            raise ValueError(f"{path.name} column contract mismatch")
        return list(reader)


def _raw_ledger(bundle: InputBundle, roots: dict[str, str]) -> list[dict]:
    ledger = []
    for logical, digest in bundle.input_hashes:
        solver, relative = logical.split(":", 1)
        root_key = "supplemental_sentaurus_root" if solver == "supplemental" else f"{solver}_root"
        ledger.append({
            "logical_id": logical, "solver": solver, "relative_path": relative,
            "path": str((Path(roots[root_key]) / relative).resolve()), "sha256": digest,
        })
    for solver, root_key in (("vela", "vela_root"), ("sentaurus", "sentaurus_root"),
                             ("supplemental", "supplemental_sentaurus_root")):
        for relative in ("manifest.json", "seal.json"):
            path = (Path(roots[root_key]) / relative).resolve()
            ledger.append({"logical_id": f"{solver}:{relative}", "solver": solver,
                           "relative_path": relative, "path": str(path),
                           "sha256": _sha256(path)})
    return sorted(ledger, key=lambda item: item["logical_id"])


def _load_raw_bundle(report: dict) -> tuple[InputBundle, dict[str, str], list[dict]]:
    try:
        provenance = report["payload"]["localization_control"]["input_provenance"]
        roots = provenance["input_roots"]
        persisted_ledger = provenance["raw_inputs"]
    except (KeyError, TypeError) as error:
        raise ValueError("authoritative raw-input provenance is missing") from error
    required = {"vela_root", "sentaurus_root", "supplemental_sentaurus_root"}
    if not isinstance(roots, dict) or set(roots) != required:
        raise ValueError("authoritative input roots mismatch")
    canonical_roots = {key: str(Path(value).resolve()) for key, value in roots.items()}
    if canonical_roots != roots:
        raise ValueError("authoritative input roots are not exact resolved paths")
    bundle = load_input_bundle(roots["vela_root"], roots["sentaurus_root"],
                               roots["supplemental_sentaurus_root"])
    expected_ledger = _raw_ledger(bundle, roots)
    if persisted_ledger != expected_ledger:
        raise ValueError("authoritative raw-input ledger semantic mismatch")
    for item in expected_ledger:
        path = Path(item["path"])
        if not path.is_file() or _sha256(path) != item["sha256"]:
            raise ValueError(f"raw input hash mismatch: {item['logical_id']}")
    return bundle, roots, expected_ledger


def _verify_report_manifest(root: Path, raw_ledger: list[dict]) -> dict:
    manifest = _load_json(root / "report_manifest.json")
    if not isinstance(manifest, dict) or manifest.get("schema") != "vela.pn2d_minimal6_inverse_report_manifest.v1":
        raise ValueError("report manifest schema mismatch")
    if set(manifest.get("exclusions", ())) != REPORT_EXCLUSIONS:
        raise ValueError("report manifest exclusions mismatch")
    if manifest.get("inputs") != raw_ledger:
        raise ValueError("report manifest raw-input ledger mismatch")
    actual = {}
    for path in sorted(item for item in root.rglob("*") if item.is_file()):
        relative = path.relative_to(root).as_posix()
        if relative not in REPORT_EXCLUSIONS:
            actual[relative] = _sha256(path)
    if manifest.get("artifacts") != actual:
        raise ValueError("report artifact hash mismatch")
    return manifest


def _verify_package_manifest(root: Path, raw_ledger: list[dict] | None = None) -> None:
    path = root / "package_manifest.json"
    if not path.exists():
        return
    manifest = _load_json(path)
    if not isinstance(manifest, dict) or manifest.get("schema") != "vela.pn2d_minimal6_inverse_package_manifest.v1":
        raise ValueError("package manifest schema mismatch")
    if set(manifest.get("exclusions", ())) != PACKAGE_EXCLUSIONS:
        raise ValueError("package manifest exclusions mismatch")
    for item in manifest.get("raw_inputs", []):
        input_path = Path(item["path"])
        if not input_path.is_file() or _sha256(input_path) != item.get("sha256"):
            raise ValueError(f"package raw input hash mismatch: {item.get('logical_id')}")
    if raw_ledger is not None and manifest.get("raw_inputs") != raw_ledger:
        raise ValueError("package raw-input ledger mismatch")
    actual = {}
    for member in sorted(item for item in root.rglob("*") if item.is_file()):
        relative = member.relative_to(root).as_posix()
        if relative not in PACKAGE_EXCLUSIONS:
            actual[relative] = _sha256(member)
    if manifest.get("artifacts") != actual:
        raise ValueError("package manifest hash mismatch")


def _observation_expected(row: Observation) -> dict[str, str]:
    values = {
        "solver": row.solver, "topology": row.topology, "bias_V": row.bias_V,
        "support_kind": row.support_kind.value, "support_id": row.support_id,
        "quantity": row.quantity, "component": row.component, "raw_value": row.raw_value,
        "raw_unit": row.raw_unit, "value_si": row.value_si, "unit_si": row.unit_si,
        "coordinate_frame": row.coordinate_frame, "orientation": row.orientation,
        "conversion": row.conversion, "status": row.status.value,
        "source_path": row.source_path, "source_sha256": row.source_sha256,
    }
    return {column: _format_csv(values[column]) for column in OBSERVATION_COLUMNS}


def _verify_observations(root: Path, rows: tuple[Observation, ...]) -> None:
    for support in SupportKind:
        expected = [_observation_expected(row) for row in rows if row.support_kind is support]
        actual = _read_csv(root / f"observations_{support.value}.csv", OBSERVATION_COLUMNS)
        if actual != expected:
            raise ValueError(f"observations_{support.value}.csv raw semantic mismatch")


def _index(rows: tuple[Observation, ...]) -> dict[tuple, Observation]:
    return {(row.solver, row.topology, row.bias_V, str(row.support_id), row.quantity,
             row.component): row for row in rows}


def _finite(row: Observation | None) -> float | None:
    if row is None or row.status is not SampleStatus.VALID or row.value_si is None:
        return None
    value = float(row.value_si)
    return value if math.isfinite(value) else None


def _candidate_csv(metrics: list[dict]) -> list[dict[str, str]]:
    return [{column: _format_csv(row.get(column)) for column in CANDIDATE_COLUMNS}
            for row in metrics]


def _expected_replacement(reported: object) -> dict:
    if not isinstance(reported, dict):
        raise ValueError("typed replacement result must be an object")
    fixed_candidates = sorted({
        "triangle_minus_grad_psi", "node_area_weighted_minus_grad_psi",
        "edge_area_weighted_minus_grad_psi", "signed_edge_minus_delta_psi_over_h",
        "triangle_qf_gradient_current", "node_area_weighted_qf_gradient_current",
        "edge_area_weighted_qf_gradient_current", "signed_edge_qf_difference_current",
        "current_inverted_qf_gradient", "signed_edge_sg_density_current",
        "signed_edge_drift_diffusion_current", "electric_field_magnitude",
        "qf_gradient_magnitude", "electric_field_current_aligned",
        "qf_gradient_current_aligned",
    })
    expected_fixed = {
        "status": SampleStatus.MISSING_FIELD.value,
        "unavailable_factor": "mobility",
        "dependency_order": list(INVERSE_DEPENDENCIES),
        "baseline": None, "one_factor": [], "forward": [], "reverse": [],
        "full_replacement": None, "adjacent_interactions": [], "closure": None,
        "evidence_source": "typed_candidate_evidence",
        "evidence_candidates": fixed_candidates,
    }
    for key, expected in expected_fixed.items():
        if reported.get(key) != expected:
            raise ValueError(f"typed replacement {key} mismatch")
    observed = reported.get("observed_prediction_dex")
    if not isinstance(observed, dict) or set(observed) != {
        "gradient_recovery", "current_semantics", "impact_driving_field",
    }:
        raise ValueError("typed replacement observed-factor mismatch")
    if not all(
        isinstance(value, (int, float)) and not isinstance(value, bool)
        and math.isfinite(float(value))
        for value in observed.values()
    ):
        raise ValueError("typed replacement evidence must be finite")
    return reported


def _replacement_csv(replacement: dict) -> list[dict[str, str]]:
    if replacement["closure"] is None:
        row = {
            "sequence": "unavailable", "step": 0,
            "factor": replacement["unavailable_factor"], "value": None,
            "incremental_dex": None, "closure_abs_dex": None,
        }
        return [{
            column: _format_csv(row.get(column))
            for column in REPLACEMENT_COLUMNS
        }]
    rows = []
    for sequence in ("one_factor", "forward", "reverse"):
        for step, item in enumerate(replacement[sequence]):
            row = {"sequence": sequence, "step": step, "factor": item["factor"],
                   "value": item["value"],
                   "incremental_dex": item.get("incremental_dex", item.get("delta_dex")),
                   "closure_abs_dex": None}
            rows.append({column: _format_csv(row.get(column)) for column in REPLACEMENT_COLUMNS})
    for name in ("forward_abs_dex", "reverse_abs_dex", "direct_abs_dex"):
        row = {"sequence": "closure", "step": 0, "factor": name, "value": None,
               "incremental_dex": None, "closure_abs_dex": replacement["closure"][name]}
        rows.append({column: _format_csv(row.get(column)) for column in REPLACEMENT_COLUMNS})
    return rows


def _expected_semantic(rows: tuple[Observation, ...], replacement: dict, context) -> dict:
    index = _index(rows)
    states = sorted({(row.topology, row.bias_V) for row in rows if row.solver == "sentaurus"})
    state = states[0]
    nodes = sorted({str(row.support_id) for row in rows if row.solver == "sentaurus"
                    and (row.topology, row.bias_V) == state and row.quantity == "coordinate"})
    topology = context["mesh_by_topology"].get(state[0], {})
    triangles = topology.get("triangles", {})
    cell_id, triangle_nodes = next(iter(triangles.items()), (None, ()))
    replay: dict[str, object] = {"state": [state[0], state[1]]}
    if len(triangle_nodes) == 3 and all(node in nodes for node in triangle_nodes):
        points, potentials = [], []
        for node in triangle_nodes:
            points.append([_finite(index.get(("sentaurus", *state, node, "coordinate", "x"))),
                           _finite(index.get(("sentaurus", *state, node, "coordinate", "y")))])
            potentials.append(_finite(index.get(("sentaurus", *state, node,
                                                  "ElectrostaticPotential", "component0"))))
        try:
            gradient = _independent_triangle_gradient(points, potentials)
            replay["triangle_gradient"] = {
                "status": "valid", "support_kind": "cell", "support_id": cell_id,
                "support_ids": list(triangle_nodes), "points_m": points,
                "values_V": potentials, "value_V_per_m": list(gradient),
            }
        except (TypeError, ValueError):
            replay["triangle_gradient"] = {"status": "incompatible_support"}
    else:
        replay["triangle_gradient"] = {"status": "insufficient_data"}
    node = triangle_nodes[0] if triangle_nodes else nodes[0]
    density = _finite(index.get(("sentaurus", *state, node, "eDensity", "component0")))
    mobility = _finite(index.get(("sentaurus", *state, node, "eMobility", "component0")))
    jx = _finite(index.get(("sentaurus", *state, node, "eCurrentDensity", "component0")))
    jy = _finite(index.get(("sentaurus", *state, node, "eCurrentDensity", "component1")))
    if None not in (density, mobility, jx, jy) and density > 0.0 and mobility > 0.0:
        gradient = _independent_current_inverted_gradient("electron", density, mobility, (jx, jy))
        replay["current_inverted_gradient"] = {
            "status": "valid", "carrier": "electron", "density_m3": density,
            "mobility_m2_per_Vs": mobility, "current_A_per_m2": [jx, jy],
            "value_V_per_m": list(gradient)}
    else:
        replay["current_inverted_gradient"] = {"status": "insufficient_data"}
    alpha_n = _finite(index.get(("sentaurus", *state, node, "eAlphaAvalanche", "component0")))
    alpha_p = _finite(index.get(("sentaurus", *state, node, "hAlphaAvalanche", "component0")))
    hjx = _finite(index.get(("sentaurus", *state, node, "hCurrentDensity", "component0")))
    hjy = _finite(index.get(("sentaurus", *state, node, "hCurrentDensity", "component1")))
    if alpha_n is not None and alpha_n > 0.0:
        parameters = context["parameters"]
        gamma = float(parameters["gamma"])
        switch_field = float(parameters["switch_field_V_m"])
        inversion = None
        for branch in ("low", "high"):
            prefactor, critical_field = parameters["electron"][branch]
            field, status = _independent_van_overstraeten_alpha(
                alpha_n, prefactor=prefactor, critical_field=critical_field,
                gamma=gamma, branch=branch, switch_field=switch_field,
            )
            candidate = {
                "status": status, "alpha_m_inv": alpha_n,
                "parameter_identity": context["provenance"]["parameter_identity"],
                "branch": branch, "prefactor_m_inv": prefactor,
                "critical_field_V_per_m": critical_field, "gamma": gamma,
                "switch_field_V_per_m": switch_field, "field_V_per_m": field,
            }
            if inversion is None or status is SampleStatus.VALID:
                inversion = candidate
            if status is SampleStatus.VALID:
                break
        replay["inverse_alpha"] = inversion
    else:
        replay["inverse_alpha"] = {"status": "insufficient_data"}
    if None not in (alpha_n, alpha_p, jx, jy, hjx, hjy):
        generation = _independent_impact_generation(alpha_n, (jx, jy), alpha_p, (hjx, hjy))
        replay["generation"] = {"status": "valid", "alpha_n_m_inv": alpha_n,
                                "alpha_p_m_inv": alpha_p, "jn_A_per_m2": [jx, jy],
                                "jp_A_per_m2": [hjx, hjy], "value_m3_s_inv": generation}
    else:
        replay["generation"] = {"status": "insufficient_data"}
    native = [_finite(index.get(("sentaurus", *state, current, "ImpactIonization", "component0")))
              for current in triangle_nodes]
    triangle = replay["triangle_gradient"]
    if triangle["status"] == "valid" and len(native) == 3 and all(value is not None for value in native):
        (x0, y0), (x1, y1), (x2, y2) = triangle["points_m"]
        area = 0.5 * abs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0))
        replay["generation_support_integral"] = {
            "status": "valid", "support_kind": "triangle_from_raw_nodes",
            "support_ids": list(triangle_nodes), "native_generation_m3_s_inv": native,
            "area_m2": area, "depth_m": 0.01,
            "integral_s_inv": statistics.mean(native) * area * 0.01}
    else:
        replay["generation_support_integral"] = {"status": "insufficient_data"}
    replay["replacement_closure"] = replacement["closure"]
    return replay


def _sentaurus_version(roots: dict[str, str]) -> str:
    versions = set()
    for key in ("sentaurus_root", "supplemental_sentaurus_root"):
        manifest = _load_json(Path(roots[key]) / "manifest.json")
        value = manifest.get("provenance", {}).get("sentaurus_version", manifest.get("sentaurus_version"))
        if isinstance(value, str) and value:
            versions.add(value)
    if len(versions) > 1:
        raise ValueError("raw inputs declare inconsistent Sentaurus versions")
    return next(iter(versions)) if versions else "not_declared_by_input_contract"


def _verify_raw_semantics(root: Path, report: dict, bundle: InputBundle,
                          roots: dict[str, str], rows: tuple[Observation, ...]) -> list[str]:
    payload = report["payload"]
    if report.get("schema") != "vela.pn2d_minimal6_physics_inverse_audit.v1" or report.get("diagnostic_only") is not True:
        raise ValueError("authoritative report schema mismatch")
    if report.get("phase_base") != "a5524cf" or payload.get("production_cpp_changed") is not False:
        raise ValueError("production baseline guard mismatch")
    if payload.get("discovery_keys") != [[a, b] for a, b in bundle.discovery_keys] or payload.get("holdout_keys") != [[a, b] for a, b in bundle.holdout_keys]:
        raise ValueError("raw discovery/holdout reconstruction mismatch")
    if payload.get("thresholds") != asdict(AcceptanceThresholds()):
        raise ValueError("raw threshold reconstruction mismatch")
    if payload.get("field_inventory") != field_inventory(bundle):
        raise ValueError("raw field inventory reconstruction mismatch")
    counts = {status.value: 0 for status in SampleStatus}
    for row in rows:
        counts[row.status.value] += 1
    if payload.get("sample_status_counts") != counts:
        raise ValueError("raw sample-status reconstruction mismatch")
    _verify_observations(root, rows)
    context = _independent_vela_context(roots)
    metrics, classifications, replacement = recompute_science(
        rows, context["mesh_by_topology"], context["parameters"],
        context["thermal_voltage_V"], bundle.discovery_keys,
        AcceptanceThresholds(),
    )
    reported_metrics = payload.get("candidate_metrics")
    metric_key = lambda row: (row["candidate"], row["split"], row["quantity"])
    if (not isinstance(reported_metrics, list)
            or sorted(reported_metrics, key=metric_key) != sorted(metrics, key=metric_key)):
        raise ValueError("independent typed candidate metric reconstruction mismatch")
    reported_classifications = payload.get("classifications")
    classification_key = lambda row: row["candidate"]
    if (not isinstance(reported_classifications, list)
            or sorted(reported_classifications, key=classification_key)
            != sorted(classifications, key=classification_key)):
        raise ValueError("independent typed candidate classification reconstruction mismatch")
    metrics = reported_metrics
    classifications = reported_classifications
    if _read_csv(root / "candidate_metrics.csv", CANDIDATE_COLUMNS) != _candidate_csv(metrics):
        raise ValueError("candidate_metrics.csv raw semantic mismatch")
    persisted = _load_json(root / "candidate_classifications.json")
    if persisted != {"classifications": classifications}:
        raise ValueError("candidate_classifications.json raw semantic mismatch")
    reported_replacement = payload.get("replacement_closure")
    if not isinstance(reported_replacement, list) or len(reported_replacement) != 1:
        raise ValueError("authoritative typed replacement payload mismatch")
    if payload.get("replacement_closure") != [replacement]:
        raise ValueError("independent replacement evidence reconstruction mismatch")
    if _read_csv(root / "replacement_matrix.csv", REPLACEMENT_COLUMNS) != _replacement_csv(replacement):
        raise ValueError("replacement_matrix.csv raw semantic mismatch")
    expected_context = {
        "provenance": context["provenance"],
        "van_overstraeten_parameters": context["parameters"],
        "thermal_voltage_V": context["thermal_voltage_V"],
    }
    reported_context = payload["localization_control"]["input_provenance"].get(
        "diagnostic_context"
    )
    if reported_context != expected_context:
        raise ValueError("Vela production parameter reconstruction mismatch")
    for carrier in ("electron", "hole"):
        for branch, field in (
            ("low", context["parameters"]["switch_field_V_m"] * 0.5),
            ("high", context["parameters"]["switch_field_V_m"] * 2.0),
        ):
            alpha, selected_branch = _independent_forward_van_overstraeten(
                field, context["parameters"], carrier
            )
            prefactor, critical_field = context["parameters"][carrier][branch]
            recovered, status = _independent_van_overstraeten_alpha(
                alpha, prefactor=prefactor, critical_field=critical_field,
                gamma=context["parameters"]["gamma"], branch=branch,
                switch_field=context["parameters"]["switch_field_V_m"],
            )
            if (
                selected_branch != branch or status is not SampleStatus.VALID
                or recovered is None
                or not math.isclose(recovered, field, rel_tol=1.0e-12)
            ):
                raise ValueError("Vela production Van Overstraeten branch replay mismatch")
    expected_replay = _expected_semantic(rows, replacement, context)
    if payload["localization_control"].get("semantic_replay") != expected_replay:
        raise ValueError("raw semantic replay mismatch")
    expected_mirror = _expected_mirror_invariance(rows)
    if expected_mirror["status"] != "pass":
        raise ValueError("raw mirror invariance failed")
    if payload["localization_control"].get("mirror_invariance") != expected_mirror:
        raise ValueError("authoritative mirror invariance reconstruction mismatch")
    if payload.get("sentaurus_version") != _sentaurus_version(roots):
        raise ValueError("raw Sentaurus version reconstruction mismatch")
    expected_input_manifest = {
        "schema": "vela.pn2d_minimal6_inverse_inputs.v1",
        "common_keys": [[a, b] for a, b in bundle.common_keys],
        "discovery_keys": [[a, b] for a, b in bundle.discovery_keys],
        "holdout_keys": [[a, b] for a, b in bundle.holdout_keys],
        "field_inventory": field_inventory(bundle),
        "executable_sha256": list(bundle.executable_hashes),
        "tracked_source_sha256": dict(bundle.tracked_source_hashes),
        "input_sha256": dict(bundle.input_hashes),
    }
    if _load_json(root / "input_manifest.json") != expected_input_manifest:
        raise ValueError("input_manifest.json raw semantic mismatch")
    if payload.get("input_manifest_sha256") != _sha256(root / "input_manifest.json"):
        raise ValueError("authoritative input-manifest binding mismatch")
    return ["raw_observations", "raw_candidate_metrics", "raw_classifications",
            "raw_replacement_sequences", "raw_triangle_gradient",
            "raw_current_inverted_gradient", "raw_inverse_alpha",
            "raw_local_generation", "raw_generation_support_integral",
            "mirror_invariance"]


def _verify_figures(root: Path) -> list[str]:
    manifest = _load_json(root / "figure_manifest.json")
    entries = manifest.get("figures") if isinstance(manifest, dict) else None
    if not isinstance(entries, list) or [item.get("name") for item in entries] != list(FIGURES):
        raise ValueError("fixed figure set mismatch")
    renderer = manifest.get("renderer", {})
    if renderer.get("backend") != "Pillow PNG + ReportLab invariant PDF" or "Agg" in json.dumps(renderer):
        raise ValueError("figure renderer metadata mismatch")
    required = {"question", "takeaway", "family", "variant", "row_grain_sufficiency",
                "fields", "palette_policy", "output_paths", "qa_surface"}
    for item in entries:
        name = item["name"]
        png, pdf = root / "figures" / f"{name}.png", root / "figures" / f"{name}.pdf"
        if _sha256(png) != item["png_sha256"] or _sha256(pdf) != item["pdf_sha256"]:
            raise ValueError("figure file hash mismatch")
        if _pixel_sha256(png) != item["png_pixel_sha256"]:
            raise ValueError("figure decoded-pixel hash mismatch")
        if set(item["chart_contract"]) != required:
            raise ValueError("figure chart-contract mismatch")
    return ["five_figure_pairs", "png_pixel_hashes", "truthful_renderer_metadata"]


def _write_or_validate(path: Path, payload: object, label: str) -> None:
    encoded = _write_json_bytes(payload)
    if path.exists():
        if path.read_bytes() != encoded:
            raise ValueError(f"{label} is not byte-stable")
    else:
        path.write_bytes(encoded)


def _write_package(root: Path, raw_ledger: list[dict]) -> None:
    artifacts = {}
    for path in sorted(item for item in root.rglob("*") if item.is_file()):
        relative = path.relative_to(root).as_posix()
        if relative not in PACKAGE_EXCLUSIONS:
            artifacts[relative] = _sha256(path)
    payload = {"schema": "vela.pn2d_minimal6_inverse_package_manifest.v1",
               "exclusions": sorted(PACKAGE_EXCLUSIONS),
               "raw_inputs": raw_ledger, "artifacts": artifacts}
    _write_or_validate(root / "package_manifest.json", payload, "package manifest")


def verify_report(root: str | Path) -> dict:
    report_root = Path(root).resolve()
    if not report_root.is_dir():
        raise ValueError("report root is not a directory")
    # Validate an existing package envelope before trusting any persisted report.
    _verify_package_manifest(report_root)
    report = _load_json(report_root / "physics_inverse_audit.json")
    if not isinstance(report, dict):
        raise ValueError("authoritative report must be an object")
    bundle, roots, raw_ledger = _load_raw_bundle(report)
    _verify_package_manifest(report_root, raw_ledger)
    report_manifest = _verify_report_manifest(report_root, raw_ledger)
    rows = tuple(bundle.observations)
    checks = _verify_raw_semantics(report_root, report, bundle, roots, rows)
    checks.extend(_verify_figures(report_root))
    checks.extend(("raw_input_hashes", "report_artifact_hashes"))
    scientific_payload = {
        "candidate_metrics": report["payload"]["candidate_metrics"],
        "classifications": report["payload"]["classifications"],
        "replacement_closure": report["payload"]["replacement_closure"],
    }
    scientific_payload_sha256 = hashlib.sha256(
        _write_json_bytes(scientific_payload)
    ).hexdigest()
    result = {
        "schema": "vela.pn2d_minimal6_inverse_verification.v1", "passed": True,
        "checks": sorted(checks),
        "report_manifest_sha256": _sha256(report_root / "report_manifest.json"),
        "verified_input_count": len(raw_ledger),
        "verified_artifact_count": len(report_manifest["artifacts"]),
        "scientific_payload_sha256": scientific_payload_sha256,
    }
    _write_or_validate(report_root / "verification.json", result, "verification result")
    _write_package(report_root, raw_ledger)
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report-root", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    result = verify_report(parse_args().report_root)
    print(json.dumps(result, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
