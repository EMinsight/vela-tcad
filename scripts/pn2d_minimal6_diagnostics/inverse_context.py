"""Ledger-bound topology and Vela production physics context."""

from __future__ import annotations

import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping


TOPOLOGIES = ("sketch", "mirror")
HEADER_RELATIVE = "source/tracked/include/vela/physics/ImpactIonizationModel.h"
DECK_PREFIX = "source/decks/"
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

MESH_RELATIVE = "source/topologies/{topology}/mesh.json"
_NUMBER = r"([0-9.eE+-]+)"
_BOLTZMANN_EV_PER_K = 8.617333262145e-5


@dataclass(frozen=True)
class DiagnosticContext:
    mesh_by_topology: dict[str, dict[str, dict[str, tuple[str, ...]]]]
    van_overstraeten_parameters: dict[str, object]
    thermal_voltage_V: float
    provenance: dict[str, object]


def _finite(value: object, label: str, *, positive: bool = False) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{label} must be finite") from error
    if not math.isfinite(result) or (positive and result <= 0.0):
        raise ValueError(f"{label} must be finite and positive")
    return result


def _sealed_path(root: Path, hashes: Mapping[str, str], relative: str) -> Path:
    if relative not in hashes:
        raise ValueError(f"Vela sealed context lacks {relative}")
    path = (root / relative).resolve()
    try:
        path.relative_to(root.resolve())
    except ValueError as error:
        raise ValueError("Vela sealed context path escapes root") from error
    if not path.is_file():
        raise ValueError(f"Vela sealed context lacks {relative}")
    return path


def _node_key(value: str) -> tuple[int, object]:
    return (0, int(value)) if value.isdecimal() else (1, value)


def _load_mesh(path: Path, topology: str) -> dict[str, dict[str, tuple[str, ...]]]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"sealed {topology} mesh is not valid JSON") from error
    if not isinstance(payload, dict) or payload.get("coordinate_unit") != "um":
        raise ValueError(f"sealed {topology} mesh coordinate contract mismatch")
    raw_nodes, raw_triangles = payload.get("nodes"), payload.get("triangles")
    if not isinstance(raw_nodes, list) or not isinstance(raw_triangles, list):
        raise ValueError(f"sealed {topology} mesh topology contract mismatch")
    coordinates: dict[str, tuple[float, float]] = {}
    for item in raw_nodes:
        if not isinstance(item, dict):
            raise ValueError(f"sealed {topology} mesh node contract mismatch")
        node = str(item.get("id"))
        if node in coordinates:
            raise ValueError(f"sealed {topology} mesh has duplicate node")
        coordinates[node] = (
            _finite(item.get("x"), f"{topology} mesh x"),
            _finite(item.get("y"), f"{topology} mesh y"),
        )
    triangles: dict[str, tuple[str, ...]] = {}
    for item in raw_triangles:
        if not isinstance(item, dict) or not isinstance(item.get("node_ids"), list):
            raise ValueError(f"sealed {topology} mesh triangle contract mismatch")
        cell = str(item.get("id"))
        nodes = tuple(str(node) for node in item["node_ids"])
        if cell in triangles or len(nodes) != 3 or len(set(nodes)) != 3:
            raise ValueError(f"sealed {topology} mesh triangle contract mismatch")
        if any(node not in coordinates for node in nodes):
            raise ValueError(f"sealed {topology} mesh triangle references unknown node")
        points = tuple(coordinates[node] for node in nodes)
        twice_area = (
            (points[1][0] - points[0][0]) * (points[2][1] - points[0][1])
            - (points[2][0] - points[0][0]) * (points[1][1] - points[0][1])
        )
        if not math.isfinite(twice_area) or abs(twice_area) <= 1.0e-300:
            raise ValueError(f"sealed {topology} mesh has degenerate triangle")
        triangles[cell] = nodes
    if not triangles:
        raise ValueError(f"sealed {topology} mesh has no triangles")
    edge_pairs = {
        tuple(sorted((nodes[index], nodes[(index + 1) % 3]), key=_node_key))
        for nodes in triangles.values()
        for index in range(3)
    }
    edges = {
        str(index): pair
        for index, pair in enumerate(sorted(edge_pairs, key=lambda pair: (
            _node_key(pair[0]), _node_key(pair[1])
        )))
    }
    return {
        "triangles": dict(sorted(triangles.items(), key=lambda item: _node_key(item[0]))),
        "edges": edges,
    }


def _header_string(text: str, member: str) -> str:
    matches = re.findall(
        rf'\bstd::string\s+{re.escape(member)}\s*=\s*"([^"]*)"\s*;', text
    )
    if len(matches) != 1:
        raise ValueError(f"Vela production header lacks unique string {member}")
    return matches[0]


def _header_number(text: str, member: str) -> float:
    matches = re.findall(
        rf"\bReal\s+{re.escape(member)}\s*=\s*{_NUMBER}\s*;", text
    )
    if len(matches) != 1:
        raise ValueError(f"Vela production header lacks unique numeric {member}")
    return _finite(matches[0], f"Vela production {member}")


def _header_bool(text: str, member: str) -> bool:
    matches = re.findall(
        rf"\bbool\s+{re.escape(member)}\s*=\s*(true|false)\s*;", text
    )
    if len(matches) != 1:
        raise ValueError(f"Vela production header lacks unique bool {member}")
    return matches[0] == "true"


def _validate_cross_field_config(config: dict[str, object], label: str) -> None:
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


def _validate_effective_config(config: dict[str, object], label: str) -> None:
    expected = (
        {key for _, key in _IMPACT_STRING_MEMBERS}
        | {key for _, key in _IMPACT_NUMBER_MEMBERS}
        | {key for _, key in _IMPACT_BOOL_MEMBERS}
    )
    if set(config) != expected or len(config) != 38:
        raise ValueError(
            f"{label} does not contain the complete canonical effective config"
        )
    for key, allowed in _ALLOWED_STRINGS.items():
        value = config[key]
        if not isinstance(value, str) or value not in allowed:
            raise ValueError(f"{label} has unsupported {key}")
    for _, key in _IMPACT_NUMBER_MEMBERS:
        value = config[key]
        if isinstance(value, bool):
            raise ValueError(f"{label} has invalid numeric {key}")
        number = _finite(
            value, f"{label} {key}", positive=key in _POSITIVE_CONFIG_KEYS
        )
        if key not in _POSITIVE_CONFIG_KEYS and number < 0.0:
            raise ValueError(f"{label} {key} must be non-negative")
        config[key] = number
    volume = float(config["source_volume_factor"])
    if volume != 0.0 and not 0.5 <= volume <= 1.0:
        raise ValueError(
            f"{label} source_volume_factor must be 0 or within [0.5, 1.0]"
        )
    if type(config["debug_raw_vanoverstraeten"]) is not bool:
        raise ValueError(f"{label} has invalid debug_raw_vanoverstraeten")
    _validate_cross_field_config(config, label)



def _load_header_effective_defaults(path: Path) -> dict[str, object]:
    text = path.read_text(encoding="utf-8")
    config: dict[str, object] = {
        key: _header_string(text, member)
        for member, key in _IMPACT_STRING_MEMBERS
    }
    config.update({
        key: _header_number(text, member)
        for member, key in _IMPACT_NUMBER_MEMBERS
    })
    config.update({
        key: _header_bool(text, member)
        for member, key in _IMPACT_BOOL_MEMBERS
    })
    if len(config) != 38:
        raise ValueError("Vela production header default config is incomplete")
    return config


def _impact_entries(value: object) -> list[object]:
    entries: list[object] = []
    if isinstance(value, dict):
        for key, child in value.items():
            if key == "impact_ionization":
                entries.append(child)
            else:
                entries.extend(_impact_entries(child))
    elif isinstance(value, list):
        for child in value:
            entries.extend(_impact_entries(child))
    return entries


def _load_deck_effective_config(
    path: Path,
    relative: str,
    defaults: dict[str, object],
) -> dict[str, object]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"sealed Vela deck {relative} is not valid JSON") from error
    entries = _impact_entries(payload)
    if len(entries) != 1:
        raise ValueError(
            f"sealed Vela deck {relative} must declare one impact_ionization config"
        )
    raw = entries[0]
    if isinstance(raw, str):
        overrides: dict[str, object] = {"model": raw}
    elif isinstance(raw, dict):
        unknown = set(raw) - _IMPACT_DECK_KEYS
        if unknown:
            raise ValueError(
                f"sealed Vela deck {relative} has unknown impact override "
                f"{sorted(unknown)[0]}"
            )
        overrides = {}
        for key in _DIRECT_STRING_KEYS:
            if key in raw:
                if not isinstance(raw[key], str):
                    raise ValueError(
                        f"sealed Vela deck {relative} has invalid {key}"
                    )
                overrides[key] = raw[key]
        for key in _DIRECT_NUMBER_KEYS:
            if key in raw:
                if isinstance(raw[key], bool):
                    raise ValueError(
                        f"sealed Vela deck {relative} has invalid {key}"
                    )
                overrides[key] = _finite(
                    raw[key], f"sealed Vela deck {relative} {key}"
                )
        for key in _DIRECT_BOOL_KEYS:
            if key in raw:
                if type(raw[key]) is not bool:
                    raise ValueError(
                        f"sealed Vela deck {relative} has invalid {key}"
                    )
                overrides[key] = raw[key]
        if "quasi_fermi_carrier_trucation" in raw:
            if "quasi_fermi_carrier_truncation" in raw:
                raise ValueError(
                    f"sealed Vela deck {relative} declares both truncation spellings"
                )
            value = raw["quasi_fermi_carrier_trucation"]
            if isinstance(value, bool):
                raise ValueError(
                    f"sealed Vela deck {relative} has invalid "
                    "quasi_fermi_carrier_trucation"
                )
            overrides["quasi_fermi_carrier_truncation"] = _finite(
                value,
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
                        f"sealed Vela deck {relative} has unknown interpolation "
                        f"override {sorted(unknown_nested)[0]}"
                    )
                if "mode" in interpolation:
                    if not isinstance(interpolation["mode"], str):
                        raise ValueError(
                            f"sealed Vela deck {relative} has invalid interpolation mode"
                        )
                    overrides["driving_force_interpolation"] = interpolation["mode"]
                for source, target in _NESTED_REF_KEYS.items():
                    if source in interpolation:
                        value = interpolation[source]
                        if isinstance(value, bool):
                            raise ValueError(
                                f"sealed Vela deck {relative} has invalid {source}"
                            )
                        overrides[target] = _finite(
                            value, f"sealed Vela deck {relative} {source}"
                        )
            else:
                raise ValueError(
                    f"sealed Vela deck {relative} has invalid "
                    "driving_force_interpolation"
                )
    else:
        raise ValueError(
            f"sealed Vela deck {relative} impact_ionization is invalid"
        )
    effective = dict(defaults)
    effective.update(overrides)
    _validate_effective_config(effective, f"sealed Vela deck {relative}")
    return effective


def _parameters_from_effective(
    config: dict[str, object],
) -> tuple[dict[str, object], float]:
    reference = float(config["reference_temperature_K"])
    temperature = float(config["temperature_K"])
    phonon = float(config["phonon_energy_eV"])
    gamma = (
        math.tanh(phonon / (2.0 * _BOLTZMANN_EV_PER_K * reference))
        / math.tanh(phonon / (2.0 * _BOLTZMANN_EV_PER_K * temperature))
    )
    a_scale = float(config["A_scale"])
    b_scale = float(config["B_scale"])
    parameters: dict[str, object] = {
        "gamma": gamma,
        "switch_field_V_m": float(config["switch_field_V_m"]),
        "electron": {
            "low": (
                a_scale * float(config["electron_a_low_m_inv"]),
                b_scale * float(config["electron_b_low_V_m"]),
            ),
            "high": (
                a_scale * float(config["electron_a_high_m_inv"]),
                b_scale * float(config["electron_b_high_V_m"]),
            ),
        },
        "hole": {
            "low": (
                a_scale * float(config["hole_a_low_m_inv"]),
                b_scale * float(config["hole_b_low_V_m"]),
            ),
            "high": (
                a_scale * float(config["hole_a_high_m_inv"]),
                b_scale * float(config["hole_b_high_V_m"]),
            ),
        },
        "phonon_energy_eV": phonon,
        "reference_temperature_K": reference,
        "temperature_K": temperature,
    }
    return parameters, _BOLTZMANN_EV_PER_K * temperature


def load_sealed_vela_context(
    root: str | Path, member_hashes: Mapping[str, str],
) -> DiagnosticContext:
    base = Path(root).resolve()
    mesh = {}
    mesh_provenance = {}
    for topology in TOPOLOGIES:
        relative = MESH_RELATIVE.format(topology=topology)
        path = _sealed_path(base, member_hashes, relative)
        mesh[topology] = _load_mesh(path, topology)
        mesh_provenance[topology] = {
            "logical_id": f"vela:{relative}", "sha256": member_hashes[relative],
        }
    header_path = _sealed_path(base, member_hashes, HEADER_RELATIVE)
    header_defaults = _load_header_effective_defaults(header_path)
    deck_sources: list[dict[str, object]] = []
    effective_configs: list[dict[str, object]] = []
    deck_relatives = sorted(
        relative for relative in member_hashes
        if relative.startswith(DECK_PREFIX) and relative.endswith(".json")
    )
    if len(deck_relatives) != 40:
        raise ValueError(
            "Vela sealed context requires exactly 40 hash-bound simulation decks"
        )
    for relative in deck_relatives:
        path = _sealed_path(base, member_hashes, relative)
        effective = _load_deck_effective_config(
            path, relative, header_defaults
        )
        effective_configs.append(effective)
        deck_sources.append({
            "logical_id": f"vela:{relative}",
            "sha256": member_hashes[relative],
            "effective_config": dict(effective),
        })
    effective_config = effective_configs[0]
    if any(config != effective_config for config in effective_configs[1:]):
        raise ValueError(
            "sealed Vela decks have inconsistent effective impact config"
        )
    parameters, thermal_voltage = _parameters_from_effective(effective_config)
    provenance = {
        "parameter_identity": "sealed_vela_production_defaults",
        "parameter_source": f"vela:{HEADER_RELATIVE}",
        "parameter_sha256": member_hashes[HEADER_RELATIVE],
        "mesh_sources": mesh_provenance,
        "deck_sources": deck_sources,
        "effective_impact_ionization_config": dict(effective_config),
    }
    return DiagnosticContext(mesh, parameters, thermal_voltage, provenance)
