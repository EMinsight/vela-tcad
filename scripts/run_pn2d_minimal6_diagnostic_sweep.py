#!/usr/bin/env python3
"""Create immutable, segmented Minimal6 Vela/Sentaurus diagnostic-sweep packages."""
from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import math
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))
from scripts.pn2d_minimal6_diagnostics.counterfactual import (
    integrate_native_nodal_per_unit_depth,
    sentaurus_alpha_current_nodal,
)
from scripts.pn2d_minimal6_diagnostics.schemas import validate_sweep_manifest_v1
DISCLAIMER = "minimal6 diagnostic sweep; not a physical BV curve"
SCHEMA = "vela.pn2d_minimal6_sweep_manifest.v1"
VELA_SOURCE_AREA_CM2_PER_UM2 = 1.0e-8
BRANCH_THRESHOLD_VERSION = "v1: multiplication=[0.1,10], leakage<=1e-3"
GEOMETRIC_ZERO_SOURCE_INTEGRAL = 1.0e-285


def integer_targets() -> tuple[float, ...]:
    return tuple(float(-value) for value in range(21))


def classify_branch(
    sentaurus_current_A_per_um: float | None,
    vela_current_A_per_um: float | None,
    *,
    geometric_zero: bool = False,
) -> str:
    if geometric_zero or sentaurus_current_A_per_um is None or vela_current_A_per_um is None or sentaurus_current_A_per_um == 0.0:
        return "unidentified"
    ratio = abs(vela_current_A_per_um / sentaurus_current_A_per_um)
    if 0.1 <= ratio <= 10.0:
        return "multiplication_like"
    if ratio <= 1.0e-3:
        return "leakage_like"
    return "unidentified"


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _strip_allowed(deck: dict[str, Any]) -> dict[str, Any]:
    value = copy.deepcopy(deck)
    for key in ("mesh_file", "node_doping_file"):
        value.pop(key, None)
    sweep = value.get("sweep", {})
    for key in ("start", "stop", "initial_state_file", "write_state_every_point_prefix", "csv_file"):
        sweep.pop(key, None)
    for cfg in sweep.get("diagnostics", {}).values():
        if isinstance(cfg, dict):
            cfg.pop("csv_file", None)
            cfg.pop("summary_file", None)
    return value


def validate_segment_deck(template: dict[str, Any], generated: dict[str, Any]) -> None:
    if _strip_allowed(template) != _strip_allowed(generated):
        raise ValueError("generated segment changed immutable physics or solver configuration")
    for key in ("mesh_file", "node_doping_file"):
        if not isinstance(generated.get(key), str) or not generated[key]:
            raise ValueError(f"segment deck lacks topology input {key}")
    sweep = generated.get("sweep")
    if not isinstance(sweep, dict) or not all(math.isfinite(float(sweep[key])) for key in ("start", "stop")):
        raise ValueError("segment deck lacks finite start/stop")
    if sweep["stop"] != sweep["start"] - 1.0:
        raise ValueError("diagnostic segments must end at the next exact integer volt")
    if not isinstance(sweep.get("write_state_every_point_prefix"), str) or not sweep["write_state_every_point_prefix"]:
        raise ValueError("segment deck lacks unique state prefix")


def bias_token(bias: float) -> str:
    return ("m" if bias < 0.0 else "") + f"{abs(bias):.6f}".replace(".", "p")


def segment_state_path(root: Path, topology: str, segment_start: float, target_bias: float) -> Path:
    return root / "vela" / topology / "states" / f"segment_{abs(int(segment_start)):02d}_bias_{bias_token(target_bias)}.csv"

def make_segment_deck(template: dict[str, Any], *, topology: str, segment_start: float, root: Path, restart: Path | None) -> dict[str, Any]:
    deck = copy.deepcopy(template)
    deck["mesh_file"] = str(root / "inputs" / topology / "mesh.json")
    deck["node_doping_file"] = str(root / "inputs" / topology / "doping.csv")
    sweep = deck["sweep"]
    sweep.update({"start": segment_start, "stop": segment_start - 1.0,
                  "csv_file": str(root / "vela" / topology / f"segment_{abs(int(segment_start)):02d}.csv"),
                  "write_state_every_point_prefix": str(root / "vela" / topology / "states" / f"segment_{abs(int(segment_start)):02d}")})
    diagnostics = sweep.get("diagnostics", {})
    for name, cfg in diagnostics.items():
        if isinstance(cfg, dict) and "csv_file" in cfg:
            cfg["csv_file"] = str(root / "vela" / topology / "diagnostics" / f"segment_{abs(int(segment_start)):02d}_{name}.csv")
    if restart is None:
        sweep.pop("initial_state_file", None)
    else:
        sweep["initial_state_file"] = str(restart)
    validate_segment_deck(template, deck)
    return deck


def read_vela_convergence_metadata(
    curve_csv: Path, target_bias_V: float,
) -> dict[str, Any] | None:
    """Read typed convergence evidence from one exact Vela curve endpoint."""
    with curve_csv.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {
            "bias_V",
            "converged",
            "iterations",
            "newton_iterations",
            "step_diagnostics",
        }
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            return None
        rows = [
            row
            for row in reader
            if abs(float(row["bias_V"]) - target_bias_V) <= 1.0e-12
        ]
    if len(rows) != 1:
        raise ValueError(
            "Vela convergence curve lacks exactly one target-bias endpoint"
        )
    row = rows[0]

    def finite_integer(name: str) -> int:
        value = float(row[name])
        if not math.isfinite(value) or not value.is_integer() or value < 0.0:
            raise ValueError(f"Vela convergence {name} must be a finite integer")
        return int(value)

    converged_token = row["converged"].strip().lower()
    if converged_token in {"1", "true"}:
        converged = True
    elif converged_token in {"0", "false"}:
        converged = False
    else:
        raise ValueError("Vela convergence flag must be boolean")

    metadata: dict[str, Any] = {
        "converged": converged,
        "iterations": finite_integer("iterations"),
        "newton_iterations": finite_integer("newton_iterations"),
        "step_diagnostics": row["step_diagnostics"],
    }
    for name in (
        "validation_diagnostics",
        "failure_reason",
        "newton_failure_class",
        "newton_failure_diagnostics_json",
    ):
        if name in row and row[name] not in (None, ""):
            metadata[name] = row[name]
    return metadata


def read_vela_endpoint(
    curve_csv: Path, terminal_csv: Path, edge_source_csv: Path, target_bias_V: float,
) -> dict[str, float]:
    """Read exact endpoint observables and independently close Vela's SG source."""
    with curve_csv.open(newline="", encoding="utf-8") as handle:
        curve_rows = [row for row in csv.DictReader(handle) if abs(float(row["bias_V"]) - target_bias_V) <= 1.0e-12]
    if len(curve_rows) != 1:
        raise ValueError("Vela curve lacks exactly one target-bias endpoint")
    with terminal_csv.open(newline="", encoding="utf-8") as handle:
        terminal_rows = [row for row in csv.DictReader(handle) if abs(float(row["bias_V"]) - target_bias_V) <= 1.0e-12]
    contacts = {row["contact"]: row for row in terminal_rows}
    if set(contacts) != {"Anode", "Cathode"}:
        raise ValueError("Vela terminal diagnostics lack both contacts at the endpoint")
    with edge_source_csv.open(newline="", encoding="utf-8") as handle:
        edge_rows = [
            row for row in csv.DictReader(handle)
            if abs(float(row["bias_V"]) - target_bias_V) <= 1.0e-12
        ]
    if not edge_rows:
        raise ValueError("Vela edge diagnostics lack the target-bias endpoint")
    native_source = float(contacts["Anode"]["sg_avalanche_source_integral_total"])
    reconstructed_source = sum(float(row["edge_source_integral"]) for row in edge_rows)
    if not math.isclose(native_source, reconstructed_source, rel_tol=1.0e-12, abs_tol=1.0e-285):
        raise ValueError(
            "Vela native/reconstructed source closure exceeds tolerance: "
            f"native={native_source:.17g}, reconstructed={reconstructed_source:.17g}"
        )
    # Native diagnostics combine cm^-1 and cm^-2 s^-1 values with mesh areas
    # stored numerically in um^2. Convert um^2 to cm^2 for the published source.
    values = {
        "anode_current_A_per_um": float(contacts["Anode"]["I_sgflux_A_per_um"]),
        "cathode_current_A_per_um": float(contacts["Cathode"]["I_sgflux_A_per_um"]),
        "max_field_V_per_m": float(curve_rows[0]["max_electric_field_V_per_m"]),
        "native_source_integral_s_inv_per_cm": native_source * VELA_SOURCE_AREA_CM2_PER_UM2,
        "reconstructed_source_integral_s_inv_per_cm": reconstructed_source * VELA_SOURCE_AREA_CM2_PER_UM2,
    }
    if not all(math.isfinite(value) for value in values.values()):
        raise ValueError("Vela endpoint contains non-finite observables")
    return values

def _read_export_scalar(path: Path) -> dict[int, float]:
    with path.open(newline="", encoding="utf-8") as handle:
        return {int(row["node_id"]): float(row["component0"]) for row in csv.DictReader(handle)}


def _read_export_vector(path: Path) -> dict[int, tuple[float, float]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return {int(row["node_id"]): (float(row["component0"]), float(row["component1"])) for row in csv.DictReader(handle)}


def _required_field_region(fields: list[dict[str, Any]], name: str, *, components: int, unit: str, region_name: str = "R.Si") -> int:
    matches = [
        field for field in fields
        if field.get("name") == name and field.get("region_name") == region_name
        and field.get("components") == components and field.get("unit") == unit
    ]
    if len(matches) != 1:
        raise ValueError(
            f"Sentaurus export lacks exactly one contract-compatible {name} field for {region_name}"
        )
    return int(matches[0]["region"])


def _contact_field_region(fields: list[dict[str, Any]], name: str, contact: str, unit: str) -> int:
    matches = [field for field in fields if field.get("name") == name and field.get("region_name") == contact]
    if len(matches) != 1 or matches[0].get("components") != 1 or matches[0].get("unit") != unit:
        raise ValueError(f"Sentaurus export lacks exactly one {name} field for {contact}")
    return int(matches[0]["region"])


def _single_contact_scalar(export_dir: Path, name: str, region: int) -> float:
    values = _read_export_scalar(export_dir / "fields" / f"{name}_region{region}.csv")
    if len(values) != 1:
        raise ValueError(f"Sentaurus contact field {name} must contain one scalar")
    value = next(iter(values.values()))
    if not math.isfinite(value):
        raise ValueError(f"Sentaurus contact field {name} is non-finite")
    return value


def read_sentaurus_endpoint(export_dir: Path, mesh_path: Path, target_bias_V: float) -> dict[str, Any]:
    """Recover exact, unit-normalized endpoint observables from one imported Sentaurus TDR."""
    field_manifest = export_dir / "field_manifest.json"
    if not field_manifest.is_file() or not mesh_path.is_file():
        raise FileNotFoundError("Sentaurus endpoint requires exported fields and the topology mesh")
    fields = json.loads(field_manifest.read_text(encoding="utf-8")).get("fields", [])
    requirements = {
        "ElectricField": (2, "V*cm^-1"), "ImpactIonization": (1, "cm^-3*s^-1"),
        "eAlphaAvalanche": (1, "cm^-1"), "hAlphaAvalanche": (1, "cm^-1"),
        "eCurrentDensity": (2, "A*cm^-2"), "hCurrentDensity": (2, "A*cm^-2"),
    }
    regions = {name: _required_field_region(fields, name, components=components, unit=unit) for name, (components, unit) in requirements.items()}
    if len(set(regions.values())) != 1:
        raise ValueError("Sentaurus volume fields do not share one semiconductor region")
    anode_current_region = _contact_field_region(fields, "ContactCurrentFlux", "Anode", "A")
    cathode_current_region = _contact_field_region(fields, "ContactCurrentFlux", "Cathode", "A")
    anode_voltage_region = _contact_field_region(fields, "ContactExternalVoltage", "Anode", "V")
    cathode_voltage_region = _contact_field_region(fields, "ContactExternalVoltage", "Cathode", "V")
    actual_bias = _single_contact_scalar(export_dir, "ContactExternalVoltage", anode_voltage_region) - _single_contact_scalar(export_dir, "ContactExternalVoltage", cathode_voltage_region)
    if abs(actual_bias - target_bias_V) > 1.0e-12:
        raise ValueError("Sentaurus endpoint does not match the exact target bias")
    region = regions["ImpactIonization"]
    field_dir = export_dir / "fields"
    mesh = json.loads(mesh_path.read_text(encoding="utf-8"))
    native = integrate_native_nodal_per_unit_depth(mesh, _read_export_scalar(field_dir / f"ImpactIonization_region{region}.csv"))
    reconstructed = integrate_native_nodal_per_unit_depth(
        mesh,
        sentaurus_alpha_current_nodal(
            _read_export_scalar(field_dir / f"eAlphaAvalanche_region{region}.csv"),
            _read_export_vector(field_dir / f"eCurrentDensity_region{region}.csv"),
            _read_export_scalar(field_dir / f"hAlphaAvalanche_region{region}.csv"),
            _read_export_vector(field_dir / f"hCurrentDensity_region{region}.csv"),
        ),
    )
    maximum_field = max(math.hypot(*value) for value in _read_export_vector(field_dir / f"ElectricField_region{region}.csv").values()) * 100.0
    observables = {
        "anode_current_A_per_um": _single_contact_scalar(export_dir, "ContactCurrentFlux", anode_current_region),
        "cathode_current_A_per_um": _single_contact_scalar(export_dir, "ContactCurrentFlux", cathode_current_region),
        "max_field_V_per_m": maximum_field,
        "native_source_integral_s_inv_per_cm": native["value_s_inv_per_unit_depth"],
        "reconstructed_source_integral_s_inv_per_cm": reconstructed["value_s_inv_per_unit_depth"],
    }
    if not all(math.isfinite(value) for value in observables.values()):
        raise ValueError("Sentaurus endpoint contains non-finite observables")
    return {"actual_bias_V": actual_bias, "observables": observables,
            "depth_convention": native["depth_convention"], "current_conversion": "Sentaurus 2-D ContactCurrentFlux A compared numerically with Vela A/um"}

def _converted_topology_mesh(source_mesh: Path) -> dict[str, Any]:
    mesh = json.loads(source_mesh.read_text(encoding="utf-8"))
    nodes = mesh.get("nodes")
    if not isinstance(nodes, list) or not nodes:
        raise ValueError("authoritative topology mesh requires a non-empty nodes list")
    for node in nodes:
        if not isinstance(node, dict) or "x" not in node or "y" not in node:
            raise ValueError("authoritative topology mesh node lacks x/y coordinates")
        x_m, y_m = float(node["x"]), float(node["y"])
        if not math.isfinite(x_m) or not math.isfinite(y_m):
            raise ValueError("authoritative topology mesh coordinates must be finite")
        node["x"] = x_m * 1.0e6
        node["y"] = y_m * 1.0e6
    mesh["coordinate_unit"] = "um"
    return mesh


def _authoritative_state_root_identity(
    authoritative_state_root: Path,
    topologies: tuple[str, ...] = ("sketch", "mirror"),
) -> dict[str, Any]:
    root = authoritative_state_root.resolve()
    manifest_path = root / "manifest.json"
    if not manifest_path.is_file():
        raise FileNotFoundError(
            f"authoritative state root manifest is missing: {manifest_path}"
        )
    topology_sources: dict[str, dict[str, str]] = {}
    for topology in topologies:
        relative_root = Path("states") / topology / "p0V" / "export"
        mesh = root / relative_root / "mesh.json"
        doping = root / relative_root / "doping.csv"
        for member in (mesh, doping):
            if not member.is_file():
                raise FileNotFoundError(
                    f"authoritative topology input is missing: {member}"
                )
        topology_sources[topology] = {
            "mesh_path": f"states/{topology}/p0V/export/mesh.json",
            "mesh_sha256": _sha(mesh),
            "doping_path": f"states/{topology}/p0V/export/doping.csv",
            "doping_sha256": _sha(doping),
        }
    return {
        "resolved_path": str(root),
        "manifest_json_sha256": _sha(manifest_path),
        "topology_sources": topology_sources,
    }


def copy_topology_inputs(authoritative_state_root: Path, destination_root: Path, topologies: tuple[str, ...] = ("sketch", "mirror")) -> dict[str, dict[str, str]]:
    """Prepare canonical 0 V topology inputs for Vela TCAD-internal units.

    Sentaurus state exports serialize mesh coordinates in metres. A Vela deck
    with scaling.mode = unit_scaling interprets JSON coordinates as
    micrometres, so the mesh must be converted by 1e6 while doping.csv remains
    in its exported cm^-3 units.
    """
    hashes: dict[str, dict[str, str]] = {}
    for topology in topologies:
        source = authoritative_state_root / "states" / topology / "p0V" / "export"
        target = destination_root / "inputs" / topology
        target.mkdir(parents=True, exist_ok=True)
        source_mesh = source / "mesh.json"
        source_doping = source / "doping.csv"
        for member in (source_mesh, source_doping):
            if not member.is_file():
                raise FileNotFoundError(f"authoritative topology input is missing: {member}")

        mesh = _converted_topology_mesh(source_mesh)

        copied_mesh = target / "mesh.json"
        _write_json(copied_mesh, mesh)
        copied_doping = target / "doping.csv"
        shutil.copyfile(source_doping, copied_doping)
        hashes[topology] = {
            "mesh.json": _sha(copied_mesh),
            "doping.csv": _sha(copied_doping),
        }
    return hashes

def _accepted_branch_groups(manifest: dict[str, Any]) -> list[list[dict[str, Any]]]:
    accepted = manifest.get("accepted_checkpoints", [])
    if not isinstance(accepted, list):
        return []
    groups: dict[tuple[str, float], list[dict[str, Any]]] = {}
    for row in accepted:
        if not isinstance(row, dict) or row.get("status") != "accepted":
            continue
        identity = (str(row.get("topology")), float(row.get("target_bias_V", math.nan)))
        groups.setdefault(identity, []).append(row)
    return list(groups.values())


def _derive_common_exact_branch_evidence(
    rows: list[dict[str, Any]],
) -> tuple[str, dict[str, Any]] | None:
    if len(rows) != 2:
        return None
    by_solver = {row.get("solver"): row for row in rows}
    if len(by_solver) != 2 or set(by_solver) != {"vela", "sentaurus"}:
        return None
    vela = by_solver["vela"]
    sentaurus = by_solver["sentaurus"]
    if vela.get("topology") != sentaurus.get("topology"):
        return None
    vela_target = float(vela.get("target_bias_V", math.nan))
    sentaurus_target = float(sentaurus.get("target_bias_V", math.nan))
    if vela_target != sentaurus_target:
        return None
    for row in (vela, sentaurus):
        actual = float(row.get("actual_bias_V", math.nan))
        if not math.isfinite(actual) or abs(actual - vela_target) > 1.0e-12:
            return None

    vela_observables = vela["observables"]
    sentaurus_observables = sentaurus["observables"]
    geometric_zero = any(
        abs(float(observables["native_source_integral_s_inv_per_cm"]))
        <= GEOMETRIC_ZERO_SOURCE_INTEGRAL
        and abs(float(observables["reconstructed_source_integral_s_inv_per_cm"]))
        <= GEOMETRIC_ZERO_SOURCE_INTEGRAL
        for observables in (vela_observables, sentaurus_observables)
    )
    vela_current = float(vela_observables["anode_current_A_per_um"])
    sentaurus_current = float(sentaurus_observables["anode_current_A_per_um"])
    ratio = (
        None
        if sentaurus_current == 0.0
        else abs(vela_current / sentaurus_current)
    )
    evidence = {
        "vela_anode_current_A_per_um": vela_current,
        "sentaurus_anode_current_A_per_um": sentaurus_current,
        "absolute_vela_over_sentaurus": ratio,
        "geometric_zero": geometric_zero,
        "threshold_version": BRANCH_THRESHOLD_VERSION,
    }
    classification = classify_branch(
        sentaurus_current,
        vela_current,
        geometric_zero=geometric_zero,
    )
    return classification, evidence


def _refresh_common_exact_branch_evidence(manifest: dict[str, Any]) -> None:
    manifest.setdefault("branch_threshold_version", BRANCH_THRESHOLD_VERSION)
    for rows in _accepted_branch_groups(manifest):
        for row in rows:
            row["branch_classification"] = "unidentified"
            row["branch_threshold_version"] = BRANCH_THRESHOLD_VERSION
            row.pop("branch_ratio_evidence", None)
        derived = _derive_common_exact_branch_evidence(rows)
        if derived is None:
            continue
        classification, evidence = derived
        for row in rows:
            row["branch_classification"] = classification
            row["branch_ratio_evidence"] = dict(evidence)


def _validate_branch_ratio_evidence_shape(
    evidence: Any,
    expected: dict[str, Any],
) -> None:
    required = {
        "vela_anode_current_A_per_um",
        "sentaurus_anode_current_A_per_um",
        "absolute_vela_over_sentaurus",
        "geometric_zero",
        "threshold_version",
    }
    if not isinstance(evidence, dict) or set(evidence) != required:
        raise ValueError("common exact branch ratio evidence has invalid fields")

    def finite_number(value: Any, name: str) -> None:
        if (
            isinstance(value, bool)
            or not isinstance(value, (int, float))
            or not math.isfinite(float(value))
        ):
            raise ValueError(f"branch ratio evidence {name} must be a finite number")

    finite_number(
        evidence["vela_anode_current_A_per_um"],
        "vela_anode_current_A_per_um",
    )
    finite_number(
        evidence["sentaurus_anode_current_A_per_um"],
        "sentaurus_anode_current_A_per_um",
    )
    if type(evidence["geometric_zero"]) is not bool:
        raise ValueError("branch ratio evidence geometric_zero must be boolean")
    if (
        type(evidence["threshold_version"]) is not str
        or evidence["threshold_version"] != BRANCH_THRESHOLD_VERSION
    ):
        raise ValueError("branch ratio evidence threshold version is noncanonical")

    ratio = evidence["absolute_vela_over_sentaurus"]
    if expected["absolute_vela_over_sentaurus"] is None:
        if ratio is not None:
            raise ValueError("zero reference current requires a null branch ratio")
    else:
        finite_number(ratio, "absolute_vela_over_sentaurus")


def _validate_serialized_branch_evidence(manifest: dict[str, Any]) -> None:
    if manifest.get("branch_threshold_version") != BRANCH_THRESHOLD_VERSION:
        raise ValueError("sweep branch threshold version is noncanonical")
    for rows in _accepted_branch_groups(manifest):
        derived = _derive_common_exact_branch_evidence(rows)
        if derived is None:
            for row in rows:
                if row.get("branch_classification") != "unidentified":
                    raise ValueError("side-only accepted branch must be unidentified")
                if row.get("branch_threshold_version") != BRANCH_THRESHOLD_VERSION:
                    raise ValueError("side-only accepted branch threshold is noncanonical")
                if "branch_ratio_evidence" in row:
                    raise ValueError("side-only accepted branch must not carry ratio evidence")
            continue
        expected_classification, expected_evidence = derived
        for row in rows:
            if row.get("branch_classification") != expected_classification:
                raise ValueError("common exact branch classification does not match currents")
            if row.get("branch_threshold_version") != BRANCH_THRESHOLD_VERSION:
                raise ValueError("common exact branch threshold is noncanonical")
            evidence = row.get("branch_ratio_evidence")
            _validate_branch_ratio_evidence_shape(evidence, expected_evidence)
            if evidence != expected_evidence:
                raise ValueError("common exact branch ratio evidence is incomplete or inconsistent")


def record_transition(
    manifest: dict[str, Any], *, solver: str, topology: str, start_bias_V: float,
    target_bias_V: float, exit_code: int, actual_bias_V: float | None,
    state_path: Path | None, observables: dict[str, float] | None, diagnostics: dict[str, str] | None = None,
    incomplete_reason: str | None = None,
) -> dict[str, Any]:
    """Append immutable evidence for one segment; rejected rows carry no observables."""
    exact = actual_bias_V is not None and abs(actual_bias_V - target_bias_V) <= 1.0e-12
    accepted = exit_code == 0 and exact and state_path is not None and state_path.is_file() and observables is not None
    supplied_convergence = None if diagnostics is None else diagnostics.get("convergence_metadata")
    convergence_metadata = (
        dict(supplied_convergence)
        if isinstance(supplied_convergence, dict) and supplied_convergence
        else {"exit_code": exit_code}
    )
    row: dict[str, Any] = {
        "solver": solver, "topology": topology, "start_bias_V": start_bias_V,
        "target_bias_V": target_bias_V, "actual_bias_V": actual_bias_V,
        "exit_code": exit_code, "status": "accepted" if accepted else "rejected",
        "state_path": None if state_path is None else str(state_path),
        "state_sha256": _sha(state_path) if accepted else None,
        "observables": dict(observables) if accepted else None,
        "stdout": "" if diagnostics is None else diagnostics.get("stdout", ""),
        "stderr": "" if diagnostics is None else diagnostics.get("stderr", ""),
        "branch_classification": "unidentified",
        "branch_threshold_version": BRANCH_THRESHOLD_VERSION,
        "convergence_metadata": convergence_metadata,
        "current_extraction_formula": (
            "signed SG contact flux in A/um"
            if solver == "vela"
            else "Sentaurus 2-D ContactCurrentFlux A compared numerically with Vela A/um"
        ),
    }
    if accepted:
        required = {"anode_current_A_per_um", "cathode_current_A_per_um", "max_field_V_per_m", "native_source_integral_s_inv_per_cm", "reconstructed_source_integral_s_inv_per_cm"}
        if observables is None or set(observables) != required or not all(math.isfinite(float(value)) for value in observables.values()):
            raise ValueError("accepted checkpoint lacks complete finite observables")
        manifest.setdefault("accepted_checkpoints", []).append(row)
        _refresh_common_exact_branch_evidence(manifest)
    else:
        if incomplete_reason:
            row["incomplete_reason"] = incomplete_reason
        manifest.setdefault("failed_transitions", []).append(row)
        if manifest.get("failed_transition") is None:
            manifest["failed_transition"] = row
    return row


def record_sentaurus_checkpoint(
    manifest: dict[str, Any], *, topology: str, start_bias_V: float, target_bias_V: float,
    state_path: Path, export_dir: Path, mesh_path: Path,
    package_root: Path | None = None,
) -> dict[str, Any]:
    """Promote one imported exact-bias Sentaurus TDR into immutable sweep evidence."""
    matches = [segment for segment in manifest.get("sentaurus_segments", [])
               if segment.get("topology") == topology and float(segment.get("target_bias_V")) == target_bias_V]
    if len(matches) != 1:
        raise ValueError("Sentaurus checkpoint does not match exactly one generated segment")
    evidence_state = state_path
    canonical_state_path: Path | None = None
    if package_root is not None:
        checkpoint_value = matches[0].get("checkpoint_tdr")
        if not isinstance(checkpoint_value, str) or not checkpoint_value:
            raise ValueError("Sentaurus segment lacks canonical checkpoint path")
        canonical_state_path = Path(checkpoint_value)
        expected_checkpoint = (
            Path("sentaurus")
            / topology
            / "checkpoints"
            / f"{topology}_{bias_token(target_bias_V)}.tdr"
        )
        if canonical_state_path != expected_checkpoint:
            raise ValueError("Sentaurus checkpoint path is not canonical")
        package_root_resolved = package_root.resolve()
        evidence_state = (package_root_resolved / canonical_state_path).resolve()
        try:
            evidence_state.relative_to(package_root_resolved)
        except ValueError as error:
            raise ValueError(
                "Sentaurus checkpoint path escapes package root"
            ) from error
    try:
        if package_root is not None and state_path.resolve() != evidence_state:
            evidence_state.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(state_path, evidence_state)
        endpoint = read_sentaurus_endpoint(export_dir, mesh_path, target_bias_V)
        row = record_transition(
            manifest, solver="sentaurus", topology=topology, start_bias_V=start_bias_V,
            target_bias_V=target_bias_V, exit_code=0, actual_bias_V=endpoint["actual_bias_V"],
            state_path=evidence_state, observables=endpoint["observables"],
            diagnostics={"stdout": "", "stderr": ""},
        )
        row["export_dir"] = str(export_dir)
        row["export_field_manifest_sha256"] = _sha(export_dir / "field_manifest.json")
        row["depth_convention"] = endpoint["depth_convention"]
        row["current_conversion"] = endpoint["current_conversion"]
    except (OSError, ValueError) as error:
        row = record_transition(
            manifest, solver="sentaurus", topology=topology, start_bias_V=start_bias_V,
            target_bias_V=target_bias_V, exit_code=1, actual_bias_V=None,
            state_path=state_path if state_path.is_file() else None, observables=None,
            diagnostics={"stdout": "", "stderr": str(error)}, incomplete_reason=f"Sentaurus checkpoint export is incomplete or invalid: {error}",
        )
    if canonical_state_path is not None and row["status"] == "accepted":
        row["state_path"] = str(canonical_state_path)
    matches[0]["status"] = row["status"]
    matches[0]["state_path"] = row["state_path"]
    matches[0]["state_sha256"] = row["state_sha256"]
    return row

def run_vela_subprocess_segment(root: Path, executable: Path, segment: dict[str, Any]) -> dict[str, Any]:
    """Run one generated Vela deck and return evidence without inventing a native source."""
    deck = root / str(segment["deck"])
    completed = subprocess.run([str(executable), "--config", str(deck)], cwd=root, text=True, capture_output=True)
    target = float(segment["target_bias_V"])
    topology = str(segment["topology"])
    start = float(segment["start_bias_V"])
    state = segment_state_path(root, topology, start, target)
    diagnostics = root / "vela" / topology / "diagnostics"
    terminal = diagnostics / f"segment_{abs(int(start)):02d}_terminal_current_method_compare.csv"
    edge_source = diagnostics / f"segment_{abs(int(start)):02d}_sg_avalanche_edges.csv"
    curve = root / "vela" / topology / f"segment_{abs(int(start)):02d}.csv"
    convergence_metadata = None
    if curve.is_file():
        try:
            convergence_metadata = read_vela_convergence_metadata(curve, target)
        except (ValueError, KeyError, csv.Error) as error:
            return {
                "exit_code": completed.returncode,
                "actual_bias_V": None,
                "state_path": None,
                "observables": None,
                "stdout": completed.stdout,
                "stderr": completed.stderr,
                "convergence_metadata": None,
                "incomplete_reason": f"Vela convergence evidence invalid: {error}",
            }
    if completed.returncode != 0 or not state.is_file() or not terminal.is_file() or not edge_source.is_file() or not curve.is_file():
        return {
            "exit_code": completed.returncode,
            "actual_bias_V": None,
            "state_path": None,
            "observables": None,
            "stdout": completed.stdout,
            "stderr": completed.stderr,
            "convergence_metadata": convergence_metadata,
        }
    if convergence_metadata is not None and not convergence_metadata["converged"]:
        return {
            "exit_code": completed.returncode,
            "actual_bias_V": None,
            "state_path": None,
            "observables": None,
            "stdout": completed.stdout,
            "stderr": completed.stderr,
            "convergence_metadata": convergence_metadata,
            "incomplete_reason": "Vela exact endpoint is not converged",
        }
    try:
        endpoint = read_vela_endpoint(curve, terminal, edge_source, target)
    except (ValueError, KeyError, csv.Error) as error:
        return {
            "exit_code": completed.returncode,
            "actual_bias_V": None,
            "state_path": None,
            "observables": None,
            "stdout": completed.stdout,
            "stderr": completed.stderr,
            "convergence_metadata": convergence_metadata,
            "incomplete_reason": f"Vela endpoint evidence invalid: {error}",
        }
    return {"exit_code": completed.returncode, "actual_bias_V": target, "state_path": state,
            "observables": endpoint, "stdout": completed.stdout, "stderr": completed.stderr,
            "convergence_metadata": convergence_metadata}

def execute_segments(
    manifest: dict[str, Any],
    root: Path,
    runner,
    *,
    package_local_state_paths: bool = False,
) -> None:
    """Run each topology independently and retain its first failed segment."""
    segments = manifest.get("segments", [])
    topologies = list(dict.fromkeys(str(segment["topology"]) for segment in segments))
    for topology in topologies:
        for segment in segments:
            if str(segment["topology"]) != topology or segment.get("status") != "pending":
                continue
            result = runner(segment)
            row = record_transition(
                manifest, solver=str(segment["solver"]), topology=topology,
                start_bias_V=float(segment["start_bias_V"]), target_bias_V=float(segment["target_bias_V"]),
                exit_code=int(result.get("exit_code", 1)), actual_bias_V=result.get("actual_bias_V"),
                state_path=result.get("state_path"), observables=result.get("observables"),
                diagnostics={
                    "stdout": str(result.get("stdout", "")),
                    "stderr": str(result.get("stderr", "")),
                    "convergence_metadata": result.get("convergence_metadata"),
                },
                incomplete_reason=None if result.get("incomplete_reason") is None else str(result["incomplete_reason"]),
            )
            if package_local_state_paths and row["status"] == "accepted":
                expected_state = segment_state_path(
                    root,
                    topology,
                    float(segment["start_bias_V"]),
                    float(segment["target_bias_V"]),
                )
                actual_state = Path(str(result["state_path"])).resolve()
                if actual_state != expected_state.resolve():
                    raise ValueError("Vela accepted state path is noncanonical")
                row["state_path"] = str(expected_state.relative_to(root))
            segment["status"] = row["status"]
            if row["status"] != "accepted":
                break

def validate_sweep_manifest(manifest: dict[str, Any], *, package_root: Path | None = None) -> None:
    if manifest.get("schema") != SCHEMA:
        raise ValueError("invalid sweep manifest schema")
    targets = tuple(manifest.get("targets_V", ()))
    if targets != integer_targets() and targets != (0.0, -1.0):
        raise ValueError("sweep manifest has non-canonical targets")
    if manifest.get("interpolation", "forbidden") != "forbidden":
        raise ValueError("sweep manifest must forbid interpolation")
    for row in manifest.get("accepted_checkpoints", []):
        if row.get("status") != "accepted" or abs(float(row["actual_bias_V"]) - float(row["target_bias_V"])) > 1.0e-12:
            raise ValueError("accepted checkpoint is not exact")
        state = Path(str(row["state_path"]))
        if not state.is_absolute():
            if package_root is None:
                raise ValueError("relative accepted checkpoint requires package_root")
            state = package_root.resolve() / state
        if not state.is_file() or _sha(state) != row.get("state_sha256"):
            raise ValueError("accepted checkpoint state is missing or hash-tampered")
        if not isinstance(row.get("observables"), dict):
            raise ValueError("accepted checkpoint lacks observables")
    failed_rows = manifest.get("failed_transitions", [])
    if any(row.get("status") != "rejected" or row.get("observables") is not None for row in failed_rows):
        raise ValueError("failed transition must preserve no fabricated observables")
    failed = manifest.get("failed_transition")
    if failed is not None and (failed.get("status") != "rejected" or failed.get("observables") is not None):
        raise ValueError("failed transition must preserve no fabricated observables")
    if failed is not None and (not failed_rows or failed != failed_rows[0]):
        raise ValueError("failed_transition must retain the first rejected transfer")
    if package_root is None:
        return

    root = package_root.resolve()
    if targets != integer_targets():
        raise ValueError("complete sweep package requires exact targets 0..-20 V")

    def checked_path(value: Any, expected_sha256: Any, label: str) -> Path:
        if not isinstance(value, str) or not value:
            raise ValueError(f"{label} lacks an artifact path")
        path = Path(value)
        if not path.is_absolute():
            path = root / path
        path = path.resolve()
        if not path.is_file():
            raise ValueError(f"{label} artifact is missing: {path}")
        if not isinstance(expected_sha256, str) or _sha(path) != expected_sha256:
            raise ValueError(f"{label} artifact is hash-tampered: {path}")
        return path

    template_record = manifest.get("template")
    if not isinstance(template_record, dict):
        raise ValueError("sweep manifest lacks immutable template provenance")
    template_path = checked_path(
        template_record.get("path"), template_record.get("sha256"), "template"
    )
    template = json.loads(template_path.read_text(encoding="utf-8"))

    authoritative = manifest.get("authoritative_state_root")
    if not isinstance(authoritative, dict):
        raise ValueError("complete package lacks authoritative state root provenance")
    source_root_value = authoritative.get("resolved_path")
    if not isinstance(source_root_value, str) or not Path(source_root_value).is_absolute():
        raise ValueError("authoritative state root path must be resolved")
    source_root = Path(source_root_value).resolve()
    expected_authoritative = _authoritative_state_root_identity(source_root)
    if authoritative != expected_authoritative:
        raise ValueError("authoritative state root provenance is missing or hash-tampered")

    topology_hashes = manifest.get("topology_input_sha256")
    if not isinstance(topology_hashes, dict) or set(topology_hashes) != {"sketch", "mirror"}:
        raise ValueError("complete package requires both topology input hash roots")
    for topology in ("sketch", "mirror"):
        hashes = topology_hashes[topology]
        if not isinstance(hashes, dict) or set(hashes) != {"mesh.json", "doping.csv"}:
            raise ValueError(f"{topology} topology hashes must cover mesh and doping")
        package_mesh = checked_path(
            str(Path("inputs") / topology / "mesh.json"),
            hashes["mesh.json"],
            f"{topology} mesh.json",
        )
        package_doping = checked_path(
            str(Path("inputs") / topology / "doping.csv"),
            hashes["doping.csv"],
            f"{topology} doping.csv",
        )
        source = expected_authoritative["topology_sources"][topology]
        source_mesh = source_root / Path(source["mesh_path"])
        source_doping = source_root / Path(source["doping_path"])
        if json.loads(package_mesh.read_text(encoding="utf-8")) != _converted_topology_mesh(source_mesh):
            raise ValueError(
                f"{topology} package mesh does not match authoritative reconstruction"
            )
        if package_doping.read_bytes() != source_doping.read_bytes():
            raise ValueError(
                f"{topology} package doping does not match authoritative source"
            )

    expected_vela = [
        (topology, start_bias, start_bias - 1.0)
        for topology in ("sketch", "mirror")
        for start_bias in integer_targets()[:-1]
    ]
    segments = manifest.get("segments")
    if not isinstance(segments, list):
        raise ValueError("Vela segments must be an array")
    actual_vela = [
        (
            row.get("topology"),
            float(row.get("start_bias_V", math.nan)),
            float(row.get("target_bias_V", math.nan)),
        )
        for row in segments
        if isinstance(row, dict)
    ]
    if len(actual_vela) != len(segments) or actual_vela != expected_vela:
        raise ValueError("Vela segments must contain each canonical transition exactly once")

    state_prefixes: set[str] = set()
    for index, (segment, identity) in enumerate(zip(segments, expected_vela, strict=True)):
        topology, start_bias, target_bias = identity
        if segment.get("solver") != "vela":
            raise ValueError("Vela segment has incorrect solver identity")
        expected_deck_path = root / "vela" / topology / "decks" / f"segment_{abs(int(start_bias)):02d}.json"
        if Path(str(segment.get("deck"))) != expected_deck_path.relative_to(root):
            raise ValueError("Vela segment deck path is noncanonical")
        deck_path = checked_path(
            segment.get("deck"), segment.get("deck_sha256"), f"Vela segment {index}"
        )
        deck = json.loads(deck_path.read_text(encoding="utf-8"))
        validate_segment_deck(template, deck)
        expected_input_root = root / "inputs" / topology
        if (
            Path(deck["mesh_file"]) != expected_input_root / "mesh.json"
            or Path(deck["node_doping_file"]) != expected_input_root / "doping.csv"
        ):
            raise ValueError("Vela segment uses topology-mismatched mesh/doping")
        sweep = deck["sweep"]
        if float(sweep["start"]) != start_bias or float(sweep["stop"]) != target_bias:
            raise ValueError("Vela deck bias does not match segment identity")
        expected_prefix = root / "vela" / topology / "states" / f"segment_{abs(int(start_bias)):02d}"
        prefix = str(sweep["write_state_every_point_prefix"])
        if Path(prefix) != expected_prefix or prefix in state_prefixes:
            raise ValueError("Vela segment state prefixes must be canonical and unique")
        state_prefixes.add(prefix)
        restart = sweep.get("initial_state_file")
        expected_restart = None if start_bias == 0.0 else str(
            segment_state_path(root, topology, start_bias + 1.0, start_bias)
        )
        if restart != expected_restart:
            raise ValueError("Vela segment restart does not match the prior exact checkpoint")

    sentaurus_segments = manifest.get("sentaurus_segments")
    if not isinstance(sentaurus_segments, list):
        raise ValueError("Sentaurus segments must be an array")
    expected_sentaurus = [
        (topology, bias)
        for topology in ("sketch", "mirror")
        for bias in integer_targets()
    ]
    actual_sentaurus = [
        (row.get("topology"), float(row.get("target_bias_V", math.nan)))
        for row in sentaurus_segments
        if isinstance(row, dict)
    ]
    if len(actual_sentaurus) != len(sentaurus_segments) or actual_sentaurus != expected_sentaurus:
        raise ValueError("Sentaurus segments must contain each canonical checkpoint exactly once")
    sentaurus_template_path = (
        REPO
        / "reference_tcad"
        / "pn2d_sentaurus2018_minimal6"
        / "source"
        / "pn2d_minimal6_sweep_sdevice.cmd"
    )
    sentaurus_template = sentaurus_template_path.read_text(encoding="utf-8")
    for index, (segment, identity) in enumerate(
        zip(sentaurus_segments, expected_sentaurus, strict=True)
    ):
        topology, bias = identity
        if segment.get("solver") != "sentaurus":
            raise ValueError("Sentaurus segment has incorrect solver identity")
        tag = f"{topology}_{bias_token(bias)}"
        expected_deck = Path("sentaurus") / topology / "decks" / f"{tag}.cmd"
        expected_checkpoint = Path("sentaurus") / topology / "checkpoints" / f"{tag}.tdr"
        if Path(str(segment.get("deck"))) != expected_deck:
            raise ValueError("Sentaurus segment deck path is noncanonical")
        if Path(str(segment.get("checkpoint_tdr"))) != expected_checkpoint:
            raise ValueError("Sentaurus checkpoint path is noncanonical")
        deck_path = checked_path(
            segment.get("deck"), segment.get("deck_sha256"), f"Sentaurus segment {index}"
        )
        expected_payload = sentaurus_template.replace("__BIAS_TAG__", tag).replace(
            "__TARGET_BIAS_V__", f"{bias:.1f}"
        )
        if deck_path.read_text(encoding="utf-8") != expected_payload:
            raise ValueError("Sentaurus deck does not match the immutable exact payload")

    for row in manifest.get("accepted_checkpoints", []):
        state_path = Path(str(row["state_path"]))
        if state_path.is_absolute():
            raise ValueError("strict package accepted state path must be relative")
        topology = str(row["topology"])
        target_bias = float(row["target_bias_V"])
        if row["solver"] == "vela":
            expected_state = segment_state_path(
                root,
                topology,
                float(row["start_bias_V"]),
                target_bias,
            ).relative_to(root)
        elif row["solver"] == "sentaurus":
            matches = [
                segment
                for segment in sentaurus_segments
                if segment["topology"] == topology
                and float(segment["target_bias_V"]) == target_bias
            ]
            if len(matches) != 1:
                raise ValueError("accepted Sentaurus state lacks canonical segment")
            expected_state = Path(str(matches[0]["checkpoint_tdr"]))
        else:
            raise ValueError("accepted state has invalid solver")
        if state_path != expected_state:
            raise ValueError(
                f"{row['solver']} accepted state path is noncanonical"
            )

    _validate_serialized_branch_evidence(manifest)
    validate_sweep_manifest_v1(manifest)


def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, sort_keys=True, indent=2) + "\n", encoding="utf-8")


def write_sentaurus_decks(root: Path, deck_template: Path, topologies: tuple[str, ...] = ("sketch", "mirror")) -> list[dict[str, Any]]:
    """Create independent exact-bias Sentaurus decks; no authoritative manifest is touched."""
    source = deck_template.read_text(encoding="utf-8")
    rows: list[dict[str, Any]] = []
    for topology in topologies:
        for bias in integer_targets():
            tag = f"{topology}_{bias_token(bias)}"
            path = root / "sentaurus" / topology / "decks" / f"{tag}.cmd"
            payload = source.replace("__BIAS_TAG__", tag).replace("__TARGET_BIAS_V__", f"{bias:.1f}")
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(payload, encoding="utf-8")
            rows.append({"solver": "sentaurus", "topology": topology, "target_bias_V": bias,
                         "deck": str(path.relative_to(root)), "deck_sha256": _sha(path), "status": "pending",
                         "checkpoint_tdr": str((root / "sentaurus" / topology / "checkpoints" / f"{tag}.tdr").relative_to(root))})
    return rows

def initialise_package(root: Path, template_path: Path, topologies: tuple[str, ...] = ("sketch", "mirror"), authoritative_state_root: Path | None = None) -> Path:
    template = json.loads(template_path.read_text(encoding="utf-8"))
    input_hashes = {} if authoritative_state_root is None else copy_topology_inputs(authoritative_state_root, root, topologies)
    authoritative_identity = (
        None if authoritative_state_root is None
        else _authoritative_state_root_identity(authoritative_state_root, topologies)
    )
    targets = integer_targets()
    segments: list[dict[str, Any]] = []
    for topology in topologies:
        restart: Path | None = None
        for start in targets[:-1]:
            deck = make_segment_deck(template, topology=topology, segment_start=start, root=root, restart=restart)
            deck_path = root / "vela" / topology / "decks" / f"segment_{abs(int(start)):02d}.json"
            _write_json(deck_path, deck)
            segments.append({"solver": "vela", "topology": topology, "start_bias_V": start, "target_bias_V": start - 1.0,
                             "deck": str(deck_path.relative_to(root)), "deck_sha256": _sha(deck_path), "status": "pending"})
            restart = segment_state_path(root, topology, start, start - 1.0)
    sentaurus_template = REPO / "reference_tcad" / "pn2d_sentaurus2018_minimal6" / "source" / "pn2d_minimal6_sweep_sdevice.cmd"
    sentaurus_segments = write_sentaurus_decks(root, sentaurus_template, topologies)
    manifest = {"schema": SCHEMA, "diagnostic_disclaimer": DISCLAIMER, "targets_V": list(targets),
                "template": {"path": str(template_path), "sha256": _sha(template_path)}, "authoritative_state_root": authoritative_identity, "topology_input_sha256": input_hashes, "segments": segments, "sentaurus_segments": sentaurus_segments,
                "accepted_checkpoints": [], "failed_transition": None, "failed_transitions": [], "interpolation": "forbidden",
                "branch_threshold_version": BRANCH_THRESHOLD_VERSION}
    manifest_path = root / "sweep_manifest.json"
    _write_json(manifest_path, manifest)
    return manifest_path


def record_sentaurus_results(
    manifest: dict[str, Any],
    results_path: Path,
    *,
    package_root: Path | None = None,
) -> list[dict[str, Any]]:
    """Load declared imported Sentaurus checkpoints without changing solver configuration."""
    payload = json.loads(results_path.read_text(encoding="utf-8-sig"))
    if not isinstance(payload, list):
        raise ValueError("Sentaurus results JSON must be a list")
    required = {"topology", "start_bias_V", "target_bias_V", "state_path", "export_dir", "mesh_path"}
    rows: list[dict[str, Any]] = []
    for item in payload:
        if not isinstance(item, dict) or set(item) != required:
            raise ValueError("Sentaurus result has an invalid contract")
        def resolve(value: Any) -> Path:
            path = Path(str(value))
            return path if path.is_absolute() else (results_path.parent / path).resolve()
        rows.append(record_sentaurus_checkpoint(
            manifest, topology=str(item["topology"]), start_bias_V=float(item["start_bias_V"]),
            target_bias_V=float(item["target_bias_V"]), state_path=resolve(item["state_path"]),
            export_dir=resolve(item["export_dir"]), mesh_path=resolve(item["mesh_path"]),
            package_root=package_root,
        ))
    return rows

def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--template", type=Path, default=REPO / "reference_tcad" / "pn2d_sentaurus2018_minimal6" / "vela" / "pn2d_minimal6_sweep_template.json")
    parser.add_argument("--state-root", type=Path)
    parser.add_argument("--vela-runner", type=Path)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--sentaurus-results-json", type=Path)
    args = parser.parse_args()
    args.out_dir = args.out_dir.resolve()
    manifest_path = args.out_dir / "sweep_manifest.json" if args.resume else initialise_package(args.out_dir, args.template, authoritative_state_root=args.state_root)
    if args.resume and not manifest_path.is_file():
        raise FileNotFoundError("--resume requires an existing sweep manifest")
    if args.sentaurus_results_json is not None or args.vela_runner is not None:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if args.sentaurus_results_json is not None:
            record_sentaurus_results(manifest, args.sentaurus_results_json.resolve(), package_root=args.out_dir)
        if args.vela_runner is not None:
            execute_segments(manifest, args.out_dir, lambda segment: run_vela_subprocess_segment(args.out_dir, args.vela_runner, segment), package_local_state_paths=True)
        validate_sweep_manifest(manifest, package_root=args.out_dir)
        _write_json(manifest_path, manifest)
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
