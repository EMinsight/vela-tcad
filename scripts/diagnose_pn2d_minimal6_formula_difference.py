#!/usr/bin/env python3
"""Create fixed-state PN2D minimal6 formula-difference diagnostics."""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.export_pn2d_minimal6_states import validate_member_hashes, validate_sealed_archive
from scripts.pn2d_minimal6_diagnostics.physics import (
    compare_van_overstraeten_parameters,
    infer_ni_eff,
    parse_van_overstraeten_de_man,
)
from scripts.pn2d_minimal6_diagnostics.counterfactual import (
    build_adjacent_interactions,
    symmetric_contributions,
    evaluate_formula_counterfactual,
    FACTOR_DEPENDENCIES,
    assert_counterfactual_closure,
    score_dominance,
    validate_field_units,
    validate_dependency_dag,
    validate_source_anchor_kind,
    integrate_native_nodal_per_unit_depth,
    sentaurus_alpha_current_nodal,
    source_log_gap,
    validate_formula_input,
)
from scripts.pn2d_minimal6_diagnostics.contracts import SourceKind
from scripts.pn2d_minimal6_diagnostics.schemas import DISCLAIMER, validate_formula_difference_v1
import scripts.audit_pn2d_minimal6_fixed_state as fixed_state_audit
from scripts.pn2d_minimal6_diagnostics.plots import render_formula_difference_figures

_NODE_STATE_FIELDS = {
    "ElectrostaticPotential": "V", "eDensity": "cm^-3", "hDensity": "cm^-3",
    "eQuasiFermiPotential": "V", "hQuasiFermiPotential": "V",
    "eMobility": "cm^2*V^-1*s^-1", "hMobility": "cm^2*V^-1*s^-1",
    "eVelocity": "cm*s^-1", "hVelocity": "cm*s^-1",
    "eIonIntegral": "1", "hIonIntegral": "1", "MeanIonIntegral": "1",

    "LatticeTemperature": "K",
    "eAlphaAvalanche": "cm^-1", "hAlphaAvalanche": "cm^-1",
}


SOURCE_FAMILIES = {
    "sentaurus_native_avalanche_generation": SourceKind.SENTAURUS.value,
    "sentaurus_alpha_current_reconstruction": SourceKind.DERIVED.value,
    "vela_alpha_flux_partial_volume_reconstruction": SourceKind.DERIVED.value,
}
SOURCE_VALUE_KEYS = {
    "sentaurus_native_avalanche_generation":
        "sentaurus_native_ImpactIonization_s_inv_per_unit_depth",
    "sentaurus_alpha_current_reconstruction":
        "sentaurus_alpha_current_reconstruction_s_inv_per_unit_depth",
    "vela_alpha_flux_partial_volume_reconstruction":
        "vela_reconstruction_s_inv_per_unit_depth",
}


def _validate_source_families() -> None:
    if set(SOURCE_FAMILIES) != set(SOURCE_VALUE_KEYS):
        raise ValueError("source family/value mappings must cover the same families")
    for family, source_kind in SOURCE_FAMILIES.items():
        validate_source_anchor_kind(
            family, source_kind,
            native=(family == "sentaurus_native_avalanche_generation"),
        )

def _read_scalar(path: Path) -> dict[int, float]:
    with path.open(newline="", encoding="utf-8") as handle:
        return {int(row["node_id"]): float(row["component0"]) for row in csv.DictReader(handle)}


def _read_vector(path: Path) -> dict[int, tuple[float, float]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return {int(row["node_id"]): (float(row["component0"]), float(row["component1"])) for row in csv.DictReader(handle)}


def _integer_csv(value: object, field: str) -> int:
    numeric = float(value)
    integer = int(numeric)
    if not math.isfinite(numeric) or numeric != integer:
        raise ValueError(f"audit identity {field} must be an exact integer")
    return integer


def _base(record: dict, kind: str, **fields) -> dict:
    row = {
        "run_id": "minimal6_fixed_state", "record_kind": kind,
        "topology": record["topology"], "bias_V": record["bias_V"],
        "node_id": "", "cell_id": "", "edge_id": "", "quantity": "", "component": "",
        "value": "", "unit": "", "source": "", "source_kind": SourceKind.DERIVED.value, "value_s_inv_per_unit_depth": "",
        "depth_convention": "", "status": "available",
    }
    row.update(fields)
    return row


def _validate_export_units(export_dir: Path) -> None:
    manifest = export_dir / "field_manifest.json"
    if not manifest.is_file():
        return
    fields = json.loads(manifest.read_text(encoding="utf-8")).get("fields", [])
    validate_field_units(fields, {
        "ImpactIonization": "cm^-3*s^-1", "eAlphaAvalanche": "cm^-1", "hAlphaAvalanche": "cm^-1",
        "eCurrentDensity": "A*cm^-2", "hCurrentDensity": "A*cm^-2",
        "ElectrostaticPotential": "V", "eDensity": "cm^-3", "hDensity": "cm^-3",
        "eQuasiFermiPotential": "V", "hQuasiFermiPotential": "V",
        "eMobility": "cm^2*V^-1*s^-1", "hMobility": "cm^2*V^-1*s^-1",
        "eVelocity": "cm*s^-1", "hVelocity": "cm*s^-1",
        "LatticeTemperature": "K",
        "eIonIntegral": "1", "hIonIntegral": "1", "MeanIonIntegral": "1",
    })

def _vela_parameter_agreement(export_dir: Path, models_par: dict | None) -> dict:
    if models_par is None:
        return {"status": "unavailable", "reason": "Sentaurus models.par is unavailable", "comparisons": []}
    audit_path = export_dir / "audit.json"
    if not audit_path.is_file():
        return compare_van_overstraeten_parameters(models_par, None)
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    impact = audit.get("impact_ionization", {})
    return compare_van_overstraeten_parameters(models_par, impact.get("van_overstraeten_parameters"))

def _tree_hashes(root: Path) -> dict[str, str]:
    return {
        path.relative_to(root).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in sorted(candidate for candidate in root.rglob("*") if candidate.is_file())
    }


def _csv_bytes(rows: list[dict]) -> bytes:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue().encode("utf-8")


def validate_audit_binding(state_root: Path, audit_root: Path) -> dict[str, object]:
    state_root = Path(state_root).resolve()
    audit_root = Path(audit_root).resolve()
    state_manifest_path = state_root / "manifest.json"
    audit_manifest_path = audit_root / "manifest.json"
    summary_path = audit_root / "summary.json"
    if not state_manifest_path.is_file() or not audit_manifest_path.is_file() or not summary_path.is_file():
        raise FileNotFoundError("state/audit binding requires both manifests and the audit summary")
    state_manifest = json.loads(state_manifest_path.read_text(encoding="utf-8"))
    audit_manifest = json.loads(audit_manifest_path.read_text(encoding="utf-8"))
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    expected_schema = "vela.pn2d_minimal6_fixed_state_audit.v1"
    if audit_manifest.get("schema") != expected_schema:
        raise ValueError("audit manifest schema mismatch")
    if audit_manifest.get("source_state_schema") != state_manifest.get("schema"):
        raise ValueError("audit source-state schema mismatch")
    if audit_manifest.get("topology_definitions") != state_manifest.get("states"):
        raise ValueError("audit six-state association does not match state manifest")
    if audit_manifest.get("task4_provenance") != state_manifest.get("task4_provenance"):
        raise ValueError("audit replay provenance does not match state manifest")
    if audit_manifest.get("gate_status") != "PASS":
        raise ValueError("audit gate status is not PASS")
    command_status = audit_manifest.get("command_status")
    if not isinstance(command_status, dict) or command_status.get("task4_replay") != "PASS" or command_status.get("report_generation") != "PASS":
        raise ValueError("audit replay/report command status is not PASS")
    if summary.get("schema") != expected_schema or summary.get("status") != "PASS":
        raise ValueError("audit summary is not a passing fixed-state report")
    gates = summary.get("gates")
    if not isinstance(gates, dict) or gates.get("passed") is not True or gates.get("provenance_replay_validated") is not True:
        raise ValueError("audit summary gate/provenance status is not validated")
    expected_hashes = _tree_hashes(state_root)
    recorded_hashes = audit_manifest.get("input_sha256")
    if not isinstance(recorded_hashes, dict):
        raise ValueError("audit manifest lacks input_sha256 binding")
    if recorded_hashes != expected_hashes:
        manifest_hash = expected_hashes.get("manifest.json")
        if recorded_hashes.get("manifest.json") != manifest_hash:
            raise ValueError("audit state manifest hash/run binding mismatch")
        raise ValueError("audit state input_sha256 member binding mismatch")
    producer = REPO / "build-release" / "pn2d_minimal6_operator_audit.exe"
    failures = fixed_state_audit.verify_task4_replay(state_root, producer)
    if failures:
        raise ValueError("audit replay provenance/hash validation failed: " + "; ".join(failures))
    expected = fixed_state_audit.build_report(state_root)
    expected_artifacts = {
        "node_state.csv": expected.node_rows,
        "edge_audit.csv": expected.edge_rows,
        "triangle_audit.csv": expected.triangle_rows,
    }
    for filename, rows in expected_artifacts.items():
        path = audit_root / filename
        if not path.is_file() or path.read_bytes() != _csv_bytes(rows):
            raise ValueError(f"audit artifact {filename} mutation or state mismatch")
    return _load_audit_root(audit_root)


def _load_audit_root(audit_root: Path) -> dict[str, object]:
    paths = {name: audit_root / name for name in (
        "node_state.csv", "edge_audit.csv", "triangle_audit.csv",
        "summary.json", "manifest.json",
    )}
    missing = [str(path) for path in paths.values() if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"audit root lacks required artifacts: {missing}")
    summary = json.loads(paths["summary.json"].read_text(encoding="utf-8"))
    manifest = json.loads(paths["manifest.json"].read_text(encoding="utf-8"))
    expected_counts = {"node_state": 36, "edge_audit": 54, "triangle_audit": 24}
    if summary.get("schema") != "vela.pn2d_minimal6_fixed_state_audit.v1":
        raise ValueError("audit root has the wrong summary schema")
    if summary.get("status") != "PASS" or summary.get("row_counts") != expected_counts:
        raise ValueError("audit root does not carry a passing exact 36/54/24 summary")
    if manifest.get("schema") != "vela.pn2d_minimal6_fixed_state_audit.v1":
        raise ValueError("audit root has the wrong manifest schema")
    rows = {}
    specifications = {
        "node": ("node_state.csv", 36, ("topology_id", "bias_V", "node_id")),
        "edge": ("edge_audit.csv", 54, ("topology_id", "bias_V", "node0", "node1")),
        "triangle": ("triangle_audit.csv", 24, ("topology_id", "bias_V", "cell_id")),
    }
    expected_states = {(topology, bias) for topology in ("sketch", "mirror")
                       for bias in (0.0, -12.0, -19.0)}
    for kind, (filename, count, key_fields) in specifications.items():
        with paths[filename].open(newline="", encoding="utf-8") as handle:
            values = list(csv.DictReader(handle))
        if len(values) != count:
            raise ValueError(f"audit {kind} rows require exact count {count}")
        keys = [tuple(row[field] for field in key_fields) for row in values]
        if len(set(keys)) != count:
            raise ValueError(f"audit {kind} rows contain duplicate identities")
        identities = {(row["topology_id"], float(row["bias_V"])) for row in values}
        if identities != expected_states:
            raise ValueError(f"audit {kind} rows do not cover the exact six states")
        rows[kind] = values
    rows["artifact_hashes"] = {
        filename: hashlib.sha256(path.read_bytes()).hexdigest()
        for filename, path in paths.items()
    }
    return rows


def _resolved_state(state: dict, state_root: Path) -> dict:
    result = dict(state)
    for key in ("export_dir", "bundle_dir"):
        if key not in result:
            continue
        path = Path(result[key])
        if not path.is_absolute():
            path = state_root / path
        result[key] = str(path.resolve())
    return result

def _state_sources(state: dict, audit_state: dict[str, list[dict]]) -> dict:
    export_dir = Path(state["export_dir"])
    if "member_sha256" in state:
        validate_member_hashes(export_dir, state["member_sha256"])
    _validate_export_units(export_dir)
    mesh = json.loads((export_dir / "mesh.json").read_text(encoding="utf-8"))
    fields = export_dir / "fields"
    native = integrate_native_nodal_per_unit_depth(
        mesh, _read_scalar(fields / "ImpactIonization_region0.csv")
    )
    sentaurus = integrate_native_nodal_per_unit_depth(
        mesh,
        sentaurus_alpha_current_nodal(
            _read_scalar(fields / "eAlphaAvalanche_region0.csv"),
            _read_vector(fields / "eCurrentDensity_region0.csv"),
            _read_scalar(fields / "hAlphaAvalanche_region0.csv"),
            _read_vector(fields / "hCurrentDensity_region0.csv"),
        ),
    )
    vela = sum(
        float(row["vela_total_source_integral_per_m_s"])
        for row in audit_state["triangle"]
    ) * 1.0e-2
    model_path = Path(state.get("bundle_dir", export_dir.parent / "source")) / "models.par"
    models_par = parse_van_overstraeten_de_man(model_path) if model_path.is_file() else None
    vela_parameter_agreement = _vela_parameter_agreement(export_dir, models_par)
    return {
        "topology": state["topology_id"],
        "bias_V": state["requested_bias_V"],
        "export_dir": str(export_dir),
        "models_par_sha256": None if models_par is None else models_par["sha256"],
        "models_par_parameters": models_par,
        "vela_parameter_agreement": vela_parameter_agreement,
        "audit_node_rows": audit_state["node"],
        "audit_edge_rows": audit_state["edge"],
        "audit_triangle_rows": audit_state["triangle"],
        "depth_convention": native["depth_convention"],
        "sentaurus_native_ImpactIonization_s_inv_per_unit_depth": native["value_s_inv_per_unit_depth"],
        "sentaurus_alpha_current_reconstruction_s_inv_per_unit_depth": sentaurus["value_s_inv_per_unit_depth"],
        "vela_reconstruction_s_inv_per_unit_depth": vela,
        "sentaurus_native_minus_reconstruction": source_log_gap(native["value_s_inv_per_unit_depth"], sentaurus["value_s_inv_per_unit_depth"]),
        "vela_native_minus_reconstruction": source_log_gap(native["value_s_inv_per_unit_depth"], vela),
    }

def _node_state_rows(record: dict):
    fields = Path(record["export_dir"]) / "fields"
    raw = {}
    for quantity, unit in _NODE_STATE_FIELDS.items():
        path = fields / f"{quantity}_region0.csv"
        if not path.is_file():
            raise FileNotFoundError(f"missing required raw node field {path}")
        raw[quantity] = _read_scalar(path)
    node_ids = set(raw["ElectrostaticPotential"])
    for quantity, values in raw.items():
        if set(values) != node_ids:
            raise ValueError(f"raw node field {quantity} has inconsistent node IDs")
    for quantity, unit in _NODE_STATE_FIELDS.items():
        for node_id in sorted(node_ids):
            yield _base(record, "node_state", node_id=node_id, quantity=quantity, value=raw[quantity][node_id], unit=unit, source="sentaurus_export", source_kind=SourceKind.SENTAURUS.value)
    for node_id in sorted(node_ids):
        thermal_voltage = 8.617333262145e-5 * raw["LatticeTemperature"][node_id]
        arguments = {
            "psi_V": raw["ElectrostaticPotential"][node_id],
            "phin_V": raw["eQuasiFermiPotential"][node_id],
            "phip_V": raw["hQuasiFermiPotential"][node_id],
            "n_cm3": raw["eDensity"][node_id],
            "p_cm3": raw["hDensity"][node_id],
            "thermal_voltage_V": thermal_voltage,
        }
        electron_log = math.log(arguments["n_cm3"]) + (arguments["psi_V"] - arguments["phin_V"]) / thermal_voltage
        hole_log = math.log(arguments["p_cm3"]) + (arguments["phip_V"] - arguments["psi_V"]) / thermal_voltage
        try:
            ni = infer_ni_eff(**arguments)
        except OverflowError:
            log_max = math.log(sys.float_info.max)
            ni = {
                "electron_cm3": math.exp(electron_log) if electron_log <= log_max else None,
                "hole_cm3": math.exp(hole_log) if hole_log <= log_max else None,
                "relative_residual": -math.expm1(-abs(electron_log - hole_log)),
            }
        for quantity, value in (("ni_eff_electron", ni["electron_cm3"]),
                                ("ni_eff_hole", ni["hole_cm3"])):
            yield _base(record, "node_state", node_id=node_id, quantity=quantity,
                        value="" if value is None else value, unit="cm^-3",
                        source="recomputed_from_sentaurus_export",
                        status="unavailable" if value is None else "available")
        yield _base(record, "node_state", node_id=node_id, quantity="ni_eff_electron_log10_cm3", value=electron_log / math.log(10.0), unit="log10(cm^-3)", source="recomputed_from_sentaurus_export")
        yield _base(record, "node_state", node_id=node_id, quantity="ni_eff_hole_log10_cm3", value=hole_log / math.log(10.0), unit="log10(cm^-3)", source="recomputed_from_sentaurus_export")
        reference_log = math.log(1.45e10)
        yield _base(record, "node_state", node_id=node_id, quantity="bgn_shift_electron_V", value=thermal_voltage * (electron_log - reference_log), unit="V", source="recomputed_from_sentaurus_export")
        yield _base(record, "node_state", node_id=node_id, quantity="bgn_shift_hole_V", value=thermal_voltage * (hole_log - reference_log), unit="V", source="recomputed_from_sentaurus_export")
        yield _base(record, "node_state", node_id=node_id, quantity="ni_eff_electron_over_reference_log10", value=(electron_log - reference_log) / math.log(10.0), unit="dex", source="recomputed_from_sentaurus_export")
        yield _base(record, "node_state", node_id=node_id, quantity="ni_eff_hole_over_reference_log10", value=(hole_log - reference_log) / math.log(10.0), unit="dex", source="recomputed_from_sentaurus_export")
        yield _base(record, "node_state", node_id=node_id, quantity="ni_eff_relative_residual", value=ni["relative_residual"], unit="1", source="recomputed_from_sentaurus_export")

def _audit_unit(quantity: str) -> str:
    suffixes = (
        ("_m2_per_V_s", "m^2*V^-1*s^-1"),
        ("_A_per_m2", "A/m^2"),
        ("_V_per_m", "V/m"),
        ("_per_m2_s", "m^-2*s^-1"),
        ("_per_m_s", "m^-1*s^-1"),
        ("_per_s", "s^-1"),
        ("_m3", "m^-3"),
        ("_per_m", "m^-1"),
        ("_m2", "m^2"),
        ("_cm3", "cm^-3"),
        ("_per_cm", "cm^-1"),
        ("_V", "V"),
    )
    return next((unit for suffix, unit in suffixes if suffix in quantity), "1")


def _audit_source(quantity: str) -> tuple[str, str]:
    if quantity.startswith(("raw_", "sentaurus_")):
        return "fixed_state_audit", SourceKind.SENTAURUS.value
    if quantity.startswith("vela_"):
        return "fixed_state_audit", SourceKind.VELA.value
    return "fixed_state_audit", SourceKind.DERIVED.value


def _numeric_audit_values(row: dict, excluded: set[str]):
    for quantity, text in row.items():
        if quantity in excluded:
            continue
        try:
            value = float(text)
        except (TypeError, ValueError):
            continue
        if not math.isfinite(value):
            raise ValueError(f"audit quantity {quantity} is non-finite")
        yield quantity, value


def _audit_replay_rows(record: dict):
    for row in record["audit_node_rows"]:
        node_id = _integer_csv(row["node_id"], "node_id")
        for quantity, value in _numeric_audit_values(
            row, {"topology_id", "bias_V", "node_id"}
        ):
            source, source_kind = _audit_source(quantity)
            yield _base(
                record, "node_replay", node_id=node_id, quantity=quantity,
                value=value, unit=_audit_unit(quantity), source=source,
                source_kind=source_kind,
            )
    for row in record["audit_edge_rows"]:
        edge_id = _integer_csv(row["edge_id"], "edge_id")
        for carrier in ("electron", "hole"):
            for component, key in (
                ("signed_projection", f"sentaurus_{carrier}_projection_A_per_m2"),
                ("magnitude", f"sentaurus_{carrier}_magnitude_A_per_m2"),
            ):
                yield _base(
                    record, "edge_raw", edge_id=edge_id,
                    quantity=f"{carrier}_current", component=component,
                    value=float(row[key]), unit="A/m^2", source="fixed_state_audit",
                    source_kind=SourceKind.SENTAURUS.value,
                )
        for quantity, value in _numeric_audit_values(
            row, {"topology_id", "bias_V", "edge_id", "node0", "node1"}
        ):
            source, source_kind = _audit_source(quantity)
            yield _base(
                record, "edge_replay", edge_id=edge_id, quantity=quantity,
                value=value, unit=_audit_unit(quantity), source=source,
                source_kind=source_kind,
            )
    for row in record["audit_triangle_rows"]:
        cell_id = _integer_csv(row["cell_id"], "cell_id")
        for quantity, prefix, sign in (
            ("minus_grad_psi", "psi", -1.0),
            ("grad_phin", "phin", 1.0),
            ("grad_phip", "phip", 1.0),
        ):
            x = sign * float(row[f"python_grad_{prefix}_x_V_per_m"])
            y = sign * float(row[f"python_grad_{prefix}_y_V_per_m"])
            for component, value in (
                ("x", x), ("y", y), ("magnitude", math.hypot(x, y)),
                ("direction_rad", math.atan2(y, x)),
            ):
                yield _base(
                    record, "cell_replay", cell_id=cell_id, quantity=quantity,
                    component=component, value=value,
                    unit="rad" if component == "direction_rad" else "V/m",
                    source="fixed_state_audit", source_kind=SourceKind.DERIVED.value,
                )
        for quantity, value in _numeric_audit_values(
            row, {"topology_id", "bias_V", "cell_id", "node0", "node1", "node2"}
        ):
            source, source_kind = _audit_source(quantity)
            node_match = re.search(r"_node(\d+)_source_partition", quantity)
            local_match = re.match(r"vela_local_edge(\d+)_", quantity)
            fields = {"cell_id": cell_id}
            kind = "cell_replay"
            if node_match:
                fields["node_id"] = int(node_match.group(1))
                kind = "node_source_mapping"
            if local_match:
                fields["edge_id"] = _integer_csv(row[f"vela_local_edge{local_match.group(1)}_edge_id"], "local edge_id")
            yield _base(
                record, kind, quantity=quantity, value=value,
                unit=_audit_unit(quantity), source=source,
                source_kind=source_kind, **fields,
            )


def _write_artifacts(out_dir: Path, records: list[dict], waterfall_paths: list[dict]) -> None:
    fields = ["run_id", "record_kind", "topology", "bias_V", "node_id", "cell_id", "edge_id", "quantity", "component", "value", "unit", "source", "source_kind", "value_s_inv_per_unit_depth", "depth_convention", "status"]
    with (out_dir / "quantity_ledger.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for record in records:
            source_families = (
                (family, SOURCE_VALUE_KEYS[family], SOURCE_FAMILIES[family])
                for family in SOURCE_VALUE_KEYS
            )
            for source, key, source_kind in source_families:
                validate_source_anchor_kind(
                    source, source_kind,
                    native=(source == "sentaurus_native_avalanche_generation"),
                )
                writer.writerow(_base(record, "source_integral", quantity="avalanche_generation", unit="s^-1 per 1 cm depth", source=source, source_kind=source_kind, value_s_inv_per_unit_depth=record[key], depth_convention=record["depth_convention"]))
            writer.writerows(_node_state_rows(record))
            writer.writerows(_audit_replay_rows(record))
    with (out_dir / "factor_waterfall.csv").open("w", newline="", encoding="utf-8") as handle:
        fields = ["topology", "bias_V", "path_identity", "order_index", "factor", "contribution_dex", "recomputed", "status"]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for path in waterfall_paths:
            availability = {row["factor"]: row["status"] for row in path["factor_availability"]}
            for path_identity in ("forward", "reverse"):
                definition = path[path_identity]
                contributions = {row["factor"]: row for row in definition["contributions"]}
                for order_index, factor in enumerate(definition["order"]):
                    contribution = contributions.get(factor, {})
                    writer.writerow({
                        "topology": path["topology"], "bias_V": path["bias_V"],
                        "path_identity": path_identity, "order_index": order_index,
                        "factor": factor, "contribution_dex": contribution.get("contribution_dex", ""),
                        "recomputed": "|".join(contribution.get("recomputed", [])),
                        "status": availability[factor],
                    })
            writer.writerow({
                "topology": path["topology"], "bias_V": path["bias_V"],
                "path_identity": "closure", "order_index": len(FACTOR_DEPENDENCIES),
                "factor": "unattributed_residual", "contribution_dex": path["residual_dex"],
                "recomputed": "", "status": path["status"],
            })
    markdown = """# PN2D minimal6 formula difference

minimal6 diagnostic sweep; not a physical BV curve

All source integrals use 1 cm out-of-plane length.

## Source labels

- Native Sentaurus anchor: `ImpactIonization`.
- Sentaurus reconstruction: `(alpha_e*|J_e| + alpha_h*|J_h|)/q`.
- Vela reconstruction: `alpha*flux*partial_volume` from the triangle audit.

## Root-cause implementation map

- Sentaurus parameter entry: `reference_tcad/pn2d_sentaurus2018_minimal6/source/models.par`; export deck: `reference_tcad/pn2d_sentaurus2018_minimal6/source/pn2d_minimal6_state_sdevice.cmd`.
- Independent formula implementation: `scripts/pn2d_minimal6_diagnostics/physics.py` (`infer_ni_eff`, `van_overstraeten_alpha`) and `scripts/pn2d_minimal6_diagnostics/counterfactual.py`.
- C++ control implementation: `src/physics/ImpactIonizationModel.cpp` (`VanOverstraetenImpactIonization`) and `src/simulation/DCSweep.cpp` (avalanche source assembly).

Parameter agreement is a control only; it does not establish a causal factor without a closed counterfactual substitution.

Counterfactual substitutions follow the explicit dependency DAG in forward and reverse order; every nonzero state closes with the named residual within 1e-10 dex.
"""
    (out_dir / "root_cause_summary.md").write_text(markdown, encoding="utf-8")


def _unavailable_counterfactual(record: dict) -> dict:
    gap = record["vela_native_minus_reconstruction"]
    path = {
        "topology": record["topology"], "bias_V": record["bias_V"],
        "dependency_order": list(FACTOR_DEPENDENCIES),
        "factor_availability": [
            {"factor": factor, "status": "unavailable", "reason": "Vela raw state/operator substitution is not exported"}
            for factor in FACTOR_DEPENDENCIES
        ],
        "forward": {"order": list(FACTOR_DEPENDENCIES), "contributions": []},
        "reverse": {"order": list(reversed(FACTOR_DEPENDENCIES)), "contributions": []},
        "interactions": [], "native_gap_dex": gap["dex"], "residual_dex": gap["dex"], "status": "insufficient_data",
    }
    if gap["classification"] == "available":
        assert_counterfactual_closure(native_gap_dex=gap["dex"], contributions_dex=[], residual_dex=gap["dex"])
    return path

def _mapping_scale(row: dict, source: str, carrier: str) -> float:
    local = sum(
        float(row[f"{source}_local_edge{index}_{carrier}_source_integral_per_m_s"])
        for index in range(3)
    )
    partition = sum(
        float(row[f"{source}_{carrier}_node{node}_source_partition_per_m_s"])
        for node in range(1, 7)
    )
    if local == 0.0:
        if partition != 0.0:
            raise ValueError(f"{source} source-to-node mapping is nonconservative")
        return 1.0
    return partition / local


def _formula_operator_inputs(record: dict):
    node_rows = {_integer_csv(row["node_id"], "node_id"): row for row in record["audit_node_rows"]}
    edge_rows = {
        tuple(sorted((_integer_csv(row["node0"], "node0"), _integer_csv(row["node1"], "node1")))): row
        for row in record["audit_edge_rows"]
    }
    baseline = {factor: [] for factor in FACTOR_DEPENDENCIES}
    replacements = {factor: [] for factor in FACTOR_DEPENDENCIES}
    replacements["current_semantics"] = []
    unavailable = {}
    for triangle in sorted(record["audit_triangle_rows"], key=lambda row: _integer_csv(row["cell_id"], "cell_id")):
        gradient_by_carrier = {
            "electron": math.hypot(float(triangle["python_grad_phin_x_V_per_m"]), float(triangle["python_grad_phin_y_V_per_m"])),
            "hole": math.hypot(float(triangle["python_grad_phip_x_V_per_m"]), float(triangle["python_grad_phip_y_V_per_m"])),
        }
        mapping = {
            (source, carrier): _mapping_scale(triangle, source, carrier)
            for source in ("vela", "python") for carrier in ("electron", "hole")
        }
        for local_index in range(3):
            prefix = f"vela_local_edge{local_index}_"
            first = _integer_csv(triangle[f"node{local_index}"], f"node{local_index}")
            second_index = (local_index + 1) % 3
            second = _integer_csv(triangle[f"node{second_index}"], f"node{second_index}")
            edge = edge_rows[tuple(sorted((first, second)))]
            for carrier, density_name, qf_name in (
                ("electron", "n", "phin"), ("hole", "p", "phip")
            ):
                density = float(triangle[prefix + f"{carrier}_midpoint_density_m3"])
                mobility = float(triangle[prefix + f"{carrier}_mobility_m2_per_V_s"])
                flux = abs(float(triangle[prefix + f"{carrier}_flux_proxy_per_m2_s"]))
                edge_field = flux / (mobility * density) if mobility > 0.0 and density > 0.0 else abs(float(edge[f"delta_{qf_name}_over_h_V_per_m"]))
                baseline["ni_eff/BGN"].append(density)
                replacements["ni_eff/BGN"].append(0.5 * (
                    float(node_rows[first][f"sentaurus_{density_name}_m3"])
                    + float(node_rows[second][f"sentaurus_{density_name}_m3"])
                ))
                baseline["gradient_recovery"].append(edge_field)
                replacements["gradient_recovery"].append(gradient_by_carrier[carrier])
                baseline["mobility"].append(mobility)
                replacements["mobility"].append(0.5 * (
                    float(node_rows[first][f"sentaurus_{carrier}_mobility_m2_per_V_s"])
                    + float(node_rows[second][f"sentaurus_{carrier}_mobility_m2_per_V_s"])
                ))
                baseline["current_semantics"].append(None)
                replacements["current_semantics"].append(
                    abs(float(edge[f"sentaurus_{carrier}_magnitude_A_per_m2"])) / 1.602176634e-19
                )
                baseline["impact_driving_field"].append(abs(float(triangle[prefix + f"{carrier}_cell_qf_field_V_per_m"])))
                electric_x = 0.5 * (
                    float(node_rows[first]["sentaurus_electric_field_x_V_per_m"])
                    + float(node_rows[second]["sentaurus_electric_field_x_V_per_m"])
                )
                electric_y = 0.5 * (
                    float(node_rows[first]["sentaurus_electric_field_y_V_per_m"])
                    + float(node_rows[second]["sentaurus_electric_field_y_V_per_m"])
                )
                replacements["impact_driving_field"].append(math.hypot(electric_x, electric_y))
                baseline["alpha_law"].append(float(triangle[prefix + f"{carrier}_alpha_per_m"]))
                replacements["alpha_law"].append(0.5 * (
                    float(node_rows[first][f"sentaurus_{carrier}_alpha_per_m"])
                    + float(node_rows[second][f"sentaurus_{carrier}_alpha_per_m"])
                ))
                baseline["partial_volume"].append(float(triangle[prefix + "truncated_partial_volume_m2"]))
                replacements["partial_volume"].append(float(triangle[f"python_local_edge{local_index}_truncated_partial_volume_m2"]))
                baseline["source_to_node_mapping"].append(1.0e-2 * mapping[("vela", carrier)])
                replacements["source_to_node_mapping"].append(1.0e-2 * mapping[("python", carrier)])
    baseline["current_semantics"] = None
    baseline = {factor: (value if value is None else tuple(value)) for factor, value in baseline.items()}
    replacements = {factor: tuple(value) for factor, value in replacements.items()}
    unavailable_inputs = {
        "ni_eff/BGN": "raw carrier-density averaging is not an independently inferred ni_eff/BGN replay",
        "impact_driving_field": "direct exported alpha lacks coefficient provenance for independent impact-field replay",
        "alpha_law": "direct exported alpha confounds driving field and law because coefficient provenance is absent",
        "source_to_node_mapping": (
            "Sentaurus mapping weights are not independently exported; conservative Vela/Python weights are controls only"
        ),
    }
    for factor, reason in unavailable_inputs.items():
        replacements.pop(factor)
        unavailable[factor] = reason
    return baseline, replacements, unavailable


def _closed_counterfactual(record: dict) -> dict:
    if record["vela_native_minus_reconstruction"]["classification"] != "available":
        return _unavailable_counterfactual(record)
    native = abs(float(record["sentaurus_native_ImpactIonization_s_inv_per_unit_depth"]))
    baseline, replacements, unavailable = _formula_operator_inputs(record)
    evaluated = evaluate_formula_counterfactual(
        native=native,
        baseline_values=baseline,
        replacement_values=replacements,
        unavailable_reasons=unavailable,
    )
    engine = evaluated.pop("engine")
    availability = {row["factor"]: row["status"] for row in evaluated["factor_availability"]}
    interactions = [
        interaction for interaction in build_adjacent_interactions(
            evaluated["forward"]["contributions"],
            evaluated["reverse"]["contributions"],
            engine.evaluate_replacements,
        )
        if availability[interaction["first_factor"]] == "available"
        and availability[interaction["second_factor"]] == "available"
    ]
    evaluated.update({
        "topology": record["topology"],
        "bias_V": record["bias_V"],
        "interactions": interactions,
        "symmetric_contributions": symmetric_contributions(
            evaluated["forward"]["contributions"],
            evaluated["reverse"]["contributions"],
        ),
        "status": "closed" if not unavailable else "insufficient_data",
    })
    return evaluated

def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--state-root", type=Path, required=True)
    parser.add_argument("--audit-root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--qa-reviewer", default="unreviewed")
    parser.add_argument("--qa-date", default="")
    parser.add_argument("--qa-status", choices=("pending_visual_inspection", "reviewed"), default="pending_visual_inspection")
    args = parser.parse_args()
    validate_dependency_dag(FACTOR_DEPENDENCIES)
    _validate_source_families()
    manifest_path = args.state_root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    base = validate_formula_input(manifest)
    if manifest.get("schema") == "vela.pn2d_minimal6_states.v2":
        validate_sealed_archive(args.state_root)
    audit = validate_audit_binding(args.state_root, args.audit_root)
    resolved_states = [_resolved_state(state, args.state_root) for state in manifest["states"]]
    records = []
    for state in resolved_states:
        topology = state["topology_id"]
        bias = float(state["requested_bias_V"])
        audit_state = {
            kind: [
                row for row in audit[kind]
                if row["topology_id"] == topology and float(row["bias_V"]) == bias
            ]
            for kind in ("node", "edge", "triangle")
        }
        if tuple(len(audit_state[kind]) for kind in ("node", "edge", "triangle")) != (6, 9, 4):
            raise ValueError(f"audit rows do not bind exact state {topology} {bias:g} V")
        records.append(_state_sources(state, audit_state))
    waterfall_paths = [_closed_counterfactual(record) for record in records]
    residuals = [
        {
            "name": "sentaurus_internal_semantics_residual",
            "topology": record["topology"],
            "bias_V": record["bias_V"],
            "classification": (
                "available" if path["residual_dex"] is not None
                else record["sentaurus_native_minus_reconstruction"]["classification"]
            ),
            "dex": path["residual_dex"],
        }
        for record, path in zip(records, waterfall_paths)
    ]
    interactions = [
        {"topology": path["topology"], "bias_V": path["bias_V"], **interaction}
        for path in waterfall_paths
        for interaction in path["interactions"]
    ]
    dominance = score_dominance([{
        "topology": path["topology"],
        "bias_V": path["bias_V"],
        "native_gap_dex": path["native_gap_dex"] or 0.0,
        "residual_dex": path["residual_dex"] or 0.0,
        "symmetric_contributions": path.get("symmetric_contributions", {}),
        "factor_availability": path["factor_availability"],
    } for path in waterfall_paths])
    state_matrix = [{
        "topology_id": state["topology_id"],
        "requested_bias_V": state["requested_bias_V"],
        "actual_bias_V": state["actual_bias_V"],
        "status": state["status"],
    } for state in manifest["states"]]
    artifact_hashes = {
        "state_manifest_sha256": hashlib.sha256(manifest_path.read_bytes()).hexdigest(),
        "audit_artifact_sha256": audit["artifact_hashes"],
        "models_par_sha256": sorted({
            item["models_par_sha256"] for item in records
            if item["models_par_sha256"] is not None
        }),
    }
    report = {
        "schema": "vela.pn2d_minimal6_formula_difference.v1",
        "diagnostic_disclaimer": DISCLAIMER,
        "input_provenance": {"state_manifest": str(manifest_path.resolve())},
        "audit_provenance": {
            "audit_root": str(args.audit_root.resolve()),
            "artifact_sha256": audit["artifact_hashes"],
        },
        "state_matrix": state_matrix,
        "row_counts": base["row_counts"],
        "waterfall_paths": waterfall_paths,
        "interactions": interactions,
        "dominance_rules": dominance,
        "sentaurus_internal_semantics_residual": residuals,
        "vela_parameter_agreement": [
            {"topology": item["topology"], "bias_V": item["bias_V"],
             **item["vela_parameter_agreement"]} for item in records
        ],
        "artifact_hashes": artifact_hashes,
        "records": records,
        "root_cause_status": dominance["status"],
        "root_cause_reason": dominance.get(
            "reason", "closed fixed-state counterfactual decomposition"
        ),
    }
    validate_formula_difference_v1(report)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "root_cause_summary.json").write_text(json.dumps(report, sort_keys=True, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    _write_artifacts(args.out_dir, records, waterfall_paths)
    figure_manifest = render_formula_difference_figures(
        ledger_path=args.out_dir / "quantity_ledger.csv",
        waterfall_path=args.out_dir / "factor_waterfall.csv",
        report_path=args.out_dir / "root_cause_summary.json",
        out_dir=args.out_dir,
        reviewer=args.qa_reviewer,
        reviewed_on=args.qa_date,
        qa_status=args.qa_status,
    )
    (args.out_dir / "figure_manifest.json").write_text(
        json.dumps(figure_manifest, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
