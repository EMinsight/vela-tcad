#!/usr/bin/env python3
"""Create fixed-state PN2D minimal6 formula-difference diagnostics."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.export_pn2d_minimal6_states import validate_member_hashes
from scripts.pn2d_minimal6_diagnostics.physics import (
    compare_van_overstraeten_parameters,
    infer_ni_eff,
    parse_van_overstraeten_de_man,
    van_overstraeten_alpha,
)
from scripts.pn2d_minimal6_diagnostics.counterfactual import (
    FACTOR_DEPENDENCIES,
    assert_counterfactual_closure,
    score_dominance,
    validate_field_units,
    validate_source_anchor_kind,
    integrate_native_nodal_per_unit_depth,
    integrate_vela_reconstructed_per_unit_depth,
    sentaurus_alpha_current_nodal,
    source_log_gap,
    validate_formula_input,
)
from scripts.pn2d_minimal6_diagnostics.schemas import DISCLAIMER, validate_formula_difference_v1
from scripts.pn2d_minimal6_diagnostics.support import project_vector_to_edge
from scripts.pn2d_minimal6_diagnostics.plots import render_formula_difference_figures

_NODE_STATE_FIELDS = {
    "ElectrostaticPotential": "V", "eDensity": "cm^-3", "hDensity": "cm^-3",
    "eQuasiFermiPotential": "V", "hQuasiFermiPotential": "V",
    "eMobility": "cm^2*V^-1*s^-1", "hMobility": "cm^2*V^-1*s^-1",
    "eVelocity": "cm*s^-1", "hVelocity": "cm*s^-1",
    "LatticeTemperature": "K",
    "eAlphaAvalanche": "cm^-1", "hAlphaAvalanche": "cm^-1",
}


def _read_scalar(path: Path) -> dict[int, float]:
    with path.open(newline="", encoding="utf-8") as handle:
        return {int(row["node_id"]): float(row["component0"]) for row in csv.DictReader(handle)}


def _read_vector(path: Path) -> dict[int, tuple[float, float]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return {int(row["node_id"]): (float(row["component0"]), float(row["component1"])) for row in csv.DictReader(handle)}


def _base(record: dict, kind: str, **fields) -> dict:
    row = {
        "run_id": "minimal6_fixed_state", "record_kind": kind,
        "topology": record["topology"], "bias_V": record["bias_V"],
        "node_id": "", "cell_id": "", "edge_id": "", "quantity": "", "component": "",
        "value": "", "unit": "", "source": "", "value_s_inv_per_unit_depth": "",
        "depth_convention": "", "status": "available",
    }
    row.update(fields)
    return row


def _vector_rows(record: dict, kind: str, *, entity: str, entity_id: int, quantity: str, x: float, y: float, unit: str):
    common = {entity: entity_id, "quantity": quantity, "unit": unit}
    return [
        _base(record, kind, **common, component="x", value=x),
        _base(record, kind, **common, component="y", value=y),
        _base(record, kind, **common, component="magnitude", value=math.hypot(x, y)),
        _base(record, kind, **common, component="direction_rad", value=math.atan2(y, x)),
    ]


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

def _state_sources(state: dict) -> dict:
    export_dir = Path(state["export_dir"])
    if "member_sha256" in state:
        validate_member_hashes(export_dir, state["member_sha256"])
    _validate_export_units(export_dir)
    mesh = json.loads((export_dir / "mesh.json").read_text(encoding="utf-8"))
    fields = export_dir / "fields"
    native = integrate_native_nodal_per_unit_depth(mesh, _read_scalar(fields / "ImpactIonization_region0.csv"))
    sentaurus = integrate_native_nodal_per_unit_depth(
        mesh,
        sentaurus_alpha_current_nodal(
            _read_scalar(fields / "eAlphaAvalanche_region0.csv"),
            _read_vector(fields / "eCurrentDensity_region0.csv"),
            _read_scalar(fields / "hAlphaAvalanche_region0.csv"),
            _read_vector(fields / "hCurrentDensity_region0.csv"),
        ),
    )
    with (export_dir / "vela_triangle_audit.csv").open(newline="", encoding="utf-8") as handle:
        vela = integrate_vela_reconstructed_per_unit_depth(csv.DictReader(handle))
    model_path = Path(state.get("bundle_dir", export_dir.parent / "source")) / "models.par"
    models_par = parse_van_overstraeten_de_man(model_path) if model_path.is_file() else None
    vela_parameter_agreement = _vela_parameter_agreement(export_dir, models_par)
    return {
        "topology": state["topology_id"], "bias_V": state["requested_bias_V"], "export_dir": str(export_dir),
        "models_par_sha256": None if models_par is None else models_par["sha256"],
        "models_par_parameters": models_par,
        "vela_parameter_agreement": vela_parameter_agreement,
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
            yield _base(record, "node_state", node_id=node_id, quantity=quantity, value=raw[quantity][node_id], unit=unit, source="sentaurus_export")
    for node_id in sorted(node_ids):
        thermal_voltage = 8.617333262145e-5 * raw["LatticeTemperature"][node_id]
        ni = infer_ni_eff(
            psi_V=raw["ElectrostaticPotential"][node_id],
            phin_V=raw["eQuasiFermiPotential"][node_id],
            phip_V=raw["hQuasiFermiPotential"][node_id],
            n_cm3=raw["eDensity"][node_id],
            p_cm3=raw["hDensity"][node_id],
            thermal_voltage_V=thermal_voltage,
        )
        yield _base(record, "node_state", node_id=node_id, quantity="ni_eff_electron", value=ni["electron_cm3"], unit="cm^-3", source="recomputed_from_sentaurus_export")
        yield _base(record, "node_state", node_id=node_id, quantity="ni_eff_hole", value=ni["hole_cm3"], unit="cm^-3", source="recomputed_from_sentaurus_export")
        yield _base(record, "node_state", node_id=node_id, quantity="ni_eff_relative_residual", value=ni["relative_residual"], unit="1", source="recomputed_from_sentaurus_export")

def _audit_replay_rows(record: dict):
    audit_path = Path(record["export_dir"]) / "vela_triangle_audit.csv"
    with audit_path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            try:
                cell_id = int(row["cell_id"])
                gradients = (("minus_grad_psi", "grad_psi_x_V_per_m", "grad_psi_y_V_per_m"),
                             ("grad_phin", "grad_phin_x_V_per_m", "grad_phin_y_V_per_m"),
                             ("grad_phip", "grad_phip_x_V_per_m", "grad_phip_y_V_per_m"))
                for quantity, x_key, y_key in gradients:
                    x, y = float(row[x_key]), float(row[y_key])
                    if quantity == "minus_grad_psi":
                        x, y = -x, -y
                    yield from _vector_rows(record, "cell_replay", entity="cell_id", entity_id=cell_id, quantity=quantity, x=x, y=y, unit="V/m")
            except KeyError as exc:
                raise ValueError(f"audit row lacks required cell field {exc.args[0]}") from exc
            local_indices = sorted({int(match.group(1)) for key in row if (match := re.fullmatch(r"local_edge(\d+)_edge_id", key))})
            for index in local_indices:
                prefix = f"local_edge{index}_"
                try:
                    edge_id = int(row[prefix + "edge_id"])
                except (KeyError, ValueError) as exc:
                    raise ValueError(f"audit row has invalid {prefix}edge_id") from exc
                quantities = (
                    ("electron_qf_field", "electron_cell_qf_field_V_per_m", "V/m"),
                    ("hole_qf_field", "hole_cell_qf_field_V_per_m", "V/m"),
                    ("electron_midpoint_density", "electron_midpoint_density_m3", "m^-3"),
                    ("hole_midpoint_density", "hole_midpoint_density_m3", "m^-3"),
                    ("electron_mobility", "electron_mobility_m2_per_V_s", "m^2*V^-1*s^-1"),
                    ("hole_mobility", "hole_mobility_m2_per_V_s", "m^2*V^-1*s^-1"),
                    ("electron_alpha", "electron_alpha_per_m", "m^-1"),
                    ("hole_alpha", "hole_alpha_per_m", "m^-1"),
                    ("electron_sg_flux", "electron_flux_proxy_per_m2_s", "m^-2*s^-1"),
                    ("hole_sg_flux", "hole_flux_proxy_per_m2_s", "m^-2*s^-1"),
                    ("electron_source_integral", "electron_source_integral_per_m_s", "m^-1*s^-1"),
                    ("hole_source_integral", "hole_source_integral_per_m_s", "m^-1*s^-1"),
                )
                for quantity, suffix, unit in quantities:
                    key = prefix + suffix
                    if key not in row:
                        raise ValueError(f"audit row lacks required edge field {key}")
                    yield _base(record, "edge_replay", cell_id=cell_id, edge_id=edge_id, quantity=quantity, value=float(row[key]), unit=unit, source="vela_triangle_audit")
                if record["models_par_parameters"] is not None:
                    for carrier in ("electron", "hole"):
                        params = record["models_par_parameters"][carrier]
                        field_v_per_cm = abs(float(row[prefix + f"{carrier}_cell_qf_field_V_per_m"])) / 100.0
                        alpha = van_overstraeten_alpha(
                            field_v_per_cm,
                            params["a_low_cm_inv"], params["b_low_v_per_cm"],
                            params["a_high_cm_inv"], params["b_high_v_per_cm"],
                            record["models_par_parameters"]["switch_field_v_per_cm"],
                        )
                        yield _base(record, "edge_replay", cell_id=cell_id, edge_id=edge_id,
                                    quantity=f"sentaurus_dem_{carrier}_alpha_recomputed", value=alpha,
                                    unit="cm^-1", source="models_par_reference")


def _sentaurus_edge_current_rows(record: dict):
    export_dir = Path(record["export_dir"])
    mesh = json.loads((export_dir / "mesh.json").read_text(encoding="utf-8"))
    coordinates = {int(node["id"]): (float(node["x"]), float(node["y"])) for node in mesh["nodes"]}
    fields = export_dir / "fields"
    carriers = {"electron": _read_vector(fields / "eCurrentDensity_region0.csv"), "hole": _read_vector(fields / "hCurrentDensity_region0.csv")}
    edges = set()
    for triangle in mesh["triangles"]:
        ids = [int(node) for node in triangle["node_ids"]]
        for first, second in ((ids[0], ids[1]), (ids[1], ids[2]), (ids[2], ids[0])):
            edges.add((min(first, second), max(first, second)))
    for edge_id, (first, second) in enumerate(sorted(edges)):
        for carrier, vectors in carriers.items():
            vector = (0.5 * (vectors[first][0] + vectors[second][0]), 0.5 * (vectors[first][1] + vectors[second][1]))
            yield _base(record, "edge_raw", edge_id=edge_id, quantity=f"{carrier}_current", component="signed_projection", value=project_vector_to_edge(vector, coordinates[first], coordinates[second]), unit="A/cm^2", source="sentaurus_export")
            yield _base(record, "edge_raw", edge_id=edge_id, quantity=f"{carrier}_current", component="magnitude", value=math.hypot(*vector), unit="A/cm^2", source="sentaurus_export")


def _write_artifacts(out_dir: Path, records: list[dict], waterfall_paths: list[dict]) -> None:
    fields = ["run_id", "record_kind", "topology", "bias_V", "node_id", "cell_id", "edge_id", "quantity", "component", "value", "unit", "source", "value_s_inv_per_unit_depth", "depth_convention", "status"]
    with (out_dir / "quantity_ledger.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for record in records:
            for source, key in (("sentaurus_native_avalanche_generation", "sentaurus_native_ImpactIonization_s_inv_per_unit_depth"), ("sentaurus_alpha_current_reconstruction", "sentaurus_alpha_current_reconstruction_s_inv_per_unit_depth"), ("vela_alpha_flux_partial_volume_reconstruction", "vela_reconstruction_s_inv_per_unit_depth")):
                validate_source_anchor_kind(source, native=(source == "sentaurus_native_avalanche_generation"))
                writer.writerow(_base(record, "source_integral", quantity="avalanche_generation", unit="s^-1 per 1 cm depth", source=source, value_s_inv_per_unit_depth=record[key], depth_convention=record["depth_convention"]))
            writer.writerows(_node_state_rows(record))
            writer.writerows(_audit_replay_rows(record))
            writer.writerows(_sentaurus_edge_current_rows(record))
    with (out_dir / "factor_waterfall.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["topology", "bias_V", "factor", "contribution_dex", "status"])
        writer.writeheader()
        for path in waterfall_paths:
            for factor in path["factor_availability"]:
                writer.writerow({"topology": path["topology"], "bias_V": path["bias_V"],
                                 "factor": factor["factor"], "contribution_dex": "", "status": factor["status"]})
            writer.writerow({"topology": path["topology"], "bias_V": path["bias_V"],
                             "factor": "unattributed_residual", "contribution_dex": path["residual_dex"], "status": path["status"]})
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

Counterfactual factor substitutions are unavailable until Vela raw nodal states and operator inputs are exported; no causal ranking is emitted.
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

def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--state-root", type=Path, required=True)
    parser.add_argument("--audit-root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--qa-reviewer", default="unreviewed")
    parser.add_argument("--qa-date", default="")
    parser.add_argument("--qa-status", choices=("pending_visual_inspection", "reviewed"), default="pending_visual_inspection")
    args = parser.parse_args()
    manifest_path = args.state_root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    base = validate_formula_input(manifest)
    records = [_state_sources(state) for state in manifest["states"]]
    residuals = [{"topology": item["topology"], "bias_V": item["bias_V"], **item["sentaurus_native_minus_reconstruction"]} for item in records]
    waterfall_paths = [_unavailable_counterfactual(record) for record in records]
    dominance = score_dominance([{
        "topology": path["topology"], "bias_V": path["bias_V"],
        "native_gap_dex": path["native_gap_dex"] or 0.0, "residual_dex": path["residual_dex"] or 0.0,
        "symmetric_contributions": {},
    } for path in waterfall_paths])
    report = {
        "schema": "vela.pn2d_minimal6_formula_difference.v1", "diagnostic_disclaimer": DISCLAIMER,
        "input_provenance": {"state_manifest": str(manifest_path)}, "audit_provenance": {"audit_root": str(args.audit_root)},
        "state_matrix": manifest["states"], "row_counts": base["row_counts"], "waterfall_paths": waterfall_paths, "interactions": [], "dominance_rules": dominance,
        "sentaurus_internal_semantics_residual": residuals,
        "vela_parameter_agreement": [{"topology": item["topology"], "bias_V": item["bias_V"], **item["vela_parameter_agreement"]} for item in records],
        "artifact_hashes": {"state_manifest_sha256": hashlib.sha256(manifest_path.read_bytes()).hexdigest(), "models_par_sha256": sorted({item["models_par_sha256"] for item in records if item["models_par_sha256"] is not None})},
        "records": records, "root_cause_status": "insufficient_data",
        "root_cause_reason": "raw Sentaurus and Vela-replay ledgers are complete, but Vela raw state needed for named counterfactual substitutions is unavailable",
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
