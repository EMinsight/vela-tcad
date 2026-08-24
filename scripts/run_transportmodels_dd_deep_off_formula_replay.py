#!/usr/bin/env python3
"""Replay three Sentaurus DD deep-off states through Vela production formulas."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[1]
REF = REPO / "build-release/reference_tcad/transportmodels_sentaurus2022"
MANIFEST = (
    REF
    / "sentaurus_vm_runs/dd_deep_off_spatial_oracles_20260824"
    / "dd_deep_off_spatial_oracles_manifest.json"
)
DD_CONFIG = (
    REF
    / "vela_baseline/dd_dg_fixed_contract_v1_2026-08-24/runs/dd"
    / "03_dd_idvg_curve.json"
)
SENTAURUS_CURVE = REF / "run02/normalized/dd_idvg.csv"
SENTAURUS_PLT = REF / "run02/full_raw/IdVgs_n6_des.plt"
VELA_CURVE = (
    REF
    / "vela_baseline/dd_dg_fixed_contract_v1_2026-08-24"
    / "dd_idvg_completed_prefix.csv"
)
VELA_SRH = (
    REF
    / "vela_baseline/dd_dg_srh_corrected_cold_regression_2026-08-23/runs/dd"
    / "dd_idvg_curve_srh_balance.csv"
)
VELA_RUN = (
    REF
    / "vela_baseline/dd_dg_fixed_contract_v1_2026-08-24/runs/dd"
)
VELA_STATE_FILES = {
    -1.00: VELA_RUN / "dd_idvg_final_bias_relax_state_bias_m1p000000.csv",
    -0.84: VELA_RUN / "dd_idvg_curve_state_bias_m0p840000.csv",
    -0.68: VELA_RUN / "dd_idvg_curve_state_bias_m0p680000.csv",
}
LEGACY_FIXED_STATE = (
    REF / "reports/idvg_dd_deep_off_fixed_state_20260823/summary.json"
)
OUTPUT = REF / "reports/transportmodels_dd_deep_off_formula_replay_20260824"
SUMMARY = OUTPUT / "summary.json"
TABLE = OUTPUT / "bias_summary.csv"
DEFAULT_RUNNER = REPO / "build-release/vela_example_runner.exe"
Q = 1.602176634e-19


def load_base() -> Any:
    path = REPO / "scripts/run_transportmodels_sentaurus_formula_replay.py"
    spec = importlib.util.spec_from_file_location("transportmodels_formula_replay", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.IDVG_CONFIG = DD_CONFIG
    return module


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def at_bias(rows: list[dict[str, str]], bias: float) -> dict[str, str] | None:
    matches = [row for row in rows if abs(float(row["bias_V"]) - bias) < 1.0e-10]
    return matches[-1] if matches else None


def sentaurus_terminal_rows() -> list[dict[str, float]]:
    path = REPO / "scripts/sentaurus_import.py"
    spec = importlib.util.spec_from_file_location("sentaurus_import", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    text = SENTAURUS_PLT.read_text(errors="ignore")
    datasets = module.parse_quoted_list(text, "datasets")
    values = module.parse_values_block(text, len(datasets))
    result: list[dict[str, float]] = []
    for values_row in values[:3]:
        row = dict(zip(datasets, values_row, strict=True))
        result.append(
            {
                "bias_V": float(row["gate OuterVoltage"]),
                "gate_A": float(row["gate TotalCurrent"]),
                "substrate_A": float(row["substrate TotalCurrent"]),
                "drain_A": float(row["drain TotalCurrent"]),
                "source_A": float(row["source TotalCurrent"]),
            }
        )
    return result


def sentaurus_srh_integral(base: Any, export_dir: Path) -> dict[str, float]:
    nodes = {
        int(row["id"]): (float(row["x_um"]), float(row["y_um"]))
        for row in read_csv(export_dir / "nodes.csv")
    }
    srh = base.field(export_dir, "srhRecombination", 3)
    signed = 0.0
    absolute = 0.0
    for row in read_csv(export_dir / "elements.csv"):
        if row["region"] != "R.Substrate":
            continue
        ids = tuple(int(row[f"node{i}"]) for i in range(3))
        (x0, y0), (x1, y1), (x2, y2) = (nodes[node] for node in ids)
        area_um2 = 0.5 * abs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0))
        rate = sum(srh[node] for node in ids) / 3.0
        signed += rate * area_um2
        absolute += abs(rate) * area_um2
    # cm^-3 s^-1 x um^2 x 1 um width -> 1e-12 carriers/s.
    return {
        "signed_A_per_um": Q * signed * 1.0e-12,
        "absolute_A_per_um": Q * absolute * 1.0e-12,
    }


def contact_scalars(export_dir: Path, field_name: str) -> dict[str, float]:
    manifest = json.loads((export_dir / "field_manifest.json").read_text(encoding="utf-8"))
    result: dict[str, float] = {}
    for item in manifest["fields"]:
        if item["name"] != field_name:
            continue
        rows = read_csv(export_dir / "fields" / item["csv_file"])
        if rows:
            result[item["region_name"]] = float(rows[0]["component0"])
    return result


def contact_nodes(export_dir: Path) -> dict[str, set[int]]:
    result: dict[str, set[int]] = {}
    for row in read_csv(export_dir / "contacts.csv"):
        result[row["name"]] = {
            int(value) for value in row["node_ids"].split(";") if value
        }
    return result


def vela_contact_cut_currents(
    export_dir: Path, edge_rows: list[dict[str, str]]
) -> dict[str, dict[str, float]]:
    result: dict[str, dict[str, float]] = {}
    for name, nodes in contact_nodes(export_dir).items():
        electron_particles = 0.0
        hole_particles = 0.0
        crossing_edges = 0
        for row in edge_rows:
            at0 = int(row["node0"]) in nodes
            at1 = int(row["node1"]) in nodes
            if at0 == at1:
                continue
            # Positive orientation is from the contact basin into the device.
            sign = 1.0 if at0 else -1.0
            electron_particles += sign * float(row["electron_particle_line_flux_per_m_s"])
            hole_particles += sign * float(row["hole_particle_line_flux_per_m_s"])
            crossing_edges += 1
        electron_A_um = -Q * electron_particles * 1.0e-6
        hole_A_um = Q * hole_particles * 1.0e-6
        result[name] = {
            "crossing_edge_count": crossing_edges,
            "electron_A_per_um": electron_A_um,
            "hole_A_per_um": hole_A_um,
            "total_A_per_um": electron_A_um + hole_A_um,
        }
    return result


def contact_cut_precision_metrics(
    export_dir: Path, edge_rows: list[dict[str, str]], contact: str
) -> dict[str, float]:
    nodes = contact_nodes(export_dir)[contact]
    contributions: list[float] = []
    qf_drops: list[float] = []
    qf_drop_ulps: list[float] = []
    for row in edge_rows:
        at0 = int(row["node0"]) in nodes
        at1 = int(row["node1"]) in nodes
        if at0 == at1:
            continue
        sign = 1.0 if at0 else -1.0
        contribution = sign * (
            -Q * float(row["electron_particle_line_flux_per_m_s"])
            + Q * float(row["hole_particle_line_flux_per_m_s"])
        ) * 1.0e-6
        phin0 = float(row["phin0_V"])
        phin1 = float(row["phin1_V"])
        drop = abs(phin1 - phin0)
        ulp = max(math.ulp(phin0), math.ulp(phin1))
        contributions.append(contribution)
        qf_drops.append(drop)
        qf_drop_ulps.append(drop / ulp if ulp > 0.0 else math.nan)
    net = sum(contributions)
    ordered_abs = sorted((abs(value) for value in contributions), reverse=True)
    ordered_drop = sorted(qf_drops)
    ordered_ulp = sorted(value for value in qf_drop_ulps if math.isfinite(value))
    median = lambda values: values[len(values) // 2] if values else math.nan
    return {
        "crossing_edge_count": len(contributions),
        "nonzero_qf_drop_edge_count": sum(drop > 0.0 for drop in qf_drops),
        "median_abs_qf_drop_V": median(ordered_drop),
        "maximum_abs_qf_drop_V": max(qf_drops, default=math.nan),
        "median_abs_qf_drop_ulp": median(ordered_ulp),
        "maximum_abs_qf_drop_ulp": max(ordered_ulp, default=math.nan),
        "net_current_A_per_um": net,
        "sum_abs_edge_current_A_per_um": sum(ordered_abs),
        "edge_cancellation_condition": sum(ordered_abs) / max(abs(net), 1.0e-300),
        "largest_edge_to_net_ratio": max(ordered_abs, default=0.0) / max(abs(net), 1.0e-300),
    }


def qf_drive_metrics(
    base: Any,
    export_dir: Path,
    mobility_rows: list[dict[str, str]],
    active_edge_ids: set[int],
) -> dict[str, Any]:
    sent_gradient = base.field(export_dir, "eGradQuasiFermi", 3, 2)
    all_errors: list[float] = []
    active_errors: list[float] = []
    pairs: list[tuple[int, float, float]] = []
    for row in mobility_rows:
        n0, n1 = int(row["node0"]), int(row["node1"])
        if n0 not in sent_gradient or n1 not in sent_gradient:
            continue
        sent = 0.5 * (
            math.hypot(*sent_gradient[n0]) + math.hypot(*sent_gradient[n1])
        )
        vela = abs(float(row["electron_mobility_field_V_m"])) / 100.0
        pairs.append((int(row["edge_id"]), sent, vela))
    floor = max((max(sent, vela) for _, sent, vela in pairs), default=0.0) * 1.0e-10
    for edge, sent, vela in pairs:
        if sent <= floor or vela <= floor:
            continue
        error = abs(math.log10(vela) - math.log10(sent))
        all_errors.append(error)
        if edge in active_edge_ids:
            active_errors.append(error)
    return {
        "sentaurus_field": "eGradQuasiFermi/Vector",
        "vela_field": "electron_mobility_field_V_m",
        "all_edges_abs_error_dex": base.error_stats(all_errors),
        "current_carrying_edges_abs_error_dex": base.error_stats(active_errors),
    }


def substitute_sentaurus_nodal_mobility(
    base: Any, export_dir: Path, edge_rows: list[dict[str, str]]
) -> list[dict[str, str]]:
    sent_n = base.field(export_dir, "eMobility", 3)
    sent_p = base.field(export_dir, "hMobility", 3)
    result: list[dict[str, str]] = []
    for source in edge_rows:
        row = dict(source)
        n0, n1 = int(row["node0"]), int(row["node1"])
        if n0 in sent_n and n1 in sent_n:
            sent_mu = 0.5 * (sent_n[n0] + sent_n[n1]) * 1.0e-4
            vela_mu = float(row["electron_mobility_m2_V_s"])
            if vela_mu > 0.0:
                row["electron_particle_line_flux_per_m_s"] = str(
                    float(row["electron_particle_line_flux_per_m_s"])
                    * sent_mu
                    / vela_mu
                )
        if n0 in sent_p and n1 in sent_p:
            sent_mu = 0.5 * (sent_p[n0] + sent_p[n1]) * 1.0e-4
            vela_mu = float(row["hole_mobility_m2_V_s"])
            if vela_mu > 0.0:
                row["hole_particle_line_flux_per_m_s"] = str(
                    float(row["hole_particle_line_flux_per_m_s"])
                    * sent_mu
                    / vela_mu
                )
        result.append(row)
    return result


def recombination_by_variant(
    rows: list[dict[str, str]], silicon_nodes: set[int], scale: float
) -> dict[str, float]:
    variants = sorted({row["variant"] for row in rows})
    return {
        variant: scale * sum(
            float(row["electron_recombination"])
            for row in rows
            if row["variant"] == variant and int(row["node_id"]) in silicon_nodes
        )
        for variant in variants
    }


def execute_case(
    base: Any, runner: Path, state: dict[str, Any], source_scale: float
) -> dict[str, Any]:
    bias = float(state["gate_bias_V"])
    export_dir = Path(state["export_dir"])
    case = {
        "group": "dd_idvg_deep_off",
        "bias_kind": "gate",
        "bias_V": bias,
        "gate_bias_V": bias,
        "drain_bias_V": 1.1,
        "export_dir": export_dir,
    }
    run_dir = OUTPUT / base.slug("gate", bias)
    run_dir.mkdir(parents=True, exist_ok=True)
    state_file = run_dir / "sentaurus_state_for_vela.csv"
    base.convert_state(export_dir, state_file)
    sg_csv, sg_status = base.run_probe(
        runner, case, run_dir, "sg_edge_flux_probe", state_file
    )
    mobility_csv, mobility_status = base.run_probe(
        runner, case, run_dir, "edge_mobility_probe", state_file
    )
    term_csv, term_status = base.run_probe(
        runner, case, run_dir, "newton_carrier_term_probe", state_file
    )
    residual_csv, residual_status = base.run_probe(
        runner, case, run_dir, "newton_residual_probe", state_file
    )
    feedback_dir = base.feedback_fields(state_file, run_dir / "feedback_fields")
    feedback_csv, feedback_status = base.run_probe(
        runner,
        case,
        run_dir,
        "newton_feedback_substitution_probe",
        state_file,
        feedback_dir,
    )
    self_dir = run_dir / "vela_self_consistent_state"
    self_dir.mkdir(parents=True, exist_ok=True)
    self_sg_csv, self_sg_status = base.run_probe(
        runner,
        case,
        self_dir,
        "sg_edge_flux_probe",
        VELA_STATE_FILES[round(bias, 2)],
    )

    sg_rows = read_csv(sg_csv)
    mobility_rows = read_csv(mobility_csv)
    term_rows = read_csv(term_csv)
    self_sg_rows = read_csv(self_sg_csv)
    silicon_nodes = set(base.field(export_dir, "srhRecombination", 3))
    fixed_srh = source_scale * sum(
        float(row["electron_recombination"])
        for row in term_rows
        if int(row["node_id"]) in silicon_nodes
    )
    sent_srh = sentaurus_srh_integral(base, export_dir)
    feedback = recombination_by_variant(
        read_csv(feedback_csv), silicon_nodes, source_scale
    )
    density_effect = (
        abs(feedback["density_only"] - feedback["baseline"])
        / max(abs(feedback["baseline"]), 1.0e-300)
    )
    tdr_current = contact_scalars(export_dir, "ContactCurrentFlux")
    tdr_voltage = contact_scalars(export_dir, "ContactExternalVoltage")
    edge_comparison = run_dir / "sentaurus_vela_edge_formula_comparison.csv"
    transport = base.edge_transport_metrics(
        export_dir, sg_rows, mobility_rows, edge_comparison
    )
    comparison_rows = read_csv(edge_comparison)
    max_sent_line = max(
        (float(row["sentaurus_eLineCurrent_A_m"]) for row in comparison_rows),
        default=0.0,
    )
    active_edge_ids = {
        int(row["edge_id"])
        for row in comparison_rows
        if float(row["sentaurus_eLineCurrent_A_m"]) > max_sent_line * 1.0e-3
    }
    transport["configured_quasi_fermi_drive"] = qf_drive_metrics(
        base, export_dir, mobility_rows, active_edge_ids
    )
    contact_cut = vela_contact_cut_currents(export_dir, sg_rows)
    mobility_substituted_cut = vela_contact_cut_currents(
        export_dir, substitute_sentaurus_nodal_mobility(base, export_dir, sg_rows)
    )
    result = {
        "gate_bias_V": bias,
        "drain_bias_V": 1.1,
        "sentaurus_tdr_contact_voltage_V": tdr_voltage,
        "sentaurus_tdr_contact_current_flux_A": tdr_current,
        "sentaurus_tdr_contact_kcl_A": sum(tdr_current.values()),
        "density": base.endpoint_density_metrics(export_dir, sg_rows),
        "transport": transport,
        "srh_shape": base.srh_shape_metrics(export_dir, term_rows),
        "srh_integral": {
            "sentaurus": sent_srh,
            "vela_formula_on_sentaurus_state_A_per_um": fixed_srh,
            "magnitude_ratio_vela_over_sentaurus": abs(fixed_srh)
            / max(abs(sent_srh["signed_A_per_um"]), 1.0e-300),
            "feedback_variants_A_per_um": feedback,
            "exact_sentaurus_density_relative_effect": density_effect,
            "carrier_term_to_A_per_um_scale": source_scale,
        },
        "vela_sg_contact_cut": contact_cut,
        "vela_sg_contact_cut_sentaurus_nodal_mobility": mobility_substituted_cut,
        "vela_sg_drain_cut_precision": contact_cut_precision_metrics(
            export_dir, sg_rows, "drain"
        ),
        "vela_self_consistent_sg_contact_cut": vela_contact_cut_currents(
            export_dir, self_sg_rows
        ),
        "vela_self_consistent_sg_drain_cut_precision": contact_cut_precision_metrics(
            export_dir, self_sg_rows, "drain"
        ),
        "statuses": {
            "sg": sg_status,
            "mobility": mobility_status,
            "carrier_terms": term_status,
            "residual": residual_status,
            "feedback": feedback_status,
            "vela_self_consistent_sg": self_sg_status,
        },
        "artifacts": {
            "sentaurus_state": str(state_file.resolve()),
            "sg_edges": str(sg_csv.resolve()),
            "edge_mobility": str(mobility_csv.resolve()),
            "carrier_terms": str(term_csv.resolve()),
            "residual": str(residual_csv.resolve()),
            "feedback": str(feedback_csv.resolve()),
            "vela_self_consistent_sg": str(self_sg_csv.resolve()),
        },
    }
    (run_dir / "summary.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    return result


def build_table(cases: list[dict[str, Any]]) -> list[dict[str, Any]]:
    sent_curve = read_csv(SENTAURUS_CURVE)
    sent_terminals = sentaurus_terminal_rows()
    vela_curve = read_csv(VELA_CURVE)
    vela_srh = read_csv(VELA_SRH)
    vela_seed_srh = read_csv(
        REF
        / "vela_baseline/dd_dg_fixed_contract_v1_2026-08-24/runs/dd"
        / "dd_idvg_curve_srh_balance.csv"
    )
    rows: list[dict[str, Any]] = []
    for case in cases:
        bias = float(case["gate_bias_V"])
        sent = at_bias(sent_curve, bias)
        sent_terminal = min(
            sent_terminals, key=lambda row: abs(float(row["bias_V"]) - bias)
        )
        vela = at_bias(vela_curve, bias)
        balance = at_bias(vela_srh, bias) or at_bias(vela_seed_srh, bias)
        sent_id = float(sent["current_total"]) if sent else math.nan
        vela_id = float(vela["current_total_A_per_um"]) if vela else math.nan
        t = case["transport"]
        rows.append(
            {
                "gate_bias_V": bias,
                "sentaurus_curve_Id_A_per_um": sent_id,
                "sentaurus_curve_source_A_per_um": sent_terminal["source_A"],
                "sentaurus_curve_substrate_A_per_um": sent_terminal["substrate_A"],
                "sentaurus_curve_gate_A_per_um": sent_terminal["gate_A"],
                "sentaurus_curve_kcl_A_per_um": sum(
                    sent_terminal[key]
                    for key in ("source_A", "substrate_A", "gate_A", "drain_A")
                ),
                "sentaurus_tdr_drain_flux_A": case["sentaurus_tdr_contact_current_flux_A"].get("drain", math.nan),
                "sentaurus_tdr_kcl_A": case["sentaurus_tdr_contact_kcl_A"],
                "sentaurus_tdr_drain_vs_curve_relative_difference": abs(
                    case["sentaurus_tdr_contact_current_flux_A"].get("drain", math.nan)
                    - sent_id
                ) / abs(sent_id),
                "vela_forward_Id_A_per_um": vela_id,
                "vela_forward_relative_error": abs(vela_id - sent_id) / abs(sent_id),
                "vela_forward_abs_error_dex": abs(math.log10(abs(vela_id)) - math.log10(abs(sent_id))),
                "vela_forward_kcl_ratio": float(balance["id_to_kcl_residual_ratio"]) if balance else math.nan,
                "sentaurus_srh_A_per_um": case["srh_integral"]["sentaurus"]["signed_A_per_um"],
                "vela_fixed_state_srh_A_per_um": case["srh_integral"]["vela_formula_on_sentaurus_state_A_per_um"],
                "vela_to_sentaurus_srh_ratio": case["srh_integral"]["magnitude_ratio_vela_over_sentaurus"],
                "exact_density_srh_relative_effect": case["srh_integral"]["exact_sentaurus_density_relative_effect"],
                "electron_density_p95_dex": case["density"]["electron_density_abs_error_dex"]["p95"],
                "electron_mobility_p95_dex": t["electron_mobility_abs_error_dex"]["p95"],
                "active_electron_mobility_p95_dex": t["current_carrying_edges"]["electron_mobility_abs_error_dex"]["p95"],
                "configured_qf_drive_p95_dex": t["configured_quasi_fermi_drive"]["all_edges_abs_error_dex"]["p95"],
                "active_configured_qf_drive_p95_dex": t["configured_quasi_fermi_drive"]["current_carrying_edges_abs_error_dex"]["p95"],
                "electron_sg_current_p95_dex": t["electron_sg_line_current_abs_error_dex"]["p95"],
                "active_electron_sg_current_p95_dex": t["current_carrying_edges"]["electron_sg_line_current_abs_error_dex"]["p95"],
                "srh_shape_tv": case["srh_shape"]["normalized_total_variation_distance"],
                "vela_sg_drain_cut_A_per_um": case["vela_sg_contact_cut"]["drain"]["total_A_per_um"],
                "vela_sg_drain_cut_sentaurus_mobility_A_per_um": case["vela_sg_contact_cut_sentaurus_nodal_mobility"]["drain"]["total_A_per_um"],
                "drain_cut_median_qf_drop_V": case["vela_sg_drain_cut_precision"]["median_abs_qf_drop_V"],
                "drain_cut_median_qf_drop_ulp": case["vela_sg_drain_cut_precision"]["median_abs_qf_drop_ulp"],
                "drain_cut_edge_cancellation_condition": case["vela_sg_drain_cut_precision"]["edge_cancellation_condition"],
                "drain_cut_largest_edge_to_net_ratio": case["vela_sg_drain_cut_precision"]["largest_edge_to_net_ratio"],
                "vela_self_consistent_sg_drain_cut_A_per_um": case["vela_self_consistent_sg_contact_cut"]["drain"]["total_A_per_um"],
                "vela_self_consistent_drain_cut_median_qf_drop_V": case["vela_self_consistent_sg_drain_cut_precision"]["median_abs_qf_drop_V"],
                "vela_self_consistent_drain_cut_median_qf_drop_ulp": case["vela_self_consistent_sg_drain_cut_precision"]["median_abs_qf_drop_ulp"],
                "vela_self_consistent_drain_cut_edge_cancellation_condition": case["vela_self_consistent_sg_drain_cut_precision"]["edge_cancellation_condition"],
            }
        )
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, default=DEFAULT_RUNNER)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.check:
        report = json.loads(SUMMARY.read_text(encoding="utf-8"))
        for case in report["cases"]:
            for artifact in case["artifacts"].values():
                if not Path(artifact).is_file():
                    raise FileNotFoundError(artifact)
        if len(read_csv(TABLE)) != 3:
            raise RuntimeError("Expected three DD deep-off summary rows")
        print("TransportModels DD deep-off formula replay check: PASS")
        return 0

    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    legacy = json.loads(LEGACY_FIXED_STATE.read_text(encoding="utf-8"))
    source_scale = float(
        legacy["srh_integral_A_per_um"]["carrier_term_to_A_per_um_scale"]
    )
    base = load_base()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    cases = [
        execute_case(base, args.runner, state, source_scale)
        for state in manifest["dd_deep_off_states"]
    ]
    table = build_table(cases)
    base.write_csv(TABLE, table)
    report = {
        "schema": "vela.transportmodels_dd_deep_off_formula_replay.v1",
        "status": "complete",
        "method": "immutable Sentaurus states replayed through Vela production density, mobility, SRH, SG, residual, and feedback probes",
        "source_scale_note": "Physical SRH current scale is inherited from the previously calibrated identical 3315-node mesh and unit-scaling contract.",
        "cases": cases,
        "summary_table": str(TABLE.resolve()),
    }
    SUMMARY.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "complete", "cases": len(cases), "table": str(TABLE)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
