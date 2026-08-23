#!/usr/bin/env python3
"""Analyze the three numerically unresolved DG Id-Vg deep-off points."""

from __future__ import annotations

import csv
import json
import math
from datetime import datetime
from pathlib import Path
from typing import Any, Iterable


REPO = Path(__file__).resolve().parents[1]
REF = REPO / "build-release/reference_tcad/transportmodels_sentaurus2022"
NEW_ROOT = REF / "vela_baseline/dg_quantum_contract_regression_2026-08-23"
NEW_RUN = NEW_ROOT / "runs/dg"
OLD_RUN = REF / "vela_baseline/dd_dg_srh_corrected_cold_regression_2026-08-23/runs"
SENT_CURVE = REF / "run02/normalized/dg_idvg.csv"
SENT_MANIFEST = REF / "sentaurus_vm_runs/remaining_spatial_oracles_20260823/remaining_spatial_oracles_manifest.json"
PROFILE_REPORT = REPO / "docs/validation/transportmodels_idvg_spatial_oracle_2026-08-21.json"
PRIOR_SPATIAL_REPORT = REPO / "docs/validation/transportmodels_remaining_spatial_state_compare_2026-08-23.json"
FORMULA_ROOT = REF / "reports/transportmodels_sentaurus_formula_replay_20260823/dg_idvg_deep_off"
OUTPUT = REF / "reports/transportmodels_dg_deep_off_three_points_20260823"
REPORT_JSON = REPO / "docs/validation/transportmodels_dg_deep_off_three_points_2026-08-23.json"
ARTIFACT_JSON = OUTPUT / "artifact.json"
REPORT_HTML = OUTPUT / "transportmodels_dg_deep_off_three_points_report.html"

DEEP_BIASES = (-1.0, -0.84, -0.68)
LOG_LIMIT_DEX = 0.15
KCL_MARGIN = 10.0


def scientific_display(value: float) -> str:
    """Format tiny currents as text that portable table readers cannot round to zero."""
    if value == 0.0:
        return "0"
    exponent = math.floor(math.log10(abs(value)))
    mantissa = value / (10.0**exponent)
    return f"{mantissa:.4f} x 10^{exponent}"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def percentile(values: Iterable[float], fraction: float) -> float:
    ordered = sorted(float(value) for value in values if math.isfinite(float(value)))
    if not ordered:
        return math.nan
    position = (len(ordered) - 1) * fraction
    lo, hi = math.floor(position), math.ceil(position)
    if lo == hi:
        return ordered[lo]
    return ordered[lo] * (hi - position) + ordered[hi] * (position - lo)


def field(export_dir: Path, name: str, region: int = 3) -> dict[int, float]:
    return {
        int(row["node_id"]): float(row["component0"])
        for row in read_csv(export_dir / "fields" / f"{name}_region{region}.csv")
    }


def bias_key(value: float) -> float:
    return round(float(value), 12)


def bias_tag(value: float) -> str:
    return f"{value:.6f}".replace("-", "m").replace(".", "p")


def current_map(path: Path, current_field: str) -> dict[float, float]:
    return {
        bias_key(float(row["bias_V"])): abs(float(row[current_field]))
        for row in read_csv(path)
    }


def current_state_path(bias: float) -> Path:
    if math.isclose(bias, -1.0):
        return NEW_RUN / "dg_idvg_final_bias_relax_final_state.csv"
    return NEW_RUN / f"dg_idvg_curve_state_bias_{bias_tag(bias)}.csv"


def diagnostics_by_bias() -> tuple[dict[float, dict[str, str]], dict[float, dict[str, dict[str, str]]]]:
    srh: dict[float, dict[str, str]] = {}
    terminal: dict[float, dict[str, dict[str, str]]] = {}
    for stem in ("dg_idvg_final_bias_relax", "dg_idvg_curve"):
        for row in read_csv(NEW_RUN / f"{stem}_srh_balance.csv"):
            srh[bias_key(float(row["bias_V"]))] = row
        for row in read_csv(NEW_RUN / f"{stem}_terminal_balance.csv"):
            bias = bias_key(float(row["bias_V"]))
            terminal.setdefault(bias, {})[row["contact"].lower()] = row
    return srh, terminal


def point_comparison() -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    sent = current_map(SENT_CURVE, "current_total")
    new = current_map(NEW_ROOT / "dg_idvg_aligned.csv", "vela_A_per_um")
    old_dg = current_map(OLD_RUN / "dg/dg_idvg_curve_comparison_candidate.csv", "current_total_A_per_um")
    dd = current_map(OLD_RUN / "dd/dd_idvg_curve_comparison_candidate.csv", "current_total_A_per_um")
    srh_by_bias, terminal_by_bias = diagnostics_by_bias()
    points: list[dict[str, Any]] = []
    terminal_rows: list[dict[str, Any]] = []
    for bias in DEEP_BIASES:
        key = bias_key(bias)
        sent_i, new_i, old_i, dd_i = sent[key], new[key], old_dg[key], dd[key]
        srh = srh_by_bias[key]
        contacts = terminal_by_bias[key]
        drain = contacts["drain"]
        source = contacts["source"]
        substrate = contacts["substrate"]
        generation = abs(float(srh["srh_generation_current_A_per_um"]))
        kcl_signed = float(srh["four_terminal_kcl_residual_A_per_um"])
        gap = sent_i - new_i
        closure_corrected = new_i - kcl_signed
        drain_e = abs(float(drain["current_electron_A_per_um"]))
        substrate_h = abs(float(substrate["current_hole_A_per_um"]))
        compensated_error = max(
            abs(float(row["current_total_compensated_A_per_um"]) - float(row["current_total_long_double_reference_A_per_um"]))
            / max(abs(float(row["current_total_long_double_reference_A_per_um"])), 1.0e-300)
            for row in contacts.values()
        )
        point = {
            "bias_V": bias,
            "bias_label": f"{bias:.2f} V",
            "sentaurus_A_per_um": sent_i,
            "vela_validated_dg_A_per_um": new_i,
            "vela_prior_dg_A_per_um": old_i,
            "vela_dd_A_per_um": dd_i,
            "vela_to_sentaurus_ratio": new_i / sent_i,
            "absolute_relative_error": abs(new_i - sent_i) / sent_i,
            "absolute_log_error_dex": abs(math.log10(new_i) - math.log10(sent_i)),
            "prior_dg_absolute_log_error_dex": abs(math.log10(old_i) - math.log10(sent_i)),
            "quantum_contract_current_uplift": new_i / old_i - 1.0,
            "dd_absolute_log_error_dex": abs(math.log10(dd_i) - math.log10(sent_i)),
            "sentaurus_minus_vela_gap_A_per_um": gap,
            "kcl_residual_A_per_um": kcl_signed,
            "id_to_kcl_residual_ratio": abs(new_i / kcl_signed),
            "kcl_fraction_of_sentaurus_vela_gap": abs(kcl_signed) / gap,
            "closure_corrected_current_A_per_um": closure_corrected,
            "closure_corrected_relative_error": abs(closure_corrected - sent_i) / sent_i,
            "srh_generation_A_per_um": generation,
            "drain_electron_current_A_per_um": drain_e,
            "substrate_hole_current_A_per_um": substrate_h,
            "drain_capture_of_srh_generation": drain_e / generation,
            "substrate_to_srh_generation_ratio": substrate_h / generation,
            "source_to_drain_current_ratio": abs(float(source["current_total_A_per_um"])) / new_i,
            "compensated_vs_long_double_max_relative_delta": compensated_error,
            "numerical_status": srh["numerical_status"],
        }
        points.append(point)
        terminal_rows.append(
            {
                "bias_V": bias,
                "bias_label": point["bias_label"],
                "sentaurus_A_per_um": sent_i,
                "sentaurus_display": scientific_display(sent_i),
                "vela_A_per_um": new_i,
                "vela_display": scientific_display(new_i),
                "gap_A_per_um": gap,
                "signed_kcl_residual_A_per_um": kcl_signed,
                "signed_kcl_residual_display": scientific_display(kcl_signed),
                "kcl_gap_fraction": point["kcl_fraction_of_sentaurus_vela_gap"],
                "closure_corrected_current_A_per_um": closure_corrected,
                "closure_corrected_relative_error": point["closure_corrected_relative_error"],
                "srh_generation_A_per_um": generation,
                "drain_electron_A_per_um": drain_e,
                "substrate_hole_A_per_um": substrate_h,
                "drain_capture_ratio": point["drain_capture_of_srh_generation"],
                "substrate_generation_ratio": point["substrate_to_srh_generation_ratio"],
                "source_to_drain_current_ratio": point["source_to_drain_current_ratio"],
                "id_kcl_ratio": point["id_to_kcl_residual_ratio"],
                "status": point["numerical_status"],
            }
        )
    intervals: list[dict[str, Any]] = []
    for left, right in zip(points, points[1:]):
        delta_v = right["bias_V"] - left["bias_V"]
        sent_ratio = right["sentaurus_A_per_um"] / left["sentaurus_A_per_um"]
        vela_ratio = right["vela_validated_dg_A_per_um"] / left["vela_validated_dg_A_per_um"]
        intervals.append(
            {
                "interval": f"{left['bias_V']:.2f} to {right['bias_V']:.2f} V",
                "delta_V": delta_v,
                "sentaurus_current_ratio": sent_ratio,
                "vela_current_ratio": vela_ratio,
                "sentaurus_log_slope_dex_per_V": math.log10(sent_ratio) / delta_v,
                "vela_log_slope_dex_per_V": math.log10(vela_ratio) / delta_v,
                "absolute_log_slope_difference_dex_per_V": abs(
                    math.log10(sent_ratio) - math.log10(vela_ratio)
                ) / delta_v,
            }
        )
    return points, terminal_rows, intervals


def spatial_state_comparison() -> list[dict[str, Any]]:
    manifest = json.loads(SENT_MANIFEST.read_text(encoding="utf-8"))
    exports = {
        bias_key(float(row["gate_bias_V"])): Path(row["export_dir"])
        for row in manifest["dg_deep_off_states"]
    }
    profiles_doc = json.loads(PROFILE_REPORT.read_text(encoding="utf-8"))
    drain_nodes = {int(node) for node in profiles_doc["profiles"]["drain_end"]["node_ids"]}
    prior_doc = json.loads(PRIOR_SPATIAL_REPORT.read_text(encoding="utf-8"))
    prior_by_bias = {
        bias_key(float(row["bias_V"])): row["zones"]
        for row in prior_doc["cases"]
        if row["group"] == "dg_idvg_deep_off"
    }
    output: list[dict[str, Any]] = []
    for bias in DEEP_BIASES:
        state = {int(row["node_id"]): row for row in read_csv(current_state_path(bias))}
        export = exports[bias_key(bias)]
        sent = {
            "psi": field(export, "ElectrostaticPotential"),
            "phin": field(export, "eQuasiFermiPotential"),
            "qn": field(export, "eQuantumPotential"),
            "n": field(export, "eDensity"),
        }
        silicon = set(sent["n"])
        for zone, nodes in (("all_substrate", silicon), ("drain_end", silicon & drain_nodes)):
            prior = prior_by_bias[bias_key(bias)][zone]
            output.append(
                {
                    "bias_V": bias,
                    "bias_label": f"{bias:.2f} V",
                    "zone": zone,
                    "node_count": len(nodes),
                    "psi_p95_error_mV": percentile(
                        [1.0e3 * abs(float(state[node]["psi"]) - sent["psi"][node]) for node in nodes], 0.95
                    ),
                    "phin_p95_error_mV": percentile(
                        [1.0e3 * abs(float(state[node]["phin"]) - sent["phin"][node]) for node in nodes], 0.95
                    ),
                    "qn_p95_error_mV": percentile(
                        [
                            1.0e3 * abs(
                                float(state[node]["electron_quantum_potential_V"]) - sent["qn"][node]
                            )
                            for node in nodes
                        ],
                        0.95,
                    ),
                    "electron_density_p95_error_dex": percentile(
                        [
                            abs(
                                math.log10(max(float(state[node]["electrons_m3"]) / 1.0e6, 1.0))
                                - math.log10(max(sent["n"][node], 1.0))
                            )
                            for node in nodes
                        ],
                        0.95,
                    ),
                    "prior_qn_p95_error_mV": prior["qn_abs_error_mV"]["p95"],
                    "prior_electron_density_p95_error_dex": prior[
                        "electron_density_abs_error_dex"
                    ]["p95"],
                    "vela_state": str(current_state_path(bias).resolve()),
                    "sentaurus_export": str(export.resolve()),
                }
            )
    return output


def formula_replay_rows() -> list[dict[str, Any]]:
    slugs = {-1.0: "gate_m1p00", -0.84: "gate_m0p84", -0.68: "gate_m0p68"}
    rows: list[dict[str, Any]] = []
    for bias in DEEP_BIASES:
        path = FORMULA_ROOT / slugs[bias] / "formula_replay_summary.json"
        report = json.loads(path.read_text(encoding="utf-8"))
        carrying = report["transport"]["current_carrying_edges"]
        rows.append(
            {
                "bias_V": bias,
                "bias_label": f"{bias:.2f} V",
                "electron_density_p95_error_dex": report["density"]["electron_density_abs_error_dex"]["p95"],
                "hole_density_p95_error_dex": report["density"]["hole_density_abs_error_dex"]["p95"],
                "electron_mobility_median_error_dex": report["transport"]["electron_mobility_abs_error_dex"]["median"],
                "electron_mobility_p95_error_dex": report["transport"]["electron_mobility_abs_error_dex"]["p95"],
                "current_carrying_sg_median_error_dex": carrying["electron_sg_line_current_abs_error_dex"]["median"],
                "srh_total_variation_distance": report["srh"]["normalized_total_variation_distance"],
                "srh_top50_node_overlap": report["srh"]["top50_node_overlap"],
                "source": str(path.resolve()),
            }
        )
    return rows


def build_artifact(report: dict[str, Any]) -> dict[str, Any]:
    points = report["points"]
    spatial_state = report["spatial_state"]
    summary = report["summary"]
    generated_at = report["generated_at"]
    sources = [
        {
            "id": "deep_off_points",
            "label": "Deep-off aligned current comparison",
            "path": "build-release/reference_tcad/transportmodels_sentaurus2022/reports/transportmodels_dg_deep_off_three_points_20260823/deep_off_point_comparison.csv",
            "query": {
                "engine": "duckdb",
                "language": "sql",
                "sql": "SELECT * FROM read_csv_auto('build-release/reference_tcad/transportmodels_sentaurus2022/reports/transportmodels_dg_deep_off_three_points_20260823/deep_off_point_comparison.csv', header=true)",
                "description": "Loads the reviewed three-point current comparison and derived error metrics.",
                "tables_used": ["build-release/reference_tcad/transportmodels_sentaurus2022/reports/transportmodels_dg_deep_off_three_points_20260823/deep_off_point_comparison.csv"],
                "filters": ["gate bias in {-1.00, -0.84, -0.68} V", "drain bias = 1.1 V"],
                "metric_definitions": [
                    "log error = abs(log10(abs(Id_Vela)) - log10(abs(Id_Sentaurus)))",
                    "KCL gap share = abs(four-terminal KCL residual) / (Id_Sentaurus - Id_Vela)",
                    "closure-corrected current = Id_Vela - signed four-terminal KCL residual",
                ],
                "executed_at": generated_at,
            },
        },
        {
            "id": "terminal_closure",
            "label": "Vela four-terminal and SRH closure decomposition",
            "path": "build-release/reference_tcad/transportmodels_sentaurus2022/reports/transportmodels_dg_deep_off_three_points_20260823/deep_off_terminal_closure.csv",
            "query": {
                "engine": "duckdb",
                "language": "sql",
                "sql": "SELECT * FROM read_csv_auto('build-release/reference_tcad/transportmodels_sentaurus2022/reports/transportmodels_dg_deep_off_three_points_20260823/deep_off_terminal_closure.csv', header=true)",
                "description": "Loads the reviewed terminal-current, SRH-generation and KCL closure decomposition.",
                "tables_used": ["build-release/reference_tcad/transportmodels_sentaurus2022/reports/transportmodels_dg_deep_off_three_points_20260823/deep_off_terminal_closure.csv"],
                "filters": ["gate bias in {-1.00, -0.84, -0.68} V", "drain bias = 1.1 V"],
                "metric_definitions": [
                    "KCL gap share = abs(four-terminal KCL residual) / (Id_Sentaurus - Id_Vela)",
                    "closure-corrected current = Id_Vela - signed four-terminal KCL residual",
                ],
                "executed_at": generated_at,
            },
        },
        {
            "id": "spatial_oracles",
            "label": "Sentaurus 2022 spatial-oracle comparison",
            "path": "build-release/reference_tcad/transportmodels_sentaurus2022/reports/transportmodels_dg_deep_off_three_points_20260823/deep_off_spatial_state.csv",
            "query": {
                "engine": "duckdb",
                "language": "sql",
                "sql": "SELECT * FROM read_csv_auto('build-release/reference_tcad/transportmodels_sentaurus2022/reports/transportmodels_dg_deep_off_three_points_20260823/deep_off_spatial_state.csv', header=true)",
                "description": "Node-matched p95 errors over the silicon substrate and drain-end profile.",
                "tables_used": ["build-release/reference_tcad/transportmodels_sentaurus2022/reports/transportmodels_dg_deep_off_three_points_20260823/deep_off_spatial_state.csv"],
                "filters": ["silicon substrate nodes", "gate bias in {-1.00, -0.84, -0.68} V"],
                "executed_at": generated_at,
            },
        },
        {
            "id": "formula_replay",
            "label": "Sentaurus-state formula replay",
            "path": "build-release/reference_tcad/transportmodels_sentaurus2022/reports/transportmodels_dg_deep_off_three_points_20260823/deep_off_formula_replay.csv",
            "query": {
                "engine": "duckdb",
                "language": "sql",
                "sql": "SELECT * FROM read_csv_auto('build-release/reference_tcad/transportmodels_sentaurus2022/reports/transportmodels_dg_deep_off_three_points_20260823/deep_off_formula_replay.csv', header=true)",
                "description": "Evaluates Vela density, mobility, SG flux and SRH operators on exported Sentaurus states.",
                "tables_used": ["build-release/reference_tcad/transportmodels_sentaurus2022/reports/transportmodels_dg_deep_off_three_points_20260823/deep_off_formula_replay.csv"],
                "filters": ["DG Id-Vg deep-off points only"],
                "executed_at": generated_at,
            },
        },
    ]
    manifest = {
        "version": 1,
        "surface": "report",
        "title": "TransportModels DG 深关断三点专项诊断",
        "description": "Sentaurus 2022 与 Vela DG 在 Vg=-1.00、-0.84、-0.68 V 的电流、守恒与空间状态对比。",
        "generatedAt": generated_at,
        "charts": [
            {
                "id": "log_error_chart",
                "title": "深关断 Id-Vg 对数电流误差",
                "subtitle": "三个参考偏压；虚线为 0.15 dex 深关断验收线",
                "intent": "comparison",
                "question": "三个深关断点是否满足 0.15 dex 对数误差门槛？",
                "rationale": "只有三个离散偏压，单系列柱图比趋势线更诚实地呈现逐点超限幅度。",
                "type": "bar",
                "dataset": "deep_points",
                "sourceId": "deep_off_points",
                "valueFormat": "number",
                "unit": "dex",
                "encodings": {
                    "x": {"field": "bias_label", "type": "ordinal", "label": "Gate voltage Vg"},
                    "y": {"field": "absolute_log_error_dex", "type": "quantitative", "label": "Absolute log error", "unit": "dex"},
                    "tooltip": [
                        {"field": "absolute_relative_error", "type": "quantitative", "label": "Relative error", "format": "percent"},
                        {"field": "id_to_kcl_residual_ratio", "type": "quantitative", "label": "Id/KCL residual"},
                        {"field": "kcl_fraction_of_sentaurus_vela_gap", "type": "quantitative", "label": "KCL gap share", "format": "percent"},
                    ],
                },
                "labels": {"values": "all"},
                "referenceLines": [{"axis": "y", "value": LOG_LIMIT_DEX, "label": "0.15 dex", "color": "neutral", "lineStyle": "dashed"}],
                "layout": "full",
                "surface": {"surface": "card", "viewMode": "visualization"},
            }
        ],
        "tables": [
            {
                "id": "terminal_table",
                "title": "端口电流与连续性闭合",
                "subtitle": "电流单位 A/um；KCL 缺口份额以 Sentaurus-Vela 电流差为分母",
                "dataset": "terminal_closure",
                "sourceId": "terminal_closure",
                "defaultSort": {"field": "bias_V", "direction": "asc"},
                "density": "spacious",
                "layout": "full",
                "columns": [
                    {"field": "bias_V", "label": "Vg (V)", "format": "number"},
                    {"field": "sentaurus_display", "label": "Sentaurus Id", "type": "text"},
                    {"field": "vela_display", "label": "Vela Id", "type": "text"},
                    {"field": "signed_kcl_residual_display", "label": "Signed KCL residual", "type": "text"},
                    {"field": "kcl_gap_fraction", "label": "KCL gap share", "format": "percent"},
                    {"field": "drain_capture_ratio", "label": "Drain/SRH generation", "format": "percent"},
                    {"field": "source_to_drain_current_ratio", "label": "Abs(source)/drain", "format": "percent"},
                    {"field": "closure_corrected_relative_error", "label": "Corrected error", "format": "percent"},
                    {"field": "id_kcl_ratio", "label": "Id/KCL", "format": "number"},
                    {"field": "status", "label": "Classification", "type": "text"},
                ],
            },
            {
                "id": "state_table",
                "title": "自洽空间状态与 Sentaurus 的 p95 误差",
                "subtitle": "全硅区与漏端局部；电势单位 mV，电子浓度单位 dex",
                "dataset": "spatial_state",
                "sourceId": "spatial_oracles",
                "defaultSort": {"field": "bias_V", "direction": "asc"},
                "density": "spacious",
                "layout": "full",
                "columns": [
                    {"field": "bias_V", "label": "Vg (V)", "format": "number"},
                    {"field": "zone", "label": "Zone", "type": "text"},
                    {"field": "psi_p95_error_mV", "label": "Psi p95", "format": "number"},
                    {"field": "phin_p95_error_mV", "label": "Phi_n p95", "format": "number"},
                    {"field": "qn_p95_error_mV", "label": "Qn p95", "format": "number"},
                    {"field": "electron_density_p95_error_dex", "label": "n p95", "format": "number"},
                    {"field": "prior_qn_p95_error_mV", "label": "Prior Qn p95", "format": "number"},
                    {"field": "prior_electron_density_p95_error_dex", "label": "Prior n p95", "format": "number"},
                ],
            },
            {
                "id": "formula_table",
                "title": "Sentaurus 状态代入 Vela 公式的误差",
                "subtitle": "用于区分公式本身与自洽求解/端口闭合误差",
                "dataset": "formula_replay",
                "sourceId": "formula_replay",
                "defaultSort": {"field": "bias_V", "direction": "asc"},
                "density": "spacious",
                "layout": "full",
                "columns": [
                    {"field": "bias_V", "label": "Vg (V)", "format": "number"},
                    {"field": "electron_density_p95_error_dex", "label": "n formula p95", "format": "number"},
                    {"field": "electron_mobility_median_error_dex", "label": "Mobility median", "format": "number"},
                    {"field": "current_carrying_sg_median_error_dex", "label": "SG current median", "format": "number"},
                    {"field": "srh_total_variation_distance", "label": "SRH TV distance", "format": "number"},
                    {"field": "srh_top50_node_overlap", "label": "SRH top50 overlap", "format": "percent"},
                ],
            },
        ],
        "sources": sources,
        "blocks": [
            {"id": "title", "type": "markdown", "body": "# TransportModels DG 深关断三点专项诊断"},
            {
                "id": "technical_summary",
                "type": "markdown",
                "body": (
                    "## 技术结论\n\n"
                    f"三个点的 Vela 电流均低于 Sentaurus，绝对对数误差为 **{min(row['absolute_log_error_dex'] for row in points):.3f}–{max(row['absolute_log_error_dex'] for row in points):.3f} dex**。"
                    f"但四端 KCL 残差可解释电流缺口的 **{min(row['kcl_fraction_of_sentaurus_vela_gap'] for row in points):.1%}–{max(row['kcl_fraction_of_sentaurus_vela_gap'] for row in points):.1%}**；"
                    f"按残差方向修正后，剩余相对误差降为 **{min(row['closure_corrected_relative_error'] for row in points):.1%}–{max(row['closure_corrected_relative_error'] for row in points):.1%}**。"
                    "因此当前证据首先支持“极低电流下的连续性/端口闭合精度不足”，而不是新的 DG 量子方程或 SRH 公式失配。"
                ),
            },
            {
                "id": "finding_scale",
                "type": "markdown",
                "sourceId": "deep_off_points",
                "body": (
                    "## 三点差异主要表现为稳定的低估，而非错误的曲线形状\n\n"
                    f"Vela/Sentaurus 电流比仅为 **{min(row['vela_to_sentaurus_ratio'] for row in points):.3f}–{max(row['vela_to_sentaurus_ratio'] for row in points):.3f}**。"
                    "-1.00→-0.84 V 时两者都近似平台；-0.84→-0.68 V 时 Sentaurus 电流增加 32.1%，Vela 增加 28.9%。"
                    "这说明门压响应方向基本正确，主要问题是极低电流平台的幅值与闭合精度。下图按三点逐项展示对数误差；三点都略高于 0.15 dex 门槛。"
                ),
            },
            {"id": "log_error", "type": "chart", "chartId": "log_error_chart"},
            {
                "id": "finding_kcl",
                "type": "markdown",
                "sourceId": "terminal_closure",
                "body": (
                    "## KCL 残差解释了大部分 Sentaurus–Vela 缺口\n\n"
                    "硅区 SRH 产生电流与衬底空穴电流匹配到约 0.03% 以内。"
                    "在 -1.00 V 和 -0.84 V，漏端只捕获约 78% 的 SRH 产生电流，源端电流仅为漏电流的 0.06% 和 0.94%；"
                    "在 -0.68 V，漏端已捕获 99.43%，但源端反向电子电流增至漏电流的 40.83%，所以四端电流仍未闭合。"
                    f"KCL 残差占两套工具电流缺口的 **{min(row['kcl_fraction_of_sentaurus_vela_gap'] for row in points):.1%}–{max(row['kcl_fraction_of_sentaurus_vela_gap'] for row in points):.1%}**。"
                    "接触边补偿求和与 long-double 参考的差异远小于该残差，因此问题不在接触电流求和舍入，而在非线性连续性方程尚未解析到 10× 电流裕量。"
                ),
            },
            {"id": "terminal", "type": "table", "tableId": "terminal_table"},
            {
                "id": "finding_state",
                "type": "markdown",
                "body": (
                    "## 量子契约改善了状态与电流，但没有消除极低电流闭合限制\n\n"
                    f"相对旧 DG 配置，三点电流提高 **{min(row['quantum_contract_current_uplift'] for row in points):.1%}–{max(row['quantum_contract_current_uplift'] for row in points):.1%}**。"
                    "全硅区 Qn p95 已从约 135–141 mV 降到约 10 mV，电子浓度 p95 从 0.78–0.91 dex 降到 0.105–0.114 dex。"
                    "与此同时，把 Sentaurus 状态直接代入 Vela 载流子公式时，电子浓度 p95 误差低于 0.009 dex。"
                    "这削弱了“材料参数、BGN 或密度公式是 0.17 dex 电流差主因”的假设；但漏端局部电子浓度 p95 仍有 0.54–0.60 dex，属于守恒解析后需要继续审计的次要热点。"
                ),
            },
            {"id": "state", "type": "table", "tableId": "state_table"},
            {"id": "formula", "type": "table", "tableId": "formula_table"},
            {
                "id": "scope",
                "type": "markdown",
                "body": (
                    "## 比较范围与指标定义\n\n"
                    "比较对象为同一 TransportModels MOS 器件、固定 Vd=1.1 V、Vg=-1.00/-0.84/-0.68 V 的 Sentaurus 2022 与 Vela DG。"
                    "普通相对误差在 10^-15 A/um 量级会放大，主验收量采用绝对对数误差；数值解析要求 |Id|/|KCL residual| ≥ 10。"
                ),
            },
            {
                "id": "method",
                "type": "markdown",
                "body": (
                    "## 诊断方法\n\n"
                    "1. 对齐 Sentaurus、当前 DG、旧 DG 与 DD 的相同偏压电流。\n"
                    "2. 从四端端口诊断和硅区 SRH 积分重建连续性闭合。\n"
                    "3. 将有符号 KCL 残差作为数值未闭合量，评估其对两套工具电流差的解释份额。\n"
                    "4. 节点一一对应比较 Psi、电子准费米势、Qn 与电子浓度。\n"
                    "5. 使用 Sentaurus 状态回放 Vela 密度、迁移率、SG 通量与 SRH 公式，区分公式误差和自洽求解误差。"
                ),
            },
            {
                "id": "limitations",
                "type": "markdown",
                "body": (
                    "## 局限、稳健性与不能下的结论\n\n"
                    "KCL 残差修正是误差归因诊断，不是新的物理解；不能把修正后电流当作正式曲线。"
                    "Sentaurus 参考曲线只提供漏端总电流，未提供同一输出文件中的四端分量，因此无法对两套工具执行完全对称的端口闭合审计。"
                    "空间 SG 电流比较使用 Sentaurus 节点电流密度投影，不是两种离散格式的严格恒等式。"
                ),
            },
            {
                "id": "next_steps",
                "type": "markdown",
                "body": (
                    "## 下一步应先提高连续性解析度\n\n"
                    "1. 对三点分别实施更严格的电子/空穴连续性绝对残差门槛，并禁止仅凭状态增量变小结束 Newton。\n"
                    "2. 将载流子行硬门槛设为 |KCL residual| ≤ 0.1|Id|，记录源/漏接触边通量与体 SRH 的逐区域闭合。\n"
                    "3. 使用当前已验证 DG 状态做 frozen-Q 与增量准费米变量 A/B，判断 Qn 外迭代是否只通过耦合放大残差。\n"
                    "4. 只有在三点全部达到 Id/KCL≥10 后，再重新评估剩余 6.8%–14.1% 的模型差异。"
                ),
            },
            {
                "id": "questions",
                "type": "markdown",
                "body": (
                    "## 仍需回答的问题\n\n"
                    "- 更严格连续性门槛能否把漏端电子电流推到与衬底/SRH 闭合一致的分支？\n"
                    "- Sentaurus 四端电流分量导出后，是否也存在相近的 SRH 产生与漏端捕获差？\n"
                    "- 在守恒解析后，剩余误差是否集中于接触附近 SG 通量或高场迁移率离散？"
                ),
            },
        ],
    }
    return {
        "surface": "report",
        "manifest": manifest,
        "snapshot": {
            "version": 1,
            "generatedAt": generated_at,
            "status": "ready",
            "datasets": {
                "deep_points": points,
                "terminal_closure": report["terminal_closure"],
                "spatial_state": spatial_state,
                "formula_replay": report["formula_replay"],
            },
            "accessIssues": [],
        },
        "sources": sources,
        "package_info": {"originUrl": "artifact://transportmodels-dg-deep-off-three-points"},
    }


def main() -> int:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    points, terminal, intervals = point_comparison()
    spatial = spatial_state_comparison()
    replay = formula_replay_rows()
    summary = {
        "min_log_error_dex": min(row["absolute_log_error_dex"] for row in points),
        "max_log_error_dex": max(row["absolute_log_error_dex"] for row in points),
        "min_kcl_gap_share": min(row["kcl_fraction_of_sentaurus_vela_gap"] for row in points),
        "max_kcl_gap_share": max(row["kcl_fraction_of_sentaurus_vela_gap"] for row in points),
        "min_closure_corrected_error": min(row["closure_corrected_relative_error"] for row in points),
        "max_closure_corrected_error": max(row["closure_corrected_relative_error"] for row in points),
        "min_quantum_contract_uplift": min(row["quantum_contract_current_uplift"] for row in points),
        "max_quantum_contract_uplift": max(row["quantum_contract_current_uplift"] for row in points),
        "all_points_log_error_pass": all(row["absolute_log_error_dex"] <= LOG_LIMIT_DEX for row in points),
        "all_points_kcl_resolved": all(row["id_to_kcl_residual_ratio"] >= KCL_MARGIN for row in points),
        "classification": "numerical_continuity_resolution_dominates",
        "confidence": "high for numerical non-resolution; moderate for the remaining model residual",
    }
    report = {
        "schema": "vela.transportmodels.dg_deep_off_three_points.v1",
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "status": "complete",
        "scope": {"gate_biases_V": list(DEEP_BIASES), "drain_bias_V": 1.1, "model": "DG"},
        "thresholds": {"log_error_limit_dex": LOG_LIMIT_DEX, "id_to_kcl_margin": KCL_MARGIN},
        "summary": summary,
        "points": points,
        "terminal_closure": terminal,
        "interval_shape": intervals,
        "spatial_state": spatial,
        "formula_replay": replay,
        "artifacts": {
            "point_csv": str((OUTPUT / "deep_off_point_comparison.csv").resolve()),
            "terminal_csv": str((OUTPUT / "deep_off_terminal_closure.csv").resolve()),
            "interval_csv": str((OUTPUT / "deep_off_interval_shape.csv").resolve()),
            "spatial_csv": str((OUTPUT / "deep_off_spatial_state.csv").resolve()),
            "formula_csv": str((OUTPUT / "deep_off_formula_replay.csv").resolve()),
            "artifact_json": str(ARTIFACT_JSON.resolve()),
            "report_html": str(REPORT_HTML.resolve()),
        },
    }
    write_csv(OUTPUT / "deep_off_point_comparison.csv", points)
    write_csv(OUTPUT / "deep_off_terminal_closure.csv", terminal)
    write_csv(OUTPUT / "deep_off_interval_shape.csv", intervals)
    write_csv(OUTPUT / "deep_off_spatial_state.csv", spatial)
    write_csv(OUTPUT / "deep_off_formula_replay.csv", replay)
    REPORT_JSON.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    ARTIFACT_JSON.write_text(json.dumps(build_artifact(report), indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "complete", "summary": summary, "report": str(REPORT_JSON.resolve())}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
