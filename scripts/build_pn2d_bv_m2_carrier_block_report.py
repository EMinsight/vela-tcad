#!/usr/bin/env python3
"""Build the canonical report artifact for the M2 carrier-block audit."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    output = args.output or args.root / "report_artifact.json"
    result = json.loads((args.root / "result.json").read_text(encoding="utf-8"))
    cases = rows(args.root / "case_summary.csv")
    modes = rows(args.root / "dominant_singular_modes.csv")
    columns = rows(args.root / "dominant_columns.csv")

    conditioning = [{
        "reverse_bias_V": abs(float(row["bias_V"])),
        "state": "baseline" if row["variant"] == "vela_baseline" else "joint QFP",
        "l2_equilibrated_condition": float(row["l2_equilibrated_condition_number"]),
    } for row in cases]
    counterfactual = []
    for row in cases:
        if abs(float(row["bias_V"])) < 19.5:
            continue
        state = "baseline" if row["variant"] == "vela_baseline" else "joint QFP"
        for operator, field in (
            ("no cross-carrier", "no_cross_carrier_relative_difference_from_full"),
            ("no avalanche", "no_avalanche_relative_difference_from_full"),
            ("no recombination", "no_recombination_relative_difference_from_full"),
        ):
            counterfactual.append({
                "reverse_bias_V": abs(float(row["bias_V"])),
                "state_bias": f"{state} {abs(float(row['bias_V'])):g} V",
                "operator": operator,
                "relative_step_change": float(row[field]),
            })
    mode_energy = []
    grouped_modes: dict[tuple[float, str], list[dict[str, str]]] = {}
    for row in modes:
        key = (float(row["bias_V"]), row["variant"])
        grouped_modes.setdefault(key, []).append(row)
    for (bias, variant), group in sorted(grouped_modes.items()):
        top_two = sorted(group, key=lambda row: int(row["step_energy_rank"]))[:2]
        mode_energy.append({
            "reverse_bias_V": abs(bias),
            "state": "baseline" if variant == "vela_baseline" else "joint QFP",
            "top_two_step_energy_share": sum(float(row["step_energy_fraction"]) for row in top_two),
            "minimum_relative_singular_value": min(float(row["relative_singular_value"]) for row in top_two),
        })
    case_table = [{
        "reverse_bias_V": abs(float(row["bias_V"])),
        "state": "baseline" if row["variant"] == "vela_baseline" else "joint QFP",
        "carrier_step_norm_V": float(row["full_physical_step_norm_V"]),
        "l2_condition": float(row["l2_equilibrated_condition_number"]),
        "phin_residual": float(row["phin_residual_norm"]),
        "phip_residual": float(row["phip_residual_norm"]),
        "no_cross_step_change": float(row["no_cross_carrier_relative_difference_from_full"]),
        "no_avalanche_step_change": float(row["no_avalanche_relative_difference_from_full"]),
    } for row in cases]

    knee_joint = [row for row in cases if row["variant"] == "sent_qfp_only" and abs(float(row["bias_V"])) >= 19.5]
    knee_base = [row for row in cases if row["variant"] == "vela_baseline" and abs(float(row["bias_V"])) >= 19.5]
    min_joint_condition = min(float(row["l2_equilibrated_condition_number"]) for row in knee_joint)
    max_joint_condition = max(float(row["l2_equilibrated_condition_number"]) for row in knee_joint)
    min_base_condition = min(float(row["l2_equilibrated_condition_number"]) for row in knee_base)
    max_base_condition = max(float(row["l2_equilibrated_condition_number"]) for row in knee_base)
    joint_top_two = [row["top_two_step_energy_share"] for row in mode_energy if row["state"] == "joint QFP" and row["reverse_bias_V"] >= 19.5]
    joint_update_columns = [row for row in columns if row["variant"] == "sent_qfp_only" and row["ranking"] == "update"]
    near_joint_updates = sum(float(row["distance_to_junction_um"]) <= 0.25 for row in joint_update_columns)

    base = "build-release/pn2d-bv-m2-carrier-block-decomposition-20260801"
    sources = [
        {"id": "result", "label": "Machine-readable verdict", "path": f"{base}/result.json"},
        {"id": "cases", "label": "Carrier-block case summary", "path": f"{base}/case_summary.csv",
         "query": {"engine": "duckdb", "language": "sql", "description": "Read conditions, residuals, and counterfactual solve differences.",
                   "sql": "SELECT * FROM read_csv_auto('case_summary.csv', header=true) ORDER BY abs(bias_V), variant",
                   "tables_used": ["case_summary.csv"], "filters": ["four frozen M2 biases", "baseline and joint-QFP states"]}},
        {"id": "modes", "label": "Dominant singular modes", "path": f"{base}/dominant_singular_modes.csv",
         "query": {"engine": "duckdb", "language": "sql", "description": "Read step-energy-ranked free carrier-block singular modes.",
                   "sql": "SELECT * FROM read_csv_auto('dominant_singular_modes.csv', header=true) ORDER BY abs(bias_V), variant, step_energy_rank",
                   "tables_used": ["dominant_singular_modes.csv"], "filters": ["top 12 modes per frozen case"]}},
        {"id": "columns", "label": "Dominant QFP columns and updates", "path": f"{base}/dominant_columns.csv"},
        {"id": "determinism", "label": "Independent-run hashes", "path": f"{base}/determinism.csv"},
    ]
    charts = [
        {"id": "condition_chart", "title": "L2-equilibrated free carrier-block condition number",
         "description": "Same M2 mesh and inputs; condition number is computed after row and column L2 equilibration on unconstrained QFP rows and columns.",
         "type": "line", "dataset": "conditioning", "sourceId": "cases",
         "encodings": {"x": {"field": "reverse_bias_V", "type": "quantitative", "title": "Reverse-bias magnitude (V)"},
                       "y": {"field": "l2_equilibrated_condition", "type": "quantitative", "title": "Condition number"},
                       "color": {"field": "state", "type": "nominal"}}},
        {"id": "mode_energy_chart", "title": "Newton-step energy captured by the two dominant singular modes",
         "description": "Projection of the production sparse-solver carrier step onto right singular vectors of the free carrier block.",
         "type": "bar", "dataset": "mode_energy", "sourceId": "modes",
         "encodings": {"x": {"field": "reverse_bias_V", "type": "nominal", "title": "Reverse-bias magnitude (V)"},
                       "y": {"field": "top_two_step_energy_share", "type": "quantitative", "title": "Step-energy share"},
                       "color": {"field": "state", "type": "nominal"}}},
        {"id": "counterfactual_chart", "title": "Carrier-step sensitivity to frozen operator removal",
         "description": "Relative change from the full carrier-only solution on the -19.5 to -20 V knee lattice; the residual is held fixed.",
         "type": "bar", "dataset": "counterfactual", "sourceId": "cases",
         "encodings": {"x": {"field": "operator", "type": "nominal", "title": "Counterfactual matrix"},
                       "y": {"field": "relative_step_change", "type": "quantitative", "title": "Relative step change"},
                       "color": {"field": "state_bias", "type": "nominal"}}},
    ]
    tables = [{
        "id": "case_table", "title": "Frozen-state carrier-block metrics",
        "description": "Exact residual, condition, step, and coupling sensitivity values.",
        "dataset": "case_table", "sourceId": "cases",
        "columns": [
            {"field": "reverse_bias_V", "label": "|V|", "format": "number"},
            {"field": "state", "label": "State", "format": "text"},
            {"field": "carrier_step_norm_V", "label": "Carrier step (V)", "format": "number"},
            {"field": "l2_condition", "label": "L2 condition", "format": "number"},
            {"field": "phin_residual", "label": "Electron residual", "format": "number"},
            {"field": "phip_residual", "label": "Hole residual", "format": "number"},
            {"field": "no_cross_step_change", "label": "No-cross change", "format": "percent"},
            {"field": "no_avalanche_step_change", "label": "No-avalanche change", "format": "percent"},
        ],
    }]
    blocks = [
        {"id": "title", "type": "markdown", "body": "# PN2D BV M2 carrier-block 线性解分解"},
        {"id": "summary", "type": "markdown", "sourceId": "result", "body": (
            "## 技术摘要\n\n"
            "冻结状态双运行通过，typed outcome 为 `carrier_block_linear_solve_decomposed`。"
            f"完整线性闭合误差最大 `{result['verdict']['maximum_full_linear_closure']:.3e}`，row-scaled 与 full 解最大相对差 `{result['verdict']['maximum_row_scaled_step_relative_difference']:.3e}`。"
            "未修改 SG/Laux、掺杂、网格、continuation、验收门限或生产默认值。"
        )},
        {"id": "finding1", "type": "markdown", "sourceId": "cases", "body": (
            "## 主要发现\n\n"
            f"在 -19.5 至 -20 V，joint-QFP 的 L2 平衡条件数为 `{min_joint_condition:.1f}`–`{max_joint_condition:.1f}`，baseline 为 `{min_base_condition:.1f}`–`{max_base_condition:.1f}`。"
            f"joint-QFP 的前两个超软模态承载 `{min(joint_top_two):.1%}`–`{max(joint_top_two):.1%}` 的实际 carrier Newton 步能量。"
        )},
        {"id": "condition", "type": "chart", "chartId": "condition_chart"},
        {"id": "modes", "type": "chart", "chartId": "mode_energy_chart"},
        {"id": "finding2", "type": "markdown", "sourceId": "modes", "body": (
            "## 耦合与空间支撑\n\n"
            "雪崩项贡献了近乎全部电子—空穴交叉块；复合项对 carrier 步的相对影响约为 1e-6。"
            "去掉交叉块会显著改变 joint-QFP 步长，但方向余弦仍高于 0.992，表明当前证据是幅值放大而不是更新方向翻转。"
            f"joint-QFP 的 `{near_joint_updates}/{len(joint_update_columns)}` 个 top-10 更新节点都位于结区 0.25 um 内。"
            "最强两个模态通常落在 x=0.75/1.25 um 的均匀掺杂结区肩部，并不直接接触 x=1.0 um 的补偿结区三角形。"
        )},
        {"id": "counterfactual", "type": "chart", "chartId": "counterfactual_chart"},
        {"id": "table", "type": "table", "tableId": "case_table"},
        {"id": "method", "type": "markdown", "body": (
            "## 范围、方法与限制\n\n"
            "SVD 与条件数仅使用去除接触约束后的 free carrier QFP 行列；六个反事实解使用同一冻结残差。"
            "raw 与 production row-scaled 矩阵因尺度跨度导致相对秩阈值只解析部分模态，因此跨状态判断采用满秩 L2 平衡条件数。"
            "本实验能定位线性反馈链的放大支撑，但不能单独证明 Sentaurus 与 Vela 的三角形内掺杂插值不同。"
        )},
        {"id": "next", "type": "markdown", "body": (
            "## 下一步\n\n"
            "保持生产物理不变，对 x=0.75–1.25 um 结区肩部执行只读的控制体积/三角形掺杂插值审计，并将两个超软模态投影到 transport 与 avalanche 分项 Jacobian。"
            "只有发现错误符号、缺失导数或与 Sentaurus 不一致的掺杂/控制体积语义后，才提出 opt-in 修正。"
        )},
    ]
    artifact = {
        "surface": "report",
        "manifest": {"version": 1, "surface": "report", "title": "M2 carrier-block 线性解分解",
                     "generatedAt": "2026-08-01T12:30:00+08:00", "cards": [], "charts": charts,
                     "tables": tables, "blocks": blocks, "sources": sources},
        "snapshot": {"version": 1, "generatedAt": "2026-08-01T12:30:00+08:00",
                     "status": "ready", "datasets": {
            "conditioning": conditioning,
            "mode_energy": mode_energy,
            "counterfactual": counterfactual,
            "case_table": case_table,
        }},
        "sources": sources,
    }
    output.write_text(json.dumps(artifact, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
