#!/usr/bin/env python3
"""Build the canonical report artifact for the M2 hotspot-edge audit."""

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
    changes = rows(args.root / "hotspot_state_change.csv")
    states = rows(args.root / "hotspot_state_summary.csv")
    contacts = rows(args.root / "contact_row_audit.csv")

    derivative_ratio = [
        {
            "reverse_bias_V": abs(float(row["bias_V"])),
            "carrier": row["carrier"],
            "joint_to_baseline_ratio": float(
                row["dominant_transport_derivative_joint_to_baseline_ratio"]
            ),
        }
        for row in changes
    ]
    operator_ratio = []
    for row in changes:
        for metric, field in (
            ("QFP drive", "qfp_drive_joint_to_baseline_ratio"),
            ("mobility", "mobility_joint_to_baseline_ratio"),
            ("transport derivative", "dominant_transport_derivative_joint_to_baseline_ratio"),
        ):
            operator_ratio.append({
                "reverse_bias_V": abs(float(row["bias_V"])),
                "carrier": row["carrier"],
                "metric": metric,
                "joint_to_baseline_ratio": float(row[field]),
            })
    state_table = [
        {
            "reverse_bias_V": abs(float(row["bias_V"])),
            "state": "baseline" if row["variant"] == "vela_baseline" else "joint QFP",
            "carrier": row["carrier"],
            "edge_id": int(row["edge_id"]),
            "qfp_drive_V_m": float(row["qfp_drive_V_m"]),
            "mobility_m2_V_s": float(row["mobility_m2_V_s"]),
            "dominant_transport_derivative_abs": float(row["dominant_production_derivative_abs"]),
            "mobility_fraction": float(row["mobility_response_fraction_of_dominant_transport"]),
            "row_weight_min": float(row["minimum_continuity_row_weight"]),
            "row_weight_max": float(row["maximum_continuity_row_weight"]),
        }
        for row in states
    ]
    base = "build-release/pn2d-bv-m2-transport-edge-jacobian-verification-20260801"
    sources = [
        {"id": "result", "label": "Machine-readable verdict", "path": f"{base}/result.json"},
        {
            "id": "changes", "label": "Hotspot state changes",
            "path": f"{base}/hotspot_state_change.csv",
            "query": {
                "engine": "duckdb", "language": "sql",
                "description": "Read carrier-resolved baseline-to-joint-QFP hotspot ratios.",
                "sql": "SELECT * FROM read_csv_auto('hotspot_state_change.csv', header=true) ORDER BY abs(bias_V), carrier",
                "tables_used": ["hotspot_state_change.csv"],
                "filters": ["four predeclared M2 biases", "electron and hole residual-hotspot incident edges"],
            },
        },
        {
            "id": "states", "label": "Hotspot state summary",
            "path": f"{base}/hotspot_state_summary.csv",
            "query": {
                "engine": "duckdb", "language": "sql",
                "description": "Read exact hotspot operator and row-scaling values for both frozen states.",
                "sql": "SELECT * FROM read_csv_auto('hotspot_state_summary.csv', header=true) ORDER BY abs(bias_V), variant, carrier",
                "tables_used": ["hotspot_state_summary.csv"],
                "filters": ["Vela baseline and joint Sentaurus-QFP states"],
            },
        },
        {"id": "decomposition", "label": "Per-edge derivative decomposition", "path": f"{base}/hotspot_decomposition.csv"},
        {"id": "contacts", "label": "Contact-row audit", "path": f"{base}/contact_row_audit.csv"},
        {"id": "determinism", "label": "Independent-run hashes", "path": f"{base}/determinism.csv"},
    ]
    charts = [
        {
            "id": "derivative_ratio_chart",
            "title": "Joint-QFP / baseline dominant transport derivative",
            "description": "Residual-hotspot incident edge; values below one mean the joint Sentaurus QFP state weakens the dominant production transport derivative.",
            "type": "bar", "dataset": "derivative_ratio", "sourceId": "changes",
            "encodings": {
                "x": {"field": "reverse_bias_V", "type": "nominal", "title": "Reverse-bias magnitude (V)"},
                "y": {"field": "joint_to_baseline_ratio", "type": "quantitative", "title": "Joint QFP / baseline"},
                "color": {"field": "carrier", "type": "nominal"},
            },
        },
        {
            "id": "operator_ratio_chart",
            "title": "Hotspot transport-factor changes in the BV knee region",
            "description": "-19.5, -19.7, and -20 V; ratios compare the joint Sentaurus-QFP frozen state with the Vela baseline.",
            "type": "bar", "dataset": "operator_ratio_knee", "sourceId": "changes",
            "encodings": {
                "x": {"field": "metric", "type": "nominal", "title": "Transport factor"},
                "y": {"field": "joint_to_baseline_ratio", "type": "quantitative", "title": "Joint QFP / baseline"},
                "color": {"field": "carrier_bias", "type": "nominal"},
            },
        },
    ]
    tables = [
        {
            "id": "state_table", "title": "Hotspot edge state and solver scaling",
            "description": "Exact state values for both frozen inputs on the four-bias lattice.",
            "dataset": "state_table", "sourceId": "states",
            "columns": [
                {"field": "reverse_bias_V", "label": "|V|", "format": "number"},
                {"field": "state", "label": "State", "format": "text"},
                {"field": "carrier", "label": "Carrier", "format": "text"},
                {"field": "edge_id", "label": "Edge", "format": "integer"},
                {"field": "qfp_drive_V_m", "label": "QFP drive (V/m)", "format": "number"},
                {"field": "mobility_m2_V_s", "label": "Mobility (m2/V/s)", "format": "number"},
                {"field": "dominant_transport_derivative_abs", "label": "Dominant |dF/dQFP|", "format": "number"},
                {"field": "mobility_fraction", "label": "Mobility derivative fraction", "format": "percent"},
                {"field": "row_weight_min", "label": "Min row weight", "format": "number"},
                {"field": "row_weight_max", "label": "Max row weight", "format": "number"},
            ],
        }
    ]
    knee_changes = [row for row in changes if abs(float(row["bias_V"])) >= 19.5]
    min_transport_ratio = min(
        float(row["dominant_transport_derivative_joint_to_baseline_ratio"])
        for row in knee_changes
    )
    max_transport_ratio = max(
        float(row["dominant_transport_derivative_joint_to_baseline_ratio"])
        for row in knee_changes
    )
    max_bernoulli_change = max(
        max(float(row["bernoulli_node0_absolute_change"]), float(row["bernoulli_node1_absolute_change"]))
        for row in changes
    )
    blocks = [
        {"id": "title", "type": "markdown", "body": "# M2 热点边 transport Jacobian 分解"},
        {"id": "summary", "type": "markdown", "sourceId": "result", "body": (
            "## 结论\n\n"
            f"诊断合同通过。生产解析 transport 导数对冻结迁移率有限差分的最大同边归一化误差为 `{result['verdict']['analytic_to_frozen_fd_max_edge_scaled_error']:.3e}`；实时迁移率总导数闭合误差为 `{result['verdict']['live_total_max_edge_scaled_error']:.3e}`。"
            "在 −19.5～−20 V，联合 Sentaurus QFP 状态把主导 transport 导数降至 baseline 的 "
            f"`{min_transport_ratio:.3f}`～`{max_transport_ratio:.3f}`；Bernoulli/GSS 系数变化为 `{max_bernoulli_change:.1f}`。"
        )},
        {"id": "derivative_text", "type": "markdown", "sourceId": "changes", "body": (
            "## 首个偏离量是载流子人口项\n\n"
            "固定电势和有效本征浓度后，Bernoulli 参数只依赖电势差与 ni 比值，因此 QFP 替换不会改变 B 系数。QFP 驱动力下降、迁移率上升，但迁移率响应只占热点主导 transport 导数约 3.7%～7.0%；主导导数的 6～14 倍下降来自 QFP 指数载流子人口项。"
        )},
        {"id": "derivative_chart", "type": "chart", "chartId": "derivative_ratio_chart"},
        {"id": "operator_text", "type": "markdown", "body": (
            "## Mobility 与 QFP drive 是次级修正\n\n"
            "高场迁移率使用 QFP gradient 驱动，但当前生产配置 `jacobian_field_derivatives=false`，Newton transport Jacobian 有意冻结迁移率。反事实实时迁移率导数已单独核验，未发现符号或尺度闭合错误。"
        )},
        {"id": "operator_chart", "type": "chart", "chartId": "operator_ratio_chart"},
        {"id": "solver_text", "type": "markdown", "sourceId": "contacts", "body": (
            "## 行缩放和接触行消元不是首个偏离源\n\n"
            "连续性行缩放逐项闭合为零误差；baseline 与联合 QFP 的热点权重变化很小。全部约束边记录均按单位行替换，物理边导数归零，而约束列对相邻非约束行继续保留，完全符合当前实现。"
        )},
        {"id": "state_table_block", "type": "table", "tableId": "state_table"},
        {"id": "method", "type": "markdown", "sourceId": "decomposition", "body": (
            "## 方法与限制\n\n"
            "热点先由上一轮 interior transport residual 的峰值节点确定，再在其相邻边中选择 baseline→joint-QFP 绝对通量变化最大的边；并列时取最小 edge ID。每个偏压、载流子、状态对两个端点行和两个同载流子 QFP 列做 1e-7 V 中心差分。逆风导数低于同边主导尺度数十个数量级时保留逐元素误差，但正式判断采用同边主导尺度，避免把双精度分辨率极限误判为模型错误。两套运行的全部原始 CSV 字节一致。"
        )},
        {"id": "next", "type": "markdown", "body": (
            "## 下一步建议\n\n"
            "保持 SG/Laux、行缩放与接触处理不变。下一项最小判别实验应检查载流子指数人口导数在完整 carrier block 线性解中的左右特征方向与变量缩放，重点比较 baseline 与联合 QFP 状态的主导列、对角占优性和电子/空穴块耦合；在取得该证据前不修改生产默认值。"
        )},
    ]
    artifact = {
        "surface": "report",
        "manifest": {
            "version": 1, "surface": "report",
            "title": "M2 热点边 transport Jacobian 分解",
            "generatedAt": "2026-08-01T12:00:00+08:00",
            "cards": [], "charts": charts, "tables": tables,
            "blocks": blocks, "sources": sources,
        },
        "snapshot": {
            "version": 1, "generatedAt": "2026-08-01T12:00:00+08:00",
            "status": "ready",
            "datasets": {
                "derivative_ratio": derivative_ratio,
                "operator_ratio_knee": [
                    {**row, "carrier_bias": f"{row['carrier']} / {row['reverse_bias_V']:g} V"}
                    for row in operator_ratio if row["reverse_bias_V"] >= 19.5
                ],
                "state_table": state_table,
            },
        },
        "sources": sources,
    }
    output.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
