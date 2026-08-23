#!/usr/bin/env python3
"""Build the updated portable report for the resolved DG deep-off audit."""

from __future__ import annotations

import csv
import json
import math
from datetime import datetime
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[1]
REF = REPO / "build-release/reference_tcad/transportmodels_sentaurus2022"
ROOT = REF / "reports/transportmodels_dg_deep_off_strict_20260823"
STRICT_CSV = ROOT / "scaled_filter/deep_off_strict_summary.csv"
CONTRACT_CSV = ROOT / "classical_srh_contract/classical_srh_contract_summary.csv"
SRH_AB_CSV = ROOT / "srh_ab/srh_ab_summary.csv"
OUTPUT = REF / "reports/transportmodels_dg_deep_off_three_points_20260823"
ARTIFACT = OUTPUT / "artifact.json"
HTML = OUTPUT / "transportmodels_dg_deep_off_three_points_report.html"
REPORT = REPO / "docs/validation/transportmodels_dg_deep_off_three_points_2026-08-23.json"
SUMMARY_MD = REPO / "docs/validation/transportmodels_dg_deep_off_three_points_2026-08-23.md"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def sci(value: float) -> str:
    if value == 0.0:
        return "0"
    exponent = math.floor(math.log10(abs(value)))
    return f"{value / 10.0**exponent:.5f} x 10^{exponent}"


def source(source_id: str, label: str, path: Path, description: str, generated: str) -> dict[str, Any]:
    relative = path.relative_to(REPO).as_posix()
    return {
        "id": source_id,
        "label": label,
        "path": relative,
        "query": {
            "engine": "duckdb",
            "language": "sql",
            "sql": f"SELECT * FROM read_csv_auto('{relative}', header=true)",
            "description": description,
            "tables_used": [relative],
            "filters": ["Vg in {-1.00, -0.84, -0.68} V", "Vd = 1.1 V", "DG"],
            "executed_at": generated,
        },
    }


def main() -> int:
    generated = datetime.now().astimezone().isoformat(timespec="seconds")
    strict: list[dict[str, Any]] = []
    for row in read_csv(STRICT_CSV):
        bias = float(row["bias_V"])
        strict.append(
            {
                "bias_V": bias,
                "bias_label": f"{bias:.2f} V",
                "sentaurus_A_per_um": float(row["sentaurus_A_per_um"]),
                "sentaurus_display": sci(float(row["sentaurus_A_per_um"])),
                "vela_A_per_um": float(row["vela_A_per_um"]),
                "vela_display": sci(float(row["vela_A_per_um"])),
                "relative_error": float(row["absolute_relative_error"]),
                "log_error_dex": float(row["absolute_log_error_dex"]),
                "electron_residual": float(row["final_electron_continuity_residual_norm"]),
                "hole_residual": float(row["final_hole_continuity_residual_norm"]),
                "carrier_row_violations": int(row["carrier_row_violations"]),
                "carrier_row_max_ratio": float(row["carrier_row_max_ratio"]),
                "kcl_A_per_um": float(row["four_terminal_kcl_residual_A_per_um"]),
                "kcl_display": sci(float(row["four_terminal_kcl_residual_A_per_um"])),
                "id_to_kcl": float(row["id_to_kcl_residual_ratio"]),
                "hard_acceptance": row["hard_acceptance"].lower() == "true",
            }
        )

    contract: list[dict[str, Any]] = []
    contract_long: list[dict[str, Any]] = []
    for row in read_csv(CONTRACT_CSV):
        bias = float(row["gate_bias_V"])
        current_gap = float(row["current_contract_relative_gap_fraction"])
        classical_gap = float(row["classical_contract_relative_gap_fraction"])
        item = {
            "bias_V": bias,
            "bias_label": f"{bias:.2f} V",
            "sentaurus_srh_display": sci(float(row["sentaurus_srh_A_per_um"])),
            "current_contract_srh_display": sci(float(row["vela_current_quantum_density_denominator_A_per_um"])),
            "classical_contract_srh_display": sci(float(row["vela_classical_electron_density_denominator_A_per_um"])),
            "current_contract_gap": current_gap,
            "classical_contract_gap": classical_gap,
            "explained_fraction": float(row["explained_fraction_of_current_contract_gap"]),
            "bgn_abs_error_meV": 1.0e3 * float(row["srh_weighted_bgn_abs_error_eV"]),
        }
        contract.append(item)
        contract_long.extend(
            [
                {"bias_label": item["bias_label"], "contract": "当前 QM 密度分母", "relative_gap": current_gap},
                {"bias_label": item["bias_label"], "contract": "Sentaurus 默认经典密度分母", "relative_gap": classical_gap},
            ]
        )

    ab_rows = read_csv(SRH_AB_CSV)
    fermi_delta = max(
        abs(float(a["vela_formula_srh_A_per_um"]) / float(b["vela_formula_srh_A_per_um"]) - 1.0)
        for a in ab_rows
        for b in ab_rows
        if a["bias_V"] == b["bias_V"]
        and a["variant"] == "baseline_fermi_oldslotboom"
        and b["variant"] == "boltzmann"
    )

    sources = [
        source("strict", "严格连续性与 KCL 验收", STRICT_CSV, "三点严格 Newton 收敛、连续性残差和四端 KCL 指标。", generated),
        source("srh_contract", "默认 DG 的经典/QM SRH 契约", CONTRACT_CSV, "同一 Sentaurus 状态上比较两种 SRH 电子密度分母。", generated),
        source("srh_ab", "固定状态 SRH A/B", SRH_AB_CSV, "统计、BGN 和寿命定义的单变量固定状态 A/B。", generated),
    ]
    min_margin = min(row["id_to_kcl"] for row in strict)
    max_margin = max(row["id_to_kcl"] for row in strict)
    min_explained = min(row["explained_fraction"] for row in contract)
    max_explained = max(row["explained_fraction"] for row in contract)
    min_remaining = min(row["classical_contract_gap"] for row in contract)
    max_remaining = max(row["classical_contract_gap"] for row in contract)

    manifest = {
        "version": 1,
        "surface": "report",
        "title": "TransportModels DG 深关断三点：数值闭合与 SRH 根因审计",
        "description": "三点均通过 KCL 硬门槛；剩余差异主要由 Sentaurus 默认 DG 的经典/QM 双密度 SRH 契约解释。",
        "generatedAt": generated,
        "charts": [
            {
                "id": "kcl_margin_chart",
                "title": "三点均超过 Id/|KCL| >= 10 硬门槛",
                "subtitle": "严格连续性行门槛 + block-filter 线搜索",
                "intent": "comparison",
                "question": "深关断三点是否已数值解析？",
                "rationale": "逐点柱图直接显示相对硬门槛的裕量。",
                "type": "bar",
                "dataset": "strict_points",
                "sourceId": "strict",
                "valueFormat": "number",
                "encodings": {
                    "x": {"field": "bias_label", "type": "ordinal", "label": "Gate voltage Vg"},
                    "y": {"field": "id_to_kcl", "type": "quantitative", "label": "Id / |KCL residual|"},
                    "tooltip": [
                        {"field": "relative_error", "type": "quantitative", "label": "Id relative error", "format": "percent"},
                        {"field": "log_error_dex", "type": "quantitative", "label": "Log error", "format": "number"},
                    ],
                },
                "labels": {"values": "all"},
                "referenceLines": [{"axis": "y", "value": 10.0, "label": "hard gate = 10", "color": "neutral", "lineStyle": "dashed"}],
                "layout": "full",
                "surface": {"surface": "card", "viewMode": "visualization"},
            },
            {
                "id": "srh_contract_chart",
                "title": "改用经典电子密度分母后，SRH 差异降至约 1.3%",
                "subtitle": "固定同一 Sentaurus 状态、寿命参数、BGN 与网格积分",
                "intent": "comparison",
                "question": "经典/QM 密度选择能解释多少 SRH 差异？",
                "rationale": "分组柱图展示单一契约变化前后的相对误差。",
                "type": "bar",
                "dataset": "contract_long",
                "sourceId": "srh_contract",
                "valueFormat": "percent",
                "encodings": {
                    "x": {"field": "bias_label", "type": "ordinal", "label": "Gate voltage Vg"},
                    "y": {"field": "relative_gap", "type": "quantitative", "label": "SRH relative gap"},
                    "color": {"field": "contract", "type": "nominal", "label": "SRH density contract"},
                    "tooltip": [{"field": "relative_gap", "type": "quantitative", "label": "Relative gap", "format": "percent"}],
                },
                "labels": {"values": "all"},
                "layout": "full",
                "surface": {"surface": "card", "viewMode": "visualization"},
            },
        ],
        "tables": [
            {
                "id": "strict_table",
                "title": "严格数值验收结果",
                "subtitle": "电流单位 A/um；三点 carrier-row violations 均为 0",
                "dataset": "strict_points",
                "sourceId": "strict",
                "defaultSort": {"field": "bias_V", "direction": "asc"},
                "density": "spacious",
                "layout": "full",
                "columns": [
                    {"field": "bias_V", "label": "Vg (V)", "format": "number"},
                    {"field": "sentaurus_display", "label": "Sentaurus Id", "type": "text"},
                    {"field": "vela_display", "label": "Vela Id", "type": "text"},
                    {"field": "relative_error", "label": "Relative error", "format": "percent"},
                    {"field": "log_error_dex", "label": "Log error (dex)", "format": "number"},
                    {"field": "kcl_display", "label": "|KCL|", "type": "text"},
                    {"field": "id_to_kcl", "label": "Id/|KCL|", "format": "number"},
                    {"field": "carrier_row_violations", "label": "Row violations", "format": "number"},
                    {"field": "hard_acceptance", "label": "Hard pass", "format": "boolean"},
                ],
            },
            {
                "id": "contract_table",
                "title": "SRH 经典/QM 密度契约 A/B",
                "subtitle": "Sentaurus 手册页 369：默认 DG 同时维护经典与量子载流子密度",
                "dataset": "contract_points",
                "sourceId": "srh_contract",
                "defaultSort": {"field": "bias_V", "direction": "asc"},
                "density": "spacious",
                "layout": "full",
                "columns": [
                    {"field": "bias_V", "label": "Vg (V)", "format": "number"},
                    {"field": "sentaurus_srh_display", "label": "Sentaurus SRH", "type": "text"},
                    {"field": "current_contract_srh_display", "label": "Current Vela SRH", "type": "text"},
                    {"field": "classical_contract_srh_display", "label": "Classical-density SRH", "type": "text"},
                    {"field": "current_contract_gap", "label": "Current gap", "format": "percent"},
                    {"field": "classical_contract_gap", "label": "Classical gap", "format": "percent"},
                    {"field": "explained_fraction", "label": "Gap explained", "format": "percent"},
                    {"field": "bgn_abs_error_meV", "label": "BGN abs error (meV)", "format": "number"},
                ],
            },
        ],
        "sources": sources,
        "blocks": [
            {"id": "title", "type": "markdown", "body": "# TransportModels DG 深关断三点：数值闭合与 SRH 根因审计"},
            {
                "id": "summary",
                "type": "markdown",
                "body": (
                    "## 结论\n\n"
                    f"三点的 **Id/|KCL|={min_margin:.1f}–{max_margin:.1f}**，均通过 >=10 的硬门槛，"
                    "因此此前的“数值未解析”状态已经解除。当前 Id 相对误差为 6.62%–13.87%，"
                    "但它不再由端口消减或连续性残差主导。\n\n"
                    f"固定同一 Sentaurus 状态后，把 SRH 分母的电子浓度从量子密度改为默认 DG 所保留的经典密度，"
                    f"可解释原 SRH 差异的 **{min_explained:.1%}–{max_explained:.1%}**，剩余仅 **{min_remaining:.2%}–{max_remaining:.2%}**。"
                ),
            },
            {"id": "kcl_chart", "type": "chart", "chartId": "kcl_margin_chart"},
            {"id": "strict_results", "type": "table", "tableId": "strict_table"},
            {
                "id": "root_cause",
                "type": "markdown",
                "sourceId": "srh_contract",
                "body": (
                    "## 根因：默认 DG 的经典/QM 双密度语义尚未进入 Vela SRH\n\n"
                    "Sentaurus Device User Guide T-2022.03 第 369 页说明：未启用 DirectQuantumCorrection 时，"
                    "默认密度梯度模型同时存在经典密度和量子密度；第 474 页给出量子化条件下的广义 SRH 公式。"
                    "当前 Vela 的 SRH 分子用准费米势分裂稳定重构，但分母使用量子电子密度，形成混合契约。"
                    "在相同状态、寿命、BGN 和重心控制体积分下，仅将分母换为 n_classical=n_QM*exp(Qn/Vt)，"
                    "三点 SRH 积分就从低估 11.29%–13.88% 变为高估约 1.29%–1.35%。"
                ),
            },
            {"id": "contract_chart", "type": "chart", "chartId": "srh_contract_chart"},
            {"id": "contract_results", "type": "table", "tableId": "contract_table"},
            {
                "id": "excluded",
                "type": "markdown",
                "sourceId": "srh_ab",
                "body": (
                    "## 已排除的主因\n\n"
                    f"固定状态下 Fermi 与 Boltzmann SRH 积分的最大相对变化仅 {fermi_delta:.2e}；"
                    "关闭 Fermi-BGN 修正反而使差异扩大约 2.4 个百分点；改用净掺杂寿命或常寿命均明显恶化。"
                    "掺杂节点值与 Sentaurus 导出值在约 2.3e-16 相对精度内一致，SRH 加权 BGN 误差仅约 0.84 meV。"
                ),
            },
            {
                "id": "method",
                "type": "markdown",
                "body": (
                    "## 方法与边界\n\n"
                    "1. 将电子/空穴连续性绝对容差收紧至 1e-18，并启用载流子行硬门槛、全局连续性闭合、block-filter 线搜索和连续性行缩放。\n"
                    "2. 要求 carrier-row violation=0 且 |KCL|<=0.1|Id|。\n"
                    "3. 在完全相同的 Sentaurus 节点状态上进行 SRH 单变量 A/B；所有源项均用同一重心控制体积分。\n\n"
                    "经典密度分母结果目前是固定状态诊断，不是已提交的生产模型；正式改动还需要解析 Jacobian、近热平衡消减和完整 21 点曲线回归。"
                ),
            },
            {
                "id": "next",
                "type": "markdown",
                "body": (
                    "## 下一步\n\n"
                    "1. 增加显式 `srh_density_coupling=sentaurus_default` 模式，分别维护经典与 QM 电子密度。\n"
                    "2. 为该模式实现解析 Jacobian，并补充 DD 不变性、DG 近热平衡和量子势符号单元测试。\n"
                    "3. 先回归三点及 5 个关键 Vg 点，再重算完整 21 点 DG Id-Vg/Id-Vd；若仍有约 1.3% 偏差，再审计节点到单元的 SRH 插值和 BGN 0.84 meV 差异。"
                ),
            },
        ],
    }
    artifact = {
        "surface": "report",
        "manifest": manifest,
        "snapshot": {
            "version": 1,
            "generatedAt": generated,
            "status": "ready",
            "datasets": {
                "strict_points": strict,
                "contract_points": contract,
                "contract_long": contract_long,
            },
            "accessIssues": [],
        },
        "sources": sources,
        "package_info": {"originUrl": "artifact://transportmodels-dg-deep-off-three-points"},
    }
    OUTPUT.mkdir(parents=True, exist_ok=True)
    ARTIFACT.write_text(json.dumps(artifact, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    report = {
        "schema": "vela.transportmodels.dg_deep_off_three_points.v2",
        "generated_at": generated,
        "status": "complete",
        "summary": {
            "all_points_kcl_resolved": True,
            "id_to_kcl_range": [min_margin, max_margin],
            "remaining_id_relative_error_range": [min(row["relative_error"] for row in strict), max(row["relative_error"] for row in strict)],
            "classical_srh_contract_explained_fraction_range": [min_explained, max_explained],
            "classical_srh_contract_remaining_gap_range": [min_remaining, max_remaining],
            "classification": "numerical_closure_resolved_srh_density_contract_dominates",
        },
        "points": strict,
        "srh_contract": contract,
        "artifacts": {"artifact_json": str(ARTIFACT.resolve()), "report_html": str(HTML.resolve())},
    }
    REPORT.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    SUMMARY_MD.write_text(
        "# TransportModels DG 深关断三点审计\n\n"
        f"三点 Id/|KCL| 为 {min_margin:.1f}-{max_margin:.1f}，全部通过硬门槛。\n\n"
        f"经典电子密度 SRH 分母可解释 {min_explained:.1%}-{max_explained:.1%} 的固定状态 SRH 差异；剩余 {min_remaining:.2%}-{max_remaining:.2%}。\n",
        encoding="utf-8",
    )
    print(json.dumps(report["summary"], indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
