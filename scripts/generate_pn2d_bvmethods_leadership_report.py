#!/usr/bin/env python3
"""Build the PN2D and BVmethods NMOS leadership progress report artifact."""

from __future__ import annotations

import argparse
import base64
import csv
import io
import json
import math
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import matplotlib.pyplot as plt


TITLE = "Vela TCAD 二维器件验证进展：PN2D 与 BVmethods NMOS"


plt.style.use("seaborn-v0_8-whitegrid")
plt.rcParams.update(
    {
        "font.family": "Microsoft YaHei",
        "axes.unicode_minus": False,
        "figure.facecolor": "white",
        "axes.facecolor": "#fbfdff",
    }
)


def png_data_uri(fig: plt.Figure) -> str:
    buffer = io.BytesIO()
    fig.savefig(buffer, format="png", dpi=150, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return "data:image/png;base64," + base64.b64encode(buffer.getvalue()).decode("ascii")


def chart_html(title: str, alt: str, caption: str, source_label: str, uri: str) -> str:
    return (
        "<figure>"
        f"<h3>{title}</h3>"
        f'<img src="{uri}" alt="{alt}" />'
        f"<figcaption>{caption}<br><small>来源：{source_label}</small></figcaption>"
        "</figure>"
    )


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def source(source_id: str, label: str, path: str, description: str) -> dict:
    if path.endswith(".csv"):
        sql = f"SELECT * FROM read_csv_auto('{path}', header=true)"
    elif path.endswith(".json"):
        sql = f"SELECT * FROM read_json_auto('{path}')"
    else:
        sql = f"SELECT '{path}' AS evidence_path"
    return {
        "id": source_id,
        "label": label,
        "path": path,
        "query": {
            "engine": "duckdb",
            "language": "sql",
            "description": description,
            "sql": sql,
            "tables_used": [path],
        },
    }


def table_spec(
    table_id: str,
    title: str,
    description: str,
    dataset: str,
    source_id: str,
    columns: list[tuple[str, str, str]],
    sort_field: str,
) -> dict:
    return {
        "id": table_id,
        "title": title,
        "description": description,
        "dataset": dataset,
        "sourceId": source_id,
        "columns": [
            {"field": field, "label": label, "format": fmt}
            for field, label, fmt in columns
        ],
        "defaultSort": {"field": sort_field, "direction": "asc"},
    }


def chart_spec(
    chart_id: str,
    title: str,
    description: str,
    chart_type: str,
    dataset: str,
    source_id: str,
    encodings: dict,
    question: str,
    rationale: str,
    max_rows: int = 500,
    reference_lines: list[dict] | None = None,
) -> dict:
    spec = {
        "id": chart_id,
        "title": title,
        "description": description,
        "type": chart_type,
        "dataset": dataset,
        "sourceId": source_id,
        "question": question,
        "rationale": rationale,
        "intent": "comparison" if chart_type in {"bar", "horizontalBar"} else "relationship",
        "encodings": encodings,
        "legend": {"position": "bottom", "interactive": True},
        "maxRows": max_rows,
    }
    if reference_lines:
        spec["referenceLines"] = reference_lines
    return spec


def build(repo: Path, artifact_relpath: str) -> dict:
    generated_at = datetime.now().astimezone().isoformat(timespec="seconds")

    pn_accept_path = repo / "build-release/pn2d-node-volume-atomic-default-acceptance-v2-20260801/acceptance.json"
    pn_accept = read_json(pn_accept_path)
    pn_mesh_path = repo / "build-release/pn2d-task10-balanced-junction-vm-v2-20260731/M2/converted/mesh.json"
    pn_mesh = read_json(pn_mesh_path)
    pn_doping_path = repo / "build-release/pn2d-task10-balanced-junction-vm-v2-20260731/M2/converted/doping.csv"
    pn_doping = {int(row["node_id"]): row for row in read_csv(pn_doping_path)}

    nmos_mesh_path = repo / "build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela/mesh.json"
    nmos_mesh = read_json(nmos_mesh_path)
    bv_summary_path = repo / "build-release/reference_tcad/bvmethods_sentaurus2018/run01/analysis/bvmethods_summary.csv"
    bv_summary_raw = read_csv(bv_summary_path)
    path_root = repo / "build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela_validation/sg_current_oriented_final_full_branch_20260805"

    pn_mesh_rows = []
    for node in pn_mesh["nodes"]:
        dop = pn_doping[node["id"]]
        net = float(dop["donors_cm3"]) - float(dop["acceptors_cm3"])
        side = "n型区" if net > 0 else "p型区" if net < 0 else "结区"
        pn_mesh_rows.append(
            {
                "node_id": node["id"],
                "x_um": node["x"],
                "y_um": node["y"],
                "doping_region": side,
            }
        )

    node_regions: dict[int, set[int]] = defaultdict(set)
    for triangle in nmos_mesh["triangles"]:
        for node_id in triangle["node_ids"]:
            node_regions[int(node_id)].add(int(triangle["region_id"]))
    region_material = {int(r["id"]): r["material"] for r in nmos_mesh["regions"]}
    nmos_mesh_rows = []
    for node in nmos_mesh["nodes"]:
        materials = {region_material[r] for r in node_regions[node["id"]]}
        # Interface nodes are colored by transport-material priority for display only.
        material = "Silicon" if "Si" in materials else "Nitride" if "Nitride" in materials else "SiO₂"
        nmos_mesh_rows.append(
            {
                "node_id": node["id"],
                "x_um": node["x"],
                "display_y_um": -node["y"],
                "material": material,
            }
        )
    if len(nmos_mesh_rows) > 2000:
        grouped: dict[str, list[dict]] = defaultdict(list)
        for row in nmos_mesh_rows:
            grouped[row["material"]].append(row)
        quotas = {"Silicon": 1400, "SiO₂": 400, "Nitride": 200}
        sampled: list[dict] = []
        sampled_ids: set[int] = set()
        for material, group in grouped.items():
            take = min(quotas.get(material, 0), len(group))
            for index in (int(i * len(group) / take) for i in range(take)):
                row = group[index]
                sampled.append(row)
                sampled_ids.add(row["node_id"])
        if len(sampled) < 2000:
            for row in nmos_mesh_rows:
                if row["node_id"] not in sampled_ids:
                    sampled.append(row)
                    sampled_ids.add(row["node_id"])
                    if len(sampled) == 2000:
                        break
        nmos_mesh_rows = sampled

    forward_rows = []
    for row in pn_accept["forward_iv"]["anchors"]["rows"]:
        forward_rows.append(
            {
                "bias_V": row["bias_V"],
                "relative_error_pct": 100.0 * row["mixed_voronoi_relative_error"],
                "sentaurus_mA_per_um": 1000.0 * row["sentaurus_A_per_um"],
                "vela_mA_per_um": 1000.0 * row["mixed_voronoi_A_per_um"],
            }
        )

    pn_error_rows = []
    pn_acceptance_rows = []
    off_rmse = {"M0": 0.0000888, "M2": 0.0000387}
    continuity = {"M0": 0.002200, "M2": 0.00000192}
    mesh_counts = {"M0": (27, 32), "M2": (115, 191)}
    for grid in ("M0", "M2"):
        metrics = pn_accept["grids"][grid]["avalanche_on_metrics"]
        pn_error_rows.extend(
            [
                {"metric": f"{grid} 全区RMSE", "error_dex": metrics["all_nonzero"]["rmse_dex"]},
                {"metric": f"{grid} 最大误差", "error_dex": metrics["all_nonzero"]["maximum_dex"]},
                {"metric": f"{grid} 拐点RMSE", "error_dex": metrics["knee"]["rmse_dex"]},
            ]
        )
        sent_bv = metrics["V_break_V"]["sentaurus"]
        vela_bv = metrics["V_break_V"]["vela"]
        pn_acceptance_rows.append(
            {
                "grid_order": 0 if grid == "M0" else 1,
                "grid": grid,
                "nodes_triangles": f"{mesh_counts[grid][0]} / {mesh_counts[grid][1]}",
                "off_rmse_dex": off_rmse[grid],
                "on_rmse_dex": metrics["all_nonzero"]["rmse_dex"],
                "max_error_dex": metrics["all_nonzero"]["maximum_dex"],
                "sentaurus_bv_V": abs(sent_bv),
                "vela_bv_V": abs(vela_bv),
                "bv_error_V": abs(sent_bv - vela_bv),
                "continuity_ratio": continuity[grid],
                "conclusion": "通过",
            }
        )

    low_bias_values = [
        (0.001, 8.716406e-11, 1.028661e-10),
        (0.002, 1.706412e-10, 2.015401e-10),
        (0.005, 4.000987e-10, 4.737424e-10),
        (0.010, 7.199892e-10, 8.563496e-10),
    ]
    nmos_low_bias_rows = []
    nmos_low_bias_table_rows = []
    for bias, vela, sent in low_bias_values:
        label = f"{bias * 1000:g} mV"
        nmos_low_bias_rows.extend(
            [
                {"bias": label, "simulator": "Sentaurus", "current_nA_per_um": sent * 1e9},
                {"bias": label, "simulator": "Vela", "current_nA_per_um": vela * 1e9},
            ]
        )
        nmos_low_bias_table_rows.append(
            {
                "bias_mV": bias * 1000,
                "sentaurus_nA_per_um": sent * 1e9,
                "vela_nA_per_um": vela * 1e9,
                "vela_sentaurus_ratio": vela / sent,
                "log_error_dex": abs(math.log10(vela / sent)),
            }
        )

    method_labels = {
        "ABA_poisson": "ABA（Poisson近似）",
        "ABA_coupled": "IIC（耦合）",
        "resistor": "外接电阻",
        "voltage2current": "电压转电流",
        "continuation": "连续法",
        "transient": "瞬态法",
    }
    bv_method_rows = []
    for order, row in enumerate(bv_summary_raw):
        bv_method_rows.append(
            {
                "method_order": order,
                "method": method_labels[row["method"]],
                "bv_V": float(row["bv_V"]),
                "criterion": row["criterion"],
                "sentaurus_node": int(row["node"]),
            }
        )

    path_rows = []
    for csv_path in path_root.rglob("path_ionization_integrals.csv"):
        for row in read_csv(csv_path):
            rank = int(row["path_rank"])
            if rank <= 3:
                path_rows.append(
                    {
                        "bias_V": float(row["bias_V"]),
                        "rank": f"第{rank}路径",
                        "rank_number": rank,
                        "seed_node": int(row["seed_node_id"]),
                        "mean_integral": float(row["mean_ionization_integral"]),
                        "electron_integral": float(row["electron_ionization_integral"]),
                        "hole_integral": float(row["hole_ionization_integral"]),
                    }
                )
    path_rows.sort(key=lambda row: (row["bias_V"], row["rank_number"]))

    endpoint_bias = max(row["bias_V"] for row in path_rows)
    endpoint_lookup = {
        row["rank_number"]: row for row in path_rows if abs(row["bias_V"] - endpoint_bias) < 1e-9
    }
    sent_endpoint = {1: 1.794564, 2: 1.544187, 3: 1.428772}
    endpoint_rows = []
    for rank in (1, 2, 3):
        row = endpoint_lookup[rank]
        sent_value = sent_endpoint[rank]
        endpoint_rows.append(
            {
                "rank": rank,
                "seed_node": row["seed_node"],
                "vela_mean": row["mean_integral"],
                "sentaurus_mean": sent_value,
                "relative_error_pct": 100.0 * (row["mean_integral"] / sent_value - 1.0),
            }
        )

    input_rows = [
        {"order": 1, "dimension": "验证职责", "pn2d": "基础方程、IV与雪崩击穿精度基准", "nmos": "复杂MOS结构与六种击穿提取方法基准"},
        {"order": 2, "dimension": "器件结构", "pn2d": "2.0 µm × 0.5 µm 突变PN结，结位置1.0 µm", "nmos": "栅长0.13 µm、2 nm栅氧、1.13 µm总宽、1.0 µm衬底"},
        {"order": 3, "dimension": "主要掺杂", "pn2d": "p/n两侧均为1×10¹⁷ cm⁻³", "nmos": "衬底1×10¹⁸；源漏峰值5×10²⁰；延伸区1×10¹⁹ cm⁻³"},
        {"order": 4, "dimension": "网格", "pn2d": "M0：27节点/32三角形；M2：115/191", "nmos": "2719节点、5210三角形、6区域、4电极"},
        {"order": 5, "dimension": "电极与扫描", "pn2d": "阳极/阴极；反向0至−20 V，另有正向IV守护", "nmos": "源/漏/栅/衬底；漏端自适应高压扫描"},
        {"order": 6, "dimension": "统计与输运", "pn2d": "Poisson + 电子/空穴连续性；SG离散", "nmos": "Fermi–Dirac密度、欧姆接触、广义Einstein/SG"},
        {"order": 7, "dimension": "关键物理模型", "pn2d": "SRH、OldSlotboom、掺杂/高场迁移率、Van Overstraeten雪崩", "nmos": "Fermi、SRH、OldSlotboom、掺杂/高场迁移率、Eparallel雪崩、E2隧穿"},
        {"order": 8, "dimension": "当前结论", "pn2d": "M0/M2均通过验收", "nmos": "基础输运与路径算法已打通，IIC/电阻/电压转电流待闭合"},
    ]

    stage_rows = [
        {"order": 1, "stage": "官方算例与网格导入", "status": "已完成", "evidence": "2719节点、5210三角形；1909个半导体节点坐标对齐", "remaining": "无"},
        {"order": 2, "stage": "Fermi–Dirac与欧姆接触", "status": "已完成", "evidence": "密度、接触平衡态及广义Einstein/SG已进入统一实现", "remaining": "继续用高压分支回归守护"},
        {"order": 3, "stage": "低偏压IV", "status": "基本对齐", "evidence": "1–10 mV电流比0.841–0.847，误差小于0.08 dex", "remaining": "消除约15%幅值偏差"},
        {"order": 4, "stage": "严格高压分支", "status": "已完成", "evidence": "7.0–10.4 V共36/36状态收敛，并达到10.448 V参考端点", "remaining": "按最终IIC定义重建闭合分支"},
        {"order": 5, "stage": "连续单元场线与路径积分", "status": "进行中", "evidence": "三条主路径端点误差3.7%–8.6%；路径回归测试通过", "remaining": "闭合第三物理路径出现电压与停止语义"},
        {"order": 6, "stage": "电流型IIC（目标6.377494 V）", "status": "待闭合", "evidence": "电场基本对齐，差异已收敛到SG电流空间支撑和α分布", "remaining": "使积分雪崩电流/漏电流在参考电压相等"},
        {"order": 7, "stage": "外接电阻（目标6.379792 V）", "status": "待完成", "evidence": "Sentaurus参考曲线和判据已提取", "remaining": "实现器件—电阻电路耦合并完成对标验收"},
        {"order": 8, "stage": "电压转电流（目标6.383184 V）", "status": "待完成", "evidence": "Sentaurus参考曲线和切换判据已提取", "remaining": "实现边界模式切换并完成对标验收"},
    ]

    root_cause_rows = [
        {"order": 1, "control": "源映射几何重构", "source_ratio": 1.000000, "reading": "数值舍入范围内闭合；不是当前主因"},
        {"order": 2, "control": "只替换为Sentaurus α", "source_ratio": 1.3733, "reading": "α空间分布可解释部分差异"},
        {"order": 3, "control": "只替换为Sentaurus边投影电流", "source_ratio": 1.02536, "reading": "单独看电流幅值接近，但空间支撑仍重要"},
        {"order": 4, "control": "同时替换Sentaurus α和矢量电流", "source_ratio": 2.11046, "reading": "接近直接Sentaurus产生率2.09583，仅差约0.7%"},
    ]

    sources = [
        source("pn2d_acceptance", "PN2D原子配置验收结果", "build-release/pn2d-node-volume-atomic-default-acceptance-v2-20260801/acceptance.json", "读取PN2D M0/M2击穿误差、正向IV锚点和验收状态。"),
        source("pn2d_acceptance_doc", "PN2D验收说明", "docs/validation/pn2d_node_volume_policy_atomic_default_acceptance_2026-08-01.md", "读取PN2D分支点数、连续性闭合和512项回归测试结论。"),
        source("pn2d_mesh", "PN2D M2网格", "build-release/pn2d-task10-balanced-junction-vm-v2-20260731/M2/converted/mesh.json", "展示M2的115个网格节点；p/n分区由同目录doping.csv计算。"),
        source("pn2d_deck", "PN2D Sentaurus输入脚本", "reference_tcad/pn2d_sentaurus2018_coarse7x3/source/pn2d_bv_sdevice.cmd", "核对偏压扫描、求解方程和物理模型。"),
        source("nmos_mesh", "BVmethods NMOS导入网格", "build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela/mesh.json", "展示2719节点、5210三角形、6区域和4电极。"),
        source("nmos_deck", "BVmethods NMOS Sentaurus器件脚本", "build-release/reference_tcad/bvmethods_sentaurus2018/run01/full_raw/pp4_des.cmd", "核对Fermi、SRH、迁移率、雪崩模型和高压求解流程。"),
        source("nmos_low_bias", "NMOS Fermi–Dirac实现验证记录", "docs/validation/bvmethods_nmos_fermi_dirac_implementation_2026-08-02.md", "读取1–10 mV低偏压Vela/Sentaurus电流对比。"),
        source("bvmethods_summary", "Sentaurus六种击穿方法汇总", "build-release/reference_tcad/bvmethods_sentaurus2018/run01/analysis/bvmethods_summary.csv", "读取六种官方方法的击穿电压和判据。"),
        source("path_followup", "NMOS路径离化积分跟踪记录", "docs/validation/bvmethods_nmos_path_ionization_followup_2026-08-04.md", "读取完整高压分支、三条主路径端点及源差异定位结果。"),
        source("report_synthesis", "本报告结构化数据", artifact_relpath, "将上述冻结证据重排为项目进展图表和里程碑表。"),
    ]

    charts = [
        chart_spec(
            "pn2d_mesh_chart", "PN2D M2网格与p/n区域", "115个节点；结附近网格加密，以捕捉高场与雪崩源。", "scatter", "pn2d_mesh", "pn2d_mesh",
            {"x": {"field": "x_um", "type": "quantitative", "title": "横向位置 x（µm）"}, "y": {"field": "y_um", "type": "quantitative", "title": "纵向位置 y（µm）"}, "color": {"field": "doping_region", "type": "nominal"}, "tooltip": [{"field": "node_id", "title": "节点"}, {"field": "doping_region", "title": "区域"}]},
            "基础PN结网格是否把分辨率放在关键结区？", "散点图直接显示节点分布和p/n区域。", 200,
        ),
        chart_spec(
            "pn2d_forward_error_chart", "PN2D正向IV锚点相对误差", "六个代表性偏压点；最大0.4066%，中位数0.2700%。", "bar", "pn2d_forward", "pn2d_acceptance",
            {"x": {"field": "bias_V", "type": "ordinal", "title": "正向偏压（V）"}, "y": {"field": "relative_error_pct", "type": "quantitative", "title": "相对误差（%）"}, "tooltip": [{"field": "sentaurus_mA_per_um", "title": "Sentaurus（mA/µm）"}, {"field": "vela_mA_per_um", "title": "Vela（mA/µm）"}]},
            "基础正向输运误差有多大？", "条形图适合比较离散偏压锚点。", 20,
        ),
        chart_spec(
            "pn2d_bv_error_chart", "PN2D雪崩电流曲线误差（dex）", "0.005 dex约等于1.16%的电流倍率差；M0/M2均低于该量级。", "bar", "pn2d_errors", "pn2d_acceptance",
            {"x": {"field": "metric", "type": "ordinal", "title": "网格与指标"}, "y": {"field": "error_dex", "type": "quantitative", "title": "log₁₀电流误差（dex）"}},
            "不同网格下的雪崩电流误差是否稳定且足够小？", "按网格和误差统计量比较验收结果。", 20,
        ),
        chart_spec(
            "nmos_mesh_chart", "BVmethods NMOS二维网格与材料", "全模型2719节点；为保证便携报告性能，图中按材料等距抽样显示2000个节点。", "scatter", "nmos_mesh", "nmos_mesh",
            {"x": {"field": "x_um", "type": "quantitative", "title": "横向位置 x（µm）"}, "y": {"field": "display_y_um", "type": "quantitative", "title": "显示高度（µm）"}, "color": {"field": "material", "type": "nominal"}, "tooltip": [{"field": "node_id", "title": "节点"}, {"field": "material", "title": "材料"}]},
            "复杂NMOS需要多大的几何和材料分辨率？", "材料着色的网格节点显示模型复杂度和局部加密。", 3000,
        ),
        chart_spec(
            "nmos_low_bias_chart", "NMOS低偏压漏电流：Vela与Sentaurus", "1–10 mV；Vela为Sentaurus的84.1%–84.7%，误差小于0.08 dex。", "bar", "nmos_low_bias", "nmos_low_bias",
            {"x": {"field": "bias", "type": "ordinal", "title": "漏端偏压"}, "y": {"field": "current_nA_per_um", "type": "quantitative", "title": "漏电流（nA/µm）"}, "color": {"field": "simulator", "type": "nominal"}},
            "复杂NMOS在进入雪崩前，基础输运是否已经同量级？", "分组条形图直接显示两套仿真器的电流幅值。", 20,
        ),
        chart_spec(
            "bvmethods_chart", "Sentaurus六种击穿提取方法", "除Poisson近似ABA外，五种自洽/电路方法集中在6.377–6.384 V。", "bar", "bv_methods", "bvmethods_summary",
            {"x": {"field": "method", "type": "ordinal", "title": "提取方法"}, "y": {"field": "bv_V", "type": "quantitative", "title": "击穿电压（V）"}, "tooltip": [{"field": "criterion", "title": "判据"}]},
            "官方参考中，不同工程提取方法给出的击穿电压是否一致？", "方法条形图突出主结果的聚集度和Poisson近似偏低。", 20,
        ),
        chart_spec(
            "path_integral_chart", "Vela高压分支前三条离化路径积分", "阈值1代表单条路径达到雪崩自持条件；当前算术均值定义下第三路径约在7.277 V过阈。", "line", "path_integrals", "path_followup",
            {"x": {"field": "bias_V", "type": "quantitative", "title": "漏端偏压（V）"}, "y": {"field": "mean_integral", "type": "quantitative", "title": "平均离化积分"}, "color": {"field": "rank", "type": "nominal"}, "tooltip": [{"field": "seed_node", "title": "种子节点"}]},
            "修正后的路径算法能否连续跟踪三条最危险雪崩通道？", "折线图展示路径随偏压的单调演化和阈值交叉。", 200,
            [{"axis": "y", "value": 1.0, "label": "路径阈值 = 1", "lineStyle": "dashed", "color": "neutral"}],
        ),
    ]

    tables = [
        table_spec("input_table", "两套算例的定位与输入参数", "同一套求解器先用简单PN结验公式，再用复杂NMOS验工程方法。", "inputs", "report_synthesis", [("order", "序号", "number"), ("dimension", "维度", "text"), ("pn2d", "PN2D", "text"), ("nmos", "BVmethods NMOS", "text")], "order"),
        table_spec("pn_acceptance_table", "PN2D验收结果", "BV按电流曲线、击穿电压和方程闭合三类指标共同验收。", "pn_acceptance", "pn2d_acceptance_doc", [("grid", "网格", "text"), ("nodes_triangles", "节点/三角形", "text"), ("off_rmse_dex", "非雪崩RMSE（dex）", "number"), ("on_rmse_dex", "雪崩RMSE（dex）", "number"), ("max_error_dex", "最大误差（dex）", "number"), ("sentaurus_bv_V", "Sentaurus |BV|（V）", "number"), ("vela_bv_V", "Vela |BV|（V）", "number"), ("bv_error_V", "BV误差（V）", "number"), ("continuity_ratio", "连续性闭合比", "number"), ("conclusion", "结论", "text")], "grid"),
        table_spec("nmos_low_bias_table", "NMOS低偏压数值明细", "该区间主要检查接触、载流子统计和SG输运，不受雪崩放大主导。", "nmos_low_bias_table", "nmos_low_bias", [("bias_mV", "偏压（mV）", "number"), ("sentaurus_nA_per_um", "Sentaurus（nA/µm）", "number"), ("vela_nA_per_um", "Vela（nA/µm）", "number"), ("vela_sentaurus_ratio", "Vela/Sentaurus", "number"), ("log_error_dex", "误差（dex）", "number")], "bias_mV"),
        table_spec("bv_method_table", "六种官方方法的精确参考值", "ABA Poisson是近似预估；其余方法用于工程击穿提取。", "bv_methods", "bvmethods_summary", [("method_order", "序号", "number"), ("method", "方法", "text"), ("bv_V", "击穿电压（V）", "number"), ("criterion", "提取判据", "text")], "method_order"),
        table_spec("path_endpoint_table", "10.448 V三条主路径端点对比", "Vela已复现三条主通道的积分量级，但路径出现顺序和停止判据仍需闭合。", "path_endpoint", "path_followup", [("rank", "路径排名", "number"), ("seed_node", "Vela种子节点", "number"), ("vela_mean", "Vela积分", "number"), ("sentaurus_mean", "Sentaurus积分", "number"), ("relative_error_pct", "相对差（%）", "number")], "rank"),
        table_spec("stage_table", "BVmethods NMOS功能与验收进度", "“已完成”表示有可重复证据；“待完成”表示尚不能宣称与Sentaurus一致。", "stages", "report_synthesis", [("order", "序号", "number"), ("stage", "工作包", "text"), ("status", "状态", "text"), ("evidence", "当前证据", "text"), ("remaining", "尚需完成", "text")], "order"),
        table_spec("root_cause_table", "6.4 V雪崩源差异的替换实验", "比值为替换后Vela源积分相对原Vela；组合替换与直接Sentaurus产生率仅差约0.7%。", "root_cause", "path_followup", [("order", "序号", "number"), ("control", "控制实验", "text"), ("source_ratio", "源积分倍率", "number"), ("reading", "说明", "text")], "order"),
    ]

    blocks = [
        {"id": "title", "type": "markdown", "body": f"# {TITLE}"},
        {"id": "executive_summary", "type": "markdown", "body": "## Executive Summary｜进展摘要\n\n- **基础能力已经闭环。** PN2D在两套网格上均通过验收：击穿电压与Sentaurus仅差**0.001 V**，雪崩电流曲线最大误差低于**0.0047 dex**，正向IV锚点最大相对误差**0.4066%**。\n- **复杂NMOS已跨过主要工程门槛。** 官方2719节点模型已经导入；Fermi–Dirac统计、欧姆接触、广义SG输运、E2隧穿和连续单元路径积分均已落地，高压分支已稳定计算到**10.448 V**。\n- **当前尚不能宣布BVmethods完全对齐。** 三条主雪崩路径在参考端点的积分误差为**3.7%–8.6%**，但第三条物理路径的出现电压和平均积分定义尚未与Sentaurus闭合；因此电流型IIC、外接电阻和电压转电流三项仍属于下一阶段。\n- **阶段性结论：**“简单器件精度已验收，复杂器件核心物理链路已打通，剩余工作聚焦于击穿判据和电路边界方法，而不是重新搭建求解器。”"},
        {"id": "background", "type": "markdown", "body": "## 1. 为什么要做两套算例\n\nTCAD可以理解为一座“虚拟电学实验室”：先把器件切成许多小网格，再在每个网格上联立求解电势、电子和空穴的守恒方程。**IV**检查加电压后流过多少电流，**BV**检查电场把载流子加速到碰撞电离、最终发生雪崩的临界电压。\n\nPN2D结构简单，适合把公式、离散、单位和雪崩源逐层验清；BVmethods NMOS包含栅氧、侧墙、源漏结和四个电极，适合验证真实器件中的高场通道与工程击穿提取方法。二者是“基础标尺 + 复杂实战”的递进关系。"},
        {"id": "input_table_block", "type": "table", "tableId": "input_table"},
        {"id": "solver", "type": "markdown", "body": "## 2. 求解方法与关键物理量\n\n每个偏压点同时求解三类方程：Poisson方程决定电势和电场，电子/空穴连续性方程保证电荷守恒，Scharfetter–Gummel（SG）格式稳定计算跨网格电流。高偏压时再计算碰撞电离系数 α，并沿最危险的载流子路径积分；积分达到1表示该路径具备自持雪崩条件。\n\n本轮重点比较六类关键量：端电流、击穿电压、电势/电场、电子/空穴密度、SG电流空间分布、α与路径离化积分。求解采用逐步升压和自适应步长；上一个偏压点的收敛状态作为下一个点的初值，以降低高压非线性求解失败率。"},
        {"id": "pn2d", "type": "markdown", "body": "## 3. PN2D：基础基准已经闭环\n\nM0用于快速回归，M2在PN结附近增加节点，用于检查网格加密后结论是否稳定。下图可见高场结区的网格更密；这正是雪崩产生率最敏感的位置。"},
        {"id": "pn_mesh", "type": "chart", "chartId": "pn2d_mesh_chart"},
        {"id": "pn_forward_text", "type": "markdown", "body": "正向IV检查普通导通输运。六个锚点的误差均低于0.41%，说明在不依赖雪崩放大的工作区，接触、迁移率和载流子输运已经与参考高度一致。"},
        {"id": "pn_forward", "type": "chart", "chartId": "pn2d_forward_error_chart"},
        {"id": "pn_bv_text", "type": "markdown", "body": "反向BV进一步放大任何局部误差。M0/M2的全区RMSE、最大误差和拐点RMSE仍保持在0.005 dex以下，击穿电压均只差0.001 V；这证明结果不是依赖单一粗网格偶然得到。"},
        {"id": "pn_bv", "type": "chart", "chartId": "pn2d_bv_error_chart"},
        {"id": "pn_table", "type": "table", "tableId": "pn_acceptance_table"},
        {"id": "pn_conclusion", "type": "markdown", "body": "**PN2D结论：** 两套网格、29/29反向偏压点、独立重复运行哈希和512/512项Release测试全部通过。该验收严格限定在合格的非钝角PN2D三角网格族，不外推为任意器件或全局默认。"},
        {"id": "nmos", "type": "markdown", "body": "## 4. BVmethods NMOS：复杂物理链路已打通\n\n该模型横向仅1.13 µm，却包含硅、2 nm栅氧和氮化硅侧墙等六个区域。源漏结、栅边缘和材料界面均需局部加密，因此网格规模约为PN2D M2的24倍。下图把器件纵向坐标翻转为常见剖面方向；共享界面节点仅为显示而按输运材料优先着色。"},
        {"id": "nmos_mesh", "type": "chart", "chartId": "nmos_mesh_chart"},
        {"id": "nmos_low_text", "type": "markdown", "body": "先在1–10 mV低偏压区验证基础输运，可避免雪崩模型掩盖接触或统计误差。Vela电流约为Sentaurus的84%，四个点的比例稳定，表明数量级和偏压响应已经一致，剩余约15%主要是幅值闭合问题。"},
        {"id": "nmos_low_chart_block", "type": "chart", "chartId": "nmos_low_bias_chart"},
        {"id": "nmos_low_table_block", "type": "table", "tableId": "nmos_low_bias_table"},
        {"id": "methods", "type": "markdown", "body": "## 5. 官方BVmethods参考：五种工程方法集中在约6.38 V\n\nSentaurus官方算例用六种方法从不同角度提取击穿。Poisson近似ABA忽略自洽载流子反馈，给出5.306 V的早期预估；耦合IIC、外接电阻、电压转电流、连续法和瞬态法集中在6.377–6.384 V，最大跨度仅约6.2 mV。因此后五种方法是Vela工程对标的主要目标。"},
        {"id": "methods_chart_block", "type": "chart", "chartId": "bvmethods_chart"},
        {"id": "methods_table_block", "type": "table", "tableId": "bv_method_table"},
        {"id": "paths", "type": "markdown", "body": "## 6. 路径离化积分：能连续跟踪，但击穿拓扑尚未闭合\n\n连续单元场线和best-vertex插值已经完成，7.0–10.448 V的36个高压状态全部收敛。前三条路径随偏压单调增强，说明算法已经能稳定跟踪危险通道，而不是在网格节点间跳变。"},
        {"id": "path_chart_block", "type": "chart", "chartId": "path_integral_chart"},
        {"id": "path_warning", "type": "markdown", "body": "图中第三路径按Vela当前“电子/空穴算术均值”在约**7.277 V**达到1，但这**不是最终闭合的击穿电压**。Sentaurus的第三条物理路径约到9.15–9.25 V才进入前三，且官方10.448 V停止点包含自适应步长越过阈值的影响；Sentaurus精确的 `I_mean` 代数定义也未在手册中公开。"},
        {"id": "path_table_block", "type": "table", "tableId": "path_endpoint_table"},
        {"id": "progress", "type": "markdown", "body": "## 7. 当前完成边界\n\n本项目已经从“能否运行NMOS”推进到“如何把击穿判据的最后几项语义与Sentaurus对齐”。以下状态表刻意区分功能存在、数值可比和正式验收，避免把局部进展误报为完整闭环。"},
        {"id": "stage_table_block", "type": "table", "tableId": "stage_table"},
        {"id": "root_cause", "type": "markdown", "body": "## 8. 剩余差异已经定位到哪里\n\n6.4 V控制实验表明，几何源映射可以精确重构；把Sentaurus的α分布和矢量电流同时代入Vela后，源积分与直接Sentaurus产生率仅差约0.7%。因此主要障碍已经收敛到**SG电流的空间支撑 + α的空间分布/路径语义**，而不是电场峰值、端电流提取器或几何权重。"},
        {"id": "root_cause_table_block", "type": "table", "tableId": "root_cause_table"},
        {"id": "recommendations", "type": "markdown", "body": "## 9. Recommended Next Steps｜建议下一步\n\n1. **先闭合路径型IIC。** 明确Sentaurus路径停止场、电子/空穴独立路径及 `I_mean` 排名语义，在冻结高压状态上复现第三路径出现电压。\n2. **再闭合电流型IIC。** 固定约6.38 V，逐边比较准费米势差、Bernoulli参数、广义Einstein因子、迁移率、几何权重和α，使积分雪崩电流与漏端传导电流相等。\n3. **实现外接电阻法。** 增加器件漏端电压与串联电阻压降的联立电路方程，以6.379792 V为正式验收目标。\n4. **实现电压转电流法。** 在高导数/负微分区从电压边界切换为电流边界，以6.383184 V为验收目标。\n5. **建立统一验收门槛。** 建议三种目标方法的BV绝对误差≤10 mV，同时要求分支连续、方程闭合和独立重复运行一致。"},
        {"id": "questions", "type": "markdown", "body": "## 10. Further Questions｜待确认事项\n\n- 下一阶段是否以“电流型IIC → 外接电阻 → 电压转电流”为唯一优先级，暂缓扩展到连续法和瞬态法？\n- 正式验收是否采用“BV误差≤10 mV + 曲线/方程闭合 + 可重复运行”三重门槛？\n- PN2D已具备独立技术评审条件，是否安排与NMOS后续开发并行的代码和科学评审？"},
        {"id": "caveats", "type": "markdown", "body": "## 11. Caveats and Assumptions｜口径与限制\n\n- Sentaurus作为本轮数值参考，不代表Vela必须复制其所有内部实现；验收看的是定义一致、守恒闭合和可重复结果。\n- PN2D结论只适用于已验收的M0/M2非钝角Tri3网格及其原子配置，不自动外推到NMOS或其他材料体系。\n- NMOS的7.277 V是Vela当前路径均值定义下的诊断交点，不应在汇报中称为器件BV。目标BV仍以官方IIC 6.377494 V、外接电阻6.379792 V和电压转电流6.383184 V为准。\n- 外接电阻和电压转电流当前有完整Sentaurus参考，但尚未完成Vela正式对标验收。\n- 便携打包器无法使用本机普通Chrome完成原生静态图表提取，因此报告按可视化技能的失败回退规则，把同一批已审阅数据生成的PNG嵌入HTML；原生图表定义和明细表仍保留在结构化artifact中。"},
    ]

    datasets = {
        "pn2d_mesh": pn_mesh_rows,
        "nmos_mesh": nmos_mesh_rows,
        "pn2d_forward": forward_rows,
        "pn2d_errors": pn_error_rows,
        "pn_acceptance": pn_acceptance_rows,
        "nmos_low_bias": nmos_low_bias_rows,
        "nmos_low_bias_table": nmos_low_bias_table_rows,
        "bv_methods": bv_method_rows,
        "path_integrals": path_rows,
        "path_endpoint": endpoint_rows,
        "inputs": input_rows,
        "stages": stage_rows,
        "root_cause": root_cause_rows,
    }

    static_charts: dict[str, str] = {}

    fig, ax = plt.subplots(figsize=(10, 3.8))
    for label, color in (("p型区", "#D64F4F"), ("n型区", "#2979B8"), ("结区", "#737B86")):
        group = [row for row in pn_mesh_rows if row["doping_region"] == label]
        if group:
            ax.scatter([r["x_um"] for r in group], [r["y_um"] for r in group], s=22, label=label, color=color)
    ax.axvline(1.0, color="#31363F", linestyle="--", linewidth=1, label="结位置")
    ax.set(xlabel="横向位置 x（µm）", ylabel="纵向位置 y（µm）", title="PN2D M2网格：结区局部加密")
    ax.set_aspect("equal", adjustable="box")
    ax.legend(ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.18), frameon=False)
    static_charts["pn2d_mesh_chart"] = chart_html(
        "PN2D M2网格与p/n区域",
        "PN2D M2网格散点图，p型区和n型区在1微米结位置相接",
        "115个节点；结附近的局部加密用于解析高场和雪崩源。",
        "PN2D M2网格与掺杂文件",
        png_data_uri(fig),
    )

    fig, ax = plt.subplots(figsize=(9.5, 4.2))
    x = list(range(len(forward_rows)))
    values = [row["relative_error_pct"] for row in forward_rows]
    bars = ax.bar(x, values, color="#2E73B8", width=0.62)
    ax.set_xticks(x, [f'{row["bias_V"]} V' for row in forward_rows])
    ax.set(ylabel="相对误差（%）", title="PN2D正向IV误差保持在0.41%以内")
    ax.set_ylim(0, max(values) * 1.28)
    ax.bar_label(bars, fmt="%.3f%%", padding=3)
    static_charts["pn2d_forward_error_chart"] = chart_html(
        "PN2D正向IV锚点相对误差",
        "六个正向偏压锚点的Vela相对Sentaurus电流误差条形图",
        "最大0.4066%，中位数0.2700%；普通导通输运已高度对齐。",
        "PN2D原子配置验收结果",
        png_data_uri(fig),
    )

    fig, ax = plt.subplots(figsize=(10, 4.2))
    labels = [row["metric"] for row in pn_error_rows]
    values = [row["error_dex"] for row in pn_error_rows]
    colors = ["#2E73B8" if label.startswith("M0") else "#F28E2B" for label in labels]
    bars = ax.bar(range(len(values)), values, color=colors, width=0.66)
    ax.set_xticks(range(len(values)), labels, rotation=18, ha="right")
    ax.set(ylabel=r"$\log_{10}$电流误差（dex）", title="两套网格的雪崩电流误差均低于0.005 dex")
    ax.axhline(0.005, color="#69717A", linestyle="--", linewidth=1, label="0.005 dex量级线")
    ax.set_ylim(0, 0.0058)
    ax.bar_label(bars, fmt="%.4f", padding=3, fontsize=8)
    ax.legend(frameon=False)
    static_charts["pn2d_bv_error_chart"] = chart_html(
        "PN2D雪崩电流曲线误差",
        "M0和M2网格的全区均方根误差、最大误差和拐点误差比较",
        "M0/M2结论一致，说明精度不依赖单一网格。",
        "PN2D原子配置验收结果",
        png_data_uri(fig),
    )

    fig, ax = plt.subplots(figsize=(10, 5.6))
    material_style = {
        "Silicon": ("#73A7D5", 5.0, 0.62),
        "SiO₂": ("#F3A64A", 8.0, 0.82),
        "Nitride": ("#8E6BBE", 8.0, 0.82),
    }
    for material in ("Silicon", "SiO₂", "Nitride"):
        group = [row for row in nmos_mesh_rows if row["material"] == material]
        color, size, alpha = material_style[material]
        ax.scatter(
            [r["x_um"] for r in group],
            [r["display_y_um"] for r in group],
            s=size,
            alpha=alpha,
            label="SiO$_2$" if material == "SiO₂" else material,
            color=color,
            linewidths=0,
        )
    ax.set(xlabel="横向位置 x（µm）", ylabel="显示高度（µm）", title="BVmethods NMOS：多材料与局部加密网格")
    ax.set_aspect("equal", adjustable="box")
    ax.legend(ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.12), frameon=False)
    static_charts["nmos_mesh_chart"] = chart_html(
        "BVmethods NMOS二维网格与材料",
        "NMOS网格散点剖面图，硅衬底上方包含栅氧和氮化硅侧墙",
        "全模型2719节点；图中按材料等距抽样2000个节点以控制便携HTML体积。",
        "BVmethods NMOS导入网格",
        png_data_uri(fig),
    )

    fig, ax = plt.subplots(figsize=(9.5, 4.4))
    labels = [f"{bias * 1000:g} mV" for bias, _, _ in low_bias_values]
    sent_values = [sent * 1e9 for _, _, sent in low_bias_values]
    vela_values = [vela * 1e9 for _, vela, _ in low_bias_values]
    x = list(range(len(labels)))
    width = 0.36
    sent_bars = ax.bar([v - width / 2 for v in x], sent_values, width, label="Sentaurus", color="#5B5F97")
    vela_bars = ax.bar([v + width / 2 for v in x], vela_values, width, label="Vela", color="#2CA58D")
    ax.set_xticks(x, labels)
    ax.set(ylabel="漏电流（nA/µm）", title="NMOS低偏压漏电流已对齐到稳定比例")
    ax.legend(frameon=False)
    ax.bar_label(sent_bars, fmt="%.3f", padding=3, fontsize=8)
    ax.bar_label(vela_bars, fmt="%.3f", padding=3, fontsize=8)
    ax.set_ylim(0, max(sent_values) * 1.25)
    static_charts["nmos_low_bias_chart"] = chart_html(
        "NMOS低偏压漏电流：Vela与Sentaurus",
        "1到10毫伏四个偏压点的Vela和Sentaurus漏电流分组条形图",
        "Vela为Sentaurus的84.1%–84.7%，误差小于0.08 dex；剩余是幅值闭合问题。",
        "NMOS Fermi–Dirac实现验证记录",
        png_data_uri(fig),
    )

    fig, ax = plt.subplots(figsize=(10, 4.8))
    ordered_methods = sorted(bv_method_rows, key=lambda row: row["method_order"], reverse=True)
    y = list(range(len(ordered_methods)))
    values = [row["bv_V"] for row in ordered_methods]
    colors = ["#B3B8BF" if "Poisson" in row["method"] else "#2E73B8" for row in ordered_methods]
    bars = ax.barh(y, values, color=colors, height=0.64)
    ax.set_yticks(y, [row["method"] for row in ordered_methods])
    ax.set(xlabel="击穿电压（V）", title="除Poisson近似外，五种工程方法集中在约6.38 V")
    ax.set_xlim(5.1, 6.55)
    for bar, value in zip(bars, values):
        ax.text(value + 0.015, bar.get_y() + bar.get_height() / 2, f"{value:.6f}", va="center", fontsize=9)
    static_charts["bvmethods_chart"] = chart_html(
        "Sentaurus六种击穿提取方法",
        "六种Sentaurus击穿提取方法的击穿电压水平条形图",
        "后五种方法跨度约6.2 mV；Poisson近似ABA用于早期预估，不能与自洽方法等同。",
        "Sentaurus六种击穿方法汇总",
        png_data_uri(fig),
    )

    fig, ax = plt.subplots(figsize=(10, 4.8))
    path_colors = {1: "#2E73B8", 2: "#F28E2B", 3: "#2CA58D"}
    for rank in (1, 2, 3):
        group = [row for row in path_rows if row["rank_number"] == rank]
        ax.plot(
            [row["bias_V"] for row in group],
            [row["mean_integral"] for row in group],
            marker="o",
            markersize=3.5,
            linewidth=2,
            label=f"第{rank}路径",
            color=path_colors[rank],
        )
    ax.axhline(1.0, color="#4B5158", linestyle="--", linewidth=1.2, label="路径阈值=1")
    ax.axvline(7.277244958, color="#B44B4B", linestyle=":", linewidth=1.2, label="当前rank-3诊断交点")
    ax.set(xlabel="漏端偏压（V）", ylabel="平均离化积分", title="Vela高压分支可连续跟踪前三条雪崩路径")
    ax.legend(ncol=3, frameon=False, loc="upper left")
    static_charts["path_integral_chart"] = chart_html(
        "Vela高压分支前三条离化路径积分",
        "7到10.448伏范围前三条路径平均离化积分随偏压变化的折线图",
        "三条路径均单调增强；7.277 V仅是当前算术均值定义的诊断交点，不是最终器件BV。",
        "NMOS路径离化积分跟踪记录",
        png_data_uri(fig),
    )

    for block in blocks:
        chart_id = block.get("chartId")
        if block.get("type") == "chart" and chart_id in static_charts:
            block.clear()
            block.update(
                {
                    "id": f"static_{chart_id}",
                    "type": "html",
                    "body": static_charts[chart_id],
                }
            )
    for block in blocks:
        if block.get("id") == "methods_table_block":
            block.clear()
            block.update(
                {
                    "id": "native_bvmethods_audit",
                    "type": "chart",
                    "chartId": "bvmethods_chart",
                }
            )
    return {
        "surface": "report",
        "manifest": {
            "version": 1,
            "surface": "report",
            "title": TITLE,
            "description": "PN2D与BVmethods NMOS验证进展、证据、风险和下一步。",
            "generatedAt": generated_at,
            "cards": [],
            "charts": charts,
            "tables": tables,
            "blocks": blocks,
            "sources": sources,
        },
        "snapshot": {"version": 1, "generatedAt": generated_at, "status": "ready", "datasets": datasets},
        "sources": sources,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    repo = args.repo.resolve()
    output = args.output.resolve()
    artifact_relpath = output.relative_to(repo).as_posix()
    artifact = build(repo, artifact_relpath)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(artifact, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
