#!/usr/bin/env python3
"""Generate compact daily-report figures for the Sentaurus BVmethods NMOS audit."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.font_manager import FontProperties, fontManager
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np


BLUE = "#2F6B9A"
GOLD = "#D39B2A"
ORANGE = "#D66A2C"
CHARCOAL = "#263238"
MID_GREY = "#7A858B"
LIGHT_GREY = "#D9E0E4"


def choose_font() -> FontProperties:
    available = {font.name for font in fontManager.ttflist}
    for name in ("Microsoft YaHei", "SimHei", "Noto Sans CJK SC", "Arial Unicode MS"):
        if name in available:
            plt.rcParams["font.sans-serif"] = [name, "DejaVu Sans"]
            plt.rcParams["axes.unicode_minus"] = False
            return FontProperties(family=name)
    return FontProperties(family="DejaVu Sans")


def configure_style() -> None:
    choose_font()
    plt.rcParams.update(
        {
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.edgecolor": CHARCOAL,
            "axes.labelcolor": CHARCOAL,
            "axes.titlecolor": CHARCOAL,
            "xtick.color": CHARCOAL,
            "ytick.color": CHARCOAL,
            "text.color": CHARCOAL,
            "grid.color": LIGHT_GREY,
            "grid.linewidth": 0.7,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.frameon": False,
            "savefig.facecolor": "white",
        }
    )


def repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def draw_mesh(mesh_path: Path, output_path: Path) -> None:
    mesh = json.loads(mesh_path.read_text(encoding="utf-8"))
    node_by_id = {int(node["id"]): (float(node["x"]), float(node["y"])) for node in mesh["nodes"]}
    triangles = sorted(mesh["triangles"], key=lambda item: int(item["id"]))
    vertices = [[node_by_id[int(node_id)] for node_id in tri["node_ids"]] for tri in triangles]

    region_by_id = {int(region["id"]): region for region in mesh["regions"]}
    region_colors = {
        0: "#BFD7EA",
        1: "#D8E6F1",
        2: "#D8E6F1",
        3: "#F1D59B",
        4: "#C8D3A5",
        5: "#C8D3A5",
    }
    face_colors = [region_colors[int(tri["region_id"])] for tri in triangles]

    fig, axes = plt.subplots(1, 2, figsize=(15.6, 7.8), gridspec_kw={"width_ratios": [1.0, 1.18]})

    for ax, zoom in zip(axes, (False, True)):
        collection = PolyCollection(
            vertices,
            facecolors=face_colors,
            edgecolors="#66757D",
            linewidths=0.08 if not zoom else 0.28,
            alpha=0.96,
        )
        ax.add_collection(collection)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel("横向位置 x (µm)")
        ax.set_ylabel("纵向位置 y (µm)")
        ax.grid(False)
        if zoom:
            ax.set_xlim(-0.27, 0.27)
            ax.set_ylim(0.12, -0.215)
            ax.set_title("栅氧—沟道局部网格", loc="left", fontsize=13, weight="bold")
        else:
            ax.set_xlim(-0.60, 0.60)
            ax.set_ylim(1.04, -0.23)
            ax.set_title("完整二维 NMOS 网格", loc="left", fontsize=13, weight="bold")

    contact_styles = {
        "gate": (BLUE, "Gate"),
        "source": (CHARCOAL, "Source"),
        "drain": (ORANGE, "Drain"),
        "substrate": (GOLD, "Substrate"),
    }
    for ax in axes:
        for contact in mesh["contacts"]:
            name = contact["name"].lower()
            color, _ = contact_styles[name]
            points = np.array([node_by_id[int(node_id)] for node_id in contact["node_ids"]])
            ax.scatter(points[:, 0], points[:, 1], s=13 if ax is axes[0] else 20, color=color, zorder=5)

    axes[0].annotate("Gate", xy=(0.0, -0.20), xytext=(0.19, -0.16), arrowprops={"arrowstyle": "->", "color": BLUE}, color=BLUE, weight="bold")
    axes[0].annotate("Source", xy=(-0.39, 0.0), xytext=(-0.51, 0.18), arrowprops={"arrowstyle": "->", "color": CHARCOAL}, weight="bold")
    axes[0].annotate("Drain", xy=(0.39, 0.0), xytext=(0.39, 0.18), arrowprops={"arrowstyle": "->", "color": ORANGE}, color=ORANGE, weight="bold")
    axes[0].annotate("Substrate", xy=(0.0, 1.0), xytext=(0.18, 0.89), arrowprops={"arrowstyle": "->", "color": GOLD}, color="#8A661B", weight="bold")

    axes[1].axhline(0.0, color=CHARCOAL, linewidth=1.2, linestyle="--", alpha=0.8)
    axes[1].annotate(
        "Si/绝缘层共享节点\nni_eff 采用 Si 输运材料",
        xy=(0.070, 0.0),
        xytext=(0.145, 0.095),
        ha="center",
        va="center",
        fontsize=10.5,
        arrowprops={"arrowstyle": "->", "color": CHARCOAL},
    )
    axes[1].text(0.0, -0.105, "Gate / SiO2", ha="center", va="center", color=BLUE, weight="bold")
    axes[1].text(0.0, 0.055, "Silicon channel", ha="center", va="center", color="#7A5813", weight="bold")

    legend_handles = [
        Patch(facecolor=region_colors[3], edgecolor="#66757D", label="Silicon"),
        Patch(facecolor=region_colors[0], edgecolor="#66757D", label="SiO2"),
        Patch(facecolor=region_colors[4], edgecolor="#66757D", label="Nitride spacer"),
        Line2D([], [], marker="o", linestyle="", color=BLUE, label="Gate contact"),
        Line2D([], [], marker="o", linestyle="", color=ORANGE, label="Drain contact"),
    ]
    fig.legend(handles=legend_handles, loc="lower center", ncol=5, bbox_to_anchor=(0.5, 0.01), fontsize=10)
    fig.suptitle("Sentaurus BVmethods 二维 NMOS 器件与网格", x=0.055, y=0.99, ha="left", fontsize=18, weight="bold")
    fig.text(0.055, 0.925, f"实际导入网格：{len(mesh['nodes'])} 个节点、{len(mesh['triangles'])} 个三角形、{len(region_by_id)} 个区域", fontsize=10.5, color=MID_GREY)
    fig.subplots_adjust(left=0.065, right=0.985, top=0.855, bottom=0.11, wspace=0.19)
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def read_sentaurus_iic(path: Path) -> dict[str, np.ndarray]:
    rows = []
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            bias = float(row["inner_voltage_V"])
            if bias >= 1.0:
                rows.append(
                    (
                        bias,
                        abs(float(row["drain_total_current_A_per_um"])),
                        abs(float(row["avalanche_current_A_per_um"])),
                    )
                )
    rows.sort()
    return {
        "bias": np.array([row[0] for row in rows]),
        "id": np.array([row[1] for row in rows]),
        "iava": np.array([row[2] for row in rows]),
    }


def write_metric_source(output_path: Path) -> None:
    rows = [
        ("low_bias_iv", "0.001 V", "drain_current_A_per_um", 8.716406e-11, 1.028661e-10, 0.847355, "strict Fermi-Dirac branch"),
        ("low_bias_iv", "0.002 V", "drain_current_A_per_um", 1.706412e-10, 2.015401e-10, 0.846686, "strict Fermi-Dirac branch"),
        ("low_bias_iv", "0.005 V", "drain_current_A_per_um", 4.000987e-10, 4.737424e-10, 0.844549, "strict Fermi-Dirac branch"),
        ("low_bias_iv", "0.010 V", "drain_current_A_per_um", 7.199892e-10, 8.563496e-10, 0.840765, "strict Fermi-Dirac branch"),
        ("peak_ratio_6p38V", "6.38 V", "electric_field", 2.445959e8, 2.279055e8, 1.07323, "old non-closed branch diagnostic"),
        ("peak_ratio_6p38V", "6.38 V", "electron_alpha", 4.599973e7, 4.400812e7, 1.04526, "old non-closed branch diagnostic"),
        ("peak_ratio_6p38V", "6.38 V", "hole_alpha", 6.663211e7, 3.498695e7, 1.90448, "old non-closed branch diagnostic"),
        ("peak_ratio_6p38V", "6.38 V", "avalanche_generation", 1.5767e32, 5.5783e36, 2.8265e-5, "corrected-unit old branch diagnostic"),
    ]
    with output_path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle)
        writer.writerow(["panel", "bias", "quantity", "vela", "sentaurus", "vela_over_sentaurus", "scope"])
        writer.writerows(rows)


def draw_physics(iic_path: Path, output_path: Path) -> None:
    # Valid low-bias results from the 2026-08-02 Fermi-Dirac implementation report.
    bias_mv = np.array([1.0, 2.0, 5.0, 10.0])
    vela_low_nA = np.array([8.716406e-11, 1.706412e-10, 4.000987e-10, 7.199892e-10]) * 1e9
    sentaurus_low_nA = np.array([1.028661e-10, 2.015401e-10, 4.737424e-10, 8.563496e-10]) * 1e9

    # Strict Vela postprocess-only branch rebuilt on 2026-08-03.
    vela_bias = np.array([1.0, 2.0, 4.0, 5.0, 5.725])
    vela_id_uA = np.array([2.3782694018477446e-9, 2.9093203249673174e-9, 3.9857869877108275e-9, 4.5900629993873126e-9, 4.9964043447318881e-9]) * 1e6

    iic = read_sentaurus_iic(iic_path)
    ratio_labels = ["电场峰值", "电子 α 峰值", "空穴 α 峰值", "雪崩生成率峰值"]
    ratios = np.array([1.07323, 1.04526, 1.90448, 2.8265e-5])

    fig = plt.figure(figsize=(16.0, 10.0))
    grid = fig.add_gridspec(2, 2, height_ratios=[1.02, 0.82], hspace=0.42, wspace=0.28)
    ax_low = fig.add_subplot(grid[0, 0])
    ax_iic = fig.add_subplot(grid[0, 1])
    ax_ratio = fig.add_subplot(grid[1, :])

    ax_low.plot(bias_mv, sentaurus_low_nA, color=CHARCOAL, linewidth=2.2, marker="s", markersize=6, label="Sentaurus")
    ax_low.plot(bias_mv, vela_low_nA, color=BLUE, linewidth=2.2, marker="o", markersize=6, label="Vela Fermi–Dirac")
    ax_low.fill_between(bias_mv, vela_low_nA, sentaurus_low_nA, color=BLUE, alpha=0.10)
    ax_low.set_xlim(0, 10.5)
    ax_low.set_ylim(0, 0.95)
    ax_low.set_xlabel("漏极偏压 Vd (mV)")
    ax_low.set_ylabel("漏极电流 Id (nA/µm)")
    ax_low.set_title("低偏压 IV 对比", loc="left", fontsize=13, weight="bold")
    ax_low.text(0.03, 0.92, "Vela/Sentaurus = 0.841–0.847\n误差 < 0.08 dex", transform=ax_low.transAxes, va="top", fontsize=10.5, color=BLUE, weight="bold")
    ax_low.legend(loc="lower right")
    ax_low.grid(axis="y")

    ax_iic.plot(iic["bias"], iic["id"] * 1e6, color=CHARCOAL, linewidth=2.2, marker="o", markersize=4.2, label="Sentaurus Id")
    ax_iic.plot(iic["bias"], iic["iava"] * 1e6, color=ORANGE, linewidth=2.2, marker="s", markersize=4.2, label="Sentaurus Iava")
    ax_iic.scatter(vela_bias, vela_id_uA, color=BLUE, marker="D", s=42, facecolors="white", linewidths=1.6, label="Vela DD（仅后处理）", zorder=5)
    ax_iic.axvline(6.377494278, color=BLUE, linestyle="--", linewidth=1.4, alpha=0.85)
    ax_iic.axvline(6.734425890, color=GOLD, linestyle=":", linewidth=2.0, alpha=0.95)
    ax_iic.text(6.31, 0.006, "官方稀疏提取\n6.377 V", rotation=90, ha="right", va="bottom", fontsize=8.8, color=BLUE)
    ax_iic.text(6.78, 0.006, "密集检查点交点\n6.734 V", rotation=90, ha="left", va="bottom", fontsize=8.8, color="#8A661B")
    ax_iic.set_yscale("log")
    ax_iic.set_xlim(0.8, 7.2)
    ax_iic.set_ylim(0.0015, 30.0)
    ax_iic.set_xlabel("漏极偏压 Vd (V)")
    ax_iic.set_ylabel("电流 (µA/µm，对数轴)")
    ax_iic.set_title("IIC 电流交点与当前 Vela 分支", loc="left", fontsize=13, weight="bold")
    ax_iic.grid(axis="y", which="both")
    ax_iic.legend(loc="upper left", fontsize=9)

    y = np.arange(len(ratio_labels))
    colors = [BLUE, BLUE, GOLD, ORANGE]
    for yi, ratio, color in zip(y, ratios, colors):
        lo, hi = sorted((1.0, ratio))
        ax_ratio.hlines(yi, lo, hi, color=color, linewidth=3.2, alpha=0.85)
        ax_ratio.scatter([ratio], [yi], s=78, color=color, edgecolor="white", linewidth=1.0, zorder=4)
        label = f"{ratio:.3f}×" if ratio >= 0.01 else f"{ratio:.2e}×"
        ax_ratio.annotate(label, xy=(ratio, yi), xytext=(7 if ratio >= 1 else -7, 0), textcoords="offset points", ha="left" if ratio >= 1 else "right", va="center", weight="bold", color=color)
    ax_ratio.axvline(1.0, color=CHARCOAL, linewidth=1.4, linestyle="--")
    ax_ratio.set_xscale("log")
    ax_ratio.set_xlim(1e-6, 5.0)
    ax_ratio.set_yticks(y, ratio_labels)
    ax_ratio.invert_yaxis()
    ax_ratio.set_xlabel("Vela / Sentaurus 峰值比（对数轴）")
    ax_ratio.set_title("6.38 V 旧分支根因诊断：场和 α 接近，雪崩源缺少电流支撑", loc="left", fontsize=13, weight="bold")
    ax_ratio.grid(axis="x", which="both")
    ax_ratio.text(0.0, -0.30, "注：该面板采用修正单位后的旧非闭合分支，仅用于定位差异源头，不作为最终验收结果。", transform=ax_ratio.transAxes, fontsize=9.5, color=MID_GREY)

    fig.suptitle("NMOS 关键物理量对比与 IIC 闭合进展", x=0.055, y=0.99, ha="left", fontsize=18, weight="bold")
    fig.text(0.055, 0.94, "2026-08-02 至 2026-08-03；电流按每微米器件深度归一化", fontsize=10.5, color=MID_GREY)
    fig.subplots_adjust(left=0.095, right=0.97, top=0.875, bottom=0.10)
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=repo_root() / "docs" / "validation" / "daily_report_2026-08-03")
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    root = repo_root()
    run_root = root / "build-release" / "reference_tcad" / "bvmethods_sentaurus2018" / "run01"
    mesh_path = run_root / "vela" / "mesh.json"
    iic_path = run_root / "vela_validation" / "iic_postprocess_20260803" / "analysis" / "multibias_sentaurus" / "sentaurus_exact_extended_curve.csv"

    draw_mesh(mesh_path, args.output_dir / "nmos_mesh_overview.png")
    draw_physics(iic_path, args.output_dir / "nmos_key_physics_comparison.png")
    write_metric_source(args.output_dir / "daily_report_key_metrics.csv")


if __name__ == "__main__":
    configure_style()
    main()
