#!/usr/bin/env python3
"""Render the TransportModels mesh and completed DG DC comparisons."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


REPO_ROOT = Path(__file__).resolve().parents[1]
BASELINE_DIR = (
    REPO_ROOT
    / "build-release/reference_tcad/transportmodels_sentaurus2022/vela_baseline"
)
DEFAULT_OUTPUT_DIR = (
    REPO_ROOT
    / "docs/progress_report_2026Q3/2026-08-19_transportmodels_dg_daily_report/figures"
)
DEFAULT_DG_IDVD_WORKFLOW_DIR = (
    BASELINE_DIR
    / "workflow_dg_idvd_resume_1p9_direct2p0_fixed1p5_outer200_run01"
)


def configure_style() -> None:
    plt.rcParams.update(
        {
            "font.family": ["Microsoft YaHei", "DejaVu Sans"],
            "axes.unicode_minus": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.edgecolor": "#343A40",
            "axes.labelcolor": "#252A31",
            "axes.titlecolor": "#20252B",
            "xtick.color": "#4E5661",
            "ytick.color": "#4E5661",
            "grid.color": "#D9DEE5",
            "grid.linewidth": 0.7,
            "grid.alpha": 0.8,
            "legend.frameon": False,
            "savefig.facecolor": "white",
            "savefig.bbox": "tight",
        }
    )


def read_current_csv(path: Path, current_field: str) -> dict[float, float]:
    rows: dict[float, float] = {}
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        for row in csv.DictReader(stream):
            rows[float(row["bias_V"])] = abs(float(row[current_field]))
    return rows


def load_dg_idvg() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    reference_path = (
        BASELINE_DIR
        / "generated/reference_curves/transportmodels_sentaurus2022_dg_idvg_reference.csv"
    )
    prefix_path = BASELINE_DIR / "dg_idvg_prefix_through_m036.csv"
    resume_path = (
        BASELINE_DIR
        / "workflow_dg_outer80_resume_m036_run01/dg_idvg_curve.csv"
    )

    reference = read_current_csv(reference_path, "current_total")
    vela = read_current_csv(prefix_path, "current_total_A_per_um")
    vela.update(read_current_csv(resume_path, "current_total_A_per_um"))

    reference_biases = sorted(reference)
    if len(reference_biases) != 21 or len(vela) != 21:
        raise RuntimeError(
            f"Expected two 21-point Id-Vg curves, got reference={len(reference)} "
            f"and Vela={len(vela)}"
        )
    for bias in reference_biases:
        if not any(math.isclose(bias, candidate, abs_tol=1.0e-12) for candidate in vela):
            raise RuntimeError(f"Vela curve is missing reference bias {bias:g} V")

    vela_by_reference_bias = []
    for bias in reference_biases:
        key = min(vela, key=lambda candidate: abs(candidate - bias))
        vela_by_reference_bias.append(vela[key])

    return (
        np.asarray(reference_biases),
        np.asarray([reference[bias] for bias in reference_biases]),
        np.asarray(vela_by_reference_bias),
    )


def load_dg_idvd(
    workflow_dir: Path,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    reference_path = (
        BASELINE_DIR
        / "generated/reference_curves/transportmodels_sentaurus2022_dg_idvd_reference.csv"
    )
    candidate_path = workflow_dir / "dg_idvd_curve_comparison_candidate.csv"

    reference = read_current_csv(reference_path, "current_total")
    vela = read_current_csv(candidate_path, "current_total_A_per_um")
    reference_biases = sorted(reference)
    if len(reference_biases) != 21 or len(vela) != 21:
        raise RuntimeError(
            f"Expected two 21-point Id-Vd curves, got reference={len(reference)} "
            f"and Vela={len(vela)}"
        )
    for bias in reference_biases:
        if not any(math.isclose(
                bias, candidate, abs_tol=1.0e-12) for candidate in vela):
            raise RuntimeError(f"Vela Id-Vd curve is missing reference bias {bias:g} V")

    vela_by_reference_bias = []
    for bias in reference_biases:
        key = min(vela, key=lambda candidate: abs(candidate - bias))
        vela_by_reference_bias.append(vela[key])

    return (
        np.asarray(reference_biases),
        np.asarray([reference[bias] for bias in reference_biases]),
        np.asarray(vela_by_reference_bias),
    )


def render_idvg(output_dir: Path) -> dict[str, float]:
    bias, sentaurus, vela = load_dg_idvg()
    signed_log_ratio = np.log10(vela / sentaurus)

    fig = plt.figure(figsize=(10.8, 8.2), constrained_layout=True)
    grid = fig.add_gridspec(2, 1, height_ratios=[3.2, 1.15])
    ax = fig.add_subplot(grid[0])
    error_ax = fig.add_subplot(grid[1], sharex=ax)

    ax.semilogy(
        bias,
        sentaurus,
        color="#245B8A",
        linewidth=2.4,
        marker="o",
        markersize=5.5,
        markerfacecolor="white",
        markeredgewidth=1.5,
        label="Sentaurus 2022",
        zorder=3,
    )
    ax.semilogy(
        bias,
        vela,
        color="#D27A1E",
        linewidth=2.1,
        linestyle="--",
        marker="s",
        markersize=4.8,
        markerfacecolor="#D27A1E",
        markeredgecolor="white",
        markeredgewidth=0.6,
        label="Vela DG",
        zorder=4,
    )
    ax.set_title("TransportModels 电子 DG MOS：Id–Vg 对比", fontsize=16, pad=34)
    ax.text(
        0.0,
        1.01,
        "Vd = 1.1 V；21 个完全一致的栅压点；漏极总电流按器件宽度归一化",
        transform=ax.transAxes,
        fontsize=10.5,
        color="#59626D",
        va="bottom",
    )
    ax.set_ylabel("|Id|  (A/μm)", fontsize=11.5)
    ax.set_ylim(1.0e-24, 1.0e-2)
    ax.grid(True, which="major")
    ax.grid(True, which="minor", linewidth=0.35, alpha=0.35)
    ax.legend(loc="lower right", fontsize=10.5)
    ax.tick_params(axis="x", labelbottom=False)

    error_ax.axhline(0.0, color="#343A40", linewidth=1.0, zorder=1)
    error_ax.plot(
        bias,
        signed_log_ratio,
        color="#8A6F24",
        linewidth=1.8,
        marker="D",
        markersize=4.0,
        markerfacecolor="white",
        markeredgewidth=1.0,
        zorder=2,
    )
    error_ax.fill_between(
        bias,
        signed_log_ratio,
        0.0,
        color="#E9D9A5",
        alpha=0.45,
        zorder=0,
    )
    error_ax.set_xlabel("栅极电压 Vg  (V)", fontsize=11.5)
    error_ax.set_ylabel("log10\n(Vela/Sentaurus)", fontsize=10)
    error_ax.set_xlim(-1.04, 2.24)
    error_ax.set_xticks(np.arange(-1.0, 2.21, 0.4))
    error_ax.set_ylim(min(-9.3, float(signed_log_ratio.min()) - 0.4), 0.9)
    error_ax.grid(True, which="major")
    error_ax.text(
        0.995,
        0.08,
        "0 表示完全一致；负值表示 Vela 电流较小",
        transform=error_ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=9.2,
        color="#626A74",
    )

    fig.text(
        0.01,
        -0.012,
        "数据：Sentaurus T-2022.03-SP2 TransportModels DG 参考曲线；Vela DG 已完成 Id–Vg 结果。",
        fontsize=9.2,
        color="#626A74",
    )

    for suffix in ("png", "svg"):
        fig.savefig(
            output_dir / f"transportmodels_dg_idvg_comparison.{suffix}",
            dpi=240 if suffix == "png" else None,
        )
    plt.close(fig)

    return {
        "points": float(len(bias)),
        "min_signed_log_ratio": float(signed_log_ratio.min()),
        "max_signed_log_ratio": float(signed_log_ratio.max()),
        "median_abs_log_error": float(np.median(np.abs(signed_log_ratio))),
    }


def render_idvd(output_dir: Path, workflow_dir: Path) -> dict[str, float]:
    bias, sentaurus, vela = load_dg_idvd(workflow_dir)
    nonzero = bias > 0.0
    relative_error = np.zeros_like(bias)
    relative_error[nonzero] = (vela[nonzero] - sentaurus[nonzero]) / sentaurus[nonzero]

    fig = plt.figure(figsize=(10.8, 8.2), constrained_layout=True)
    grid = fig.add_gridspec(2, 1, height_ratios=[3.2, 1.15])
    ax = fig.add_subplot(grid[0])
    error_ax = fig.add_subplot(grid[1], sharex=ax)

    ax.plot(
        bias,
        sentaurus,
        color="#245B8A",
        linewidth=2.4,
        marker="o",
        markersize=5.5,
        markerfacecolor="white",
        markeredgewidth=1.5,
        label="Sentaurus 2022",
        zorder=3,
    )
    ax.plot(
        bias,
        vela,
        color="#D27A1E",
        linewidth=2.1,
        linestyle="--",
        marker="s",
        markersize=4.8,
        markerfacecolor="#D27A1E",
        markeredgecolor="white",
        markeredgewidth=0.6,
        label="Vela DG",
        zorder=4,
    )
    ax.set_title("TransportModels 电子 DG MOS：Id–Vd 对比", fontsize=16, pad=34)
    ax.text(
        0.0,
        1.01,
        "Vg = 1.0 V；21 个完全一致的漏压点；漏极总电流按器件宽度归一化",
        transform=ax.transAxes,
        fontsize=10.5,
        color="#59626D",
        va="bottom",
    )
    ax.set_ylabel("|Id|  (A/μm)", fontsize=11.5)
    ax.tick_params(axis="x", labelbottom=False)
    ax.grid(True, which="major")
    ax.legend(loc="lower right", fontsize=10.5)

    error_ax.axhline(0.0, color="#343A40", linewidth=1.0, zorder=1)
    error_ax.plot(
        bias[nonzero],
        100.0 * relative_error[nonzero],
        color="#8A6F24",
        linewidth=1.8,
        marker="D",
        markersize=4.0,
        markerfacecolor="white",
        markeredgewidth=1.0,
        zorder=2,
    )
    error_ax.fill_between(
        bias[nonzero],
        100.0 * relative_error[nonzero],
        0.0,
        color="#E9D9A5",
        alpha=0.45,
        zorder=0,
    )
    error_ax.set_xlabel("漏极电压 Vd  (V)", fontsize=11.5)
    error_ax.set_ylabel("相对误差\n(%)", fontsize=10)
    error_ax.set_xlim(-0.04, 2.04)
    error_ax.set_xticks(np.arange(0.0, 2.01, 0.2))
    error_ax.grid(True, which="major")
    error_ax.text(
        0.995,
        0.08,
        "正值表示 Vela 电流较大；0 V 点不计算相对误差",
        transform=error_ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=9.2,
        color="#626A74",
    )

    fig.text(
        0.01,
        -0.012,
        "数据：Sentaurus T-2022.03-SP2 TransportModels DG 参考曲线；Vela DG 完整 Id–Vd 结果。",
        fontsize=9.2,
        color="#626A74",
    )

    for suffix in ("png", "svg"):
        fig.savefig(
            output_dir / f"transportmodels_dg_idvd_comparison.{suffix}",
            dpi=240 if suffix == "png" else None,
        )
    plt.close(fig)

    return {
        "points": float(len(bias)),
        "median_abs_relative_error": float(
            np.median(np.abs(relative_error[nonzero]))
        ),
        "max_abs_relative_error": float(np.max(np.abs(relative_error[nonzero]))),
        "endpoint_relative_error": float(relative_error[-1]),
    }


def render_mesh(output_dir: Path) -> dict[str, int]:
    mesh_path = BASELINE_DIR / "generated/vela/mesh.json"
    with mesh_path.open("r", encoding="utf-8") as stream:
        mesh = json.load(stream)

    node_by_id = {node["id"]: node for node in mesh["nodes"]}
    node_ids = [node["id"] for node in mesh["nodes"]]
    node_index = {node_id: index for index, node_id in enumerate(node_ids)}

    # Sentaurus coordinates use x as depth and y as the lateral direction.
    lateral = np.asarray([node_by_id[node_id]["y"] for node_id in node_ids])
    depth = np.asarray([node_by_id[node_id]["x"] for node_id in node_ids])
    triangles = np.asarray(
        [
            [node_index[node_id] for node_id in triangle["node_ids"]]
            for triangle in mesh["triangles"]
        ],
        dtype=int,
    )
    region_id = np.asarray([triangle["region_id"] for triangle in mesh["triangles"]])
    material_by_region = {
        region["id"]: region["material"] for region in mesh["regions"]
    }
    materials = ["Si", "SiO2", "PolySilicon", "Nitride"]
    material_code = {material: index for index, material in enumerate(materials)}
    face_values = np.asarray(
        [material_code[material_by_region[int(identifier)]] for identifier in region_id]
    )

    material_colors = {
        "Si": "#DCE8F2",
        "SiO2": "#F3E4B5",
        "PolySilicon": "#E4C8D8",
        "Nitride": "#D7E3C5",
    }
    cmap = ListedColormap([material_colors[material] for material in materials])
    norm = BoundaryNorm(np.arange(-0.5, len(materials) + 0.5), cmap.N)
    triangulation = mtri.Triangulation(lateral, depth, triangles)

    fig, axes = plt.subplots(1, 2, figsize=(13.4, 7.8), constrained_layout=True)
    contact_colors = {
        "gate": "#245B8A",
        "source": "#C96A1B",
        "drain": "#C96A1B",
        "substrate": "#343A40",
    }

    for ax, title, y_limits in (
        (axes[0], "全器件网格", (1.04, -0.12)),
        (axes[1], "表面沟道与源漏区放大", (0.22, -0.12)),
    ):
        ax.tripcolor(
            triangulation,
            facecolors=face_values,
            cmap=cmap,
            norm=norm,
            edgecolors="none",
            shading="flat",
            zorder=0,
        )
        ax.triplot(
            triangulation,
            color="#626B75",
            linewidth=0.22,
            alpha=0.72,
            zorder=1,
        )
        for contact in mesh["contacts"]:
            coordinates = sorted(
                (
                    node_by_id[node_id]["y"],
                    node_by_id[node_id]["x"],
                )
                for node_id in contact["node_ids"]
            )
            x_contact, y_contact = zip(*coordinates)
            ax.plot(
                x_contact,
                y_contact,
                color=contact_colors[contact["name"]],
                linewidth=3.0,
                solid_capstyle="round",
                zorder=4,
            )
        ax.set_title(title, fontsize=14, pad=10)
        ax.set_xlabel("横向位置  (μm)", fontsize=11)
        ax.set_ylabel("表面法向深度  (μm)", fontsize=11)
        ax.set_xlim(-0.2, 0.2)
        ax.set_ylim(*y_limits)
        ax.set_aspect("equal", adjustable="box")
        ax.grid(False)

    axes[0].annotate("Gate", (0.0, -0.1012), xytext=(0.0, -0.075), ha="center", fontsize=9.5, color="#245B8A")
    axes[0].annotate("Source", (-0.135, 0.0), xytext=(-0.135, 0.045), ha="center", fontsize=9.5, color="#A65316")
    axes[0].annotate("Drain", (0.135, 0.0), xytext=(0.135, 0.045), ha="center", fontsize=9.5, color="#A65316")
    axes[0].annotate("Substrate", (0.0, 1.0), xytext=(0.0, 0.965), ha="center", fontsize=9.5, color="#343A40")
    axes[1].annotate("Gate", (0.0, -0.1012), xytext=(0.0, -0.073), ha="center", fontsize=9.5, color="#245B8A")
    axes[1].annotate("Source", (-0.135, 0.0), xytext=(-0.135, 0.037), ha="center", fontsize=9.5, color="#A65316")
    axes[1].annotate("Drain", (0.135, 0.0), xytext=(0.135, 0.037), ha="center", fontsize=9.5, color="#A65316")

    material_legend = [
        Patch(
            facecolor=material_colors[material],
            edgecolor="#59626D",
            linewidth=0.5,
            label=material,
        )
        for material in materials
    ]
    contact_legend = [
        Line2D([0], [0], color="#245B8A", linewidth=3.0, label="Gate contact"),
        Line2D([0], [0], color="#C96A1B", linewidth=3.0, label="Source/Drain contacts"),
        Line2D([0], [0], color="#343A40", linewidth=3.0, label="Substrate contact"),
    ]
    fig.legend(
        handles=material_legend + contact_legend,
        loc="outside lower center",
        ncol=7,
        fontsize=9.3,
        handlelength=1.8,
        columnspacing=1.25,
    )
    fig.suptitle(
        "TransportModels MOS 器件有限体积网格（3315 节点 / 6456 三角形）",
        fontsize=16,
        y=1.025,
    )

    for suffix in ("png", "svg"):
        fig.savefig(
            output_dir / f"transportmodels_device_mesh.{suffix}",
            dpi=240 if suffix == "png" else None,
        )
    plt.close(fig)

    return {
        "nodes": len(mesh["nodes"]),
        "triangles": len(mesh["triangles"]),
        "regions": len(mesh["regions"]),
        "contacts": len(mesh["contacts"]),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--idvd-workflow-dir",
        type=Path,
        default=DEFAULT_DG_IDVD_WORKFLOW_DIR,
    )
    parser.add_argument(
        "--summary-json",
        type=Path,
        help=("figure QA summary and chart map; defaults to "
              "OUTPUT_DIR/transportmodels_dg_figure_summary.json"),
    )
    args = parser.parse_args()

    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    configure_style()

    mesh_summary = render_mesh(output_dir)
    idvg_summary = render_idvg(output_dir)
    idvd_summary = render_idvd(output_dir, args.idvd_workflow_dir.resolve())
    summary = {
        "output_dir": str(output_dir),
        "sources": {
            "idvg_reference": str((
                BASELINE_DIR
                / "generated/reference_curves/transportmodels_sentaurus2022_dg_idvg_reference.csv"
            ).resolve()),
            "idvg_candidate": str((
                BASELINE_DIR
                / "workflow_dg_outer80_resume_m036_run01/dg_idvg_curve_comparison_candidate.csv"
            ).resolve()),
            "idvd_reference": str((
                BASELINE_DIR
                / "generated/reference_curves/transportmodels_sentaurus2022_dg_idvd_reference.csv"
            ).resolve()),
            "idvd_candidate": str((
                args.idvd_workflow_dir.resolve()
                / "dg_idvd_curve_comparison_candidate.csv"
            )),
        },
        "mesh": mesh_summary,
        "idvg": idvg_summary,
        "idvd": idvd_summary,
        "chart_map": [
            {
                "segment": "device mesh",
                "question": "What device geometry and mesh support the comparison?",
                "family": "matrix and structure",
                "type": "categorical triangular mesh with overview and surface zoom",
                "fields": ["node coordinates", "triangle region", "contact"],
                "supported_claim": "DD and DG comparisons share the same imported 2-D mesh.",
                "palette_policy": "relaxed multi-category materials plus neutral contacts",
                "outputs": ["transportmodels_device_mesh.png", "transportmodels_device_mesh.svg"],
            },
            {
                "segment": "DG Id-Vg validation",
                "question": "How closely does Vela reproduce the Sentaurus gate sweep?",
                "family": "trend and benchmark",
                "type": "highlighted two-series semilog line with signed log-error panel",
                "fields": ["gate bias", "drain current", "signed log10 ratio"],
                "supported_claim": "The on-state shape is close, while deep-off current is below Vela's numerical floor.",
                "palette_policy": "hard two-root cap with marker and line-style redundancy",
                "outputs": ["transportmodels_dg_idvg_comparison.png", "transportmodels_dg_idvg_comparison.svg"],
            },
            {
                "segment": "DG Id-Vd validation",
                "question": "How closely does Vela reproduce the Sentaurus drain sweep?",
                "family": "trend and benchmark",
                "type": "highlighted two-series line with nonzero-bias relative-error panel",
                "fields": ["drain bias", "drain current", "relative error"],
                "supported_claim": "The complete nonzero-bias drain sweep is monotonic and remains within 10.72% relative error.",
                "palette_policy": "hard two-root cap with marker and line-style redundancy",
                "outputs": ["transportmodels_dg_idvd_comparison.png", "transportmodels_dg_idvd_comparison.svg"],
            },
        ],
    }
    summary_path = (args.summary_json or (
        output_dir / "transportmodels_dg_figure_summary.json")).resolve()
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
