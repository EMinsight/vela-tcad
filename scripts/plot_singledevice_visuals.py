#!/usr/bin/env python3
"""Render SingleDevice mesh and validation comparison as static PNG files."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PolyCollection
from matplotlib.patches import Patch
import numpy as np


ROOT = Path(__file__).resolve().parents[1]
BUILD = ROOT / "build-release/reference_tcad/singledevice_sentaurus2018"
OUT = Path.home() / ".codex/visualizations/2026/08/11/019fefda-670a-7aa0-86f2-4e6a8dda6c61"


def load_curve(path: Path, current_column: str) -> tuple[np.ndarray, np.ndarray]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    return (
        np.asarray([float(row["bias_V"]) for row in rows]),
        np.asarray([float(row[current_column]) for row in rows]),
    )


def plot_mesh() -> Path:
    mesh_path = BUILD / "vela_import_fixedmaterials/vela/mesh.json"
    mesh = json.loads(mesh_path.read_text(encoding="utf-8"))
    nodes = {int(node["id"]): (float(node["y"]), float(node["x"])) for node in mesh["nodes"]}
    triangles = mesh["triangles"]
    region_names = {int(r["id"]): (r["name"], r["material"]) for r in mesh["regions"]}
    colors = {
        0: "#79a9dc", 1: "#79a9dc", 2: "#79a9dc",
        3: "#7fc38b", 4: "#db9a59", 5: "#b190d4", 6: "#b190d4",
    }
    polys_by_region: dict[int, list[list[tuple[float, float]]]] = {i: [] for i in range(7)}
    segments: list[list[tuple[float, float]]] = []
    for tri in triangles:
        poly = [nodes[int(node_id)] for node_id in tri["node_ids"]]
        polys_by_region[int(tri["region_id"])].append(poly)
        segments.extend([[poly[0], poly[1]], [poly[1], poly[2]], [poly[2], poly[0]]])

    fig, axes = plt.subplots(1, 2, figsize=(14, 6.8))
    views = [((-0.20, 0.20), (1.03, -0.12), "Full device"),
             ((-0.11, 0.11), (0.14, -0.112), "Gate-channel mesh zoom")]
    for ax, (xlim, ylim, title) in zip(axes, views):
        for region_id, polys in polys_by_region.items():
            ax.add_collection(PolyCollection(
                polys, facecolor=colors[region_id], edgecolor="none", alpha=0.42))
        ax.add_collection(LineCollection(segments, colors="#314252", linewidths=0.22, alpha=0.45))
        for contact in mesh["contacts"]:
            points = np.asarray([nodes[int(i)] for i in contact["node_ids"]])
            order = np.argsort(points[:, 0] + points[:, 1] * 1e-5)
            ax.plot(points[order, 0], points[order, 1], color="#111111", lw=3.2)
            center = points.mean(axis=0)
            offset = -0.018 if contact["name"] == "substrate" else 0.012
            ax.text(center[0], center[1] + offset, contact["name"].capitalize(),
                    ha="center", va="center", fontsize=9, weight="bold")
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.set_xlabel("Lateral coordinate y (um)")
        ax.set_ylabel("Depth coordinate x (um)")
        ax.set_title(title, weight="bold")
        ax.grid(False)
        ax.set_facecolor("#fbfcfe")
    axes[0].text(0, 0.52, "R.Substrate / Si", ha="center", fontsize=11)
    axes[1].text(0, -0.055, "Poly-Si", ha="center", fontsize=10)
    axes[1].text(0, 0.065, "Channel / Si", ha="center", fontsize=10)
    axes[1].axhline(0.0, color="#1b1b1b", lw=1.25)
    handles = [
        Patch(color=colors[3], alpha=.55, label="Si substrate"),
        Patch(color=colors[4], alpha=.55, label="Poly-Si gate"),
        Patch(color=colors[0], alpha=.55, label="SiO2"),
        Patch(color=colors[5], alpha=.55, label="Si3N4 spacer"),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=4, frameon=False,
               bbox_to_anchor=(0.5, 0.925))
    fig.suptitle("SingleDevice 2-D structure and imported triangular mesh",
                 fontsize=15, weight="bold", y=0.98)
    fig.subplots_adjust(left=.065, right=.985, bottom=.09, top=.82, wspace=.22)
    path = OUT / "singledevice-mesh.png"
    fig.savefig(path, dpi=190, facecolor="white")
    plt.close(fig)
    return path


def plot_results() -> Path:
    ref_dir = BUILD / "vela_import_fixedmaterials/reference_curves"
    candidate_dir = BUILD / "fixed_state_curve"
    lin_s = load_curve(ref_dir / "singledevice_sentaurus2018_idvg_lin_reference.csv", "current_total")
    sat_s = load_curve(ref_dir / "singledevice_sentaurus2018_idvg_sat_reference.csv", "current_total")
    lin_v = load_curve(candidate_dir / "lin_candidate.csv", "current_total_A_per_um")
    sat_v = load_curve(candidate_dir / "sat_candidate.csv", "current_total_A_per_um")
    lin_err = 100 * np.abs(lin_v[1] - lin_s[1]) / np.abs(lin_s[1])
    sat_err = 100 * np.abs(sat_v[1] - sat_s[1]) / np.abs(sat_s[1])

    fig = plt.figure(figsize=(14, 7.2))
    grid = fig.add_gridspec(2, 2, width_ratios=[1.75, 0.85], height_ratios=[1.45, .75])
    ax_iv = fig.add_subplot(grid[0, 0])
    ax_err = fig.add_subplot(grid[1, 0], sharex=ax_iv)
    ax_end = fig.add_subplot(grid[:, 1])
    blue, orange = "#3078bd", "#d36f22"
    ax_iv.semilogy(*lin_s, "o-", color=blue, ms=3.4, lw=1.8, label="Linear Sentaurus")
    ax_iv.semilogy(*lin_v, "o--", color=blue, mfc="white", ms=3.4, lw=1.8, label="Linear Vela")
    ax_iv.semilogy(*sat_s, "s-", color=orange, ms=3.4, lw=1.8, label="Saturation Sentaurus")
    ax_iv.semilogy(*sat_v, "s--", color=orange, mfc="white", ms=3.4, lw=1.8, label="Saturation Vela")
    ax_iv.set_ylabel("Drain current Id (A/um)")
    ax_iv.set_ylim(7e-15, 3e-3)
    ax_iv.grid(True, which="both", alpha=.22)
    ax_iv.legend(frameon=False, ncol=2, loc="upper left")
    ax_iv.set_title("Id-Vg: frozen imported quantum-potential states", weight="bold")

    ax_err.plot(lin_s[0], lin_err, "o--", color=blue, ms=3.2, label="Linear")
    ax_err.plot(sat_s[0], sat_err, "s--", color=orange, ms=3.2, label="Saturation")
    ax_err.set_xlabel("Gate voltage Vg (V)")
    ax_err.set_ylabel("Absolute relative error (%)")
    ax_err.set_ylim(0, 10)
    ax_err.grid(True, alpha=.25)
    ax_err.text(2.18, lin_err.max() + .35, f"max {lin_err.max():.2f}%", color=blue, ha="right")
    ax_err.text(2.18, sat_err.max() + .35, f"max {sat_err.max():.2f}%", color=orange, ha="right")

    labels = ["Linear\nVd=0.1 V", "Saturation\nVd=1.1 V"]
    changes = [2.0993, 3.5866]
    residuals = [1.1854e6, 1.0538e6]
    bars = ax_end.bar(labels, changes, color=[blue, orange], alpha=.72, width=.62)
    ax_end.axhline(5e-4, color="#333333", ls=":", lw=1.3,
                   label="Acceptance 5e-4 V")
    ax_end.set_ylabel("Raw quantum-potential change (V)")
    ax_end.set_ylim(0, 4.15)
    ax_end.set_title("Self-consistent Eq. 231 endpoints", weight="bold")
    ax_end.grid(True, axis="y", alpha=.22)
    for bar, value, residual in zip(bars, changes, residuals):
        ax_end.text(bar.get_x() + bar.get_width()/2, value + .10, f"{value:.3f} V",
                    ha="center", weight="bold")
        ax_end.text(bar.get_x() + bar.get_width()/2, value * .58,
                    f"NOT CONVERGED\n500 iterations\nresidual {residual:.3e}",
                    ha="center", va="center", color="white", fontsize=9, weight="bold")
    ax_end.legend(frameon=False, loc="upper left")
    fig.suptitle("SingleDevice Sentaurus-Vela validation comparison", fontsize=15,
                 weight="bold", y=.98)
    fig.subplots_adjust(left=.065, right=.985, bottom=.11, top=.88,
                        wspace=.14, hspace=.13)
    path = OUT / "singledevice-results-comparison.png"
    fig.savefig(path, dpi=190, facecolor="white")
    plt.close(fig)
    return path


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    print(plot_mesh())
    print(plot_results())
