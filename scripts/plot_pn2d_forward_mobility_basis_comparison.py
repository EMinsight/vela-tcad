#!/usr/bin/env python3
"""Plot Sentaurus and Vela mobility-doping-basis forward-IV curves."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


CURVES = (
    (
        "net_doping",
        "Vela net doping (production default)",
        "#4C78A8",
        "--",
    ),
    (
        "total_impurity",
        "Vela total impurity",
        "#54A24B",
        "-.",
    ),
    (
        "cell_reconstructed_total_impurity",
        "Vela cell-reconstructed total impurity",
        "#E45756",
        "-",
    ),
)


def read_sentaurus(path: Path) -> tuple[np.ndarray, np.ndarray]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.DictReader(stream))
    by_bias: dict[float, float] = {}
    for row in rows:
        bias = float(row["bias_V"])
        current = abs(float(row["current_total"]))
        if bias not in by_bias or current > by_bias[bias]:
            by_bias[bias] = current
    biases = np.asarray(sorted(by_bias))
    return biases, np.asarray([by_bias[bias] for bias in biases])


def read_vela(path: Path) -> tuple[np.ndarray, np.ndarray]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        rows = [
            row
            for row in csv.DictReader(stream)
            if row["converged"] == "1"
        ]
    by_bias = {
        round(float(row["bias_V"]), 12):
        abs(float(row["current_total_A_per_um"]))
        for row in rows
    }
    biases = np.asarray(sorted(by_bias))
    return biases, np.asarray([by_bias[bias] for bias in biases])


def read_exact_anchors(
    path: Path,
) -> tuple[np.ndarray, np.ndarray, dict[str, np.ndarray]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.DictReader(stream))
    bases = [basis for basis, *_ in CURVES]
    by_basis = {
        basis: {
            float(row["bias_V"]): float(row["vela_A_per_um"])
            for row in rows
            if row["basis"] == basis
        }
        for basis in bases
    }
    sentaurus_by_bias = {
        float(row["bias_V"]): float(row["sentaurus_A_per_um"])
        for row in rows
    }
    biases = np.asarray(sorted(sentaurus_by_bias))
    sentaurus = np.asarray([sentaurus_by_bias[bias] for bias in biases])
    vela = {
        basis: np.asarray([by_basis[basis][bias] for bias in biases])
        for basis in bases
    }
    return biases, sentaurus, vela


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sentaurus", type=Path, required=True)
    parser.add_argument("--candidate-root", type=Path, required=True)
    parser.add_argument("--anchor-comparison", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()

    sx, sy = read_sentaurus(args.sentaurus)
    anchor_x, anchor_sentaurus, anchor_vela = read_exact_anchors(
        args.anchor_comparison
    )
    vela_curves: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for basis, _, _, _ in CURVES:
        path = (
            args.candidate_root
            / basis
            / f"vela_{basis}.csv"
        )
        vela_curves[basis] = read_vela(path)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    comparison_path = (
        args.out_dir / "forward_iv_mobility_basis_comparison.csv"
    )
    grid = vela_curves["net_doping"][0]
    with comparison_path.open("w", newline="", encoding="utf-8") as stream:
        fields = [
            "bias_V",
            "sentaurus_A_per_um",
            *[f"{basis}_A_per_um" for basis, *_ in CURVES],
            *[f"{basis}_relative_error_percent" for basis, *_ in CURVES],
        ]
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        sent_interp = np.interp(grid, sx, sy)
        for index, bias in enumerate(grid):
            row: dict[str, float] = {
                "bias_V": float(bias),
                "sentaurus_A_per_um": float(sent_interp[index]),
            }
            for basis, *_ in CURVES:
                vx, vy = vela_curves[basis]
                value = float(np.interp(bias, vx, vy))
                row[f"{basis}_A_per_um"] = value
                row[f"{basis}_relative_error_percent"] = (
                    100.0 * (value - sent_interp[index])
                    / sent_interp[index]
                    if sent_interp[index] > 1e-30
                    else float("nan")
                )
            writer.writerow(row)

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9.5,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "legend.fontsize": 8.5,
        }
    )
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.4), dpi=180)
    sentaurus_color = "#B279A2"

    log_ax, linear_ax, error_ax = axes
    for ax in (log_ax, linear_ax):
        plot = ax.semilogy if ax is log_ax else ax.plot
        plot(
            sx,
            np.maximum(sy, 1e-30) if ax is log_ax else sy * 1e3,
            color=sentaurus_color,
            linewidth=2.4,
            label="Sentaurus 2018 (low-field)",
            zorder=5,
        )
        for basis, label, color, linestyle in CURVES:
            vx, vy = vela_curves[basis]
            plot(
                vx,
                np.maximum(vy, 1e-30) if ax is log_ax else vy * 1e3,
                color=color,
                linestyle=linestyle,
                linewidth=1.8,
                label=label,
            )
        plot(
            anchor_x,
            (
                np.maximum(anchor_sentaurus, 1e-30)
                if ax is log_ax
                else anchor_sentaurus * 1e3
            ),
            linestyle="none",
            marker="o",
            markersize=4.5,
            markerfacecolor="none",
            markeredgecolor=sentaurus_color,
            markeredgewidth=1.2,
            label="Sentaurus exact-bias states",
            zorder=6,
        )
        ax.set_xlim(0, 20)
        ax.set_xlabel("Anode voltage (V)")
        ax.grid(True, which="major", color="#D8D8D8", linewidth=0.65)
        if ax is log_ax:
            ax.grid(True, which="minor", color="#EEEEEE", linewidth=0.35)

    log_ax.set_title("Absolute current — log scale")
    log_ax.set_ylabel("|Anode current| (A/µm)")
    log_ax.set_ylim(
        1e-20,
        max(
            float(sy.max()),
            *(float(values.max()) for _, values in vela_curves.values()),
        )
        * 1.8,
    )
    log_ax.legend(frameon=False, loc="lower right")

    linear_ax.set_title("Absolute current — linear scale")
    linear_ax.set_ylabel("|Anode current| (mA/µm)")
    linear_ax.set_ylim(bottom=0)

    for basis, label, color, linestyle in CURVES:
        relative_error = (
            100.0
            * (anchor_vela[basis] - anchor_sentaurus)
            / anchor_sentaurus
        )
        error_ax.plot(
            anchor_x,
            relative_error,
            color=color,
            linestyle=linestyle,
            marker="o",
            markersize=3.5,
            linewidth=1.8,
            label=label,
        )
    error_ax.axhline(0.0, color="#777777", linewidth=0.9)
    error_ax.set_title("Relative difference from Sentaurus")
    error_ax.set_xlabel("Anode voltage (V)")
    error_ax.set_ylabel("(Vela − Sentaurus) / Sentaurus (%)")
    error_ax.set_xlim(1, 20)
    error_ax.grid(True, color="#E0E0E0", linewidth=0.65)

    fig.suptitle(
        "PN2D coarse7x3 forward I–V: mobility doping-basis comparison",
        fontsize=14,
        y=0.98,
    )
    fig.text(
        0.5,
        0.015,
        "Avalanche off; SRH on; low-field Masetti mobility; "
        "same mesh and current normalization.",
        ha="center",
        color="#555555",
        fontsize=8.5,
    )
    fig.tight_layout(rect=(0, 0.05, 1, 0.94))
    stem = args.out_dir / "pn2d_coarse7x3_forward_iv_mobility_basis_0v20v"
    fig.savefig(stem.with_suffix(".png"))
    fig.savefig(stem.with_suffix(".svg"))
    plt.close(fig)


if __name__ == "__main__":
    main()
