#!/usr/bin/env python3
"""Merge segmented Vela forward-IV results and compare against Sentaurus."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def read_vela(path: Path) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    with path.open(newline="", encoding="utf-8-sig") as stream:
        for row in csv.DictReader(stream):
            if row["converged"] != "1":
                continue
            rows.append(
                {
                    "bias_V": float(row["bias_V"]),
                    "current_A_per_um": float(row["current_total_A_per_um"]),
                }
            )
    return rows


def read_sentaurus(path: Path) -> list[dict[str, float]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return [
            {
                "bias_V": float(row["bias_V"]),
                "current_A_per_um": float(row["current_total"]),
            }
            for row in csv.DictReader(stream)
        ]


def unique_sorted(rows: list[dict[str, float]]) -> list[dict[str, float]]:
    by_bias: dict[float, dict[str, float]] = {}
    for row in rows:
        by_bias[round(row["bias_V"], 12)] = row
    return [by_bias[key] for key in sorted(by_bias)]


def write_csv(path: Path, rows: list[dict[str, float]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=["bias_V", "current_A_per_um"])
        writer.writeheader()
        writer.writerows(rows)


def interpolate_abs(rows: list[dict[str, float]], biases: np.ndarray) -> np.ndarray:
    x = np.asarray([row["bias_V"] for row in rows])
    y = np.abs(np.asarray([row["current_A_per_um"] for row in rows]))
    return np.interp(biases, x, y)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vela", type=Path, nargs="+", required=True)
    parser.add_argument("--sentaurus", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()

    vela = unique_sorted(
        [row for path in args.vela for row in read_vela(path)]
    )
    sentaurus = unique_sorted(read_sentaurus(args.sentaurus))
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "vela_forward_iv_merged_0v20v.csv", vela)

    vx = np.asarray([row["bias_V"] for row in vela])
    vy = np.abs(np.asarray([row["current_A_per_um"] for row in vela]))
    sx = np.asarray([row["bias_V"] for row in sentaurus])
    sy = np.abs(np.asarray([row["current_A_per_um"] for row in sentaurus]))

    compare_biases = np.asarray([0.5, 0.7, 1.0, 2.0, 5.0, 10.0, 15.0, 20.0])
    vela_i = interpolate_abs(vela, compare_biases)
    sentaurus_i = interpolate_abs(sentaurus, compare_biases)
    comparison = [
        {
            "bias_V": float(bias),
            "vela_abs_current_A_per_um": float(vi),
            "sentaurus_abs_current_A_per_um": float(si),
            "vela_over_sentaurus": float(vi / si) if si else None,
            "relative_difference": float((vi - si) / si) if si else None,
        }
        for bias, vi, si in zip(compare_biases, vela_i, sentaurus_i)
    ]
    with (args.out_dir / "forward_iv_selected_bias_comparison.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=comparison[0].keys())
        writer.writeheader()
        writer.writerows(comparison)

    stitch_left = next(
        (row for row in vela if abs(row["bias_V"] - 9.9) < 1e-10), None
    )
    stitch_right = next(
        (row for row in vela if abs(row["bias_V"] - 9.901) < 1e-10), None
    )
    stitch_relative_jump = (
        abs(stitch_right["current_A_per_um"] - stitch_left["current_A_per_um"])
        / max(abs(stitch_left["current_A_per_um"]), 1e-300)
        if stitch_left is not None and stitch_right is not None
        else None
    )
    summary = {
        "vela_points": len(vela),
        "sentaurus_points": len(sentaurus),
        "vela_bias_range_V": [float(vx.min()), float(vx.max())],
        "sentaurus_bias_range_V": [float(sx.min()), float(sx.max())],
        "stitch_9p9_to_9p901_relative_current_jump": stitch_relative_jump,
        "comparison": comparison,
    }
    (args.out_dir / "forward_iv_comparison_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.labelsize": 10,
            "legend.fontsize": 9,
        }
    )
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 5.4), dpi=160)
    vela_color = "#2F6B9A"
    sentaurus_color = "#C47A2C"

    ax = axes[0]
    ax.semilogy(
        sx,
        np.maximum(sy, 1e-30),
        color=sentaurus_color,
        linewidth=2.0,
        label="Sentaurus 2018 legacy IV (low-field)",
    )
    ax.semilogy(
        vx,
        np.maximum(vy, 1e-30),
        color=vela_color,
        linestyle="--",
        linewidth=2.0,
        label="Vela low-field (Masetti)",
    )
    ax.set_title("Forward I-V comparison (log scale)")
    ax.set_xlabel("Anode voltage (V)")
    ax.set_ylabel("|Anode current| (A/um)")
    ax.set_xlim(0, 20)
    ax.set_ylim(1e-20, max(float(vy.max()), float(sy.max())) * 1.8)
    ax.grid(True, which="major", color="#D8D8D8", linewidth=0.7)
    ax.grid(True, which="minor", color="#EEEEEE", linewidth=0.4)
    ax.legend(frameon=False, loc="lower right")

    ax = axes[1]
    ax.plot(
        sx,
        sy * 1e3,
        color=sentaurus_color,
        linewidth=2.0,
        label="Sentaurus 2018 legacy IV (low-field)",
    )
    ax.plot(
        vx,
        vy * 1e3,
        color=vela_color,
        linestyle="--",
        linewidth=2.0,
        label="Vela low-field (Masetti)",
    )
    ax.set_title("Forward I-V comparison (linear scale)")
    ax.set_xlabel("Anode voltage (V)")
    ax.set_ylabel("|Anode current| (mA/um)")
    ax.set_xlim(0, 20)
    ax.set_ylim(bottom=0)
    ax.grid(True, color="#E0E0E0", linewidth=0.7)
    ax.legend(frameon=False, loc="upper left")

    fig.suptitle("PN2D coarse7x3: avalanche off, SRH on", fontsize=14, y=0.98)
    fig.text(
        0.5,
        0.015,
        "Same coarse7x3 mesh; avalanche off and SRH on. Both curves use "
        "doping-dependent low-field mobility without high-field saturation.",
        ha="center",
        color="#555555",
        fontsize=8.5,
    )
    fig.tight_layout(rect=(0, 0.05, 1, 0.94))
    fig.savefig(args.out_dir / "pn2d_coarse7x3_forward_iv_avalanche_off_0v20v.png")
    fig.savefig(args.out_dir / "pn2d_coarse7x3_forward_iv_avalanche_off_0v20v.svg")


if __name__ == "__main__":
    main()
