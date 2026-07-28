#!/usr/bin/env python3
"""Compare PN2D forward IV before and after the 2-D charge-area fix."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def read_curve(path: Path, current_column: str) -> tuple[np.ndarray, np.ndarray]:
    points: dict[float, float] = {}
    with path.open(newline="", encoding="utf-8-sig") as stream:
        for row in csv.DictReader(stream):
            if row.get("converged", "1") != "1":
                continue
            bias = float(row["bias_V"])
            points[round(bias, 12)] = abs(float(row[current_column]))
    biases = np.asarray(sorted(points))
    currents = np.asarray([points[bias] for bias in biases])
    return biases, currents


def exact_sentaurus_current(fields_root: Path, bias: int) -> float:
    path = fields_root / f"{bias}v" / "fields" / "ContactCurrentFlux_region2.csv"
    with path.open(newline="", encoding="utf-8-sig") as stream:
        row = next(csv.DictReader(stream))
    return abs(float(row["component0"]))


def exact_curve_current(
    biases: np.ndarray, currents: np.ndarray, bias: int
) -> float:
    matches = np.flatnonzero(np.isclose(biases, float(bias), atol=1e-10))
    if len(matches) != 1:
        raise ValueError(f"Expected one exact curve point at {bias} V")
    return float(currents[matches[0]])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sentaurus", type=Path, required=True)
    parser.add_argument("--sentaurus-fields", type=Path, required=True)
    parser.add_argument("--before", type=Path, required=True)
    parser.add_argument("--after", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()

    sx, sy = read_curve(args.sentaurus, "current_total")
    bx, by = read_curve(args.before, "current_total_A_per_um")
    ax, ay = read_curve(args.after, "current_total_A_per_um")
    args.out_dir.mkdir(parents=True, exist_ok=True)

    common = np.asarray([1, 2, 5, 10, 15, 20], dtype=float)
    sent = np.asarray(
        [exact_sentaurus_current(args.sentaurus_fields, int(bias)) for bias in common]
    )
    before = np.asarray(
        [exact_curve_current(bx, by, int(bias)) for bias in common]
    )
    after = np.asarray(
        [exact_curve_current(ax, ay, int(bias)) for bias in common]
    )
    before_rel = (before - sent) / sent
    after_rel = (after - sent) / sent

    rows = []
    for bias, s_cur, b_cur, a_cur in zip(common, sent, before, after):
        rows.append(
            {
                "bias_V": bias,
                "sentaurus_A_per_um": s_cur,
                "vela_before_A_per_um": b_cur,
                "vela_after_A_per_um": a_cur,
                "before_over_sentaurus": b_cur / s_cur,
                "after_over_sentaurus": a_cur / s_cur,
                "after_relative_difference_percent": 100.0 * (a_cur - s_cur) / s_cur,
            }
        )
    with (args.out_dir / "charge_area_fix_iv_selected_biases.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)

    summary = {
        "comparison_basis": "exact same-bias Sentaurus TDR contact-current anchors",
        "comparison_biases_V": common.tolist(),
        "comparison_points": len(common),
        "log10_rmse_before_decade": float(
            np.sqrt(np.mean(np.square(np.log10(before / sent))))
        ),
        "log10_rmse_after_decade": float(
            np.sqrt(np.mean(np.square(np.log10(after / sent))))
        ),
        "median_abs_relative_error_before_percent": float(
            100.0 * np.median(np.abs(before_rel))
        ),
        "median_abs_relative_error_after_percent": float(
            100.0 * np.median(np.abs(after_rel))
        ),
        "endpoint_20V": rows[-1],
    }
    (args.out_dir / "charge_area_fix_iv_summary.json").write_text(
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
    fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.4), dpi=170)
    sent_color = "#C56A21"
    before_color = "#7A7A7A"
    after_color = "#246B9C"

    left = axes[0]
    left.semilogy(sx, np.maximum(sy, 1e-30), color=sent_color, linewidth=2.3,
                  label="Sentaurus")
    left.semilogy(bx, np.maximum(by, 1e-30), color=before_color, linewidth=1.7,
                  linestyle=":", label="Vela before fix")
    left.semilogy(ax, np.maximum(ay, 1e-30), color=after_color, linewidth=2.1,
                  linestyle="--", label="Vela after charge-area fix")
    left.set(
        title="Forward I-V",
        xlabel="Anode voltage (V)",
        ylabel="|Anode current| (A/um)",
        xlim=(0, 20),
    )
    left.set_ylim(1e-20, 3e-2)
    left.grid(True, which="major", color="#D8D8D8", linewidth=0.7)
    left.grid(True, which="minor", color="#EEEEEE", linewidth=0.4)
    left.legend(frameon=False, loc="lower right")

    right = axes[1]
    right.plot(common, 100.0 * before_rel, color=before_color, linewidth=1.7,
               linestyle=":", marker="o", label="Before fix")
    right.plot(common, 100.0 * after_rel, color=after_color, linewidth=2.1,
               linestyle="--", marker="s", label="After fix")
    right.axhline(0.0, color="#333333", linewidth=0.8)
    right.set(
        title="Relative current difference at exact TDR anchors",
        xlabel="Anode voltage (V)",
        ylabel="(Vela - Sentaurus) / Sentaurus (%)",
        xlim=(0, 20),
    )
    right.grid(True, color="#E0E0E0", linewidth=0.7)
    right.legend(frameon=False, loc="best")

    fig.suptitle("PN2D coarse7x3: 2-D charge-area scaling correction", fontsize=14)
    fig.text(
        0.5,
        0.012,
        "Avalanche off; SRH on; matched low-field mobility configuration.",
        ha="center",
        color="#555555",
        fontsize=8.5,
    )
    fig.tight_layout(rect=(0, 0.045, 1, 0.94))
    fig.savefig(args.out_dir / "pn2d_charge_area_fix_forward_iv_0v20v.png")
    fig.savefig(args.out_dir / "pn2d_charge_area_fix_forward_iv_0v20v.svg")


if __name__ == "__main__":
    main()
