#!/usr/bin/env python3
"""Compare scheduled Vela BV controls with Sentaurus on/off terminal currents."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def log_interpolate(rows: list[dict[str, str]], bias: float) -> float:
    points = sorted(
        (abs(float(row["bias_V"])), abs(float(row["current_total_A_per_um"])))
        for row in rows
        if row["converged"] == "1"
    )
    for (x0, y0), (x1, y1) in zip(points, points[1:]):
        if x0 <= bias <= x1:
            if abs(x1 - x0) < 1.0e-15:
                return y1
            weight = (bias - x0) / (x1 - x0)
            log0 = math.log(max(y0, 1.0e-300))
            log1 = math.log(max(y1, 1.0e-300))
            return math.exp(log0 + weight * (log1 - log0))
    if abs(points[-1][0] - bias) < 1.0e-12:
        return points[-1][1]
    raise ValueError(f"bias {bias} is outside Vela curve")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vela-dir", type=Path, required=True)
    parser.add_argument("--sentaurus-control", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    sent_rows = read_csv(args.sentaurus_control)
    sent = {
        int(round(abs(float(row["bias_V"])))): {
            "on": abs(float(row["sentaurus_avalanche_on_A_per_um"])),
            "off": abs(float(row["sentaurus_avalanche_off_A_per_um"])),
        }
        for row in sent_rows
    }
    curves: dict[str, dict[str, list[dict[str, str]]]] = {}
    for basis in ("net_doping", "cell_reconstructed_total_impurity"):
        curves[basis] = {}
        for impact in ("on", "off"):
            path = args.vela_dir / f"{basis}_{impact}" / f"{basis}_{impact}.csv"
            curves[basis][impact] = read_csv(path)

    comparison: list[dict[str, object]] = []
    metrics: dict[str, object] = {}
    for basis, impact_curves in curves.items():
        squared_on: list[float] = []
        squared_off: list[float] = []
        for bias in range(1, 21):
            vela_on = log_interpolate(impact_curves["on"], bias)
            vela_off = log_interpolate(impact_curves["off"], bias)
            sent_on = sent[bias]["on"]
            sent_off = sent[bias]["off"]
            squared_on.append(math.log10(vela_on / sent_on) ** 2)
            squared_off.append(math.log10(vela_off / sent_off) ** 2)
            comparison.append(
                {
                    "basis": basis,
                    "bias_V": -bias,
                    "vela_on_A_per_um": vela_on,
                    "vela_off_A_per_um": vela_off,
                    "vela_gain": vela_on / vela_off,
                    "sentaurus_on_A_per_um": sent_on,
                    "sentaurus_off_A_per_um": sent_off,
                    "sentaurus_gain": sent_on / sent_off,
                    "on_ratio_vela_to_sentaurus": vela_on / sent_on,
                    "off_ratio_vela_to_sentaurus": vela_off / sent_off,
                    "gain_ratio_vela_to_sentaurus": (vela_on / vela_off) / (sent_on / sent_off),
                }
            )
        metrics[basis] = {
            "log10_rmse_on": math.sqrt(sum(squared_on) / len(squared_on)),
            "log10_rmse_off": math.sqrt(sum(squared_off) / len(squared_off)),
            "minus20": comparison[-1],
        }

    comparison_path = args.out_dir / "integer_bias_comparison.csv"
    with comparison_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=comparison[0].keys())
        writer.writeheader()
        writer.writerows(comparison)
    (args.out_dir / "metrics.json").write_text(
        json.dumps(metrics, indent=2), encoding="utf-8"
    )

    fig, axes = plt.subplots(1, 2, figsize=(12.4, 5.0), dpi=180)
    sent_biases = sorted(sent)
    axes[0].semilogy(
        sent_biases, [sent[b]["on"] for b in sent_biases],
        "o-", label="Sentaurus avalanche on", color="#7A5195"
    )
    axes[0].semilogy(
        sent_biases, [sent[b]["off"] for b in sent_biases],
        "o--", label="Sentaurus avalanche off", color="#EF5675"
    )
    colors = {"net_doping": "#2F4B7C", "cell_reconstructed_total_impurity": "#FFA600"}
    for basis, impact_curves in curves.items():
        biases = list(range(1, 21))
        on_values = [log_interpolate(impact_curves["on"], b) for b in biases]
        off_values = [log_interpolate(impact_curves["off"], b) for b in biases]
        axes[0].semilogy(biases, on_values, "-", label=f"Vela {basis} on", color=colors[basis])
        axes[0].semilogy(biases, off_values, "--", label=f"Vela {basis} off", color=colors[basis])
        axes[1].semilogy(
            biases,
            [on / off for on, off in zip(on_values, off_values)],
            "-",
            label=f"Vela {basis}",
            color=colors[basis],
        )
    axes[1].semilogy(
        sent_biases,
        [sent[b]["on"] / sent[b]["off"] for b in sent_biases],
        "o-",
        label="Sentaurus",
        color="#7A5195",
    )
    axes[0].set_ylabel("|Anode current| (A/um)")
    axes[1].set_ylabel("Avalanche gain |I(on)| / |I(off)|")
    for axis in axes:
        axis.set_xlabel("Reverse bias magnitude (V)")
        axis.set_xlim(1, 20)
        axis.grid(True, which="both", alpha=0.3)
        axis.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(args.out_dir / "reverse_bv_comparison.png")
    fig.savefig(args.out_dir / "reverse_bv_comparison.svg")
    plt.close(fig)


if __name__ == "__main__":
    main()
