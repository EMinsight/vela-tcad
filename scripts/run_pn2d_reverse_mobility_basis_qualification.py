#!/usr/bin/env python3
"""Qualify PN2D reverse BV behavior for mobility doping-basis candidates."""

from __future__ import annotations

import argparse
import csv
import json
import math
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


BASES = ("net_doping", "cell_reconstructed_total_impurity")


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def run_case(
    runner: Path,
    base: dict,
    out_dir: Path,
    basis: str,
    avalanche_enabled: bool,
) -> dict[str, Path]:
    tag = f"{basis}_{'on' if avalanche_enabled else 'off'}"
    case_dir = out_dir / tag
    case_dir.mkdir(parents=True, exist_ok=True)
    config = json.loads(json.dumps(base))
    config["solver"]["mobility"]["doping_concentration_basis"] = basis
    if not avalanche_enabled:
        config["solver"]["impact_ionization"] = {"model": "none"}
    output_csv = case_dir / f"{tag}.csv"
    state_csv = case_dir / f"{tag}_last_state.csv"
    newton_csv = case_dir / f"{tag}_newton_history.csv"
    continuity_csv = case_dir / f"{tag}_continuity_balance.csv"
    source_csv = case_dir / f"{tag}_triangle_gss_sources.csv"
    config["output_csv"] = str(output_csv.resolve())
    sweep = config["sweep"]
    sweep["start"] = 0.0
    sweep["stop"] = -20.0
    sweep["step"] = -0.05
    sweep.pop("bias_points", None)
    sweep["write_state_file"] = str(state_csv.resolve())
    sweep["write_vtk"] = False
    sweep["vtk_prefix"] = str((case_dir / "vtk" / tag).resolve())
    sweep["diagnostics"] = {
        "continuity_balance": {
            "enabled": True,
            "contacts": ["Anode", "Cathode"],
            "csv_file": str(continuity_csv.resolve()),
        },
        "newton_history": {
            "enabled": True,
            "csv_file": str(newton_csv.resolve()),
        },
    }
    if avalanche_enabled:
        sweep["diagnostics"]["triangle_gss_sources"] = {
            "enabled": True,
            "csv_file": str(source_csv.resolve()),
        }
    config_path = case_dir / f"simulation_{tag}.json"
    config_path.write_text(json.dumps(config, indent=2), encoding="utf-8")
    completed = subprocess.run(
        [str(runner), "--config", str(config_path)],
        cwd=case_dir,
        text=True,
        capture_output=True,
        check=False,
    )
    (case_dir / "runner.stdout.log").write_text(
        completed.stdout, encoding="utf-8"
    )
    (case_dir / "runner.stderr.log").write_text(
        completed.stderr, encoding="utf-8"
    )
    if completed.returncode:
        raise RuntimeError(
            f"{tag} failed with {completed.returncode}: "
            f"{completed.stderr[-2000:]}"
        )
    return {
        "config": config_path,
        "iv": output_csv,
        "source": source_csv,
        "vtk_dir": case_dir / "vtk",
    }


def current_by_bias(path: Path) -> dict[int, float]:
    result: dict[int, float] = {}
    for row in read_rows(path):
        if row.get("converged") != "1":
            continue
        raw_bias = abs(float(row["bias_V"])); bias = int(round(raw_bias))
        if abs(raw_bias - bias) > 1.0e-8:
            continue
        result[bias] = abs(float(row["current_total_A_per_um"]))
    return result


def sentaurus_integer_curves(
    path: Path,
) -> tuple[dict[int, float], dict[int, float]]:
    on: dict[int, float] = {}
    off: dict[int, float] = {}
    for row in read_rows(path):
        bias = int(round(abs(float(row["bias_V"]))))
        on[bias] = abs(float(row["sentaurus_avalanche_on_A_per_um"]))
        off[bias] = abs(float(row["sentaurus_avalanche_off_A_per_um"]))
    return on, off


def gain_threshold_voltage(
    current_on: dict[int, float],
    current_off: dict[int, float],
    threshold: float,
) -> float | None:
    points = [
        (bias, current_on[bias] / current_off[bias])
        for bias in sorted(current_on)
        if bias > 0 and bias in current_off and current_off[bias] > 0.0
    ]
    for (v0, gain0), (v1, gain1) in zip(points, points[1:]):
        if gain0 < threshold <= gain1:
            log0 = math.log(max(gain0, 1.0e-300))
            log1 = math.log(max(gain1, 1.0e-300))
            target = math.log(threshold)
            fraction = (target - log0) / (log1 - log0)
            return v0 + fraction * (v1 - v0)
    return None


def source_summary(path: Path, basis: str) -> list[dict[str, object]]:
    grouped: dict[int, list[dict[str, str]]] = {}
    for row in read_rows(path):
        raw_bias = abs(float(row["bias_V"])); bias = int(round(raw_bias))
        if abs(raw_bias - bias) > 1.0e-8:
            continue
        grouped.setdefault(bias, []).append(row)
    summaries: list[dict[str, object]] = []
    for bias, rows in sorted(grouped.items()):
        cell_sources: dict[int, float] = {}
        total = 0.0
        weighted_x = 0.0
        weighted_y = 0.0
        for row in rows:
            source = max(0.0, float(row["edge_source_integral"]))
            total += source
            cell = int(row["cell_id"])
            cell_sources[cell] = cell_sources.get(cell, 0.0) + source
            midpoint_x = 0.5 * (float(row["x0_um"]) + float(row["x1_um"]))
            midpoint_y = 0.5 * (float(row["y0_um"]) + float(row["y1_um"]))
            weighted_x += source * midpoint_x
            weighted_y += source * midpoint_y
        peak_cell, peak_source = max(
            cell_sources.items(), key=lambda item: item[1]
        )
        summaries.append(
            {
                "basis": basis,
                "bias_V": -bias,
                "total_source_integral": total,
                "peak_cell_id": peak_cell,
                "peak_cell_source_integral": peak_source,
                "source_centroid_x_um": weighted_x / total if total else 0.0,
                "source_centroid_y_um": weighted_y / total if total else 0.0,
            }
        )
    return summaries


def make_plot(
    path: Path,
    curves: dict[str, tuple[dict[int, float], dict[int, float]]],
) -> None:
    colors = {
        "Sentaurus": "#B279A2",
        "net_doping": "#4C78A8",
        "cell_reconstructed_total_impurity": "#E45756",
    }
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 5.2), dpi=180)
    for name, (on, off) in curves.items():
        biases = np.asarray(sorted(on))
        on_values = np.asarray([on[value] for value in biases])
        gains = np.asarray([on[value] / off[value] for value in biases])
        axes[0].semilogy(
            biases,
            np.maximum(on_values, 1.0e-30),
            marker="o",
            linewidth=1.8,
            markersize=3.5,
            color=colors[name],
            label=name,
        )
        axes[1].semilogy(
            biases,
            np.maximum(gains, 1.0e-8),
            marker="o",
            linewidth=1.8,
            markersize=3.5,
            color=colors[name],
            label=name,
        )
    axes[0].set_title("Reverse current")
    axes[0].set_xlabel("Reverse bias |V| (V)")
    axes[0].set_ylabel("|Anode current| (A/µm)")
    axes[1].set_title("Avalanche multiplication")
    axes[1].set_xlabel("Reverse bias |V| (V)")
    axes[1].set_ylabel("|I(on)| / |I(off)|")
    axes[1].axhline(2.0, color="#777777", linestyle=":", linewidth=1.0)
    for axis in axes:
        axis.set_xlim(1, 20)
        axis.grid(True, which="both", alpha=0.3, linewidth=0.6)
        axis.legend(frameon=False)
    fig.suptitle(
        "PN2D coarse7x3 reverse-BV mobility-basis qualification",
        fontsize=14,
    )
    fig.tight_layout()
    fig.savefig(path)
    fig.savefig(path.with_suffix(".svg"))
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--base-config", type=Path, required=True)
    parser.add_argument("--sentaurus-control", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()

    base = json.loads(args.base_config.read_text(encoding="utf-8-sig"))
    args.out_dir.mkdir(parents=True, exist_ok=True)
    artifacts: dict[str, dict[str, Path]] = {}
    for basis in BASES:
        for enabled in (True, False):
            key = f"{basis}_{'on' if enabled else 'off'}"
            artifacts[key] = run_case(
                args.runner.resolve(),
                base,
                args.out_dir.resolve(),
                basis,
                enabled,
            )

    sent_on, sent_off = sentaurus_integer_curves(args.sentaurus_control)
    comparison_rows: list[dict[str, object]] = []
    curves: dict[str, tuple[dict[int, float], dict[int, float]]] = {
        "Sentaurus": (sent_on, sent_off)
    }
    for basis in BASES:
        on = current_by_bias(artifacts[f"{basis}_on"]["iv"])
        off = current_by_bias(artifacts[f"{basis}_off"]["iv"])
        curves[basis] = (on, off)
        for bias in range(1, 21):
            comparison_rows.append(
                {
                    "basis": basis,
                    "bias_V": -bias,
                    "sentaurus_on_A_per_um": sent_on[bias],
                    "sentaurus_off_A_per_um": sent_off[bias],
                    "sentaurus_gain": sent_on[bias] / sent_off[bias],
                    "vela_on_A_per_um": on[bias],
                    "vela_off_A_per_um": off[bias],
                    "vela_gain": on[bias] / off[bias],
                    "on_current_ratio_to_sentaurus": on[bias] / sent_on[bias],
                    "gain_ratio_to_sentaurus": (
                        (on[bias] / off[bias])
                        / (sent_on[bias] / sent_off[bias])
                    ),
                }
            )
    write_rows(args.out_dir / "reverse_iv_gain_comparison.csv", comparison_rows)

    source_rows: list[dict[str, object]] = []
    for basis in BASES:
        source_rows.extend(
            source_summary(artifacts[f"{basis}_on"]["source"], basis)
        )
    write_rows(args.out_dir / "avalanche_source_summary.csv", source_rows)

    threshold_rows: list[dict[str, object]] = []
    for name, (on, off) in curves.items():
        threshold_rows.append(
            {
                "basis": name,
                "gain_1p5_voltage_V": gain_threshold_voltage(on, off, 1.5),
                "gain_2_voltage_V": gain_threshold_voltage(on, off, 2.0),
                "gain_5_voltage_V": gain_threshold_voltage(on, off, 5.0),
                "gain_20_voltage_V": gain_threshold_voltage(on, off, 20.0),
            }
        )
    write_rows(args.out_dir / "breakdown_gain_thresholds.csv", threshold_rows)
    make_plot(args.out_dir / "reverse_iv_gain_comparison.png", curves)
    (args.out_dir / "qualification_manifest.json").write_text(
        json.dumps(
            {
                "base_config": str(args.base_config.resolve()),
                "sentaurus_control": str(args.sentaurus_control.resolve()),
                "artifacts": {
                    key: {name: str(path.resolve()) for name, path in value.items()}
                    for key, value in artifacts.items()
                },
            },
            indent=2,
        ),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
