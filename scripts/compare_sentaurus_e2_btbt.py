#!/usr/bin/env python3
"""Verify the documented Sentaurus E2 law against exported nodal fields."""

from __future__ import annotations

import argparse
import csv
import math
import statistics
from pathlib import Path


DEFAULT_ROOT = Path(
    "build-release/reference_tcad/bvmethods_sentaurus2018/run01/"
    "sentaurus_iic_multibias_exact_20260803/imported"
)


def read_scalar(path: Path) -> dict[int, float]:
    with path.open(newline="", encoding="utf-8") as handle:
        return {
            int(row["node_id"]): float(row["component0"])
            for row in csv.DictReader(handle)
        }


def read_vector_magnitude(path: Path) -> dict[int, float]:
    with path.open(newline="", encoding="utf-8") as handle:
        return {
            int(row["node_id"]): math.hypot(
                float(row["component0"]), float(row["component1"])
            )
            for row in csv.DictReader(handle)
        }


def token(bias: float) -> str:
    return f"iic_v{bias:.6f}".replace(".", "p")


def percentile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    position = fraction * (len(ordered) - 1)
    lo = int(math.floor(position))
    hi = int(math.ceil(position))
    weight = position - lo
    return ordered[lo] * (1.0 - weight) + ordered[hi] * weight


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=DEFAULT_ROOT)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--biases", nargs="+", type=float,
                        default=[1.0, 2.0, 4.0, 5.0, 6.0, 6.4])
    parser.add_argument("--A", type=float, default=3.4e21,
                        help="E2 A in cm^-1 s^-1 V^-2")
    parser.add_argument("--B", type=float, default=22.6e6,
                        help="E2 B in V/cm")
    args = parser.parse_args()

    args.output.parent.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, float]] = []
    with args.output.open("w", newline="", encoding="utf-8") as output:
        fieldnames = [
            "bias_V", "selected_nodes", "median_predicted_over_exported",
            "p05_predicted_over_exported", "p95_predicted_over_exported",
            "median_abs_log10_error", "max_exported_generation_cm3_s",
        ]
        writer = csv.DictWriter(output, fieldnames=fieldnames)
        writer.writeheader()
        for bias in args.biases:
            fields_dir = args.root / token(bias) / "fields"
            field = read_vector_magnitude(fields_dir / "ElectricField_region3.csv")
            exported = read_scalar(fields_dir / "Band2BandGeneration_region3.csv")
            max_exported = max(exported.values())
            threshold = max(1.0e5, max_exported * 1.0e-12)
            ratios: list[float] = []
            log_errors: list[float] = []
            for node, generation in exported.items():
                magnitude = field.get(node, 0.0)
                if generation <= threshold or magnitude <= 0.0:
                    continue
                predicted = args.A * magnitude * magnitude * math.exp(-args.B / magnitude)
                ratio = predicted / generation
                ratios.append(ratio)
                log_errors.append(abs(math.log10(ratio)))
            row = {
                "bias_V": bias,
                "selected_nodes": len(ratios),
                "median_predicted_over_exported": statistics.median(ratios),
                "p05_predicted_over_exported": percentile(ratios, 0.05),
                "p95_predicted_over_exported": percentile(ratios, 0.95),
                "median_abs_log10_error": statistics.median(log_errors),
                "max_exported_generation_cm3_s": max_exported,
            }
            rows.append(row)
            writer.writerow(row)

    for row in rows:
        print(
            f"{row['bias_V']:.3f} V: median E2/export="
            f"{row['median_predicted_over_exported']:.6g}, "
            f"nodes={int(row['selected_nodes'])}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
