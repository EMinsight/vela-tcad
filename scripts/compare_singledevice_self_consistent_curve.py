#!/usr/bin/env python3
"""Compare a descending self-consistent Vela SingleDevice sweep to Sentaurus."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


def rows(path: Path) -> list[tuple[float, float]]:
    with path.open(newline="") as handle:
        parsed = list(csv.DictReader(handle))
    bias_key = "bias_V" if "bias_V" in parsed[0] else "gate_voltage_V"
    current_key = next(
        key for key in (
            "current_total_A_per_um", "current_total", "drain_current_A_per_um")
        if key in parsed[0]
    )
    return [(float(row[bias_key]), abs(float(row[current_key]))) for row in parsed]


def log_interpolate(points: list[tuple[float, float]], bias: float) -> float:
    ordered = sorted(points)
    for candidate_bias, current in ordered:
        if math.isclose(candidate_bias, bias, rel_tol=0.0, abs_tol=1.0e-10):
            return current
    for (left_bias, left_current), (right_bias, right_current) in zip(ordered, ordered[1:]):
        if left_bias <= bias <= right_bias:
            weight = (bias - left_bias) / (right_bias - left_bias)
            return math.exp(
                math.log(max(left_current, 1.0e-300)) * (1.0 - weight)
                + math.log(max(right_current, 1.0e-300)) * weight
            )
    raise ValueError(f"bias {bias} is outside the Vela sweep")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vela", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    args = parser.parse_args()

    vela = rows(args.vela)
    reference = rows(args.reference)
    output: list[dict[str, float]] = []
    for bias, expected in reference:
        actual = log_interpolate(vela, bias)
        output.append({
            "gate_voltage_V": bias,
            "sentaurus_current_A_per_um": expected,
            "vela_current_A_per_um": actual,
            "relative_error": abs(actual - expected) / max(expected, 1.0e-300),
            "orders_of_magnitude_error": abs(
                math.log10(max(actual, 1.0e-300) / max(expected, 1.0e-300))),
        })

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output[0]))
        writer.writeheader()
        writer.writerows(output)
    summary = {
        "points": len(output),
        "max_relative_error": max(row["relative_error"] for row in output),
        "max_orders_of_magnitude": max(
            row["orders_of_magnitude_error"] for row in output),
        "trend_match": all(
            output[i]["vela_current_A_per_um"] <= output[i + 1]["vela_current_A_per_um"]
            for i in range(len(output) - 1)
        ),
    }
    args.summary.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
