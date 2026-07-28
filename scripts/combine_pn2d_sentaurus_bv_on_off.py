#!/usr/bin/env python3
"""Combine refreshed Sentaurus BV on/off exports at integer reverse biases."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def read_curve(path: Path) -> list[tuple[float, float]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.DictReader(stream))
    return sorted(
        (abs(float(row["bias_V"])), abs(float(row["current_total"])))
        for row in rows
    )


def log_interpolate(points: list[tuple[float, float]], bias: float) -> float:
    for (x0, y0), (x1, y1) in zip(points, points[1:]):
        if x0 <= bias <= x1:
            if abs(x1 - x0) < 1.0e-15:
                return y1
            weight = (bias - x0) / (x1 - x0)
            return math.exp(
                math.log(max(y0, 1.0e-300))
                + weight * (math.log(max(y1, 1.0e-300)) - math.log(max(y0, 1.0e-300)))
            )
    raise ValueError(f"bias {bias} is outside curve")


parser = argparse.ArgumentParser()
parser.add_argument("--on", type=Path, required=True)
parser.add_argument("--off", type=Path, required=True)
parser.add_argument("--out", type=Path, required=True)
args = parser.parse_args()
on = read_curve(args.on)
off = read_curve(args.off)
rows = []
for bias in range(1, 21):
    current_on = log_interpolate(on, bias)
    current_off = log_interpolate(off, bias)
    rows.append(
        {
            "bias_V": -bias,
            "sentaurus_avalanche_on_A_per_um": current_on,
            "sentaurus_avalanche_off_A_per_um": current_off,
            "sentaurus_gain_on_over_off": current_on / current_off,
        }
    )
args.out.parent.mkdir(parents=True, exist_ok=True)
with args.out.open("w", newline="", encoding="utf-8") as stream:
    writer = csv.DictWriter(stream, fieldnames=rows[0].keys())
    writer.writeheader()
    writer.writerows(rows)
