#!/usr/bin/env python3
"""Summarize corrected SLOT-LDMOS IALMob voltage branches and BVDS."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any


CASES = ("ialmob_off", "ialmob_on")


class BranchAnalysisError(RuntimeError):
    """Raised when an IV branch cannot be interpreted safely."""


def read_converged_points(path: Path) -> list[dict[str, float]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    points: list[dict[str, float]] = []
    for row in rows:
        if row.get("converged") not in {"1", "true", "True"}:
            continue
        voltage_text = row.get("inner_voltage_V") or row.get("bias_V")
        current_text = row.get("current_total_A_per_um")
        if not voltage_text or not current_text:
            raise BranchAnalysisError(f"missing voltage/current columns in {path}")
        voltage = float(voltage_text)
        current = abs(float(current_text))
        if not math.isfinite(voltage) or not math.isfinite(current):
            raise BranchAnalysisError(f"non-finite IV point in {path}")
        points.append({"voltage_V": voltage, "current_A_per_um": current})
    if not points:
        raise BranchAnalysisError(f"no converged IV points in {path}")
    points.sort(key=lambda point: point["voltage_V"])
    return points


def locate_bvds(
    points: list[dict[str, float]], criterion_A_per_um: float
) -> dict[str, Any]:
    if criterion_A_per_um <= 0.0:
        raise BranchAnalysisError("BVDS current criterion must be positive")
    for index, upper in enumerate(points):
        if upper["current_A_per_um"] < criterion_A_per_um:
            continue
        if index == 0:
            return {
                "status": "criterion_at_or_below_first_point",
                "bvds_V": upper["voltage_V"],
                "bracket": [upper, upper],
            }
        lower = points[index - 1]
        if lower["current_A_per_um"] <= 0.0:
            fraction = (
                (criterion_A_per_um - lower["current_A_per_um"])
                / (upper["current_A_per_um"] - lower["current_A_per_um"])
            )
            method = "linear_current"
        else:
            log_lower = math.log(lower["current_A_per_um"])
            log_upper = math.log(upper["current_A_per_um"])
            fraction = (
                (math.log(criterion_A_per_um) - log_lower)
                / (log_upper - log_lower)
            )
            method = "log_current"
        bvds = lower["voltage_V"] + fraction * (
            upper["voltage_V"] - lower["voltage_V"]
        )
        return {
            "status": "located",
            "bvds_V": bvds,
            "interpolation": method,
            "bracket": [lower, upper],
        }
    return {
        "status": "criterion_not_reached",
        "last_point": points[-1],
        "maximum_current_A_per_um": max(
            point["current_A_per_um"] for point in points
        ),
    }


def analyze(
    root: Path,
    output_dir: Path,
    criterion_A_per_um: float = 1.0e-7,
) -> dict[str, Any]:
    cases: dict[str, Any] = {}
    for case in CASES:
        points = read_converged_points(root / case / "iv.csv")
        cases[case] = {
            "point_count": len(points),
            "first_point": points[0],
            "last_point": points[-1],
            "bvds": locate_bvds(points, criterion_A_per_um),
        }
    off_bv = cases["ialmob_off"]["bvds"]
    on_bv = cases["ialmob_on"]["bvds"]
    both_located = off_bv["status"] == "located" and on_bv["status"] == "located"
    summary = {
        "schema": "vela.slot_ldmos.corrected_ialmob_branch_result.v1",
        "criterion_A_per_um": criterion_A_per_um,
        "controlled_delta": "Enhanced Lombardi Enormal mobility only",
        "cases": cases,
        "ialmob_bvds_shift_V": (
            on_bv["bvds_V"] - off_bv["bvds_V"] if both_located else None
        ),
        "requires_post_fold_continuation": not both_located,
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "corrected_ialmob_branch_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    lines = [
        "# Corrected SLOT-LDMOS IALMob branches",
        "",
        "| Case | Points | Last V (V) | Last Id (A/um) | BVDS status | BVDS (V) |",
        "|---|---:|---:|---:|---|---:|",
    ]
    for case in CASES:
        metrics = cases[case]
        result = metrics["bvds"]
        lines.append(
            f"| {case} | {metrics['point_count']} | "
            f"{metrics['last_point']['voltage_V']:.9g} | "
            f"{metrics['last_point']['current_A_per_um']:.9e} | "
            f"{result['status']} | "
            f"{result.get('bvds_V', float('nan')):.9g} |"
        )
    lines.extend([
        "",
        (
            f"IALMob BVDS shift (on - off): {summary['ialmob_bvds_shift_V']:.9g} V"
            if both_located
            else "At least one fixed-voltage branch has not bracketed the current criterion."
        ),
        "",
    ])
    (output_dir / "corrected_ialmob_branch_summary.md").write_text(
        "\n".join(lines), encoding="utf-8"
    )
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--criterion", type=float, default=1.0e-7)
    args = parser.parse_args()
    print(json.dumps(
        analyze(args.root, args.output_dir, args.criterion),
        indent=2,
        sort_keys=True,
    ))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
