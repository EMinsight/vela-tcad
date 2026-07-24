#!/usr/bin/env python3
"""Independently verify the Minimal6 high-field drive evidence."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from pathlib import Path


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def stats(values: list[float]) -> dict[str, float]:
    ordered = sorted(values)
    return {
        "median": statistics.median(values),
        "p95": ordered[math.ceil(0.95 * len(ordered)) - 1],
        "maximum": ordered[-1],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    manifest = json.loads(
        (root / "manifest.json").read_text(encoding="utf-8")
    )
    samples = rows(root / "highfield_drive_samples.csv")
    summary = {
        row["carrier"]: row
        for row in rows(root / "highfield_drive_summary.csv")
    }
    failures: list[str] = []
    if manifest.get("status") != "valid":
        failures.append("manifest status is not valid")
    if len(samples) != 320:
        failures.append(f"sample count is {len(samples)}, expected 320")
    gates = {
        "median": 0.005,
        "p95": 0.03,
        "maximum": 0.05,
    }
    recomputed: dict[str, dict[str, float]] = {}
    for carrier in ("electron", "hole"):
        carrier_rows = [
            row for row in samples if row["carrier"] == carrier
        ]
        if len(carrier_rows) != 160:
            failures.append(
                f"{carrier} sample count is {len(carrier_rows)}, expected 160"
            )
            continue
        result = stats(
            [
                float(row["electric_replay_abs_error_dex"])
                for row in carrier_rows
            ]
        )
        recomputed[carrier] = result
        recorded = summary.get(carrier)
        if recorded is None:
            failures.append(f"missing {carrier} summary")
            continue
        for name, threshold in gates.items():
            key = f"electric_replay_{name}_dex"
            difference = abs(float(recorded[key]) - result[name])
            if difference > 1.0e-15:
                failures.append(
                    f"{carrier} {name} summary mismatch: {difference}"
                )
            if result[name] > threshold:
                failures.append(
                    f"{carrier} {name} {result[name]} exceeds {threshold}"
                )
        native = stats(
            [
                float(row["native_qf_replay_abs_error_dex"])
                for row in carrier_rows
            ]
        )
        triangle = stats(
            [
                float(row["triangle_qf_replay_abs_error_dex"])
                for row in carrier_rows
            ]
        )
        if not (
            result["median"] < native["median"]
            and result["p95"] < native["p95"]
            and result["median"] < triangle["median"]
            and result["p95"] < triangle["p95"]
        ):
            failures.append(
                f"{carrier} electric-field branch is not best on median/P95"
            )
    verification = {
        "status": "passed" if not failures else "failed",
        "failure_count": len(failures),
        "failures": failures,
        "sample_count": len(samples),
        "recomputed_electric_replay": recomputed,
    }
    (root / "independent_verification.json").write_text(
        json.dumps(verification, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    print(json.dumps(verification, indent=2))
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
