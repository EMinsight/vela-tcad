#!/usr/bin/env python3
"""Independently verify the Minimal6 high-field box-current evidence."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from pathlib import Path


BRANCH = "sentaurus_lowfield_element_electric_field"


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def quantile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * fraction
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    return ordered[lower] + (
        ordered[upper] - ordered[lower]
    ) * (position - lower)


def close(left: float, right: float, tolerance: float = 1.0e-15) -> bool:
    return abs(left - right) <= tolerance * max(1.0, abs(left), abs(right))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    manifest = json.loads(
        (root / "manifest.json").read_text(encoding="utf-8")
    )
    element = rows(root / "element_edge_samples.csv")
    edge = rows(root / "box_edge_samples.csv")
    summary = rows(root / "box_edge_summary.csv")
    terminal = rows(root / "terminal_comparison.csv")
    kcl = rows(root / "total_current_kcl.csv")
    failures: list[str] = []

    if manifest.get("status") != "bounded_gate_failure":
        failures.append("manifest does not preserve bounded KCL failure")
    if len(element) != 4800:
        failures.append(f"element sample count is {len(element)}, expected 4800")
    if len(edge) != 4320:
        failures.append(f"edge sample count is {len(edge)}, expected 4320")

    recomputed: dict[str, dict[str, float]] = {}
    for carrier in ("electron", "hole"):
        selected = [
            row
            for row in edge
            if row["branch"] == BRANCH
            and row["carrier"] == carrier
            and row["status"] == "valid"
        ]
        if len(selected) != 200:
            failures.append(
                f"{carrier} valid edge count is {len(selected)}, expected 200"
            )
            continue
        errors: list[float] = []
        signs: list[float] = []
        for row in selected:
            reference = float(row["reference_A_per_um"])
            candidate = float(row["candidate_A_per_um"])
            reference_mobility = float(
                row["reference_mobility_m2_per_Vs"]
            )
            candidate_mobility = float(
                row["candidate_mobility_m2_per_Vs"]
            )
            replay = reference * candidate_mobility / reference_mobility
            if not close(candidate, replay):
                failures.append(
                    f"{carrier} edge {row['edge_id']} current scaling mismatch"
                )
                break
            errors.append(abs(math.log10(abs(candidate) / abs(reference))))
            signs.append(
                float(
                    math.copysign(1.0, candidate)
                    == math.copysign(1.0, reference)
                )
            )
        result = {
            "median_abs_dex": statistics.median(errors),
            "p95_abs_dex": quantile(errors, 0.95),
            "maximum_abs_dex": max(errors),
            "sign_agreement_fraction": statistics.mean(signs),
        }
        recomputed[carrier] = result
        recorded = next(
            (
                row
                for row in summary
                if row["branch"] == BRANCH
                and row["carrier"] == carrier
            ),
            None,
        )
        if recorded is None:
            failures.append(f"missing {carrier} summary")
            continue
        for name, value in result.items():
            if not close(float(recorded[name]), value):
                failures.append(f"{carrier} {name} summary mismatch")
        if result["median_abs_dex"] > 0.01:
            failures.append(f"{carrier} median edge gate failed")
        if result["p95_abs_dex"] > 0.05:
            failures.append(f"{carrier} P95 edge gate failed")
        if result["sign_agreement_fraction"] != 1.0:
            failures.append(f"{carrier} sign gate failed")

    branch_terminal = [
        row for row in terminal if row["branch"] == BRANCH
    ]
    maximum_terminal = max(
        float(row["relative_error"]) for row in branch_terminal
    )
    if maximum_terminal > 0.02:
        failures.append("terminal-current gate failed")
    if not close(
        maximum_terminal,
        float(manifest["maximum_terminal_relative_error"]),
    ):
        failures.append("terminal-current manifest maximum mismatch")

    branch_kcl = [row for row in kcl if row["branch"] == BRANCH]
    maximum_kcl = max(
        float(row["relative_to_terminal_total"]) for row in branch_kcl
    )
    if maximum_kcl <= 1.0e-8:
        failures.append("expected fixed-state KCL gate failure was not retained")
    if maximum_kcl > 1.0e-3:
        failures.append("fixed-state KCL perturbation exceeds bounded limit")
    if not close(
        maximum_kcl,
        float(manifest["maximum_kcl_relative_to_terminal"]),
    ):
        failures.append("KCL manifest maximum mismatch")
    expected_gates = {
        "edge_count": True,
        "edge_median": True,
        "edge_p95": True,
        "edge_sign": True,
        "terminal": True,
        "kcl": False,
    }
    if manifest.get("gates") != expected_gates:
        failures.append("manifest gate vector differs from expected result")

    verification = {
        "status": (
            "passed_expected_bounded_gate_failure"
            if not failures
            else "failed"
        ),
        "failure_count": len(failures),
        "failures": failures,
        "recomputed_edge_summary": recomputed,
        "maximum_terminal_relative_error": maximum_terminal,
        "maximum_kcl_relative_to_terminal": maximum_kcl,
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
