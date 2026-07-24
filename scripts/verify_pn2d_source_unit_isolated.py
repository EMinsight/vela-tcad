#!/usr/bin/env python3
"""Independently verify the isolated continuity source-unit audit."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    manifest = json.loads(
        (root / "manifest.json").read_text(encoding="utf-8")
    )
    terms = rows(root / "term_comparison.csv")
    edges = rows(root / "edge_flux_control.csv")
    jacobian = rows(root / "jacobian_comparison.csv")
    updates = rows(root / "first_update_comparison.csv")
    failures: list[str] = []
    expected = {
        "terms": (len(terms), 72),
        "edges": (len(edges), 108),
        "jacobian": (len(jacobian), 30),
        "updates": (len(updates), 48),
    }
    for name, (observed, target) in expected.items():
        if observed != target:
            failures.append(f"{name} count is {observed}, expected {target}")
    maximum_edge = max(
        float(row["relative_difference"]) for row in edges
    )
    maximum_flux = max(
        float(row["flux_relative_difference"]) for row in terms
    )
    srh = [
        float(row["srh_scaled_over_one"])
        for row in terms
        if row["srh_scaled_over_one"] != ""
    ]
    impact = [
        float(row["impact_scaled_over_one"])
        for row in terms
        if row["impact_scaled_over_one"] != ""
    ]
    maximum_srh_error = max(abs(value - 1.0e-8) for value in srh)
    maximum_impact_error = max(
        abs(value - 1.0e-8) for value in impact
    )
    maximum_factor_one_jacobian = max(
        float(row["factor_one_rel_diff"]) for row in jacobian
    )
    maximum_factor_scaled_jacobian = max(
        float(row["factor_scaled_rel_diff"]) for row in jacobian
    )
    first_update_reduced_fraction = sum(
        float(float(row["scaled_over_one_abs_delta"]) < 1.0)
        for row in updates
    ) / len(updates)
    checks = {
        "edge_flux_unchanged": maximum_edge <= 1.0e-14,
        "term_flux_unchanged": maximum_flux <= 1.0e-14,
        "srh_ratio": maximum_srh_error <= 1.0e-14,
        "impact_ratio": maximum_impact_error <= 1.0e-14,
        "factor_scaled_jacobian": maximum_factor_scaled_jacobian <= 1.0e-8,
        "forward_iv_insensitive": float(
            manifest["maximum_forward_iv_relative_change"]
        )
        <= 1.0e-6,
        "dimensional_factor": math.isclose(
            float(manifest["dimensional_source_integral_factor"]),
            1.0e-8,
            rel_tol=0.0,
            abs_tol=1.0e-20,
        ),
    }
    failures.extend(
        name for name, passed in checks.items() if not passed
    )
    if (
        manifest.get("typed_outcome")
        != "source_factor_dimensionally_required_forward_iv_insensitive"
    ):
        failures.append("typed outcome mismatch")
    verification = {
        "status": "passed" if not failures else "failed",
        "failure_count": len(failures),
        "failures": failures,
        "checks": checks,
        "maximum_edge_flux_relative_difference": maximum_edge,
        "maximum_term_flux_relative_difference": maximum_flux,
        "maximum_srh_ratio_error_from_1e_8": maximum_srh_error,
        "maximum_impact_ratio_error_from_1e_8": maximum_impact_error,
        "maximum_factor_one_analytic_fd_jacobian_difference": (
            maximum_factor_one_jacobian
        ),
        "maximum_factor_scaled_analytic_fd_jacobian_difference": (
            maximum_factor_scaled_jacobian
        ),
        "first_update_reduced_fraction": first_update_reduced_fraction,
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
