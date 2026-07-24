#!/usr/bin/env python3
"""Independently verify the high-field residual/first-step follow-up."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def update_key(row: dict[str, str]) -> tuple[str, str, str, str, str]:
    return (
        row["topology"],
        row["bias_V"],
        row["mode"],
        row["carrier"],
        row["node_id"],
    )


def close(left: float, right: float, tolerance: float = 1.0e-14) -> bool:
    return abs(left - right) <= tolerance * max(1.0, abs(left), abs(right))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    manifest = json.loads(
        (root / "manifest.json").read_text(encoding="utf-8")
    )
    residual = rows(root / "residual_waterfall.csv")
    norms = rows(root / "state_residual_norms.csv")
    updates = rows(root / "first_update.csv")
    jacobian = rows(root / "jacobian_audit.csv")
    failures: list[str] = []
    expected = {
        "residual": (len(residual), 800),
        "norms": (len(norms), 400),
        "updates": (len(updates), 144),
        "jacobian": (len(jacobian), 90),
    }
    for name, (observed, target) in expected.items():
        if observed != target:
            failures.append(f"{name} count is {observed}, expected {target}")

    production = [
        row for row in residual if row["branch"] == "vela_production"
    ]
    closure = max(
        abs(
            float(row["final_residual_normalized_units"])
            - float(row["recorded_production_residual"])
        )
        / max(
            1.0,
            abs(float(row["final_residual_normalized_units"])),
            abs(float(row["recorded_production_residual"])),
        )
        for row in production
    )
    if closure > 1.0e-12:
        failures.append("production edge-to-node replay does not close")
    if not close(
        closure,
        float(manifest["maximum_production_edge_to_node_closure"]),
    ):
        failures.append("production closure manifest mismatch")
    if {
        row["source_unit_policy"] for row in residual
    } != {"unchanged_production_snapshot"}:
        failures.append("mobility branches changed source-unit policy")

    element_electric = [
        row
        for row in norms
        if row["branch"]
        == "sentaurus_lowfield_element_electric_field"
    ]
    element_residual_all_improved = all(
        float(row["ratio_to_production"]) < 1.0
        for row in element_electric
    )
    if len(element_electric) != 80 or not element_residual_all_improved:
        failures.append(
            "element-electric residual does not improve all 80 pairs"
        )

    baseline = {
        update_key(row): row
        for row in updates
        if row["branch"] == "vela_global_qfp_config"
    }
    electric = {
        update_key(row): row
        for row in updates
        if row["branch"] == "vela_global_electric_field_config"
    }
    if set(baseline) != set(electric) or len(electric) != 48:
        failures.append("first-update paired lattice is incomplete")
    ratios = {
        key: float(electric[key]["absolute_delta_qfp_V"])
        / max(float(baseline[key]["absolute_delta_qfp_V"]), 1.0e-300)
        for key in set(baseline) & set(electric)
    }
    first_update_all_improved = all(value < 1.0 for value in ratios.values())
    if not first_update_all_improved:
        failures.append("electric first update does not improve every pair")

    max_jacobian = max(float(row["rel_diff"]) for row in jacobian)
    if max_jacobian > 1.0e-8:
        failures.append("analytic/FD Jacobian gate failed")
    if not close(
        max_jacobian,
        float(manifest["maximum_jacobian_relative_difference"]),
    ):
        failures.append("Jacobian manifest maximum mismatch")

    contacts = {0, 2, 3, 4}
    maximum_boundary_difference = 0.0
    for topology in ("mirror", "sketch"):
        for bias in (1, 10, 20):
            work = root / "raw" / topology / f"m{bias}V"
            for suffix in ("carrier", "coupled"):
                qfp = {
                    int(row["node_id"]): row
                    for row in rows(
                        work / f"vela_global_qfp_config_{suffix}.csv"
                    )
                }
                electric_rows = {
                    int(row["node_id"]): row
                    for row in rows(
                        work
                        / f"vela_global_electric_field_config_{suffix}.csv"
                    )
                }
                for node in contacts:
                    for field in (
                        "phin_residual",
                        "phip_residual",
                        "trial_phin_residual",
                        "trial_phip_residual",
                    ):
                        difference = abs(
                            float(qfp[node][field])
                            - float(electric_rows[node][field])
                        )
                        maximum_boundary_difference = max(
                            maximum_boundary_difference, difference
                        )
    if maximum_boundary_difference > 1.0e-14:
        failures.append("contact boundary rows changed")

    if manifest.get("typed_outcome") != "mobility_candidate_causal":
        failures.append("typed outcome is not causal")
    if manifest.get("production_formula_modified") is not False:
        failures.append("follow-up modified a production formula")

    verification = {
        "status": "passed" if not failures else "failed",
        "failure_count": len(failures),
        "failures": failures,
        "production_edge_to_node_closure": closure,
        "element_electric_residual_all_80_pairs_improved": (
            element_residual_all_improved
        ),
        "electric_first_update_all_48_pairs_improved": (
            first_update_all_improved
        ),
        "first_update_ratio_minimum": min(ratios.values()),
        "first_update_ratio_maximum": max(ratios.values()),
        "maximum_jacobian_relative_difference": max_jacobian,
        "maximum_contact_boundary_difference": (
            maximum_boundary_difference
        ),
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
