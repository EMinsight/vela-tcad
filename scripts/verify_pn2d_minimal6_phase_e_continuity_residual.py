#!/usr/bin/env python3
"""Independently verify a PN2D Minimal6 Phase E evidence root."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("root")
    args = parser.parse_args()
    root = Path(args.root).resolve()
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    failures: list[str] = []
    for name, expected in manifest["outputs"].items():
        path = root / name
        if not path.is_file():
            failures.append(f"missing output {name}")
        elif sha256(path) != expected:
            failures.append(f"hash mismatch {name}")

    waterfall = rows(root / "residual_waterfall.csv")
    boundary = rows(root / "boundary_audit.csv")
    jacobian = rows(root / "jacobian_audit.csv")
    updates = rows(root / "first_update.csv")
    states = rows(root / "state_metrics.csv")
    source_controls = rows(root / "source_unit_scaling_control.csv")
    controls = rows(root / "controlled_branch_experiments.csv")
    expected_counts = {
        "waterfall": (len(waterfall), 480),
        "boundary": (len(boundary), 160),
        "jacobian": (len(jacobian), 200),
        "updates": (len(updates), 80),
        "states": (len(states), 40),
        "source_controls": (len(source_controls), 160),
    }
    for label, (actual, expected) in expected_counts.items():
        if actual != expected:
            failures.append(f"{label} count {actual} != {expected}")
    if len(controls) != manifest["contracts"]["controlled_experiment_count"]:
        failures.append("controlled experiment count mismatch")

    keys = {
        (row["topology"], row["bias_V"], row["branch"], row["carrier"], row["node_id"])
        for row in waterfall
    }
    if len(keys) != len(waterfall):
        failures.append("duplicate waterfall identity")
    for row in waterfall:
        terms = [
            float(row["sg_divergence_normalized"]),
            float(row["srh_normalized"]),
            float(row["impact_normalized"]),
            float(row["gauge_normalized"]),
            float(row["contact_boundary_flux_normalized"]),
        ]
        final = float(row["final_residual_normalized_units"])
        scale = max(1.0, abs(final), *(abs(value) for value in terms))
        if abs(final - sum(terms)) > 2.0e-13 * scale:
            failures.append("waterfall algebraic closure failed")
            break
        physical = float(row["final_unconverted_SI_equivalent_per_m_s"])
        expected = final * float(row["continuity_scale_SI_per_m_s"])
        if abs(physical - expected) > 2.0e-13 * max(1.0, abs(expected)):
            failures.append("SI residual scaling closure failed")
            break

    for row in source_controls:
        factor = float(row["source_unit_factor"])
        expected = (
            float(row["sg_divergence_normalized"])
            + factor * float(row["original_srh_normalized"])
            + factor * float(row["original_impact_normalized"])
        )
        actual = float(row["source_corrected_residual_normalized_units"])
        if abs(actual - expected) > 2.0e-13 * max(1.0, abs(expected)):
            failures.append("source-unit control closure failed")
            break
        if factor != 1.0e-8:
            failures.append("unexpected source-unit factor")
            break

    if max(float(row["electron_density_abs_dex"]) for row in boundary) > 1.0e-4:
        failures.append("electron boundary density gate failed")
    if max(float(row["hole_density_abs_dex"]) for row in boundary) > 1.0e-4:
        failures.append("hole boundary density gate failed")
    if max(
        max(abs(float(row["phin_error_V"])), abs(float(row["phip_error_V"])))
        for row in boundary
    ) > 1.0e-10:
        failures.append("contact QFP gate failed")
    if not all(math.isfinite(float(row["rel_diff"])) for row in jacobian):
        failures.append("non-finite Jacobian audit")
    if not all(
        math.isfinite(float(row["max_abs_production_residual_normalized_units"]))
        for row in states
    ):
        failures.append("non-finite state metric")

    result = {
        "status": "passed" if not failures else "failed",
        "failures": failures,
        "verified_output_count": len(manifest["outputs"]),
        "waterfall_row_count": len(waterfall),
        "state_count": len(states),
        "controlled_experiment_count": len(controls),
    }
    (root / "independent_verification.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8", newline="\n"
    )
    print(json.dumps(result, sort_keys=True))
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())

