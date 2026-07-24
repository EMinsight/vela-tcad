#!/usr/bin/env python3
"""Independently verify Minimal6 current-factor attribution."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


TARGET = (
    "sentaurus",
    "native_element_electric_field",
    "native_elements",
)
BASELINE = ("vela", "global_edge_qfp", "global_edge")


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def close(left: float, right: float, tolerance: float = 2.0e-14) -> bool:
    return abs(left - right) <= tolerance * max(1.0, abs(left), abs(right))


def key(row: dict[str, str]) -> tuple[str, float, str, int]:
    return (
        row["topology"],
        float(row["bias_V"]),
        row["carrier"],
        int(row["edge_id"]),
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    manifest = json.loads(
        (root / "manifest.json").read_text(encoding="utf-8")
    )
    factorial = rows(root / "factorial_samples.csv")
    orders = rows(root / "replacement_orders.csv")
    shapley = rows(root / "shapley_samples.csv")
    central = rows(root / "central_tail.csv")
    source_box = rows(Path(manifest["inputs"]["box_current"]["path"]))
    source_target = {
        key(row): row
        for row in source_box
        if row["branch"]
        == "sentaurus_lowfield_element_electric_field"
        and row["status"] == "valid"
    }
    source_baseline = {
        key(row): row
        for row in source_box
        if row["branch"] == "vela_imported_state_production_mobility"
        and row["status"] == "valid"
    }
    failures: list[str] = []
    if manifest.get("status") != "valid":
        failures.append("manifest status is not valid")
    expected_counts = {
        "factorial": (len(factorial), 4800),
        "orders": (len(orders), 7200),
        "shapley": (len(shapley), 400),
        "central": (len(central), 480),
    }
    for name, (observed, expected) in expected_counts.items():
        if observed != expected:
            failures.append(
                f"{name} count is {observed}, expected {expected}"
            )

    lattice: dict[
        tuple[str, float, str, int],
        dict[tuple[str, str, str], dict[str, str]],
    ] = {}
    for row in factorial:
        lattice.setdefault(key(row), {})[
            (
                row["low_field_source"],
                row["drive"],
                row["support"],
            )
        ] = row
    if any(len(value) != 12 for value in lattice.values()):
        failures.append("at least one active edge lacks the 2x3x2 lattice")

    maximum_closure = 0.0
    for row in shapley:
        local_key = key(row)
        closure = (
            float(row["target_signed_dex"])
            - float(row["baseline_signed_dex"])
            - float(row["low_field_shapley_dex"])
            - float(row["drive_shapley_dex"])
            - float(row["support_shapley_dex"])
        )
        maximum_closure = max(maximum_closure, abs(closure))
        if abs(closure) > 1.0e-12:
            failures.append(f"Shapley closure failed at {local_key}")
            break
        target = lattice[local_key][TARGET]
        baseline = lattice[local_key][BASELINE]
        recorded_target = source_target.get(local_key)
        recorded_baseline = source_baseline.get(local_key)
        if recorded_target is None or recorded_baseline is None:
            failures.append(f"source branch missing at {local_key}")
            break
        if not close(
            float(target["candidate_A_per_um"]),
            float(recorded_target["candidate_A_per_um"]),
        ):
            failures.append(f"target branch mismatch at {local_key}")
            break
        if not close(
            float(baseline["candidate_A_per_um"]),
            float(recorded_baseline["candidate_A_per_um"]),
        ):
            failures.append(f"baseline branch mismatch at {local_key}")
            break

    paired = all(
        row["same_edge_paired_baseline"] == "True" for row in orders
    )
    if not paired:
        failures.append("replacement order lacks same-edge paired baseline")
    if manifest.get("fitted_parameter_count") != 0:
        failures.append("factorial introduced a fitted parameter")
    if not close(
        maximum_closure,
        float(manifest["maximum_shapley_closure_dex"]),
    ):
        failures.append("manifest closure maximum mismatch")
    if manifest.get("ranking_stable") is not False:
        failures.append("rank disagreement was not retained")
    if manifest.get("typed_outcome") != "interaction_dominant":
        failures.append("typed outcome differs from recomputed evidence")
    if not all(
        {int(row["node0"]), int(row["node1"])} == {1, 5}
        and math.isfinite(float(row["reference_abs_A_per_um"]))
        and math.isfinite(float(row["candidate_abs_A_per_um"]))
        for row in central
    ):
        failures.append("central 1-5 tail is not separately preserved")

    verification = {
        "status": "passed" if not failures else "failed",
        "failure_count": len(failures),
        "failures": failures,
        "active_edge_count": len(lattice),
        "maximum_recomputed_shapley_closure_dex": maximum_closure,
        "baseline_and_target_source_branches_close": not any(
            "branch mismatch" in failure for failure in failures
        ),
        "central_tail_separately_preserved": not any(
            "central 1-5 tail" in failure for failure in failures
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
