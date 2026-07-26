#!/usr/bin/env python3
"""Independently verify PN2D imported-state QFP-update evidence."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any

FILES = (
    "topology_gate.csv",
    "residual_decomposition.csv",
    "first_qfp_updates.csv",
    "jacobian_blocks.csv",
    "disabled_controls.csv",
)
EXPECTED_COUNTS = {
    "topology_gate.csv": 9,
    "residual_decomposition.csv": 468,
    "first_qfp_updates.csv": 600,
    "jacobian_blocks.csv": 90,
    "disabled_controls.csv": 468,
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="ascii") as stream:
        return list(csv.DictReader(stream))


def finite(row: dict[str, str], keys: tuple[str, ...]) -> bool:
    return all(math.isfinite(float(row[key])) for key in keys)


def verify(root_a: Path, root_b: Path) -> dict[str, Any]:
    manifests = [json.loads((root / "manifest.json").read_text(encoding="ascii")) for root in (root_a, root_b)]
    errors: list[str] = []
    for index, (root, manifest) in enumerate(zip((root_a, root_b), manifests), start=1):
        if manifest.get("schema") != "pn2d_imported_state_qfp_update_v1":
            errors.append(f"root {index}: schema mismatch")
        if manifest.get("typed_outcome") != "operator_improvement_without_qfp_causality":
            errors.append(f"root {index}: unexpected typed outcome")
        if manifest.get("gates", {}).get("task8_authorized") is not False:
            errors.append(f"root {index}: Task 8 must remain unauthorized")
        for name in FILES:
            actual = sha256(root / name)
            if actual != manifest.get("hashes", {}).get(name):
                errors.append(f"root {index}: manifest hash mismatch for {name}")
            count = len(rows(root / name))
            if count != EXPECTED_COUNTS[name]:
                errors.append(f"root {index}: {name} row count {count}")
    for name in FILES:
        if sha256(root_a / name) != sha256(root_b / name):
            errors.append(f"A/B deterministic hash mismatch for {name}")

    residuals = rows(root_a / "residual_decomposition.csv")
    updates = rows(root_a / "first_qfp_updates.csv")
    jacobians = rows(root_a / "jacobian_blocks.csv")
    controls = rows(root_a / "disabled_controls.csv")
    topology = rows(root_a / "topology_gate.csv")
    if {float(row["bias_V"]) for row in updates} != {-1.0, -10.0, -20.0}:
        errors.append("first-update bias lattice mismatch")
    if {row["mode"] for row in updates} != {"carrier_only", "coupled"}:
        errors.append("first-update mode lattice mismatch")
    if {row["variant"] for row in updates} != {"production_triangle", "element_edge_opt_in"}:
        errors.append("first-update variant lattice mismatch")
    if not all(finite(row, ("delta_qfp_V", "residual_before", "residual_after")) for row in updates):
        errors.append("non-finite first update")
    max_coordinate_error = max(float(row["max_coordinate_error_um"]) for row in topology)
    if max_coordinate_error > 1.0e-12:
        errors.append(f"topology coordinate error {max_coordinate_error}")
    max_closure = max(float(row["closure_relative"]) for row in residuals)
    if max_closure > 1.0e-12:
        errors.append(f"residual decomposition closure {max_closure}")

    boundary: dict[tuple[Any, ...], dict[str, float]] = {}
    for row in residuals:
        if int(row["is_boundary"]) != 1:
            continue
        key = (row["topology"], row["bias_V"], row["carrier"], row["node_id"])
        boundary.setdefault(key, {})[row["variant"]] = float(row["final_residual_normalized"])
    max_boundary_difference = max(abs(pair["production_triangle"] - pair["element_edge_opt_in"]) for pair in boundary.values())
    if max_boundary_difference != 0.0:
        errors.append(f"boundary difference {max_boundary_difference}")
    if any(float(row["avalanche"]) != 0.0 for row in controls if row["control"] == "avalanche_off"):
        errors.append("avalanche-off control contains avalanche source")
    if any(float(row["srh"]) != 0.0 for row in controls if row["control"] == "srh_off"):
        errors.append("SRH-off control contains SRH source")

    grouped: dict[tuple[str, float, str, str], dict[str, float]] = {}
    for row in updates:
        key = (row["topology"], float(row["bias_V"]), row["mode"], row["carrier"])
        grouped.setdefault(key, {}).setdefault(row["variant"], 0.0)
        grouped[key][row["variant"]] += float(row["delta_qfp_V"]) ** 2
    non_improved = [key for key, pair in grouped.items() if math.sqrt(pair["element_edge_opt_in"]) >= math.sqrt(pair["production_triangle"])]
    if not any(key[0].startswith("minimal6_") for key in non_improved):
        errors.append("missing Minimal6 non-improvement evidence")
    if not any(key[0] == "coarse7x3" for key in non_improved):
        errors.append("missing coarse7x3 non-improvement evidence")

    nonzero_jacobian = [row for row in jacobians if row["block"] == "sg_avalanche" and max(abs(float(row["analytic_norm"])), abs(float(row["fd_norm"]))) >= 1.0e-8]
    imported_jacobian_max_relative = max(float(row["rel_diff"]) for row in nonzero_jacobian)
    imported_jacobian_gate = imported_jacobian_max_relative <= 1.0e-8

    for root in (root_a, root_b):
        for path in root.glob("raw/*/m*V/*_terms.json"):
            cfg = json.loads(path.read_text(encoding="ascii"))
            impact = cfg.get("solver", {}).get("impact_ionization")
            if impact is None:
                continue
            if impact.get("model") != "van_overstraeten" or impact.get("driving_force") != "quasi_fermi_gradient":
                errors.append(f"default model/driver changed in {path}")
            for fitted in ("A_scale", "B_scale"):
                if fitted in impact:
                    errors.append(f"fitted {fitted} present in {path}")

    return {
        "schema": "pn2d_imported_state_qfp_update_independent_verification_v1",
        "pass": not errors,
        "errors": errors,
        "sealed_hashes_match": all(sha256(root_a / name) == sha256(root_b / name) for name in FILES),
        "residual_closure_max_relative": max_closure,
        "boundary_max_abs_difference": max_boundary_difference,
        "topology_max_coordinate_error_um": max_coordinate_error,
        "first_update_non_improved_group_count": len(non_improved),
        "imported_state_sg_avalanche_jacobian_max_relative_nonzero": imported_jacobian_max_relative,
        "imported_state_sg_avalanche_jacobian_gate_1e_8": imported_jacobian_gate,
        "task8_authorized": False,
        "typed_outcome": "operator_improvement_without_qfp_causality",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root-a", type=Path, required=True)
    parser.add_argument("--root-b", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = verify(args.root_a.resolve(), args.root_b.resolve())
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="ascii", newline="\n")
    print(json.dumps({"pass": result["pass"], "typed_outcome": result["typed_outcome"]}, sort_keys=True))
    return 0 if result["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
