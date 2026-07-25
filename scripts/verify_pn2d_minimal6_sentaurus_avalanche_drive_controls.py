#!/usr/bin/env python3
"""Independently verify the Minimal6 avalanche drive-control comparison."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any


EXPECTED_VARIANTS = {
    "implicit_default",
    "explicit_electric_field",
    "grad_qf_aval_dens_grad_qf",
}
EXPECTED_TOPOLOGIES = {"mirror", "sketch"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-root", type=Path, required=True)
    parser.add_argument("--comparison-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="ascii") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"empty CSV: {path}")
    return rows


def numeric_equal(left: float, right: float) -> bool:
    return math.isclose(left, right, rel_tol=1.0e-15, abs_tol=1.0e-300)


def verify(
    raw_root: Path,
    comparison_root: Path,
    output: Path,
) -> dict[str, Any]:
    raw_root = raw_root.resolve()
    comparison_root = comparison_root.resolve()
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=True)

    manifest_path = comparison_root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="ascii"))
    if manifest["status"] != "valid_sentaurus_avalanche_drive_comparison":
        raise ValueError("comparison manifest status is not valid")
    if manifest["reference_variant"] != "explicit_grad_qf":
        raise ValueError("comparison reference is not explicit_grad_qf")
    if set(manifest["candidate_variants"]) != EXPECTED_VARIANTS:
        raise ValueError("candidate variant contract mismatch")
    if set(manifest["topologies"]) != EXPECTED_TOPOLOGIES:
        raise ValueError("topology contract mismatch")

    state_path = comparison_root / "state_summary.csv"
    quantity_path = comparison_root / "quantity_comparison.csv"
    for name, expected in manifest["output_sha256"].items():
        actual = sha256(comparison_root / name)
        if actual != expected:
            raise ValueError(f"output hash mismatch: {name}")
    for name, expected in manifest["input_sha256"].items():
        actual = sha256(raw_root / Path(name))
        if actual != expected:
            raise ValueError(f"input hash mismatch: {name}")

    states = load_csv(state_path)
    quantities = load_csv(quantity_path)
    if len(states) != int(manifest["state_count"]):
        raise ValueError("state count mismatch")
    if len(quantities) != int(manifest["quantity_comparison_count"]):
        raise ValueError("quantity comparison count mismatch")

    biases = {float(value) for value in manifest["biases_V"]}
    observed_states = {
        (
            row["topology"],
            float(row["bias_V"]),
            row["candidate_variant"],
        )
        for row in states
    }
    expected_states = {
        (topology, bias, variant)
        for topology in EXPECTED_TOPOLOGIES
        for bias in biases
        for variant in EXPECTED_VARIANTS
    }
    if observed_states != expected_states:
        raise ValueError("state matrix mismatch")

    implicit_quantities = [
        row
        for row in quantities
        if row["candidate_variant"] == "implicit_default"
    ]
    if not implicit_quantities:
        raise ValueError("implicit-default comparison is missing")
    if any(row["exact"] != "1" for row in implicit_quantities):
        raise ValueError("implicit default is not an exact parsed-value match")
    if any(
        float(row["reference_value"]) != float(row["candidate_value"])
        for row in implicit_quantities
    ):
        raise ValueError("implicit default reference/candidate values differ")

    metric_names = list(next(iter(manifest["maxima_by_variant"].values())))
    recomputed_maxima = {}
    for variant in EXPECTED_VARIANTS:
        selected = [
            row for row in states if row["candidate_variant"] == variant
        ]
        recomputed_maxima[variant] = {}
        for metric in metric_names:
            actual = max(abs(float(row[metric])) for row in selected)
            expected = float(manifest["maxima_by_variant"][variant][metric])
            if not numeric_equal(actual, expected):
                raise ValueError(
                    f"manifest maximum mismatch: {variant}/{metric}"
                )
            recomputed_maxima[variant][metric] = actual

    electric_distinct = any(
        row["exact"] == "0"
        for row in quantities
        if row["candidate_variant"] == "explicit_electric_field"
        and (
            row["field"] in {"alpha_n_cm_inv", "alpha_p_cm_inv"}
            or row["field"].startswith("generation_")
        )
    )
    aval_dens_distinct = any(
        row["exact"] == "0"
        for row in quantities
        if row["candidate_variant"] == "grad_qf_aval_dens_grad_qf"
        and (
            row["field"].startswith("generation_")
            or row["field"].startswith("qg_")
            or "AvalancheIntegral" in row["field"]
        )
    )
    if electric_distinct != bool(manifest["explicit_electric_field_distinct"]):
        raise ValueError("electric-field distinctness mismatch")
    if aval_dens_distinct != bool(manifest["aval_dens_grad_qf_distinct"]):
        raise ValueError("AvalDensGradQF distinctness mismatch")
    if not bool(manifest["implicit_default_exact_match"]):
        raise ValueError("manifest did not pass implicit-default exact gate")

    verification = {
        "schema_version": 1,
        "status": "independently_verified",
        "experiment": manifest["experiment"],
        "reference_variant": manifest["reference_variant"],
        "state_count": len(states),
        "quantity_comparison_count": len(quantities),
        "implicit_default_exact_match": True,
        "explicit_electric_field_distinct": electric_distinct,
        "aval_dens_grad_qf_distinct": aval_dens_distinct,
        "recomputed_maxima_by_variant": recomputed_maxima,
        "input_sha256": {
            manifest_path.relative_to(comparison_root).as_posix(): sha256(
                manifest_path
            ),
            state_path.relative_to(comparison_root).as_posix(): sha256(
                state_path
            ),
            quantity_path.relative_to(comparison_root).as_posix(): sha256(
                quantity_path
            ),
        },
    }
    verification_path = output / "verification.json"
    verification_path.write_text(
        json.dumps(verification, indent=2, sort_keys=True) + "\n",
        encoding="ascii",
        newline="\n",
    )
    return verification


def main() -> int:
    args = parse_args()
    verification = verify(
        args.raw_root,
        args.comparison_root,
        args.output,
    )
    print(json.dumps(verification, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
