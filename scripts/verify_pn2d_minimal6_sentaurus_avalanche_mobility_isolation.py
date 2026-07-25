#!/usr/bin/env python3
"""Independently verify the avalanche mobility-isolation comparison."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any


REFERENCE = "lowfield_mobility_avalanche_electric_field"
CANDIDATE = "lowfield_mobility_avalanche_grad_qf"


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
    if manifest["status"] != "valid_sentaurus_avalanche_mobility_isolation":
        raise ValueError("comparison manifest status is invalid")
    if (
        manifest["reference_variant"] != REFERENCE
        or manifest["candidate_variant"] != CANDIDATE
    ):
        raise ValueError("variant contract mismatch")
    for name, expected in manifest["input_sha256"].items():
        if sha256(raw_root / Path(name)) != expected:
            raise ValueError(f"input hash mismatch: {name}")
    for name, expected in manifest["output_sha256"].items():
        if sha256(comparison_root / name) != expected:
            raise ValueError(f"output hash mismatch: {name}")

    states = load_csv(comparison_root / "state_summary.csv")
    quantities = load_csv(comparison_root / "quantity_comparison.csv")
    if len(states) != int(manifest["state_count"]):
        raise ValueError("state count mismatch")
    if len(quantities) != int(manifest["quantity_comparison_count"]):
        raise ValueError("quantity count mismatch")
    biases = {float(value) for value in manifest["biases_V"]}
    observed = {
        (
            row["topology"],
            float(row["bias_V"]),
            row["candidate_variant"],
            row["reference_variant"],
        )
        for row in states
    }
    expected = {
        (topology, bias, CANDIDATE, REFERENCE)
        for topology in ("mirror", "sketch")
        for bias in biases
    }
    if observed != expected:
        raise ValueError("state lattice mismatch")
    distinct = any(row["exact"] == "0" for row in quantities)
    if distinct != bool(manifest["candidate_distinct"]) or not distinct:
        raise ValueError("candidate distinctness mismatch")

    metrics = [key for key in states[0] if key.startswith("max_")]
    maxima = {}
    for metric in metrics:
        actual = max(abs(float(row[metric])) for row in states)
        expected_value = float(manifest["maxima"][metric])
        if not math.isclose(
            actual,
            expected_value,
            rel_tol=1.0e-15,
            abs_tol=1.0e-300,
        ):
            raise ValueError(f"maximum mismatch: {metric}")
        maxima[metric] = actual
    result = {
        "schema_version": 1,
        "status": "independently_verified",
        "experiment": manifest["experiment"],
        "sentaurus_release": manifest["sentaurus_release"],
        "state_count": len(states),
        "quantity_comparison_count": len(quantities),
        "candidate_distinct": distinct,
        "recomputed_maxima": maxima,
        "input_sha256": {
            "manifest.json": sha256(manifest_path),
            "state_summary.csv": sha256(
                comparison_root / "state_summary.csv"
            ),
            "quantity_comparison.csv": sha256(
                comparison_root / "quantity_comparison.csv"
            ),
        },
    }
    (output / "verification.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n",
        encoding="ascii",
        newline="\n",
    )
    return result


def main() -> int:
    args = parse_args()
    result = verify(args.raw_root, args.comparison_root, args.output)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
