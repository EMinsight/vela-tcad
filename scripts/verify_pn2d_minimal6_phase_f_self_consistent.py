#!/usr/bin/env python3
"""Independently verify one or two Phase F evidence roots."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def row_count(path: Path) -> int:
    with path.open(newline="", encoding="ascii") as handle:
        return sum(1 for _ in csv.DictReader(handle))


def verify(root: Path) -> dict:
    root = root.resolve()
    manifest = json.loads((root / "manifest.json").read_text(encoding="ascii"))
    if manifest["experiment"] != "pn2d_minimal6_phase_f_self_consistent":
        raise ValueError("unexpected Phase F experiment identity")
    for name, expected in manifest["outputs"].items():
        if sha256(root / name) != expected:
            raise ValueError(f"hash mismatch: {name}")
    expected_counts = {
        "state_node_comparison.csv": 80,
        "mobility_element_comparison.csv": 320,
        "directed_edge_current_comparison.csv": 400,
        "terminal_source_comparison.csv": 40,
        "summary.csv": 11,
    }
    for name, expected in expected_counts.items():
        if row_count(root / name) != expected:
            raise ValueError(f"row-count mismatch: {name}")
    if manifest["contracts"]["sentaurus_current_support"] != "box_operator_reconstruction":
        raise ValueError("Sentaurus current support is mislabeled")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("--compare-root", type=Path)
    args = parser.parse_args()
    first = verify(args.root)
    result = {
        "status": "passed",
        "outcome": first["outcome"]["status"],
        "verified_output_count": len(first["outputs"]),
    }
    if args.compare_root is not None:
        second = verify(args.compare_root)
        if first != second:
            raise ValueError("Phase F manifests are not byte-equivalent in content")
        result["deterministic_pair"] = True
    (args.root / "independent_verification.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n",
        encoding="ascii",
        newline="\n",
    )
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
