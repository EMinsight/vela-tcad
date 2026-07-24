#!/usr/bin/env python3
"""Independently verify a PN2D Minimal6 Phase D evidence root."""

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

    parameters = rows(root / "parameter_comparison.csv")
    native = rows(root / "native_element_decomposition.csv")
    edges = rows(root / "box_edge_mobility_decomposition.csv")
    adjacent = rows(root / "box_edge_adjacent_elements.csv")
    central = rows(root / "central_edge_1_5_decomposition.csv")
    controls = rows(root / "parameter_substitution_controls.csv")
    summary = rows(root / "summary.csv")
    expected_counts = {
        "parameters": (len(parameters), 25),
        "native": (len(native), 320),
        "edges": (len(edges), 720),
        "adjacent": (len(adjacent), 960),
        "central": (len(central), 160),
        "controls": (len(controls), 4),
        "summary": (len(summary), 10),
    }
    for label, (actual, expected) in expected_counts.items():
        if actual != expected:
            failures.append(f"{label} count {actual} != {expected}")
    if any(row["status"] == "mismatch" for row in parameters):
        failures.append("documented numeric parameter mismatch")
    active = [row for row in edges if row["status"] == "valid"]
    if len(active) != 400:
        failures.append(f"active edge count {len(active)} != 400")
    for row in active:
        model = float(row["model_signed_log10_ratio_dex"])
        support = float(row["support_signed_log10_ratio_dex"])
        total = float(row["total_signed_log10_ratio_dex"])
        closure = float(row["decomposition_closure_dex"])
        if abs(total - model - support) > 2.0e-14 or abs(closure) > 2.0e-14:
            failures.append("edge log decomposition closure failed")
            break
    replay = [
        float(row["inferred_replay_relative_error"])
        for row in native if row["inferred_replay_relative_error"]
    ]
    if len(replay) != 320 or max(replay) > 2.0e-14:
        failures.append("native inferred mobility replay failed")
    documented = [
        row for row in controls
        if row["control"] == "sentaurus_documented_parameter_substitution"
    ]
    if len(documented) != 2 or any(
        float(row["maximum_change_from_vela_current_parameters"]) != 0.0
        for row in documented
    ):
        failures.append("documented parameter substitution is not a no-op")
    if not all(
        math.isfinite(float(row["median_abs_dex"]))
        and math.isfinite(float(row["p95_abs_dex"]))
        for row in summary
    ):
        failures.append("non-finite summary metric")
    result = {
        "status": "passed" if not failures else "failed",
        "failures": failures,
        "verified_output_count": len(manifest["outputs"]),
        "native_count": len(native),
        "active_edge_count": len(active),
        "maximum_inferred_replay_relative_error": max(replay) if replay else None,
    }
    (root / "independent_verification.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8", newline="\n"
    )
    print(json.dumps(result, sort_keys=True))
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
