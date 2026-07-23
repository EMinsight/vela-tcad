#!/usr/bin/env python3
"""Attach production C++ baseline edge mobilities to the current-proxy table."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def triangle_mobility(path: Path) -> dict[tuple[int, int, str], float]:
    values: dict[tuple[int, int, str], list[float]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            for local in range(3):
                prefix = f"local_edge{local}"
                pair = tuple(
                    sorted(
                        (
                            int(row[f"{prefix}_node0"]),
                            int(row[f"{prefix}_node1"]),
                        )
                    )
                )
                for carrier in ("electron", "hole"):
                    key = (pair[0], pair[1], carrier)
                    values.setdefault(key, []).append(
                        float(row[f"{prefix}_{carrier}_mobility_m2_per_V_s"])
                    )
    result: dict[tuple[int, int, str], float] = {}
    for key, samples in values.items():
        scale = max(max(abs(value) for value in samples), 1.0e-300)
        if max(samples) - min(samples) > 2.0e-14 * scale:
            raise ValueError(f"cell-local mobility disagrees for edge {key}: {samples}")
        result[key] = samples[0]
    if len(result) != 18:
        raise ValueError(f"expected 18 carrier-edge mobilities, got {len(result)}")
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--current-edges", type=Path, required=True)
    parser.add_argument("--baseline-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    cache: dict[tuple[str, float], dict[tuple[int, int, str], float]] = {}
    rows: list[dict[str, str]] = []
    with args.current_edges.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError("current-edge input has no header")
        fieldnames = list(reader.fieldnames)
        for row in reader:
            topology, bias = row["topology"], float(row["bias_V"])
            state_key = (topology, bias)
            if state_key not in cache:
                triangle = (
                    args.baseline_root
                    / topology
                    / f"m{abs(int(bias))}V"
                    / "triangles.csv"
                )
                cache[state_key] = triangle_mobility(triangle)
            node0, node1 = sorted((int(row["node0"]), int(row["node1"])))
            mobility = cache[state_key][(node0, node1, row["carrier"])]
            row["vela_masetti_native_state_mobility_m2_per_Vs"] = format(
                mobility, ".17g"
            )
            rows.append(row)
    if len(cache) != 40 or len(rows) != 720:
        raise ValueError(f"expected 40 states/720 rows, got {len(cache)}/{len(rows)}")
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    manifest = {
        "schema_version": 1,
        "status": "valid",
        "semantics": (
            "vela_masetti_native_state_mobility_m2_per_Vs is replaced by the "
            "actual production C++ baseline edge mobility from triangle audit output"
        ),
        "state_count": len(cache),
        "row_count": len(rows),
        "input_current_edges": str(args.current_edges.resolve()),
        "input_current_edges_sha256": sha256(args.current_edges),
        "baseline_root": str(args.baseline_root.resolve()),
        "output": str(args.output.resolve()),
        "output_sha256": sha256(args.output),
    }
    manifest_path = args.output.with_suffix(".manifest.json")
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"state_count": len(cache), "row_count": len(rows)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
