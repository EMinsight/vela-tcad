#!/usr/bin/env python3
"""Merge a multi-region Sentaurus TDR CSV export into a Vela restart state."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


FIELDS = {
    "psi": "ElectrostaticPotential",
    "phin": "eQuasiFermiPotential",
    "phip": "hQuasiFermiPotential",
    "electrons_m3": "eDensity",
    "holes_m3": "hDensity",
    "electron_quantum_potential_V": "eQuantumPotential",
}


def read_field(path: Path) -> dict[int, float]:
    values: dict[int, float] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            node = int(row["node_id"])
            value = float(row["component0"])
            if node in values and not math.isclose(values[node], value, rel_tol=1e-8, abs_tol=1e-10):
                raise ValueError(f"inconsistent duplicate node {node} in {path}")
            values[node] = value
    return values


def merge_field(export_dir: Path, manifest: dict, name: str) -> tuple[dict[int, float], str]:
    merged: dict[int, float] = {}
    unit = ""
    for entry in manifest.get("fields", []):
        if entry.get("name") != name or entry.get("mapping_status") != "complete":
            continue
        unit = str(entry.get("unit", unit))
        values = read_field(export_dir / "fields" / str(entry["csv_file"]))
        for node, value in values.items():
            if node in merged and not math.isclose(merged[node], value, rel_tol=1e-8, abs_tol=1e-10):
                raise ValueError(f"field {name} differs across regions at shared node {node}")
            merged[node] = value
    return merged, unit


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--export-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    manifest = json.loads((args.export_dir / "field_manifest.json").read_text())
    node_count = sum(1 for _ in (args.export_dir / "nodes.csv").open()) - 1
    columns: dict[str, dict[int, float]] = {}
    units: dict[str, str] = {}
    for output_name, sentaurus_name in FIELDS.items():
        columns[output_name], units[output_name] = merge_field(
            args.export_dir, manifest, sentaurus_name)

    required_all_nodes = ("psi", "phin", "phip", "electron_quantum_potential_V")
    for name in required_all_nodes:
        missing = set(range(node_count)) - columns[name].keys()
        if missing:
            raise ValueError(f"{name} is missing {len(missing)} nodes")

    # Densities exist only on carrier-supporting regions. Insulator values are
    # unused by the DD assembler and are written as zero.
    density_scale = {
        name: 1.0e6 if units[name].replace(" ", "") in {"cm^-3", "1/cm^3"} else 1.0
        for name in ("electrons_m3", "holes_m3")
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["node_id", *FIELDS])
        for node in range(node_count):
            # Sentaurus may carry arbitrary quasi-Fermi placeholders through
            # insulators.  Vela represents those non-transport rows as pinned
            # algebraic unknowns, so normalize nodes that have no density
            # support instead of feeding the placeholders to Newton.
            has_carrier_support = (
                node in columns["electrons_m3"] or
                node in columns["holes_m3"]
            )
            writer.writerow([
                node,
                columns["psi"][node],
                columns["phin"][node] if has_carrier_support else 0.0,
                columns["phip"][node] if has_carrier_support else 0.0,
                columns["electrons_m3"].get(node, 0.0) * density_scale["electrons_m3"],
                columns["holes_m3"].get(node, 0.0) * density_scale["holes_m3"],
                columns["electron_quantum_potential_V"][node],
            ])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
