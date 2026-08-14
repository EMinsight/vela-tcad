#!/usr/bin/env python3
"""Audit Sentaurus Eq. 231 potential-like traces and element jumps."""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


def merge_with_regions(export_dir: Path, field_name: str) -> dict[int, list[tuple[int, float]]]:
    manifest = json.loads((export_dir / "field_manifest.json").read_text())
    values: dict[int, list[tuple[int, float]]] = defaultdict(list)
    for entry in manifest.get("fields", []):
        if (entry.get("name") != field_name or
                entry.get("mapping_status") != "complete"):
            continue
        with (export_dir / "fields" / entry["csv_file"]).open(newline="") as stream:
            for row in csv.DictReader(stream):
                values[int(row["node_id"])].append(
                    (int(entry["region"]), float(row["component0"])))
    return values


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--export-dir", type=Path, required=True)
    parser.add_argument("--nodes", type=int, nargs="*", default=[])
    args = parser.parse_args()
    fields = {
        name: merge_with_regions(args.export_dir, name)
        for name in ("ConductionBandEnergy", "eQuantumPotential",
                     "ElectrostaticPotential", "eQuasiFermiPotential",
                     "ElectronAffinity", "BandgapNarrowing")
    }
    potential_like = {
        node: fields["ConductionBandEnergy"][node][0][1] +
        fields["eQuantumPotential"][node][0][1]
        for node in fields["ConductionBandEnergy"]
    }
    reference_offsets = {
        node: fields["ConductionBandEnergy"][node][0][1] +
        fields["ElectrostaticPotential"][node][0][1] +
        fields["ElectronAffinity"][node][0][1]
        for node in fields["ConductionBandEnergy"]
    }
    duplicate_audit = {}
    for name, values in fields.items():
        differences = [
            max(value for _, value in occurrences) -
            min(value for _, value in occurrences)
            for occurrences in values.values() if len(occurrences) > 1
        ]
        duplicate_audit[name] = {
            "shared_node_occurrences": len(differences),
            "max_cross_region_difference": max(differences, default=0.0),
        }

    jumps = []
    material_nodes: dict[str, set[int]] = defaultdict(set)
    with (args.export_dir / "elements.csv").open(newline="") as stream:
        for row in csv.DictReader(stream):
            nodes = [int(row[f"node{i}"]) for i in range(3)]
            material_nodes[row["material"]].update(nodes)
            for first, second in ((0, 1), (1, 2), (2, 0)):
                n0, n1 = nodes[first], nodes[second]
                transport = row["material"] in {"Silicon", "PolySilicon"}
                drive = fields[
                    "eQuasiFermiPotential" if transport else
                    "ElectrostaticPotential"]
                w0 = (-drive[n0][0][1] - potential_like[n0]) / 0.025851999786
                w1 = (-drive[n1][0][1] - potential_like[n1]) / 0.025851999786
                jumps.append({
                    "jump_V": abs(potential_like[n0] - potential_like[n1]),
                    "drive_jump_V": abs(
                        drive[n0][0][1] - drive[n1][0][1]),
                    "w_jump": abs(w0 - w1),
                    "node0": n0,
                    "node1": n1,
                    "region": row["region"],
                    "material": row["material"],
                })
    potential_jumps = sorted(
        jumps, key=lambda row: row["jump_V"], reverse=True)
    w_jumps = sorted(jumps, key=lambda row: row["w_jump"], reverse=True)
    material_parameter_audit = {}
    for material, nodes in material_nodes.items():
        base_affinities = [
            fields["ElectronAffinity"][node][0][1] -
            0.5 * fields["BandgapNarrowing"][node][0][1]
            for node in nodes
        ]
        material_parameter_audit[material] = {
            "node_count": len(nodes),
            "inferred_base_affinity_min_eV": min(base_affinities),
            "inferred_base_affinity_max_eV": max(base_affinities),
            "inferred_base_affinity_mean_eV": (
                sum(base_affinities) / len(base_affinities)),
        }
    selected = {}
    for node in args.nodes:
        selected[str(node)] = {
            name: occurrences.get(node, [])
            for name, occurrences in fields.items()
        }
        selected[str(node)]["potential_like_V"] = potential_like.get(node)
    print(json.dumps({
        "potential_like_min_V": min(potential_like.values()),
        "potential_like_max_V": max(potential_like.values()),
        "band_energy_reference_min_V": min(reference_offsets.values()),
        "band_energy_reference_max_V": max(reference_offsets.values()),
        "duplicate_audit": duplicate_audit,
        "material_parameter_audit": material_parameter_audit,
        "largest_potential_like_element_edge_jumps": potential_jumps[:20],
        "largest_w_element_edge_jumps": w_jumps[:20],
        "selected_nodes": selected,
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
