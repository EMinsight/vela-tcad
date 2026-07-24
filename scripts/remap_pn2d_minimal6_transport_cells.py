#!/usr/bin/env python3
"""Map Sentaurus region-cell order to the Minimal6 triangle ids."""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
from pathlib import Path


def load_psi(path: Path, topology: str) -> dict[int, float]:
    output: dict[int, float] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if (
                row["solver"] == "sentaurus"
                and row["topology"] == topology
                and float(row["bias_V"]) == -1.0
                and row["support_kind"] == "node"
                and row["quantity"] == "ElectrostaticPotential"
                and row["component"] == "component0"
                and row["status"] == "valid"
            ):
                output[int(row["support_id"])] = float(row["value_si"])
    if set(output) != set(range(6)):
        raise ValueError(f"missing Sentaurus psi for {topology}")
    return output


def p1_field(
    nodes: list[int],
    coordinates: dict[int, tuple[float, float]],
    psi: dict[int, float],
) -> tuple[float, float]:
    node0, node1, node2 = nodes
    x0, y0 = coordinates[node0]
    x1, y1 = coordinates[node1]
    x2, y2 = coordinates[node2]
    a11, a12 = x1 - x0, y1 - y0
    a21, a22 = x2 - x0, y2 - y0
    b1 = psi[node1] - psi[node0]
    b2 = psi[node2] - psi[node0]
    determinant = a11 * a22 - a12 * a21
    grad_x = (b1 * a22 - a12 * b2) / determinant
    grad_y = (a11 * b2 - b1 * a21) / determinant
    return -grad_x, -grad_y


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--observations", type=Path, required=True)
    parser.add_argument("--transport", type=Path, required=True)
    parser.add_argument("--mesh-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--mapping-output", type=Path, required=True)
    args = parser.parse_args()

    with args.transport.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    header = list(rows[0])
    mappings: dict[str, dict[int, int]] = {}
    mapping_rows: list[dict[str, object]] = []
    for topology in ("mirror", "sketch"):
        mesh = json.loads(
            (args.mesh_root / topology / "mesh.json").read_text(encoding="utf-8")
        )
        scale = 1.0e-6 if mesh["coordinate_unit"] == "um" else 1.0
        coordinates = {
            int(node["id"]): (float(node["x"]) * scale, float(node["y"]) * scale)
            for node in mesh["nodes"]
        }
        psi = load_psi(args.observations, topology)
        predicted = {
            int(triangle["id"]): p1_field(
                [int(node) for node in triangle["node_ids"]], coordinates, psi
            )
            for triangle in mesh["triangles"]
        }
        observed = {
            int(row["cell_id"]): (
                float(row["electric_field_x_V_per_m"]),
                float(row["electric_field_y_V_per_m"]),
            )
            for row in rows
            if row["topology"] == topology and float(row["bias_V"]) == -1.0
        }
        best: tuple[float, tuple[int, ...]] | None = None
        for permutation in itertools.permutations(range(4)):
            score = 0.0
            for triangle, sent_cell in enumerate(permutation):
                px, py = predicted[triangle]
                ox, oy = observed[sent_cell]
                scale2 = max(px * px + py * py, ox * ox + oy * oy, 1.0)
                score += ((px - ox) ** 2 + (py - oy) ** 2) / scale2
            candidate = (score, permutation)
            if best is None or candidate < best:
                best = candidate
        assert best is not None
        score, permutation = best
        mappings[topology] = {
            sent_cell: triangle
            for triangle, sent_cell in enumerate(permutation)
        }
        for triangle, sent_cell in enumerate(permutation):
            px, py = predicted[triangle]
            ox, oy = observed[sent_cell]
            residual = math.hypot(px - ox, py - oy) / max(
                math.hypot(px, py), math.hypot(ox, oy), 1.0
            )
            mapping_rows.append(
                {
                    "topology": topology,
                    "vela_triangle_id": triangle,
                    "sentaurus_region_cell_index": sent_cell,
                    "electric_field_relative_residual": residual,
                    "permutation_rms_relative_residual": math.sqrt(score / 4.0),
                }
            )

    for row in rows:
        row["cell_id"] = str(mappings[row["topology"]][int(row["cell_id"])])
    rows.sort(key=lambda row: (row["topology"], float(row["bias_V"]), int(row["cell_id"])))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=header)
        writer.writeheader()
        writer.writerows(rows)
    with args.mapping_output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(mapping_rows[0]))
        writer.writeheader()
        writer.writerows(mapping_rows)
    print(json.dumps({"status": "valid", "mappings": mappings}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
