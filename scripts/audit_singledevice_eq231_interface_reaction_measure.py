#!/usr/bin/env python3
"""Compare mixed-Voronoi and barycentric Eq. 231 interface reactions."""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


def mixed_measures(points: list[tuple[float, float]], area: float) -> list[float]:
    result = [0.0, 0.0, 0.0]
    for i in range(3):
        j = (i + 1) % 3
        k = (i + 2) % 3
        uj = (points[j][0] - points[i][0], points[j][1] - points[i][1])
        uk = (points[k][0] - points[i][0], points[k][1] - points[i][1])
        if uj[0] * uk[0] + uj[1] * uk[1] <= 0.0:
            result[:] = [0.25 * area] * 3
            result[i] = 0.5 * area
            return result

    def squared_distance(a: int, b: int) -> float:
        return sum((points[b][d] - points[a][d]) ** 2 for d in range(2))

    def cotangent(vertex: int, a: int, b: int) -> float:
        ua = (points[a][0] - points[vertex][0], points[a][1] - points[vertex][1])
        ub = (points[b][0] - points[vertex][0], points[b][1] - points[vertex][1])
        cross = abs(ua[0] * ub[1] - ua[1] * ub[0])
        return (ua[0] * ub[0] + ua[1] * ub[1]) / cross

    for i in range(3):
        j = (i + 1) % 3
        k = (i + 2) % 3
        result[i] = 0.125 * (
            squared_distance(i, k) * cotangent(j, i, k)
            + squared_distance(i, j) * cotangent(k, i, j)
        )
    scale = area / sum(result)
    return [value * scale for value in result]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cells", type=Path, required=True)
    parser.add_argument("--nodes", type=Path, required=True)
    parser.add_argument("--limit", type=int, default=40)
    args = parser.parse_args()

    with args.cells.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    by_cell: dict[int, list[dict[str, str]]] = defaultdict(list)
    materials: dict[int, set[str]] = defaultdict(set)
    for row in rows:
        by_cell[int(row["cell_id"])].append(row)
        materials[int(row["node_id"])].add(row["material"])

    mixed_reaction: dict[int, float] = defaultdict(float)
    barycentric_reaction: dict[int, float] = defaultdict(float)
    for cell_rows in by_cell.values():
        cell_rows.sort(key=lambda row: int(row["local_node"]))
        points = [
            (float(row["x_internal"]), float(row["y_internal"]))
            for row in cell_rows
        ]
        area = float(cell_rows[0]["area_m2"])
        # Coordinate units cancel in the barycentric/mixed ratio.
        coordinate_area = 0.5 * abs(
            (points[1][0] - points[0][0]) * (points[2][1] - points[0][1])
            - (points[1][1] - points[0][1]) * (points[2][0] - points[0][0])
        )
        measures = mixed_measures(points, coordinate_area)
        for local, row in enumerate(cell_rows):
            node = int(row["node_id"])
            reaction = float(row["reaction"])
            mixed_reaction[node] += reaction
            barycentric_reaction[node] += reaction * (
                coordinate_area / 3.0 / measures[local]
            )

    with args.nodes.open(newline="", encoding="utf-8") as stream:
        nodes = {int(row["node_id"]): row for row in csv.DictReader(stream)}
    output = []
    for node, material_set in materials.items():
        row = nodes[node]
        if len(material_set) < 2 or row["is_active"] != "1" or row["is_dirichlet"] != "0":
            continue
        stiffness = float(row["stiffness"])
        mixed = mixed_reaction[node]
        barycentric = barycentric_reaction[node]
        residual = float(row["raw_total"])
        output.append({
            "node": node,
            "materials": "/".join(sorted(material_set)),
            "abs_residual": abs(residual),
            "required_multiplier": -stiffness / mixed if mixed else math.nan,
            "barycentric_multiplier": barycentric / mixed if mixed else math.nan,
            "barycentric_residual": stiffness + barycentric,
        })
    output.sort(key=lambda row: row["abs_residual"], reverse=True)
    writer = csv.DictWriter(
        __import__("sys").stdout,
        fieldnames=list(output[0]) if output else ["node"],
    )
    writer.writeheader()
    writer.writerows(output[: args.limit])


if __name__ == "__main__":
    main()
