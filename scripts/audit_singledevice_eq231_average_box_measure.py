#!/usr/bin/env python3
"""Replace Vela's local mixed measure by a Sentaurus Measure debug oracle.

This is intentionally a fixed-state audit.  It answers whether the apparent
material-interface "reaction multiplier" is actually the local
AverageBoxMethod element-vertex measure on a non-Delaunay mesh.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from collections import defaultdict
from pathlib import Path


MEASURE_LINE = re.compile(
    r"^\s*(?P<grd>\d+)\s+(?P<des>-?\d+)\s+(?P<type>\d+)\s+"
    r"(?P<values>[^#]+?)\s*$"
)


def read_measures(path: Path) -> dict[int, list[float]]:
    text = path.read_text(encoding="utf-8")
    match = re.search(r"\n\s*Measure\s*\{.*?\n(.*?)\n\s*\}", text, re.DOTALL)
    if match is None:
        raise ValueError(f"Measure block not found in {path}")
    result: dict[int, list[float]] = {}
    for line in match.group(1).splitlines():
        parsed = MEASURE_LINE.match(line)
        if parsed is None or int(parsed.group("type")) != 2:
            continue
        result[int(parsed.group("des"))] = [
            float(value) for value in parsed.group("values").split()
        ]
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--cells", type=Path, required=True)
    parser.add_argument("--nodes", type=Path, required=True)
    parser.add_argument("--measure-debug", type=Path, required=True)
    parser.add_argument("--current-interface-multiplier", type=float, default=1.0)
    parser.add_argument("--limit", type=int, default=40)
    parser.add_argument(
        "--inspect-node",
        type=int,
        help="Report this active free node even when it is single-material.",
    )
    args = parser.parse_args()

    mesh = json.loads(args.mesh.read_text(encoding="utf-8"))
    triangle_nodes = {
        int(cell["id"]): [int(node) for node in cell.get("node_ids", cell.get("nodes"))]
        for cell in mesh["triangles"]
    }
    sentaurus = read_measures(args.measure_debug)
    with args.cells.open(newline="", encoding="utf-8") as stream:
        cells = list(csv.DictReader(stream))
    with args.nodes.open(newline="", encoding="utf-8") as stream:
        nodes = {int(row["node_id"]): row for row in csv.DictReader(stream)}

    material_sets: dict[int, set[str]] = defaultdict(set)
    current_reaction: dict[int, float] = defaultdict(float)
    average_reaction: dict[int, float] = defaultdict(float)
    measure_ratios: dict[int, list[float]] = defaultdict(list)
    rows_by_cell: dict[int, list[dict[str, str]]] = defaultdict(list)
    for row in cells:
        rows_by_cell[int(row["cell_id"])].append(row)
    for rows in rows_by_cell.values():
        rows.sort(key=lambda item: int(item["local_node"]))
    for row in cells:
        cell = int(row["cell_id"])
        local = int(row["local_node"])
        node = int(row["node_id"])
        if triangle_nodes[cell][local] != node:
            raise ValueError(
                f"cell {cell} local {local}: mesh node {triangle_nodes[cell][local]} != CSV node {node}"
            )
        if cell not in sentaurus:
            raise ValueError(f"cell {cell} missing from Sentaurus Measure block")
        area_m2 = float(row["area_m2"])
        # The diagnostic is from the current sentaurus_box assembly and hence
        # the current local measure is reaction-proportional.  Recompute that
        # measure geometrically from reaction only through the known
        # mixed-measure sum: each cell's three shares sum to its area.
        # We recover each local mixed share below from lambda/coefficient is
        # cumbersome, so use the same acute/obtuse geometry routine directly.
        points = [
            (
                float(rows_by_cell[cell][k]["x_internal"]),
                float(rows_by_cell[cell][k]["y_internal"]),
            )
            for k in range(3)
        ]
        vectors = []
        non_acute = -1
        for k in range(3):
            j, m = (k + 1) % 3, (k + 2) % 3
            uj = (points[j][0] - points[k][0], points[j][1] - points[k][1])
            um = (points[m][0] - points[k][0], points[m][1] - points[k][1])
            vectors.append((uj, um))
            if uj[0] * um[0] + uj[1] * um[1] <= 0.0 and non_acute < 0:
                non_acute = k
        if non_acute >= 0:
            shares = [0.25 * area_m2] * 3
            shares[non_acute] = 0.5 * area_m2
        else:
            def dist2(a: int, b: int) -> float:
                return sum((points[b][d] - points[a][d]) ** 2 for d in range(2))

            def cot(vertex: int, a: int, b: int) -> float:
                ua = (points[a][0] - points[vertex][0], points[a][1] - points[vertex][1])
                ub = (points[b][0] - points[vertex][0], points[b][1] - points[vertex][1])
                return (ua[0] * ub[0] + ua[1] * ub[1]) / abs(ua[0] * ub[1] - ua[1] * ub[0])

            shares_internal = []
            for k in range(3):
                j, m = (k + 1) % 3, (k + 2) % 3
                shares_internal.append(0.125 * (dist2(k, m) * cot(j, k, m) + dist2(k, j) * cot(m, k, j)))
            coordinate_area = 0.5 * abs(
                (points[1][0] - points[0][0]) * (points[2][1] - points[0][1])
                - (points[1][1] - points[0][1]) * (points[2][0] - points[0][0])
            )
            scale = area_m2 / coordinate_area
            shares = [share * scale for share in shares_internal]

        sent_measure_m2 = sentaurus[cell][local] * 1.0e-12
        ratio = sent_measure_m2 / shares[local]
        base_reaction = float(row["reaction"]) / args.current_interface_multiplier
        current_reaction[node] += base_reaction
        average_reaction[node] += base_reaction * ratio
        measure_ratios[node].append(ratio)
        material_sets[node].add(row["material"])

    output = []
    for node, materials in material_sets.items():
        record = nodes[node]
        selected = args.inspect_node is not None and node == args.inspect_node
        if (
            (len(materials) < 2 and not selected)
            or record["is_active"] != "1"
            or record["is_dirichlet"] != "0"
        ):
            continue
        stiffness = float(record["stiffness"])
        output.append({
            "node": node,
            "materials": "/".join(sorted(materials)),
            "stiffness": stiffness,
            "current_base_reaction": current_reaction[node],
            "sentaurus_average_reaction": average_reaction[node],
            "sentaurus_average_residual": stiffness + average_reaction[node],
            "min_local_measure_ratio": min(measure_ratios[node]),
            "max_local_measure_ratio": max(measure_ratios[node]),
        })
    output.sort(key=lambda item: abs(item["sentaurus_average_residual"]), reverse=True)
    if args.inspect_node is not None:
        output = [item for item in output if item["node"] == args.inspect_node]
    if not output:
        raise ValueError(f"no active free diagnostics for node {args.inspect_node}")
    writer = csv.DictWriter(__import__("sys").stdout, fieldnames=list(output[0]))
    writer.writeheader()
    writer.writerows(output[: args.limit])


if __name__ == "__main__":
    main()
