#!/usr/bin/env python3
"""Fit a material-pair half-jump trace from adjacent oxide rows."""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


def branch(h: float) -> float:
    return __import__("math").expm1(h) if h < 0.0 else h + 0.5 * h * h


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--edges", type=Path, required=True)
    parser.add_argument("--cells", type=Path, required=True)
    parser.add_argument("--nodes", type=Path, required=True)
    parser.add_argument("--transport-material", required=True)
    parser.add_argument("--nontransport-material", required=True)
    parser.add_argument(
        "--side-material",
        help="material whose pure adjacent rows are fitted; defaults to nontransport",
    )
    parser.add_argument("--current-offset", type=float, default=0.0)
    args = parser.parse_args()

    materials: dict[int, set[str]] = defaultdict(set)
    with args.cells.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            materials[int(row["node_id"])].add(row["material"])
    pair = {args.transport_material, args.nontransport_material}
    side_material = args.side_material or args.nontransport_material
    pair_nodes = {node for node, values in materials.items() if values == pair}
    pure_side = {
        node for node, values in materials.items()
        if values == {side_material}
    }
    with args.nodes.open(newline="", encoding="utf-8") as stream:
        nodes = {int(row["node_id"]): row for row in csv.DictReader(stream)}
    affected: dict[int, list[dict[str, str]]] = defaultdict(list)
    with args.edges.open(newline="", encoding="utf-8") as stream:
        for edge in csv.DictReader(stream):
            row = int(edge["row_node"])
            column = int(edge["column_node"])
            if (
                edge["material"] == side_material
                and row in pure_side
                and column in pair_nodes
                and nodes[row]["is_active"] == "1"
                and nodes[row]["is_dirichlet"] == "0"
            ):
                affected[row].append(edge)

    def evaluate(offset: float) -> tuple[float, list[tuple[int, float]]]:
        delta = offset - args.current_offset
        residuals = []
        for row, edges in affected.items():
            correction = 0.0
            for edge in edges:
                h = float(edge["half_jump"])
                stiffness = float(edge["stiffness"])
                correction += -2.0 * stiffness * (branch(h + delta) - branch(h))
            residual = float(nodes[row]["raw_total"]) + correction
            residuals.append((row, residual))
        rms = (
            sum(residual * residual for _, residual in residuals) / len(residuals)
        ) ** 0.5
        return rms, residuals

    left, right = -0.25, 0.25
    phi = (5.0**0.5 - 1.0) / 2.0
    c = right - phi * (right - left)
    d = left + phi * (right - left)
    for _ in range(100):
        if evaluate(c)[0] < evaluate(d)[0]:
            right, d = d, c
            c = right - phi * (right - left)
        else:
            left, c = c, d
            d = left + phi * (right - left)
    offset = 0.5 * (left + right)
    rms, residuals = evaluate(offset)
    residuals.sort(key=lambda item: abs(item[1]), reverse=True)
    print(json.dumps({
        "material_pair": sorted(pair),
        "side_material": side_material,
        "adjacent_row_count": len(affected),
        "current_offset": args.current_offset,
        "fitted_offset": offset,
        "rms_residual": rms,
        "max_abs_residual": max(abs(value) for _, value in residuals),
        "worst": [
            {"node": node, "residual": residual}
            for node, residual in residuals[:10]
        ],
    }, indent=2))


if __name__ == "__main__":
    main()
