#!/usr/bin/env python3
"""Fit a cross-state affine Si/oxide Eq. 231 reaction trace.

The fixed-state rows constrain reaction residuals while a direct Sentaurus
NewtonPlot probe constrains the reaction Jacobian.  Combining both avoids the
false solution obtained by fitting a single multiplicative weight.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

import numpy as np


def endpoint_rows(
    label: str,
    cells_path: Path,
    nodes_path: Path,
    silicon_weight: float,
    oxide_weight: float,
) -> tuple[list[list[float]], list[float], list[dict[str, float | int | str]]]:
    with nodes_path.open(newline="", encoding="utf-8") as stream:
        nodes = {int(row["node_id"]): row for row in csv.DictReader(stream)}
    reactions: dict[int, dict[str, float]] = defaultdict(lambda: defaultdict(float))
    materials: dict[int, set[str]] = defaultdict(set)
    with cells_path.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            node = int(row["node_id"])
            material = row["material"]
            reactions[node][material] += float(row["reaction"])
            materials[node].add(material)

    matrix: list[list[float]] = []
    target: list[float] = []
    metadata: list[dict[str, float | int | str]] = []
    for node, material_set in materials.items():
        if material_set != {"Si", "SiO2"}:
            continue
        record = nodes[node]
        if record["is_active"] != "1" or record["is_dirichlet"] != "0":
            continue
        output_lambda = float(record["output_lambda_V"])
        if abs(output_lambda) <= 1.0e-14:
            continue
        silicon_coefficient = (
            reactions[node]["Si"] / silicon_weight / output_lambda
        )
        oxide_coefficient = (
            reactions[node]["SiO2"] / oxide_weight / output_lambda
        )
        matrix.append([
            silicon_coefficient * output_lambda,
            oxide_coefficient * output_lambda,
            silicon_coefficient,
            oxide_coefficient,
        ])
        target.append(-float(record["stiffness"]))
        metadata.append({
            "kind": "residual",
            "endpoint": label,
            "node": node,
            "lambda_V": output_lambda,
        })
    return matrix, target, metadata


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lin-cells", type=Path, required=True)
    parser.add_argument("--lin-nodes", type=Path, required=True)
    parser.add_argument("--sat-cells", type=Path, required=True)
    parser.add_argument("--sat-nodes", type=Path, required=True)
    parser.add_argument("--current-silicon-weight", type=float, required=True)
    parser.add_argument("--current-oxide-weight", type=float, required=True)
    parser.add_argument("--jacobian-node", type=int, required=True)
    parser.add_argument("--jacobian-target", type=float, required=True)
    parser.add_argument("--jacobian-weight", type=float, default=100.0)
    args = parser.parse_args()

    matrix: list[list[float]] = []
    target: list[float] = []
    metadata: list[dict[str, float | int | str]] = []
    endpoint_data = {}
    for label, cells, nodes in (
        ("lin", args.lin_cells, args.lin_nodes),
        ("sat", args.sat_cells, args.sat_nodes),
    ):
        rows, values, row_metadata = endpoint_rows(
            label, cells, nodes,
            args.current_silicon_weight,
            args.current_oxide_weight,
        )
        endpoint_data[label] = (rows, values, row_metadata)
        matrix.extend(rows)
        target.extend(values)
        metadata.extend(row_metadata)

    residual_design = np.asarray(
        [[row[0], row[1]] for row in matrix], dtype=float
    )
    residual_rhs = np.asarray(target, dtype=float)
    residual_weights, _, _, _ = np.linalg.lstsq(
        residual_design, residual_rhs, rcond=None
    )
    residual_error = residual_design @ residual_weights - residual_rhs

    lin_rows, _, lin_metadata = endpoint_data["lin"]
    match = [
        row for row, meta in zip(lin_rows, lin_metadata)
        if meta["node"] == args.jacobian_node
    ]
    if len(match) != 1:
        raise ValueError("Jacobian node is not one free Si/SiO2 interface row")
    jacobian_row = [match[0][2], match[0][3], 0.0, 0.0]
    matrix.append([args.jacobian_weight * value for value in jacobian_row])
    target.append(args.jacobian_weight * args.jacobian_target)
    metadata.append({"kind": "jacobian", "endpoint": "lin", "node": args.jacobian_node})

    design = np.asarray(matrix, dtype=float)
    rhs = np.asarray(target, dtype=float)
    coefficients, _, _, _ = np.linalg.lstsq(design, rhs, rcond=None)
    prediction = design @ coefficients
    error = prediction - rhs

    endpoint_summary = {}
    offset = 0
    for label in ("lin", "sat"):
        rows, values, row_metadata = endpoint_data[label]
        count = len(rows)
        local_error = error[offset: offset + count]
        order = np.argsort(np.abs(local_error))[::-1][:10]
        endpoint_summary[label] = {
            "node_count": count,
            "max_abs_residual": float(np.max(np.abs(local_error))),
            "rms_residual": float(np.sqrt(np.mean(local_error * local_error))),
            "worst": [
                {
                    "node": int(row_metadata[index]["node"]),
                    "residual": float(local_error[index]),
                }
                for index in order
            ],
        }
        offset += count

    print(json.dumps({
        "schema": "vela.singledevice.eq231.si_oxide_affine_closure.v1",
        "cross_state_multiplicative_fit": {
            "silicon_reaction_weight": float(residual_weights[0]),
            "oxide_reaction_weight": float(residual_weights[1]),
            "max_abs_residual": float(np.max(np.abs(residual_error))),
            "rms_residual": float(
                np.sqrt(np.mean(residual_error * residual_error))
            ),
        },
        "silicon_reaction_weight": float(coefficients[0]),
        "oxide_reaction_weight": float(coefficients[1]),
        "silicon_reaction_offset_V": float(coefficients[2]),
        "oxide_reaction_offset_V": float(coefficients[3]),
        "jacobian": {
            "node": args.jacobian_node,
            "target": args.jacobian_target,
            "predicted": float(np.dot(jacobian_row, coefficients)),
            "error": float(np.dot(jacobian_row, coefficients) - args.jacobian_target),
        },
        "endpoints": endpoint_summary,
    }, indent=2))


if __name__ == "__main__":
    main()
