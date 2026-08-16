#!/usr/bin/env python3
"""Fit the PolySilicon/SiO2 affine Eq. 231 trace with a central Jacobian oracle."""

from __future__ import annotations

import argparse
import csv
import json
import statistics
from collections import defaultdict
from pathlib import Path

import numpy as np

from analyze_singledevice_eq231_sentaurus_jacobian import (
    read_scalar,
    read_scalar_occurrences,
)


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def central_derivative(
    root: Path, row: int, region: int, column: int
) -> float:
    positive = root / f"import_n{column}_p"
    negative = root / f"import_n{column}_m"
    rhs_positive = read_scalar(
        positive, "eQuantumPotentialRhs", region
    )[row]
    rhs_negative = read_scalar(
        negative, "eQuantumPotentialRhs", region
    )[row]
    state_positive = read_scalar_occurrences(
        positive, "eQuantumPotential"
    )[column]
    state_negative = read_scalar_occurrences(
        negative, "eQuantumPotential"
    )[column]
    deltas = [
        plus - minus
        for plus, minus in zip(state_positive, state_negative)
        if plus != minus
    ]
    if not deltas:
        raise ValueError(f"no central state delta for node {column}")
    return (rhs_positive - rhs_negative) / max(deltas, key=abs)


def state_equations(
    label: str,
    cells_path: Path,
    nodes_path: Path,
    current_poly_weight: float,
    current_oxide_weight: float,
    current_poly_offset: float,
    current_oxide_offset: float,
) -> list[dict[str, object]]:
    nodes = {int(row["node_id"]): row for row in read_rows(nodes_path)}
    reaction = defaultdict(lambda: defaultdict(float))
    materials = defaultdict(set)
    for row in read_rows(cells_path):
        node = int(row["node_id"])
        material = row["material"]
        if material not in {"PolySilicon", "SiO2"}:
            continue
        reaction[node][material] += float(row["reaction"])
        materials[node].add(material)

    result = []
    for node, material_set in materials.items():
        record = nodes[node]
        if (
            material_set != {"PolySilicon", "SiO2"}
            or record["is_active"] != "1"
            or record["is_dirichlet"] != "0"
        ):
            continue
        output_lambda = float(record["output_lambda_V"])
        poly_denominator = (
            current_poly_weight * output_lambda + current_poly_offset
        )
        oxide_denominator = (
            current_oxide_weight * output_lambda + current_oxide_offset
        )
        poly_coefficient = reaction[node]["PolySilicon"] / poly_denominator
        oxide_coefficient = reaction[node]["SiO2"] / oxide_denominator
        result.append({
            "state": label,
            "node": node,
            "lambda_V": output_lambda,
            "features": np.asarray([
                poly_coefficient * output_lambda,
                oxide_coefficient * output_lambda,
                poly_coefficient,
                oxide_coefficient,
            ]),
            "target": -float(record["stiffness"]),
            "poly_coefficient": poly_coefficient,
            "oxide_coefficient": oxide_coefficient,
        })
    return result


def constrained_least_squares(
    matrix: np.ndarray,
    target: np.ndarray,
    constraint: np.ndarray,
    value: float,
) -> np.ndarray:
    norm_squared = float(constraint @ constraint)
    if norm_squared == 0.0:
        raise ValueError("zero Jacobian constraint")
    particular = constraint * (value / norm_squared)
    _, _, vh = np.linalg.svd(constraint.reshape(1, -1), full_matrices=True)
    nullspace = vh[1:].T
    reduced, _, _, _ = np.linalg.lstsq(
        matrix @ nullspace,
        target - matrix @ particular,
        rcond=None,
    )
    return particular + nullspace @ reduced


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--state", action="append", nargs=3, required=True,
        metavar=("LABEL", "CELLS", "NODES"),
    )
    parser.add_argument("--fitted-edges", type=Path, required=True)
    parser.add_argument("--central-root", type=Path, required=True)
    parser.add_argument("--row-node", type=int, default=2075)
    parser.add_argument("--region", type=int, default=0)
    parser.add_argument(
        "--columns", type=int, nargs="+",
        default=[2075, 2072, 2073, 2074, 2120],
    )
    parser.add_argument("--current-poly-weight", type=float, required=True)
    parser.add_argument("--current-oxide-weight", type=float, required=True)
    parser.add_argument("--current-poly-offset", type=float, default=0.0)
    parser.add_argument("--current-oxide-offset", type=float, default=0.0)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    fitted_by_column: dict[int, float] = defaultdict(float)
    for edge in read_rows(args.fitted_edges):
        if int(edge["row_node"]) == args.row_node:
            fitted_by_column[int(edge["column_node"])] += float(
                edge["jacobian_contribution"]
            )

    central = {
        column: central_derivative(
            args.central_root, args.row_node, args.region, column
        )
        for column in args.columns
    }
    controls = [column for column in args.columns if column != args.row_node]
    scales = [
        fitted_by_column[column] / central[column]
        for column in controls
        if central[column] != 0.0 and fitted_by_column[column] != 0.0
    ]
    scale = statistics.median(scales)
    stiffness_diagonal = -sum(fitted_by_column.values())
    sentaurus_total_diagonal = central[args.row_node] * scale
    sentaurus_reaction_diagonal = (
        sentaurus_total_diagonal - stiffness_diagonal
    )

    equations = []
    for label, cells, nodes in args.state:
        equations.extend(state_equations(
            label,
            Path(cells),
            Path(nodes),
            args.current_poly_weight,
            args.current_oxide_weight,
            args.current_poly_offset,
            args.current_oxide_offset,
        ))
    matrix = np.vstack([row["features"] for row in equations])
    target = np.asarray([row["target"] for row in equations])

    oracle_rows = [
        row for row in equations
        if row["state"] == args.state[0][0]
        and row["node"] == args.row_node
    ]
    if len(oracle_rows) != 1:
        raise ValueError("central oracle row is absent or ambiguous")
    oracle = oracle_rows[0]
    constraint = np.asarray([
        oracle["poly_coefficient"],
        oracle["oxide_coefficient"],
        0.0,
        0.0,
    ])
    fitted = constrained_least_squares(
        matrix, target, constraint, sentaurus_reaction_diagonal
    )
    predicted = matrix @ fitted
    errors = predicted - target

    by_state = {}
    for label, _, _ in args.state:
        indices = [
            index for index, row in enumerate(equations)
            if row["state"] == label
        ]
        state_errors = errors[indices]
        worst = sorted(indices, key=lambda index: abs(errors[index]), reverse=True)[:10]
        by_state[label] = {
            "node_count": len(indices),
            "max_abs_residual": float(np.max(np.abs(state_errors))),
            "rms_residual": float(np.sqrt(np.mean(state_errors ** 2))),
            "worst": [
                {
                    "node": equations[index]["node"],
                    "residual": float(errors[index]),
                }
                for index in worst
            ],
        }

    output = {
        "schema": "vela.singledevice_eq231_poly_oxide_central_closure.v1",
        "central_jacobian": {
            "row_node": args.row_node,
            "region": args.region,
            "raw_derivatives": central,
            "vela_by_column": dict(fitted_by_column),
            "control_scales": scales,
            "median_newtonplot_to_vela_scale": scale,
            "vela_stiffness_diagonal": stiffness_diagonal,
            "sentaurus_total_diagonal": sentaurus_total_diagonal,
            "sentaurus_reaction_diagonal": sentaurus_reaction_diagonal,
            "fitted_reaction_diagonal": float(constraint @ fitted),
        },
        "recommended_config": {
            "sentaurus_interface_polysilicon_reaction_weight": float(fitted[0]),
            "sentaurus_interface_insulator_at_polysilicon_reaction_weight": float(fitted[1]),
            "sentaurus_interface_polysilicon_reaction_offset_V": float(fitted[2]),
            "sentaurus_interface_insulator_at_polysilicon_reaction_offset_V": float(fitted[3]),
        },
        "global_max_abs_residual": float(np.max(np.abs(errors))),
        "global_rms_residual": float(np.sqrt(np.mean(errors ** 2))),
        "states": by_state,
    }
    text = json.dumps(output, indent=2) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    print(text, end="")


if __name__ == "__main__":
    main()
