#!/usr/bin/env python3
"""Recover the fitted-flux and local interface terms in one Eq. 231 row."""

from __future__ import annotations

import argparse
import csv
import json
import statistics
from collections import defaultdict
from pathlib import Path

from analyze_singledevice_eq231_sentaurus_jacobian import (
    read_scalar,
    read_scalar_occurrences,
)


def derivative(
    baseline_export: Path,
    perturbation_root: Path,
    row: int,
    region: int,
    column: int,
) -> float:
    baseline_rhs = read_scalar(baseline_export, "eQuantumPotentialRhs", region)[row]
    baseline_state = read_scalar_occurrences(baseline_export, "eQuantumPotential")
    probe_dir = perturbation_root / f"perturbation_{column}_export"
    probe_rhs = read_scalar(probe_dir, "eQuantumPotentialRhs", region)[row]
    probe_state = read_scalar_occurrences(probe_dir, "eQuantumPotential")
    deltas = [
        value - reference
        for value in probe_state[column]
        for reference in baseline_state[column]
        if value != reference
    ]
    if not deltas:
        raise ValueError(f"no observed perturbation for node {column}")
    return (probe_rhs - baseline_rhs) / max(deltas, key=abs)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--edges", type=Path, required=True)
    parser.add_argument("--nodes", type=Path, required=True)
    parser.add_argument("--baseline-export", type=Path, required=True)
    parser.add_argument("--perturbation-root", type=Path, required=True)
    parser.add_argument("--row", type=int, required=True)
    parser.add_argument("--region", type=int, required=True)
    parser.add_argument("--controls", type=int, nargs="+", required=True)
    parser.add_argument("--current-interface-multiplier", type=float, default=1.0)
    args = parser.parse_args()

    with args.edges.open(newline="", encoding="utf-8") as stream:
        edge_rows = [
            edge for edge in csv.DictReader(stream)
            if int(edge["row_node"]) == args.row
        ]
    vela_offdiag: dict[int, float] = defaultdict(float)
    for edge in edge_rows:
        vela_offdiag[int(edge["column_node"])] += float(
            edge["jacobian_contribution"]
        )
    raw_controls = {
        column: derivative(
            args.baseline_export, args.perturbation_root,
            args.row, args.region, column,
        )
        for column in args.controls
    }
    scales = [
        vela_offdiag[column] / raw
        for column, raw in raw_controls.items()
        if raw != 0.0 and vela_offdiag[column] != 0.0
    ]
    scale = statistics.median(scales)
    raw_diagonal = derivative(
        args.baseline_export, args.perturbation_root,
        args.row, args.region, args.row,
    )
    sentaurus_diagonal = scale * raw_diagonal
    fitted_flux_diagonal = -sum(vela_offdiag.values())
    sentaurus_local_diagonal = sentaurus_diagonal - fitted_flux_diagonal

    with args.nodes.open(newline="", encoding="utf-8") as stream:
        node = next(
            row for row in csv.DictReader(stream)
            if int(row["node_id"]) == args.row
        )
    output_lambda = float(node["output_lambda_V"])
    stiffness = float(node["stiffness"])
    current_reaction = (
        float(node["reaction"]) / args.current_interface_multiplier
    )
    current_reaction_diagonal = current_reaction / output_lambda
    required_local_residual = -stiffness
    sentaurus_linear_reaction_at_state = sentaurus_local_diagonal * output_lambda
    sentaurus_affine_source = (
        required_local_residual - sentaurus_linear_reaction_at_state
    )

    print(json.dumps({
        "row": args.row,
        "region": args.region,
        "controls": {
            str(column): {
                "raw_sentaurus_derivative": raw_controls[column],
                "vela_shifted_fitted_derivative": vela_offdiag[column],
                "scale": (
                    vela_offdiag[column] / raw_controls[column]
                    if raw_controls[column] != 0.0 else None
                ),
            }
            for column in args.controls
        },
        "newtonplot_to_vela_scale": scale,
        "sentaurus_diagonal": sentaurus_diagonal,
        "fitted_flux_diagonal": fitted_flux_diagonal,
        "sentaurus_local_diagonal": sentaurus_local_diagonal,
        "current_reaction_diagonal": current_reaction_diagonal,
        "output_lambda_V": output_lambda,
        "required_local_residual": required_local_residual,
        "sentaurus_linear_reaction_at_state": sentaurus_linear_reaction_at_state,
        "sentaurus_affine_source": sentaurus_affine_source,
    }, indent=2))


if __name__ == "__main__":
    main()
