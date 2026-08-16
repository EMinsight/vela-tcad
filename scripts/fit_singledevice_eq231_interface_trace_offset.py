#!/usr/bin/env python3
"""Fit the Sentaurus default material-interface half-jump trace offset.

The scale of NewtonPlot's eQuantumPotentialRhs is arbitrary.  It is recovered
from non-interface columns whose Vela fitted-edge Jacobian is already known to
match.  The remaining interface-column derivative then determines the single
dimensionless half-jump offset without fitting a complete device curve.
"""

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


def fitted_edges(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def observed_derivative(
    baseline_export: Path,
    perturbation_root: Path,
    row_node: int,
    region: int,
    column_node: int,
) -> float:
    baseline_rhs = read_scalar(
        baseline_export, "eQuantumPotentialRhs", region
    )[row_node]
    baseline_state = read_scalar_occurrences(
        baseline_export, "eQuantumPotential"
    )
    export_dir = perturbation_root / f"perturbation_{column_node}_export"
    perturbed_rhs = read_scalar(
        export_dir, "eQuantumPotentialRhs", region
    )[row_node]
    perturbed_state = read_scalar_occurrences(
        export_dir, "eQuantumPotential"
    )
    deltas = [
        probe - reference
        for probe in perturbed_state[column_node]
        for reference in baseline_state[column_node]
        if probe != reference
    ]
    if not deltas:
        raise ValueError(f"no observed state delta for node {column_node}")
    delta = max(deltas, key=abs)
    return (perturbed_rhs - baseline_rhs) / delta


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--fitted-edges", type=Path, required=True)
    parser.add_argument("--nodes", type=Path, required=True)
    parser.add_argument("--baseline-export", type=Path, required=True)
    parser.add_argument("--perturbation-root", type=Path, required=True)
    parser.add_argument("--thermal-voltage", type=float, default=0.025851999786)
    parser.add_argument(
        "--probe", action="append", nargs=4,
        metavar=("ROW", "REGION", "INTERFACE_COLUMN", "CONTROL_COLUMNS"),
        required=True,
        help="CONTROL_COLUMNS is a comma-separated list",
    )
    args = parser.parse_args()

    edges = fitted_edges(args.fitted_edges)
    node_rows = {}
    with args.nodes.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            node_rows[int(row["node_id"])] = row

    results = []
    for row_text, region_text, interface_text, controls_text in args.probe:
        row_node = int(row_text)
        region = int(region_text)
        interface_column = int(interface_text)
        controls = [int(value) for value in controls_text.split(",")]
        row_edges = [
            edge for edge in edges if int(edge["row_node"]) == row_node
        ]
        vela_by_column: dict[int, float] = defaultdict(float)
        for edge in row_edges:
            vela_by_column[int(edge["column_node"])] += float(
                edge["jacobian_contribution"]
            )

        scales = []
        for column in controls:
            sentaurus = observed_derivative(
                args.baseline_export, args.perturbation_root,
                row_node, region, column,
            )
            if sentaurus != 0.0 and vela_by_column[column] != 0.0:
                scales.append(vela_by_column[column] / sentaurus)
        if not scales:
            raise ValueError(f"no usable controls for row {row_node}")
        scale = statistics.median(scales)

        sentaurus_interface = scale * observed_derivative(
            args.baseline_export, args.perturbation_root,
            row_node, region, interface_column,
        )
        vela_interface = vela_by_column[interface_column]
        interface_edges = [
            edge for edge in row_edges
            if int(edge["column_node"]) == interface_column
        ]
        stiffness_sum = sum(float(edge["stiffness"]) for edge in interface_edges)
        if stiffness_sum == 0.0:
            raise ValueError(f"zero interface stiffness for row {row_node}")
        half_jump_offset = (
            (sentaurus_interface - vela_interface)
            * args.thermal_voltage / stiffness_sum
        )

        residual = float(node_rows[row_node]["raw_total"])
        correction = 0.0
        for edge in interface_edges:
            h = float(edge["half_jump"])
            stiffness = float(edge["stiffness"])
            old = h + 0.5 * h * h
            shifted = (
                h + half_jump_offset
                + 0.5 * (h + half_jump_offset) ** 2
            )
            correction += -2.0 * stiffness * (shifted - old)

        results.append({
            "row_node": row_node,
            "region": region,
            "interface_column": interface_column,
            "control_columns": controls,
            "newtonplot_to_vela_scale": scale,
            "vela_interface_derivative": vela_interface,
            "sentaurus_interface_derivative": sentaurus_interface,
            "half_jump_offset": half_jump_offset,
            "potential_like_trace_shift_V": (
                -2.0 * args.thermal_voltage * half_jump_offset
            ),
            "current_residual": residual,
            "predicted_flux_correction": correction,
            "predicted_corrected_residual": residual + correction,
        })

    print(json.dumps({
        "probes": results,
        "median_half_jump_offset": statistics.median(
            result["half_jump_offset"] for result in results
        ),
    }, indent=2))


if __name__ == "__main__":
    main()
