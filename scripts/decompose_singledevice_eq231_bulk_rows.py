#!/usr/bin/env python3
"""Decompose ordinary-Silicon SingleDevice Eq. 231 rows.

The report separates the Vela Formula-0 row into its linear P1 Laplacian,
fitted nonlinear correction, and reaction terms.  When Sentaurus NewtonPlot
perturbation exports are supplied, it also recovers a row-local unit scale
from the theta=0 off-diagonal Jacobian and compares the full Formula-0
Jacobian without fitting that scale to the nonlinear result.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path


def read_field(export_dir: Path, name: str, region: int) -> dict[int, float]:
    manifest = json.loads(
        (export_dir / "field_manifest.json").read_text(encoding="utf-8")
    )
    matches = [
        item for item in manifest["fields"]
        if item["name"] == name and int(item["region"]) == region
    ]
    if len(matches) != 1:
        raise ValueError(
            f"expected one {name!r} field for region {region} in "
            f"{export_dir}, got {len(matches)}"
        )
    with (export_dir / "fields" / matches[0]["csv_file"]).open(
        newline="", encoding="utf-8"
    ) as stream:
        return {
            int(row["node_id"]): float(row["component0"])
            for row in csv.DictReader(stream)
        }


def triangle_geometry(points: list[tuple[float, float]]) -> tuple[float, list[float], list[float]]:
    x = [point[0] for point in points]
    y = [point[1] for point in points]
    twice_area = abs(
        (x[1] - x[0]) * (y[2] - y[0])
        - (x[2] - x[0]) * (y[1] - y[0])
    )
    if twice_area == 0.0:
        raise ValueError("degenerate triangle")
    return (
        0.5 * twice_area,
        [y[1] - y[2], y[2] - y[0], y[0] - y[1]],
        [x[2] - x[1], x[0] - x[2], x[1] - x[0]],
    )


def fitted_jump(value: float) -> float:
    return math.expm1(value) if value < 0.0 else value + 0.5 * value * value


def fitted_jump_derivative(value: float) -> float:
    return math.exp(value) if value < 0.0 else 1.0 + value


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--cells", type=Path, required=True)
    parser.add_argument("--nodes", type=Path, required=True)
    parser.add_argument("--baseline-export", type=Path)
    parser.add_argument("--theta-zero-export", type=Path)
    parser.add_argument("--reference-lambda-csv", type=Path)
    parser.add_argument("--reaction-coefficient-scale", type=float, default=1.0)
    parser.add_argument("--full-perturb-root", type=Path)
    parser.add_argument("--theta-zero-perturb-root", type=Path)
    parser.add_argument("--row-nodes", type=int, nargs="+", required=True)
    parser.add_argument("--region", type=int, default=3)
    parser.add_argument("--thermal-voltage", type=float, default=0.025851999786)
    parser.add_argument("--delta", type=float, default=1.0e-5)
    parser.add_argument("--output-json", type=Path)
    parser.add_argument("--output-csv", type=Path)
    args = parser.parse_args()

    mesh = json.loads(args.mesh.read_text(encoding="utf-8"))
    coordinates = {
        int(node["id"]): (float(node["x"]), float(node["y"]))
        for node in mesh["nodes"]
    }
    triangles = {
        int(cell["id"]): [int(value) for value in cell.get("node_ids", cell.get("nodes"))]
        for cell in mesh["triangles"]
    }

    rows_by_cell: dict[int, list[dict[str, str]]] = defaultdict(list)
    with args.cells.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            rows_by_cell[int(row["cell_id"])].append(row)
    for rows in rows_by_cell.values():
        rows.sort(key=lambda row: int(row["local_node"]))

    node_diagnostics: dict[int, dict[str, str]] = {}
    with args.nodes.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            node_diagnostics[int(row["node_id"])] = row

    if (args.baseline_export is None) != (args.theta_zero_export is None):
        raise ValueError(
            "--baseline-export and --theta-zero-export must be supplied together"
        )
    baseline_rhs: dict[int, float] = {}
    theta_zero_rhs: dict[int, float] = {}
    baseline_state: dict[int, float] = {}
    if args.baseline_export is not None:
        baseline_rhs = read_field(
            args.baseline_export, "eQuantumPotentialRhs", args.region
        )
        theta_zero_rhs = read_field(
            args.theta_zero_export, "eQuantumPotentialRhs", args.region
        )
        baseline_state = read_field(
            args.baseline_export, "eQuantumPotential", args.region
        )
    if args.reference_lambda_csv is not None:
        with args.reference_lambda_csv.open(newline="", encoding="utf-8") as stream:
            restart_lambda = {
                int(row["node_id"]): float(row["electron_quantum_potential_V"])
                for row in csv.DictReader(stream)
            }
        if baseline_state:
            for node in args.row_nodes:
                if not math.isclose(
                    baseline_state[node], restart_lambda[node],
                    rel_tol=0.0, abs_tol=1.0e-12,
                ):
                    raise ValueError(
                        f"reference lambda mismatch at node {node}: "
                        f"{baseline_state[node]} versus {restart_lambda[node]}"
                    )
        baseline_state = restart_lambda
    if not baseline_state:
        raise ValueError(
            "one of --baseline-export or --reference-lambda-csv is required"
        )

    reports = []
    flat_cells = []
    for row_node in args.row_nodes:
        per_cell = []
        p1_jacobian: dict[int, float] = defaultdict(float)
        fitted_jacobian: dict[int, float] = defaultdict(float)
        fitted_stiffness = 0.0
        p1_stiffness = 0.0
        reaction = 0.0
        reaction_diagonal = 0.0

        for cell_id, records in rows_by_cell.items():
            row_records = [
                row for row in records
                if int(row["node_id"]) == row_node
                and row["is_active"] == "1"
                and row["is_dirichlet"] == "0"
            ]
            if not row_records:
                continue
            row_record = row_records[0]
            nodes = triangles[cell_id]
            local_row = int(row_record["local_node"])
            records_by_local = {int(row["local_node"]): row for row in records}
            points = [coordinates[node] for node in nodes]
            area, b, c = triangle_geometry(points)
            four_area_squared = 4.0 * area * area
            w = [float(records_by_local[local]["w"]) for local in range(3)]

            local_p1 = 0.0
            local_fitted = 0.0
            edge_rows = []
            for local, node in enumerate(nodes):
                stiffness = area * (
                    b[local_row] * b[local]
                    + c[local_row] * c[local]
                ) / four_area_squared
                local_p1 += -stiffness * w[local]
                p1_jacobian[node] += stiffness / args.thermal_voltage
                if local == local_row:
                    continue
                half_jump = 0.5 * (w[local] - w[local_row])
                contribution = -2.0 * stiffness * fitted_jump(half_jump)
                derivative = (
                    stiffness * fitted_jump_derivative(half_jump)
                    / args.thermal_voltage
                )
                local_fitted += contribution
                fitted_jacobian[node] += derivative
                fitted_jacobian[row_node] -= derivative
                edge_rows.append({
                    "column_node": node,
                    "stiffness": stiffness,
                    "half_jump": half_jump,
                    "p1_contribution": -stiffness * (w[local] - w[local_row]),
                    "formula0_contribution": contribution,
                    "formula0_minus_p1": (
                        contribution + stiffness * (w[local] - w[local_row])
                    ),
                })

            recorded_fitted = float(row_record["stiffness"])
            local_reaction = float(row_record["reaction"])
            local_lambda = float(row_record["lambda_V"])
            local_reaction_diagonal = (
                local_reaction / local_lambda
                if abs(local_lambda) > 1.0e-30 else float("nan")
            )
            fitted_stiffness += recorded_fitted
            p1_stiffness += local_p1
            reaction += local_reaction
            if math.isfinite(local_reaction_diagonal):
                reaction_diagonal += local_reaction_diagonal
            per_cell.append({
                "cell_id": cell_id,
                "nodes": nodes,
                "area_m2": float(row_record["area_m2"]),
                "lambda_V": local_lambda,
                "w": w,
                "p1_stiffness": local_p1,
                "formula0_stiffness": recorded_fitted,
                "formula0_recomputed": local_fitted,
                "formula0_minus_p1": recorded_fitted - local_p1,
                "reaction": local_reaction,
                "reaction_diagonal_per_V": local_reaction_diagonal,
                "total": recorded_fitted + local_reaction,
                "edges": edge_rows,
            })

        fitted_jacobian[row_node] += reaction_diagonal
        p1_jacobian[row_node] += reaction_diagonal
        diagnostic = node_diagnostics[row_node]
        closure_lambda = (
            -fitted_stiffness / reaction_diagonal
            if abs(reaction_diagonal) > 1.0e-30 else float("nan")
        )
        reference_lambda = baseline_state[row_node]
        corrected_reaction_diagonal = (
            reaction_diagonal * args.reaction_coefficient_scale
        )
        corrected_reaction = corrected_reaction_diagonal * reference_lambda
        report = {
            "node_id": row_node,
            "coordinate_internal": [
                float(diagnostic["x_internal"]),
                float(diagnostic["y_internal"]),
            ],
            "potential_like_V": float(diagnostic["potential_like_V"]),
            "output_lambda_V": float(diagnostic["output_lambda_V"]),
            "vela": {
                "p1_stiffness": p1_stiffness,
                "formula0_stiffness": fitted_stiffness,
                "formula0_nonlinear_correction": fitted_stiffness - p1_stiffness,
                "reaction": reaction,
                "reaction_diagonal_per_V": reaction_diagonal,
                "p1_total": p1_stiffness + reaction,
                "formula0_total": fitted_stiffness + reaction,
                "closure_lambda_V": closure_lambda,
                "closure_minus_output_lambda_V": (
                    closure_lambda - float(diagnostic["output_lambda_V"])
                ),
                "reference_lambda_V": reference_lambda,
                "reference_minus_output_lambda_V": (
                    reference_lambda - float(diagnostic["output_lambda_V"])
                ),
                "reaction_coefficient_scale": args.reaction_coefficient_scale,
                "corrected_reaction_diagonal_per_V": corrected_reaction_diagonal,
                "corrected_reference_reaction": corrected_reaction,
                "corrected_reference_formula0_total": (
                    fitted_stiffness + corrected_reaction
                ),
            },
            "cells": per_cell,
        }
        if baseline_rhs:
            report["sentaurus_rhs"] = {
                "formula0": baseline_rhs[row_node],
                "theta_zero": theta_zero_rhs[row_node],
            }

        neighbors = sorted(set(p1_jacobian) | set(fitted_jacobian))
        if (
            baseline_rhs
            and args.full_perturb_root
            and args.theta_zero_perturb_root
        ):
            full_derivatives = {}
            theta_derivatives = {}
            observed_deltas = {}
            for node in neighbors:
                full_dir = args.full_perturb_root / f"full_{node}_export"
                theta_dir = args.theta_zero_perturb_root / f"theta0_{node}_export"
                if not full_dir.is_dir() or not theta_dir.is_dir():
                    continue
                perturbed_full_rhs = read_field(
                    full_dir, "eQuantumPotentialRhs", args.region
                )
                perturbed_theta_rhs = read_field(
                    theta_dir, "eQuantumPotentialRhs", args.region
                )
                perturbed_state = read_field(
                    full_dir, "eQuantumPotential", args.region
                )
                observed_delta = perturbed_state[node] - baseline_state[node]
                if abs(observed_delta) < 1.0e-15:
                    observed_delta = args.delta
                observed_deltas[node] = observed_delta
                full_derivatives[node] = (
                    perturbed_full_rhs[row_node] - baseline_rhs[row_node]
                ) / observed_delta
                theta_derivatives[node] = (
                    perturbed_theta_rhs[row_node] - theta_zero_rhs[row_node]
                ) / observed_delta

            scale_samples = [
                p1_jacobian[node] / theta_derivatives[node]
                for node in neighbors
                if node != row_node
                and node in theta_derivatives
                and abs(theta_derivatives[node]) > 1.0e-20
                and abs(p1_jacobian[node]) > 1.0e-20
            ]
            scale = statistics.median(scale_samples)
            jacobian_rows = []
            for node in neighbors:
                if node not in full_derivatives or node not in theta_derivatives:
                    continue
                sentaurus_full_scaled = full_derivatives[node] * scale
                sentaurus_theta_scaled = theta_derivatives[node] * scale
                jacobian_rows.append({
                    "column_node": node,
                    "observed_delta_V": observed_deltas[node],
                    "sentaurus_formula0_scaled": sentaurus_full_scaled,
                    "sentaurus_theta_zero_scaled": sentaurus_theta_scaled,
                    "vela_formula0": fitted_jacobian.get(node, 0.0),
                    "vela_theta_zero": p1_jacobian.get(node, 0.0),
                    "formula0_difference": (
                        fitted_jacobian.get(node, 0.0) - sentaurus_full_scaled
                    ),
                    "theta_zero_difference": (
                        p1_jacobian.get(node, 0.0) - sentaurus_theta_scaled
                    ),
                })
            sentaurus_theta_residual_scaled = theta_zero_rhs[row_node] * scale
            sentaurus_theta_diagonal = theta_derivatives[row_node] * scale
            sentaurus_theta_offdiagonal_sum = sum(
                theta_derivatives[node] * scale
                for node in neighbors
                if node != row_node and node in theta_derivatives
            )
            sentaurus_reaction_diagonal = (
                sentaurus_theta_diagonal + sentaurus_theta_offdiagonal_sum
            )
            report["sentaurus_scaled"] = {
                "theta_zero_to_vela_scale": scale,
                "scale_samples": scale_samples,
                "formula0_residual": baseline_rhs[row_node] * scale,
                "theta_zero_residual": sentaurus_theta_residual_scaled,
                "inferred_nonlinear_correction": (
                    baseline_rhs[row_node] - theta_zero_rhs[row_node]
                ) * scale,
                "vela_p1_minus_sentaurus_theta_zero": (
                    p1_stiffness + reaction - sentaurus_theta_residual_scaled
                ),
                "vela_nonlinear_minus_sentaurus": (
                    fitted_stiffness - p1_stiffness
                    - (baseline_rhs[row_node] - theta_zero_rhs[row_node]) * scale
                ),
                "reaction_diagonal_per_V": sentaurus_reaction_diagonal,
                "reaction_diagonal_ratio_to_vela": (
                    sentaurus_reaction_diagonal / reaction_diagonal
                ),
                "reference_lambda_reaction": (
                    sentaurus_reaction_diagonal * reference_lambda
                ),
                "reference_lambda_formula0_total": (
                    fitted_stiffness
                    + sentaurus_reaction_diagonal * reference_lambda
                ),
                "jacobian": jacobian_rows,
            }

        reports.append(report)
        for cell in per_cell:
            flat_cells.append({
                "node_id": row_node,
                "cell_id": cell["cell_id"],
                "p1_stiffness": cell["p1_stiffness"],
                "formula0_stiffness": cell["formula0_stiffness"],
                "formula0_minus_p1": cell["formula0_minus_p1"],
                "reaction": cell["reaction"],
                "total": cell["total"],
            })

    output = {
        "schema": "vela.singledevice.eq231.bulk_row_decomposition.v1",
        "region": args.region,
        "thermal_voltage_V": args.thermal_voltage,
        "delta_V": args.delta,
        "rows": reports,
    }
    payload = json.dumps(output, indent=2) + "\n"
    if args.output_json:
        args.output_json.parent.mkdir(parents=True, exist_ok=True)
        args.output_json.write_text(payload, encoding="utf-8")
    else:
        print(payload, end="")
    if args.output_csv:
        args.output_csv.parent.mkdir(parents=True, exist_ok=True)
        with args.output_csv.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(flat_cells[0]))
            writer.writeheader()
            writer.writerows(flat_cells)


if __name__ == "__main__":
    main()
