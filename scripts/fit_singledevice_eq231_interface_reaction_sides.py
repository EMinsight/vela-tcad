#!/usr/bin/env python3
"""Fit material-side Eq. 231 reaction weights at shared interface nodes."""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

import numpy as np


TRANSPORT = {"Si", "Silicon", "PolySilicon", "Germanium", "GaAs", "SiC"}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cells", type=Path, required=True)
    parser.add_argument("--nodes", type=Path, required=True)
    parser.add_argument("--current-interface-multiplier", type=float, default=1.0)
    parser.add_argument(
        "--current-weight", action="append", default=[],
        metavar="MATERIAL=WEIGHT",
    )
    parser.add_argument("--group-by-region", action="store_true")
    args = parser.parse_args()
    current_weights = {
        key: float(value)
        for item in args.current_weight
        for key, value in [item.split("=", 1)]
    }

    with args.nodes.open(newline="", encoding="utf-8") as stream:
        nodes = {int(row["node_id"]): row for row in csv.DictReader(stream)}
    material_reaction: dict[int, dict[str, float]] = defaultdict(lambda: defaultdict(float))
    materials: dict[int, set[str]] = defaultdict(set)
    region_material: dict[str, str] = {}
    side_lambda: dict[int, dict[str, float]] = defaultdict(dict)
    with args.cells.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            node = int(row["node_id"])
            material = row["material"]
            key = row["region_name"] if args.group_by_region else material
            material_reaction[node][key] += (
                float(row["reaction"]) /
                current_weights.get(material, args.current_interface_multiplier)
            )
            materials[node].add(key)
            region_material[key] = material
            side_lambda[node][key] = float(row["lambda_V"])

    grouped: dict[tuple[str, str], list[tuple[int, float, float, float, float, float, float]]] = defaultdict(list)
    for node, material_set in materials.items():
        record = nodes[node]
        if record["is_active"] != "1" or record["is_dirichlet"] != "0":
            continue
        transport = sorted(
            value for value in material_set
            if region_material[value] in TRANSPORT
        )
        nontransport = sorted(
            value for value in material_set
            if region_material[value] not in TRANSPORT
        )
        if len(transport) != 1 or len(nontransport) != 1:
            continue
        left, right = transport[0], nontransport[0]
        grouped[(left, right)].append((
            node,
            material_reaction[node][left],
            material_reaction[node][right],
            -float(record["stiffness"]),
            float(record["output_lambda_V"]),
            side_lambda[node][left],
            side_lambda[node][right],
        ))

    output = {}
    for pair, rows in sorted(grouped.items()):
        matrix = np.asarray([[row[1], row[2]] for row in rows], dtype=float)
        target = np.asarray([row[3] for row in rows], dtype=float)
        weights, _, _, _ = np.linalg.lstsq(matrix, target, rcond=None)
        predicted = matrix @ weights
        errors = predicted - target
        worst = np.argsort(np.abs(errors))[::-1][:10]
        nonzero_rows = [row for row in rows if abs(row[4]) > 1.0e-12]
        affine_matrix = np.asarray([
            [row[1], row[2], row[1] / row[4], row[2] / row[4]]
            for row in nonzero_rows
        ], dtype=float)
        affine_target = np.asarray([row[3] for row in nonzero_rows], dtype=float)
        affine_weights, _, _, _ = np.linalg.lstsq(
            affine_matrix, affine_target, rcond=None
        )
        affine_predicted = affine_matrix @ affine_weights
        affine_errors = affine_predicted - affine_target
        affine_worst = np.argsort(np.abs(affine_errors))[::-1][:10]
        quadratic_matrix = np.asarray([
            [row[1], row[2], row[1] * row[4], row[2] * row[4]]
            for row in rows
        ], dtype=float)
        quadratic_weights, _, _, _ = np.linalg.lstsq(
            quadratic_matrix, target, rcond=None
        )
        quadratic_predicted = quadratic_matrix @ quadratic_weights
        quadratic_errors = quadratic_predicted - target
        quadratic_worst = np.argsort(np.abs(quadratic_errors))[::-1][:10]
        native_matrix = np.asarray([
            [row[1] / row[4] * row[5], row[2] / row[4] * row[6]]
            for row in rows
        ], dtype=float)
        native_weights, _, _, _ = np.linalg.lstsq(native_matrix, target, rcond=None)
        native_predicted = native_matrix @ native_weights
        native_errors = native_predicted - target
        trace_lambda = np.asarray([row[4] for row in rows], dtype=float)
        trace_coefficient = np.asarray([
            (row[1] + row[2]) / row[4] for row in rows
        ], dtype=float)
        target_trace = target / trace_coefficient
        polynomial_traces = {}
        for degree in range(1, 6):
            coefficients = np.polynomial.polynomial.polyfit(
                trace_lambda, target_trace, degree
            )
            fitted_trace = np.polynomial.polynomial.polyval(
                trace_lambda, coefficients
            )
            fitted_reaction = trace_coefficient * fitted_trace
            fit_errors = fitted_reaction - target
            polynomial_traces[str(degree)] = {
                "coefficients_low_to_high": [float(value) for value in coefficients],
                "max_abs_residual": float(np.max(np.abs(fit_errors))),
                "rms_residual": float(np.sqrt(np.mean(fit_errors * fit_errors))),
            }
        output["/".join(pair)] = {
            "node_count": len(rows),
            "transport_weight": float(weights[0]),
            "nontransport_weight": float(weights[1]),
            "max_abs_residual": float(np.max(np.abs(errors))),
            "rms_residual": float(np.sqrt(np.mean(errors * errors))),
            "worst": [
                {
                    "node": rows[index][0],
                    "target_reaction": rows[index][3],
                    "predicted_reaction": float(predicted[index]),
                    "residual": float(errors[index]),
                }
                for index in worst
            ],
            "affine_side_trace": {
                "transport_slope": float(affine_weights[0]),
                "nontransport_slope": float(affine_weights[1]),
                "transport_offset_V": float(affine_weights[2]),
                "nontransport_offset_V": float(affine_weights[3]),
                "max_abs_residual": float(np.max(np.abs(affine_errors))),
                "rms_residual": float(np.sqrt(np.mean(affine_errors * affine_errors))),
                "worst": [
                    {
                        "node": nonzero_rows[index][0],
                        "target_reaction": nonzero_rows[index][3],
                        "predicted_reaction": float(affine_predicted[index]),
                        "residual": float(affine_errors[index]),
                    }
                    for index in affine_worst
                ],
            },
            "quadratic_side_trace": {
                "transport_linear": float(quadratic_weights[0]),
                "nontransport_linear": float(quadratic_weights[1]),
                "transport_quadratic_per_V": float(quadratic_weights[2]),
                "nontransport_quadratic_per_V": float(quadratic_weights[3]),
                "max_abs_residual": float(np.max(np.abs(quadratic_errors))),
                "rms_residual": float(np.sqrt(np.mean(quadratic_errors * quadratic_errors))),
                "worst": [
                    {
                        "node": rows[index][0],
                        "target_reaction": rows[index][3],
                        "predicted_reaction": float(quadratic_predicted[index]),
                        "residual": float(quadratic_errors[index]),
                    }
                    for index in quadratic_worst
                ],
            },
            "single_eliminated_trace_polynomial": polynomial_traces,
            "native_material_side_trace": {
                "transport_weight": float(native_weights[0]),
                "nontransport_weight": float(native_weights[1]),
                "max_abs_residual": float(np.max(np.abs(native_errors))),
                "rms_residual": float(np.sqrt(np.mean(native_errors * native_errors))),
            },
        }
    if not args.group_by_region:
        columns = ["Si", "PolySilicon", "SiO2"]
        global_rows = []
        global_affine_rows = []
        global_targets = []
        global_nodes = []
        for node, material_set in materials.items():
            record = nodes[node]
            if (
                record["is_active"] != "1"
                or record["is_dirichlet"] != "0"
                or "SiO2" not in material_set
                or not (material_set & {"Si", "PolySilicon"})
            ):
                continue
            global_rows.append([
                material_reaction[node].get(column, 0.0) for column in columns
            ])
            output_lambda = float(record["output_lambda_V"])
            global_affine_rows.append(
                global_rows[-1] + [
                    value / output_lambda for value in global_rows[-1]
                ]
            )
            global_targets.append(-float(record["stiffness"]))
            global_nodes.append(node)
        global_matrix = np.asarray(global_rows, dtype=float)
        global_target = np.asarray(global_targets, dtype=float)
        global_weights, _, _, _ = np.linalg.lstsq(
            global_matrix, global_target, rcond=None
        )
        global_errors = global_matrix @ global_weights - global_target
        global_affine_matrix = np.asarray(global_affine_rows, dtype=float)
        global_affine_weights, _, _, _ = np.linalg.lstsq(
            global_affine_matrix, global_target, rcond=None
        )
        global_affine_errors = (
            global_affine_matrix @ global_affine_weights - global_target
        )
        output["_global_material_weights"] = {
            "columns": columns,
            "weights": [float(value) for value in global_weights],
            "max_abs_residual": float(np.max(np.abs(global_errors))),
            "rms_residual": float(np.sqrt(np.mean(global_errors * global_errors))),
            "affine_slopes": [float(value) for value in global_affine_weights[:3]],
            "affine_offsets_V": [float(value) for value in global_affine_weights[3:]],
            "affine_max_abs_residual": float(np.max(np.abs(global_affine_errors))),
            "affine_rms_residual": float(
                np.sqrt(np.mean(global_affine_errors * global_affine_errors))
            ),
        }
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
