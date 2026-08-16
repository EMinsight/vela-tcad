#!/usr/bin/env python3
"""Recover one Sentaurus Eq. 231 Jacobian row from NewtonPlot probes."""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


def read_scalar(export_dir: Path, field: str, region: int) -> dict[int, float]:
    manifest = json.loads(
        (export_dir / "field_manifest.json").read_text(encoding="utf-8")
    )
    matches = [
        entry for entry in manifest["fields"]
        if entry["name"] == field and int(entry["region"]) == region
    ]
    if len(matches) != 1:
        raise ValueError(
            f"expected one {field!r} entry for region {region}, got {len(matches)}"
        )
    with (export_dir / "fields" / matches[0]["csv_file"]).open(
        newline="", encoding="utf-8"
    ) as stream:
        return {
            int(row["node_id"]): float(row["component0"])
            for row in csv.DictReader(stream)
        }


def read_scalar_occurrences(export_dir: Path, field: str) -> dict[int, list[float]]:
    manifest = json.loads(
        (export_dir / "field_manifest.json").read_text(encoding="utf-8")
    )
    result: dict[int, list[float]] = defaultdict(list)
    for entry in manifest["fields"]:
        if entry["name"] != field:
            continue
        with (export_dir / "fields" / entry["csv_file"]).open(
            newline="", encoding="utf-8"
        ) as stream:
            for row in csv.DictReader(stream):
                result[int(row["node_id"])].append(float(row["component0"]))
    return dict(result)


def current_p1_row(
    mesh_path: Path,
    cells_path: Path,
    row_node: int,
    thermal_voltage: float,
    theta: float,
) -> tuple[dict[int, float], float, dict[str, float]]:
    mesh = json.loads(mesh_path.read_text(encoding="utf-8"))
    coordinates = {
        int(node["id"]): (float(node["x"]), float(node["y"]))
        for node in mesh["nodes"]
    }
    triangles = {
        int(cell["id"]): [int(node) for node in cell.get("node_ids", cell.get("nodes"))]
        for cell in mesh["triangles"]
    }
    by_cell: dict[int, list[dict[str, str]]] = defaultdict(list)
    with cells_path.open(newline="", encoding="utf-8") as stream:
        for record in csv.DictReader(stream):
            by_cell[int(record["cell_id"])].append(record)

    jacobian: dict[int, float] = defaultdict(float)
    residual_parts = defaultdict(float)
    for cell_id, records in by_cell.items():
        row_records = [
            record for record in records
            if int(record["node_id"]) == row_node
            and record["is_active"] == "1"
            and record["is_dirichlet"] == "0"
        ]
        if not row_records:
            continue
        row_record = row_records[0]
        nodes = triangles[cell_id]
        local_row = int(row_record["local_node"])
        records_by_local = {int(record["local_node"]): record for record in records}
        points = [coordinates[node] for node in nodes]
        x = [point[0] for point in points]
        y = [point[1] for point in points]
        twice_area = abs(
            (x[1] - x[0]) * (y[2] - y[0])
            - (x[2] - x[0]) * (y[1] - y[0])
        )
        area = 0.5 * twice_area
        b = [y[1] - y[2], y[2] - y[0], y[0] - y[1]]
        c = [x[2] - x[1], x[0] - x[2], x[1] - x[0]]
        w = [float(records_by_local[local]["w"]) for local in range(3)]
        gradient_x = sum(w[local] * b[local] for local in range(3)) / twice_area
        gradient_y = sum(w[local] * c[local] for local in range(3)) / twice_area

        for local, node in enumerate(nodes):
            stiffness = (
                b[local_row] * b[local] + c[local_row] * c[local]
            ) / (2.0 * twice_area)
            gradient_dot_shape = (
                gradient_x * b[local] + gradient_y * c[local]
            ) / twice_area
            jacobian[node] += stiffness / thermal_voltage
            jacobian[node] += (
                theta * (-2.0 / thermal_voltage)
                * gradient_dot_shape * area / 3.0
            )
        reaction_diagonal = (
            float(row_record["reaction"]) / float(row_record["lambda_V"])
        )
        jacobian[row_node] += reaction_diagonal
        for key in ("stiffness", "gradient_squared", "reaction", "total"):
            residual_parts[key] += float(row_record[key])

    return dict(jacobian), sum(residual_parts.values()) - residual_parts["total"], dict(residual_parts)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--cells", type=Path, required=True)
    parser.add_argument("--baseline-export", type=Path, required=True)
    parser.add_argument("--theta-zero-export", type=Path, required=True)
    parser.add_argument("--perturbation-root", type=Path, required=True)
    parser.add_argument("--perturbation-prefix", default="perturbation_")
    parser.add_argument("--row-node", type=int, required=True)
    parser.add_argument("--region", type=int, required=True)
    parser.add_argument("--perturb-nodes", type=int, nargs="+", required=True)
    parser.add_argument("--delta", type=float, default=1.0e-5)
    parser.add_argument("--thermal-voltage", type=float, default=0.025851999786)
    parser.add_argument("--theta", type=float, default=0.5)
    args = parser.parse_args()

    baseline = read_scalar(
        args.baseline_export, "eQuantumPotentialRhs", args.region
    )[args.row_node]
    baseline_state = read_scalar_occurrences(
        args.baseline_export, "eQuantumPotential"
    )
    theta_zero_rhs = read_scalar(
        args.theta_zero_export, "eQuantumPotentialRhs", args.region
    )[args.row_node]
    vela_jacobian, _, residual_parts = current_p1_row(
        args.mesh, args.cells, args.row_node,
        args.thermal_voltage, args.theta,
    )
    theta_zero_residual = residual_parts["stiffness"] + residual_parts["reaction"]
    sentaurus_to_vela_scale = theta_zero_residual / theta_zero_rhs

    rows = []
    for node in args.perturb_nodes:
        export_dir = (
            args.perturbation_root
            / f"{args.perturbation_prefix}{node}_export"
        )
        perturbed = read_scalar(
            export_dir, "eQuantumPotentialRhs", args.region
        )[args.row_node]
        perturbed_state = read_scalar_occurrences(
            export_dir, "eQuantumPotential"
        )
        state_deltas = [
            probe - reference
            for probe in perturbed_state[node]
            for reference in baseline_state[node]
            if abs(probe - reference) > 0.0
        ]
        observed_delta = max(state_deltas, key=abs) if state_deltas else args.delta
        sentaurus_normalized = (perturbed - baseline) / observed_delta
        sentaurus_scaled = sentaurus_normalized * sentaurus_to_vela_scale
        vela = vela_jacobian.get(node, 0.0)
        rows.append({
            "column_node": node,
            "observed_delta_V": observed_delta,
            "sentaurus_normalized_derivative": sentaurus_normalized,
            "sentaurus_scaled_derivative": sentaurus_scaled,
            "vela_p1_derivative": vela,
            "vela_minus_sentaurus": vela - sentaurus_scaled,
        })

    print(json.dumps({
        "row_node": args.row_node,
        "region": args.region,
        "delta_V": args.delta,
        "baseline_rhs": baseline,
        "theta_zero_rhs": theta_zero_rhs,
        "theta_zero_vela_residual": theta_zero_residual,
        "sentaurus_to_vela_scale": sentaurus_to_vela_scale,
        "fixed_state_vela_residual_parts": residual_parts,
        "jacobian_row": rows,
    }, indent=2))


if __name__ == "__main__":
    main()
