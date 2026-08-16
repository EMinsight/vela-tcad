#!/usr/bin/env python3
"""Replay Eq. 231 source terms with element-vertex box measures."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path


def triangle_area(a: dict, b: dict, c: dict) -> float:
    return 0.5 * abs(
        (b["x"] - a["x"]) * (c["y"] - a["y"])
        - (b["y"] - a["y"]) * (c["x"] - a["x"])
    )


def angle_degrees(a: dict, b: dict, c: dict) -> float:
    ux, uy = b["x"] - a["x"], b["y"] - a["y"]
    vx, vy = c["x"] - a["x"], c["y"] - a["y"]
    norm = math.hypot(ux, uy) * math.hypot(vx, vy)
    if norm == 0.0:
        return 0.0
    cosine = max(-1.0, min(1.0, (ux * vx + uy * vy) / norm))
    return math.degrees(math.acos(cosine))


def cotangent(a: dict, b: dict, opposite: dict) -> float:
    ux, uy = a["x"] - opposite["x"], a["y"] - opposite["y"]
    vx, vy = b["x"] - opposite["x"], b["y"] - opposite["y"]
    cross = ux * vy - uy * vx
    if abs(cross) < 1.0e-30:
        return 0.0
    return (ux * vx + uy * vy) / abs(cross)


def distance_squared(a: dict, b: dict) -> float:
    return (b["x"] - a["x"]) ** 2 + (b["y"] - a["y"]) ** 2


def mixed_voronoi_shares(nodes: list[dict]) -> list[float]:
    area = triangle_area(*nodes)
    angles = [
        angle_degrees(nodes[0], nodes[1], nodes[2]),
        angle_degrees(nodes[1], nodes[2], nodes[0]),
        angle_degrees(nodes[2], nodes[0], nodes[1]),
    ]
    maximum = max(range(3), key=angles.__getitem__)
    if angles[maximum] > 90.0:
        shares = [0.25 * area] * 3
        shares[maximum] = 0.5 * area
        return shares

    shares: list[float] = []
    for i in range(3):
        j, k = (i + 1) % 3, (i + 2) % 3
        cot_j = cotangent(nodes[i], nodes[k], nodes[j])
        cot_k = cotangent(nodes[i], nodes[j], nodes[k])
        shares.append(
            0.125
            * (
                distance_squared(nodes[i], nodes[k]) * cot_j
                + distance_squared(nodes[i], nodes[j]) * cot_k
            )
        )
    return shares


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--cells", type=Path, required=True)
    parser.add_argument("--node-diagnostics", type=Path)
    parser.add_argument("--top", type=int, default=20)
    parser.add_argument("--inspect-node", type=int)
    parser.add_argument(
        "--quadratic",
        choices=(
            "cell_p1",
            "edge_energy_raw",
            "edge_energy_positive",
            "node_region_gradient_barycentric",
            "node_region_gradient_box",
            "node_region_gradient_component_minmod",
            "node_region_gradient_min_norm",
            "box_vector_scheme_a",
            "box_vector_scheme_b",
        ),
        default="cell_p1",
    )
    parser.add_argument("--theta", type=float, default=0.5)
    parser.add_argument(
        "--reaction-trace",
        choices=(
            "material_side",
            "owner_nodal",
            "owner_transport_interface",
            "maximum_material_side",
        ),
        default="material_side",
    )
    args = parser.parse_args()

    owner_lambda: dict[int, float] = {}
    if args.reaction_trace in ("owner_nodal", "owner_transport_interface"):
        if args.node_diagnostics is None:
            raise ValueError("--node-diagnostics is required for owner_nodal")
        with args.node_diagnostics.open(newline="", encoding="utf-8") as stream:
            owner_lambda = {
                int(row["node_id"]): float(row["output_lambda_V"])
                for row in csv.DictReader(stream)
            }

    mesh = json.loads(args.mesh.read_text(encoding="utf-8"))
    nodes = {int(item["id"]): item for item in mesh["nodes"]}
    triangles = {int(item["id"]): item for item in mesh["triangles"]}
    shares_by_cell: dict[int, list[float]] = {}
    maximum_angle_by_cell: dict[int, float] = {}
    for cell_id, cell in triangles.items():
        cell_nodes = cell.get("nodes", cell.get("node_ids"))
        cell_points = [nodes[int(node)] for node in cell_nodes]
        shares_by_cell[cell_id] = mixed_voronoi_shares(cell_points)
        maximum_angle_by_cell[cell_id] = max(
            angle_degrees(cell_points[i], cell_points[(i + 1) % 3],
                          cell_points[(i + 2) % 3])
            for i in range(3)
        )

    assembled = defaultdict(
        lambda: {
            "stiffness": 0.0,
            "quadratic": 0.0,
            "reaction": 0.0,
            "baseline_quadratic": 0.0,
            "baseline_reaction": 0.0,
        }
    )
    metadata: dict[int, dict[str, str]] = {}
    rows_by_cell: dict[int, list[dict[str, str]]] = defaultdict(list)
    with args.cells.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            rows_by_cell[int(row["cell_id"])].append(row)
    node_transport_flags: dict[int, set[bool]] = defaultdict(set)
    node_region_lambdas: dict[int, dict[int, float]] = defaultdict(dict)
    for cell_rows in rows_by_cell.values():
        for row in cell_rows:
            node = int(row["node_id"])
            node_transport_flags[node].add(
                row["is_transport"] == "1"
            )
            node_region_lambdas[node][int(row["region_id"])] = float(
                row["lambda_V"]
            )

    # Reconstruct a material-side nodal gradient from adjacent element
    # gradients. Keeping region_id in the key prevents averaging a material
    # driving field across an Si/SiO2 Formula 0 interface.
    gradient_sum = defaultdict(lambda: [0.0, 0.0, 0.0])
    gradient_samples: dict[tuple[int, int], list[tuple[float, float]]] = defaultdict(list)
    box_volume = defaultdict(float)
    scheme_a_numerator = defaultdict(lambda: [0.0, 0.0])
    scheme_b_matrix = defaultdict(lambda: [0.0, 0.0, 0.0])
    scheme_b_rhs = defaultdict(lambda: [0.0, 0.0])
    for cell_id, cell_rows in rows_by_cell.items():
        cell_rows.sort(key=lambda row: int(row["local_node"]))
        triangle = triangles[cell_id]
        triangle_nodes = triangle.get("nodes", triangle.get("node_ids"))
        points = [nodes[int(index)] for index in triangle_nodes]
        area_internal = triangle_area(*points)
        x = [point["x"] for point in points]
        y = [point["y"] for point in points]
        b = [y[1] - y[2], y[2] - y[0], y[0] - y[1]]
        c = [x[2] - x[1], x[0] - x[2], x[1] - x[0]]
        w = [float(row["w"]) for row in cell_rows]
        four_area_squared = 4.0 * area_internal * area_internal
        local_stiffness = [
            [
                area_internal
                * (b[a] * b[d] + c[a] * c[d])
                / four_area_squared
                for d in range(3)
            ]
            for a in range(3)
        ]
        gx = sum(w[local] * b[local] for local in range(3)) / (2.0 * area_internal)
        gy = sum(w[local] * c[local] for local in range(3)) / (2.0 * area_internal)
        for local, row in enumerate(cell_rows):
            if args.quadratic == "node_region_gradient_box":
                weight = shares_by_cell[cell_id][local]
            else:
                weight = area_internal / 3.0
            key = (int(row["node_id"]), int(row["region_id"]))
            gradient_sum[key][0] += weight * gx
            gradient_sum[key][1] += weight * gy
            gradient_sum[key][2] += weight
            gradient_samples[key].append((gx, gy))
            box_volume[key] += shares_by_cell[cell_id][local]
        for first in range(3):
            for second in range(first + 1, 3):
                conductance = max(-local_stiffness[first][second], 0.0)
                if conductance == 0.0:
                    continue
                dx = x[second] - x[first]
                dy = y[second] - y[first]
                distance = math.hypot(dx, dy)
                if distance == 0.0:
                    continue
                xi_x, xi_y = dx / distance, dy / distance
                delta = w[second] - w[first]
                for local, sign in ((first, 1.0), (second, -1.0)):
                    row = cell_rows[local]
                    key = (int(row["node_id"]), int(row["region_id"]))
                    local_dx = sign * dx
                    local_dy = sign * dy
                    local_delta = sign * delta
                    scheme_a_numerator[key][0] += (
                        conductance * local_dx * local_delta
                    )
                    scheme_a_numerator[key][1] += (
                        conductance * local_dy * local_delta
                    )
                    scheme_b_matrix[key][0] += conductance * xi_x * xi_x
                    scheme_b_matrix[key][1] += conductance * xi_x * xi_y
                    scheme_b_matrix[key][2] += conductance * xi_y * xi_y
                    projected = local_delta / distance
                    scheme_b_rhs[key][0] += conductance * (sign * xi_x) * projected
                    scheme_b_rhs[key][1] += conductance * (sign * xi_y) * projected

    node_region_gradient = {
        key: (value[0] / value[2], value[1] / value[2])
        for key, value in gradient_sum.items()
        if value[2] > 0.0
    }
    if args.quadratic == "node_region_gradient_component_minmod":
        for key, samples in gradient_samples.items():
            components = []
            for axis in range(2):
                values = [sample[axis] for sample in samples]
                if all(value > 0.0 for value in values):
                    components.append(min(values))
                elif all(value < 0.0 for value in values):
                    components.append(max(values))
                else:
                    components.append(0.0)
            node_region_gradient[key] = (components[0], components[1])
    elif args.quadratic == "node_region_gradient_min_norm":
        node_region_gradient.update({
            key: min(samples, key=lambda value: value[0] ** 2 + value[1] ** 2)
            for key, samples in gradient_samples.items()
        })
    elif args.quadratic == "box_vector_scheme_a":
        node_region_gradient.update({
            key: (
                numerator[0] / (2.0 * box_volume[key]),
                numerator[1] / (2.0 * box_volume[key]),
            )
            for key, numerator in scheme_a_numerator.items()
            if box_volume[key] > 0.0
        })
    elif args.quadratic == "box_vector_scheme_b":
        for key, matrix in scheme_b_matrix.items():
            xx, xy, yy = matrix
            determinant = xx * yy - xy * xy
            if determinant <= 1.0e-30:
                continue
            rhs_x, rhs_y = scheme_b_rhs[key]
            node_region_gradient[key] = (
                (yy * rhs_x - xy * rhs_y) / determinant,
                (-xy * rhs_x + xx * rhs_y) / determinant,
            )

    for cell_id, cell_rows in rows_by_cell.items():
        cell_rows.sort(key=lambda row: int(row["local_node"]))
        triangle = triangles[cell_id]
        triangle_nodes = triangle.get("nodes", triangle.get("node_ids"))
        points = [nodes[int(index)] for index in triangle_nodes]
        area_internal = triangle_area(*points)
        x = [point["x"] for point in points]
        y = [point["y"] for point in points]
        b = [y[1] - y[2], y[2] - y[0], y[0] - y[1]]
        c = [x[2] - x[1], x[0] - x[2], x[1] - x[0]]
        four_area_squared = 4.0 * area_internal * area_internal
        stiffness = [
            [
                area_internal
                * (b[a] * b[d] + c[a] * c[d])
                / four_area_squared
                for d in range(3)
            ]
            for a in range(3)
        ]
        w = [float(row["w"]) for row in cell_rows]
        for row in cell_rows:
            if row["is_active"] != "1" or row["is_dirichlet"] == "1":
                continue
            node = int(row["node_id"])
            local = int(row["local_node"])
            share_ratio = shares_by_cell[cell_id][local] / (area_internal / 3.0)
            if args.quadratic == "cell_p1":
                quadratic = float(row["gradient_squared"]) * share_ratio
            elif (
                args.quadratic.startswith("node_region_gradient_")
                or args.quadratic.startswith("box_vector_scheme_")
            ):
                gx, gy = node_region_gradient[(node, int(row["region_id"]))]
                quadratic = (
                    args.theta
                    * (gx * gx + gy * gy)
                    * shares_by_cell[cell_id][local]
                )
            else:
                quadratic = 0.0
                for other in range(3):
                    if other == local:
                        continue
                    conductance = -stiffness[local][other]
                    if args.quadratic == "edge_energy_positive":
                        conductance = max(conductance, 0.0)
                    quadratic += (
                        args.theta
                        * 0.5
                        * conductance
                        * (w[other] - w[local]) ** 2
                    )
            assembled[node]["stiffness"] += float(row["stiffness"])
            assembled[node]["quadratic"] += quadratic
            reaction = float(row["reaction"])
            use_owner_trace = args.reaction_trace == "owner_nodal" or (
                args.reaction_trace == "owner_transport_interface"
                and node_transport_flags[node] == {False, True}
            )
            if use_owner_trace:
                material_lambda = float(row["lambda_V"])
                if abs(material_lambda) > 1.0e-30:
                    reaction *= owner_lambda[node] / material_lambda
            elif (
                args.reaction_trace == "maximum_material_side"
                and len(node_region_lambdas[node]) > 1
            ):
                material_lambda = float(row["lambda_V"])
                if abs(material_lambda) > 1.0e-30:
                    reaction *= (
                        max(node_region_lambdas[node].values())
                        / material_lambda
                    )
            assembled[node]["reaction"] += reaction * share_ratio
            assembled[node]["baseline_quadratic"] += float(row["gradient_squared"])
            assembled[node]["baseline_reaction"] += float(row["reaction"])
            metadata[node] = {
                "region_name": row["region_name"],
                "material": row["material"],
            }

    ranked = []
    for node, terms in assembled.items():
        total = terms["stiffness"] + terms["quadratic"] + terms["reaction"]
        baseline_total = (
            terms["stiffness"]
            + terms["baseline_quadratic"]
            + terms["baseline_reaction"]
        )
        ranked.append(
            {
                "node_id": node,
                **metadata[node],
                **terms,
                "total": total,
                "absolute_total": abs(total),
                "baseline_total": baseline_total,
            }
        )
    ranked.sort(key=lambda row: row["absolute_total"], reverse=True)
    output = {
        "geometry": {
            "triangle_count": len(triangles),
            "obtuse_triangle_count": sum(
                angle > 90.0 for angle in maximum_angle_by_cell.values()
            ),
            "maximum_angle_deg": max(maximum_angle_by_cell.values()),
        },
        "maximum": ranked[0],
        "top": ranked[: args.top],
    }
    if args.inspect_node is not None:
        inspected_terms = assembled.get(args.inspect_node)
        if inspected_terms is not None:
            inspected_total = (
                inspected_terms["stiffness"]
                + inspected_terms["quadratic"]
                + inspected_terms["reaction"]
            )
            output["inspection_assembled"] = {
                **inspected_terms,
                "total": inspected_total,
            }
        output["inspection"] = {
            str(key): {
                "cell_gradients": gradient_samples[key],
                "selected_gradient": node_region_gradient[key],
            }
            for key in gradient_samples
            if key[0] == args.inspect_node
        }
        output["inspection_cells"] = [
            {
                "cell_id": cell_id,
                "maximum_angle_deg": maximum_angle_by_cell[cell_id],
                "mixed_shares": shares_by_cell[cell_id],
            }
            for cell_id, rows in rows_by_cell.items()
            if any(int(row["node_id"]) == args.inspect_node for row in rows)
        ]
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
