#!/usr/bin/env python3
"""Audit carrier-specific Eparallel recovery on selected BVmethods NMOS edges.

The script reproduces Vela's area-weighted adjacent-cell electric field and
length-weighted nodal least-squares SG current direction.  It then compares
the resulting edge Eparallel with Sentaurus nodal vectors and with the field
inverted from Sentaurus eAlphaAvalanche using the configured van Overstraeten
electron coefficient.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Any

import audit_bvmethods_nmos_hotspot_slope_and_criteria as base


TARGET_EDGES = (2663, 2666, 2672, 2675, 2681)
DEFAULT_MESH = base.RUN / "vela/mesh.json"
DEFAULT_VELA_RUNS = {
    6.4: base.RUN / "vela_validation/btbt_e2_iic_qf_vector_fixed6p4_rerun_20260805",
    7.0: base.RUN / "vela_validation/btbt_e2_iic_qf_vector_fixed7p0_20260805",
}
DEFAULT_IMPLEMENTED_RUNS = {
    6.4: base.RUN / "vela_validation/btbt_e2_iic_qf_vector_nodal_vertex_star_fixed6p4_20260805",
    7.0: base.RUN / "vela_validation/btbt_e2_iic_qf_vector_nodal_vertex_star_fixed7p0_20260805",
}
DEFAULT_OUT = base.RUN / "vela_validation/qf_vector_eparallel_edge_interpolation_20260805"

# ImpactIonizationModel defaults used by the validation configuration at 300 K.
VAN_OVERSTRAETEN_ELECTRON_A_M_INV = 7.03e7
VAN_OVERSTRAETEN_ELECTRON_B_V_M = 1.231e8


Vector = tuple[float, float]


def add(a: Vector, b: Vector) -> Vector:
    return a[0] + b[0], a[1] + b[1]


def scale(factor: float, value: Vector) -> Vector:
    return factor * value[0], factor * value[1]


def dot(a: Vector, b: Vector) -> float:
    return a[0] * b[0] + a[1] * b[1]


def norm(value: Vector) -> float:
    return math.hypot(value[0], value[1])


def eparallel(electric_field: Vector, current: Vector) -> float:
    magnitude = norm(current)
    return max(dot(electric_field, current) / magnitude, 0.0) if magnitude > 0.0 else 0.0


def alpha_from_field(field_v_m: float) -> float:
    if field_v_m <= 0.0:
        return 0.0
    return VAN_OVERSTRAETEN_ELECTRON_A_M_INV * math.exp(
        -VAN_OVERSTRAETEN_ELECTRON_B_V_M / field_v_m
    )


def field_from_alpha(alpha_m_inv: float) -> float:
    if not 0.0 < alpha_m_inv < VAN_OVERSTRAETEN_ELECTRON_A_M_INV:
        return math.nan
    return -VAN_OVERSTRAETEN_ELECTRON_B_V_M / math.log(
        alpha_m_inv / VAN_OVERSTRAETEN_ELECTRON_A_M_INV
    )


def number(row: dict[str, str], key: str) -> float:
    return float(row.get(key, "0") or 0.0)


def read_state(path: Path) -> dict[int, float]:
    return {
        int(row["node_id"]): float(row["psi"])
        for row in base.csv_rows(path)
    }


def triangle_gradient(
    node_ids: list[int], nodes: dict[int, Vector], values: dict[int, float]
) -> tuple[Vector, float]:
    n0, n1, n2 = node_ids
    x0, y0 = nodes[n0]
    x1, y1 = nodes[n1]
    x2, y2 = nodes[n2]
    dx10, dy10 = x1 - x0, y1 - y0
    dx20, dy20 = x2 - x0, y2 - y0
    determinant = dx10 * dy20 - dx20 * dy10
    if abs(determinant) <= 1.0e-30:
        return (0.0, 0.0), 0.0
    dv10 = values[n1] - values[n0]
    dv20 = values[n2] - values[n0]
    gradient = (
        (dv10 * dy20 - dv20 * dy10) / determinant,
        (dx10 * dv20 - dx20 * dv10) / determinant,
    )
    return gradient, 0.5 * abs(determinant)


def nodal_least_squares_current(
    node: int,
    node_edges: dict[int, list[dict[str, str]]],
) -> Vector:
    a00 = a01 = a11 = b0 = b1 = 0.0
    fallback = (0.0, 0.0)
    fallback_weight = 0.0
    used = 0
    for row in node_edges.get(node, []):
        if number(row, "electron_mobility_m2_V_s") <= 0.0:
            continue
        length = math.hypot(
            number(row, "x1_um") - number(row, "x0_um"),
            number(row, "y1_um") - number(row, "y0_um"),
        )
        if length <= 1.0e-30:
            continue
        tangent = (
            (number(row, "x1_um") - number(row, "x0_um")) / length,
            (number(row, "y1_um") - number(row, "y0_um")) / length,
        )
        flux = number(row, "electron_raw_signed_flux_proxy")
        weight = length
        a00 += weight * tangent[0] * tangent[0]
        a01 += weight * tangent[0] * tangent[1]
        a11 += weight * tangent[1] * tangent[1]
        b0 += weight * tangent[0] * flux
        b1 += weight * tangent[1] * flux
        fallback = add(fallback, scale(weight * flux, tangent))
        fallback_weight += weight
        used += 1
    determinant = a00 * a11 - a01 * a01
    determinant_scale = max(abs(a00 * a11), abs(a01 * a01), 1.0e-300)
    if used >= 2 and abs(determinant) > 1.0e-24 * determinant_scale:
        return (
            (b0 * a11 - b1 * a01) / determinant,
            (a00 * b1 - a01 * b0) / determinant,
        )
    return scale(1.0 / fallback_weight, fallback) if fallback_weight > 0.0 else (0.0, 0.0)


def mesh_context(mesh_path: Path) -> dict[str, Any]:
    mesh = json.loads(mesh_path.read_text(encoding="utf-8"))
    nodes = {int(row["id"]): (float(row["x"]), float(row["y"])) for row in mesh["nodes"]}
    regions = {int(row["id"]): row for row in mesh["regions"]}
    triangles = {
        int(row["id"]): {
            "nodes": [int(node) for node in row["node_ids"]],
            "region_id": int(row["region_id"]),
            "material": regions[int(row["region_id"])]["material"],
        }
        for row in mesh["triangles"]
    }
    pair_cells: dict[tuple[int, int], list[int]] = defaultdict(list)
    node_cells: dict[int, list[int]] = defaultdict(list)
    for cell_id, cell in triangles.items():
        ids = cell["nodes"]
        for node in ids:
            node_cells[node].append(cell_id)
        for left, right in ((ids[0], ids[1]), (ids[1], ids[2]), (ids[2], ids[0])):
            pair_cells[tuple(sorted((left, right)))].append(cell_id)
    return {
        "nodes": nodes,
        "triangles": triangles,
        "pair_cells": pair_cells,
        "node_cells": node_cells,
    }


def sentaurus_vectors(state: Path) -> dict[str, dict[int, Vector] | dict[int, float]]:
    electric_raw = base.field(state, "ElectricField")
    current_raw = base.field(state, "eCurrentDensity")
    alpha_raw = base.field(state, "eAlphaAvalanche")
    potential_raw = base.field(state, "ElectrostaticPotential")
    electric = {
        node: (values[0] * 100.0, values[1] * 100.0)
        for node, values in electric_raw.items()
    }
    current = {
        node: (values[0] * 1.0e4, values[1] * 1.0e4)
        for node, values in current_raw.items()
    }
    alpha = {node: abs(values[0]) * 100.0 for node, values in alpha_raw.items()}
    potential = {node: values[0] for node, values in potential_raw.items()}
    return {
        "electric": electric,
        "current": current,
        "alpha": alpha,
        "potential": potential,
    }


def audit_bias(
    bias: float,
    run: Path,
    sent_state: Path,
    mesh: dict[str, Any],
    target_edges: tuple[int, ...],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, Any]]]:
    rows = [
        row for row in base.csv_rows(run / "sg_avalanche_edges.csv")
        if math.isclose(number(row, "bias_V"), bias, abs_tol=1.0e-9)
    ]
    edge_by_id = {int(row["edge_id"]): row for row in rows}
    node_edges: dict[int, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        node_edges[int(row["node0"])].append(row)
        node_edges[int(row["node1"])].append(row)
    nodal_particle_current = {
        node: nodal_least_squares_current(node, node_edges)
        for node in mesh["nodes"]
    }
    psi = read_state(run / "last_state.csv")
    cell_data: dict[int, tuple[Vector, float]] = {
        cell_id: triangle_gradient(cell["nodes"], mesh["nodes"], psi)
        for cell_id, cell in mesh["triangles"].items()
    }
    sent = sentaurus_vectors(sent_state)
    sent_cell_data: dict[int, tuple[Vector, float]] = {}
    for cell_id, cell in mesh["triangles"].items():
        if all(node in sent["potential"] for node in cell["nodes"]):
            sent_cell_data[cell_id] = triangle_gradient(
                cell["nodes"], mesh["nodes"], sent["potential"]
            )

    def nodal_transport_electric(
        node: int, data: dict[int, tuple[Vector, float]]
    ) -> Vector:
        weighted = (0.0, 0.0)
        total_area = 0.0
        for cell_id in mesh["node_cells"][node]:
            if mesh["triangles"][cell_id]["material"].lower() not in {"si", "silicon"}:
                continue
            if cell_id not in data:
                continue
            gradient, area_um2 = data[cell_id]
            weighted = add(weighted, scale(area_um2, gradient))
            total_area += area_um2
        return scale(-1.0e6 / total_area, weighted) if total_area > 0.0 else (0.0, 0.0)
    edge_records: list[dict[str, Any]] = []
    cell_records: list[dict[str, Any]] = []
    stencil_records: list[dict[str, Any]] = []

    for edge_id in target_edges:
        row = edge_by_id[edge_id]
        n0, n1 = int(row["node0"]), int(row["node1"])
        cells = mesh["pair_cells"][tuple(sorted((n0, n1)))]
        weighted_gradient = (0.0, 0.0)
        total_area = 0.0
        for cell_id in cells:
            gradient, area_um2 = cell_data[cell_id]
            weighted_gradient = add(weighted_gradient, scale(area_um2, gradient))
            total_area += area_um2
        weighted_gradient = scale(1.0 / total_area, weighted_gradient)
        vela_electric = scale(-1.0e6, weighted_gradient)
        vela_j0 = scale(-1.0, nodal_particle_current[n0])
        vela_j1 = scale(-1.0, nodal_particle_current[n1])
        vela_current = scale(0.5, add(vela_j0, vela_j1))
        vela_recovered_eparallel = eparallel(vela_electric, vela_current)
        vela_nodal_e0 = nodal_transport_electric(n0, cell_data)
        vela_nodal_e1 = nodal_transport_electric(n1, cell_data)
        vela_nodal_electric = scale(0.5, add(vela_nodal_e0, vela_nodal_e1))
        vela_nodal_eparallel = eparallel(vela_nodal_electric, vela_current)
        vela_cos = (
            dot(vela_electric, vela_current) / (norm(vela_electric) * norm(vela_current))
            if norm(vela_electric) > 0.0 and norm(vela_current) > 0.0 else 0.0
        )

        sent_e0 = sent["electric"][n0]
        sent_e1 = sent["electric"][n1]
        sent_j0 = sent["current"][n0]
        sent_j1 = sent["current"][n1]
        sent_alpha0 = sent["alpha"][n0]
        sent_alpha1 = sent["alpha"][n1]
        sent_direct0 = eparallel(sent_e0, sent_j0)
        sent_direct1 = eparallel(sent_e1, sent_j1)
        sent_inverse0 = field_from_alpha(sent_alpha0)
        sent_inverse1 = field_from_alpha(sent_alpha1)
        sent_electric = scale(0.5, add(sent_e0, sent_e1))
        sent_current = scale(0.5, add(sent_j0, sent_j1))
        sent_vector_first = eparallel(sent_electric, sent_current)
        sent_p1_weighted_gradient = (0.0, 0.0)
        sent_p1_total_area = 0.0
        for cell_id in cells:
            gradient, area_um2 = sent_cell_data[cell_id]
            sent_p1_weighted_gradient = add(
                sent_p1_weighted_gradient, scale(area_um2, gradient)
            )
            sent_p1_total_area += area_um2
        sent_p1_electric = scale(
            -1.0e6 / sent_p1_total_area, sent_p1_weighted_gradient
        )
        sent_p1_eparallel = eparallel(sent_p1_electric, sent_current)
        sent_p1_nodal_e0 = nodal_transport_electric(n0, sent_cell_data)
        sent_p1_nodal_e1 = nodal_transport_electric(n1, sent_cell_data)
        sent_p1_nodal_electric = scale(0.5, add(sent_p1_nodal_e0, sent_p1_nodal_e1))
        sent_p1_nodal_eparallel = eparallel(sent_p1_nodal_electric, sent_current)
        sent_cos = (
            dot(sent_electric, sent_current) / (norm(sent_electric) * norm(sent_current))
            if norm(sent_electric) > 0.0 and norm(sent_current) > 0.0 else 0.0
        )
        sent_alpha_mean = 0.5 * (sent_alpha0 + sent_alpha1)

        projected_cells: list[float] = []
        for cell_id in cells:
            gradient, area_um2 = cell_data[cell_id]
            cell_electric = scale(-1.0e6, gradient)
            cell_eparallel = eparallel(cell_electric, vela_current)
            sent_gradient, _ = sent_cell_data[cell_id]
            sent_cell_electric = scale(-1.0e6, sent_gradient)
            sent_cell_eparallel = eparallel(sent_cell_electric, sent_current)
            projected_cells.append(cell_eparallel)
            cell_records.append({
                "bias_V": bias,
                "edge_id": edge_id,
                "cell_id": cell_id,
                "material": mesh["triangles"][cell_id]["material"],
                "node_ids": ";".join(str(node) for node in mesh["triangles"][cell_id]["nodes"]),
                "area_um2": area_um2,
                "electric_x_V_m": cell_electric[0],
                "electric_y_V_m": cell_electric[1],
                "electric_magnitude_V_m": norm(cell_electric),
                "projected_eparallel_V_m": cell_eparallel,
                "sentaurus_p1_electric_x_V_m": sent_cell_electric[0],
                "sentaurus_p1_electric_y_V_m": sent_cell_electric[1],
                "sentaurus_p1_electric_magnitude_V_m": norm(sent_cell_electric),
                "sentaurus_p1_projected_eparallel_V_m": sent_cell_eparallel,
            })

        for endpoint in (n0, n1):
            for incident in node_edges[endpoint]:
                stencil_records.append({
                    "bias_V": bias,
                    "target_edge_id": edge_id,
                    "endpoint_node": endpoint,
                    "incident_edge_id": int(incident["edge_id"]),
                    "incident_node0": int(incident["node0"]),
                    "incident_node1": int(incident["node1"]),
                    "active_electron_edge": number(incident, "electron_mobility_m2_V_s") > 0.0,
                    "edge_length_m": number(incident, "edge_length_m"),
                    "electron_raw_signed_flux_proxy": number(incident, "electron_raw_signed_flux_proxy"),
                })

        csv_eparallel = number(row, "electron_impact_field_V_per_m")
        vela_alpha = number(row, "electron_alpha_m_inv")
        edge_records.append({
            "bias_V": bias,
            "edge_id": edge_id,
            "node0": n0,
            "node1": n1,
            "x_mid_um": 0.5 * (number(row, "x0_um") + number(row, "x1_um")),
            "y_mid_um": 0.5 * (number(row, "y0_um") + number(row, "y1_um")),
            "adjacent_cell_ids": ";".join(str(cell) for cell in cells),
            "vela_csv_eparallel_V_m": csv_eparallel,
            "vela_recovered_eparallel_V_m": vela_recovered_eparallel,
            "vela_recovery_relative_error": (
                vela_recovered_eparallel / csv_eparallel - 1.0 if csv_eparallel > 0.0 else math.nan
            ),
            "vela_electric_magnitude_V_m": norm(vela_electric),
            "vela_electric_x_V_m": vela_electric[0],
            "vela_electric_y_V_m": vela_electric[1],
            "vela_current_x_proxy": vela_current[0],
            "vela_current_y_proxy": vela_current[1],
            "vela_electric_current_cosine": vela_cos,
            "vela_nodal_p1_electric_magnitude_V_m": norm(vela_nodal_electric),
            "vela_nodal_p1_eparallel_V_m": vela_nodal_eparallel,
            "vela_nodal_p1_alpha_m_inv": alpha_from_field(vela_nodal_eparallel),
            "vela_edge_over_nodal_p1_eparallel": (
                vela_recovered_eparallel / vela_nodal_eparallel
                if vela_nodal_eparallel > 0.0 else math.nan
            ),
            "vela_adjacent_cell_projected_min_V_m": min(projected_cells),
            "vela_adjacent_cell_projected_max_V_m": max(projected_cells),
            "vela_alpha_csv_m_inv": vela_alpha,
            "vela_alpha_from_recovered_eparallel_m_inv": alpha_from_field(vela_recovered_eparallel),
            "vela_alpha_formula_relative_error": (
                alpha_from_field(vela_recovered_eparallel) / vela_alpha - 1.0
                if vela_alpha > 0.0 else math.nan
            ),
            "vela_electron_flux_proxy": number(row, "electron_flux_proxy"),
            "vela_electron_current_vector_magnitude_proxy": norm(vela_current),
            "sentaurus_node0_direct_eparallel_V_m": sent_direct0,
            "sentaurus_node1_direct_eparallel_V_m": sent_direct1,
            "sentaurus_endpoint_mean_direct_eparallel_V_m": 0.5 * (sent_direct0 + sent_direct1),
            "sentaurus_vector_first_eparallel_V_m": sent_vector_first,
            "sentaurus_electric_magnitude_V_m": norm(sent_electric),
            "sentaurus_electric_x_V_m": sent_electric[0],
            "sentaurus_electric_y_V_m": sent_electric[1],
            "sentaurus_current_x_A_m2": sent_current[0],
            "sentaurus_current_y_A_m2": sent_current[1],
            "sentaurus_electric_current_cosine": sent_cos,
            "sentaurus_p1_adjacent_cell_electric_magnitude_V_m": norm(sent_p1_electric),
            "sentaurus_p1_adjacent_cell_eparallel_V_m": sent_p1_eparallel,
            "sentaurus_p1_over_exported_vector_eparallel": (
                sent_p1_eparallel / sent_vector_first if sent_vector_first > 0.0 else math.nan
            ),
            "sentaurus_p1_nodal_recovered_electric_magnitude_V_m": norm(
                sent_p1_nodal_electric
            ),
            "sentaurus_p1_nodal_recovered_eparallel_V_m": sent_p1_nodal_eparallel,
            "sentaurus_p1_nodal_over_exported_vector_eparallel": (
                sent_p1_nodal_eparallel / sent_vector_first
                if sent_vector_first > 0.0 else math.nan
            ),
            "sentaurus_node0_p1_over_exported_electric_magnitude": (
                norm(sent_p1_nodal_e0) / norm(sent_e0) if norm(sent_e0) > 0.0 else math.nan
            ),
            "sentaurus_node1_p1_over_exported_electric_magnitude": (
                norm(sent_p1_nodal_e1) / norm(sent_e1) if norm(sent_e1) > 0.0 else math.nan
            ),
            "sentaurus_node0_alpha_inverted_eparallel_V_m": sent_inverse0,
            "sentaurus_node1_alpha_inverted_eparallel_V_m": sent_inverse1,
            "sentaurus_endpoint_mean_alpha_inverted_eparallel_V_m": 0.5 * (sent_inverse0 + sent_inverse1),
            "sentaurus_mean_alpha_inverted_eparallel_V_m": field_from_alpha(sent_alpha_mean),
            "sentaurus_endpoint_mean_alpha_m_inv": sent_alpha_mean,
            "vela_nodal_p1_over_sentaurus_endpoint_mean_alpha": (
                alpha_from_field(vela_nodal_eparallel) / sent_alpha_mean
                if sent_alpha_mean > 0.0 else math.nan
            ),
            "sentaurus_alpha_from_direct_endpoint_mean_m_inv": 0.5 * (
                alpha_from_field(sent_direct0) + alpha_from_field(sent_direct1)
            ),
            "vela_over_sentaurus_alpha_inverted_eparallel": (
                csv_eparallel / (0.5 * (sent_inverse0 + sent_inverse1))
            ),
            "vela_over_sentaurus_vector_first_eparallel": csv_eparallel / sent_vector_first,
        })
    return edge_records, cell_records, stencil_records


def growth_records(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    index = {(row["bias_V"], row["edge_id"]): row for row in records}
    output: list[dict[str, Any]] = []
    for edge_id in TARGET_EDGES:
        left, right = index[(6.4, edge_id)], index[(7.0, edge_id)]
        def growth(key: str) -> float:
            return right[key] / left[key] if left[key] != 0.0 else math.nan
        vela_e0 = (left["vela_electric_x_V_m"], left["vela_electric_y_V_m"])
        vela_e1 = (right["vela_electric_x_V_m"], right["vela_electric_y_V_m"])
        vela_j0 = (left["vela_current_x_proxy"], left["vela_current_y_proxy"])
        vela_j1 = (right["vela_current_x_proxy"], right["vela_current_y_proxy"])
        sent_e0 = (left["sentaurus_electric_x_V_m"], left["sentaurus_electric_y_V_m"])
        sent_e1 = (right["sentaurus_electric_x_V_m"], right["sentaurus_electric_y_V_m"])
        sent_j0 = (left["sentaurus_current_x_A_m2"], left["sentaurus_current_y_A_m2"])
        sent_j1 = (right["sentaurus_current_x_A_m2"], right["sentaurus_current_y_A_m2"])
        vela_base = eparallel(vela_e0, vela_j0)
        sent_base = eparallel(sent_e0, sent_j0)
        output.append({
            "edge_id": edge_id,
            "node0": left["node0"],
            "node1": left["node1"],
            "vela_eparallel_growth_factor": growth("vela_csv_eparallel_V_m"),
            "vela_electric_magnitude_growth_factor": growth("vela_electric_magnitude_V_m"),
            "vela_alignment_cosine_growth_factor": growth("vela_electric_current_cosine"),
            "vela_nodal_p1_eparallel_growth_factor": growth("vela_nodal_p1_eparallel_V_m"),
            "vela_nodal_p1_alpha_growth_factor": growth("vela_nodal_p1_alpha_m_inv"),
            "vela_field_only_growth_with_6p4_current_direction": (
                eparallel(vela_e1, vela_j0) / vela_base
            ),
            "vela_direction_only_factor_with_6p4_electric_field": (
                eparallel(vela_e0, vela_j1) / vela_base
            ),
            "vela_alpha_growth_factor": growth("vela_alpha_csv_m_inv"),
            "sentaurus_direct_eparallel_growth_factor": growth("sentaurus_endpoint_mean_direct_eparallel_V_m"),
            "sentaurus_vector_first_eparallel_growth_factor": growth("sentaurus_vector_first_eparallel_V_m"),
            "sentaurus_alpha_inverted_eparallel_growth_factor": growth("sentaurus_endpoint_mean_alpha_inverted_eparallel_V_m"),
            "sentaurus_electric_magnitude_growth_factor": growth("sentaurus_electric_magnitude_V_m"),
            "sentaurus_p1_adjacent_cell_electric_growth_factor": growth(
                "sentaurus_p1_adjacent_cell_electric_magnitude_V_m"
            ),
            "sentaurus_p1_adjacent_cell_eparallel_growth_factor": growth(
                "sentaurus_p1_adjacent_cell_eparallel_V_m"
            ),
            "sentaurus_p1_nodal_recovered_eparallel_growth_factor": growth(
                "sentaurus_p1_nodal_recovered_eparallel_V_m"
            ),
            "sentaurus_alignment_cosine_growth_factor": growth("sentaurus_electric_current_cosine"),
            "sentaurus_field_only_growth_with_6p4_current_direction": (
                eparallel(sent_e1, sent_j0) / sent_base
            ),
            "sentaurus_direction_only_factor_with_6p4_electric_field": (
                eparallel(sent_e0, sent_j1) / sent_base
            ),
            "sentaurus_alpha_growth_factor": growth("sentaurus_endpoint_mean_alpha_m_inv"),
            "vela_minus_sentaurus_alpha_inverted_eparallel_growth": (
                growth("vela_csv_eparallel_V_m")
                - growth("sentaurus_endpoint_mean_alpha_inverted_eparallel_V_m")
            ),
            "vela_minus_sentaurus_alpha_growth": (
                growth("vela_alpha_csv_m_inv")
                - growth("sentaurus_endpoint_mean_alpha_m_inv")
            ),
        })
    return output


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=Path, default=DEFAULT_MESH)
    parser.add_argument("--sentaurus-states", type=Path, default=base.DEFAULT_SENT_STATES)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()

    states = {
        round(base.state_bias(path), 9): path
        for path in args.sentaurus_states.iterdir()
        if path.is_dir() and path.name.startswith("iic_v")
    }
    mesh = mesh_context(args.mesh)
    edges: list[dict[str, Any]] = []
    cells: list[dict[str, Any]] = []
    stencils: list[dict[str, Any]] = []
    for bias, run in DEFAULT_VELA_RUNS.items():
        edge_rows, cell_rows, stencil_rows = audit_bias(
            bias, run, states[bias], mesh, TARGET_EDGES
        )
        edges += edge_rows
        cells += cell_rows
        stencils += stencil_rows
    growth = growth_records(edges)
    edge_index = {(row["bias_V"], row["edge_id"]): row for row in edges}
    implemented: list[dict[str, Any]] = []
    for bias, run in DEFAULT_IMPLEMENTED_RUNS.items():
        actual_rows = {
            int(row["edge_id"]): row
            for row in base.csv_rows(run / "sg_avalanche_edges.csv")
            if math.isclose(number(row, "bias_V"), bias, abs_tol=1.0e-9)
        }
        for edge_id in TARGET_EDGES:
            reference = edge_index[(bias, edge_id)]
            actual = actual_rows[edge_id]
            actual_field = number(actual, "electron_impact_field_V_per_m")
            actual_alpha = number(actual, "electron_alpha_m_inv")
            implemented.append({
                "bias_V": bias,
                "edge_id": edge_id,
                "actual_nodal_vertex_star_eparallel_V_m": actual_field,
                "predicted_nodal_vertex_star_eparallel_V_m": reference["vela_nodal_p1_eparallel_V_m"],
                "actual_vs_predicted_eparallel_relative_error": (
                    actual_field / reference["vela_nodal_p1_eparallel_V_m"] - 1.0
                ),
                "sentaurus_exported_vector_eparallel_V_m": reference["sentaurus_vector_first_eparallel_V_m"],
                "actual_over_sentaurus_eparallel": (
                    actual_field / reference["sentaurus_vector_first_eparallel_V_m"]
                ),
                "actual_nodal_vertex_star_alpha_m_inv": actual_alpha,
                "predicted_nodal_vertex_star_alpha_m_inv": reference["vela_nodal_p1_alpha_m_inv"],
                "actual_vs_predicted_alpha_relative_error": (
                    actual_alpha / reference["vela_nodal_p1_alpha_m_inv"] - 1.0
                ),
                "sentaurus_endpoint_mean_alpha_m_inv": reference["sentaurus_endpoint_mean_alpha_m_inv"],
                "actual_over_sentaurus_alpha": (
                    actual_alpha / reference["sentaurus_endpoint_mean_alpha_m_inv"]
                ),
                "actual_electron_flux_proxy": number(actual, "electron_flux_proxy"),
                "legacy_electron_flux_proxy": reference["vela_electron_flux_proxy"],
                "actual_over_legacy_electron_flux": (
                    number(actual, "electron_flux_proxy")
                    / reference["vela_electron_flux_proxy"]
                ),
            })
    implemented_index = {
        (row["bias_V"], row["edge_id"]): row for row in implemented
    }
    implemented_growth: list[dict[str, Any]] = []
    for edge_id in TARGET_EDGES:
        left = implemented_index[(6.4, edge_id)]
        right = implemented_index[(7.0, edge_id)]
        implemented_growth.append({
            "edge_id": edge_id,
            "actual_eparallel_growth_factor": (
                right["actual_nodal_vertex_star_eparallel_V_m"]
                / left["actual_nodal_vertex_star_eparallel_V_m"]
            ),
            "sentaurus_eparallel_growth_factor": next(
                row["sentaurus_vector_first_eparallel_growth_factor"]
                for row in growth if row["edge_id"] == edge_id
            ),
            "actual_alpha_growth_factor": (
                right["actual_nodal_vertex_star_alpha_m_inv"]
                / left["actual_nodal_vertex_star_alpha_m_inv"]
            ),
            "sentaurus_alpha_growth_factor": next(
                row["sentaurus_alpha_growth_factor"]
                for row in growth if row["edge_id"] == edge_id
            ),
        })
    max_recovery_error = max(abs(row["vela_recovery_relative_error"]) for row in edges)
    max_alpha_formula_error = max(abs(row["vela_alpha_formula_relative_error"]) for row in edges)
    sent_p1_export_errors = [
        row["sentaurus_p1_nodal_over_exported_vector_eparallel"] - 1.0
        for row in edges
    ]
    vela_nodal_sent_export_errors = [
        row["vela_nodal_p1_eparallel_V_m"]
        / row["sentaurus_vector_first_eparallel_V_m"] - 1.0
        for row in edges
    ]
    direction_only_changes = [
        row["vela_direction_only_factor_with_6p4_electric_field"] - 1.0
        for row in growth
    ]
    summary = {
        "target_edges": list(TARGET_EDGES),
        "biases_V": sorted(DEFAULT_VELA_RUNS),
        "vela_recovery_max_abs_relative_error": max_recovery_error,
        "vela_alpha_formula_max_abs_relative_error": max_alpha_formula_error,
        "sentaurus_p1_nodal_recovery_vs_exported_eparallel_relative_error_range": [
            min(sent_p1_export_errors), max(sent_p1_export_errors)
        ],
        "vela_nodal_p1_vs_sentaurus_exported_eparallel_relative_error_range": [
            min(vela_nodal_sent_export_errors), max(vela_nodal_sent_export_errors)
        ],
        "vela_fixed_field_current_direction_only_relative_change_range": [
            min(direction_only_changes), max(direction_only_changes)
        ],
        "diagnosis": (
            "The current edge-only adjacent-cell field stencil is the first differing "
            "operator. A full vertex-star area-weighted P1 gradient followed by endpoint "
            "averaging reproduces Sentaurus exported Eparallel growth; SG current direction "
            "is not the limiting term."
        ),
        "implemented_mode_max_abs_eparallel_prediction_relative_error": max(
            abs(row["actual_vs_predicted_eparallel_relative_error"])
            for row in implemented
        ),
        "implemented_mode_max_abs_alpha_prediction_relative_error": max(
            abs(row["actual_vs_predicted_alpha_relative_error"])
            for row in implemented
        ),
        "implemented_mode_electron_flux_relative_change_range": [
            min(row["actual_over_legacy_electron_flux"] - 1.0 for row in implemented),
            max(row["actual_over_legacy_electron_flux"] - 1.0 for row in implemented),
        ],
        "growth": growth,
        "implemented_mode_growth": implemented_growth,
    }
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "edge_eparallel_audit.csv", edges)
    write_csv(args.out_dir / "adjacent_cell_projection_audit.csv", cells)
    write_csv(args.out_dir / "incident_sg_stencil_audit.csv", stencils)
    write_csv(args.out_dir / "edge_growth_decomposition.csv", growth)
    write_csv(args.out_dir / "implemented_mode_validation.csv", implemented)
    write_csv(args.out_dir / "implemented_mode_growth.csv", implemented_growth)
    (args.out_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, allow_nan=True) + "\n", encoding="utf-8"
    )
    print(args.out_dir / "summary.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
