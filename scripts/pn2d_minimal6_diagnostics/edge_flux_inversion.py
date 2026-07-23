"""Deterministic edge-support and SG inversion helpers for Minimal6."""

from __future__ import annotations

import math
from collections import defaultdict
from collections.abc import Mapping, Sequence

from .qfp_sg_replacement import (
    ELEMENTARY_CHARGE_C,
    density_sg_flux,
    qf_sg_flux,
)


def _carrier(value: str) -> str:
    carrier = value.strip().lower()
    if carrier not in {"electron", "hole"}:
        raise ValueError(f"unsupported carrier {value!r}")
    return carrier


def canonical_edges(
    triangles: Sequence[Sequence[int]],
) -> tuple[tuple[int, int], ...]:
    edges: set[tuple[int, int]] = set()
    for triangle in triangles:
        if len(triangle) != 3:
            raise ValueError("every Minimal6 cell must be a triangle")
        a, b, c = (int(value) for value in triangle)
        edges.update(
            {
                tuple(sorted((a, b))),
                tuple(sorted((b, c))),
                tuple(sorted((c, a))),
            }
        )
    return tuple(sorted(edges))


def _unit_tangent(
    coordinates: Mapping[int, tuple[float, float]],
    edge: tuple[int, int],
) -> tuple[float, float]:
    x0, y0 = coordinates[edge[0]]
    x1, y1 = coordinates[edge[1]]
    length = math.hypot(x1 - x0, y1 - y0)
    if length <= 0.0:
        raise ValueError(f"edge {edge} has zero length")
    return (x1 - x0) / length, (y1 - y0) / length


def _project(vector: tuple[float, float], tangent: tuple[float, float]) -> float:
    return vector[0] * tangent[0] + vector[1] * tangent[1]


def edge_current_supports(
    coordinates: Mapping[int, tuple[float, float]],
    triangles: Sequence[Sequence[int]],
    node_vectors: Mapping[int, tuple[float, float]],
) -> dict[tuple[int, int], dict[str, float]]:
    """Reconstruct declared scalar edge supports from nodal vector current.

    The P1 line integral is listed explicitly to prove that it is exactly the
    endpoint mean, not an independent edge-current observation.
    """

    triangle_vectors: list[tuple[float, float]] = []
    adjacent: dict[tuple[int, int], list[int]] = defaultdict(list)
    for triangle_index, raw_triangle in enumerate(triangles):
        triangle = tuple(int(value) for value in raw_triangle)
        if len(triangle) != 3:
            raise ValueError("every Minimal6 cell must be a triangle")
        missing = [node for node in triangle if node not in node_vectors]
        if missing:
            raise ValueError(f"node current is missing nodes {missing}")
        triangle_vectors.append(
            (
                sum(node_vectors[node][0] for node in triangle) / 3.0,
                sum(node_vectors[node][1] for node in triangle) / 3.0,
            )
        )
        for edge in canonical_edges((triangle,)):
            adjacent[edge].append(triangle_index)

    result: dict[tuple[int, int], dict[str, float]] = {}
    for edge in canonical_edges(triangles):
        tangent = _unit_tangent(coordinates, edge)
        endpoint_vector = (
            (node_vectors[edge[0]][0] + node_vectors[edge[1]][0]) * 0.5,
            (node_vectors[edge[0]][1] + node_vectors[edge[1]][1]) * 0.5,
        )
        cell_ids = adjacent[edge]
        cell_vector = (
            sum(triangle_vectors[index][0] for index in cell_ids) / len(cell_ids),
            sum(triangle_vectors[index][1] for index in cell_ids) / len(cell_ids),
        )
        endpoint_tangent = _project(endpoint_vector, tangent)
        result[edge] = {
            "endpoint_mean_tangent": endpoint_tangent,
            "p1_line_mean_tangent": endpoint_tangent,
            "adjacent_cell_mean_tangent": _project(cell_vector, tangent),
            "endpoint_mean_magnitude": math.hypot(*endpoint_vector),
            "adjacent_cell_mean_magnitude": math.hypot(*cell_vector),
        }
    return result


def continuity_flux_from_current(carrier: str, current_A_per_m2: float) -> float:
    sign = -1.0 if _carrier(carrier) == "electron" else 1.0
    return sign * float(current_A_per_m2) / ELEMENTARY_CHARGE_C


def required_positive_mobility(
    *, reference_flux: float, unit_mobility_flux: float
) -> dict[str, str | float | None]:
    reference = float(reference_flux)
    unit = float(unit_mobility_flux)
    if not math.isfinite(reference) or not math.isfinite(unit):
        return {"classification": "nonfinite", "mobility_m2_per_Vs": None}
    if unit == 0.0:
        return {"classification": "zero_operator", "mobility_m2_per_Vs": None}
    mobility = reference / unit
    if mobility < 0.0:
        return {
            "classification": "sign_incompatible",
            "mobility_m2_per_Vs": None,
        }
    return {
        "classification": "available",
        "mobility_m2_per_Vs": mobility,
    }


def staged_sg_flux(
    formulation: str,
    carrier: str,
    state: Mapping[str, float],
) -> float:
    length = float(state["length_m"])
    if length <= 0.0:
        raise ValueError("edge length must be positive")
    thermal_voltage = float(state["thermal_voltage_V"])
    coefficient = (
        float(state["mobility_m2_per_Vs"]) * thermal_voltage / length
    )
    if formulation == "qf_sg":
        return qf_sg_flux(
            carrier,
            state["ni0_m3"],
            state["ni1_m3"],
            state["psi0_V"],
            state["psi1_V"],
            state["qf0_V"],
            state["qf1_V"],
            thermal_voltage,
            coefficient,
        )
    if formulation == "density_sg":
        return density_sg_flux(
            carrier,
            state["density0_m3"],
            state["density1_m3"],
            state["psi0_V"],
            state["psi1_V"],
            thermal_voltage,
            coefficient,
        )
    raise ValueError(f"unsupported SG formulation {formulation!r}")
