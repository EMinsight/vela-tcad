#!/usr/bin/env python3
"""Audit same-state PN2D edge-current semantics without changing solver code."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
import statistics
import sys
from pathlib import Path
from typing import Any, Iterable


ELEMENTARY_CHARGE_C = 1.602176634e-19
A_CM2_TO_PARTICLE_FLUX_M2_S = 1.0e4 / ELEMENTARY_CHARGE_C
AREA_RELATIVE_TOLERANCE = 1.0e-12
COORDINATE_MATCH_TOLERANCE_UM = 1.0e-8

WEIGHTED_TRIANGLE_FIELDS = (
    "electron_flux_proxy",
    "hole_flux_proxy",
    "electron_alpha_m_inv",
    "hole_alpha_m_inv",
    "electron_cell_qf_field_V_per_m",
    "hole_cell_qf_field_V_per_m",
    "electron_edge_qf_field_V_per_m",
    "hole_edge_qf_field_V_per_m",
    "electron_midpoint_density_m3",
    "hole_midpoint_density_m3",
    "electron_mobility_m2_V_s",
    "hole_mobility_m2_V_s",
)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--triangle-csv", type=Path, required=True)
    parser.add_argument("--sg-edge-csv", type=Path, required=True)
    parser.add_argument("--vtk-root", type=Path, required=True)
    parser.add_argument("--sentaurus-root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--biases", default="-12,-19,-19.4")
    parser.add_argument("--focus-edge", type=int, default=50)
    parser.add_argument("--top-n", type=int, default=20)
    args = parser.parse_args(argv)
    try:
        args.biases = [float(item.strip()) for item in args.biases.split(",") if item.strip()]
    except ValueError as exc:
        parser.error(f"invalid --biases value: {exc}")
    if not args.biases:
        parser.error("--biases must contain at least one bias")
    if args.top_n <= 0:
        parser.error("--top-n must be positive")
    return args


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def finite_float(value: Any, default: float = math.nan) -> float:
    if value in (None, ""):
        return default
    try:
        result = float(value)
    except (TypeError, ValueError):
        return default
    return result if math.isfinite(result) else default


def optional_float(value: Any) -> float | None:
    result = finite_float(value)
    return result if math.isfinite(result) else None


def _edge_key(row: dict[str, Any]) -> tuple[float, int]:
    return (float(row["bias_V"]), int(row["edge_id"]))


def _weighted_average(rows: list[dict[str, Any]], field: str, weights: list[float]) -> float:
    pairs = [
        (finite_float(row.get(field)), weight)
        for row, weight in zip(rows, weights)
        if math.isfinite(finite_float(row.get(field))) and weight > 0.0
    ]
    denominator = sum(weight for _, weight in pairs)
    if denominator > 0.0:
        return sum(value * weight for value, weight in pairs) / denominator
    values = [finite_float(row.get(field)) for row in rows]
    values = [value for value in values if math.isfinite(value)]
    return statistics.fmean(values) if values else math.nan


def aggregate_triangle_rows(
    rows: Iterable[dict[str, Any]],
) -> dict[tuple[float, int], dict[str, Any]]:
    grouped: dict[tuple[float, int], list[dict[str, Any]]] = {}
    for row in rows:
        grouped.setdefault(_edge_key(row), []).append(row)

    result: dict[tuple[float, int], dict[str, Any]] = {}
    for key, adjacent in grouped.items():
        weights = [max(0.0, finite_float(row.get("truncated_partial_volume_m2"), 0.0)) for row in adjacent]
        area = sum(weights)
        aggregate: dict[str, Any] = {
            "bias_V": key[0],
            "edge_id": key[1],
            "adjacent_cell_row_count": len(adjacent),
            "partial_volume_sum_m2": area,
        }
        for name in ("node0", "node1", "x0_um", "y0_um", "x1_um", "y1_um", "edge_length_m"):
            if name in adjacent[0]:
                aggregate[name] = adjacent[0][name]
        for field in WEIGHTED_TRIANGLE_FIELDS:
            aggregate[field] = _weighted_average(adjacent, field, weights)

        aggregate["electron_pdf_gradqf_flux_m2_s"] = aggregate["electron_flux_proxy"]
        aggregate["hole_pdf_gradqf_flux_m2_s"] = aggregate["hole_flux_proxy"]
        aggregate["electron_qf_source_proxy"] = sum(
            finite_float(row.get("electron_alpha_m_inv"), 0.0)
            * finite_float(row.get("electron_flux_proxy"), 0.0)
            * weight
            for row, weight in zip(adjacent, weights)
        )
        aggregate["hole_qf_source_proxy"] = sum(
            finite_float(row.get("hole_alpha_m_inv"), 0.0)
            * finite_float(row.get("hole_flux_proxy"), 0.0)
            * weight
            for row, weight in zip(adjacent, weights)
        )
        result[key] = aggregate
    return result


def enforce_shared_edge_area_gate(
    triangle_by_edge: dict[tuple[float, int], dict[str, Any]],
    sg_rows: Iterable[dict[str, Any]],
    tolerance: float = AREA_RELATIVE_TOLERANCE,
) -> dict[str, Any]:
    sg_by_edge = {_edge_key(row): row for row in sg_rows}
    common = sorted(set(triangle_by_edge).intersection(sg_by_edge))
    if not common:
        raise ValueError("shared edge area gate has no common (bias, edge_id) rows")
    errors: list[float] = []
    mismatches: list[dict[str, Any]] = []
    for key in common:
        triangle_area = finite_float(triangle_by_edge[key].get("partial_volume_sum_m2"), 0.0)
        sg_area = finite_float(sg_by_edge[key].get("edge_area_proxy_m2"), 0.0)
        scale = max(abs(triangle_area), abs(sg_area), 1.0e-300)
        relative_error = abs(triangle_area - sg_area) / scale
        errors.append(relative_error)
        if not relative_error < tolerance:
            mismatches.append({
                "bias_V": key[0], "edge_id": key[1],
                "triangle_partial_volume_sum_m2": triangle_area,
                "sg_edge_area_proxy_m2": sg_area,
                "relative_error": relative_error,
            })
    if mismatches:
        first = mismatches[0]
        raise ValueError(
            "shared edge area mismatch: "
            f"bias={first['bias_V']:g}, edge={first['edge_id']}, "
            f"relative_error={first['relative_error']:.17g}, tolerance={tolerance:.17g}"
        )
    return {
        "passed": True,
        "tolerance": tolerance,
        "common_edge_count": len(common),
        "max_relative_error": max(errors),
    }


def _canonical_endpoint_order(
    point0: tuple[float, float], point1: tuple[float, float]
) -> tuple[int, int]:
    return (0, 1) if point0 <= point1 else (1, 0)


def project_endpoint_vector(
    point0: tuple[float, float],
    point1: tuple[float, float],
    vector0: tuple[float, ...],
    vector1: tuple[float, ...],
    *,
    electron_continuity: bool = False,
) -> dict[str, float]:
    first, second = _canonical_endpoint_order(point0, point1)
    points = (point0, point1)
    start, end = points[first], points[second]
    dx, dy = end[0] - start[0], end[1] - start[1]
    length = math.hypot(dx, dy)
    if length == 0.0:
        raise ValueError("cannot project a vector on a zero-length edge")
    mean_x = 0.5 * (float(vector0[0]) + float(vector1[0]))
    mean_y = 0.5 * (float(vector0[1]) + float(vector1[1]))
    projection = mean_x * dx / length + mean_y * dy / length
    if electron_continuity:
        projection = -projection
    return {
        "canonical_projection": projection,
        "magnitude": math.hypot(mean_x, mean_y),
        "mean_x": mean_x,
        "mean_y": mean_y,
    }


def _take_numeric_tokens(lines: list[str], index: int, count: int) -> tuple[list[float], int]:
    values: list[float] = []
    while len(values) < count and index < len(lines):
        values.extend(float(token) for token in lines[index].split())
        index += 1
    if len(values) != count:
        raise ValueError(f"legacy VTK expected {count} numeric values, got {len(values)}")
    return values, index


def parse_legacy_ascii_vtk(path: Path) -> dict[str, Any]:
    lines = [line.strip() for line in path.read_text(encoding="ascii").splitlines()]
    if len(lines) < 4 or lines[2].upper() != "ASCII":
        raise ValueError(f"only legacy ASCII VTK is supported: {path}")
    result: dict[str, Any] = {"points": [], "cells": [], "point_data": {}, "cell_data": {}}
    context: str | None = None
    context_count = 0
    index = 0
    while index < len(lines):
        line = lines[index]
        parts = line.split()
        if not parts:
            index += 1
            continue
        keyword = parts[0].upper()
        if keyword == "POINTS":
            count = int(parts[1])
            values, index = _take_numeric_tokens(lines, index + 1, count * 3)
            result["points"] = [tuple(values[offset:offset + 3]) for offset in range(0, len(values), 3)]
            continue
        if keyword == "CELLS":
            count = int(parts[1])
            index += 1
            cells = []
            for _ in range(count):
                tokens = [int(token) for token in lines[index].split()]
                if not tokens or len(tokens) != tokens[0] + 1:
                    raise ValueError(f"invalid legacy VTK cell row in {path}")
                cells.append(tuple(tokens[1:]))
                index += 1
            result["cells"] = cells
            continue
        if keyword == "CELL_TYPES":
            count = int(parts[1])
            _, index = _take_numeric_tokens(lines, index + 1, count)
            continue
        if keyword in ("POINT_DATA", "CELL_DATA"):
            context = "point_data" if keyword == "POINT_DATA" else "cell_data"
            context_count = int(parts[1])
            index += 1
            continue
        if keyword == "SCALARS" and context is not None:
            name = parts[1]
            components = int(parts[3]) if len(parts) > 3 else 1
            index += 1
            if index >= len(lines) or not lines[index].upper().startswith("LOOKUP_TABLE"):
                raise ValueError(f"legacy VTK SCALARS {name} lacks LOOKUP_TABLE")
            values, index = _take_numeric_tokens(lines, index + 1, context_count * components)
            if components == 1:
                result[context][name] = values
            else:
                result[context][name] = [
                    tuple(values[offset:offset + components])
                    for offset in range(0, len(values), components)
                ]
            continue
        if keyword == "VECTORS" and context is not None:
            name = parts[1]
            values, index = _take_numeric_tokens(lines, index + 1, context_count * 3)
            result[context][name] = [
                tuple(values[offset:offset + 3]) for offset in range(0, len(values), 3)
            ]
            continue
        index += 1
    return result


_SENTAURUS_BIAS_RE = re.compile(r"^sentaurus_([+-]?[0-9]+(?:\.[0-9]+)?)v$", re.IGNORECASE)


def select_sentaurus_export(root: Path, requested_bias_V: float) -> dict[str, Any]:
    candidates: list[tuple[float, Path]] = []
    for path in root.iterdir():
        if not path.is_dir():
            continue
        match = _SENTAURUS_BIAS_RE.match(path.name)
        if match:
            candidates.append((float(match.group(1)), path))
    if not candidates:
        raise FileNotFoundError(f"no sentaurus_<bias>v exports under {root}")
    selected_bias, selected_path = min(
        candidates,
        key=lambda item: (abs(item[0] - requested_bias_V), abs(item[0]), item[0]),
    )
    return {
        "requested_bias_V": requested_bias_V,
        "selected_bias_V": selected_bias,
        "path": selected_path,
        "exact_match": math.isclose(selected_bias, requested_bias_V, abs_tol=1.0e-12),
    }


_VELA_VTK_BIAS_RE = re.compile(r"_([+-]?[0-9]+(?:\.[0-9]+)?)V\.vtk$", re.IGNORECASE)


def select_exact_vela_vtk(root: Path, bias_V: float) -> Path:
    matches = []
    for path in root.glob("dc_sweep_*_*V.vtk"):
        match = _VELA_VTK_BIAS_RE.search(path.name)
        if match and math.isclose(float(match.group(1)), bias_V, abs_tol=1.0e-12):
            matches.append(path)
    if len(matches) != 1:
        raise FileNotFoundError(
            f"expected exactly one dc_sweep_*_{bias_V:g}V.vtk under {root}, got {len(matches)}"
        )
    return matches[0]


def percentile(values: Iterable[float], fraction: float) -> float | None:
    ordered = sorted(float(value) for value in values if math.isfinite(float(value)))
    if not ordered:
        return None
    if len(ordered) == 1:
        return ordered[0]
    position = (len(ordered) - 1) * fraction
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def summarize_active_support(
    rows: list[dict[str, Any]], candidates: Iterable[str]
) -> dict[str, Any]:
    exact = [row for row in rows if bool(row.get("exact_match"))]
    positive = [
        finite_float(row.get("sentaurus_vector_magnitude_flux_m2_s"))
        for row in exact
        if finite_float(row.get("sentaurus_vector_magnitude_flux_m2_s")) > 0.0
    ]
    threshold = percentile(positive, 0.80)
    active = [
        row for row in exact
        if threshold is not None
        and finite_float(row.get("sentaurus_vector_magnitude_flux_m2_s")) > threshold
    ]
    result = {
        "exact_row_count": len(exact),
        "positive_reference_count": len(positive),
        "positive_p80_threshold": threshold,
        "active_row_count": len(active),
        "active_coverage": len(active) / len(exact) if exact else 0.0,
        "candidates": {},
    }
    for candidate in candidates:
        values = [
            abs(finite_float(row.get(f"{candidate}_log10_abs_error")))
            for row in active
            if math.isfinite(finite_float(row.get(f"{candidate}_log10_abs_error")))
        ]
        result["candidates"][candidate] = {
            "count": len(values),
            "coverage": len(values) / len(active) if active else 0.0,
            "median_abs_log10_error": statistics.median(values) if values else None,
            "p95_abs_log10_error": percentile(values, 0.95),
        }
    return result


def evaluate_contact_policy_gate(rows: list[dict[str, Any]]) -> dict[str, Any]:
    contact = [
        row for row in rows
        if bool(row.get("exact_match"))
        and row.get("edge_class") == "contact_edge"
        and bool(row.get("active_support"))
    ]
    improvements = [
        finite_float(row.get("qf_log10_abs_error"))
        - finite_float(row.get("fallback_log10_abs_error"))
        for row in contact
        if math.isfinite(finite_float(row.get("qf_log10_abs_error")))
        and math.isfinite(finite_float(row.get("fallback_log10_abs_error")))
    ]
    improved = [value for value in improvements if value >= 0.3]
    coverage = len(improved) / len(improvements) if improvements else 0.0
    interior = [
        finite_float(row.get("fallback_log10_abs_error"))
        - finite_float(row.get("qf_log10_abs_error"))
        for row in rows
        if bool(row.get("exact_match"))
        and row.get("edge_class") != "contact_edge"
        and math.isfinite(finite_float(row.get("qf_log10_abs_error")))
        and math.isfinite(finite_float(row.get("fallback_log10_abs_error")))
    ]
    interior_worsening = statistics.median(interior) if interior else None
    return {
        "exact_contact_active_count": len(contact),
        "comparable_contact_count": len(improvements),
        "improved_by_at_least_0_3_dex_count": len(improved),
        "improvement_coverage": coverage,
        "interior_count": len(interior),
        "interior_median_worsening_dex": interior_worsening,
        "recommend_explicit_contact_policy": bool(
            improvements and coverage >= 0.5
            and interior_worsening is not None and interior_worsening > 0.0
        ),
    }


def van_overstraeten_alpha(field_V_per_m: float, carrier: str) -> float:
    field = abs(float(field_V_per_m))
    if field <= 0.0 or not math.isfinite(field):
        return 0.0
    if carrier == "electron":
        prefactor, critical = 7.03e7, 1.231e8
    elif carrier == "hole" and field < 4.0e7:
        prefactor, critical = 1.582e8, 2.036e8
    elif carrier == "hole":
        prefactor, critical = 6.71e7, 1.693e8
    else:
        raise ValueError(f"unknown carrier: {carrier}")
    return prefactor * math.exp(-critical / field)


def log10_abs_error(candidate: float | None, reference: float | None) -> float | None:
    if candidate is None or reference is None:
        return None
    if not math.isfinite(candidate) or not math.isfinite(reference) or candidate == 0.0 or reference == 0.0:
        return None
    return abs(math.log10(abs(candidate) / abs(reference)))


def _nearest_point_index(
    points: list[tuple[float, ...]], target: tuple[float, float]
) -> int:
    distances = [math.hypot(point[0] - target[0], point[1] - target[1]) for point in points]
    index = min(range(len(points)), key=distances.__getitem__)
    if distances[index] > COORDINATE_MATCH_TOLERANCE_UM:
        raise ValueError(f"no coordinate match for endpoint {target}; nearest distance={distances[index]:.6g} um")
    return index


def _read_sentaurus_field(export: Path, stem: str) -> dict[int, Any]:
    files = sorted((export / "fields").glob(f"{stem}_region*.csv"))
    result: dict[int, Any] = {}
    for path in files:
        for row in read_csv(path):
            node_id = int(row["node_id"])
            components = [
                finite_float(row[name]) for name in sorted(row) if name.startswith("component")
            ]
            if len(components) == 1:
                result[node_id] = components[0]
            elif components:
                result[node_id] = tuple(components)
            elif "value" in row:
                result[node_id] = finite_float(row["value"])
    return result


def load_sentaurus_export(export: Path) -> dict[str, Any]:
    nodes = read_csv(export / "nodes.csv")
    points = {int(row["id"]): (float(row["x_um"]), float(row["y_um"])) for row in nodes}
    return {
        "points": points,
        "electron_current": _read_sentaurus_field(export, "eCurrentDensity"),
        "hole_current": _read_sentaurus_field(export, "hCurrentDensity"),
        "electron_alpha": _read_sentaurus_field(export, "eAlphaAvalanche"),
        "hole_alpha": _read_sentaurus_field(export, "hAlphaAvalanche"),
    }


def _sentaurus_node_for_coordinate(data: dict[str, Any], point: tuple[float, float]) -> int:
    node_ids = list(data["points"])
    coordinates = [data["points"][node_id] for node_id in node_ids]
    return node_ids[_nearest_point_index(coordinates, point)]


def _sentaurus_edge_semantics(
    data: dict[str, Any], point0: tuple[float, float], point1: tuple[float, float], carrier: str
) -> dict[str, float | None]:
    node0 = _sentaurus_node_for_coordinate(data, point0)
    node1 = _sentaurus_node_for_coordinate(data, point1)
    vectors = data[f"{carrier}_current"]
    if node0 not in vectors or node1 not in vectors:
        raise ValueError(f"Sentaurus {carrier} current lacks endpoint nodes {node0}/{node1}")
    projected = project_endpoint_vector(
        point0, point1, vectors[node0], vectors[node1],
        electron_continuity=(carrier == "electron"),
    )
    alpha_values = data[f"{carrier}_alpha"]
    alpha = None
    if node0 in alpha_values and node1 in alpha_values:
        # Sentaurus avalanche alpha exports are cm^-1.
        alpha = 100.0 * 0.5 * (float(alpha_values[node0]) + float(alpha_values[node1]))
    return {
        "projection_flux_m2_s": projected["canonical_projection"] * A_CM2_TO_PARTICLE_FLUX_M2_S,
        "magnitude_flux_m2_s": projected["magnitude"] * A_CM2_TO_PARTICLE_FLUX_M2_S,
        "alpha_m_inv": alpha,
    }


def _vtk_edge_semantics(
    data: dict[str, Any], point0: tuple[float, float], point1: tuple[float, float], carrier: str
) -> dict[str, float]:
    points = data["points"]
    node0 = _nearest_point_index(points, point0)
    node1 = _nearest_point_index(points, point1)
    field = "ElectronCurrentDensityVector" if carrier == "electron" else "HoleCurrentDensityVector"
    vectors = data["point_data"].get(field)
    if vectors is None:
        raise ValueError(f"VTK POINT_DATA lacks {field}")
    projected = project_endpoint_vector(
        point0, point1, vectors[node0], vectors[node1],
        electron_continuity=(carrier == "electron"),
    )
    return {
        "projection_flux_m2_s": projected["canonical_projection"] * A_CM2_TO_PARTICLE_FLUX_M2_S,
        "magnitude_flux_m2_s": projected["magnitude"] * A_CM2_TO_PARTICLE_FLUX_M2_S,
    }


def _reorder_endpoint_values(sg: dict[str, Any], first: int, second: int) -> dict[str, Any]:
    def value(base: str, endpoint: int) -> float | None:
        return optional_float(sg.get(f"{base}{endpoint}"))

    return {
        "endpoint_n0_m3": value("electron_sg_n", first),
        "endpoint_n1_m3": value("electron_sg_n", second),
        "endpoint_p0_m3": value("hole_sg_p", first),
        "endpoint_p1_m3": value("hole_sg_p", second),
        "endpoint_psi0_V": value("electron_sg_psi", first),
        "endpoint_psi1_V": value("electron_sg_psi", second),
        "endpoint_phin0_V": value("electron_sg_phin", first),
        "endpoint_phin1_V": value("electron_sg_phin", second),
        "endpoint_phip0_V": value("hole_sg_phip", first),
        "endpoint_phip1_V": value("hole_sg_phip", second),
    }


def build_audit_rows(
    triangle_rows: list[dict[str, Any]],
    sg_rows: list[dict[str, Any]],
    vtk_root: Path,
    sentaurus_root: Path,
    biases: Iterable[float],
) -> tuple[list[dict[str, Any]], dict[str, Any], list[str]]:
    requested = list(biases)
    triangle_all = aggregate_triangle_rows(triangle_rows)
    area_gate = enforce_shared_edge_area_gate(triangle_all, sg_rows)
    sg_by_edge = {_edge_key(row): row for row in sg_rows}
    rows: list[dict[str, Any]] = []
    issues: list[str] = []
    if any("hole_sg_p0" not in row or "hole_sg_phip0" not in row for row in sg_rows):
        issues.append(
            "SG edge CSV does not export hole endpoint p/phip fields; endpoint_p* and endpoint_phip* are blank."
        )

    for bias in requested:
        vtk_path = select_exact_vela_vtk(vtk_root, bias)
        vtk = parse_legacy_ascii_vtk(vtk_path)
        selection = select_sentaurus_export(sentaurus_root, bias)
        sentaurus = load_sentaurus_export(selection["path"])
        if not selection["exact_match"]:
            issues.append(
                f"Requested Sentaurus bias {bias:g} V used nearest export "
                f"{selection['selected_bias_V']:g} V; exact_match=false and excluded from accuracy summaries."
            )
        keys = sorted(key for key in set(triangle_all).intersection(sg_by_edge) if math.isclose(key[0], bias, abs_tol=1e-12))
        for key in keys:
            triangle = triangle_all[key]
            sg = sg_by_edge[key]
            point0 = (float(sg["x0_um"]), float(sg["y0_um"]))
            point1 = (float(sg["x1_um"]), float(sg["y1_um"]))
            first, second = _canonical_endpoint_order(point0, point1)
            ordered_points = (point0, point1)
            start, end = ordered_points[first], ordered_points[second]
            vela = {
                carrier: _vtk_edge_semantics(vtk, point0, point1, carrier)
                for carrier in ("electron", "hole")
            }
            sent = {
                carrier: _sentaurus_edge_semantics(sentaurus, point0, point1, carrier)
                for carrier in ("electron", "hole")
            }
            area = finite_float(sg.get("edge_area_proxy_m2"), 0.0)
            triangle_area = finite_float(triangle.get("partial_volume_sum_m2"), 0.0)
            area_error = abs(area - triangle_area) / max(abs(area), abs(triangle_area), 1.0e-300)
            electric_field = finite_float(sg.get("electric_field_V_per_m"), 0.0)
            electron_genius = optional_float(sg.get("electron_sg_production_abs_continuity_particle_flux_m2_s"))
            electron_source = "electron_sg_production_abs_continuity_particle_flux_m2_s"
            if electron_genius is None:
                electron_genius = optional_float(sg.get("electron_flux_proxy"))
                electron_source = "electron_flux_proxy_fallback"
                issues.append("SG electron production continuity flux is absent on at least one edge; electron_flux_proxy fallback used.")
            hole_genius = optional_float(sg.get("hole_raw_flux_proxy"))
            hole_source = "hole_raw_flux_proxy"
            row: dict[str, Any] = {
                "bias_V": bias,
                "sentaurus_bias_V": selection["selected_bias_V"],
                "exact_match": selection["exact_match"],
                "edge_id": key[1],
                "edge_class": sg.get("edge_class", ""),
                "contact": sg.get("edge_class") == "contact_edge",
                "source_rank": "",
                "node0": sg.get(f"node{first}"),
                "node1": sg.get(f"node{second}"),
                "x0_um": start[0], "y0_um": start[1], "x1_um": end[0], "y1_um": end[1],
                **_reorder_endpoint_values(sg, first, second),
                "electron_gss_midpoint_density_m3": triangle.get("electron_midpoint_density_m3"),
                "hole_gss_midpoint_density_m3": triangle.get("hole_midpoint_density_m3"),
                "electron_mobility_m2_V_s": triangle.get("electron_mobility_m2_V_s"),
                "hole_mobility_m2_V_s": triangle.get("hole_mobility_m2_V_s"),
                "electron_cell_qf_field_V_per_m": triangle.get("electron_cell_qf_field_V_per_m"),
                "hole_cell_qf_field_V_per_m": triangle.get("hole_cell_qf_field_V_per_m"),
                "electron_edge_qf_field_V_per_m": triangle.get("electron_edge_qf_field_V_per_m"),
                "hole_edge_qf_field_V_per_m": triangle.get("hole_edge_qf_field_V_per_m"),
                "electric_field_V_per_m": electric_field,
                "electron_qf_alpha_m_inv": triangle.get("electron_alpha_m_inv"),
                "hole_qf_alpha_m_inv": triangle.get("hole_alpha_m_inv"),
                "electron_electric_fallback_alpha_m_inv": van_overstraeten_alpha(electric_field, "electron"),
                "hole_electric_fallback_alpha_m_inv": van_overstraeten_alpha(electric_field, "hole"),
                "electron_pdf_gradqf_flux_m2_s": triangle.get("electron_pdf_gradqf_flux_m2_s"),
                "hole_pdf_gradqf_flux_m2_s": triangle.get("hole_pdf_gradqf_flux_m2_s"),
                "electron_genius_sg_flux_m2_s": electron_genius,
                "hole_genius_sg_flux_m2_s": hole_genius,
                "electron_genius_flux_source": electron_source,
                "hole_genius_flux_source": hole_source,
                "electron_vela_vector_projection_flux_m2_s": vela["electron"]["projection_flux_m2_s"],
                "electron_vela_vector_magnitude_flux_m2_s": vela["electron"]["magnitude_flux_m2_s"],
                "hole_vela_vector_projection_flux_m2_s": vela["hole"]["projection_flux_m2_s"],
                "hole_vela_vector_magnitude_flux_m2_s": vela["hole"]["magnitude_flux_m2_s"],
                "electron_sentaurus_vector_projection_flux_m2_s": sent["electron"]["projection_flux_m2_s"],
                "electron_sentaurus_vector_magnitude_flux_m2_s": sent["electron"]["magnitude_flux_m2_s"],
                "hole_sentaurus_vector_projection_flux_m2_s": sent["hole"]["projection_flux_m2_s"],
                "hole_sentaurus_vector_magnitude_flux_m2_s": sent["hole"]["magnitude_flux_m2_s"],
                "partial_volume_sum_m2": triangle_area,
                "edge_area_proxy_m2": area,
                "area_relative_error": area_error,
                "adjacent_cell_row_count": triangle.get("adjacent_cell_row_count"),
                "electron_pdf_qf_source_proxy": triangle.get("electron_qf_source_proxy"),
                "hole_pdf_qf_source_proxy": triangle.get("hole_qf_source_proxy"),
            }
            for carrier in ("electron", "hole"):
                projection_ref = sent[carrier]["projection_flux_m2_s"]
                magnitude_ref = sent[carrier]["magnitude_flux_m2_s"]
                row[f"{carrier}_pdf_log10_abs_error"] = log10_abs_error(
                    optional_float(row[f"{carrier}_pdf_gradqf_flux_m2_s"]), magnitude_ref
                )
                row[f"{carrier}_genius_log10_abs_error"] = log10_abs_error(
                    optional_float(row[f"{carrier}_genius_sg_flux_m2_s"]), magnitude_ref
                )
                row[f"{carrier}_vela_projection_log10_abs_error"] = log10_abs_error(
                    optional_float(row[f"{carrier}_vela_vector_projection_flux_m2_s"]), projection_ref
                )
                row[f"{carrier}_vela_magnitude_log10_abs_error"] = log10_abs_error(
                    optional_float(row[f"{carrier}_vela_vector_magnitude_flux_m2_s"]), magnitude_ref
                )
                qf_alpha = optional_float(row[f"{carrier}_qf_alpha_m_inv"])
                fallback_alpha = optional_float(row[f"{carrier}_electric_fallback_alpha_m_inv"])
                genius_flux = optional_float(row[f"{carrier}_genius_sg_flux_m2_s"])
                qf_source = qf_alpha * genius_flux * area if qf_alpha is not None and genius_flux is not None else None
                fallback_source = fallback_alpha * genius_flux * area if fallback_alpha is not None and genius_flux is not None else None
                sent_alpha = sent[carrier]["alpha_m_inv"]
                sent_source = sent_alpha * magnitude_ref * area if sent_alpha is not None else None
                row[f"{carrier}_qf_source_proxy"] = qf_source
                row[f"{carrier}_fallback_source_proxy"] = fallback_source
                row[f"{carrier}_sentaurus_source_proxy"] = sent_source
                row[f"{carrier}_qf_source_log10_abs_error"] = log10_abs_error(qf_source, sent_source)
                row[f"{carrier}_fallback_source_log10_abs_error"] = log10_abs_error(fallback_source, sent_source)
            ranked = []
            for name in (
                "electron_pdf", "electron_genius", "electron_vela_projection", "electron_vela_magnitude",
                "hole_pdf", "hole_genius", "hole_vela_projection", "hole_vela_magnitude",
            ):
                error = optional_float(row.get(f"{name}_log10_abs_error"))
                if error is not None:
                    ranked.append((error, name))
            ranked.sort()
            row["source_rank"] = ";".join(f"{rank + 1}:{name}" for rank, (_, name) in enumerate(ranked))
            rows.append(row)
    return rows, area_gate, sorted(set(issues))


def _candidate_summary_rows(rows: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], list[str]]:
    candidates = ["pdf", "genius", "vela_projection", "vela_magnitude"]
    long_rows = []
    for row in rows:
        for carrier in ("electron", "hole"):
            item = {
                "exact_match": row["exact_match"],
                "sentaurus_vector_magnitude_flux_m2_s": row[f"{carrier}_sentaurus_vector_magnitude_flux_m2_s"],
            }
            for candidate in candidates:
                item[f"{candidate}_log10_abs_error"] = row.get(f"{carrier}_{candidate}_log10_abs_error")
            long_rows.append(item)
    return long_rows, candidates


def _mark_active_and_contact_rows(rows: list[dict[str, Any]]) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    long_rows, candidates = _candidate_summary_rows(rows)
    summary = summarize_active_support(long_rows, candidates)
    threshold = summary["positive_p80_threshold"]
    contact_rows = []
    for row in rows:
        for carrier in ("electron", "hole"):
            magnitude = finite_float(row[f"{carrier}_sentaurus_vector_magnitude_flux_m2_s"])
            active = bool(row["exact_match"] and threshold is not None and magnitude > threshold)
            contact_rows.append({
                "exact_match": row["exact_match"],
                "edge_class": row["edge_class"],
                "active_support": active,
                "qf_log10_abs_error": row.get(f"{carrier}_qf_source_log10_abs_error"),
                "fallback_log10_abs_error": row.get(f"{carrier}_fallback_source_log10_abs_error"),
                "carrier": carrier,
                "edge_id": row["edge_id"],
                "bias_V": row["bias_V"],
            })
    return summary, contact_rows


def clean_json(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, float) and not math.isfinite(value):
        return None
    if isinstance(value, dict):
        return {key: clean_json(item) for key, item in value.items()}
    if isinstance(value, list):
        return [clean_json(item) for item in value]
    return value


def write_outputs(
    out_dir: Path,
    rows: list[dict[str, Any]],
    area_gate: dict[str, Any],
    active_summary: dict[str, Any],
    contact_gate: dict[str, Any],
    issues: list[str],
    focus_edge: int,
    top_n: int,
) -> dict[str, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / "same_state_edge_current_semantics.csv"
    json_path = out_dir / "same_state_edge_current_semantics.json"
    markdown_path = out_dir / "same_state_edge_current_semantics.md"
    if not rows:
        raise ValueError("audit produced no detail rows")
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(clean_json(rows))

    ranked = sorted(
        rows,
        key=lambda row: max(
            [finite_float(row.get(name), -1.0) for name in row if name.endswith("_log10_abs_error")],
            default=-1.0,
        ),
        reverse=True,
    )[:top_n]
    focus_rows = [row for row in rows if int(row["edge_id"]) == focus_edge]
    payload = {
        "schema": "vela.pn2d_bv_same_state_edge_current_semantics.v1",
        "row_count": len(rows),
        "area_gate": area_gate,
        "active_support": active_summary,
        "contact_gate": contact_gate,
        "focus_edge": focus_edge,
        "focus_rows": focus_rows,
        "top_n": top_n,
        "top_error_rows": ranked,
        "data_contract_issues": issues,
        "rows": rows,
    }
    json_path.write_text(json.dumps(clean_json(payload), indent=2) + "\n", encoding="utf-8")

    lines = [
        "# PN2D BV Same-State Edge Current Semantics Audit", "",
        "## Gates", "",
        f"- Shared-edge area gate: `passed` ({area_gate['common_edge_count']} common edges, "
        f"max relative error `{area_gate['max_relative_error']:.6g}`).",
        f"- Exact-bias active-support rows: `{active_summary['active_row_count']}` / "
        f"`{active_summary['exact_row_count']}`; positive p80 threshold "
        f"`{active_summary['positive_p80_threshold']}` particle m^-2 s^-1.",
        f"- Contact fallback improvement coverage: `{contact_gate['improvement_coverage']:.6g}`.",
        f"- Interior median worsening: `{contact_gate['interior_median_worsening_dex']}` dex.",
        f"- `recommend_explicit_contact_policy={str(contact_gate['recommend_explicit_contact_policy']).lower()}`.",
        "", "## Active-Support Accuracy", "",
        "| candidate | coverage | median abs log10 error | p95 abs log10 error |", "|---|---:|---:|---:|",
    ]
    for name, item in active_summary["candidates"].items():
        lines.append(
            f"| {name} | {item['coverage']:.6g} | {item['median_abs_log10_error']} | "
            f"{item['p95_abs_log10_error']} |"
        )
    lines.extend(["", "## Data Contract Issues", ""])
    if issues:
        lines.extend(f"- {issue}" for issue in issues)
    else:
        lines.append("- None.")
    lines.extend([
        "", "## Scope", "",
        "This is an offline classification audit. It does not modify production current, "
        "impact-ionization, contact, or solver policies.", "",
    ])
    markdown_path.write_text("\n".join(lines), encoding="utf-8")
    return {"csv": csv_path, "json": json_path, "markdown": markdown_path}


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    triangle_rows = read_csv(args.triangle_csv)
    sg_rows = read_csv(args.sg_edge_csv)
    rows, area_gate, issues = build_audit_rows(
        triangle_rows, sg_rows, args.vtk_root, args.sentaurus_root, args.biases
    )
    active_summary, contact_rows = _mark_active_and_contact_rows(rows)
    contact_gate = evaluate_contact_policy_gate(contact_rows)
    paths = write_outputs(
        args.out_dir, rows, area_gate, active_summary, contact_gate, issues,
        args.focus_edge, args.top_n,
    )
    print(json.dumps(clean_json({
        "rows": len(rows),
        "outputs": paths,
        "area_gate": area_gate,
        "active_support": active_summary,
        "contact_gate": contact_gate,
        "data_contract_issues": issues,
    }), indent=2))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (FileNotFoundError, ValueError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(2)
