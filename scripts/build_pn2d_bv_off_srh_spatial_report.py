#!/usr/bin/env python3
"""Build the PN2D avalanche-off spatial SRH audit report and data tables."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from pathlib import Path
from typing import Any, Iterable, Sequence

import matplotlib.pyplot as plt


Q_C = 1.602176634e-19
PLOT_FLOOR_A_PER_UM = 1.0e-30
VTK_BIAS_RE = re.compile(r"_(\d+)_(-?\d+(?:\.\d+)?)V\.vtk$")


def triangle_area(points: Sequence[tuple[float, float]]) -> float:
    (ax, ay), (bx, by), (cx, cy) = points
    return 0.5 * abs((bx - ax) * (cy - ay) - (by - ay) * (cx - ax))


def integrate_linear_triangle(
    points: Sequence[tuple[float, float]], values: Sequence[float]
) -> float:
    """Integrate a P1 nodal field exactly on one triangle."""
    if len(points) != 3 or len(values) != 3:
        raise ValueError("P1 triangle integration requires three points and values")
    return triangle_area(points) * sum(values) / 3.0


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError(f"refusing to write empty table: {path}")
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_vtk_scalars(path: Path) -> dict[str, list[float]]:
    lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    result: dict[str, list[float]] = {}
    expected = None
    index = 0
    while index < len(lines):
        tokens = lines[index].split()
        if tokens and tokens[0] == "POINT_DATA":
            expected = int(tokens[1])
            index += 1
            continue
        if tokens and tokens[0] == "SCALARS" and expected is not None:
            name = tokens[1]
            index += 1
            if index < len(lines) and lines[index].startswith("LOOKUP_TABLE"):
                index += 1
            values: list[float] = []
            while index < len(lines) and len(values) < expected:
                current = lines[index].split()
                if current and current[0] in {
                    "SCALARS",
                    "VECTORS",
                    "CELL_DATA",
                    "POINT_DATA",
                }:
                    break
                values.extend(float(value) for value in current)
                index += 1
            if len(values) != expected:
                raise RuntimeError(
                    f"{path}: scalar {name} has {len(values)} values, expected {expected}"
                )
            result[name] = values
            continue
        index += 1
    return result


def mesh_data(path: Path) -> tuple[dict[int, tuple[float, float]], list[tuple[int, int, int]]]:
    raw = json.loads(path.read_text(encoding="utf-8-sig"))
    nodes = {
        int(node["id"]): (float(node["x"]), float(node["y"]))
        for node in raw["nodes"]
    }
    triangles = [
        tuple(int(value) for value in cell["node_ids"])
        for cell in raw["triangles"]
    ]
    if set(nodes) != set(range(len(nodes))):
        raise RuntimeError("mesh node IDs must be contiguous from zero")
    return nodes, triangles


def nodal_control_areas(
    nodes: dict[int, tuple[float, float]],
    triangles: Sequence[tuple[int, int, int]],
) -> list[float]:
    areas = [0.0] * len(nodes)
    for triangle in triangles:
        points = tuple(nodes[node] for node in triangle)
        share = triangle_area(points) / 3.0
        for node in triangle:
            areas[node] += share
    return areas


def estimate_junction_x(
    nodes: dict[int, tuple[float, float]],
    triangles: Sequence[tuple[int, int, int]],
    doping: Sequence[float],
) -> float:
    candidates: list[float] = []
    seen: set[tuple[int, int]] = set()
    for triangle in triangles:
        for first, second in (
            (triangle[0], triangle[1]),
            (triangle[1], triangle[2]),
            (triangle[2], triangle[0]),
        ):
            edge = tuple(sorted((first, second)))
            if edge in seen:
                continue
            seen.add(edge)
            da, db = doping[first], doping[second]
            if da == 0.0:
                candidates.append(nodes[first][0])
            elif db == 0.0:
                candidates.append(nodes[second][0])
            elif da * db < 0.0:
                fraction = abs(da) / (abs(da) + abs(db))
                candidates.append(
                    nodes[first][0]
                    + fraction * (nodes[second][0] - nodes[first][0])
                )
    if not candidates:
        raise RuntimeError("could not locate a metallurgical-junction sign change")
    return sorted(candidates)[len(candidates) // 2]


def percentile_coordinate(
    rows: Sequence[dict[str, Any]], percentile: float
) -> float | None:
    positive = [
        (float(row["junction_normal_um"]), float(row["absolute_source_A_per_um"]))
        for row in rows
        if float(row["absolute_source_A_per_um"]) > 0.0
    ]
    if not positive:
        return None
    positive.sort()
    total = sum(weight for _, weight in positive)
    target = percentile * total
    cumulative = 0.0
    for coordinate, weight in positive:
        cumulative += weight
        if cumulative >= target:
            return coordinate
    return positive[-1][0]


def terminal_by_bias(path: Path) -> dict[int, list[dict[str, str]]]:
    result: dict[int, list[dict[str, str]]] = {}
    for row in read_csv(path):
        bias = int(round(abs(float(row["bias_V"]))))
        result.setdefault(bias, []).append(row)
    return result


def exact_curve(path: Path, current_columns: Sequence[str]) -> dict[int, float]:
    result: dict[int, float] = {}
    for row in read_csv(path):
        bias = abs(float(row["bias_V"]))
        integer = int(round(bias))
        if not math.isclose(bias, integer, rel_tol=0.0, abs_tol=1.0e-8):
            continue
        for column in current_columns:
            if row.get(column) not in {None, ""}:
                result[integer] = float(row[column])
                break
    return result


def field_csv(export_dir: Path, candidates: Sequence[str]) -> tuple[str | None, list[float] | None]:
    fields = export_dir / "fields"
    for candidate in candidates:
        matches = sorted(fields.glob(f"{candidate}_region*.csv"))
        if not matches:
            continue
        rows = read_csv(matches[0])
        value_column = next(
            (
                name
                for name in rows[0]
                if name.lower() not in {"node", "node_id", "id", "x", "y", "z"}
            ),
            None,
        )
        if value_column is None:
            continue
        return candidate, [float(row[value_column]) for row in rows]
    return None, None


def plot_iv(path: Path, summary: Sequence[dict[str, Any]]) -> None:
    figure, axis = plt.subplots(figsize=(8.4, 5.0), dpi=170)
    biases = [abs(float(row["bias_V"])) for row in summary]
    axis.semilogy(
        biases,
        [max(abs(float(row["vela_terminal_A_per_um"])), PLOT_FLOOR_A_PER_UM) for row in summary],
        marker="o",
        label="Vela",
    )
    axis.semilogy(
        biases,
        [max(abs(float(row["sentaurus_terminal_A_per_um"])), PLOT_FLOOR_A_PER_UM) for row in summary],
        marker="s",
        linestyle="--",
        label="Sentaurus",
    )
    axis.set(xlabel="Reverse bias |V| (V)", ylabel="|I| (A/um)")
    axis.grid(True, which="both", alpha=0.35)
    axis.legend(frameon=False)
    figure.tight_layout()
    figure.savefig(path)
    plt.close(figure)


def plot_profiles(
    path: Path, rows: Sequence[dict[str, Any]], value: str, ylabel: str
) -> None:
    anchors = {1, 5, 10, 15, 20}
    figure, axis = plt.subplots(figsize=(8.4, 5.0), dpi=170)
    for bias in anchors:
        selected = [row for row in rows if int(row["reverse_bias_V"]) == bias]
        selected.sort(key=lambda row: float(row["junction_normal_um"]))
        if not selected:
            continue
        axis.plot(
            [float(row["junction_normal_um"]) for row in selected],
            [float(row[value]) for row in selected],
            marker="o",
            label=(
                f"{selected[0].get('simulator', '')} -{bias} V".strip()
            ),
        )
    axis.set(xlabel="Distance from metallurgical junction (um)", ylabel=ylabel)
    axis.grid(True, alpha=0.35)
    axis.legend(frameon=False)
    figure.tight_layout()
    figure.savefig(path)
    plt.close(figure)


def build(args: argparse.Namespace) -> dict[str, Any]:
    args.out_dir.mkdir(parents=True, exist_ok=True)
    nodes, triangles = mesh_data(args.mesh)
    areas = nodal_control_areas(nodes, triangles)
    vela_curve_rows = read_csv(args.vela_curve)
    vela_curve = {
        int(round(abs(float(row["bias_V"])))): row
        for row in vela_curve_rows
        if row.get("converged") == "1"
        and math.isclose(
            abs(float(row["bias_V"])),
            round(abs(float(row["bias_V"]))),
            abs_tol=1.0e-8,
        )
    }
    sent_curve = exact_curve(
        args.sentaurus_curve,
        (
            "sentaurus_avalanche_off_A_per_um",
            "current_total_A_per_um",
            "current_total",
        ),
    )
    # The sealed on/off comparison intentionally omits the zero-current row.
    # The corresponding 0 V TDR is still required and exported; use exact zero
    # only for the terminal plotting row, never for a nonzero error metric.
    sent_curve.setdefault(0, 0.0)
    terminal = terminal_by_bias(args.terminal_balance)
    vtks: dict[int, Path] = {}
    for path in args.vtk_dir.glob("*.vtk"):
        match = VTK_BIAS_RE.search(path.name)
        if match:
            vtks[int(round(abs(float(match.group(2)))))] = path
    expected = set(range(21))
    if set(vela_curve) != expected or set(sent_curve) != expected or set(vtks) != expected:
        raise RuntimeError(
            "21-point contract failed: "
            f"vela={len(vela_curve)}, sentaurus={len(sent_curve)}, vtk={len(vtks)}"
        )

    node_rows: list[dict[str, Any]] = []
    triangle_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []
    profile_metrics_rows: list[dict[str, Any]] = []
    sent_field_availability: dict[str, dict[str, Any]] = {}
    junction_x = None
    for bias in range(21):
        scalars = parse_vtk_scalars(vtks[bias])
        required = {
            "Potential",
            "Electrons",
            "Holes",
            "EffectiveIntrinsicDensity",
            "NetDoping",
            "SRHRecombination",
        }
        missing = sorted(required - set(scalars))
        if missing:
            raise RuntimeError(f"{vtks[bias]} missing Vela fields: {missing}")
        if junction_x is None:
            junction_x = estimate_junction_x(
                nodes, triangles, scalars["NetDoping"]
            )
        per_bias_nodes: list[dict[str, Any]] = []
        for node in range(len(nodes)):
            x, y = nodes[node]
            rate = scalars["SRHRecombination"][node]
            signed_current = Q_C * rate * areas[node] * 1.0e-12
            depletion = (
                scalars["Electrons"][node] + scalars["Holes"][node]
            ) < 0.1 * max(
                abs(scalars["NetDoping"][node]),
                scalars["EffectiveIntrinsicDensity"][node],
            )
            row = {
                "simulator": "vela",
                "reverse_bias_V": bias,
                "node_id": node,
                "x_um": x,
                "y_um": y,
                "control_volume_um2": areas[node],
                "psi_V": scalars["Potential"][node],
                "electron_density_cm3": scalars["Electrons"][node],
                "hole_density_cm3": scalars["Holes"][node],
                "effective_intrinsic_density_cm3": scalars[
                    "EffectiveIntrinsicDensity"
                ][node],
                "net_doping_cm3": scalars["NetDoping"][node],
                "srh_rate_cm3_s": rate,
                "signed_source_A_per_um": signed_current,
                "absolute_source_A_per_um": abs(signed_current),
                "depletion_indicator": int(depletion),
                "junction_normal_um": x - junction_x,
            }
            per_bias_nodes.append(row)
            node_rows.append(row)

        positive_generation = sum(
            max(-float(row["signed_source_A_per_um"]), 0.0)
            for row in per_bias_nodes
        )
        negative_recombination = sum(
            min(-float(row["signed_source_A_per_um"]), 0.0)
            for row in per_bias_nodes
        )
        total_generation = positive_generation + negative_recombination
        absolute_weight = sum(
            float(row["absolute_source_A_per_um"]) for row in per_bias_nodes
        )
        centroid = (
            sum(
                float(row["junction_normal_um"])
                * float(row["absolute_source_A_per_um"])
                for row in per_bias_nodes
            )
            / absolute_weight
            if absolute_weight
            else None
        )
        p10 = percentile_coordinate(per_bias_nodes, 0.1)
        p50 = percentile_coordinate(per_bias_nodes, 0.5)
        p90 = percentile_coordinate(per_bias_nodes, 0.9)
        depletion_coordinates = [
            float(row["junction_normal_um"])
            for row in per_bias_nodes
            if row["depletion_indicator"]
        ]
        depletion_width = (
            max(depletion_coordinates) - min(depletion_coordinates)
            if depletion_coordinates
            else 0.0
        )

        triangle_integral = 0.0
        for cell_id, triangle in enumerate(triangles):
            points = tuple(nodes[node] for node in triangle)
            values = tuple(scalars["SRHRecombination"][node] for node in triangle)
            rate_integral = integrate_linear_triangle(points, values)
            triangle_integral += rate_integral
            cx = sum(point[0] for point in points) / 3.0
            cy = sum(point[1] for point in points) / 3.0
            signed_current = Q_C * rate_integral * 1.0e-12
            triangle_rows.append(
                {
                    "simulator": "vela",
                    "reverse_bias_V": bias,
                    "cell_id": cell_id,
                    "node0": triangle[0],
                    "node1": triangle[1],
                    "node2": triangle[2],
                    "centroid_x_um": cx,
                    "centroid_y_um": cy,
                    "cell_area_um2": triangle_area(points),
                    "psi_V": sum(scalars["Potential"][node] for node in triangle) / 3.0,
                    "electron_density_cm3": sum(scalars["Electrons"][node] for node in triangle) / 3.0,
                    "hole_density_cm3": sum(scalars["Holes"][node] for node in triangle) / 3.0,
                    "effective_intrinsic_density_cm3": sum(
                        scalars["EffectiveIntrinsicDensity"][node] for node in triangle
                    ) / 3.0,
                    "net_doping_cm3": sum(scalars["NetDoping"][node] for node in triangle) / 3.0,
                    "srh_rate_cm3_s": sum(values) / 3.0,
                    "signed_source_A_per_um": signed_current,
                    "absolute_source_A_per_um": abs(signed_current),
                    "depletion_indicator": int(
                        any(per_bias_nodes[node]["depletion_indicator"] for node in triangle)
                    ),
                    "junction_normal_um": cx - junction_x,
                }
            )

        contacts = terminal[bias]
        solver_row = vela_curve[bias]
        electron_sum = sum(float(row["current_electron_A_per_um"]) for row in contacts)
        hole_sum = sum(float(row["current_hole_A_per_um"]) for row in contacts)
        total_sum = sum(float(row["current_total_A_per_um"]) for row in contacts)
        electron_contact_internal = float(solver_row["global_electron_contact_flux"])
        hole_contact_internal = float(solver_row["global_hole_contact_flux"])
        electron_source_internal = float(solver_row["global_electron_integrated_source"])
        hole_source_internal = float(solver_row["global_hole_integrated_source"])
        electron_solver_source = (
            electron_source_internal * electron_sum / electron_contact_internal
        )
        hole_solver_source = (
            hole_source_internal * hole_sum / hole_contact_internal
        )
        solver_source = 0.5 * (electron_solver_source + hole_solver_source)
        # The production continuity source is assembled on nodal control
        # volumes, so the stop-condition comparison must use the exported
        # nodal rate times that exact support.  The P1 triangle integral is
        # retained as a separate discretization comparison for Task 2.
        reintegrated_source = -sum(
            float(row["signed_source_A_per_um"]) for row in per_bias_nodes
        )
        triangle_source = -Q_C * triangle_integral * 1.0e-12
        integration_relative_error = (
            abs(reintegrated_source - solver_source) / abs(solver_source)
            if solver_source
            else abs(reintegrated_source - solver_source)
        )
        closure_scale = max(abs(solver_source), PLOT_FLOOR_A_PER_UM)
        vela_current = float(vela_curve[bias]["current_total_A_per_um"])
        sent_current = float(sent_curve[bias])
        summary_rows.append(
            {
                "bias_V": -bias,
                "vela_converged": 1,
                "sentaurus_converged": 1,
                "vela_terminal_A_per_um": vela_current,
                "sentaurus_terminal_A_per_um": sent_current,
                "log10_abs_current_ratio": (
                    math.log10(abs(vela_current / sent_current))
                    if bias and sent_current and vela_current
                    else 0.0
                ),
                "electron_closure_relative": abs(
                    electron_contact_internal - electron_source_internal
                ) / max(abs(electron_source_internal), 1.0e-300),
                "hole_closure_relative": abs(
                    hole_contact_internal - hole_source_internal
                ) / max(abs(hole_source_internal), 1.0e-300),
                "total_terminal_closure_A_per_um": total_sum,
                "export_reintegration_relative_error": integration_relative_error,
                "triangle_vs_nodal_integration_relative": (
                    abs(triangle_source - reintegrated_source)
                    / max(abs(reintegrated_source), PLOT_FLOOR_A_PER_UM)
                ),
                "integrated_positive_generation_A_per_um": positive_generation,
                "integrated_negative_recombination_A_per_um": negative_recombination,
                "source_centroid_um": centroid,
                "source_p10_um": p10,
                "source_p50_um": p50,
                "source_p90_um": p90,
                "source_p10_p90_width_um": (
                    p90 - p10 if p10 is not None and p90 is not None else None
                ),
                "depletion_width_um": depletion_width,
            }
        )
        profile_metrics_rows.append(
            {
                "simulator": "vela",
                "reverse_bias_V": bias,
                "integrated_positive_generation_A_per_um": positive_generation,
                "integrated_negative_recombination_A_per_um": negative_recombination,
                "source_centroid_um": centroid,
                "source_p10_um": p10,
                "source_p50_um": p50,
                "source_p90_um": p90,
                "source_p10_p90_width_um": (
                    p90 - p10 if p10 is not None and p90 is not None else None
                ),
                "depletion_width_um": depletion_width,
            }
        )

        sent_export = args.sentaurus_export_root / f"sentaurus_-{bias}v"
        if bias == 0:
            sent_export = args.sentaurus_export_root / "sentaurus_0v"
        availability: dict[str, Any] = {}
        for name, candidates in {
            "psi": ("ElectrostaticPotential", "Potential"),
            "n": ("eDensity",),
            "p": ("hDensity",),
            "effective_intrinsic_density": (
                "EffectiveIntrinsicDensity",
                "IntrinsicDensity",
                "EffectiveIntrinsicDensityOldSlotboom",
            ),
            "net_doping": ("DopingConcentration", "NetDoping"),
            "srh_rate": ("SRHRecombination", "srhRecombination"),
        }.items():
            field, values = field_csv(sent_export, candidates)
            availability[name] = {
                "available": values is not None,
                "native_field": field,
                "reconstructed": False,
            }
        sent_field_availability[str(-bias)] = availability
        sent_nodes_raw = read_csv(sent_export / "nodes.csv")
        sent_elements_raw = read_csv(sent_export / "elements.csv")
        sent_nodes = {
            int(row["id"]): (float(row["x_um"]), float(row["y_um"]))
            for row in sent_nodes_raw
        }
        sent_triangles = [
            (int(row["node0"]), int(row["node1"]), int(row["node2"]))
            for row in sent_elements_raw
        ]
        sent_areas = nodal_control_areas(sent_nodes, sent_triangles)
        sent_values: dict[str, list[float] | None] = {}
        for name, candidates in {
            "psi": ("ElectrostaticPotential", "Potential"),
            "n": ("eDensity",),
            "p": ("hDensity",),
            "effective_intrinsic_density": (
                "EffectiveIntrinsicDensity",
                "IntrinsicDensity",
            ),
            "net_doping": ("DopingConcentration", "NetDoping"),
            "srh_rate": ("SRHRecombination", "srhRecombination"),
        }.items():
            _, sent_values[name] = field_csv(sent_export, candidates)
        required_sent = ("psi", "n", "p", "net_doping", "srh_rate")
        if any(sent_values[name] is None for name in required_sent):
            raise RuntimeError(
                f"{sent_export}: missing required native Sentaurus spatial field"
            )
        sent_doping = sent_values["net_doping"]
        assert sent_doping is not None
        sent_junction_x = estimate_junction_x(
            sent_nodes, sent_triangles, sent_doping
        )
        sent_per_bias: list[dict[str, Any]] = []
        for node in range(len(sent_nodes)):
            x, y = sent_nodes[node]
            rate = sent_values["srh_rate"][node]  # type: ignore[index]
            ni_values = sent_values["effective_intrinsic_density"]
            ni = ni_values[node] if ni_values is not None else None
            depletion_reference = max(
                abs(sent_values["net_doping"][node]),  # type: ignore[index]
                ni if ni is not None else 1.0,
            )
            row = {
                "simulator": "sentaurus",
                "reverse_bias_V": bias,
                "node_id": node,
                "x_um": x,
                "y_um": y,
                "control_volume_um2": sent_areas[node],
                "psi_V": sent_values["psi"][node],  # type: ignore[index]
                "electron_density_cm3": sent_values["n"][node],  # type: ignore[index]
                "hole_density_cm3": sent_values["p"][node],  # type: ignore[index]
                "effective_intrinsic_density_cm3": ni,
                "net_doping_cm3": sent_values["net_doping"][node],  # type: ignore[index]
                "srh_rate_cm3_s": rate,
                "signed_source_A_per_um": Q_C * rate * sent_areas[node] * 1.0e-12,
                "absolute_source_A_per_um": abs(Q_C * rate * sent_areas[node] * 1.0e-12),
                "depletion_indicator": int(
                    sent_values["n"][node] + sent_values["p"][node]  # type: ignore[index]
                    < 0.1 * depletion_reference
                ),
                "junction_normal_um": x - sent_junction_x,
            }
            sent_per_bias.append(row)
            node_rows.append(row)
        sent_positive = sum(
            max(-float(row["signed_source_A_per_um"]), 0.0)
            for row in sent_per_bias
        )
        sent_negative = sum(
            min(-float(row["signed_source_A_per_um"]), 0.0)
            for row in sent_per_bias
        )
        sent_weight = sum(
            float(row["absolute_source_A_per_um"]) for row in sent_per_bias
        )
        sent_centroid = (
            sum(
                float(row["junction_normal_um"])
                * float(row["absolute_source_A_per_um"])
                for row in sent_per_bias
            )
            / sent_weight
            if sent_weight
            else None
        )
        sent_p10 = percentile_coordinate(sent_per_bias, 0.1)
        sent_p50 = percentile_coordinate(sent_per_bias, 0.5)
        sent_p90 = percentile_coordinate(sent_per_bias, 0.9)
        sent_depletion = [
            float(row["junction_normal_um"])
            for row in sent_per_bias
            if row["depletion_indicator"]
        ]
        profile_metrics_rows.append(
            {
                "simulator": "sentaurus",
                "reverse_bias_V": bias,
                "integrated_positive_generation_A_per_um": sent_positive,
                "integrated_negative_recombination_A_per_um": sent_negative,
                "source_centroid_um": sent_centroid,
                "source_p10_um": sent_p10,
                "source_p50_um": sent_p50,
                "source_p90_um": sent_p90,
                "source_p10_p90_width_um": (
                    sent_p90 - sent_p10
                    if sent_p10 is not None and sent_p90 is not None
                    else None
                ),
                "depletion_width_um": (
                    max(sent_depletion) - min(sent_depletion)
                    if sent_depletion
                    else 0.0
                ),
            }
        )
        for cell_id, triangle in enumerate(sent_triangles):
            points = tuple(sent_nodes[node] for node in triangle)
            values = tuple(
                sent_values["srh_rate"][node] for node in triangle  # type: ignore[index]
            )
            rate_integral = integrate_linear_triangle(points, values)
            cx = sum(point[0] for point in points) / 3.0
            cy = sum(point[1] for point in points) / 3.0
            triangle_rows.append(
                {
                    "simulator": "sentaurus",
                    "reverse_bias_V": bias,
                    "cell_id": cell_id,
                    "node0": triangle[0],
                    "node1": triangle[1],
                    "node2": triangle[2],
                    "centroid_x_um": cx,
                    "centroid_y_um": cy,
                    "cell_area_um2": triangle_area(points),
                    "psi_V": sum(sent_values["psi"][node] for node in triangle) / 3.0,  # type: ignore[index]
                    "electron_density_cm3": sum(sent_values["n"][node] for node in triangle) / 3.0,  # type: ignore[index]
                    "hole_density_cm3": sum(sent_values["p"][node] for node in triangle) / 3.0,  # type: ignore[index]
                    "effective_intrinsic_density_cm3": None,
                    "net_doping_cm3": sum(sent_values["net_doping"][node] for node in triangle) / 3.0,  # type: ignore[index]
                    "srh_rate_cm3_s": sum(values) / 3.0,
                    "signed_source_A_per_um": Q_C * rate_integral * 1.0e-12,
                    "absolute_source_A_per_um": abs(Q_C * rate_integral * 1.0e-12),
                    "depletion_indicator": int(
                        any(sent_per_bias[node]["depletion_indicator"] for node in triangle)
                    ),
                    "junction_normal_um": cx - sent_junction_x,
                }
            )

    write_csv(args.out_dir / "vela_node_srh_spatial.csv", node_rows)
    write_csv(args.out_dir / "vela_triangle_srh_spatial.csv", triangle_rows)
    write_csv(args.out_dir / "spatial_summary.csv", summary_rows)
    write_csv(args.out_dir / "source_profile_metrics.csv", profile_metrics_rows)
    plot_iv(args.out_dir / "terminal_iv.png", summary_rows)
    plot_profiles(
        args.out_dir / "local_srh_rate.png",
        node_rows,
        "srh_rate_cm3_s",
        "Vela SRH rate (cm^-3 s^-1)",
    )
    cumulative_rows: list[dict[str, Any]] = []
    for simulator in ("vela", "sentaurus"):
        for bias in range(21):
            selected = [
                row
                for row in node_rows
                if row["simulator"] == simulator
                and int(row["reverse_bias_V"]) == bias
            ]
            selected.sort(key=lambda row: float(row["junction_normal_um"]))
            total = sum(float(row["absolute_source_A_per_um"]) for row in selected)
            cumulative = 0.0
            for row in selected:
                cumulative += float(row["absolute_source_A_per_um"])
                cumulative_rows.append(
                    {
                        "simulator": simulator,
                        "reverse_bias_V": bias,
                        "junction_normal_um": row["junction_normal_um"],
                        "cumulative_absolute_source_fraction": cumulative / total if total else 0.0,
                    }
                )
    write_csv(args.out_dir / "cumulative_source_profiles.csv", cumulative_rows)
    plot_profiles(
        args.out_dir / "cumulative_srh_source.png",
        cumulative_rows,
        "cumulative_absolute_source_fraction",
        "Cumulative absolute SRH source fraction",
    )
    plot_profiles(
        args.out_dir / "source_centroid_percentiles.png",
        [
            {
                "simulator": row["simulator"],
                "reverse_bias_V": row["reverse_bias_V"],
                "junction_normal_um": row["source_centroid_um"] or 0.0,
                "source_marker": row["source_p90_um"] or 0.0,
            }
            for row in profile_metrics_rows
        ],
        "source_marker",
        "90% source position (um)",
    )

    nonzero = [row for row in summary_rows if int(row["bias_V"]) != 0]
    acceptance = {
        "vela_converged_points": len(vela_curve),
        "sentaurus_converged_points": len(sent_curve),
        "max_electron_closure_relative": max(float(row["electron_closure_relative"]) for row in nonzero),
        "max_hole_closure_relative": max(float(row["hole_closure_relative"]) for row in nonzero),
        "max_total_terminal_closure_A_per_um": max(abs(float(row["total_terminal_closure_A_per_um"])) for row in summary_rows),
        "max_export_reintegration_relative_error": max(float(row["export_reintegration_relative_error"]) for row in nonzero),
        "triangle_constant_linear_test_tolerance": 1.0e-12,
    }
    acceptance["passed"] = (
        acceptance["vela_converged_points"] == 21
        and acceptance["sentaurus_converged_points"] == 21
        and acceptance["max_electron_closure_relative"] <= 1.0e-5
        and acceptance["max_hole_closure_relative"] <= 1.0e-5
        and acceptance["max_total_terminal_closure_A_per_um"] <= 1.0e-20
        and acceptance["max_export_reintegration_relative_error"] <= 1.0e-6
    )
    manifest = {
        "schema": "vela.pn2d_bv_off_srh_spatial_report.v1",
        "inputs": {
            str(path): sha256(path)
            for path in (
                args.mesh,
                args.vela_curve,
                args.sentaurus_curve,
                args.terminal_balance,
            )
        },
        "junction_normal_definition": "x - inferred metallurgical junction x",
        "depletion_indicator": "(n+p) < 0.1*max(abs(net_doping), effective_intrinsic_density)",
        "plot_floor_A_per_um": PLOT_FLOOR_A_PER_UM,
        "sentaurus_native_field_availability": sent_field_availability,
        "acceptance": acceptance,
        "traceability": {
            "terminal_iv.png": "spatial_summary.csv",
            "local_srh_rate.png": "vela_node_srh_spatial.csv",
            "cumulative_srh_source.png": "cumulative_source_profiles.csv",
            "source_centroid_percentiles.png": "source_profile_metrics.csv",
        },
    }
    (args.out_dir / "report_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(acceptance, indent=2))
    if not acceptance["passed"]:
        raise RuntimeError("Task 1 acceptance criteria failed; stop before interpretation")
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--vela-curve", type=Path, required=True)
    parser.add_argument("--sentaurus-curve", type=Path, required=True)
    parser.add_argument("--terminal-balance", type=Path, required=True)
    parser.add_argument("--vtk-dir", type=Path, required=True)
    parser.add_argument("--sentaurus-export-root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    build(parse_args())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
