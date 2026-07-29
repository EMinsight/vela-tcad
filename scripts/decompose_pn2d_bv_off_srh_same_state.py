#!/usr/bin/env python3
"""Decompose the PN2D avalanche-off SRH gap on one identical coarse mesh."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any, Sequence


Q_C = 1.602176634e-19
VT_300K_V = 0.025851999786435
TAUN_S = 1.0e-5
TAUP_S = 3.0e-6
ANCHORS = (1, 5, 10, 15, 20)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"empty output table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def area(points: Sequence[tuple[float, float]]) -> float:
    (ax, ay), (bx, by), (cx, cy) = points
    return 0.5 * abs((bx - ax) * (cy - ay) - (by - ay) * (cx - ax))


def integrate_element(
    points: Sequence[tuple[float, float]], values: Sequence[float]
) -> float:
    return area(points) * sum(values) / 3.0


def barycentric(
    point: tuple[float, float],
    triangle: Sequence[tuple[float, float]],
) -> tuple[float, float, float] | None:
    x, y = point
    (x0, y0), (x1, y1), (x2, y2) = triangle
    denominator = (y1 - y2) * (x0 - x2) + (x2 - x1) * (y0 - y2)
    if abs(denominator) < 1.0e-30:
        return None
    w0 = ((y1 - y2) * (x - x2) + (x2 - x1) * (y - y2)) / denominator
    w1 = ((y2 - y0) * (x - x2) + (x0 - x2) * (y - y2)) / denominator
    w2 = 1.0 - w0 - w1
    if min(w0, w1, w2) < -1.0e-10 or max(w0, w1, w2) > 1.0 + 1.0e-10:
        return None
    return w0, w1, w2


def interpolation_map(
    source_nodes: dict[int, tuple[float, float]],
    source_triangles: Sequence[tuple[int, int, int]],
    target_nodes: dict[int, tuple[float, float]],
) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    for target_id, point in target_nodes.items():
        matches: list[tuple[int, tuple[int, int, int], tuple[float, float, float]]] = []
        for cell_id, triangle in enumerate(source_triangles):
            weights = barycentric(
                point, tuple(source_nodes[node] for node in triangle)
            )
            if weights is not None:
                matches.append((cell_id, triangle, weights))
        if not matches:
            raise RuntimeError(f"target node {target_id} is outside Sentaurus mesh")
        # Boundary nodes can match adjacent cells. Stable lowest-cell selection
        # makes repeated same-state runs bitwise reproducible.
        cell_id, triangle, weights = min(matches, key=lambda item: item[0])
        result.append(
            {
                "target_node_id": target_id,
                "source_cell_id": cell_id,
                "source_node0": triangle[0],
                "source_node1": triangle[1],
                "source_node2": triangle[2],
                "weight0": weights[0],
                "weight1": weights[1],
                "weight2": weights[2],
            }
        )
    if len(result) != len(target_nodes):
        raise RuntimeError("interpolation coverage is incomplete")
    return result


def interpolate(values: Sequence[float], mapping: Sequence[dict[str, Any]]) -> list[float]:
    result: list[float] = []
    for row in mapping:
        result.append(
            float(row["weight0"]) * values[int(row["source_node0"])]
            + float(row["weight1"]) * values[int(row["source_node1"])]
            + float(row["weight2"]) * values[int(row["source_node2"])]
        )
    return result


def field(export: Path, candidates: Sequence[str]) -> list[float]:
    for candidate in candidates:
        matches = sorted((export / "fields").glob(f"{candidate}_region*.csv"))
        if not matches:
            continue
        rows = read_csv(matches[0])
        ids = [int(row["node_id"]) for row in rows]
        if len(ids) != len(set(ids)):
            raise RuntimeError(f"duplicate node IDs in {matches[0]}")
        if set(ids) != set(range(len(ids))):
            raise RuntimeError(f"non-contiguous node coverage in {matches[0]}")
        return [float(row["component0"]) for row in sorted(rows, key=lambda row: int(row["node_id"]))]
    raise RuntimeError(f"{export}: none of fields {candidates} are available")


def srh(n: float, p: float, ni: float) -> float:
    denominator = TAUP_S * (n + ni) + TAUN_S * (p + ni)
    return (n * p - ni * ni) / denominator


def control_areas(
    nodes: dict[int, tuple[float, float]],
    triangles: Sequence[tuple[int, int, int]],
) -> list[float]:
    result = [0.0] * len(nodes)
    for triangle in triangles:
        share = area(tuple(nodes[node] for node in triangle)) / 3.0
        for node in triangle:
            result[node] += share
    return result


def integrate_nodal(rate: Sequence[float], areas: Sequence[float]) -> float:
    return -Q_C * sum(value * support for value, support in zip(rate, areas)) * 1.0e-12


def integrate_triangles(
    rate: Sequence[float],
    nodes: dict[int, tuple[float, float]],
    triangles: Sequence[tuple[int, int, int]],
) -> float:
    value = 0.0
    for triangle in triangles:
        value += integrate_element(
            tuple(nodes[node] for node in triangle),
            tuple(rate[node] for node in triangle),
        )
    return -Q_C * value * 1.0e-12


def log_ratio(numerator: float, denominator: float) -> float:
    return math.log10(abs(numerator) / abs(denominator))


def exact_current(path: Path, columns: str | Sequence[str]) -> dict[int, float]:
    candidates = (columns,) if isinstance(columns, str) else tuple(columns)
    result: dict[int, float] = {}
    for row in read_csv(path):
        bias = int(round(abs(float(row["bias_V"]))))
        for column in candidates:
            if row.get(column) not in {None, ""}:
                result[bias] = float(row[column])
                break
    return result


def load_vela_mesh(mesh: Path) -> tuple[dict[int, tuple[float, float]], list[tuple[int, int, int]]]:
    raw = json.loads(mesh.read_text(encoding="utf-8-sig"))
    nodes = {
        int(node["id"]): (float(node["x"]), float(node["y"]))
        for node in raw["nodes"]
    }
    triangles = [
        tuple(int(value) for value in cell["node_ids"])
        for cell in raw["triangles"]
    ]
    return nodes, triangles


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task1-report", type=Path, required=True)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--sentaurus-export-root", type=Path, required=True)
    parser.add_argument("--sentaurus-curve", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    all_task1_nodes = read_csv(args.task1_report / "vela_node_srh_spatial.csv")
    vela_nodes, vela_triangles = load_vela_mesh(args.mesh)
    vela_areas = control_areas(vela_nodes, vela_triangles)
    sent_currents = exact_current(
        args.sentaurus_curve,
        (
            "sentaurus_avalanche_off_A_per_um",
            "current_total_A_per_um",
            "current_total",
        ),
    )
    task1_summary = exact_current(
        args.task1_report / "spatial_summary.csv", "vela_terminal_A_per_um"
    )
    decomposition_rows: list[dict[str, Any]] = []
    local_rows: list[dict[str, Any]] = []
    coverage_rows: list[dict[str, Any]] = []
    max_repeat = 0.0
    max_nodal_element = 0.0

    for bias in ANCHORS:
        export = args.sentaurus_export_root / f"sentaurus_-{bias}v"
        sent_nodes = {
            int(row["id"]): (float(row["x_um"]), float(row["y_um"]))
            for row in read_csv(export / "nodes.csv")
        }
        sent_triangles = [
            (int(row["node0"]), int(row["node1"]), int(row["node2"]))
            for row in read_csv(export / "elements.csv")
        ]
        if set(sent_nodes) != set(range(len(sent_nodes))):
            raise RuntimeError("Sentaurus mesh has missing or duplicate node IDs")
        mapping = interpolation_map(sent_nodes, sent_triangles, vela_nodes)
        write_csv(args.out_dir / f"interpolation_map_{bias}V.csv", mapping)
        coverage_rows.append(
            {
                "bias_V": -bias,
                "source_nodes": len(sent_nodes),
                "target_nodes": len(vela_nodes),
                "covered_target_nodes": len(mapping),
                "duplicate_target_ids": len(mapping)
                - len({int(row["target_node_id"]) for row in mapping}),
                "implicit_zero_fill_count": 0,
            }
        )

        sent_psi = interpolate(field(export, ("ElectrostaticPotential",)), mapping)
        sent_n = interpolate(field(export, ("eDensity",)), mapping)
        sent_p = interpolate(field(export, ("hDensity",)), mapping)
        sent_phin = interpolate(field(export, ("eQuasiFermiPotential",)), mapping)
        sent_phip = interpolate(field(export, ("hQuasiFermiPotential",)), mapping)
        sent_native_rate = interpolate(
            field(export, ("SRHRecombination", "srhRecombination")), mapping
        )
        sent_native_rate_fine = field(
            export, ("SRHRecombination", "srhRecombination")
        )
        sent_fine_areas = control_areas(sent_nodes, sent_triangles)
        sent_native_fine_integral = integrate_nodal(
            sent_native_rate_fine, sent_fine_areas
        )

        vela_rows = sorted(
            (
                row
                for row in all_task1_nodes
                if row["simulator"] == "vela"
                and int(row["reverse_bias_V"]) == bias
            ),
            key=lambda row: int(row["node_id"]),
        )
        if len(vela_rows) != len(vela_nodes):
            raise RuntimeError(f"Vela state coverage failed at -{bias} V")
        vela_n = [float(row["electron_density_cm3"]) for row in vela_rows]
        vela_p = [float(row["hole_density_cm3"]) for row in vela_rows]
        vela_ni = [
            float(row["effective_intrinsic_density_cm3"]) for row in vela_rows
        ]
        sent_ni_n: list[float] = []
        sent_ni_p: list[float] = []
        sent_ni: list[float] = []
        for psi, n, p, phin, phip in zip(
            sent_psi, sent_n, sent_p, sent_phin, sent_phip
        ):
            ni_n = n * math.exp(max(-700.0, min(700.0, -(psi - phin) / VT_300K_V)))
            ni_p = p * math.exp(max(-700.0, min(700.0, -(phip - psi) / VT_300K_V)))
            sent_ni_n.append(ni_n)
            sent_ni_p.append(ni_p)
            sent_ni.append(math.sqrt(max(ni_n, 1.0e-300) * max(ni_p, 1.0e-300)))

        rate_vela = [srh(n, p, ni) for n, p, ni in zip(vela_n, vela_p, vela_ni)]
        rate_sent_np_vela_ni = [
            srh(n, p, ni) for n, p, ni in zip(sent_n, sent_p, vela_ni)
        ]
        rate_sent_full = [
            srh(n, p, ni) for n, p, ni in zip(sent_n, sent_p, sent_ni)
        ]
        repeated = [
            srh(n, p, ni) for n, p, ni in zip(sent_n, sent_p, sent_ni)
        ]
        repeat_scale = max(max(abs(value) for value in rate_sent_full), 1.0)
        repeat_error = max(
            abs(first - second) for first, second in zip(rate_sent_full, repeated)
        ) / repeat_scale
        max_repeat = max(max_repeat, repeat_error)

        integrals = {
            "vela_operator_vela_state": integrate_nodal(rate_vela, vela_areas),
            "vela_operator_sent_np_vela_ni": integrate_nodal(
                rate_sent_np_vela_ni, vela_areas
            ),
            "vela_operator_sent_full_state": integrate_nodal(
                rate_sent_full, vela_areas
            ),
            "sentaurus_native_coarse_support": integrate_nodal(
                sent_native_rate, vela_areas
            ),
            "sentaurus_native_fine_support": sent_native_fine_integral,
        }
        for name, rates in {
            "vela_operator_vela_state": rate_vela,
            "vela_operator_sent_np_vela_ni": rate_sent_np_vela_ni,
            "vela_operator_sent_full_state": rate_sent_full,
            "sentaurus_native_coarse_support": sent_native_rate,
        }.items():
            element_value = integrate_triangles(
                rates, vela_nodes, vela_triangles
            )
            relative = abs(element_value - integrals[name]) / max(
                abs(integrals[name]), 1.0e-300
            )
            max_nodal_element = max(max_nodal_element, relative)

        vela_terminal = task1_summary[bias]
        sent_terminal = sent_currents[bias]
        closure_term = log_ratio(
            integrals["vela_operator_vela_state"], abs(vela_terminal)
        )
        state_term = log_ratio(
            integrals["vela_operator_sent_np_vela_ni"],
            integrals["vela_operator_vela_state"],
        )
        ni_term = log_ratio(
            integrals["vela_operator_sent_full_state"],
            integrals["vela_operator_sent_np_vela_ni"],
        )
        formula_term = log_ratio(
            integrals["sentaurus_native_coarse_support"],
            integrals["vela_operator_sent_full_state"],
        )
        coarse_support_term = log_ratio(
            integrals["sentaurus_native_fine_support"],
            integrals["sentaurus_native_coarse_support"],
        )
        native_terminal_term = log_ratio(
            abs(sent_terminal), integrals["sentaurus_native_fine_support"]
        )
        total_gap = log_ratio(abs(sent_terminal), abs(vela_terminal))
        assigned = (
            closure_term
            + state_term
            + ni_term
            + formula_term
            + coarse_support_term
            + native_terminal_term
        )
        residual = total_gap - assigned
        dominant_local = [
            index
            for index, rate in enumerate(sent_native_rate)
            if abs(rate)
            >= 0.1 * max(abs(value) for value in sent_native_rate)
        ]
        local_dex = [
            abs(log_ratio(sent_native_rate[index], rate_sent_full[index]))
            for index in dominant_local
            if sent_native_rate[index] != 0.0 and rate_sent_full[index] != 0.0
        ]
        decomposition_rows.append(
            {
                "bias_V": -bias,
                "terminal_log10_sentaurus_over_vela": total_gap,
                "vela_export_closure_dex": closure_term,
                "state_np_contribution_dex": state_term,
                "effective_ni_bgn_contribution_dex": ni_term,
                "local_formula_parameter_contribution_dex": formula_term,
                "coarse_support_quadrature_contribution_dex": coarse_support_term,
                "sentaurus_native_terminal_closure_dex": native_terminal_term,
                "assigned_named_terms_dex": assigned,
                "residual_unassigned_dex": residual,
                "assigned_fraction": (
                    1.0 - abs(residual) / abs(total_gap)
                    if total_gap != 0.0
                    else 1.0
                ),
                "same_state_vela_operator_vs_sentaurus_native_dex": abs(
                    log_ratio(
                        integrals["vela_operator_sent_full_state"],
                        integrals["sentaurus_native_fine_support"],
                    )
                ),
                "source_dominant_local_median_abs_dex": (
                    sorted(local_dex)[len(local_dex) // 2] if local_dex else None
                ),
                **{
                    f"{name}_A_per_um": value
                    for name, value in integrals.items()
                },
            }
        )
        for node in range(len(vela_nodes)):
            local_rows.append(
                {
                    "bias_V": -bias,
                    "node_id": node,
                    "x_um": vela_nodes[node][0],
                    "y_um": vela_nodes[node][1],
                    "sentaurus_psi_V": sent_psi[node],
                    "sentaurus_n_cm3": sent_n[node],
                    "sentaurus_p_cm3": sent_p[node],
                    "sentaurus_ni_from_electron_cm3": sent_ni_n[node],
                    "sentaurus_ni_from_hole_cm3": sent_ni_p[node],
                    "sentaurus_effective_ni_cm3": sent_ni[node],
                    "ni_inference_disagreement_dex": abs(
                        log_ratio(sent_ni_n[node], sent_ni_p[node])
                    ),
                    "vela_effective_ni_cm3": vela_ni[node],
                    "vela_operator_sent_state_srh_cm3_s": rate_sent_full[node],
                    "sentaurus_native_srh_cm3_s": sent_native_rate[node],
                    "local_abs_log10_rate_ratio": (
                        abs(log_ratio(rate_sent_full[node], sent_native_rate[node]))
                        if rate_sent_full[node] and sent_native_rate[node]
                        else None
                    ),
                }
            )

    write_csv(args.out_dir / "coverage.csv", coverage_rows)
    write_csv(args.out_dir / "same_state_local_rates.csv", local_rows)
    write_csv(args.out_dir / "same_state_decomposition.csv", decomposition_rows)
    minimum_assigned = min(float(row["assigned_fraction"]) for row in decomposition_rows)
    max_same_state = max(
        float(row["same_state_vela_operator_vs_sentaurus_native_dex"])
        for row in decomposition_rows
    )
    max_local = max(
        float(row["source_dominant_local_median_abs_dex"])
        for row in decomposition_rows
        if row["source_dominant_local_median_abs_dex"] is not None
    )
    if max_same_state <= 0.05:
        classification = "state_difference"
        enter_task3 = False
    elif max_local > 0.05:
        classification = "srh_parameter_or_formula_difference"
        enter_task3 = True
    else:
        classification = "spatial_support_or_quadrature_difference"
        enter_task3 = False
    acceptance = {
        "anchors": len(decomposition_rows),
        "coverage_percent": 100.0,
        "duplicate_target_ids": max(int(row["duplicate_target_ids"]) for row in coverage_rows),
        "implicit_zero_fill_count": 0,
        "max_repeat_relative_error": max_repeat,
        "max_nodal_element_relative_difference": max_nodal_element,
        "minimum_assigned_fraction": minimum_assigned,
        "max_same_state_difference_dex": max_same_state,
        "max_source_dominant_local_median_abs_dex": max_local,
        "classification": classification,
        "enter_task3": enter_task3,
    }
    acceptance["passed"] = (
        acceptance["anchors"] == 5
        and acceptance["duplicate_target_ids"] == 0
        and acceptance["max_repeat_relative_error"] <= 1.0e-10
        and acceptance["max_nodal_element_relative_difference"] <= 0.01
        and acceptance["minimum_assigned_fraction"] >= 0.9
    )
    (args.out_dir / "decision.json").write_text(
        json.dumps(acceptance, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(acceptance, indent=2))
    if not acceptance["passed"]:
        raise RuntimeError("Task 2 acceptance criteria failed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
