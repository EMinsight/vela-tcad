#!/usr/bin/env python3
"""Locate BVmethods NMOS electron-generation support growth differences."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

import audit_bvmethods_nmos_hotspot_slope_and_criteria as base


THRESHOLDS = (0.10, 0.30, 0.50, 0.80)
RADII_UM = (0.0025, 0.005, 0.010, 0.020, 0.040, 0.080, 0.120, 0.200)
DEFAULT_OUT = (
    base.RUN / "vela_validation/qf_vector_generation_support_growth_20260805"
)


def number(row: dict[str, str], key: str) -> float:
    return float(row.get(key, "0") or 0.0)


def mesh_data(state: Path) -> tuple[
    dict[int, tuple[float, float]], dict[int, float], list[dict[str, str]]
]:
    nodes = {
        int(row["id"]): (float(row["x_um"]), float(row["y_um"]))
        for row in base.csv_rows(state / "nodes.csv")
    }
    return nodes, base.p1_measures(state, nodes), base.csv_rows(state / "elements.csv")


def sentaurus_distribution(
    state: Path, measures: dict[int, float]
) -> dict[str, Any]:
    generation = {
        node: abs(values[0]) * 1.0e6
        for node, values in base.field(state, "eImpactIonization").items()
    }
    alpha = {
        node: abs(values[0]) * 100.0
        for node, values in base.field(state, "eAlphaAvalanche").items()
    }
    current = {
        node: math.sqrt(sum(component * component for component in values)) * 1.0e4
        for node, values in base.field(state, "eCurrentDensity").items()
    }
    source = {
        node: value * measures.get(node, 0.0)
        for node, value in generation.items()
    }
    peak_node = max(generation, key=generation.get)
    return {
        "generation": generation,
        "alpha": alpha,
        "current": current,
        "source": source,
        "peak_node": peak_node,
        "peak_generation": generation[peak_node],
        "total_source": sum(source.values()),
    }


def vela_distribution(
    edges: list[dict[str, str]], measures: dict[int, float]
) -> dict[str, Any]:
    # nodal_eparallel_p1 exports the exact carrier-specific endpoint ledger.
    # Retain the symmetric projection only for legacy diagnostic files.
    exact_nodal_ledger = bool(edges) and all(
        key in edges[0]
        for key in (
            "electron_node0_source_integral",
            "electron_node1_source_integral",
        )
    )
    source = {node: 0.0 for node in measures}
    for row in edges:
        for endpoint in (0, 1):
            node = int(row[f"node{endpoint}"])
            if node not in source:
                continue
            if exact_nodal_ledger:
                contribution = number(
                    row, f"electron_node{endpoint}_source_integral"
                ) * 1.0e-6
            else:
                contribution = 0.5 * number(
                    row, "electron_source_integral"
                ) * 1.0e-6
            source[node] += contribution
    generation = {
        node: source_value / measures[node] if measures[node] > 0.0 else 0.0
        for node, source_value in source.items()
    }
    peak_node = max(generation, key=generation.get)
    return {
        "generation": generation,
        "source": source,
        "peak_node": peak_node,
        "peak_generation": generation[peak_node],
        "total_source": sum(source.values()),
        "projection": (
            "exact_carrier_specific_endpoint_ledger"
            if exact_nodal_ledger
            else "symmetric_edge_projection"
        ),
    }


def threshold_rows(
    simulator: str,
    bias: float,
    distribution: dict[str, Any],
    measures: dict[int, float],
) -> list[dict[str, Any]]:
    total_area = sum(measures.values())
    output: list[dict[str, Any]] = []
    for threshold in THRESHOLDS:
        selected = [
            node for node, value in distribution["generation"].items()
            if value >= threshold * distribution["peak_generation"]
        ]
        area = sum(measures[node] for node in selected)
        source = sum(distribution["source"][node] for node in selected)
        output.append({
            "simulator": simulator,
            "bias_V": bias,
            "threshold_fraction_of_peak": threshold,
            "selected_nodes": len(selected),
            "support_area_m2": area,
            "support_area_um2": area * 1.0e12,
            "support_area_fraction_of_silicon": area / total_area,
            "source_per_m_s": source,
            "source_fraction": source / distribution["total_source"],
            "peak_generation_m3_s": distribution["peak_generation"],
        })
    return output


def radial_rows(
    simulator: str,
    bias: float,
    distribution: dict[str, Any],
    nodes: dict[int, tuple[float, float]],
    center: tuple[float, float],
) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    for radius in RADII_UM:
        selected = [
            node for node, point in nodes.items()
            if math.hypot(point[0] - center[0], point[1] - center[1]) <= radius
        ]
        source = sum(distribution["source"].get(node, 0.0) for node in selected)
        output.append({
            "simulator": simulator,
            "bias_V": bias,
            "radius_um": radius,
            "selected_nodes": len(selected),
            "cumulative_source_per_m_s": source,
            "cumulative_source_fraction": source / distribution["total_source"],
        })
    return output


def native_edge_threshold_rows(
    bias: float, edges: list[dict[str, str]]
) -> list[dict[str, Any]]:
    samples: list[tuple[float, float, float]] = []
    for row in edges:
        generation = (
            number(row, "electron_alpha_m_inv")
            * number(row, "electron_flux_proxy") * 1.0e4
        )
        source = number(row, "electron_source_integral") * 1.0e-6
        area = source / generation if generation > 0.0 else 0.0
        samples.append((generation, source, area))
    peak = max(sample[0] for sample in samples)
    total_source = sum(sample[1] for sample in samples)
    output: list[dict[str, Any]] = []
    for threshold in THRESHOLDS:
        selected = [sample for sample in samples if sample[0] >= threshold * peak]
        output.append({
            "bias_V": bias,
            "threshold_fraction_of_peak": threshold,
            "selected_edges": len(selected),
            "native_support_area_m2": sum(sample[2] for sample in selected),
            "native_support_area_um2": sum(sample[2] for sample in selected) * 1.0e12,
            "native_source_fraction": (
                sum(sample[1] for sample in selected) / total_source
            ),
        })
    return output


def compare_thresholds(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    indexed = {
        (row["simulator"], row["bias_V"], row["threshold_fraction_of_peak"]): row
        for row in records
    }
    output: list[dict[str, Any]] = []
    biases = sorted({row["bias_V"] for row in records})
    for bias in biases:
        for threshold in THRESHOLDS:
            sent = indexed[("sentaurus", bias, threshold)]
            vela = indexed[("vela", bias, threshold)]
            output.append({
                "bias_V": bias,
                "threshold_fraction_of_peak": threshold,
                "sentaurus_support_area_um2": sent["support_area_um2"],
                "vela_support_area_um2": vela["support_area_um2"],
                "vela_over_sentaurus_support_area": (
                    vela["support_area_m2"] / sent["support_area_m2"]
                    if sent["support_area_m2"] > 0.0 else math.nan
                ),
                "sentaurus_source_fraction": sent["source_fraction"],
                "vela_source_fraction": vela["source_fraction"],
                "vela_minus_sentaurus_source_fraction": (
                    vela["source_fraction"] - sent["source_fraction"]
                ),
            })
    return output


def compare_radial(records: list[dict[str, Any]]) -> list[dict[str, Any]]:
    indexed = {
        (row["simulator"], row["bias_V"], row["radius_um"]): row
        for row in records
    }
    output: list[dict[str, Any]] = []
    biases = sorted({row["bias_V"] for row in records})
    for bias in biases:
        for radius in RADII_UM:
            sent = indexed[("sentaurus", bias, radius)]
            vela = indexed[("vela", bias, radius)]
            output.append({
                "bias_V": bias,
                "radius_um": radius,
                "sentaurus_cumulative_source_fraction": sent["cumulative_source_fraction"],
                "vela_cumulative_source_fraction": vela["cumulative_source_fraction"],
                "vela_minus_sentaurus_cumulative_source_fraction": (
                    vela["cumulative_source_fraction"]
                    - sent["cumulative_source_fraction"]
                ),
            })
    return output


def node_growth_rows(
    nodes: dict[int, tuple[float, float]],
    sent0: dict[str, Any],
    sent1: dict[str, Any],
    vela0: dict[str, Any],
    vela1: dict[str, Any],
) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    sent_floor = max(sent0["source"].values()) * 1.0e-8
    vela_floor = max(vela0["source"].values()) * 1.0e-8
    for node, point in nodes.items():
        sent_source0 = sent0["source"].get(node, 0.0)
        sent_source1 = sent1["source"].get(node, 0.0)
        vela_source0 = vela0["source"].get(node, 0.0)
        vela_source1 = vela1["source"].get(node, 0.0)
        if sent_source0 < sent_floor or vela_source0 < vela_floor:
            continue
        sent_growth = sent_source1 / sent_source0 if sent_source0 > 0.0 else math.nan
        expected = vela_source0 * sent_growth if math.isfinite(sent_growth) else 0.0
        output.append({
            "node_id": node,
            "x_um": point[0],
            "y_um": point[1],
            "sentaurus_source_6p4_per_m_s": sent_source0,
            "sentaurus_source_7p0_per_m_s": sent_source1,
            "vela_source_6p4_per_m_s": vela_source0,
            "vela_source_7p0_per_m_s": vela_source1,
            "sentaurus_growth_factor": sent_growth,
            "vela_growth_factor": (
                vela_source1 / vela_source0 if vela_source0 > 0.0 else math.nan
            ),
            "expected_vela_7p0_if_sentaurus_growth_per_m_s": expected,
            "vela_growth_deficit_per_m_s": vela_source1 - expected,
            "vela_over_sentaurus_source_at_7p0": (
                vela_source1 / sent_source1 if sent_source1 > 0.0 else math.nan
            ),
        })
    return sorted(output, key=lambda row: row["vela_growth_deficit_per_m_s"])


def edge_growth_rows(
    edges0: list[dict[str, str]],
    edges1: list[dict[str, str]],
    sent0: dict[str, Any],
    sent1: dict[str, Any],
) -> list[dict[str, Any]]:
    by_id0 = {int(row["edge_id"]): row for row in edges0}
    by_id1 = {int(row["edge_id"]): row for row in edges1}
    total_source0 = sum(
        number(row, "electron_source_integral") * 1.0e-6 for row in edges0
    )
    sent_peak0 = max(sent0["generation"].values())
    output: list[dict[str, Any]] = []
    for edge_id in sorted(by_id0.keys() & by_id1.keys()):
        left, right = by_id0[edge_id], by_id1[edge_id]
        n0, n1 = int(left["node0"]), int(left["node1"])
        sent_g0 = 0.5 * (
            sent0["generation"].get(n0, 0.0) + sent0["generation"].get(n1, 0.0)
        )
        sent_g1 = 0.5 * (
            sent1["generation"].get(n0, 0.0) + sent1["generation"].get(n1, 0.0)
        )
        sent_growth = sent_g1 / sent_g0 if sent_g0 > 0.0 else math.nan
        sent_alpha0 = 0.5 * (
            sent0["alpha"].get(n0, 0.0) + sent0["alpha"].get(n1, 0.0)
        )
        sent_alpha1 = 0.5 * (
            sent1["alpha"].get(n0, 0.0) + sent1["alpha"].get(n1, 0.0)
        )
        sent_current0 = 0.5 * (
            sent0["current"].get(n0, 0.0) + sent0["current"].get(n1, 0.0)
        )
        sent_current1 = 0.5 * (
            sent1["current"].get(n0, 0.0) + sent1["current"].get(n1, 0.0)
        )
        source0 = number(left, "electron_source_integral") * 1.0e-6
        source1 = number(right, "electron_source_integral") * 1.0e-6
        if source0 < total_source0 * 1.0e-8 or sent_g0 < sent_peak0 * 1.0e-8:
            continue
        expected = source0 * sent_growth if math.isfinite(sent_growth) else 0.0
        output.append({
            "edge_id": edge_id,
            "node0": n0,
            "node1": n1,
            "x_mid_um": 0.5 * (number(left, "x0_um") + number(left, "x1_um")),
            "y_mid_um": 0.5 * (number(left, "y0_um") + number(left, "y1_um")),
            "edge_class": left["edge_class"],
            "sentaurus_endpoint_generation_growth_factor": sent_growth,
            "sentaurus_endpoint_alpha_growth_factor": (
                sent_alpha1 / sent_alpha0 if sent_alpha0 > 0.0 else math.nan
            ),
            "sentaurus_endpoint_current_growth_factor": (
                sent_current1 / sent_current0 if sent_current0 > 0.0 else math.nan
            ),
            "vela_source_6p4_per_m_s": source0,
            "vela_source_7p0_per_m_s": source1,
            "vela_source_growth_factor": source1 / source0 if source0 > 0.0 else math.nan,
            "vela_alpha_growth_factor": (
                number(right, "electron_alpha_m_inv")
                / number(left, "electron_alpha_m_inv")
                if number(left, "electron_alpha_m_inv") > 0.0 else math.nan
            ),
            "vela_current_growth_factor": (
                number(right, "electron_flux_proxy")
                / number(left, "electron_flux_proxy")
                if number(left, "electron_flux_proxy") > 0.0 else math.nan
            ),
            "vela_electron_impact_field_growth_factor": (
                number(right, "electron_impact_field_V_per_m")
                / number(left, "electron_impact_field_V_per_m")
                if number(left, "electron_impact_field_V_per_m") > 0.0 else math.nan
            ),
            "expected_vela_7p0_if_sentaurus_growth_per_m_s": expected,
            "vela_growth_deficit_per_m_s": source1 - expected,
            "electron_alpha_7p0_m_inv": number(right, "electron_alpha_m_inv"),
            "electron_current_density_7p0_A_m2": (
                base.Q_C * number(right, "electron_flux_proxy") * 1.0e4
            ),
        })
    return sorted(output, key=lambda row: row["vela_growth_deficit_per_m_s"])


def cell_growth_rows(
    elements: list[dict[str, str]],
    edge_rows: list[dict[str, Any]],
    nodes: dict[int, tuple[float, float]],
) -> list[dict[str, Any]]:
    edge_by_pair = {
        tuple(sorted((row["node0"], row["node1"]))): row for row in edge_rows
    }
    edge_cell_count: dict[int, int] = {}
    cell_edges: list[tuple[dict[str, str], list[dict[str, Any]]]] = []
    for cell in elements:
        if cell["material"].lower() not in {"si", "silicon"}:
            continue
        ids = [int(cell[f"node{i}"]) for i in range(3)]
        edges = [
            edge_by_pair[pair] for pair in (
                tuple(sorted((ids[0], ids[1]))),
                tuple(sorted((ids[1], ids[2]))),
                tuple(sorted((ids[2], ids[0]))),
            ) if pair in edge_by_pair
        ]
        for edge in edges:
            edge_cell_count[edge["edge_id"]] = edge_cell_count.get(edge["edge_id"], 0) + 1
        cell_edges.append((cell, edges))

    output: list[dict[str, Any]] = []
    for cell, edges in cell_edges:
        ids = [int(cell[f"node{i}"]) for i in range(3)]
        deficit = sum(
            edge["vela_growth_deficit_per_m_s"]
            / edge_cell_count[edge["edge_id"]]
            for edge in edges
        )
        output.append({
            "cell_id": int(cell["id"]),
            "node0": ids[0],
            "node1": ids[1],
            "node2": ids[2],
            "x_centroid_um": sum(nodes[node][0] for node in ids) / 3.0,
            "y_centroid_um": sum(nodes[node][1] for node in ids) / 3.0,
            "vela_growth_deficit_per_m_s": deficit,
            "edge_ids": ";".join(str(edge["edge_id"]) for edge in edges),
        })
    return sorted(output, key=lambda row: row["vela_growth_deficit_per_m_s"])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sentaurus-states", type=Path, default=base.DEFAULT_SENT_STATES)
    parser.add_argument("--vela-edges", type=Path, action="append")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--top-deficits", type=int, default=100)
    args = parser.parse_args()

    states = {
        round(base.state_bias(path), 9): path
        for path in args.sentaurus_states.iterdir()
        if path.is_dir() and path.name.startswith("iic_v")
    }
    nodes, measures, elements = mesh_data(states[6.4])
    vela_edges = base.load_vela_edges(args.vela_edges or base.DEFAULT_VELA_EDGES)
    common_biases = sorted(
        bias for bias in states.keys() & vela_edges.keys() if 6.4 <= bias <= 7.0
    )

    sent_distributions: dict[float, dict[str, Any]] = {}
    vela_distributions: dict[float, dict[str, Any]] = {}
    thresholds: list[dict[str, Any]] = []
    radial: list[dict[str, Any]] = []
    native_edge_thresholds: list[dict[str, Any]] = []
    for bias in common_biases:
        sent = sentaurus_distribution(states[bias], measures)
        vela = vela_distribution(vela_edges[bias], measures)
        raw_vela_total = sum(
            number(row, "electron_source_integral") * 1.0e-6
            for row in vela_edges[bias]
        )
        if not math.isclose(
            vela["total_source"], raw_vela_total, rel_tol=1.0e-12, abs_tol=1.0e-6
        ):
            raise ValueError(f"Vela edge-to-node source projection is not conservative at {bias:g} V")
        sent_distributions[bias] = sent
        vela_distributions[bias] = vela
        center = nodes[sent["peak_node"]]
        thresholds += threshold_rows("sentaurus", bias, sent, measures)
        thresholds += threshold_rows("vela", bias, vela, measures)
        radial += radial_rows("sentaurus", bias, sent, nodes, center)
        radial += radial_rows("vela", bias, vela, nodes, center)
        native_edge_thresholds += native_edge_threshold_rows(bias, vela_edges[bias])

    node_deficits = node_growth_rows(
        nodes,
        sent_distributions[6.4], sent_distributions[7.0],
        vela_distributions[6.4], vela_distributions[7.0],
    )
    edge_deficits = edge_growth_rows(
        vela_edges[6.4], vela_edges[7.0],
        sent_distributions[6.4], sent_distributions[7.0],
    )
    cell_deficits = cell_growth_rows(elements, edge_deficits, nodes)
    threshold_compare = compare_thresholds(thresholds)
    radial_compare = compare_radial(radial)
    threshold_index = {
        (row["bias_V"], row["threshold_fraction_of_peak"]): row
        for row in threshold_compare
    }
    native_threshold_index = {
        (row["bias_V"], row["threshold_fraction_of_peak"]): row
        for row in native_edge_thresholds
    }
    negative_edges = [
        row for row in edge_deficits if row["vela_growth_deficit_per_m_s"] < 0.0
    ]
    negative_total = sum(
        row["vela_growth_deficit_per_m_s"] for row in negative_edges
    )
    top10 = negative_edges[:10]

    summary = {
        "common_biases_V": common_biases,
        "projection": vela_distributions[common_biases[0]]["projection"],
        "integrated_electron_source_ratio_vela_over_sentaurus": {
            str(bias): (
                vela_distributions[bias]["total_source"]
                / sent_distributions[bias]["total_source"]
            )
            for bias in common_biases
        },
        "peak_electron_generation_ratio_vela_over_sentaurus": {
            str(bias): (
                vela_distributions[bias]["peak_generation"]
                / sent_distributions[bias]["peak_generation"]
            )
            for bias in common_biases
        },
        "threshold_area_ratio_6p4": {
            str(row["threshold_fraction_of_peak"]): row["vela_over_sentaurus_support_area"]
            for row in threshold_compare if row["bias_V"] == 6.4
        },
        "threshold_area_ratio_7p0": {
            str(row["threshold_fraction_of_peak"]): row["vela_over_sentaurus_support_area"]
            for row in threshold_compare if row["bias_V"] == 7.0
        },
        "threshold_support_area_growth_6p4_to_7p0": {
            str(threshold): {
                "sentaurus_factor": (
                    threshold_index[(7.0, threshold)]["sentaurus_support_area_um2"]
                    / threshold_index[(6.4, threshold)]["sentaurus_support_area_um2"]
                ),
                "vela_factor": (
                    threshold_index[(7.0, threshold)]["vela_support_area_um2"]
                    / threshold_index[(6.4, threshold)]["vela_support_area_um2"]
                ),
            }
            for threshold in THRESHOLDS
        },
        "vela_native_edge_support_area_growth_6p4_to_7p0": {
            str(threshold): (
                native_threshold_index[(7.0, threshold)]["native_support_area_um2"]
                / native_threshold_index[(6.4, threshold)]["native_support_area_um2"]
            )
            for threshold in THRESHOLDS
        },
        "top_growth_deficit_edge": edge_deficits[0],
        "top_growth_deficit_cell": cell_deficits[0],
        "top_growth_deficit_node": node_deficits[0],
        "negative_edge_growth_deficit_total_per_m_s": negative_total,
        "top10_negative_edge_deficit_fraction": (
            sum(row["vela_growth_deficit_per_m_s"] for row in top10)
            / negative_total if negative_total < 0.0 else math.nan
        ),
        "top10_negative_edge_corridor_um": {
            "x_min": min(row["x_mid_um"] for row in top10),
            "x_max": max(row["x_mid_um"] for row in top10),
            "y_min": min(row["y_mid_um"] for row in top10),
            "y_max": max(row["y_mid_um"] for row in top10),
        },
    }

    args.out_dir.mkdir(parents=True, exist_ok=True)
    base.write_csv(args.out_dir / "threshold_support_all.csv", thresholds)
    base.write_csv(args.out_dir / "threshold_support_compare.csv", threshold_compare)
    base.write_csv(args.out_dir / "radial_cumulative_source_all.csv", radial)
    base.write_csv(args.out_dir / "radial_cumulative_source_compare.csv", radial_compare)
    base.write_csv(
        args.out_dir / "vela_native_edge_threshold_support.csv",
        native_edge_thresholds,
    )
    base.write_csv(
        args.out_dir / "top_node_growth_deficits_6p4_7p0.csv",
        node_deficits[:args.top_deficits],
    )
    base.write_csv(
        args.out_dir / "top_edge_growth_deficits_6p4_7p0.csv",
        edge_deficits[:args.top_deficits],
    )
    base.write_csv(
        args.out_dir / "top_cell_growth_deficits_6p4_7p0.csv",
        cell_deficits[:args.top_deficits],
    )
    (args.out_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, allow_nan=True) + "\n", encoding="utf-8"
    )
    print(args.out_dir / "summary.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
