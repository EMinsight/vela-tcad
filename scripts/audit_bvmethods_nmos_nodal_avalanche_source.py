#!/usr/bin/env python3
"""Compare Vela mapped nodal avalanche generation with Sentaurus P1 data."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

from analyze_sentaurus_bvmethods import read_plt


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_rows(path: Path, records: list[dict[str, Any]]) -> None:
    if not records:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def percentile(values: list[float], fraction: float) -> float:
    finite = sorted(value for value in values if math.isfinite(value))
    if not finite:
        return math.nan
    position = fraction * (len(finite) - 1)
    low = int(math.floor(position))
    high = int(math.ceil(position))
    weight = position - low
    return finite[low] * (1.0 - weight) + finite[high] * weight


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vela-edges", type=Path, required=True)
    parser.add_argument("--sentaurus-state", type=Path, required=True)
    parser.add_argument("--sentaurus-plot", type=Path)
    parser.add_argument("--bias", type=float, default=6.4)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()

    nodes = {
        int(row["id"]): (float(row["x_um"]), float(row["y_um"]))
        for row in rows(args.sentaurus_state / "nodes.csv")
    }
    sent_generation = {
        int(row["node_id"]): float(row["component0"]) * 1.0e6
        for row in rows(
            args.sentaurus_state / "fields" / "ImpactIonization_region3.csv"
        )
    }
    sent_electron_generation = {
        int(row["node_id"]): float(row["component0"]) * 1.0e6
        for row in rows(
            args.sentaurus_state / "fields" / "eImpactIonization_region3.csv"
        )
    }
    sent_hole_generation = {
        int(row["node_id"]): float(row["component0"]) * 1.0e6
        for row in rows(
            args.sentaurus_state / "fields" / "hImpactIonization_region3.csv"
        )
    }
    node_measure = {node: 0.0 for node in sent_generation}
    silicon_cells = 0
    silicon_area = 0.0
    for cell in rows(args.sentaurus_state / "elements.csv"):
        if cell["material"].lower() not in {"si", "silicon"}:
            continue
        ids = [int(cell[f"node{i}"]) for i in range(3)]
        if not all(node in sent_generation for node in ids):
            continue
        points = [nodes[node] for node in ids]
        area_um2 = 0.5 * abs(
            (points[1][0] - points[0][0]) * (points[2][1] - points[0][1])
            - (points[1][1] - points[0][1]) * (points[2][0] - points[0][0])
        )
        area_m2 = area_um2 * 1.0e-12
        silicon_area += area_m2
        silicon_cells += 1
        for node in ids:
            node_measure[node] += area_m2 / 3.0

    electron_source = {node: 0.0 for node in sent_generation}
    hole_source = {node: 0.0 for node in sent_generation}
    selected_edges = 0
    vela_electron_integral = 0.0
    vela_hole_integral = 0.0
    for edge in rows(args.vela_edges):
        if not math.isclose(float(edge["bias_V"]), args.bias, abs_tol=1.0e-10):
            continue
        node0, node1 = int(edge["node0"]), int(edge["node1"])
        if node0 not in sent_generation or node1 not in sent_generation:
            continue
        # Native alpha[cm^-1]*flux[cm^-2/s]*area[um^2] converts to 1/(m*s)
        # with 1e-6. The mapped endpoint columns already contain the split.
        electron_source[node0] += float(edge.get("electron_node0_source_integral", 0.0) or 0.0)
        electron_source[node1] += float(edge.get("electron_node1_source_integral", 0.0) or 0.0)
        hole_source[node0] += float(edge.get("hole_node0_source_integral", 0.0) or 0.0)
        hole_source[node1] += float(edge.get("hole_node1_source_integral", 0.0) or 0.0)
        vela_electron_integral += float(edge["electron_source_integral"]) * 1.0e-6
        vela_hole_integral += float(edge["hole_source_integral"]) * 1.0e-6
        selected_edges += 1

    # Older dumps expose only combined node columns. Retain an exact combined
    # comparison even when component endpoint columns are unavailable.
    combined_source = {node: 0.0 for node in sent_generation}
    for edge in rows(args.vela_edges):
        if not math.isclose(float(edge["bias_V"]), args.bias, abs_tol=1.0e-10):
            continue
        node0, node1 = int(edge["node0"]), int(edge["node1"])
        if node0 not in combined_source or node1 not in combined_source:
            continue
        combined_source[node0] += float(edge["node0_source_integral"])
        combined_source[node1] += float(edge["node1_source_integral"])

    details: list[dict[str, Any]] = []
    ratios: list[float] = []
    sent_integral = 0.0
    sent_electron_integral = 0.0
    sent_hole_integral = 0.0
    vela_integral = 0.0
    for node in sorted(sent_generation):
        measure = node_measure[node]
        sent = sent_generation[node]
        vela = combined_source[node] * 1.0e-6 / measure if measure > 0.0 else 0.0
        sent_integral += sent * measure
        sent_electron_integral += sent_electron_generation[node] * measure
        sent_hole_integral += sent_hole_generation[node] * measure
        vela_integral += combined_source[node] * 1.0e-6
        ratio = vela / sent if sent > 0.0 else math.nan
        if sent > 0.0:
            ratios.append(ratio)
        x_um, y_um = nodes[node]
        details.append({
            "node_id": node,
            "x_um": x_um,
            "y_um": y_um,
            "p1_measure_m2": measure,
            "sentaurus_generation_m3_s": sent,
            "vela_generation_m3_s": vela,
            "vela_over_sentaurus": ratio,
            "sentaurus_source_per_m_s": sent * measure,
            "vela_source_per_m_s": combined_source[node] * 1.0e-6,
        })

    details.sort(key=lambda row: abs(float(row["sentaurus_source_per_m_s"])), reverse=True)
    summary = {
        "bias_V": args.bias,
        "silicon_cells": silicon_cells,
        "silicon_area_m2": silicon_area,
        "selected_edges": selected_edges,
        "sentaurus_p1_source_per_m_s": sent_integral,
        "sentaurus_p1_electron_source_per_m_s": sent_electron_integral,
        "sentaurus_p1_hole_source_per_m_s": sent_hole_integral,
        "vela_mapped_source_per_m_s": vela_integral,
        "vela_electron_source_per_m_s": vela_electron_integral,
        "vela_hole_source_per_m_s": vela_hole_integral,
        "vela_over_sentaurus_electron_integral": (
            vela_electron_integral / sent_electron_integral
        ),
        "vela_over_sentaurus_hole_integral": (
            vela_hole_integral / sent_hole_integral
        ),
        "vela_over_sentaurus_integral": vela_integral / sent_integral,
        "positive_node_ratio_p05": percentile(ratios, 0.05),
        "positive_node_ratio_p50": percentile(ratios, 0.50),
        "positive_node_ratio_p95": percentile(ratios, 0.95),
    }
    if args.sentaurus_plot is not None:
        datasets, plot_rows = read_plt(args.sentaurus_plot)
        bias_index = datasets.index("drain InnerVoltage")
        source_index = datasets.index("IntegrSemiconductor AvalancheGeneration")
        current_index = datasets.index("drain TotalCurrent")
        plot_row = min(
            plot_rows,
            key=lambda row: abs(row[bias_index] - args.bias),
        )
        plot_bias = plot_row[bias_index]
        if not math.isclose(plot_bias, args.bias, abs_tol=1.0e-8):
            raise ValueError(
                f"Sentaurus plot has no row at {args.bias:g} V; nearest is "
                f"{plot_bias:g} V"
            )
        internal_native = plot_row[source_index]
        internal_per_m_s = internal_native * 1.0e-6
        summary.update({
            "sentaurus_internal_plot_bias_V": plot_bias,
            "sentaurus_internal_source_native": internal_native,
            "sentaurus_internal_source_per_m_s": internal_per_m_s,
            "sentaurus_internal_drain_current_A_per_um": abs(plot_row[current_index]),
            "sentaurus_p1_over_internal_integral": (
                sent_integral / internal_per_m_s
            ),
            "vela_over_sentaurus_internal_integral": (
                vela_integral / internal_per_m_s
            ),
        })
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_rows(args.out_dir / "nodal_source_ledger.csv", details)
    write_rows(args.out_dir / "top100_sentaurus_source_nodes.csv", details[:100])
    (args.out_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(args.out_dir / "summary.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
