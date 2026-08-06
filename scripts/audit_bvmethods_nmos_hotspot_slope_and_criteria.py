#!/usr/bin/env python3
"""Audit BVmethods NMOS hotspot growth and path/global avalanche criteria."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path
from typing import Any

from analyze_sentaurus_bvmethods import read_plt


Q_C = 1.602176634e-19
REPO = Path(__file__).resolve().parents[1]
RUN = REPO / "build-release/reference_tcad/bvmethods_sentaurus2018/run01"
DEFAULT_SENT_STATES = RUN / "sentaurus_iic_multibias_exact_extended_20260803/imported"
DEFAULT_SENT_PLOT = RUN / "sentaurus_iic_multibias_exact_extended_20260803/raw/iic_multibias_des.plt"
DEFAULT_ABA_CURVE = RUN / "analysis/curves/ABA_coupled.csv"
DEFAULT_OUT = RUN / "vela_validation/qf_vector_hotspot_slope_criteria_audit_20260805"
DEFAULT_VELA_EDGES = [
    RUN / "vela_validation/btbt_e2_iic_qf_vector_fixed6p4_rerun_20260805/sg_avalanche_edges.csv",
    RUN / "vela_validation/btbt_e2_iic_qf_vector_branch_6p5_7p1_20260805/sg_avalanche_edges.csv",
    RUN / "vela_validation/btbt_e2_iic_qf_vector_fixed6p8_20260805/sg_avalanche_edges.csv",
    RUN / "vela_validation/btbt_e2_iic_qf_vector_fixed6p9_20260805/sg_avalanche_edges.csv",
    RUN / "vela_validation/btbt_e2_iic_qf_vector_fixed7p0_20260805/sg_avalanche_edges.csv",
]


def csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, records: list[dict[str, Any]]) -> None:
    if not records:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def field(path: Path, name: str) -> dict[int, tuple[float, ...]]:
    output: dict[int, tuple[float, ...]] = {}
    for row in csv_rows(path / "fields" / f"{name}_region3.csv"):
        components = tuple(
            float(row[key]) for key in sorted(row) if key.startswith("component")
        )
        output[int(row["node_id"])] = components
    return output


def scalar(data: dict[int, tuple[float, ...]], node: int) -> float:
    return data.get(node, (0.0,))[0]


def magnitude(data: dict[int, tuple[float, ...]], node: int) -> float:
    return math.sqrt(sum(value * value for value in data.get(node, (0.0,))))


def state_bias(path: Path) -> float:
    match = re.fullmatch(r"iic_v(\d+)p(\d+)", path.name)
    if not match:
        raise ValueError(f"unsupported Sentaurus state directory: {path.name}")
    return float(f"{match.group(1)}.{match.group(2)}")


def distinct_positive_plateaus(values: list[float]) -> list[float]:
    distinct: list[float] = []
    for value in sorted((v for v in values if v > 0.0 and math.isfinite(v)), reverse=True):
        if not distinct or abs(value - distinct[-1]) > max(1.0e-10, 1.0e-7 * value):
            distinct.append(value)
    return distinct


def p1_measures(state: Path, nodes: dict[int, tuple[float, float]]) -> dict[int, float]:
    measures = {node: 0.0 for node in nodes}
    for cell in csv_rows(state / "elements.csv"):
        if cell["material"].lower() not in {"si", "silicon"}:
            continue
        ids = [int(cell[f"node{i}"]) for i in range(3)]
        p0, p1, p2 = (nodes[node] for node in ids)
        area_m2 = 0.5e-12 * abs(
            (p1[0] - p0[0]) * (p2[1] - p0[1])
            - (p1[1] - p0[1]) * (p2[0] - p0[0])
        )
        for node in ids:
            measures[node] += area_m2 / 3.0
    return measures


def plot_index(path: Path) -> dict[float, dict[str, float]]:
    datasets, values = read_plt(path)
    wanted = {
        name: datasets.index(name)
        for name in (
            "drain InnerVoltage",
            "drain TotalCurrent",
            "IntegrSemiconductor AvalancheGeneration",
            "PhiElectron",
            "PhiHole",
        )
    }
    return {
        round(row[wanted["drain InnerVoltage"]], 9): {
            name: row[index] for name, index in wanted.items()
        }
        for row in values
    }


def sentaurus_record(state: Path, plot: dict[float, dict[str, float]]) -> dict[str, Any]:
    bias = state_bias(state)
    nodes = {
        int(row["id"]): (float(row["x_um"]), float(row["y_um"]))
        for row in csv_rows(state / "nodes.csv")
    }
    measures = p1_measures(state, nodes)
    eg = field(state, "eImpactIonization")
    total_g = field(state, "ImpactIonization")
    alpha = field(state, "eAlphaAvalanche")
    current = field(state, "eCurrentDensity")
    electric = field(state, "ElectricField")
    mean_integral = field(state, "MeanIonIntegral")
    peak_node = max(eg, key=lambda node: abs(scalar(eg, node)))
    plateaus = distinct_positive_plateaus([scalar(mean_integral, node) for node in mean_integral])
    plot_row = plot.get(round(bias, 9), {})
    source_native = abs(plot_row.get("IntegrSemiconductor AvalancheGeneration", math.nan))
    drain_current = abs(plot_row.get("drain TotalCurrent", math.nan))
    return {
        "bias_V": bias,
        "peak_node": peak_node,
        "peak_x_um": nodes[peak_node][0],
        "peak_y_um": nodes[peak_node][1],
        "peak_electron_generation_m3_s": abs(scalar(eg, peak_node)) * 1.0e6,
        "peak_electron_alpha_m_inv": abs(scalar(alpha, peak_node)) * 100.0,
        "peak_electron_current_density_A_m2": magnitude(current, peak_node) * 1.0e4,
        "peak_electric_field_V_m": magnitude(electric, peak_node) * 100.0,
        "peak_node_mean_ion_integral": scalar(mean_integral, peak_node),
        "electron_p1_source_per_m_s": sum(
            scalar(eg, node) * 1.0e6 * measures.get(node, 0.0) for node in eg
        ),
        "total_p1_source_per_m_s": sum(
            scalar(total_g, node) * 1.0e6 * measures.get(node, 0.0)
            for node in total_g
        ),
        "mean_integral_distinct_path_count": len(plateaus),
        "mean_integral_rank1": plateaus[0] if len(plateaus) >= 1 else 0.0,
        "mean_integral_rank2": plateaus[1] if len(plateaus) >= 2 else 0.0,
        "mean_integral_rank3": plateaus[2] if len(plateaus) >= 3 else 0.0,
        "plot_phi_electron": plot_row.get("PhiElectron", math.nan),
        "plot_phi_hole": plot_row.get("PhiHole", math.nan),
        "internal_total_source_native": source_native,
        "internal_total_source_per_m_s": source_native * 1.0e-6,
        "drain_current_A_per_um": drain_current,
        "internal_iava_A_per_um": Q_C * source_native * 1.0e-12,
        "internal_iava_minus_id_A_per_um": (
            Q_C * source_native * 1.0e-12 - drain_current
        ),
        "internal_iava_over_id": (
            Q_C * source_native * 1.0e-12 / drain_current
            if drain_current > 0.0 else math.nan
        ),
    }


def load_vela_edges(paths: list[Path]) -> dict[float, list[dict[str, str]]]:
    grouped: dict[float, list[dict[str, str]]] = {}
    # Later fixed-point files intentionally replace partial-branch dumps.
    for path in paths:
        local: dict[float, list[dict[str, str]]] = {}
        for row in csv_rows(path):
            local.setdefault(round(float(row["bias_V"]), 9), []).append(row)
        grouped.update(local)
    return grouped


def vela_record(bias: float, edges: list[dict[str, str]]) -> dict[str, Any]:
    def value(row: dict[str, str], key: str) -> float:
        return float(row.get(key, "0") or 0.0)

    def electron_generation(row: dict[str, str]) -> float:
        return value(row, "electron_alpha_m_inv") * value(row, "electron_flux_proxy") * 1.0e4

    peak = max(edges, key=electron_generation)
    return {
        "bias_V": bias,
        "peak_edge": int(peak["edge_id"]),
        "peak_x_um": 0.5 * (value(peak, "x0_um") + value(peak, "x1_um")),
        "peak_y_um": 0.5 * (value(peak, "y0_um") + value(peak, "y1_um")),
        "peak_electron_generation_m3_s": electron_generation(peak),
        "peak_electron_alpha_m_inv": value(peak, "electron_alpha_m_inv"),
        "peak_electron_current_density_A_m2": (
            Q_C * value(peak, "electron_flux_proxy") * 1.0e4
        ),
        "peak_electric_field_V_m": value(peak, "electric_field_V_per_m"),
        "peak_electron_impact_field_V_m": value(peak, "electron_impact_field_V_per_m"),
        "electron_source_per_m_s": sum(
            value(row, "electron_source_integral") for row in edges
        ) * 1.0e-6,
        "total_source_per_m_s": sum(
            value(row, "edge_source_integral") for row in edges
        ) * 1.0e-6,
    }


def comparison(sent: list[dict[str, Any]], vela: list[dict[str, Any]]) -> list[dict[str, Any]]:
    sent_by_bias = {round(row["bias_V"], 9): row for row in sent}
    output: list[dict[str, Any]] = []
    for v in vela:
        s = sent_by_bias.get(round(v["bias_V"], 9))
        if s is None:
            continue
        row: dict[str, Any] = {
            "bias_V": v["bias_V"],
            "sent_peak_node": s["peak_node"],
            "vela_peak_edge": v["peak_edge"],
            "hotspot_distance_um": math.hypot(
                v["peak_x_um"] - s["peak_x_um"],
                v["peak_y_um"] - s["peak_y_um"],
            ),
        }
        for key in (
            "peak_electron_generation_m3_s",
            "peak_electron_alpha_m_inv",
            "peak_electron_current_density_A_m2",
            "peak_electric_field_V_m",
        ):
            row[f"sentaurus_{key}"] = s[key]
            row[f"vela_{key}"] = v[key]
            row[f"vela_over_sentaurus_{key}"] = v[key] / s[key] if s[key] else math.nan
        row["sentaurus_electron_source_per_m_s"] = s["electron_p1_source_per_m_s"]
        row["vela_electron_source_per_m_s"] = v["electron_source_per_m_s"]
        row["vela_over_sentaurus_electron_source"] = (
            v["electron_source_per_m_s"] / s["electron_p1_source_per_m_s"]
        )
        row["sentaurus_effective_hotspot_area_m2"] = (
            s["electron_p1_source_per_m_s"] / s["peak_electron_generation_m3_s"]
        )
        row["vela_effective_hotspot_area_m2"] = (
            v["electron_source_per_m_s"] / v["peak_electron_generation_m3_s"]
        )
        row["vela_over_sentaurus_effective_hotspot_area"] = (
            row["vela_effective_hotspot_area_m2"]
            / row["sentaurus_effective_hotspot_area_m2"]
        )
        output.append(row)
    return output


def slopes(compare: list[dict[str, Any]]) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    keys = (
        "peak_electron_generation_m3_s",
        "peak_electron_alpha_m_inv",
        "peak_electron_current_density_A_m2",
        "electron_source_per_m_s",
    )
    for left, right in zip(compare, compare[1:]):
        delta = right["bias_V"] - left["bias_V"]
        row: dict[str, Any] = {
            "left_bias_V": left["bias_V"],
            "right_bias_V": right["bias_V"],
            "midpoint_bias_V": 0.5 * (left["bias_V"] + right["bias_V"]),
        }
        for key in keys:
            sent_key = f"sentaurus_{key}"
            vela_key = f"vela_{key}"
            sent_slope = math.log10(right[sent_key] / left[sent_key]) / delta
            vela_slope = math.log10(right[vela_key] / left[vela_key]) / delta
            row[f"sentaurus_{key}_slope_dex_V"] = sent_slope
            row[f"vela_{key}_slope_dex_V"] = vela_slope
            row[f"vela_minus_sentaurus_{key}_slope_dex_V"] = vela_slope - sent_slope
        output.append(row)
    return output


def linear_crossing(records: list[dict[str, Any]], key: str, threshold: float) -> float | None:
    for left, right in zip(records, records[1:]):
        y0, y1 = left[key], right[key]
        if y0 < threshold <= y1 and y1 != y0:
            return left["bias_V"] + (right["bias_V"] - left["bias_V"]) * (
                threshold - y0
            ) / (y1 - y0)
    return None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sentaurus-states", type=Path, default=DEFAULT_SENT_STATES)
    parser.add_argument("--sentaurus-plot", type=Path, default=DEFAULT_SENT_PLOT)
    parser.add_argument("--aba-curve", type=Path, default=DEFAULT_ABA_CURVE)
    parser.add_argument("--vela-edges", type=Path, action="append")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()

    plot = plot_index(args.sentaurus_plot)
    sent = sorted(
        (
            sentaurus_record(path, plot)
            for path in args.sentaurus_states.iterdir()
            if path.is_dir() and path.name.startswith("iic_v")
        ),
        key=lambda row: row["bias_V"],
    )
    edge_paths = args.vela_edges or DEFAULT_VELA_EDGES
    vela = [
        vela_record(bias, edges)
        for bias, edges in sorted(load_vela_edges(edge_paths).items())
        if 6.3 <= bias <= 7.0
    ]
    compared = comparison(sent, vela)
    slope_rows = slopes(compared)
    sent_dense = [row for row in sent if 6.3 <= row["bias_V"] <= 7.0]
    sparse_curve = []
    for row in csv_rows(args.aba_curve):
        current = abs(float(row["drain_total_current_A_per_um"]))
        iava = abs(float(row["avalanche_current_A_per_um"]))
        sparse_curve.append({
            "bias_V": float(row["inner_voltage_V"]),
            "iava_minus_id_A_per_um": iava - current,
        })
    summary = {
        "common_biases_V": [row["bias_V"] for row in compared],
        "sentaurus_rank3_mean_crossing_V": linear_crossing(
            sent_dense, "mean_integral_rank3", 1.0
        ),
        "sentaurus_current_source_crossing_V": linear_crossing(
            sent_dense, "internal_iava_minus_id_A_per_um", 0.0
        ),
        "sentaurus_sparse_linear_current_source_crossing_V": linear_crossing(
            sparse_curve, "iava_minus_id_A_per_um", 0.0
        ),
        "sentaurus_rank3_at_6p4": next(
            row["mean_integral_rank3"] for row in sent if row["bias_V"] == 6.4
        ),
        "sentaurus_iava_over_id_at_6p4": next(
            row["internal_iava_over_id"] for row in sent if row["bias_V"] == 6.4
        ),
        "sentaurus_rank1_at_6p4": next(
            row["mean_integral_rank1"] for row in sent if row["bias_V"] == 6.4
        ),
        "sentaurus_rank2_at_6p4": next(
            row["mean_integral_rank2"] for row in sent if row["bias_V"] == 6.4
        ),
        "sentaurus_rank1_at_7p0": next(
            row["mean_integral_rank1"] for row in sent if row["bias_V"] == 7.0
        ),
        "sentaurus_rank2_at_7p0": next(
            row["mean_integral_rank2"] for row in sent if row["bias_V"] == 7.0
        ),
        "sentaurus_rank3_at_7p0": next(
            row["mean_integral_rank3"] for row in sent if row["bias_V"] == 7.0
        ),
        "hotspot_generation_ratio_at_6p4": next(
            row["vela_over_sentaurus_peak_electron_generation_m3_s"]
            for row in compared if row["bias_V"] == 6.4
        ),
        "hotspot_generation_ratio_at_7p0": next(
            row["vela_over_sentaurus_peak_electron_generation_m3_s"]
            for row in compared if row["bias_V"] == 7.0
        ),
        "effective_hotspot_area_ratio_at_6p4": next(
            row["vela_over_sentaurus_effective_hotspot_area"]
            for row in compared if row["bias_V"] == 6.4
        ),
        "effective_hotspot_area_ratio_at_7p0": next(
            row["vela_over_sentaurus_effective_hotspot_area"]
            for row in compared if row["bias_V"] == 7.0
        ),
    }
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "sentaurus_criteria_vs_bias.csv", sent_dense)
    write_csv(args.out_dir / "vela_hotspot_vs_bias.csv", vela)
    write_csv(args.out_dir / "hotspot_same_bias_compare.csv", compared)
    write_csv(args.out_dir / "hotspot_slope_compare.csv", slope_rows)
    (args.out_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, allow_nan=True) + "\n", encoding="utf-8"
    )
    print(args.out_dir / "summary.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
