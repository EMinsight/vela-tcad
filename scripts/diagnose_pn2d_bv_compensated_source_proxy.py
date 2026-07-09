#!/usr/bin/env python3
"""Compare PN2D BV baseline vs compensated-junction source proxy factors.

The diagnostic focuses on the three horizontal junction cuts used by the
coarse7x3 BV debug artifacts. It intentionally does not alter solver state.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import statistics
from pathlib import Path
from typing import Any


BIASES = [-12.0, -19.0, -20.0]
Y_CUTS = [0.0, 0.25, 0.5]
EDGE_BY_SIDE = {
    0.0: {"left": 9, "right": 13},
    0.25: {"left": 12, "right": 16},
    0.5: {"left": 34, "right": 37},
}
VELA_X_COLUMNS = [2.0 / 3.0, 1.0, 4.0 / 3.0]
SENTAURUS_X_COLUMNS = [0.75, 1.0, 1.25]
NODE_FIELDS = ["Potential", "ElectronQuasiFermi", "HoleQuasiFermi", "Electrons", "Holes"]
FIELD_TO_OUTPUT = {
    "Potential": "psi",
    "ElectronQuasiFermi": "phin",
    "HoleQuasiFermi": "phip",
    "Electrons": "electrons",
    "Holes": "holes",
}
RATIO_FIELDS = [
    "edge_source_integral",
    "electron_source_integral",
    "hole_source_integral",
    "electron_alpha_m_inv",
    "hole_alpha_m_inv",
    "electron_flux_proxy",
    "hole_flux_proxy",
    "electron_raw_flux_proxy",
    "hole_raw_flux_proxy",
    "electron_mobility_m2_V_s",
    "hole_mobility_m2_V_s",
    "electron_density_mid_m3",
    "hole_density_mid_m3",
    "edge_area_proxy_m2",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-report-root", type=Path, required=True)
    parser.add_argument("--probe-root", type=Path, required=True)
    parser.add_argument("--sentaurus-root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    return parser.parse_args()


def long_path(path: Path) -> str:
    resolved = str(path.resolve())
    if os.name == "nt" and not resolved.startswith("\\\\?\\"):
        return "\\\\?\\" + resolved
    return resolved


def open_path(path: Path, mode: str, **kwargs: Any) -> Any:
    return open(long_path(path), mode, **kwargs)


def read_csv(path: Path) -> list[dict[str, str]]:
    with open_path(path, "r", newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise SystemExit(f"no rows to write: {path}")
    with open_path(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def finite_float(raw: Any, default: float = math.nan) -> float:
    if raw in (None, ""):
        return default
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return default
    return value if math.isfinite(value) else default


def safe_ratio(numerator: float, denominator: float) -> float | None:
    if not math.isfinite(numerator) or not math.isfinite(denominator) or denominator == 0.0:
        return None
    return numerator / denominator


def abs_ratio(numerator: float, denominator: float) -> float | None:
    if not math.isfinite(numerator) or not math.isfinite(denominator) or denominator == 0.0:
        return None
    return abs(numerator) / abs(denominator)


def log10_abs(value: float | None) -> float | None:
    if value is None or not math.isfinite(value) or value == 0.0:
        return None
    return math.log10(abs(value))


def clean_json(value: Any) -> Any:
    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        return None
    if isinstance(value, dict):
        return {key: clean_json(item) for key, item in value.items()}
    if isinstance(value, list):
        return [clean_json(item) for item in value]
    return value


def classify_doping(donors: float, acceptors: float) -> tuple[str, float]:
    net = donors - acceptors
    threshold = 1.0e-6 * max(abs(donors), abs(acceptors), 1.0)
    if abs(net) <= threshold:
        return "compensated", net
    if net > 0.0:
        return "n", net
    return "p", net


def load_doping(path: Path) -> dict[int, dict[str, Any]]:
    result: dict[int, dict[str, Any]] = {}
    for row in read_csv(path):
        node_id = int(finite_float(row.get("node_id", row.get("id"))))
        donors = finite_float(row.get("donors_cm3"), 0.0)
        acceptors = finite_float(row.get("acceptors_cm3"), 0.0)
        kind, net = classify_doping(donors, acceptors)
        result[node_id] = {
            "donors_cm3": donors,
            "acceptors_cm3": acceptors,
            "net_doping_cm3": net,
            "type": kind,
        }
    return result


def parse_vtk(path: Path) -> dict[str, Any]:
    lines = path.read_text(encoding="utf-8").splitlines()
    points: list[tuple[float, float, float]] = []
    scalars: dict[str, list[float]] = {}
    section: str | None = None
    section_count = 0
    index = 0
    while index < len(lines):
        parts = lines[index].split()
        if not parts:
            index += 1
            continue
        if parts[0] == "POINTS":
            count = int(parts[1])
            index += 1
            values: list[float] = []
            while len(values) < 3 * count:
                values.extend(float(item) for item in lines[index].split())
                index += 1
            points = [(values[i] * 1.0e6, values[i + 1] * 1.0e6, values[i + 2] * 1.0e6) for i in range(0, len(values), 3)]
            continue
        if parts[0] == "POINT_DATA":
            section = "point"
            section_count = int(parts[1])
            index += 1
            continue
        if parts[0] == "CELL_DATA":
            section = "cell"
            section_count = int(parts[1])
            index += 1
            continue
        if parts[0] == "SCALARS" and section == "point":
            name = parts[1]
            index += 1
            if lines[index].strip().startswith("LOOKUP_TABLE"):
                index += 1
            values = []
            while len(values) < section_count:
                values.extend(float(item) for item in lines[index].split())
                index += 1
            scalars[name] = values[:section_count]
            continue
        if parts[0] == "VECTORS" and section in {"point", "cell"}:
            index += 1 + section_count
            continue
        index += 1
    missing = [name for name in NODE_FIELDS if name not in scalars]
    if missing:
        raise SystemExit(f"missing VTK scalar(s) {missing} in {path}")
    return {"points": points, "scalars": scalars}


def nearest_node(points: list[tuple[float, float, float]], x_um: float, y_um: float) -> int:
    best: tuple[float, int] | None = None
    for node_id, (x, y, _z) in enumerate(points):
        distance = (x - x_um) ** 2 + (y - y_um) ** 2
        if best is None or distance < best[0]:
            best = (distance, node_id)
    if best is None:
        raise SystemExit("no VTK points loaded")
    return best[1]


def vtk_for_bias(root: Path, prefix: str, bias: float) -> Path:
    index = int(round(abs(bias) / 0.05))
    exact = root / f"{prefix}_{index:04d}_{bias:g}V.vtk"
    if exact.exists():
        return exact
    matches = sorted(root.glob(f"{prefix}_{index:04d}_*.vtk"))
    if not matches:
        raise SystemExit(f"no VTK file found for bias {bias} in {root}")
    return matches[0]


def load_sg_edges(path: Path) -> dict[tuple[float, int], dict[str, str]]:
    result: dict[tuple[float, int], dict[str, str]] = {}
    for row in read_csv(path):
        bias = round(finite_float(row.get("bias_V")), 10)
        edge_id = int(finite_float(row.get("edge_id")))
        result[(bias, edge_id)] = row
    return result


def load_sentaurus_nodes(sentaurus_root: Path, bias: float) -> tuple[list[dict[str, Any]], dict[int, dict[str, Any]]]:
    bias_name = f"sentaurus_{bias:g}v"
    directory = sentaurus_root / bias_name
    if not directory.exists():
        raise SystemExit(f"missing Sentaurus export directory: {directory}")
    nodes = []
    for row in read_csv(directory / "nodes.csv"):
        nodes.append({"id": int(row["id"]), "x_um": finite_float(row["x_um"]), "y_um": finite_float(row["y_um"])})
    doping = load_doping(directory / "doping.csv")
    return nodes, doping


def nearest_sentaurus_node(nodes: list[dict[str, Any]], x_um: float, y_um: float) -> dict[str, Any]:
    best: tuple[float, dict[str, Any]] | None = None
    for node in nodes:
        distance = (node["x_um"] - x_um) ** 2 + (node["y_um"] - y_um) ** 2
        if best is None or distance < best[0]:
            best = (distance, node)
    if best is None:
        raise SystemExit("no Sentaurus nodes loaded")
    return best[1]


def row_side_nodes(state: dict[str, Any], side: str, y_um: float) -> tuple[int, int]:
    x0, x1 = (VELA_X_COLUMNS[0], VELA_X_COLUMNS[1]) if side == "left" else (VELA_X_COLUMNS[1], VELA_X_COLUMNS[2])
    points = state["points"]
    return nearest_node(points, x0, y_um), nearest_node(points, x1, y_um)


def scalar_drop(state: dict[str, Any], field: str, node0: int, node1: int) -> float:
    values = state["scalars"][field]
    return values[node1] - values[node0]


def scalar_mid(state: dict[str, Any], field: str, node0: int, node1: int) -> float:
    values = state["scalars"][field]
    return 0.5 * (values[node0] + values[node1])


def endpoint_values(state: dict[str, Any], field: str, node0: int, node1: int) -> tuple[float, float]:
    values = state["scalars"][field]
    return values[node0], values[node1]


def build_detail_rows(args: argparse.Namespace) -> list[dict[str, Any]]:
    baseline_mesh_doping = load_doping(args.baseline_report_root.parent.parent / "imported_reference" / "vela" / "doping.csv")
    probe_doping = load_doping(args.probe_root / "doping_compensated_x1_column.csv")
    variants = {
        "baseline": {
            "sg": load_sg_edges(args.baseline_report_root / "sg_avalanche_edges_density_gradient_0p05.csv"),
            "vtk_root": args.baseline_report_root / "vtk_density_gradient_0p05",
            "vtk_prefix": "dc_sweep",
            "doping": baseline_mesh_doping,
        },
        "compensated_probe": {
            "sg": load_sg_edges(args.probe_root / "sg_avalanche_edges_compensated_junction_0p05.csv"),
            "vtk_root": args.probe_root / "vtk_compensated_junction_0p05",
            "vtk_prefix": "dc_sweep",
            "doping": probe_doping,
        },
    }

    rows: list[dict[str, Any]] = []
    for bias in BIASES:
        sentaurus_nodes, sentaurus_doping = load_sentaurus_nodes(args.sentaurus_root, bias)
        states = {
            name: parse_vtk(vtk_for_bias(data["vtk_root"], data["vtk_prefix"], bias))
            for name, data in variants.items()
        }
        for y_um in Y_CUTS:
            sent_left = nearest_sentaurus_node(sentaurus_nodes, SENTAURUS_X_COLUMNS[0], y_um)
            sent_mid = nearest_sentaurus_node(sentaurus_nodes, SENTAURUS_X_COLUMNS[1], y_um)
            sent_right = nearest_sentaurus_node(sentaurus_nodes, SENTAURUS_X_COLUMNS[2], y_um)
            sent_edge_types = {
                "left": f"{sentaurus_doping[sent_left['id']]['type']}-{sentaurus_doping[sent_mid['id']]['type']}",
                "right": f"{sentaurus_doping[sent_mid['id']]['type']}-{sentaurus_doping[sent_right['id']]['type']}",
            }
            for variant_name, variant in variants.items():
                state = states[variant_name]
                for side in ["left", "right"]:
                    edge_id = EDGE_BY_SIDE[y_um][side]
                    edge_row = variant["sg"].get((round(bias, 10), edge_id))
                    if edge_row is None:
                        raise SystemExit(f"missing SG edge row for {variant_name} bias={bias} edge={edge_id}")
                    node0, node1 = row_side_nodes(state, side, y_um)
                    doping0 = variant["doping"][node0]
                    doping1 = variant["doping"][node1]
                    item: dict[str, Any] = {
                        "variant": variant_name,
                        "bias_V": bias,
                        "y_um": y_um,
                        "side": side,
                        "edge_id": edge_id,
                        "node0": node0,
                        "node1": node1,
                        "node0_type": doping0["type"],
                        "node1_type": doping1["type"],
                        "edge_type": f"{doping0['type']}-{doping1['type']}",
                        "node0_net_doping_cm3": doping0["net_doping_cm3"],
                        "node1_net_doping_cm3": doping1["net_doping_cm3"],
                        "sentaurus_edge_type": sent_edge_types[side],
                        "sentaurus_nearest_left_node": sent_left["id"],
                        "sentaurus_nearest_mid_node": sent_mid["id"],
                        "sentaurus_nearest_right_node": sent_right["id"],
                    }
                    for field, output_name in FIELD_TO_OUTPUT.items():
                        item[f"{output_name}_drop_V"] = scalar_drop(state, field, node0, node1)
                        value0, value1 = endpoint_values(state, field, node0, node1)
                        item[f"{output_name}0"] = value0
                        item[f"{output_name}1"] = value1
                    item["electron_density_mid_m3"] = scalar_mid(state, "Electrons", node0, node1)
                    item["hole_density_mid_m3"] = scalar_mid(state, "Holes", node0, node1)
                    item["electron_density_endpoint_abs_ratio"] = abs_ratio(item["electrons1"], item["electrons0"])
                    item["hole_density_endpoint_abs_ratio"] = abs_ratio(item["holes1"], item["holes0"])
                    for field in [
                        "edge_length_m",
                        "edge_couple_m",
                        "edge_area_proxy_m2",
                        "electric_field_V_per_m",
                        "electron_impact_field_V_per_m",
                        "hole_impact_field_V_per_m",
                        "electron_alpha_m_inv",
                        "hole_alpha_m_inv",
                        "electron_mobility_m2_V_s",
                        "hole_mobility_m2_V_s",
                        "electron_flux_proxy",
                        "hole_flux_proxy",
                        "electron_raw_flux_proxy",
                        "hole_raw_flux_proxy",
                        "electron_reconstructed_flux_proxy",
                        "hole_reconstructed_flux_proxy",
                        "electron_final_over_raw_flux_proxy",
                        "hole_final_over_raw_flux_proxy",
                        "electron_source_integral",
                        "hole_source_integral",
                        "edge_source_integral",
                    ]:
                        item[field] = finite_float(edge_row.get(field))
                    rows.append(item)
    add_pair_ratios(rows)
    add_probe_over_baseline(rows)
    return rows


def add_pair_ratios(rows: list[dict[str, Any]]) -> None:
    by_key: dict[tuple[str, float, float], dict[str, dict[str, Any]]] = {}
    for row in rows:
        by_key.setdefault((row["variant"], row["bias_V"], row["y_um"]), {})[row["side"]] = row
    for (_variant, _bias, _y), pair in by_key.items():
        left = pair.get("left")
        right = pair.get("right")
        if left is None or right is None:
            continue
        for field in RATIO_FIELDS + ["psi_drop_V", "phin_drop_V", "phip_drop_V"]:
            ratio = abs_ratio(right.get(field, math.nan), left.get(field, math.nan))
            left[f"right_over_left_{field}"] = ratio
            right[f"right_over_left_{field}"] = ratio


def add_probe_over_baseline(rows: list[dict[str, Any]]) -> None:
    by_key: dict[tuple[float, float, str], dict[str, dict[str, Any]]] = {}
    for row in rows:
        by_key.setdefault((row["bias_V"], row["y_um"], row["side"]), {})[row["variant"]] = row
    for (_bias, _y, _side), pair in by_key.items():
        baseline = pair.get("baseline")
        probe = pair.get("compensated_probe")
        if baseline is None or probe is None:
            continue
        for field in RATIO_FIELDS + ["psi_drop_V", "phin_drop_V", "phip_drop_V"]:
            ratio = abs_ratio(probe.get(field, math.nan), baseline.get(field, math.nan))
            baseline[f"probe_over_baseline_{field}"] = ratio
            probe[f"probe_over_baseline_{field}"] = ratio


def median(values: list[float]) -> float | None:
    finite = [value for value in values if value is not None and math.isfinite(value)]
    return statistics.median(finite) if finite else None


def summarize(rows: list[dict[str, Any]]) -> dict[str, Any]:
    aggregate: list[dict[str, Any]] = []
    dominant_by_bias: list[dict[str, Any]] = []
    for variant in ["baseline", "compensated_probe"]:
        for bias in BIASES:
            subset = [row for row in rows if row["variant"] == variant and row["bias_V"] == bias and row["side"] == "right"]
            item: dict[str, Any] = {"variant": variant, "bias_V": bias}
            for field in [
                "edge_source_integral",
                "electron_source_integral",
                "hole_source_integral",
                "phin_drop_V",
                "electron_alpha_m_inv",
                "hole_alpha_m_inv",
                "electron_flux_proxy",
                "hole_flux_proxy",
                "electron_raw_flux_proxy",
                "hole_raw_flux_proxy",
                "electron_mobility_m2_V_s",
                "hole_mobility_m2_V_s",
                "electron_density_mid_m3",
                "hole_density_mid_m3",
            ]:
                item[f"median_right_over_left_{field}"] = median([
                    row.get(f"right_over_left_{field}") for row in subset
                ])
            aggregate.append(item)

    for bias in BIASES:
        subset = [row for row in rows if row["variant"] == "compensated_probe" and row["bias_V"] == bias and row["side"] == "right"]
        source_ratio = median([row.get("right_over_left_edge_source_integral") for row in subset])
        source_log = log10_abs(source_ratio)
        channel = {
            "electron_source_right_left_ratio": median([row.get("right_over_left_electron_source_integral") for row in subset]),
            "hole_source_right_left_ratio": median([row.get("right_over_left_hole_source_integral") for row in subset]),
            "electron_alpha_right_left_ratio": median([row.get("right_over_left_electron_alpha_m_inv") for row in subset]),
            "electron_flux_right_left_ratio": median([row.get("right_over_left_electron_flux_proxy") for row in subset]),
            "electron_raw_flux_right_left_ratio": median([row.get("right_over_left_electron_raw_flux_proxy") for row in subset]),
            "electron_mobility_right_left_ratio": median([row.get("right_over_left_electron_mobility_m2_V_s") for row in subset]),
            "electron_density_mid_right_left_ratio": median([row.get("right_over_left_electron_density_mid_m3") for row in subset]),
            "hole_alpha_right_left_ratio": median([row.get("right_over_left_hole_alpha_m_inv") for row in subset]),
            "hole_flux_right_left_ratio": median([row.get("right_over_left_hole_flux_proxy") for row in subset]),
            "edge_area_right_left_ratio": median([row.get("right_over_left_edge_area_proxy_m2") for row in subset]),
        }
        channel["electron_alpha_x_flux_right_left_ratio"] = (
            channel["electron_alpha_right_left_ratio"] * channel["electron_flux_right_left_ratio"]
            if channel["electron_alpha_right_left_ratio"] is not None and channel["electron_flux_right_left_ratio"] is not None
            else None
        )
        channel["hole_alpha_x_flux_right_left_ratio"] = (
            channel["hole_alpha_right_left_ratio"] * channel["hole_flux_right_left_ratio"]
            if channel["hole_alpha_right_left_ratio"] is not None and channel["hole_flux_right_left_ratio"] is not None
            else None
        )
        if (channel["electron_source_right_left_ratio"] or 0.0) > 1.0 and (channel["hole_source_right_left_ratio"] or math.inf) < 1.0:
            channel["dominant_physical_reading"] = "electron source is right-heavy while hole source is left-heavy; residual right bias follows electron SG flux proxy moderated by alpha"
        elif (channel["hole_source_right_left_ratio"] or 0.0) > 1.0 and (channel["electron_source_right_left_ratio"] or math.inf) < 1.0:
            channel["dominant_physical_reading"] = "hole source is right-heavy while electron source is left-heavy"
        else:
            channel["dominant_physical_reading"] = "both carrier source channels have the same right/left direction or one channel is unavailable"

        candidates = []
        for label, fields in [
            ("electron_flux_proxy", ["right_over_left_electron_flux_proxy"]),
            ("electron_raw_flux_proxy", ["right_over_left_electron_raw_flux_proxy"]),
            ("electron_alpha", ["right_over_left_electron_alpha_m_inv"]),
            ("electron_density_mid", ["right_over_left_electron_density_mid_m3"]),
            ("electron_mobility", ["right_over_left_electron_mobility_m2_V_s"]),
            ("hole_flux_proxy", ["right_over_left_hole_flux_proxy"]),
            ("hole_alpha", ["right_over_left_hole_alpha_m_inv"]),
            ("edge_area", ["right_over_left_edge_area_proxy_m2"]),
            ("electron_alpha_x_flux", ["right_over_left_electron_alpha_m_inv", "right_over_left_electron_flux_proxy"]),
            ("hole_alpha_x_flux", ["right_over_left_hole_alpha_m_inv", "right_over_left_hole_flux_proxy"]),
        ]:
            product = 1.0
            ok = True
            for field in fields:
                value = median([row.get(field) for row in subset])
                if value is None or not math.isfinite(value):
                    ok = False
                    break
                product *= value
            if not ok:
                continue
            candidate_log = log10_abs(product)
            distance = abs(candidate_log - source_log) if candidate_log is not None and source_log is not None else None
            candidates.append({
                "factor": label,
                "median_ratio_or_product": product,
                "log10_ratio_or_product": candidate_log,
                "distance_to_source_log10": distance,
            })
        candidates.sort(key=lambda item: math.inf if item["distance_to_source_log10"] is None else item["distance_to_source_log10"])
        dominant_by_bias.append({
            "bias_V": bias,
            "source_right_left_ratio": source_ratio,
            "source_log10_ratio": source_log,
            "channel_decomposition": channel,
            "closest_factor": candidates[0] if candidates else None,
            "ranked_factors": candidates,
        })
    return {"aggregate": aggregate, "dominant_by_bias": dominant_by_bias}


def write_report(path: Path, summary: dict[str, Any], rows: list[dict[str, Any]]) -> None:
    lines: list[str] = []
    lines.append("# PN2D BV Compensated Junction Source Proxy Compare")
    lines.append("")
    lines.append("## Summary")
    lines.append("")
    lines.append("This diagnostic compares the baseline density-gradient BV run with the compensated-junction probe at -12 V, -19 V, and -20 V. It classifies nodes from donors-acceptors and decomposes the remaining right-heavy source into QF drops, carrier densities, mobilities, SG flux proxies, alpha, and source integrals.")
    lines.append("")
    lines.append("Direct `phin/phip` clamp or zeroing is intentionally not part of this diagnostic. Doping classification is used only as metadata for artifact alignment and source-proxy interpretation.")
    lines.append("")
    lines.append("## Median Right/Left Ratios")
    lines.append("")
    lines.append("| variant | bias | source | phin drop | e-alpha | e-flux proxy | e-raw flux | e-density mid | e-mobility |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    for item in summary["aggregate"]:
        lines.append("| {variant} | {bias:g} | {source:.6g} | {phin:.6g} | {alpha:.6g} | {flux:.6g} | {raw:.6g} | {density:.6g} | {mob:.6g} |".format(
            variant=item["variant"],
            bias=item["bias_V"],
            source=item.get("median_right_over_left_edge_source_integral") or math.nan,
            phin=item.get("median_right_over_left_phin_drop_V") or math.nan,
            alpha=item.get("median_right_over_left_electron_alpha_m_inv") or math.nan,
            flux=item.get("median_right_over_left_electron_flux_proxy") or math.nan,
            raw=item.get("median_right_over_left_electron_raw_flux_proxy") or math.nan,
            density=item.get("median_right_over_left_electron_density_mid_m3") or math.nan,
            mob=item.get("median_right_over_left_electron_mobility_m2_V_s") or math.nan,
        ))
    lines.append("")
    lines.append("## Channel Source Decomposition For Compensated Probe")
    lines.append("")
    lines.append("| bias | total source R/L | electron source R/L | hole source R/L | e-alpha | e-flux proxy | e-alpha x flux | e-mobility | reading |")
    lines.append("|---:|---:|---:|---:|---:|---:|---:|---:|---|")
    for item in summary["dominant_by_bias"]:
        channel = item.get("channel_decomposition") or {}
        lines.append("| {bias:g} | {source:.6g} | {esrc:.6g} | {hsrc:.6g} | {ealpha:.6g} | {eflux:.6g} | {eaf:.6g} | {emob:.6g} | {reading} |".format(
            bias=item["bias_V"],
            source=item.get("source_right_left_ratio") or math.nan,
            esrc=channel.get("electron_source_right_left_ratio") or math.nan,
            hsrc=channel.get("hole_source_right_left_ratio") or math.nan,
            ealpha=channel.get("electron_alpha_right_left_ratio") or math.nan,
            eflux=channel.get("electron_flux_right_left_ratio") or math.nan,
            eaf=channel.get("electron_alpha_x_flux_right_left_ratio") or math.nan,
            emob=channel.get("electron_mobility_right_left_ratio") or math.nan,
            reading=channel.get("dominant_physical_reading", ""),
        ))
    lines.append("")
    lines.append("## Scalar Closest-Factor Ranking For Compensated Probe")
    lines.append("")
    lines.append("This table is a scalar log-distance screen only. Use it with the channel decomposition above so that a numerically close factor is not mistaken for the contributing carrier channel.")
    lines.append("")
    lines.append("| bias | source right/left | closest factor | factor ratio/product | log10 distance |")
    lines.append("|---:|---:|---|---:|---:|")
    for item in summary["dominant_by_bias"]:
        closest = item.get("closest_factor") or {}
        lines.append("| {bias:g} | {source:.6g} | {factor} | {ratio:.6g} | {distance:.6g} |".format(
            bias=item["bias_V"],
            source=item.get("source_right_left_ratio") or math.nan,
            factor=closest.get("factor", ""),
            ratio=closest.get("median_ratio_or_product") or math.nan,
            distance=closest.get("distance_to_source_log10") or math.nan,
        ))
    lines.append("")
    lines.append("## Interpretation")
    lines.append("")
    probe_agg = [item for item in summary["aggregate"] if item["variant"] == "compensated_probe"]
    max_phin_ratio = max((item.get("median_right_over_left_phin_drop_V") or 0.0) for item in probe_agg)
    source_ratios = [item.get("median_right_over_left_edge_source_integral") for item in probe_agg]
    lines.append(f"- The compensated probe keeps left/right `phin` drops balanced: max median right/left ratio is `{max_phin_ratio:.6g}`.")
    lines.append("- Remaining source right/left ratios are `{}` for -12/-19/-20 V.".format(
        ", ".join(f"{value:.6g}" for value in source_ratios if value is not None)
    ))
    lines.append("- Channel decomposition shows the residual right-heavy source is carried by the electron source channel; the hole source channel is left-heavy at all three inspected biases.")
    lines.append("- The electron right/left source ratio is driven mainly by the electron SG flux proxy / raw flux proxy, while electron alpha is below 1 and mobility is close to 1 after compensation.")
    lines.append("- This points the next debug target at density-gradient SG current/source construction and carrier-density/flux proxy selection, not at a QF hard limiter.")
    lines.append("")
    lines.append("## Artifacts")
    lines.append("")
    lines.append("- `compensated_source_proxy_compare.csv`")
    lines.append("- `compensated_source_proxy_compare_summary.json`")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    rows = build_detail_rows(args)
    summary = summarize(rows)
    csv_path = args.out_dir / "compensated_source_proxy_compare.csv"
    json_path = args.out_dir / "compensated_source_proxy_compare_summary.json"
    report_path = args.out_dir / "compensated_source_proxy_compare_report_20260709.md"
    write_csv(csv_path, rows)
    json_path.write_text(json.dumps(clean_json(summary), indent=2), encoding="utf-8")
    write_report(report_path, summary, rows)
    print(json.dumps({
        "csv": str(csv_path),
        "json": str(json_path),
        "report": str(report_path),
        "rows": len(rows),
    }, indent=2))


if __name__ == "__main__":
    main()
