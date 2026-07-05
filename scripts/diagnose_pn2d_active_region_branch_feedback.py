#!/usr/bin/env python3
"""Diagnose PN2D active-region branch/state feedback against Sentaurus support."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
from typing import Any, Iterable


ELEMENTARY_CHARGE_C = 1.602176634e-19

DEFAULT_REPORT_DIR = (
    Path("build-release")
    / "reference_tcad"
    / "pn2d_sentaurus2018_coarse7x3"
    / "reports"
    / "active_region_branch_feedback_20260705"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_REPORT_DIR)
    parser.add_argument("--out-dir", type=Path, default=None)
    parser.add_argument("--support-csv", type=Path, default=None)
    parser.add_argument("--node-compare-csv", type=Path, default=None)
    parser.add_argument("--biases", type=float, nargs="+", default=[-18.0, -20.0])
    parser.add_argument("--top-limit", type=int, default=20)
    parser.add_argument("--minimum-field-V-m", type=float, default=0.0)
    parser.add_argument("--cutline-biases", default="-1,-5,-10,-18")
    parser.add_argument("--cutline-axis", choices=["horizontal_y", "vertical_x"], default="horizontal_y")
    parser.add_argument("--cutline-value-um", type=float, default=0.25)
    parser.add_argument("--cutline-tolerance-um", type=float, default=1.0e-6)
    parser.add_argument("--disable-plots", action="store_true")
    return parser.parse_args()


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        return [
            {
                (key.strip() if key is not None else ""): (value.strip() if isinstance(value, str) else "")
                for key, value in row.items()
            }
            for row in reader
        ]


def write_csv_rows(path: Path, rows: list[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def optional_float(value: Any) -> float | None:
    if value is None:
        return None
    if isinstance(value, str) and value.strip() == "":
        return None
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def finite_float(value: Any, default: float = 0.0) -> float:
    parsed = optional_float(value)
    return default if parsed is None else parsed


def optional_int(value: Any) -> int | None:
    parsed = optional_float(value)
    return None if parsed is None else int(parsed)


def close(lhs: float | None, rhs: float, tol: float = 1.0e-9) -> bool:
    return lhs is not None and abs(lhs - rhs) <= tol


def median(values: Iterable[float | None]) -> float | None:
    clean = sorted(value for value in values if value is not None and math.isfinite(value))
    if not clean:
        return None
    mid = len(clean) // 2
    if len(clean) % 2:
        return clean[mid]
    return 0.5 * (clean[mid - 1] + clean[mid])


def ratio(numerator: float, denominator: float) -> float | None:
    if not math.isfinite(numerator) or not math.isfinite(denominator) or denominator == 0.0:
        return None
    return numerator / denominator


def log10_ratio(numerator: float | None, denominator: float | None) -> float | None:
    if numerator is None or denominator is None or numerator <= 0.0 or denominator <= 0.0:
        return None
    return math.log10(numerator / denominator)


def alpha_flux(electron_flux: float, hole_flux: float, electron_alpha: float, hole_alpha: float) -> float:
    return abs(electron_flux) * electron_alpha + abs(hole_flux) * hole_alpha


def branch_delta(
    vela_psi: float,
    vela_phin: float,
    vela_phip: float,
    sent_psi: float,
    sent_phin: float,
    sent_phip: float,
) -> dict[str, float]:
    vela_psi_minus_phin = vela_psi - vela_phin
    sent_psi_minus_phin = sent_psi - sent_phin
    vela_phip_minus_psi = vela_phip - vela_psi
    sent_phip_minus_psi = sent_phip - sent_psi
    return {
        "vela_psi_minus_phin_V": vela_psi_minus_phin,
        "sent_psi_minus_phin_V": sent_psi_minus_phin,
        "delta_psi_minus_phin_V": vela_psi_minus_phin - sent_psi_minus_phin,
        "vela_phip_minus_psi_V": vela_phip_minus_psi,
        "sent_phip_minus_psi_V": sent_phip_minus_psi,
        "delta_phip_minus_psi_V": vela_phip_minus_psi - sent_phip_minus_psi,
    }


def sentaurus_flux_from_current(current_density: Any) -> float:
    return abs(finite_float(current_density)) * 1.0e4 / ELEMENTARY_CHARGE_C


def edge_scores(row: dict[str, Any]) -> dict[str, float]:
    sent_e_flux = sentaurus_flux_from_current(row.get("sent_e_current"))
    sent_h_flux = sentaurus_flux_from_current(row.get("sent_h_current"))
    sent_e_alpha = finite_float(row.get("sent_e_alpha"))
    sent_h_alpha = finite_float(row.get("sent_h_alpha"))
    vela_e_flux = abs(finite_float(row.get("electron_flux_abs")))
    vela_h_flux = abs(finite_float(row.get("hole_flux_abs")))
    vela_e_alpha = finite_float(row.get("electron_alpha_m_inv"))
    vela_h_alpha = finite_float(row.get("hole_alpha_m_inv"))
    return {
        "sent_e_flux": sent_e_flux,
        "sent_h_flux": sent_h_flux,
        "sent_e_alpha": sent_e_alpha,
        "sent_h_alpha": sent_h_alpha,
        "vela_e_flux": vela_e_flux,
        "vela_h_flux": vela_h_flux,
        "vela_e_alpha": vela_e_alpha,
        "vela_h_alpha": vela_h_alpha,
        "sent_e_alpha_flux": sent_e_alpha * sent_e_flux,
        "sent_h_alpha_flux": sent_h_alpha * sent_h_flux,
        "vela_e_alpha_flux": vela_e_alpha * vela_e_flux,
        "vela_h_alpha_flux": vela_h_alpha * vela_h_flux,
        "mixed_sent_alpha_vela_flux": alpha_flux(vela_e_flux, vela_h_flux, sent_e_alpha, sent_h_alpha),
        "mixed_vela_alpha_sent_flux": alpha_flux(sent_e_flux, sent_h_flux, vela_e_alpha, vela_h_alpha),
    }


def nearest_vela_bias_for_target(rows: list[dict[str, str]], target_bias: float) -> float:
    candidates: dict[float, int] = {}
    for row in rows:
        if not close(optional_float(row.get("nearest_sentaurus_bias_V")), target_bias):
            continue
        bias = optional_float(row.get("bias_V"))
        if bias is not None:
            candidates[bias] = candidates.get(bias, 0) + 1
    if not candidates:
        raise ValueError(f"no rows mapped to nearest_sentaurus_bias_V={target_bias}")
    return min(candidates, key=lambda bias: (abs(bias - target_bias), -candidates[bias], bias))


def rows_for_target(rows: list[dict[str, str]], target_bias: float) -> tuple[float, list[dict[str, str]]]:
    vela_bias = nearest_vela_bias_for_target(rows, target_bias)
    selected = [
        row for row in rows
        if close(optional_float(row.get("nearest_sentaurus_bias_V")), target_bias)
        and close(optional_float(row.get("bias_V")), vela_bias)
    ]
    return vela_bias, selected


def select_support_rows_for_sentaurus_bias(
    rows: list[dict[str, str]],
    target_sentaurus_bias: float,
) -> tuple[float, list[dict[str, str]]]:
    return rows_for_target(rows, target_sentaurus_bias)


def ranked_edge_ids(rows: list[dict[str, str]], score_key: str, limit: int) -> list[int]:
    scored: list[tuple[float, int]] = []
    for row in rows:
        edge_id = optional_int(row.get("edge_id"))
        if edge_id is None:
            continue
        scores = edge_scores(row)
        scored.append((scores[score_key], edge_id))
    scored.sort(key=lambda item: (-item[0], item[1]))
    return [edge_id for score, edge_id in scored[:limit] if score > 0.0]


def support_label(in_sentaurus_top: bool, in_vela_top: bool) -> str:
    if in_sentaurus_top and in_vela_top:
        return "overlap"
    if in_sentaurus_top:
        return "sentaurus_only"
    if in_vela_top:
        return "vela_only"
    return "inactive"


def replay_decomposition(row: dict[str, str]) -> dict[str, float | None]:
    scores = edge_scores(row)
    sent_full = scores["sent_e_alpha_flux"] + scores["sent_h_alpha_flux"]
    vela_full = scores["vela_e_alpha_flux"] + scores["vela_h_alpha_flux"]
    return {
        "sent_full": sent_full,
        "baseline_vela": vela_full,
        "sent_alpha_vela_flux": scores["mixed_sent_alpha_vela_flux"],
        "vela_alpha_sent_flux": scores["mixed_vela_alpha_sent_flux"],
        "baseline_vela_over_sent": ratio(vela_full, sent_full),
        "sent_alpha_vela_flux_over_sent": ratio(scores["mixed_sent_alpha_vela_flux"], sent_full),
        "vela_alpha_sent_flux_over_sent": ratio(scores["mixed_vela_alpha_sent_flux"], sent_full),
    }


def limiting_factor(
    sent_alpha_vela_flux_over_sent_full: float | None,
    vela_alpha_sent_flux_over_sent_full: float | None,
    *,
    threshold: float = 1.0e-6,
) -> str:
    flux_limited = sent_alpha_vela_flux_over_sent_full is not None and sent_alpha_vela_flux_over_sent_full < threshold
    alpha_limited = vela_alpha_sent_flux_over_sent_full is not None and vela_alpha_sent_flux_over_sent_full < threshold
    if flux_limited and alpha_limited:
        return "mixed_alpha_and_flux"
    if flux_limited:
        return "flux_limited"
    if alpha_limited:
        return "alpha_limited"
    return "not_limited"


FIELD_RATIO_BINS = [
    ("lt_1e-6", -math.inf, -6.0),
    ("1e-6_to_1e-3", -6.0, -3.0),
    ("1e-3_to_1e-1", -3.0, -1.0),
    ("1e-1_to_0p5", -1.0, math.log10(0.5)),
    ("0p5_to_2", math.log10(0.5), math.log10(2.0)),
    ("2_to_10", math.log10(2.0), 1.0),
    ("10_to_1e3", 1.0, 3.0),
    ("gte_1e3", 3.0, math.inf),
]

QUANTITY_ORDER = {"psi": 0, "phin": 1, "phip": 2, "n": 3, "p": 4}


def parse_bias_csv(raw: str) -> list[float]:
    return [float(item.strip()) for item in raw.split(",") if item.strip()]


def first_present_float(row: dict[str, Any], names: list[str]) -> float | None:
    for name in names:
        value = optional_float(row.get(name))
        if value is not None:
            return abs(value)
    return None


def electric_field_value(row: dict[str, Any]) -> float | None:
    return first_present_float(row, [
        "electric_field_V_m",
        "electric_field_V_per_m",
        "electric_field_abs_V_m",
        "electric_field_abs_V_per_m",
    ])


def carrier_qf_field_value(row: dict[str, Any], carrier: str) -> float | None:
    if carrier == "electron":
        return first_present_float(row, [
            "electron_qf_field_V_m",
            "electron_qf_field_V_per_m",
            "electron_field_V_m",
            "electron_impact_field_V_per_m",
            "electron_impact_field_V_m",
        ])
    return first_present_float(row, [
        "hole_qf_field_V_m",
        "hole_qf_field_V_per_m",
        "hole_field_V_m",
        "hole_impact_field_V_per_m",
        "hole_impact_field_V_m",
    ])


def bin_label_for_ratio(value: float | None) -> str:
    if value is None or value <= 0.0 or not math.isfinite(value):
        return "invalid"
    log_value = math.log10(value)
    for label, low, high in FIELD_RATIO_BINS:
        if low <= log_value < high:
            return label
    return "invalid"


def active_edges_for_bias(membership_rows: list[dict[str, Any]], target_bias: float) -> list[dict[str, Any]]:
    return [
        row for row in membership_rows
        if close(optional_float(row.get("target_sentaurus_bias_V")), target_bias)
        and (int(row.get("in_sentaurus_top") or 0) or int(row.get("in_vela_top") or 0))
    ]


def build_active_edge_field_ratio_diagnostics(
    membership_rows: list[dict[str, Any]],
    raw_by_bias_edge: dict[tuple[float, int], dict[str, str]],
    target_biases: list[float],
    *,
    minimum_field_v_m: float = 0.0,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    field_rows: list[dict[str, Any]] = []
    histogram_rows: list[dict[str, Any]] = []
    for target_bias in target_biases:
        active = active_edges_for_bias(membership_rows, target_bias)
        for edge in active:
            edge_id = int(edge["edge_id"])
            raw = raw_by_bias_edge.get((target_bias, edge_id), {})
            efield = electric_field_value(raw)
            for carrier in ["electron", "hole"]:
                qf_field = carrier_qf_field_value(raw, carrier)
                ratio_value = ratio(qf_field or 0.0, efield or 0.0) if qf_field is not None and efield is not None else None
                hit = int(
                    qf_field is not None
                    and minimum_field_v_m > 0.0
                    and qf_field < minimum_field_v_m
                )
                field_rows.append({
                    "target_sentaurus_bias_V": target_bias,
                    "vela_bias_V": edge.get("vela_bias_V"),
                    "edge_id": edge_id,
                    "carrier": carrier,
                    "electric_field_V_m": efield,
                    "qf_field_V_m": qf_field,
                    "qf_over_electric_field": ratio_value,
                    "log10_qf_over_electric_field": None if ratio_value is None or ratio_value <= 0.0 else math.log10(ratio_value),
                    "ratio_bin": bin_label_for_ratio(ratio_value),
                    "minimum_field_V_m": minimum_field_v_m,
                    "minimum_field_cutoff_hit": hit,
                    "support_set": edge.get("support_set", ""),
                    "sentaurus_rank": edge.get("sentaurus_rank"),
                    "vela_rank": edge.get("vela_rank"),
                })
        for carrier in ["electron", "hole"]:
            subset = [
                row for row in field_rows
                if close(optional_float(row.get("target_sentaurus_bias_V")), target_bias)
                and row["carrier"] == carrier
            ]
            valid = [row for row in subset if optional_float(row.get("qf_over_electric_field")) is not None]
            hit_count = sum(int(row["minimum_field_cutoff_hit"]) for row in subset)
            denominator = len([row for row in subset if optional_float(row.get("qf_field_V_m")) is not None])
            histogram_rows.append({
                "target_sentaurus_bias_V": target_bias,
                "carrier": carrier,
                "bin": "summary",
                "bin_low_log10": None,
                "bin_high_log10": None,
                "count": len(valid),
                "fraction": 1.0 if valid else None,
                "active_edge_count": len(subset),
                "minimum_field_V_m": minimum_field_v_m,
                "minimum_field_cutoff_hits": hit_count,
                "minimum_field_cutoff_hit_rate": ratio(hit_count, denominator) if denominator else None,
            })
            for label, low, high in FIELD_RATIO_BINS:
                count = sum(1 for row in valid if row.get("ratio_bin") == label)
                histogram_rows.append({
                    "target_sentaurus_bias_V": target_bias,
                    "carrier": carrier,
                    "bin": label,
                    "bin_low_log10": None if not math.isfinite(low) else low,
                    "bin_high_log10": None if not math.isfinite(high) else high,
                    "count": count,
                    "fraction": ratio(count, len(valid)) if valid else None,
                    "active_edge_count": len(subset),
                    "minimum_field_V_m": minimum_field_v_m,
                    "minimum_field_cutoff_hits": hit_count,
                    "minimum_field_cutoff_hit_rate": ratio(hit_count, denominator) if denominator else None,
                })
    return field_rows, histogram_rows


def pearson(xs: list[float], ys: list[float]) -> float | None:
    if len(xs) != len(ys) or len(xs) < 2:
        return None
    mx = sum(xs) / len(xs)
    my = sum(ys) / len(ys)
    cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    vx = sum((x - mx) * (x - mx) for x in xs)
    vy = sum((y - my) * (y - my) for y in ys)
    if vx <= 0.0 or vy <= 0.0:
        return None
    return cov / math.sqrt(vx * vy)


def build_alpha_flux_scatter_diagnostics(
    membership_rows: list[dict[str, Any]],
    raw_by_bias_edge: dict[tuple[float, int], dict[str, str]],
    target_biases: list[float],
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    scatter_rows: list[dict[str, Any]] = []
    summary_rows: list[dict[str, Any]] = []
    for target_bias in target_biases:
        active = active_edges_for_bias(membership_rows, target_bias)
        for edge in active:
            edge_id = int(edge["edge_id"])
            raw = raw_by_bias_edge.get((target_bias, edge_id), {})
            scores = edge_scores(raw)
            combined_flux = scores["vela_e_flux"] + scores["vela_h_flux"]
            combined_alpha = scores["vela_e_alpha"] + scores["vela_h_alpha"]
            scatter_rows.append({
                "target_sentaurus_bias_V": target_bias,
                "vela_bias_V": edge.get("vela_bias_V"),
                "edge_id": edge_id,
                "support_set": edge.get("support_set", ""),
                "combined_flux_m2_s": combined_flux,
                "combined_alpha_m_inv": combined_alpha,
                "combined_alpha_flux": combined_flux * combined_alpha,
                "electron_flux_abs": scores["vela_e_flux"],
                "hole_flux_abs": scores["vela_h_flux"],
                "electron_alpha_m_inv": scores["vela_e_alpha"],
                "hole_alpha_m_inv": scores["vela_h_alpha"],
                "log10_combined_flux": math.log10(combined_flux) if combined_flux > 0.0 else None,
                "log10_combined_alpha": math.log10(combined_alpha) if combined_alpha > 0.0 else None,
            })
        subset = [
            row for row in scatter_rows
            if close(optional_float(row.get("target_sentaurus_bias_V")), target_bias)
            and optional_float(row.get("log10_combined_flux")) is not None
            and optional_float(row.get("log10_combined_alpha")) is not None
        ]
        xs = [float(row["log10_combined_flux"]) for row in subset]
        ys = [float(row["log10_combined_alpha"]) for row in subset]
        corr = pearson(xs, ys)
        summary_rows.append({
            "target_sentaurus_bias_V": target_bias,
            "edge_count": len(subset),
            "pearson_log_alpha_log_flux": corr,
            "anti_correlated": None if corr is None else int(corr < 0.0),
        })
    return scatter_rows, summary_rows


def load_node_compare_records(node_rows: list[dict[str, str]]) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for row in node_rows:
        quantity = normalized_quantity(row.get("quantity", ""))
        bias = optional_float(row.get("bias_V"))
        node_id = optional_int(row.get("node_id"))
        if quantity is None or bias is None or node_id is None:
            continue
        records.append({
            "bias_V": bias,
            "vela_bias_V": optional_float(row.get("vela_bias_V")),
            "node_id": node_id,
            "x_um": optional_float(row.get("x_um")),
            "y_um": optional_float(row.get("y_um")),
            "quantity": quantity,
            "sentaurus_value": optional_float(row.get("sentaurus_value")),
            "vela_value": optional_float(row.get("vela_value_scaled_to_sentaurus_units")),
        })
    return records


def cutline_state_rows(
    records: list[dict[str, Any]],
    target_biases: list[float],
    *,
    axis: str = "horizontal_y",
    axis_value_um: float = 0.25,
    tolerance_um: float = 1.0e-6,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for target_bias in target_biases:
        for record in records:
            if not close(optional_float(record.get("bias_V")), target_bias):
                continue
            x = optional_float(record.get("x_um"))
            y = optional_float(record.get("y_um"))
            if x is None or y is None:
                continue
            if axis == "horizontal_y":
                if abs(y - axis_value_um) > tolerance_um:
                    continue
                coordinate = x
            elif axis == "vertical_x":
                if abs(x - axis_value_um) > tolerance_um:
                    continue
                coordinate = y
            else:
                raise ValueError(f"unknown cutline axis {axis}")
            rows.append({
                "bias_V": target_bias,
                "vela_bias_V": record.get("vela_bias_V"),
                "axis": axis,
                "axis_value_um": axis_value_um,
                "coordinate_um": coordinate,
                "node_id": record["node_id"],
                "x_um": x,
                "y_um": y,
                "quantity": record["quantity"],
                "sentaurus_value": record.get("sentaurus_value"),
                "vela_value": record.get("vela_value"),
                "diff": None if record.get("sentaurus_value") is None or record.get("vela_value") is None else record["vela_value"] - record["sentaurus_value"],
            })
    rows.sort(key=lambda row: (
        float(row["bias_V"]),
        float(row["coordinate_um"]),
        QUANTITY_ORDER.get(str(row["quantity"]), 99),
        int(row["node_id"]),
    ))
    return rows


def plot_bias_token(bias: float) -> str:
    return f"{bias:g}".replace("-", "m").replace(".", "p")


def setup_matplotlib(out_dir: Path):
    mpl_config = out_dir / ".matplotlib"
    mpl_config.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_config))
    import matplotlib.pyplot as plt
    return plt


def write_field_ratio_histogram_plot(hist_rows: list[dict[str, Any]], out_dir: Path) -> dict[str, str]:
    if not hist_rows:
        return {}
    plt = setup_matplotlib(out_dir)
    plots: dict[str, str] = {}
    biases = sorted({float(row["target_sentaurus_bias_V"]) for row in hist_rows})
    labels = [label for label, _, _ in FIELD_RATIO_BINS]
    for bias in biases:
        fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2), constrained_layout=True)
        for ax, carrier in zip(axes, ["electron", "hole"]):
            values = []
            for label in labels:
                row = next((
                    item for item in hist_rows
                    if close(optional_float(item.get("target_sentaurus_bias_V")), bias)
                    and item.get("carrier") == carrier
                    and item.get("bin") == label
                ), None)
                values.append(0 if row is None else int(row.get("count") or 0))
            summary = next((
                item for item in hist_rows
                if close(optional_float(item.get("target_sentaurus_bias_V")), bias)
                and item.get("carrier") == carrier
                and item.get("bin") == "summary"
            ), {})
            hit_rate = summary.get("minimum_field_cutoff_hit_rate")
            ax.bar(range(len(labels)), values, color="#3A6EA5")
            ax.set_xticks(range(len(labels)))
            ax.set_xticklabels(labels, rotation=35, ha="right", fontsize=8)
            ax.set_ylabel("active edge count")
            ax.set_title(f"{carrier} F_qf/E, minField hit={hit_rate if hit_rate is not None else 'n/a'}")
            ax.grid(axis="y", linewidth=0.45, alpha=0.35)
        fig.suptitle(f"Active-edge F_qf/E histogram at {bias:g} V")
        path = out_dir / f"active_edge_fqf_over_e_hist_{plot_bias_token(bias)}V.png"
        fig.savefig(path, dpi=220)
        plt.close(fig)
        plots[f"field_ratio_hist_{bias:g}V_png"] = str(path)
    return plots


def write_alpha_flux_scatter_plots(
    scatter_rows: list[dict[str, Any]],
    corr_rows: list[dict[str, Any]],
    out_dir: Path,
) -> dict[str, str]:
    positive = [
        row for row in scatter_rows
        if optional_float(row.get("combined_flux_m2_s")) not in (None, 0.0)
        and optional_float(row.get("combined_alpha_m_inv")) not in (None, 0.0)
    ]
    if not positive:
        return {}
    plt = setup_matplotlib(out_dir)
    plots: dict[str, str] = {}
    biases = sorted({float(row["target_sentaurus_bias_V"]) for row in positive})
    for bias in biases:
        subset = [row for row in positive if close(optional_float(row.get("target_sentaurus_bias_V")), bias)]
        corr = next((row for row in corr_rows if close(optional_float(row.get("target_sentaurus_bias_V")), bias)), {})
        fig, ax = plt.subplots(figsize=(6.2, 5.0), constrained_layout=True)
        ax.loglog(
            [max(float(row["combined_flux_m2_s"]), 1.0e-300) for row in subset],
            [max(float(row["combined_alpha_m_inv"]), 1.0e-300) for row in subset],
            "o",
            markersize=4,
            alpha=0.78,
            color="#9B3A4D",
        )
        ax.set_xlabel("Vela edge flux support |Gamma_e|+|Gamma_h| (m^-2 s^-1)")
        ax.set_ylabel("Vela edge alpha support alpha_e+alpha_h (m^-1)")
        ax.set_title(f"alpha vs flux at {bias:g} V, r={corr.get('pearson_log_alpha_log_flux')}")
        ax.grid(True, which="both", linewidth=0.45, alpha=0.35)
        path = out_dir / f"active_edge_alpha_flux_scatter_{plot_bias_token(bias)}V.png"
        fig.savefig(path, dpi=220)
        plt.close(fig)
        plots[f"alpha_flux_scatter_{bias:g}V_png"] = str(path)
    return plots


def write_cutline_state_plots(cutline_rows: list[dict[str, Any]], out_dir: Path) -> dict[str, str]:
    if not cutline_rows:
        return {}
    plt = setup_matplotlib(out_dir)
    plots: dict[str, str] = {}
    quantities = ["psi", "phin", "phip", "n", "p"]
    biases = sorted({float(row["bias_V"]) for row in cutline_rows})
    for bias in biases:
        fig, axes = plt.subplots(len(quantities), 1, figsize=(8.0, 10.5), sharex=True, constrained_layout=True)
        for ax, quantity in zip(axes, quantities):
            subset = [
                row for row in cutline_rows
                if close(optional_float(row.get("bias_V")), bias)
                and row.get("quantity") == quantity
            ]
            xs = [float(row["coordinate_um"]) for row in subset]
            sent = [row.get("sentaurus_value") for row in subset]
            vela = [row.get("vela_value") for row in subset]
            if quantity in {"n", "p"}:
                ax.semilogy(xs, [max(float(v or 0.0), 1.0e-300) for v in sent], "-o", label="Sentaurus", markersize=3)
                ax.semilogy(xs, [max(float(v or 0.0), 1.0e-300) for v in vela], "--s", label="Vela", markersize=3)
                ax.set_ylabel(f"{quantity} (cm^-3)")
            else:
                ax.plot(xs, sent, "-o", label="Sentaurus", markersize=3)
                ax.plot(xs, vela, "--s", label="Vela", markersize=3)
                ax.set_ylabel(f"{quantity} (V)")
            ax.grid(True, linewidth=0.45, alpha=0.35)
            ax.legend(loc="best", fontsize=8)
        axes[-1].set_xlabel("cutline coordinate (um)")
        fig.suptitle(f"PN2D cutline state overlay at {bias:g} V")
        path = out_dir / f"cutline_state_overlay_{plot_bias_token(bias)}V.png"
        fig.savefig(path, dpi=220)
        plt.close(fig)
        plots[f"cutline_state_overlay_{bias:g}V_png"] = str(path)
    return plots

def build_support_membership(
    support_rows: list[dict[str, str]],
    target_biases: list[float],
    top_limit: int,
) -> tuple[list[dict[str, Any]], dict[tuple[float, int], dict[str, str]]]:
    output: list[dict[str, Any]] = []
    raw_by_bias_edge: dict[tuple[float, int], dict[str, str]] = {}
    for target_bias in target_biases:
        vela_bias, rows = rows_for_target(support_rows, target_bias)
        sent_top = set(ranked_edge_ids(rows, "sent_e_alpha_flux", top_limit))
        sent_top.update(ranked_edge_ids(rows, "sent_h_alpha_flux", top_limit))
        sent_combined = sorted(
            (
                (edge_scores(row)["sent_e_alpha_flux"] + edge_scores(row)["sent_h_alpha_flux"], optional_int(row.get("edge_id")), row)
                for row in rows
            ),
            key=lambda item: (-(item[0]), item[1] if item[1] is not None else 10**9),
        )
        vela_combined = sorted(
            (
                (edge_scores(row)["vela_e_alpha_flux"] + edge_scores(row)["vela_h_alpha_flux"], optional_int(row.get("edge_id")), row)
                for row in rows
            ),
            key=lambda item: (-(item[0]), item[1] if item[1] is not None else 10**9),
        )
        sent_top = {edge_id for _, edge_id, _ in sent_combined[:top_limit] if edge_id is not None}
        vela_top = {edge_id for _, edge_id, _ in vela_combined[:top_limit] if edge_id is not None}
        sent_rank = {edge_id: rank for rank, (_, edge_id, _) in enumerate(sent_combined, start=1) if edge_id is not None}
        vela_rank = {edge_id: rank for rank, (_, edge_id, _) in enumerate(vela_combined, start=1) if edge_id is not None}
        for _, edge_id, row in sent_combined:
            if edge_id is None:
                continue
            raw_by_bias_edge[(target_bias, edge_id)] = row
            if edge_id not in sent_top and edge_id not in vela_top:
                continue
            scores = edge_scores(row)
            sent_score = scores["sent_e_alpha_flux"] + scores["sent_h_alpha_flux"]
            vela_score = scores["vela_e_alpha_flux"] + scores["vela_h_alpha_flux"]
            output.append({
                "target_sentaurus_bias_V": target_bias,
                "vela_bias_V": vela_bias,
                "edge_id": edge_id,
                "node0": optional_int(row.get("node0")),
                "node1": optional_int(row.get("node1")),
                "edge_class": row.get("edge_class", ""),
                "in_sentaurus_top": int(edge_id in sent_top),
                "in_vela_top": int(edge_id in vela_top),
                "support_set": support_label(edge_id in sent_top, edge_id in vela_top),
                "sentaurus_rank": sent_rank.get(edge_id),
                "vela_rank": vela_rank.get(edge_id),
                "sent_combined_alpha_flux": sent_score,
                "vela_combined_alpha_flux": vela_score,
                "vela_over_sent_combined_alpha_flux": ratio(vela_score, sent_score),
                "source_integral_total": finite_float(row.get("source_integral_total")),
            })
    return output, raw_by_bias_edge


def normalized_quantity(name: str) -> str | None:
    compact = "".join(ch for ch in name.lower() if ch.isalnum())
    if compact in {"potential", "electrostaticpotential"}:
        return "psi"
    if compact in {"edensity", "electrondensity", "electronconcentration", "electrondensity"}:
        return "n"
    if compact in {"hdensity", "holedensity", "holeconcentration"}:
        return "p"
    if compact in {"equasifermi", "electronqf", "electronquasifermi", "electronquasifermipotential"}:
        return "phin"
    if compact in {"hquasifermi", "holeqf", "holequasifermi", "holequasifermipotential"}:
        return "phip"
    return None


def load_node_states(node_rows: list[dict[str, str]]) -> dict[tuple[float, int], dict[str, dict[str, float | None]]]:
    states: dict[tuple[float, int], dict[str, dict[str, float | None]]] = {}
    for row in node_rows:
        bias = optional_float(row.get("bias_V"))
        node_id = optional_int(row.get("node_id"))
        quantity = normalized_quantity(row.get("quantity", ""))
        if bias is None or node_id is None or quantity is None:
            continue
        states.setdefault((bias, node_id), {})[quantity] = {
            "sent": optional_float(row.get("sentaurus_value")),
            "vela": optional_float(row.get("vela_value_scaled_to_sentaurus_units")),
        }
    return states


def node_metric(states: dict[tuple[float, int], dict[str, dict[str, float | None]]], bias: float, node_id: int, quantity: str, side: str) -> float | None:
    return states.get((bias, node_id), {}).get(quantity, {}).get(side)


def support_edges_for_set(membership_rows: list[dict[str, Any]], target_bias: float, support_set: str) -> list[dict[str, Any]]:
    rows = [row for row in membership_rows if close(optional_float(row.get("target_sentaurus_bias_V")), target_bias)]
    if support_set == "sentaurus_top":
        return [row for row in rows if int(row["in_sentaurus_top"])]
    if support_set == "vela_top":
        return [row for row in rows if int(row["in_vela_top"])]
    if support_set == "overlap":
        return [row for row in rows if row["support_set"] == "overlap"]
    if support_set == "sentaurus_only":
        return [row for row in rows if row["support_set"] == "sentaurus_only"]
    if support_set == "vela_only":
        return [row for row in rows if row["support_set"] == "vela_only"]
    if support_set == "union":
        return rows
    raise ValueError(f"unknown support_set {support_set}")


def endpoint_nodes(edges: list[dict[str, Any]]) -> list[int]:
    nodes: set[int] = set()
    for edge in edges:
        for key in ("node0", "node1"):
            node = optional_int(edge.get(key))
            if node is not None:
                nodes.add(node)
    return sorted(nodes)


def branch_state_summary(
    membership_rows: list[dict[str, Any]],
    node_states: dict[tuple[float, int], dict[str, dict[str, float | None]]],
    target_biases: list[float],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for target_bias in target_biases:
        for support_set in ["sentaurus_top", "vela_top", "overlap", "sentaurus_only", "vela_only", "union"]:
            edges = support_edges_for_set(membership_rows, target_bias, support_set)
            nodes = endpoint_nodes(edges)
            vela_psi_minus_phin: list[float | None] = []
            sent_psi_minus_phin: list[float | None] = []
            vela_phip_minus_psi: list[float | None] = []
            sent_phip_minus_psi: list[float | None] = []
            n_log_ratios: list[float | None] = []
            p_log_ratios: list[float | None] = []
            for node_id in nodes:
                sent_psi = node_metric(node_states, target_bias, node_id, "psi", "sent")
                vela_psi = node_metric(node_states, target_bias, node_id, "psi", "vela")
                sent_phin = node_metric(node_states, target_bias, node_id, "phin", "sent")
                vela_phin = node_metric(node_states, target_bias, node_id, "phin", "vela")
                sent_phip = node_metric(node_states, target_bias, node_id, "phip", "sent")
                vela_phip = node_metric(node_states, target_bias, node_id, "phip", "vela")
                sent_n = node_metric(node_states, target_bias, node_id, "n", "sent")
                vela_n = node_metric(node_states, target_bias, node_id, "n", "vela")
                sent_p = node_metric(node_states, target_bias, node_id, "p", "sent")
                vela_p = node_metric(node_states, target_bias, node_id, "p", "vela")
                vela_psi_minus_phin.append(None if vela_psi is None or vela_phin is None else vela_psi - vela_phin)
                sent_psi_minus_phin.append(None if sent_psi is None or sent_phin is None else sent_psi - sent_phin)
                vela_phip_minus_psi.append(None if vela_phip is None or vela_psi is None else vela_phip - vela_psi)
                sent_phip_minus_psi.append(None if sent_phip is None or sent_psi is None else sent_phip - sent_psi)
                n_log_ratios.append(log10_ratio(vela_n, sent_n))
                p_log_ratios.append(log10_ratio(vela_p, sent_p))
            vela_branch_n = median(vela_psi_minus_phin)
            sent_branch_n = median(sent_psi_minus_phin)
            vela_branch_p = median(vela_phip_minus_psi)
            sent_branch_p = median(sent_phip_minus_psi)
            rows.append({
                "target_sentaurus_bias_V": target_bias,
                "support_set": support_set,
                "edge_count": len(edges),
                "endpoint_node_count": len(nodes),
                "edge_ids": ";".join(str(row["edge_id"]) for row in edges),
                "endpoint_nodes": ";".join(str(node_id) for node_id in nodes),
                "median_vela_psi_minus_phin_V": vela_branch_n,
                "median_sent_psi_minus_phin_V": sent_branch_n,
                "median_delta_psi_minus_phin_V": None if vela_branch_n is None or sent_branch_n is None else vela_branch_n - sent_branch_n,
                "median_vela_phip_minus_psi_V": vela_branch_p,
                "median_sent_phip_minus_psi_V": sent_branch_p,
                "median_delta_phip_minus_psi_V": None if vela_branch_p is None or sent_branch_p is None else vela_branch_p - sent_branch_p,
                "median_log10_n_vela_over_sentaurus": median(n_log_ratios),
                "median_log10_p_vela_over_sentaurus": median(p_log_ratios),
            })
    return rows


def replay_summary(
    membership_rows: list[dict[str, Any]],
    raw_by_bias_edge: dict[tuple[float, int], dict[str, str]],
    target_biases: list[float],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for target_bias in target_biases:
        for support_set in ["sentaurus_top", "vela_top", "overlap", "sentaurus_only", "vela_only", "union"]:
            edges = support_edges_for_set(membership_rows, target_bias, support_set)
            totals = {
                "sent_e_alpha_flux": 0.0,
                "sent_h_alpha_flux": 0.0,
                "vela_e_alpha_flux": 0.0,
                "vela_h_alpha_flux": 0.0,
                "mixed_sent_alpha_vela_flux": 0.0,
                "mixed_vela_alpha_sent_flux": 0.0,
                "source_integral_total": 0.0,
            }
            for edge in edges:
                raw = raw_by_bias_edge.get((target_bias, int(edge["edge_id"])))
                if raw is None:
                    continue
                scores = edge_scores(raw)
                for key in totals:
                    if key == "source_integral_total":
                        totals[key] += finite_float(raw.get("source_integral_total"))
                    else:
                        totals[key] += scores[key]
            sent_total = totals["sent_e_alpha_flux"] + totals["sent_h_alpha_flux"]
            vela_total = totals["vela_e_alpha_flux"] + totals["vela_h_alpha_flux"]
            rows.append({
                "target_sentaurus_bias_V": target_bias,
                "support_set": support_set,
                "edge_count": len(edges),
                "sent_combined_alpha_flux": sent_total,
                "vela_combined_alpha_flux": vela_total,
                "mixed_sent_alpha_vela_flux": totals["mixed_sent_alpha_vela_flux"],
                "mixed_vela_alpha_sent_flux": totals["mixed_vela_alpha_sent_flux"],
                "source_integral_total": totals["source_integral_total"],
                "vela_full_over_sent": ratio(vela_total, sent_total),
                "sent_alpha_vela_flux_over_sent": ratio(totals["mixed_sent_alpha_vela_flux"], sent_total),
                "vela_alpha_sent_flux_over_sent": ratio(totals["mixed_vela_alpha_sent_flux"], sent_total),
                "source_integral_over_sent": ratio(totals["source_integral_total"], sent_total),
            })
    return rows


def summarize_json(
    membership_rows: list[dict[str, Any]],
    branch_rows: list[dict[str, Any]],
    replay_rows: list[dict[str, Any]],
    top_limit: int,
) -> dict[str, Any]:
    by_bias: dict[str, dict[str, Any]] = {}
    for row in replay_rows:
        if row["support_set"] != "sentaurus_top":
            continue
        key = f"{float(row['target_sentaurus_bias_V']):.6g}"
        by_bias[key] = {
            "sentaurus_top_edges": row["edge_count"],
            "vela_full_over_sent": row["vela_full_over_sent"],
            "sent_alpha_vela_flux_over_sent": row["sent_alpha_vela_flux_over_sent"],
            "vela_alpha_sent_flux_over_sent": row["vela_alpha_sent_flux_over_sent"],
        }
    return {
        "top_limit": top_limit,
        "support_membership_rows": len(membership_rows),
        "branch_state_rows": len(branch_rows),
        "replay_rows": len(replay_rows),
        "sentaurus_top_replay_by_bias": by_bias,
    }


def main() -> int:
    args = parse_args()
    input_dir = args.input_dir
    out_dir = args.out_dir or input_dir
    support_csv = args.support_csv or input_dir / "input_coarse_current_support_compare_psi_gradient_proxy_1vgrid.csv"
    node_csv = args.node_compare_csv or input_dir / "input_coarse_node_field_compare_psi_gradient_proxy_1vgrid.csv"

    support_rows = read_csv_rows(support_csv)
    node_rows = read_csv_rows(node_csv)
    membership_rows, raw_by_bias_edge = build_support_membership(support_rows, args.biases, args.top_limit)
    node_states = load_node_states(node_rows)
    branch_rows = branch_state_summary(membership_rows, node_states, args.biases)
    replay_rows = replay_summary(membership_rows, raw_by_bias_edge, args.biases)
    field_ratio_rows, field_ratio_hist_rows = build_active_edge_field_ratio_diagnostics(
        membership_rows,
        raw_by_bias_edge,
        args.biases,
        minimum_field_v_m=args.minimum_field_V_m,
    )
    alpha_flux_scatter_rows, alpha_flux_corr_rows = build_alpha_flux_scatter_diagnostics(
        membership_rows,
        raw_by_bias_edge,
        args.biases,
    )
    cutline_biases = parse_bias_csv(args.cutline_biases)
    node_compare_records = load_node_compare_records(node_rows)
    cutline_rows = cutline_state_rows(
        node_compare_records,
        cutline_biases,
        axis=args.cutline_axis,
        axis_value_um=args.cutline_value_um,
        tolerance_um=args.cutline_tolerance_um,
    )

    write_csv_rows(out_dir / "support_membership.csv", membership_rows, [
        "target_sentaurus_bias_V", "vela_bias_V", "edge_id", "node0", "node1", "edge_class",
        "in_sentaurus_top", "in_vela_top", "support_set", "sentaurus_rank", "vela_rank",
        "sent_combined_alpha_flux", "vela_combined_alpha_flux", "vela_over_sent_combined_alpha_flux",
        "source_integral_total",
    ])
    write_csv_rows(out_dir / "support_set_branch_state_summary.csv", branch_rows, [
        "target_sentaurus_bias_V", "support_set", "edge_count", "endpoint_node_count",
        "edge_ids", "endpoint_nodes",
        "median_vela_psi_minus_phin_V", "median_sent_psi_minus_phin_V", "median_delta_psi_minus_phin_V",
        "median_vela_phip_minus_psi_V", "median_sent_phip_minus_psi_V", "median_delta_phip_minus_psi_V",
        "median_log10_n_vela_over_sentaurus", "median_log10_p_vela_over_sentaurus",
    ])
    write_csv_rows(out_dir / "alpha_flux_replay_summary.csv", replay_rows, [
        "target_sentaurus_bias_V", "support_set", "edge_count",
        "sent_combined_alpha_flux", "vela_combined_alpha_flux",
        "mixed_sent_alpha_vela_flux", "mixed_vela_alpha_sent_flux", "source_integral_total",
        "vela_full_over_sent", "sent_alpha_vela_flux_over_sent",
        "vela_alpha_sent_flux_over_sent", "source_integral_over_sent",
    ])
    write_csv_rows(out_dir / "active_edge_field_ratio_edges.csv", field_ratio_rows, [
        "target_sentaurus_bias_V", "vela_bias_V", "edge_id", "carrier",
        "electric_field_V_m", "qf_field_V_m", "qf_over_electric_field",
        "log10_qf_over_electric_field", "ratio_bin", "minimum_field_V_m",
        "minimum_field_cutoff_hit", "support_set", "sentaurus_rank", "vela_rank",
    ])
    write_csv_rows(out_dir / "active_edge_field_ratio_histogram.csv", field_ratio_hist_rows, [
        "target_sentaurus_bias_V", "carrier", "bin", "bin_low_log10", "bin_high_log10",
        "count", "fraction", "active_edge_count", "minimum_field_V_m",
        "minimum_field_cutoff_hits", "minimum_field_cutoff_hit_rate",
    ])
    write_csv_rows(out_dir / "active_edge_alpha_flux_scatter.csv", alpha_flux_scatter_rows, [
        "target_sentaurus_bias_V", "vela_bias_V", "edge_id", "support_set",
        "combined_flux_m2_s", "combined_alpha_m_inv", "combined_alpha_flux",
        "electron_flux_abs", "hole_flux_abs", "electron_alpha_m_inv", "hole_alpha_m_inv",
        "log10_combined_flux", "log10_combined_alpha",
    ])
    write_csv_rows(out_dir / "active_edge_alpha_flux_correlation.csv", alpha_flux_corr_rows, [
        "target_sentaurus_bias_V", "edge_count", "pearson_log_alpha_log_flux", "anti_correlated",
    ])
    write_csv_rows(out_dir / "cutline_state_overlay_samples.csv", cutline_rows, [
        "bias_V", "vela_bias_V", "axis", "axis_value_um", "coordinate_um", "node_id",
        "x_um", "y_um", "quantity", "sentaurus_value", "vela_value", "diff",
    ])
    plots: dict[str, str] = {}
    plot_error = None
    if not args.disable_plots:
        try:
            plots.update(write_field_ratio_histogram_plot(field_ratio_hist_rows, out_dir))
            plots.update(write_alpha_flux_scatter_plots(alpha_flux_scatter_rows, alpha_flux_corr_rows, out_dir))
            plots.update(write_cutline_state_plots(cutline_rows, out_dir))
        except Exception as exc:  # pragma: no cover - plotting availability is environment-dependent
            plot_error = str(exc)
    summary = summarize_json(membership_rows, branch_rows, replay_rows, args.top_limit)
    summary.update({
        "active_edge_field_ratio_rows": len(field_ratio_rows),
        "active_edge_field_ratio_histogram_rows": len(field_ratio_hist_rows),
        "active_edge_alpha_flux_scatter_rows": len(alpha_flux_scatter_rows),
        "active_edge_alpha_flux_correlation_rows": len(alpha_flux_corr_rows),
        "cutline_state_overlay_rows": len(cutline_rows),
        "minimum_field_V_m": args.minimum_field_V_m,
        "cutline": {
            "biases": cutline_biases,
            "axis": args.cutline_axis,
            "axis_value_um": args.cutline_value_um,
            "tolerance_um": args.cutline_tolerance_um,
        },
        "plots": plots,
        "plot_error": plot_error,
    })
    (out_dir / "active_region_branch_feedback_summary.json").write_text(
        json.dumps(summary, indent=2),
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
