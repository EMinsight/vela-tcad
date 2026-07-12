#!/usr/bin/env python3
"""Summarize PN2D full-range IV debug evidence against Sentaurus reference."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any


def read_csv(path: Path | None) -> list[dict[str, str]]:
    if path is None:
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def finite_float(text: str | None) -> float | None:
    if text is None or text == "":
        return None
    try:
        value = float(text)
    except ValueError:
        return None
    return value if math.isfinite(value) else None


def truthy(text: str | None) -> bool:
    if text is None:
        return False
    return text.strip().lower() in {"1", "true", "yes", "y", "converged", "ok"}

def in_bias_window(bias: float, bias_min: float, bias_max: float) -> bool:
    tol = max(abs(bias), abs(bias_min), abs(bias_max), 1.0) * 1.0e-9
    return bias_min - tol <= bias <= bias_max + tol


def finite_pairs(rows: list[dict[str, str]], column: str, scale: float = 1.0) -> list[tuple[float, float]]:
    pairs: list[tuple[float, float]] = []
    for row in rows:
        bias = finite_float(row.get("bias_V"))
        value = finite_float(row.get(column))
        if bias is None or value is None:
            continue
        pairs.append((bias, value * scale))
    return sorted(pairs)


def interpolate_at(pairs: list[tuple[float, float]], bias: float, mode: str) -> float | None:
    if not pairs or bias < pairs[0][0] or bias > pairs[-1][0]:
        return None
    tol = max(abs(bias), 1.0) * 1.0e-12
    for existing_bias, value in pairs:
        if abs(existing_bias - bias) <= tol:
            return value
    for (b0, v0), (b1, v1) in zip(pairs, pairs[1:]):
        if b0 <= bias <= b1 and b1 != b0:
            t = (bias - b0) / (b1 - b0)
            if (
                mode == "log_current"
                and v0 != 0.0
                and v1 != 0.0
                and math.copysign(1.0, v0) == math.copysign(1.0, v1)
            ):
                magnitude = math.exp(math.log(abs(v0)) + t * (math.log(abs(v1)) - math.log(abs(v0))))
                return math.copysign(magnitude, v0)
            return v0 + t * (v1 - v0)
    return None


def nearest_row(rows: list[dict[str, str]], bias: float) -> dict[str, str] | None:
    best: tuple[float, dict[str, str]] | None = None
    for row in rows:
        row_bias = finite_float(row.get("bias_V"))
        if row_bias is None:
            continue
        delta = abs(row_bias - bias)
        if best is None or delta < best[0]:
            best = (delta, row)
    if best is None:
        return None
    return best[1] if best[0] <= max(abs(bias), 1.0) * 1.0e-9 else None


def aligned_points(args: argparse.Namespace, reference_rows: list[dict[str, str]], candidate_rows: list[dict[str, str]]) -> list[dict[str, Any]]:
    reference_pairs = finite_pairs(reference_rows, args.reference_column, 1.0)
    candidate_pairs = finite_pairs(candidate_rows, args.candidate_column, args.candidate_scale)
    points: list[dict[str, Any]] = []
    for bias, reference_value in reference_pairs:
        if not in_bias_window(bias, args.bias_min, args.bias_max):
            continue
        candidate_value = interpolate_at(candidate_pairs, bias, args.interpolation)
        if candidate_value is None:
            continue
        if reference_value == 0.0 or candidate_value == 0.0:
            continue
        order_error = abs(math.log10(abs(candidate_value / reference_value)))
        rel_error = abs(candidate_value - reference_value) / max(abs(reference_value), 1.0e-300)
        points.append({
            "bias_V": bias,
            "reference_current": reference_value,
            "candidate_current": candidate_value,
            "orders_of_magnitude_error": order_error,
            "relative_error": rel_error,
            "candidate_row": nearest_row(candidate_rows, bias),
        })
    return points


def first_trend_mismatch(points: list[dict[str, Any]]) -> dict[str, Any] | None:
    for left, right in zip(points, points[1:]):
        ref_delta = right["reference_current"] - left["reference_current"]
        cand_delta = right["candidate_current"] - left["candidate_current"]
        eps = 1.0e-300
        if abs(ref_delta) <= eps or abs(cand_delta) <= eps:
            continue
        if math.copysign(1.0, ref_delta) != math.copysign(1.0, cand_delta):
            return {
                "bias_V": right["bias_V"],
                "reference_delta": ref_delta,
                "candidate_delta": cand_delta,
            }
    return None


def first_non_reltol(candidate_rows: list[dict[str, str]], bias_min: float, bias_max: float) -> dict[str, Any] | None:
    out: list[tuple[float, dict[str, str]]] = []
    for row in candidate_rows:
        bias = finite_float(row.get("bias_V"))
        if bias is None or not in_bias_window(bias, bias_min, bias_max):
            continue
        reason = row.get("newton_convergence_reason", "").strip()
        if reason and reason != "reltol":
            out.append((bias, row))
    if not out:
        return None
    bias, row = sorted(out)[0]
    return {
        "bias_V": bias,
        "newton_convergence_reason": row.get("newton_convergence_reason", ""),
        "newton_iterations": finite_float(row.get("newton_iterations")),
    }


def last_converged_bias(candidate_rows: list[dict[str, str]], bias_min: float, bias_max: float) -> float | None:
    biases: list[float] = []
    for row in candidate_rows:
        bias = finite_float(row.get("bias_V"))
        if bias is None or not in_bias_window(bias, bias_min, bias_max):
            continue
        if truthy(row.get("converged")):
            biases.append(bias)
    return max(biases) if biases else None


def terminal_current_consistency(rows: list[dict[str, str]], bias_min: float, bias_max: float) -> dict[str, Any]:
    groups: dict[float, list[dict[str, str]]] = {}
    for row in rows:
        bias = finite_float(row.get("bias_V"))
        current = finite_float(row.get("current_total_A_per_um"))
        if bias is None or current is None:
            continue
        if not in_bias_window(bias, bias_min, bias_max):
            continue
        groups.setdefault(bias, []).append(row)
    worst: dict[str, Any] | None = None
    for bias, group in groups.items():
        currents = [finite_float(row.get("current_total_A_per_um")) for row in group]
        finite = [value for value in currents if value is not None]
        if not finite:
            continue
        abs_sum = abs(sum(finite))
        max_abs = max(abs(value) for value in finite)
        ratio = abs_sum / max(max_abs, 1.0e-300)
        item = {
            "bias_V": bias,
            "worst_bias_V": bias,
            "contacts": [row.get("contact", "") for row in group],
            "abs_current_sum_A_per_um": abs_sum,
            "max_abs_current_A_per_um": max_abs,
            "abs_sum_over_max_abs": ratio,
        }
        if worst is None or ratio > worst["abs_sum_over_max_abs"]:
            worst = item
    return worst or {"available": False}


def top_contact_edges(rows: list[dict[str, str]], focus_bias: float | None, limit: int) -> list[dict[str, Any]]:
    if focus_bias is None:
        return []
    candidates: list[dict[str, Any]] = []
    for row in rows:
        bias = finite_float(row.get("bias_V"))
        current = finite_float(row.get("current_total_A_per_um") or row.get("current_total"))
        if bias is None or current is None:
            continue
        if abs(bias - focus_bias) > max(abs(focus_bias), 1.0) * 1.0e-9:
            continue
        phin0 = finite_float(row.get("phin0_V") or row.get("phin0"))
        phin1 = finite_float(row.get("phin1_V") or row.get("phin1"))
        item = {
            "bias_V": bias,
            "contact": row.get("current_contact") or row.get("contact") or "",
            "edge_id": row.get("edge_id", ""),
            "node0": row.get("node0", ""),
            "node1": row.get("node1", ""),
            "current_total_A_per_um": current,
            "abs_current_total_A_per_um": abs(current),
        }
        if phin0 is not None and phin1 is not None:
            item["phin_drop_V"] = phin0 - phin1
        candidates.append(item)
    candidates.sort(key=lambda item: item["abs_current_total_A_per_um"], reverse=True)
    return candidates[:limit]


def markdown_report(summary: dict[str, Any]) -> str:
    first_high = summary.get("first_high_error") or {}
    first_non = summary.get("first_non_reltol") or {}
    lines = [
        "# PN2D IV Full-Range Debug Summary",
        "",
        f"Points compared: {summary.get('points_compared', 0)}",
        f"Max order error: {summary.get('max_orders_of_magnitude_error')}",
        f"Max relative error: {summary.get('max_relative_error')}",
        f"First high-error bias V: {first_high.get('bias_V')}",
        f"First non-reltol bias V: {first_non.get('bias_V')}",
        f"Last converged bias V: {summary.get('last_converged_bias_V')}",
        "",
        "Top contact-edge contributors:",
    ]
    for edge in summary.get("top_contact_edges_by_abs_current", []):
        lines.append(
            "- bias={bias_V} contact={contact} edge={edge_id} current={current_total_A_per_um}".format(**edge)
        )
    lines.append("")
    return "\n".join(lines)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--terminal-balance", type=Path)
    parser.add_argument("--contact-edge", type=Path)
    parser.add_argument("--newton-history", type=Path)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-md", type=Path)
    parser.add_argument("--reference-column", default="current_total")
    parser.add_argument("--candidate-column", default="current_total_A_per_um")
    parser.add_argument("--candidate-scale", type=float, default=-1.0)
    parser.add_argument("--bias-min", type=float, default=0.0)
    parser.add_argument("--bias-max", type=float, default=10.0)
    parser.add_argument("--interpolation", choices=["linear", "log_current"], default="log_current")
    parser.add_argument("--max-orders-threshold", type=float, default=0.3)
    parser.add_argument("--top-edges", type=int, default=5)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    reference_rows = read_csv(args.reference)
    candidate_rows = read_csv(args.candidate)
    terminal_rows = read_csv(args.terminal_balance)
    contact_rows = read_csv(args.contact_edge)
    points = aligned_points(args, reference_rows, candidate_rows)
    first_high = next(
        (point for point in points if point["orders_of_magnitude_error"] > args.max_orders_threshold),
        None,
    )
    focus_bias = first_high["bias_V"] if first_high is not None else (points[0]["bias_V"] if points else None)
    summary = {
        "schema": "vela.pn2d_iv_full_range_debug.v1",
        "reference": str(args.reference),
        "candidate": str(args.candidate),
        "bias_range_V": [args.bias_min, args.bias_max],
        "points_compared": len(points),
        "max_orders_of_magnitude_error": max((p["orders_of_magnitude_error"] for p in points), default=None),
        "max_relative_error": max((p["relative_error"] for p in points), default=None),
        "first_high_error": first_high,
        "first_trend_mismatch": first_trend_mismatch(points),
        "first_non_reltol": first_non_reltol(candidate_rows, args.bias_min, args.bias_max),
        "last_converged_bias_V": last_converged_bias(candidate_rows, args.bias_min, args.bias_max),
        "terminal_current_consistency": terminal_current_consistency(terminal_rows, args.bias_min, args.bias_max),
        "top_contact_edges_by_abs_current": top_contact_edges(contact_rows, focus_bias, args.top_edges),
        "diagnostic_inputs": {
            "terminal_balance": str(args.terminal_balance) if args.terminal_balance else None,
            "contact_edge": str(args.contact_edge) if args.contact_edge else None,
            "newton_history": str(args.newton_history) if args.newton_history else None,
        },
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(summary, indent=2) + "\n")
    if args.output_md is not None:
        args.output_md.parent.mkdir(parents=True, exist_ok=True)
        args.output_md.write_text(markdown_report(summary))
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
