#!/usr/bin/env python3
"""Summarize the BVmethods NMOS 0--7 V E2 and IIC follow-up runs."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any


Q_C = 1.602176634e-19
REPO = Path(__file__).resolve().parents[1]
RUN_ROOT = REPO / "build-release/reference_tcad/bvmethods_sentaurus2018/run01"
DEFAULT_BRANCH = (
    RUN_ROOT / "vela_validation/btbt_e2_adaptive_0_7_20260804/branch_0_7.csv"
)
DEFAULT_NO_E2 = (
    RUN_ROOT
    / "vela_validation/iic_rebuild_20260803/trunk_1p9_2p0_earlyfloor/"
      "postprocess_only/sweep.csv"
)
DEFAULT_SENT = RUN_ROOT / "analysis/curves/ABA_coupled.csv"
DEFAULT_SENT_EXACT = (
    RUN_ROOT
    / "vela_validation/iic_postprocess_20260803/analysis/multibias_sentaurus/"
      "sentaurus_exact_extended_curve.csv"
)
DEFAULT_IIC = (
    RUN_ROOT / "vela_validation/btbt_e2_iic_reclosure_20260804/postprocess_only"
)
DEFAULT_OUTPUT = (
    RUN_ROOT / "vela_validation/btbt_e2_followup_20260804/analysis"
)


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def number(row: dict[str, str], key: str, default: float = math.nan) -> float:
    value = row.get(key, "")
    return float(value) if value not in (None, "") else default


def write_csv(path: Path, records: list[dict[str, Any]]) -> None:
    if not records:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def interpolate(
    points: list[dict[str, str]], bias: float, key: str, x_key: str = "inner_voltage_V"
) -> float:
    ordered = sorted(points, key=lambda row: number(row, x_key))
    for row in ordered:
        if math.isclose(number(row, x_key), bias, rel_tol=0.0, abs_tol=1.0e-12):
            return number(row, key)
    for left, right in zip(ordered, ordered[1:]):
        x0, x1 = number(left, x_key), number(right, x_key)
        if x0 <= bias <= x1:
            y0, y1 = number(left, key), number(right, key)
            weight = (bias - x0) / (x1 - x0)
            if y0 * y1 > 0.0:
                sign = math.copysign(1.0, y0)
                return sign * 10.0 ** (
                    math.log10(abs(y0))
                    + weight * (math.log10(abs(y1)) - math.log10(abs(y0)))
                )
            return y0 + weight * (y1 - y0)
    return math.nan


def reference_value(
    curve: list[dict[str, str]], exact: list[dict[str, str]], bias: float, key: str
) -> float:
    matches = [
        row for row in exact
        if math.isclose(number(row, "inner_voltage_V"), bias, rel_tol=0.0, abs_tol=1.0e-12)
    ]
    if matches:
        # The extended file contains two equilibrium states.  Positive-bias
        # checkpoints are unique; for 0 V use the final imported equilibrium.
        return number(matches[-1], key)
    return interpolate(curve, bias, key)


def source_to_current(source_integral: float) -> float:
    """Convert alpha[cm^-1]*flux[cm^-2/s]*area[um^2] to A/um."""
    return Q_C * source_integral * 1.0e-12


def row_at(rows: list[dict[str, str]], bias: float, key: str = "bias_V") -> dict[str, str]:
    return min(rows, key=lambda row: abs(number(row, key) - bias))


def branch_compare(
    branch: list[dict[str, str]], sentaurus: list[dict[str, str]], exact: list[dict[str, str]]
) -> list[dict[str, Any]]:
    records: list[dict[str, Any]] = []
    for row in branch:
        bias = number(row, "bias_V")
        vela = number(row, "current_total_A_per_um")
        sent = reference_value(
            sentaurus, exact, bias, "drain_total_current_A_per_um"
        )
        records.append({
            "bias_V": bias,
            "vela_Id_A_per_um": vela,
            "sentaurus_Id_A_per_um_log_interp": sent,
            "signed_relative_error": (vela - sent) / sent if sent else math.nan,
            "abs_vela_over_sentaurus": abs(vela / sent) if sent else math.nan,
            "max_electric_field_V_per_m": number(row, "max_electric_field_V_per_m"),
            "newton_iterations": int(number(row, "iterations", 0.0)),
            "carrier_row_max_ratio": number(row, "carrier_row_max_ratio"),
            "convergence_reason": row.get("newton_convergence_reason", ""),
        })
    return records


def low_bias_decomposition(
    branch: list[dict[str, str]], no_e2: list[dict[str, str]],
    sentaurus: list[dict[str, str]], exact: list[dict[str, str]]
) -> dict[str, Any]:
    e2_row = row_at(branch, 2.0)
    no_e2_row = row_at(no_e2, 2.0)
    e2 = number(e2_row, "current_total_A_per_um")
    baseline = number(no_e2_row, "current_total_A_per_um")
    sent = reference_value(
        sentaurus, exact, 2.0, "drain_total_current_A_per_um"
    )
    increment = e2 - baseline
    remaining = sent - e2
    field_e2 = number(e2_row, "max_electric_field_V_per_m")
    field_baseline = number(no_e2_row, "max_electric_field_V_per_m")
    return {
        "bias_V": 2.0,
        "vela_no_E2_Id_A_per_um": baseline,
        "vela_E2_Id_A_per_um": e2,
        "sentaurus_Id_A_per_um": sent,
        "E2_current_increment_A_per_um": increment,
        "E2_increment_over_no_E2": increment / baseline,
        "E2_fraction_of_full_current": increment / e2,
        "remaining_sentaurus_minus_vela_A_per_um": remaining,
        "remaining_relative_to_sentaurus": remaining / sent,
        "no_E2_max_electric_field_V_per_m": field_baseline,
        "E2_max_electric_field_V_per_m": field_e2,
        "electric_field_relative_change": (field_e2 - field_baseline) / field_baseline,
        "interpretation": (
            "The 2 V discrepancy is dominated by the contact/SG transport baseline, "
            "not by missing E2 generation."
        ),
    }


def iic_compare(
    iic_root: Path, sentaurus: list[dict[str, str]], exact: list[dict[str, str]]
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    sweep = read_rows(iic_root / "sweep.csv")
    avalanche = read_rows(iic_root / "avalanche_summary.csv")
    records: list[dict[str, Any]] = []
    for source_row in avalanche:
        bias = number(source_row, "bias_V")
        sweep_row = row_at(sweep, bias)
        vela_id = number(sweep_row, "current_total_A_per_um")
        vela_iava = source_to_current(number(source_row, "sum_edge_source_integral"))
        sent_id = reference_value(
            sentaurus, exact, bias, "drain_total_current_A_per_um"
        )
        sent_iava = reference_value(
            sentaurus, exact, bias, "avalanche_current_A_per_um"
        )
        records.append({
            "bias_V": bias,
            "vela_Id_A_per_um": vela_id,
            "vela_Iava_A_per_um": vela_iava,
            "vela_Iava_over_abs_Id": vela_iava / abs(vela_id),
            "sentaurus_Id_A_per_um_log_interp": sent_id,
            "sentaurus_Iava_A_per_um_log_interp": sent_iava,
            "sentaurus_Iava_over_abs_Id": abs(sent_iava / sent_id),
            "abs_vela_over_sentaurus_Id": abs(vela_id / sent_id),
            "abs_vela_over_sentaurus_Iava": abs(vela_iava / sent_iava),
            "max_electron_alpha_m_inv": number(source_row, "max_electron_alpha_m_inv"),
            "max_hole_alpha_m_inv": number(source_row, "max_hole_alpha_m_inv"),
        })
    records.sort(key=lambda row: row["bias_V"])
    def current_crossing(
        iava_key: str, id_key: str
    ) -> tuple[tuple[dict[str, Any], dict[str, Any]] | None, float]:
        for left, right in zip(records, records[1:]):
            y0 = left[iava_key] - abs(left[id_key])
            y1 = right[iava_key] - abs(right[id_key])
            if y0 * y1 <= 0.0:
                x0, x1 = left["bias_V"], right["bias_V"]
                return (left, right), x0 - y0 * (x1 - x0) / (y1 - y0)
        return None, math.nan

    bracket, estimate = current_crossing("vela_Iava_A_per_um", "vela_Id_A_per_um")
    sent_bracket, sent_estimate = current_crossing(
        "sentaurus_Iava_A_per_um_log_interp",
        "sentaurus_Id_A_per_um_log_interp",
    )
    near = min(records, key=lambda row: abs(row["bias_V"] - 6.377494277837012))
    summary = {
        "sentaurus_iic_BV_V": 6.377494277837012,
        "vela_iic_crossing_bracketed": bracket is not None,
        "vela_iic_crossing_estimate_V": estimate,
        "vela_iic_crossing_bracket_V": (
            [bracket[0]["bias_V"], bracket[1]["bias_V"]] if bracket else []
        ),
        "vela_max_tested_Iava_over_abs_Id": max(
            row["vela_Iava_over_abs_Id"] for row in records
        ),
        "vela_max_tested_bias_V": max(row["bias_V"] for row in records),
        "sentaurus_dense_current_crossing_estimate_V": sent_estimate,
        "sentaurus_dense_current_crossing_bracket_V": (
            [sent_bracket[0]["bias_V"], sent_bracket[1]["bias_V"]]
            if sent_bracket else []
        ),
        "near_sentaurus_iic": near,
    }
    return records, summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--branch", type=Path, default=DEFAULT_BRANCH)
    parser.add_argument("--no-e2", type=Path, default=DEFAULT_NO_E2)
    parser.add_argument("--sentaurus", type=Path, default=DEFAULT_SENT)
    parser.add_argument("--sentaurus-exact", type=Path, default=DEFAULT_SENT_EXACT)
    parser.add_argument("--iic-root", type=Path, default=DEFAULT_IIC)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    branch = read_rows(args.branch)
    no_e2 = read_rows(args.no_e2)
    sentaurus = read_rows(args.sentaurus)
    exact = read_rows(args.sentaurus_exact)
    comparison = branch_compare(branch, sentaurus, exact)
    low_bias = low_bias_decomposition(branch, no_e2, sentaurus, exact)
    write_csv(args.out_dir / "e2_branch_compare.csv", comparison)
    write_csv(args.out_dir / "e2_2v_baseline.csv", [low_bias])

    summary: dict[str, Any] = {
        "branch_bias_range_V": [number(branch[0], "bias_V"), number(branch[-1], "bias_V")],
        "branch_point_count": len(branch),
        "two_volt_decomposition": low_bias,
        "iic_available": False,
    }
    if (args.iic_root / "avalanche_summary.csv").exists():
        iic_records, iic_summary = iic_compare(args.iic_root, sentaurus, exact)
        write_csv(args.out_dir / "e2_iic_compare.csv", iic_records)
        summary["iic_available"] = True
        summary["iic"] = iic_summary

    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "e2_followup_summary.json").write_text(
        json.dumps(summary, indent=2, allow_nan=True) + "\n", encoding="utf-8"
    )
    print(args.out_dir / "e2_followup_summary.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
