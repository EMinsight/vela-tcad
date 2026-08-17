#!/usr/bin/env python3
"""Close the Sentaurus 2018 SingleDevice two-branch Id-Vg contract.

The closeout is deliberately independent of the simulator.  It consumes the
same-bias comparison CSVs produced after a Vela run, applies the documented
low-current hybrid rule, and computes the device-level metrics frozen by the
reference manifest.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from pathlib import Path
from typing import Any


DEFAULT_POLICY = {
    "relative_error_max": 0.10,
    "log_error_max_dex": 0.20,
    "low_current_reference_max_A_per_um": 1.0e-13,
    "low_current_absolute_error_max_A_per_um": 1.0e-14,
    "threshold_current_A_per_um": 1.0e-7,
    "subthreshold_fit_current_min_A_per_um": 1.0e-13,
    "subthreshold_fit_current_max_A_per_um": 1.0e-8,
}


def percentile(values: list[float], fraction: float) -> float:
    if not values:
        raise ValueError("percentile requires at least one value")
    ordered = sorted(values)
    position = fraction * (len(ordered) - 1)
    left = int(math.floor(position))
    right = int(math.ceil(position))
    if left == right:
        return ordered[left]
    weight = position - left
    return ordered[left] * (1.0 - weight) + ordered[right] * weight


def load_comparison(path: Path) -> list[dict[str, float]]:
    with path.open(newline="", encoding="utf-8") as handle:
        raw = list(csv.DictReader(handle))
    required = {
        "gate_voltage_V", "sentaurus_current_A_per_um", "vela_current_A_per_um"
    }
    if not raw or not required.issubset(raw[0]):
        raise ValueError(f"{path}: comparison CSV is empty or lacks {sorted(required)}")
    rows: list[dict[str, float]] = []
    for item in raw:
        bias = float(item["gate_voltage_V"])
        expected = abs(float(item["sentaurus_current_A_per_um"]))
        actual = abs(float(item["vela_current_A_per_um"]))
        if not all(math.isfinite(value) for value in (bias, expected, actual)):
            raise ValueError(f"{path}: non-finite comparison row at {bias}")
        if expected <= 0.0 or actual <= 0.0:
            raise ValueError(f"{path}: currents must be positive magnitudes at {bias}")
        absolute = abs(actual - expected)
        relative = absolute / expected
        log_error = abs(math.log10(actual / expected))
        rows.append({
            "gate_voltage_V": bias,
            "sentaurus_current_A_per_um": expected,
            "vela_current_A_per_um": actual,
            "absolute_error_A_per_um": absolute,
            "relative_error": relative,
            "log_error_dex": log_error,
        })
    rows.sort(key=lambda row: row["gate_voltage_V"])
    return rows


def threshold_voltage(rows: list[dict[str, float]], current_key: str,
                      target: float) -> float:
    for left, right in zip(rows, rows[1:]):
        il = left[current_key]
        ir = right[current_key]
        if (il <= target <= ir) or (ir <= target <= il):
            if il == ir:
                return 0.5 * (left["gate_voltage_V"] + right["gate_voltage_V"])
            fraction = (math.log10(target) - math.log10(il)) / (
                math.log10(ir) - math.log10(il))
            return left["gate_voltage_V"] + fraction * (
                right["gate_voltage_V"] - left["gate_voltage_V"])
    raise ValueError(f"threshold current {target} is outside the curve")


def subthreshold_swing(rows: list[dict[str, float]], current_key: str,
                       low: float, high: float) -> dict[str, float | int]:
    selected = [row for row in rows if low <= row[current_key] <= high]
    if len(selected) < 2:
        raise ValueError(
            f"need at least two points in the subthreshold window [{low}, {high}]")
    xs = [row["gate_voltage_V"] for row in selected]
    ys = [math.log10(row[current_key]) for row in selected]
    xmean = statistics.fmean(xs)
    ymean = statistics.fmean(ys)
    denominator = sum((x - xmean) ** 2 for x in xs)
    if denominator == 0.0:
        raise ValueError("degenerate subthreshold fit")
    slope = sum((x - xmean) * (y - ymean) for x, y in zip(xs, ys)) / denominator
    if slope <= 0.0:
        raise ValueError("subthreshold current is not increasing with gate voltage")
    return {"mV_per_dec": 1000.0 / slope, "fit_points": len(selected)}


def point_gate(row: dict[str, float], policy: dict[str, float]) -> dict[str, Any]:
    relative_pass = row["relative_error"] <= policy["relative_error_max"]
    low_current = (
        row["sentaurus_current_A_per_um"]
        < policy["low_current_reference_max_A_per_um"]
    )
    fallback_pass = (
        low_current
        and row["log_error_dex"] <= policy["log_error_max_dex"]
        and row["absolute_error_A_per_um"]
        <= policy["low_current_absolute_error_max_A_per_um"]
    )
    return {
        "gate_voltage_V": row["gate_voltage_V"],
        "status": "pass" if relative_pass or fallback_pass else "fail",
        "acceptance_path": (
            "relative" if relative_pass else
            "low_current_hybrid" if fallback_pass else "none"
        ),
        "low_current": low_current,
        "relative_error": row["relative_error"],
        "log_error_dex": row["log_error_dex"],
        "absolute_error_A_per_um": row["absolute_error_A_per_um"],
    }


def branch_metrics(rows: list[dict[str, float]], policy: dict[str, float]) -> dict[str, Any]:
    point_results = [point_gate(row, policy) for row in rows]
    log_errors = [row["log_error_dex"] for row in rows]
    reference_vth = threshold_voltage(
        rows, "sentaurus_current_A_per_um", policy["threshold_current_A_per_um"])
    vela_vth = threshold_voltage(
        rows, "vela_current_A_per_um", policy["threshold_current_A_per_um"])
    reference_ss = subthreshold_swing(
        rows, "sentaurus_current_A_per_um",
        policy["subthreshold_fit_current_min_A_per_um"],
        policy["subthreshold_fit_current_max_A_per_um"])
    vela_ss = subthreshold_swing(
        rows, "vela_current_A_per_um",
        policy["subthreshold_fit_current_min_A_per_um"],
        policy["subthreshold_fit_current_max_A_per_um"])
    ss_relative_error = abs(
        float(vela_ss["mV_per_dec"]) - float(reference_ss["mV_per_dec"])
    ) / float(reference_ss["mV_per_dec"])
    return {
        "status": "pass" if all(row["status"] == "pass" for row in point_results) else "fail",
        "points": len(rows),
        "direct_bias_match": len({row["gate_voltage_V"] for row in rows}) == len(rows),
        "trend_match": all(
            rows[index]["vela_current_A_per_um"]
            <= rows[index + 1]["vela_current_A_per_um"]
            for index in range(len(rows) - 1)
        ),
        "median_abs_log_error_dex": statistics.median(log_errors),
        "p95_abs_log_error_dex": percentile(log_errors, 0.95),
        "max_abs_log_error_dex": max(log_errors),
        "max_relative_error": max(row["relative_error"] for row in rows),
        "max_absolute_error_A_per_um": max(
            row["absolute_error_A_per_um"] for row in rows),
        "ion": {
            "gate_voltage_V": rows[-1]["gate_voltage_V"],
            "sentaurus_A_per_um": rows[-1]["sentaurus_current_A_per_um"],
            "vela_A_per_um": rows[-1]["vela_current_A_per_um"],
            "relative_error": rows[-1]["relative_error"],
        },
        "ioff": {
            "gate_voltage_V": rows[0]["gate_voltage_V"],
            "sentaurus_A_per_um": rows[0]["sentaurus_current_A_per_um"],
            "vela_A_per_um": rows[0]["vela_current_A_per_um"],
            "relative_error": rows[0]["relative_error"],
            "log_error_dex": rows[0]["log_error_dex"],
            "absolute_error_A_per_um": rows[0]["absolute_error_A_per_um"],
            "acceptance_path": point_results[0]["acceptance_path"],
        },
        "threshold": {
            "target_current_A_per_um": policy["threshold_current_A_per_um"],
            "sentaurus_V": reference_vth,
            "vela_V": vela_vth,
            "delta_V": vela_vth - reference_vth,
        },
        "subthreshold_swing": {
            "fit_current_min_A_per_um": policy["subthreshold_fit_current_min_A_per_um"],
            "fit_current_max_A_per_um": policy["subthreshold_fit_current_max_A_per_um"],
            "sentaurus_mV_per_dec": reference_ss["mV_per_dec"],
            "vela_mV_per_dec": vela_ss["mV_per_dec"],
            "sentaurus_fit_points": reference_ss["fit_points"],
            "vela_fit_points": vela_ss["fit_points"],
            "relative_error": ss_relative_error,
        },
        "point_acceptance": point_results,
    }


def closeout(linear: list[dict[str, float]], saturation: list[dict[str, float]],
             reference: dict[str, Any], policy: dict[str, float]) -> dict[str, Any]:
    lin = branch_metrics(linear, policy)
    sat = branch_metrics(saturation, policy)
    drain_delta = 1.1 - 0.1
    sentaurus_dibl = (
        lin["threshold"]["sentaurus_V"] - sat["threshold"]["sentaurus_V"]
    ) / drain_delta
    vela_dibl = (
        lin["threshold"]["vela_V"] - sat["threshold"]["vela_V"]
    ) / drain_delta
    limits = reference["acceptance"]
    derived_checks = {
        "linear_branch": lin["status"] == "pass",
        "saturation_branch": sat["status"] == "pass",
        "minimum_points": lin["points"] >= 21 and sat["points"] >= 21,
        "trend": lin["trend_match"] and sat["trend_match"],
        "median_log_error": max(
            lin["median_abs_log_error_dex"], sat["median_abs_log_error_dex"]
        ) <= float(limits["median_log10_current_error_max_dex"]),
        "p95_log_error": max(
            lin["p95_abs_log_error_dex"], sat["p95_abs_log_error_dex"]
        ) <= float(limits["p95_log10_current_error_max_dex"]),
        "ion": max(
            lin["ion"]["relative_error"], sat["ion"]["relative_error"]
        ) <= float(limits["ion_relative_error_max"]),
        "threshold": max(
            abs(lin["threshold"]["delta_V"]), abs(sat["threshold"]["delta_V"])
        ) <= float(limits["threshold_voltage_error_max_V"]),
        "dibl": abs(vela_dibl - sentaurus_dibl)
        <= float(limits["dibl_error_max_V_per_V"]),
        "subthreshold_swing": max(
            lin["subthreshold_swing"]["relative_error"],
            sat["subthreshold_swing"]["relative_error"],
        ) <= float(limits["subthreshold_swing_relative_error_max"]),
    }
    return {
        "schema": "vela.singledevice.closeout.v1",
        "status": "pass" if all(derived_checks.values()) else "fail",
        "case": reference.get("case", "singledevice_sentaurus2018"),
        "sentaurus_version": reference.get("sentaurus_version"),
        "acceptance_policy": policy,
        "branches": {"linear": lin, "saturation": sat},
        "dibl": {
            "drain_voltage_delta_V": drain_delta,
            "sentaurus_V_per_V": sentaurus_dibl,
            "vela_V_per_V": vela_dibl,
            "absolute_error_V_per_V": abs(vela_dibl - sentaurus_dibl),
        },
        "checks": derived_checks,
    }


def markdown(report: dict[str, Any]) -> str:
    lin = report["branches"]["linear"]
    sat = report["branches"]["saturation"]
    return "\n".join([
        "# SingleDevice closeout summary",
        "",
        f"Overall status: **{report['status']}**.",
        "",
        "| Metric | Linear | Saturation |",
        "|---|---:|---:|",
        f"| Points | {lin['points']} | {sat['points']} |",
        f"| Median abs log error (dex) | {lin['median_abs_log_error_dex']:.6g} | {sat['median_abs_log_error_dex']:.6g} |",
        f"| P95 abs log error (dex) | {lin['p95_abs_log_error_dex']:.6g} | {sat['p95_abs_log_error_dex']:.6g} |",
        f"| Ion relative error | {lin['ion']['relative_error']:.6%} | {sat['ion']['relative_error']:.6%} |",
        f"| Vth delta (mV) | {1000.0 * lin['threshold']['delta_V']:.6g} | {1000.0 * sat['threshold']['delta_V']:.6g} |",
        f"| SS relative error | {lin['subthreshold_swing']['relative_error']:.6%} | {sat['subthreshold_swing']['relative_error']:.6%} |",
        "",
        f"DIBL absolute error: `{report['dibl']['absolute_error_V_per_V']:.6g} V/V`.",
        "",
        "The linear deep-off point is accepted only through the documented "
        f"`{lin['ioff']['acceptance_path']}` path; raw relative, log, and absolute "
        f"errors remain recorded as `{lin['ioff']['relative_error']:.6%}`, "
        f"`{lin['ioff']['log_error_dex']:.6g} dex`, and "
        f"`{lin['ioff']['absolute_error_A_per_um']:.6g} A/um`.",
        "",
    ])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--linear", type=Path, required=True)
    parser.add_argument("--saturation", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-md", type=Path)
    args = parser.parse_args()

    reference = json.loads(args.reference.read_text(encoding="utf-8"))
    policy = dict(DEFAULT_POLICY)
    policy.update(reference.get("low_current_acceptance", {}))
    report = closeout(
        load_comparison(args.linear), load_comparison(args.saturation),
        reference, policy)
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    if args.output_md:
        args.output_md.parent.mkdir(parents=True, exist_ok=True)
        args.output_md.write_text(markdown(report), encoding="utf-8")
    print(json.dumps({"status": report["status"], "checks": report["checks"]}))
    return 0 if report["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
