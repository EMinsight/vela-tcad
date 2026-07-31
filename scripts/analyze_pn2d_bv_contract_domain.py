#!/usr/bin/env python3
"""Analyze PN2D BV parity only on a contract's predeclared effective domain."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from analyze_pn2d_avalanche_on_bv_parity import (
    CURVE_NAMES,
    KNEE_BIASES_V,
    adjacent_slopes,
    comparison_rows,
    continuous_breakpoint,
    curvature_knee,
    curve_error_metrics,
    exact_curve_index,
    gain_error_metrics,
    load_curve,
    sha256,
    slope_knee,
    validate_curve_points,
)


def analyze(raw_curves, effective_biases):
    required = tuple(sorted(set(effective_biases) | set(KNEE_BIASES_V), reverse=True))
    indexed = {}
    for name, points in raw_curves.items():
        validate_curve_points(points, name)
        rows, missing = exact_curve_index(points, required, curve_name=name)
        if missing:
            raise ValueError(f"{name} is missing exact contract rows: {missing}")
        indexed[name] = rows
    effective_rows, effective_metrics = curve_error_metrics(
        indexed["vela_on"], indexed["sentaurus_on"], effective_biases
    )
    knee_rows, knee_metrics = curve_error_metrics(
        indexed["vela_on"], indexed["sentaurus_on"], KNEE_BIASES_V
    )
    vela_slopes = adjacent_slopes(
        indexed["vela_on"], KNEE_BIASES_V, "vela_on"
    )
    sentaurus_slopes = adjacent_slopes(
        indexed["sentaurus_on"], KNEE_BIASES_V, "sentaurus_on"
    )
    slope_errors = [
        float(vela["slope_dex_per_V"]) - float(sentaurus["slope_dex_per_V"])
        for vela, sentaurus in zip(vela_slopes, sentaurus_slopes, strict=True)
    ]
    slope_rmse = math.sqrt(
        sum(error * error for error in slope_errors) / len(slope_errors)
    )
    knee_estimators = {}
    for name, slopes in (
        ("vela", vela_slopes),
        ("sentaurus", sentaurus_slopes),
    ):
        curve = indexed[f"{name}_on"]
        knee_estimators[name] = {
            "V_slope": slope_knee(slopes),
            "V_break": continuous_breakpoint(
                curve, KNEE_BIASES_V, f"{name}_on"
            ),
            "V_curvature": curvature_knee(slopes),
        }
    return {
        "schema": "vela.pn2d_bv_contract_domain_parity.v1",
        "outcome": "contract_domain_metrics_complete",
        "effective_biases_V": list(effective_biases),
        "effective_metrics": effective_metrics,
        "knee_metrics": knee_metrics,
        "gain_metrics": gain_error_metrics(indexed, effective_biases),
        "knee_estimators": knee_estimators,
        "adjacent_slope_rmse_dex_per_V": slope_rmse,
        "curve_rows": comparison_rows(indexed, required),
        "effective_error_rows": effective_rows,
        "knee_error_rows": knee_rows,
        "vela_slopes": vela_slopes,
        "sentaurus_slopes": sentaurus_slopes,
        "zero_current_outside_contract_domain_is_not_evaluated": True,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, required=True)
    for name in CURVE_NAMES:
        parser.add_argument(
            f"--{name.replace('_', '-')}-csv", type=Path, required=True
        )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    contract = json.loads(args.contract.read_text(encoding="utf-8"))
    effective_biases = tuple(
        float(value) for value in contract["bv_domain"]["exact_biases_V"]
    )
    curves = {}
    inputs = {}
    for name in CURVE_NAMES:
        path = getattr(args, f"{name}_csv").resolve()
        curves[name] = load_curve(path)
        inputs[name] = {"path": str(path), "sha256": sha256(path)}
    result = analyze(curves, effective_biases)
    result["inputs"] = inputs
    result["contract"] = {
        "path": str(args.contract.resolve()),
        "sha256": sha256(args.contract),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
