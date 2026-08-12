#!/usr/bin/env python3
"""Summarize the Vela A/B/C/D full-physics BV ablation outputs."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


VARIANTS = (
    "a_constant_no_enormal",
    "b_doping_no_enormal",
    "c_constant_enormal",
    "d_doping_enormal",
)
SENTAURUS_VOLTAGE_TO_CURRENT_V = 6.38318420057198
SENTAURUS_EXTERNAL_RESISTOR_V = 6.379791636301563


def read_target_row(path: Path, target: float) -> dict[str, str]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream))
    matches = [
        row
        for row in rows
        if row.get("boundary_control_mode") == "current"
        and abs(float(row["target_current_A_per_um"]) - target) <= 1.0e-12
    ]
    if len(matches) != 1:
        raise ValueError(f"expected one {target:g} A/um row in {path}, got {len(matches)}")
    if matches[0].get("converged") != "1":
        raise ValueError(f"target row did not converge in {path}")
    return matches[0]


def read_external_resistor_row(path: Path) -> dict[str, str]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        rows = list(csv.DictReader(stream))
    matches = [
        row
        for row in rows
        if row.get("boundary_control_mode") == "external_resistor"
        and row.get("converged") == "1"
    ]
    if len(matches) != 1:
        raise ValueError(f"expected one converged external-resistor row in {path}")
    return matches[0]


def summarize(
    root: Path,
    target: float = 1.0e-4,
    external_resistor_csv: Path | None = None,
) -> dict:
    variants = {}
    for name in VARIANTS:
        row = read_target_row(root / name / "sweep.csv", target)
        variants[name] = {
            "voltage_V": float(row["inner_voltage_V"]),
            "current_A_per_um": float(row["current_total_A_per_um"]),
            "current_residual_A_per_um": float(
                row["current_boundary_residual_A_per_um"]
            ),
            "boundary_control_evaluations": int(row["boundary_control_evaluations"]),
            "newton_iterations": int(row["newton_iterations"]),
            "global_electron_continuity_closure_ratio": float(
                row["global_electron_continuity_closure_ratio"]
            ),
            "global_hole_continuity_closure_ratio": float(
                row["global_hole_continuity_closure_ratio"]
            ),
        }
    a = variants[VARIANTS[0]]["voltage_V"]
    b = variants[VARIANTS[1]]["voltage_V"]
    c = variants[VARIANTS[2]]["voltage_V"]
    d = variants[VARIANTS[3]]["voltage_V"]
    result = {
        "schema": "vela.bvmethods_full_physics_ablation_summary.v1",
        "status": "PASS",
        "target_current_A_per_um": target,
        "variants": variants,
        "voltage_increments_V": {
            "srh_doping_dependence_B_minus_A": b - a,
            "enormal_C_minus_A": c - a,
            "combined_D_minus_A": d - a,
            "srh_with_enormal_D_minus_C": d - c,
            "enormal_with_doping_D_minus_B": d - b,
            "interaction_D_minus_B_minus_C_plus_A": d - b - c + a,
        },
        "sentaurus_full_model_acceptance": {
            "method": "voltage_to_current",
            "sentaurus_voltage_V": SENTAURUS_VOLTAGE_TO_CURRENT_V,
            "vela_voltage_V": d,
            "absolute_error_V": abs(d - SENTAURUS_VOLTAGE_TO_CURRENT_V),
            "relative_error": abs(d - SENTAURUS_VOLTAGE_TO_CURRENT_V)
            / SENTAURUS_VOLTAGE_TO_CURRENT_V,
            "relative_error_limit": 0.02,
            "status": (
                "PASS"
                if abs(d - SENTAURUS_VOLTAGE_TO_CURRENT_V)
                / SENTAURUS_VOLTAGE_TO_CURRENT_V
                <= 0.02
                else "FAIL"
            ),
        },
    }
    if external_resistor_csv is not None:
        row = read_external_resistor_row(external_resistor_csv)
        voltage = float(row["inner_voltage_V"])
        current = float(row["current_total_A_per_um"])
        result["external_resistor_cross_check"] = {
            "sentaurus_voltage_V": SENTAURUS_EXTERNAL_RESISTOR_V,
            "vela_voltage_V": voltage,
            "vela_current_A_per_um": current,
            "load_line_residual_V": float(row["load_line_residual_V"]),
            "absolute_error_vs_sentaurus_V": abs(
                voltage - SENTAURUS_EXTERNAL_RESISTOR_V
            ),
            "relative_error_vs_sentaurus": abs(
                voltage - SENTAURUS_EXTERNAL_RESISTOR_V
            )
            / SENTAURUS_EXTERNAL_RESISTOR_V,
            "absolute_difference_vs_voltage_to_current_V": abs(voltage - d),
            "relative_difference_vs_voltage_to_current": abs(voltage - d) / d,
            "status": "PASS"
            if abs(float(row["load_line_residual_V"])) <= 0.1
            and abs(voltage - SENTAURUS_EXTERNAL_RESISTOR_V)
            / SENTAURUS_EXTERNAL_RESISTOR_V
            <= 0.02
            else "FAIL",
        }
    result["status"] = (
        "PASS"
        if result["sentaurus_full_model_acceptance"]["status"] == "PASS"
        and result.get("external_resistor_cross_check", {}).get("status", "PASS")
        == "PASS"
        else "FAIL"
    )
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--target-current", type=float, default=1.0e-4)
    parser.add_argument("--external-resistor-csv", type=Path)
    args = parser.parse_args()
    result = summarize(
        args.root.resolve(),
        args.target_current,
        args.external_resistor_csv.resolve()
        if args.external_resistor_csv is not None
        else None,
    )
    output = args.output or args.root / "ablation_summary.json"
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(output.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
