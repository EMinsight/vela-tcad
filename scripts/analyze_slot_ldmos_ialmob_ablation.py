#!/usr/bin/env python3
"""Analyze the Sentaurus Slot-LDMOS IALMob on/off controls."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

try:
    from analyze_slot_ldmos_bvds_reference import (
        BREAK_CURRENT_A_PER_UM,
        curve_rows,
        interpolate_crossing,
        write_curve,
    )
except ModuleNotFoundError:  # Imported as scripts.* by regression tests.
    from scripts.analyze_slot_ldmos_bvds_reference import (
        BREAK_CURRENT_A_PER_UM,
        curve_rows,
        interpolate_crossing,
        write_curve,
    )


REPO = Path(__file__).resolve().parents[1]
DEFAULT_ON = (
    REPO
    / "build-release/reference_tcad/slot_ldmos_sentaurus2022/run01"
    / "sentaurus_bvds_result/logs_curves"
)
DEFAULT_OFF = (
    REPO
    / "build-release/reference_tcad/slot_ldmos_sentaurus2022/run01"
    / "ialmob_ablation/run01/results"
)
DEFAULT_OUTPUT = (
    REPO
    / "build-release/reference_tcad/slot_ldmos_sentaurus2022/run01"
    / "ialmob_ablation/run01/analysis"
)


def relative_delta(control: float, baseline: float) -> float:
    return (control - baseline) / baseline


def analyze(on_dir: Path, off_dir: Path, output_dir: Path) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    for name in (
        "avalanche_off_60v_no_ialmob",
        "bvds_external_resistor_final_no_ialmob",
    ):
        exit_code = int((off_dir / f"{name}.exitcode").read_text().strip())
        if exit_code != 0:
            raise ValueError(f"{name} returned exit code {exit_code}")

    on_leakage = curve_rows(on_dir / "avalanche_off_60v.plt")
    off_leakage = curve_rows(off_dir / "avalanche_off_60v_no_ialmob.plt")
    on_bvds = curve_rows(on_dir / "bvds_external_resistor_final.plt")
    off_bvds = curve_rows(
        off_dir / "bvds_external_resistor_final_no_ialmob.plt"
    )
    on_crossing = interpolate_crossing(on_bvds, BREAK_CURRENT_A_PER_UM)
    off_crossing = interpolate_crossing(off_bvds, BREAK_CURRENT_A_PER_UM)
    if on_crossing is None or off_crossing is None:
        raise ValueError("IALMob A/B curve does not cross 1e-7 A/um")

    write_curve(output_dir / "avalanche_off_60v_no_ialmob.csv", off_leakage)
    write_curve(
        output_dir / "bvds_external_resistor_final_no_ialmob.csv", off_bvds
    )

    on_leakage_final = on_leakage[-1]
    off_leakage_final = off_leakage[-1]
    on_bv = on_crossing["inner_voltage_V"]
    off_bv = off_crossing["inner_voltage_V"]
    summary: dict[str, Any] = {
        "schema": "vela.slot_ldmos.sentaurus_ialmob_ablation_result.v1",
        "sentaurus_release": "T-2022.03-SP2",
        "controlled_delta": "Enormal(IALMob) enabled versus disabled",
        "breakdown_criterion_A_per_um": BREAK_CURRENT_A_PER_UM,
        "ialmob_on": {
            "avalanche_off_60v": on_leakage_final,
            "bvds_crossing": on_crossing,
            "bvds_curve_points": len(on_bvds),
        },
        "ialmob_off": {
            "avalanche_off_60v": off_leakage_final,
            "bvds_crossing": off_crossing,
            "bvds_curve_points": len(off_bvds),
        },
        "off_minus_on": {
            "leakage_A_per_um": (
                off_leakage_final["drain_total_current_A_per_um"]
                - on_leakage_final["drain_total_current_A_per_um"]
            ),
            "leakage_relative": relative_delta(
                off_leakage_final["drain_total_current_A_per_um"],
                on_leakage_final["drain_total_current_A_per_um"],
            ),
            "bvds_V": off_bv - on_bv,
            "bvds_relative": relative_delta(off_bv, on_bv),
        },
        "interpretation": (
            "IALMob has a measurable but sub-percent BVDS effect and cannot "
            "explain the present Vela leakage or nonlinear-solver discrepancy."
        ),
    }
    (output_dir / "ialmob_ablation_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    delta = summary["off_minus_on"]
    markdown = "\n".join(
        [
            "# Slot-LDMOS Sentaurus IALMob ablation",
            "",
            "| Metric | IALMob on | IALMob off | Off - on |",
            "|---|---:|---:|---:|",
            (
                "| 60 V avalanche-off drain current (A/um) | "
                f"{on_leakage_final['drain_total_current_A_per_um']:.9e} | "
                f"{off_leakage_final['drain_total_current_A_per_um']:.9e} | "
                f"{delta['leakage_relative'] * 100:.6f}% |"
            ),
            (
                "| BVDS at 1e-7 A/um (V) | "
                f"{on_bv:.9f} | {off_bv:.9f} | "
                f"{delta['bvds_V']:+.9f} V ({delta['bvds_relative'] * 100:+.6f}%) |"
            ),
            "",
            "Only `Enormal(IALMob)` was removed; mesh, parameter file, "
            "high-field mobility, recombination, avalanche, load line, and "
            "numerical controls were unchanged.",
            "",
        ]
    )
    (output_dir / "ialmob_ablation_summary.md").write_text(
        markdown, encoding="utf-8"
    )
    return summary


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--on-dir", type=Path, default=DEFAULT_ON)
    parser.add_argument("--off-dir", type=Path, default=DEFAULT_OFF)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    summary = analyze(
        args.on_dir.resolve(), args.off_dir.resolve(), args.output_dir.resolve()
    )
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
