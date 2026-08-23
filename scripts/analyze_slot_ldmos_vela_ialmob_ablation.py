#!/usr/bin/env python3
"""Analyze Vela Slot-LDMOS Enhanced-Lombardi IALMob A/B outputs."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any


CASES = ("ialmob_off", "ialmob_on")


class IALMobAnalysisError(ValueError):
    """Raised when a completed A/B output is missing or inconsistent."""


def last_csv_row(path: Path) -> dict[str, str]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise IALMobAnalysisError(f"empty CSV: {path}")
    return rows[-1]


def read_vtk_node_scalar_peak(
    path: Path, scalar_name: str
) -> dict[str, float | int]:
    points: list[tuple[float, float, float]] = []
    point_count: int | None = None
    values: list[float] | None = None
    with path.open(encoding="utf-8") as handle:
        iterator = iter(handle)
        for line in iterator:
            words = line.split()
            if len(words) >= 3 and words[0] == "POINTS":
                point_count = int(words[1])
                while len(points) < point_count:
                    coordinates = [float(value) for value in next(iterator).split()]
                    for index in range(0, len(coordinates), 3):
                        points.append(tuple(coordinates[index : index + 3]))
                points = points[:point_count]
            if (
                len(words) >= 3
                and words[0] == "SCALARS"
                and words[1] == scalar_name
            ):
                if point_count is None:
                    raise IALMobAnalysisError("VTK scalar appears before POINTS")
                lookup = next(iterator).strip()
                if not lookup.startswith("LOOKUP_TABLE"):
                    raise IALMobAnalysisError(
                        f"missing LOOKUP_TABLE for {scalar_name}"
                    )
                values = []
                while len(values) < point_count:
                    values.extend(float(value) for value in next(iterator).split())
                values = values[:point_count]
                break
    if point_count is None or len(points) != point_count:
        raise IALMobAnalysisError(f"invalid VTK points in {path}")
    if values is None:
        raise IALMobAnalysisError(f"missing scalar {scalar_name} in {path}")
    index = max(range(point_count), key=values.__getitem__)
    x, y, z = points[index]
    return {
        "node_index": index,
        "value": values[index],
        "x_um": x,
        "y_um": y,
        "z_um": z,
    }


def case_metrics(case_dir: Path) -> dict[str, Any]:
    iv = last_csv_row(case_dir / "iv.csv")
    avalanche = last_csv_row(case_dir / "avalanche_summary.csv")
    vtk_files = sorted((case_dir / "vtk").glob("state_*.vtk"))
    if len(vtk_files) != 1:
        raise IALMobAnalysisError(
            f"expected one VTK in {case_dir / 'vtk'}, found {len(vtk_files)}"
        )
    return {
        "outer_voltage_V": float(iv["outer_voltage_V"]),
        "inner_voltage_V": float(iv["inner_voltage_V"]),
        "drain_current_A_per_um": float(iv["current_total_A_per_um"]),
        "load_line_residual_V": float(iv["load_line_residual_V"]),
        "qG_full_A_per_um": float(avalanche["qG_full"]),
        "max_electric_field_V_per_cm": float(avalanche["max_E"]),
        "max_avalanche_generation_cm3_s": float(avalanche["max_Gava"]),
        "avalanche_peak": read_vtk_node_scalar_peak(
            vtk_files[0], "AvalancheGeneration"
        ),
        "vtk": str(vtk_files[0]),
    }


def relative_delta(on: float, off: float) -> float:
    if off == 0.0:
        return math.nan
    return (on - off) / off


def analyze(root: Path, output_dir: Path) -> dict[str, Any]:
    metrics = {case: case_metrics(root / case) for case in CASES}
    off = metrics["ialmob_off"]
    on = metrics["ialmob_on"]
    off_peak = off["avalanche_peak"]
    on_peak = on["avalanche_peak"]
    peak_shift = math.hypot(
        float(on_peak["x_um"]) - float(off_peak["x_um"]),
        float(on_peak["y_um"]) - float(off_peak["y_um"]),
    )
    summary = {
        "schema": "vela.slot_ldmos.ialmob_ablation_result.v1",
        "controlled_delta": (
            "masetti_field_lombardi minus masetti_field; all other physics "
            "and numerical controls shared"
        ),
        **metrics,
        "on_minus_off": {
            "inner_voltage_V": on["inner_voltage_V"] - off["inner_voltage_V"],
            "drain_current_A_per_um": (
                on["drain_current_A_per_um"] - off["drain_current_A_per_um"]
            ),
            "drain_current_relative": relative_delta(
                on["drain_current_A_per_um"], off["drain_current_A_per_um"]
            ),
            "qG_full_A_per_um": on["qG_full_A_per_um"] - off["qG_full_A_per_um"],
            "qG_full_relative": relative_delta(
                on["qG_full_A_per_um"], off["qG_full_A_per_um"]
            ),
            "max_electric_field_relative": relative_delta(
                on["max_electric_field_V_per_cm"],
                off["max_electric_field_V_per_cm"],
            ),
            "max_avalanche_generation_relative": relative_delta(
                on["max_avalanche_generation_cm3_s"],
                off["max_avalanche_generation_cm3_s"],
            ),
            "avalanche_peak_shift_um": peak_shift,
        },
        "bvds": {
            "status": "not_reached_by_60V_outer_probe",
            "criterion_A_per_um": 1.0e-7,
            "maximum_probe_current_A_per_um": max(
                off["drain_current_A_per_um"], on["drain_current_A_per_um"]
            ),
        },
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "ialmob_probe_60v_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    delta = summary["on_minus_off"]
    markdown = "\n".join(
        [
            "# Vela Slot-LDMOS IALMob 60 V probe",
            "",
            "| Metric | IALMob off | IALMob on | On - off |",
            "|---|---:|---:|---:|",
            (
                "| Inner voltage (V) | "
                f"{off['inner_voltage_V']:.9f} | {on['inner_voltage_V']:.9f} | "
                f"{delta['inner_voltage_V']:+.9f} |"
            ),
            (
                "| Drain current (A/um) | "
                f"{off['drain_current_A_per_um']:.9e} | "
                f"{on['drain_current_A_per_um']:.9e} | "
                f"{delta['drain_current_relative'] * 100:+.6f}% |"
            ),
            (
                "| Integrated avalanche qG (A/um) | "
                f"{off['qG_full_A_per_um']:.9e} | "
                f"{on['qG_full_A_per_um']:.9e} | "
                f"{delta['qG_full_relative'] * 100:+.6f}% |"
            ),
            (
                "| Avalanche peak (x, y) um | "
                f"({off_peak['x_um']:.6g}, {off_peak['y_um']:.6g}) | "
                f"({on_peak['x_um']:.6g}, {on_peak['y_um']:.6g}) | "
                f"{peak_shift:.6g} um shift |"
            ),
            "",
            "BVDS is not inferred from this 60 V outer-load-line probe because "
            "neither case reaches 1e-7 A/um.",
            "",
        ]
    )
    (output_dir / "ialmob_probe_60v_summary.md").write_text(
        markdown, encoding="utf-8"
    )
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(analyze(args.root, args.output_dir), indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
