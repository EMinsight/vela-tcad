#!/usr/bin/env python3
"""Summarize the staged Slot-LDMOS Sentaurus BVDS reference run."""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))

from sentaurus_import import parse_quoted_list, parse_values_block  # noqa: E402


STAGES = (
    "equilibrium",
    "unit_resistor_1v_r1e12",
    "unit_resistor_1v_direct",
    "avalanche_off_60v",
    "iic_postprocess_60v",
    "avalanche_on_60v",
    "bvds_external_resistor_final",
)
RESISTANCE_OHM_UM = 1.0e12
BREAK_CURRENT_A_PER_UM = 1.0e-7


def read_plt(path: Path) -> tuple[list[str], list[list[float]]]:
    text = path.read_text(errors="ignore")
    datasets = parse_quoted_list(text, "datasets")
    if not datasets:
        raise ValueError(f"{path} does not declare datasets")
    rows = parse_values_block(text, len(datasets))
    if not rows:
        raise ValueError(f"{path} has no complete data rows")
    return datasets, rows


def value(row: list[float], datasets: list[str], name: str) -> float:
    try:
        return row[datasets.index(name)]
    except ValueError as error:
        raise ValueError(f"missing dataset {name!r}") from error


def curve_rows(path: Path) -> list[dict[str, float]]:
    datasets, rows = read_plt(path)
    return [
        {
            "time": value(row, datasets, "time"),
            "outer_voltage_V": value(row, datasets, "drain OuterVoltage"),
            "inner_voltage_V": value(row, datasets, "drain InnerVoltage"),
            "drain_total_current_A_per_um": value(
                row, datasets, "drain TotalCurrent"
            ),
        }
        for row in rows
    ]


def interpolate_crossing(
    rows: list[dict[str, float]], threshold: float
) -> dict[str, float] | None:
    for before, after in zip(rows, rows[1:]):
        y0 = abs(before["drain_total_current_A_per_um"])
        y1 = abs(after["drain_total_current_A_per_um"])
        if y0 < threshold <= y1:
            fraction = 1.0 if y1 == y0 else (threshold - y0) / (y1 - y0)
            return {
                key: before[key] + fraction * (after[key] - before[key])
                for key in (
                    "time",
                    "outer_voltage_V",
                    "inner_voltage_V",
                    "drain_total_current_A_per_um",
                )
            }
    return None


def load_line(row: dict[str, float]) -> dict[str, float]:
    drop = row["outer_voltage_V"] - row["inner_voltage_V"]
    predicted = row["drain_total_current_A_per_um"] * RESISTANCE_OHM_UM
    return {
        "voltage_drop_V": drop,
        "current_times_resistance_V": predicted,
        "residual_V": drop - predicted,
    }


def parse_iic_log(path: Path) -> dict[str, float] | None:
    text = path.read_text(errors="ignore")
    matches = re.findall(
        r"Maximum electric field:\s*([0-9.eE+-]+) V/cm at "
        r"\(([0-9.eE+-]+),([0-9.eE+-]+)\) um.*?"
        r"Electron:\s*([0-9.eE+-]+)\s*Hole:\s*([0-9.eE+-]+)",
        text,
        flags=re.DOTALL,
    )
    if not matches:
        return None
    field, x_um, y_um, electron, hole = matches[-1]
    return {
        "maximum_electric_field_V_per_cm": float(field),
        "maximum_field_x_um": float(x_um),
        "maximum_field_y_um": float(y_um),
        "electron_ionization_integral": float(electron),
        "hole_ionization_integral": float(hole),
    }


def parse_mesh_log(path: Path) -> dict[str, float | int] | None:
    text = path.read_text(errors="ignore")
    match = re.search(
        r"^\s*Total\s+([0-9.eE+-]+)\s+([0-9.eE+-]+)\s+"
        r"([0-9.eE+-]+)\s+(\d+)\s+(\d+)\s*\([^)]*\)\s+"
        r"([0-9.eE+-]+)\s*\(\s*([0-9.eE+-]+)\s*\)",
        text,
        flags=re.MULTILINE,
    )
    if not match:
        return None
    volume, box_volume, delta, cells, non_delaunay, nd_volume, nd_percent = (
        match.groups()
    )
    return {
        "geometric_volume_um2": float(volume),
        "box_method_volume_um2": float(box_volume),
        "delta_volume_percent": float(delta),
        "cell_count": int(cells),
        "non_delaunay_cell_count": int(non_delaunay),
        "non_delaunay_volume_um2": float(nd_volume),
        "non_delaunay_volume_percent": float(nd_percent),
    }


def write_curve(path: Path, rows: list[dict[str, float]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def markdown(summary: dict[str, Any]) -> str:
    lines = [
        "# Slot-LDMOS Sentaurus BVDS staged reference summary",
        "",
        "| Stage | Points | Outer V | Inner V | Drain current (A/um) |",
        "|---|---:|---:|---:|---:|",
    ]
    for name in STAGES:
        item = summary["stages"][name]
        final = item["final"]
        lines.append(
            f"| `{name}` | {item['point_count']} | "
            f"{final['outer_voltage_V']:.9g} | {final['inner_voltage_V']:.9g} | "
            f"{final['drain_total_current_A_per_um']:.9g} |"
        )
    calibration = summary["unit_calibration"]
    endpoint = summary["breakdown_current_interpolation"]
    mesh = summary["sentaurus_mesh_statistics"]
    iic = summary["iic_diagnostics"]
    lines.extend(
        [
            "",
            "## Acceptance metrics",
            "",
            f"- 1 Tohm load-line residual: "
            f"{calibration['resistor_load_line']['residual_V']:.6e} V.",
            f"- Direct/resistor 1 V current relative difference: "
            f"{calibration['direct_resistor_current_relative_difference']:.6e}.",
            f"- Interpolated inner BVDS at 1e-7 A/um: "
            f"{endpoint['inner_voltage_V']:.9g} V.",
            f"- Interpolated outer voltage at 1e-7 A/um: "
            f"{endpoint['outer_voltage_V']:.9g} V.",
            f"- Sentaurus non-Delaunay cells: {mesh['non_delaunay_cell_count']} "
            f"of {mesh['cell_count']}; DeltaVolume={mesh['delta_volume_percent']:.6e}%.",
            f"- IIC best-path maximum field: "
            f"{iic['maximum_electric_field_V_per_cm']:.6e} V/cm at "
            f"({iic['maximum_field_x_um']:.6g}, {iic['maximum_field_y_um']:.6g}) um.",
            f"- IIC electron/hole integrals: "
            f"{iic['electron_ionization_integral']:.6g} / "
            f"{iic['hole_ionization_integral']:.6g}.",
            "",
        ]
    )
    return "\n".join(lines)


def analyze(input_dir: Path, output_dir: Path) -> dict[str, Any]:
    output_dir.mkdir(parents=True, exist_ok=True)
    stages: dict[str, Any] = {}
    curves: dict[str, list[dict[str, float]]] = {}
    for name in STAGES:
        rows = curve_rows(input_dir / f"{name}.plt")
        curves[name] = rows
        write_curve(output_dir / f"{name}.csv", rows)
        stages[name] = {
            "point_count": len(rows),
            "final": rows[-1],
            "exit_code": int((input_dir / f"{name}.exitcode").read_text().strip()),
        }

    resistor = curves["unit_resistor_1v_r1e12"][-1]
    direct = curves["unit_resistor_1v_direct"][-1]
    relative_difference = abs(
        direct["drain_total_current_A_per_um"]
        - resistor["drain_total_current_A_per_um"]
    ) / max(abs(direct["drain_total_current_A_per_um"]), sys.float_info.min)
    endpoint = interpolate_crossing(
        curves["bvds_external_resistor_final"], BREAK_CURRENT_A_PER_UM
    )
    if endpoint is None:
        raise ValueError("final BVDS curve does not cross 1e-7 A/um")

    summary = {
        "schema": "vela.slot_ldmos.sentaurus_bvds_reference.v1",
        "current_unit": "A/um",
        "external_resistance_unit": "ohm*um",
        "external_resistance_ohm_um": RESISTANCE_OHM_UM,
        "stages": stages,
        "unit_calibration": {
            "resistor_load_line": load_line(resistor),
            "direct_resistor_current_relative_difference": relative_difference,
        },
        "breakdown_current_interpolation": endpoint,
        "sentaurus_mesh_statistics": parse_mesh_log(
            input_dir / "equilibrium.log_des.log"
        ),
        "iic_diagnostics": parse_iic_log(
            input_dir / "iic_postprocess_60v.log_des.log"
        ),
    }
    (output_dir / "slot_ldmos_bvds_reference_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    (output_dir / "slot_ldmos_bvds_reference_summary.md").write_text(
        markdown(summary), encoding="utf-8"
    )
    return summary


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    summary = analyze(args.input_dir.resolve(), args.output_dir.resolve())
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
