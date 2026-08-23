#!/usr/bin/env python3
"""Extract comparable BV results from the Sentaurus BVmethods example.

The official Workbench project maps its six SDevice methods to nodes 3--8:
ABA Poisson, ABA with coupled carriers (IIC), external resistor,
voltage-to-current switching, continuation, and transient.  This tool reads
the generated DF-ISE ``n*_des.plt`` files without requiring Inspect, writes
normalized per-method curves, and reproduces the extraction criteria used by
the example's ``inspect_ins.cmd``.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any, Callable

from sentaurus_import import parse_quoted_list, parse_values_block


ELEMENTARY_CHARGE_C = 1.602176634e-19
DEFAULT_CURRENT_THRESHOLD_A_PER_UM = 1.0e-4
DEFAULT_ION_INTEGRAL_THRESHOLD = 1.05

METHOD_NODES = {
    "ABA_poisson": 3,
    "ABA_coupled": 4,
    "resistor": 5,
    "voltage2current": 6,
    "continuation": 7,
    "transient": 8,
}


def interpolate_x(x0: float, y0: float, x1: float, y1: float, target: float) -> float:
    if y1 == y0:
        return x1
    fraction = (target - y0) / (y1 - y0)
    return x0 + fraction * (x1 - x0)


def first_upward_crossing(
    xs: list[float], ys: list[float], target: float,
    *, minimum_x: float = -math.inf,
) -> float | None:
    for index in range(1, len(xs)):
        if xs[index] < minimum_x:
            continue
        y0 = ys[index - 1]
        y1 = ys[index]
        if y0 < target <= y1:
            return interpolate_x(xs[index - 1], y0, xs[index], y1, target)
    return None


def read_plt(path: Path) -> tuple[list[str], list[list[float]]]:
    text = path.read_text(errors="ignore")
    datasets = parse_quoted_list(text, "datasets")
    if not datasets:
        raise ValueError(f"{path} does not declare datasets")
    return datasets, parse_values_block(text, len(datasets))


def column(
    datasets: list[str], rows: list[list[float]], name: str,
) -> list[float]:
    try:
        index = datasets.index(name)
    except ValueError as exc:
        raise ValueError(f"missing dataset {name!r}") from exc
    return [row[index] for row in rows]


def optional_column(
    datasets: list[str], rows: list[list[float]], name: str,
) -> list[float | None]:
    if name not in datasets:
        return [None] * len(rows)
    index = datasets.index(name)
    return [row[index] for row in rows]


def column_alias(
    datasets: list[str], rows: list[list[float]], names: tuple[str, ...],
) -> list[float]:
    for name in names:
        if name in datasets:
            return column(datasets, rows, name)
    raise ValueError(f"missing all dataset aliases {names!r}")


def optional_column_alias(
    datasets: list[str], rows: list[list[float]], names: tuple[str, ...],
) -> list[float | None]:
    for name in names:
        if name in datasets:
            return optional_column(datasets, rows, name)
    return [None] * len(rows)


IMPACT_GENERATION_DATASETS = (
    "IntegrSemiconductor ImpactIonization",
    "IntegrSemiconductor AvalancheGeneration",
)


def write_curve(
    path: Path,
    datasets: list[str],
    rows: list[list[float]],
) -> None:
    inner = column(datasets, rows, "drain InnerVoltage")
    outer = column(datasets, rows, "drain OuterVoltage")
    current = column(datasets, rows, "drain TotalCurrent")
    phi_e = optional_column(datasets, rows, "PhiElectron")
    phi_h = optional_column(datasets, rows, "PhiHole")
    integrated = optional_column_alias(
        datasets, rows, IMPACT_GENERATION_DATASETS)

    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow([
            "inner_voltage_V",
            "outer_voltage_V",
            "drain_total_current_A_per_um",
            "phi_electron",
            "phi_hole",
            "integrated_avalanche_generation_per_s_um",
            "avalanche_current_A_per_um",
        ])
        for values in zip(inner, outer, current, phi_e, phi_h, integrated):
            generation = values[5]
            avalanche_current = (
                generation * 1.0e-12 * ELEMENTARY_CHARGE_C
                if generation is not None else None
            )
            writer.writerow([*values, avalanche_current])


def analyze_method(
    method: str,
    datasets: list[str],
    rows: list[list[float]],
    current_threshold: float,
    ion_integral_threshold: float,
) -> dict[str, Any]:
    inner = column(datasets, rows, "drain InnerVoltage")
    outer = column(datasets, rows, "drain OuterVoltage")
    current = column(datasets, rows, "drain TotalCurrent")

    result: dict[str, Any] = {
        "method": method,
        "rows": len(rows),
        "last_inner_voltage_V": inner[-1],
        "last_outer_voltage_V": outer[-1],
        "last_drain_total_current_A_per_um": current[-1],
    }

    if method == "ABA_poisson":
        phi_e = column(datasets, rows, "PhiElectron")
        phi_h = column(datasets, rows, "PhiHole")
        electron_bv = first_upward_crossing(inner, phi_e, ion_integral_threshold)
        hole_bv = first_upward_crossing(inner, phi_h, ion_integral_threshold)
        candidates = [value for value in (electron_bv, hole_bv) if value is not None]
        result.update({
            "criterion": f"min(PhiElectron,PhiHole) crossing {ion_integral_threshold}",
            "electron_bv_V": electron_bv,
            "hole_bv_V": hole_bv,
            "bv_V": min(candidates) if candidates else None,
        })
        return result

    if method == "ABA_coupled":
        generation = column_alias(
            datasets, rows, IMPACT_GENERATION_DATASETS)
        avalanche_current = [
            value * 1.0e-12 * ELEMENTARY_CHARGE_C for value in generation
        ]
        difference = [
            avalanche - terminal
            for avalanche, terminal in zip(avalanche_current, current)
        ]
        result.update({
            "criterion": "integrated avalanche current equals drain conduction current",
            "bv_V": first_upward_crossing(inner, difference, 0.0, minimum_x=1.0),
            "last_avalanche_current_A_per_um": avalanche_current[-1],
        })
        return result

    result.update({
        "criterion": f"abs(drain current) crossing {current_threshold} A/um",
        "bv_V": first_upward_crossing(
            inner, [abs(value) for value in current], current_threshold),
    })
    return result


def write_summary_csv(path: Path, results: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "method", "node", "rows", "bv_V", "criterion",
        "last_inner_voltage_V", "last_outer_voltage_V",
        "last_drain_total_current_A_per_um",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(results)


def write_summary_markdown(path: Path, summary: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Sentaurus BVmethods extraction",
        "",
        f"Input directory: `{summary['input_dir']}`",
        "",
        "| Method | Node | Rows | Extracted BV (V) | Criterion |",
        "| --- | ---: | ---: | ---: | --- |",
    ]
    for result in summary["results"]:
        value = result.get("bv_V")
        rendered = "n/a" if value is None else f"{value:.6g}"
        lines.append(
            f"| {result['method']} | {result['node']} | {result['rows']} | "
            f"{rendered} | {result['criterion']} |")
    lines.extend([
        "",
        "The voltage reported here is the drain `InnerVoltage`. For the external",
        "resistor and transient methods, `OuterVoltage` includes the voltage drop",
        "across the series resistor and is not the device breakdown voltage.",
        "",
    ])
    path.write_text("\n".join(lines), encoding="utf-8")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--current-threshold-A-per-um", type=float,
        default=DEFAULT_CURRENT_THRESHOLD_A_PER_UM)
    parser.add_argument(
        "--ion-integral-threshold", type=float,
        default=DEFAULT_ION_INTEGRAL_THRESHOLD)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    results: list[dict[str, Any]] = []
    curves_dir = args.output_dir / "curves"

    for method, node in METHOD_NODES.items():
        path = args.input_dir / f"n{node}_des.plt"
        datasets, rows = read_plt(path)
        write_curve(curves_dir / f"{method}.csv", datasets, rows)
        result = analyze_method(
            method,
            datasets,
            rows,
            args.current_threshold_A_per_um,
            args.ion_integral_threshold,
        )
        result["node"] = node
        result["source_plt"] = str(path)
        results.append(result)

    summary = {
        "schema": "vela.sentaurus_bvmethods_analysis.v1",
        "input_dir": str(args.input_dir.resolve()),
        "current_threshold_A_per_um": args.current_threshold_A_per_um,
        "ion_integral_threshold": args.ion_integral_threshold,
        "results": results,
    }
    args.output_dir.mkdir(parents=True, exist_ok=True)
    (args.output_dir / "bvmethods_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    write_summary_csv(args.output_dir / "bvmethods_summary.csv", results)
    write_summary_markdown(args.output_dir / "bvmethods_summary.md", summary)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
