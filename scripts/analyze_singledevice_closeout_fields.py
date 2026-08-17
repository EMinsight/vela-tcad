#!/usr/bin/env python3
"""Compare three exact-bias SingleDevice states and audit terminal KCL."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


BIAS_SPECS = (
    (-0.5, "m0p500000", "vg_m0p5"),
    (0.31, "0p310000", "vg_0p31"),
    (2.2, "2p200000", "vg_2p2"),
)


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def percentile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    if not ordered:
        return 0.0
    position = fraction * (len(ordered) - 1)
    lo, hi = math.floor(position), math.ceil(position)
    if lo == hi:
        return ordered[lo]
    return ordered[lo] * (hi - position) + ordered[hi] * (position - lo)


def abs_summary(reference: dict[int, float], candidate: dict[int, float]) -> dict:
    errors = [abs(candidate[node] - value) for node, value in reference.items()
              if node in candidate and math.isfinite(value) and math.isfinite(candidate[node])]
    return {
        "count": len(errors),
        "median_abs_error": percentile(errors, 0.5),
        "p95_abs_error": percentile(errors, 0.95),
        "max_abs_error": max(errors, default=0.0),
    }


def log_summary(reference: dict[int, float], candidate: dict[int, float]) -> dict:
    errors = [abs(math.log10(candidate[node] / value))
              for node, value in reference.items()
              if node in candidate and value > 0.0 and candidate[node] > 0.0]
    return {
        "count": len(errors),
        "median_abs_log10_error_dex": percentile(errors, 0.5),
        "p95_abs_log10_error_dex": percentile(errors, 0.95),
        "max_abs_log10_error_dex": max(errors, default=0.0),
    }


def scalar_csv(path: Path) -> dict[int, float]:
    return {int(row["node_id"]): float(row["component0"]) for row in rows(path)}


def sentaurus_field(export: Path, name: str, regions: tuple[int, ...]) -> dict[int, float]:
    result: dict[int, float] = {}
    for region in regions:
        path = export / "fields" / f"{name}_region{region}.csv"
        if path.is_file():
            result.update(scalar_csv(path))
    if not result:
        raise FileNotFoundError(f"no {name} fields below {export}")
    return result


def state_fields(path: Path) -> dict[str, dict[int, float]]:
    result = {name: {} for name in (
        "psi", "phin", "phip", "electrons_m3", "holes_m3",
        "electron_quantum_potential_V")}
    for row in rows(path):
        node = int(row["node_id"])
        for name in result:
            if name in row and row[name] not in (None, ""):
                result[name][node] = float(row[name])
    return result


def kcl_at(path: Path, bias: float) -> dict:
    selected = [row for row in rows(path)
                if abs(float(row["bias_V"]) - bias) <= 1.0e-10]
    if len(selected) != 4:
        raise ValueError(f"expected four terminals at {bias} V in {path}")
    currents = {row["contact"]: float(row["current_total_A_per_um"])
                for row in selected}
    residual = abs(sum(currents.values()))
    drain = abs(currents["drain"])
    ratio = residual / max(drain, 1.0e-300)
    return {
        "terminal_currents_A_per_um": currents,
        "absolute_residual_A_per_um": residual,
        "residual_over_drain_current": ratio,
        "pass": residual <= 1.0e-14 or ratio <= 0.01,
        "rule": "absolute_residual <= 1e-14 A/um OR residual/drain <= 1%",
    }


def analyze_point(state: Path, export: Path, terminal: Path, bias: float) -> dict:
    vela = state_fields(state)
    sent = {
        "psi": sentaurus_field(export, "ElectrostaticPotential", tuple(range(7))),
        # Sentaurus writes placeholder quasi-Fermi values in insulators while
        # Vela pins those algebraic rows.  Compare transport regions only.
        "phin": sentaurus_field(export, "eQuasiFermiPotential", (3, 4)),
        "phip": sentaurus_field(export, "hQuasiFermiPotential", (3, 4)),
        "electrons_cm3": sentaurus_field(export, "eDensity", (3, 4)),
        "holes_cm3": sentaurus_field(export, "hDensity", (3, 4)),
        "electron_quantum_potential_V": sentaurus_field(
            export, "eQuantumPotential", tuple(range(7))),
    }
    vela_e_cm3 = {node: value / 1.0e6 for node, value in vela["electrons_m3"].items()}
    vela_h_cm3 = {node: value / 1.0e6 for node, value in vela["holes_m3"].items()}
    return {
        "bias": {"Vg_V": bias, "Vds_V": 0.1},
        "state_file": str(state.resolve()),
        "sentaurus_export": str(export.resolve()),
        "fields": {
            "electrostatic_potential_V": abs_summary(sent["psi"], vela["psi"]),
            "electron_quasi_fermi_V": abs_summary(sent["phin"], vela["phin"]),
            "hole_quasi_fermi_V": abs_summary(sent["phip"], vela["phip"]),
            "electron_quantum_potential_V": abs_summary(
                sent["electron_quantum_potential_V"],
                vela["electron_quantum_potential_V"]),
            "electron_density": log_summary(sent["electrons_cm3"], vela_e_cm3),
            "hole_density": log_summary(sent["holes_cm3"], vela_h_cm3),
        },
        "kcl": kcl_at(terminal, bias),
    }


def markdown(report: dict) -> str:
    lines = [
        "# SingleDevice three-bias field and KCL audit", "",
        "All field comparisons use exact Sentaurus checkpoints on the original mesh.", "",
        "| Vg (V) | psi median/p95 (V) | eDensity median/p95 (dex) | eQP median/p95 (V) | KCL abs (A/um) | KCL ratio |", "| ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for point in report["points"]:
        field = point["fields"]
        psi = field["electrostatic_potential_V"]
        density = field["electron_density"]
        quantum = field["electron_quantum_potential_V"]
        kcl = point["kcl"]
        quantum_text = (f"{quantum['median_abs_error']:.6g} / "
                        f"{quantum['p95_abs_error']:.6g}"
                        if quantum["count"] else "n/a")
        lines.append(
            f"| {point['bias']['Vg_V']:.2f} | {psi['median_abs_error']:.6g} / {psi['p95_abs_error']:.6g} | "
            f"{density['median_abs_log10_error_dex']:.6g} / {density['p95_abs_log10_error_dex']:.6g} | "
            f"{quantum_text} | "
            f"{kcl['absolute_residual_A_per_um']:.6g} | {kcl['residual_over_drain_current']:.6g} |")
    lines += ["", f"KCL gate: **{'PASS' if report['kcl_all_pass'] else 'FAIL'}**.",
              "Field metrics are diagnostic and do not authorize a new physics model.", ""]
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workflow-dir", required=True, type=Path)
    parser.add_argument("--sentaurus-root", required=True, type=Path)
    parser.add_argument("--output-json", required=True, type=Path)
    parser.add_argument("--output-md", required=True, type=Path)
    args = parser.parse_args()
    terminal = args.workflow_dir / "11_linear_idvg_terminal_balance.csv"
    points = [analyze_point(
        args.workflow_dir / f"11_linear_idvg_state_bias_{state_suffix}.csv",
        args.sentaurus_root / sentaurus_suffix, terminal, bias)
              for bias, state_suffix, sentaurus_suffix in BIAS_SPECS]
    report = {
        "schema": "vela.singledevice.three_bias_field_kcl.v1",
        "points": points,
        "kcl_all_pass": all(point["kcl"]["pass"] for point in points),
        "field_acceptance": "diagnostic_only",
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    args.output_md.parent.mkdir(parents=True, exist_ok=True)
    args.output_md.write_text(markdown(report), encoding="utf-8")
    print(json.dumps({"kcl_all_pass": report["kcl_all_pass"], "points": 3}))
    return 0 if report["kcl_all_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
