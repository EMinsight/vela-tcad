#!/usr/bin/env python3
"""Re-run the 21-point DD Id-Vg curve with the phase-7 DG shared settings.

The controlled model delta is deliberately limited to disabling the electron
quantum-potential equation.  Mesh, doping, material parameters, carrier
statistics, BGN, mobility, high-field driving force, contacts, bias points,
and nonlinear controls are inherited from the phase-7 regression helper.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[1]
BASELINE = REPO_ROOT / "build-release/reference_tcad/transportmodels_sentaurus2022/vela_baseline"
OUTPUT_ROOT = BASELINE / "dd_phase7_shared_baseline_2026-08-21"
REFERENCE = BASELINE / "generated/reference_curves/transportmodels_sentaurus2022_dd_idvg_reference.csv"
PHASE7_SCRIPT = REPO_ROOT / "scripts/run_transportmodels_dg_phase7_regression.py"
PHASE7_CONFIG = BASELINE / "dg_phase7_regression_2026-08-21/idvg/config.json"
PHASE7_STATE = BASELINE / "dg_phase7_regression_2026-08-21/idvg/state_bias_2p200000.csv"
REPORT_JSON = REPO_ROOT / "docs/validation/transportmodels_dd_phase7_shared_baseline_2026-08-21.json"
REPORT_MD = REPO_ROOT / "docs/validation/transportmodels_dd_phase7_shared_baseline_2026-08-21.md"


def load_phase7_module():
    spec = importlib.util.spec_from_file_location("transportmodels_phase7", PHASE7_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load {PHASE7_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sha256(path: Path) -> str:
    import hashlib

    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def finite_relative_error(candidate: float, reference: float) -> float | None:
    return abs(candidate - reference) / reference if reference > 1.0e-16 else None


def percentile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * fraction
    low = int(math.floor(position))
    high = int(math.ceil(position))
    if low == high:
        return ordered[low]
    return ordered[low] * (high - position) + ordered[high] * (position - low)


def aligned_rows(candidate_path: Path) -> list[dict[str, Any]]:
    with REFERENCE.open(newline="", encoding="utf-8-sig") as handle:
        reference = {
            round(float(row["bias_V"]), 12): abs(float(row["current_total"]))
            for row in csv.DictReader(handle)
        }
    with candidate_path.open(newline="", encoding="utf-8") as handle:
        candidate = list(csv.DictReader(handle))
    rows = []
    for row in candidate:
        bias = round(float(row["bias_V"]), 12)
        if bias not in reference:
            continue
        vela = abs(float(row["current_total_A_per_um"]))
        sentaurus = reference[bias]
        rows.append(
            {
                "bias_V": float(row["bias_V"]),
                "vela_A_per_um": vela,
                "sentaurus_A_per_um": sentaurus,
                "converged": row["converged"] == "1",
                "iterations": int(row["iterations"]),
                "absolute_log_error_dex": abs(
                    math.log10(max(vela, 1.0e-30))
                    - math.log10(max(sentaurus, 1.0e-30))
                ),
                "absolute_relative_error": finite_relative_error(vela, sentaurus),
            }
        )
    rows.sort(key=lambda item: item["bias_V"])
    return rows


def regime_metrics(rows: list[dict[str, Any]]) -> dict[str, Any]:
    regimes = {"off": rows[:3], "transition": rows[3:8], "on": rows[8:]}
    result: dict[str, Any] = {}
    for name, selected in regimes.items():
        relative = [
            float(row["absolute_relative_error"])
            for row in selected
            if row["absolute_relative_error"] is not None
        ]
        logarithmic = [float(row["absolute_log_error_dex"]) for row in selected]
        result[name] = {
            "points": len(selected),
            "max_absolute_log_error_dex": max(logarithmic),
            "median_absolute_log_error_dex": percentile(logarithmic, 0.5),
            "max_relative_error": max(relative) if relative else None,
            "median_relative_error": percentile(relative, 0.5) if relative else None,
        }
    return result


def make_plot(rows: list[dict[str, Any]]) -> tuple[Path, Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axis = plt.subplots(figsize=(7.6, 5.4))
    axis.semilogy(
        [row["bias_V"] for row in rows],
        [row["sentaurus_A_per_um"] for row in rows],
        "o-",
        label="Sentaurus 2022 DD",
    )
    axis.semilogy(
        [row["bias_V"] for row in rows],
        [row["vela_A_per_um"] for row in rows],
        "s--",
        label="Vela DD, phase-7 shared settings",
    )
    axis.set_xlabel("Gate voltage Vg (V)")
    axis.set_ylabel("Drain current Id (A/µm)")
    axis.set_title("TransportModels DD Id-Vg controlled baseline")
    axis.grid(True, which="both", alpha=0.25)
    axis.legend()
    fig.tight_layout()
    png = OUTPUT_ROOT / "dd_idvg_phase7_shared_comparison.png"
    svg = OUTPUT_ROOT / "dd_idvg_phase7_shared_comparison.svg"
    fig.savefig(png, dpi=180)
    fig.savefig(svg)
    plt.close(fig)
    return png, svg


def markdown(report: dict[str, Any]) -> str:
    metrics = report["metrics"]
    return "\n".join(
        [
            "# TransportModels DD Id-Vg phase-7 shared-settings baseline",
            "",
            f"Status: **{report['status']}**; completed `{report['completed_points']}/21` points.",
            "",
            "## Controlled delta",
            "",
            "The configuration is generated by the phase-7 DG configuration helper. "
            "Only `solver.electron_quantum_potential.enabled` is changed from `true` to `false`; "
            "all shared mesh, doping, material, statistics, BGN, mobility, high-field, contact, "
            "sweep, and nonlinear-solver settings remain identical.",
            "",
            "## Error summary",
            "",
            "| Regime | Points | Max log error (dex) | Median log error (dex) | Max relative error |",
            "|---|---:|---:|---:|---:|",
            *[
                f"| {name} | {row['points']} | {row['max_absolute_log_error_dex']:.6g} | "
                f"{row['median_absolute_log_error_dex']:.6g} | "
                + (f"{row['max_relative_error']:.6%} |" if row["max_relative_error"] is not None else "n/a |")
                for name, row in metrics.items()
            ],
            "",
            f"Figure: `{report['artifacts']['png']}`",
            "",
        ]
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--report-only", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.check:
        report = json.loads(REPORT_JSON.read_text(encoding="utf-8"))
        assert sha256(Path(report["candidate_csv"])) == report["candidate_sha256"]
        assert sha256(Path(report["reference_csv"])) == report["reference_sha256"]
        assert report["controlled_delta"]["electron_quantum_potential_enabled"] == [True, False]
        print("TransportModels DD phase-7 shared baseline check: PASS")
        return 0

    phase7 = load_phase7_module()
    phase7.OUTPUT_ROOT = OUTPUT_ROOT
    curve = {
        "name": "idvg",
        "label": "DD Id-Vg with phase-7 shared settings",
        "config": phase7.GENERATED / "simulation_dg_idvg.json",
        "reference": REFERENCE,
        "contact": "gate",
        "current_contact": "drain",
        "points": [2.2 - 0.16 * index for index in range(21)],
        "initial_state": PHASE7_STATE,
    }
    original_make_config = phase7.make_config

    def make_dd_config(active_curve: dict[str, Any]) -> Path:
        path = original_make_config(active_curve)
        config = json.loads(path.read_text(encoding="utf-8"))
        config["_comment"] = "TransportModels DD 21-point Id-Vg with phase-7 shared settings"
        config["solver"]["electron_quantum_potential"]["enabled"] = False
        path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
        return path

    phase7.make_config = make_dd_config
    execution = None
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    if not args.report_only:
        execution = phase7.execute_curve(curve)

    candidate = OUTPUT_ROOT / "idvg/curve_combined.csv"
    if not candidate.exists():
        candidate = OUTPUT_ROOT / "idvg/curve.csv"
    rows = aligned_rows(candidate)
    complete = len(rows) == 21 and all(row["converged"] for row in rows)
    config_path = OUTPUT_ROOT / "idvg/config.json"
    config = json.loads(config_path.read_text(encoding="utf-8"))
    phase7_config = json.loads(PHASE7_CONFIG.read_text(encoding="utf-8"))
    png, svg = make_plot(rows)
    report = {
        "schema": "vela.transportmodels.dd_phase7_shared_baseline.v1",
        "as_of": "2026-08-21",
        "status": "complete" if complete else "partial",
        "completed_points": len(rows),
        "all_converged": complete,
        "execution": execution,
        "controlled_delta": {
            "source_phase7_config": str(PHASE7_CONFIG.resolve()),
            "dd_config": str(config_path.resolve()),
            "electron_quantum_potential_enabled": [
                phase7_config["solver"]["electron_quantum_potential"]["enabled"],
                config["solver"]["electron_quantum_potential"]["enabled"],
            ],
            "materials_identical": phase7_config["materials_file"] == config["materials_file"],
            "mobility_identical": phase7_config["solver"]["mobility"] == config["solver"]["mobility"],
            "carrier_statistics_identical": phase7_config["solver"]["carrier_statistics"] == config["solver"]["carrier_statistics"],
            "bandgap_narrowing_identical": phase7_config["solver"]["bandgap_narrowing"] == config["solver"]["bandgap_narrowing"],
            "bias_points_identical": phase7_config["sweep"]["bias_points"] == config["sweep"]["bias_points"],
        },
        "metrics": regime_metrics(rows) if len(rows) == 21 else None,
        "aligned": rows,
        "candidate_csv": str(candidate.resolve()),
        "candidate_sha256": sha256(candidate),
        "reference_csv": str(REFERENCE.resolve()),
        "reference_sha256": sha256(REFERENCE),
        "artifacts": {"png": str(png.resolve()), "svg": str(svg.resolve())},
    }
    REPORT_JSON.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    REPORT_MD.write_text(markdown(report), encoding="utf-8")
    print(
        json.dumps(
            {
                "status": report["status"],
                "completed_points": report["completed_points"],
                "controlled_delta": report["controlled_delta"],
                "metrics": report["metrics"],
            },
            indent=2,
        )
    )
    return 0 if complete else 1


if __name__ == "__main__":
    raise SystemExit(main())
