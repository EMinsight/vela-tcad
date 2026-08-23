#!/usr/bin/env python3
"""A/B the current DG contract against the earlier Sentaurus-box contract."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import subprocess
from pathlib import Path
from typing import Any, Iterable


REPO = Path(__file__).resolve().parents[1]
REF = REPO / "build-release/reference_tcad/transportmodels_sentaurus2022"
BASE = REF / "vela_baseline/dd_dg_srh_corrected_cold_regression_2026-08-23/runs/dg"
OLD_CONFIG = REF / "vela_baseline/vela_fermi_bgn_ab_2026-08-21/dg_on/config.json"
SENT = REF / "sentaurus_vm_runs/remaining_spatial_oracles_20260823/exports"
OUTPUT = REF / "reports/transportmodels_quantum_contract_ab_20260823"
DEFAULT_RUNNER = REPO / "build-release/vela_example_runner_quantum_ab.exe"
REPORT_JSON = REPO / "docs/validation/transportmodels_quantum_contract_ab_2026-08-23.json"
REPORT_MD = REPO / "docs/validation/transportmodels_quantum_contract_ab_2026-08-23.md"

POINTS = (
    {
        "name": "idvd_vg1_vd2",
        "base_config": BASE / "05_dg_idvd_curve.json",
        "baseline_state": BASE / "dg_idvd_curve_state_bias_2p000000.csv",
        "export": SENT / "vd_p2p00",
        "contact": "drain",
        "bias": 2.0,
        "gate": 1.0,
        "drain": 2.0,
    },
    {
        "name": "idvg_vgm1_vd1p1",
        "base_config": BASE / "03_dg_idvg_curve.json",
        "baseline_state": BASE / "dg_idvg_final_bias_relax_final_state.csv",
        "export": SENT / "vg_m1p00",
        "contact": "gate",
        "bias": -1.0,
        "gate": -1.0,
        "drain": 1.1,
    },
)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def percentile(values: Iterable[float], fraction: float) -> float:
    ordered = sorted(float(value) for value in values if math.isfinite(float(value)))
    if not ordered:
        return math.nan
    position = (len(ordered) - 1) * fraction
    lo, hi = math.floor(position), math.ceil(position)
    if lo == hi:
        return ordered[lo]
    return ordered[lo] * (hi - position) + ordered[hi] * (position - lo)


def stats(values: Iterable[float]) -> dict[str, float]:
    data = list(values)
    return {"median": percentile(data, 0.5), "p95": percentile(data, 0.95), "maximum": max(data)}


def sent_field(export: Path, name: str) -> dict[int, float]:
    return {
        int(row["node_id"]): float(row["component0"])
        for row in read_csv(export / "fields" / f"{name}_region3.csv")
    }


def state_metrics(state_path: Path, export: Path) -> dict[str, Any]:
    state = {int(row["node_id"]): row for row in read_csv(state_path)}
    sent_qn = sent_field(export, "eQuantumPotential")
    sent_n = sent_field(export, "eDensity")
    sent_phin = sent_field(export, "eQuasiFermiPotential")
    nodes = set(sent_n)
    return {
        "qn_abs_error_mV": stats(
            1.0e3 * abs(float(state[node]["electron_quantum_potential_V"]) - sent_qn[node])
            for node in nodes
        ),
        "electron_density_abs_error_dex": stats(
            abs(
                math.log10(max(float(state[node]["electrons_m3"]) / 1.0e6, 1.0))
                - math.log10(max(sent_n[node], 1.0))
            )
            for node in nodes
        ),
        "phin_abs_error_mV": stats(
            1.0e3 * abs(float(state[node]["phin"]) - sent_phin[node]) for node in nodes
        ),
    }


def execute(point: dict[str, Any], quantum: dict[str, Any], runner: Path) -> dict[str, Any]:
    run_dir = OUTPUT / point["name"] / "sentaurus_box_insulator"
    run_dir.mkdir(parents=True, exist_ok=True)
    config = json.loads(Path(point["base_config"]).read_text(encoding="utf-8"))
    config["solver"]["electron_quantum_potential"] = quantum
    config["solver"]["verbose"] = False
    for contact in config["contacts"]:
        if contact["name"].lower() == "gate":
            contact["bias"] = point["gate"]
        elif contact["name"].lower() == "drain":
            contact["bias"] = point["drain"]
    config["output_csv"] = str((run_dir / "curve.csv").resolve())
    config["log_file"] = str((run_dir / "runner.log").resolve())
    sweep = config["sweep"]
    sweep.update(
        {
            "contact": point["contact"],
            "start": point["bias"],
            "stop": point["bias"],
            "step": 0.1,
            "bias_points": [point["bias"]],
            "initial_state_file": str(Path(point["baseline_state"]).resolve()),
            "write_vtk": False,
            "write_state_file": str((run_dir / "final_state.csv").resolve()),
            "write_state_every_point_prefix": str((run_dir / "state").resolve()),
        }
    )
    sweep.pop("diagnostics", None)
    config["_comment"] = (
        "Corrected material/SRH contract with the earlier include-insulators "
        "Sentaurus-box DG discretization"
    )
    config_path = run_dir / "config.json"
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    env = os.environ.copy()
    env["PATH"] = r"D:\msys64\ucrt64\bin" + os.pathsep + env.get("PATH", "")
    completed = subprocess.run(
        [str(runner), "--config", str(config_path)],
        cwd=REPO,
        text=True,
        capture_output=True,
        env=env,
        check=False,
    )
    (run_dir / "stdout.txt").write_text(completed.stdout, encoding="utf-8")
    (run_dir / "stderr.txt").write_text(completed.stderr, encoding="utf-8")
    curve = read_csv(run_dir / "curve.csv") if (run_dir / "curve.csv").is_file() else []
    state_path = run_dir / "final_state.csv"
    return {
        "returncode": completed.returncode,
        "converged": bool(curve) and curve[-1].get("converged") == "1",
        "current_A_per_um": float(curve[-1]["current_total_A_per_um"]) if curve else math.nan,
        "state_metrics": state_metrics(state_path, Path(point["export"])) if state_path.is_file() else None,
        "config": str(config_path.resolve()),
        "state": str(state_path.resolve()),
    }


def existing_result(point: dict[str, Any]) -> dict[str, Any]:
    run_dir = OUTPUT / point["name"] / "sentaurus_box_insulator"
    curve = read_csv(run_dir / "curve.csv")
    state_path = run_dir / "final_state.csv"
    return {
        "returncode": 0,
        "converged": curve[-1].get("converged") == "1",
        "current_A_per_um": float(curve[-1]["current_total_A_per_um"]),
        "state_metrics": state_metrics(state_path, Path(point["export"])),
        "config": str((run_dir / "config.json").resolve()),
        "state": str(state_path.resolve()),
    }


def markdown(report: dict[str, Any]) -> str:
    lines = [
        "# TransportModels DG quantum-contract A/B",
        "",
        "| Point | Variant | Converged | Current (A/um) | Id error | Qn p95 (mV) | n p95 (dex) | phin p95 (mV) |",
        "|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for point in report["points"]:
        for name in ("current_baseline", "sentaurus_box_insulator"):
            row = point[name]
            metrics = row["state_metrics"]
            lines.append(
                f"| {point['name']} | {name} | {row['converged']} | {row['current_A_per_um']:.8g} | "
                f"{row['current_absolute_relative_error']:.3%} | "
                f"{metrics['qn_abs_error_mV']['p95']:.5g} | "
                f"{metrics['electron_density_abs_error_dex']['p95']:.5g} | "
                f"{metrics['phin_abs_error_mV']['p95']:.5g} |"
            )
    lines.append("")
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, default=DEFAULT_RUNNER)
    parser.add_argument("--report-only", action="store_true")
    args = parser.parse_args()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    old_quantum = json.loads(OLD_CONFIG.read_text(encoding="utf-8"))["solver"]["electron_quantum_potential"]
    results = []
    for point in POINTS:
        baseline_curve = read_csv(
            BASE / ("dg_idvd_curve_comparison_candidate.csv" if point["contact"] == "drain" else "dg_idvg_curve_comparison_candidate.csv")
        )
        current = next(float(row["current_total_A_per_um"]) for row in baseline_curve if math.isclose(float(row["bias_V"]), point["bias"]))
        sentaurus_curve = read_csv(
            REF / "run02/normalized" / ("dg_idvd.csv" if point["contact"] == "drain" else "dg_idvg.csv")
        )
        sentaurus_current = next(
            float(row["current_total"])
            for row in sentaurus_curve
            if math.isclose(float(row["bias_V"]), point["bias"])
        )
        baseline = {
            "converged": True,
            "current_A_per_um": current,
            "current_absolute_relative_error": abs(current - sentaurus_current) / abs(sentaurus_current),
            "state_metrics": state_metrics(Path(point["baseline_state"]), Path(point["export"])),
            "state": str(Path(point["baseline_state"]).resolve()),
        }
        candidate = existing_result(point) if args.report_only else execute(
            point, old_quantum, args.runner.resolve()
        )
        candidate["current_absolute_relative_error"] = abs(
            candidate["current_A_per_um"] - sentaurus_current
        ) / abs(sentaurus_current)
        results.append(
            {
                "name": point["name"],
                "gate_bias_V": point["gate"],
                "drain_bias_V": point["drain"],
                "sentaurus_current_A_per_um": sentaurus_current,
                "current_baseline": baseline,
                "sentaurus_box_insulator": candidate,
            }
        )
    report = {
        "schema": "vela.transportmodels.quantum_contract_ab.v1",
        "as_of": "2026-08-23",
        "status": "complete",
        "points": results,
    }
    REPORT_JSON.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    REPORT_MD.write_text(markdown(report), encoding="utf-8")
    print(json.dumps({"status": "complete", "report": str(REPORT_JSON)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
