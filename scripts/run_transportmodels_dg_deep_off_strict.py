#!/usr/bin/env python3
"""Reclose the three DG Id-Vg deep-off points under hard continuity/KCL gates."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import subprocess
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[1]
REF = REPO / "build-release/reference_tcad/transportmodels_sentaurus2022"
BASE_RUN = REF / "vela_baseline/dg_quantum_contract_regression_2026-08-23/runs/dg"
BASE_CONFIG = BASE_RUN / "03_dg_idvg_curve.json"
SENT_CURVE = REF / "run02/normalized/dg_idvg.csv"
OUTPUT = REF / "reports/transportmodels_dg_deep_off_strict_20260823"
DEFAULT_RUNNER = REPO / "build-release/vela_example_runner.exe"
BIAS_STATES = {
    -1.0: BASE_RUN / "dg_idvg_final_bias_relax_final_state.csv",
    -0.84: BASE_RUN / "dg_idvg_curve_state_bias_m0p840000.csv",
    -0.68: BASE_RUN / "dg_idvg_curve_state_bias_m0p680000.csv",
}
KCL_MARGIN = 10.0


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def bias_tag(value: float) -> str:
    return ("m" if value < 0.0 else "p") + f"{abs(value):.2f}".replace(".", "p")


def runner_environment() -> dict[str, str]:
    env = os.environ.copy()
    prefixes = [r"D:\msys64\ucrt64\bin", r"D:\msys64\usr\bin"]
    env["PATH"] = os.pathsep.join(prefixes + [env.get("PATH", "")])
    return env


def sentaurus_reference() -> dict[float, float]:
    return {
        round(float(row["bias_V"]), 12): abs(float(row["current_total"]))
        for row in read_csv(SENT_CURVE)
    }


def strict_config(bias: float, run_dir: Path, variant: str) -> dict[str, Any]:
    config = json.loads(BASE_CONFIG.read_text(encoding="utf-8"))
    state = BIAS_STATES[bias]
    if not state.is_file():
        raise FileNotFoundError(state)

    config["_comment"] = (
        "DG deep-off reclosure with per-carrier absolute residual diagnostics, "
        "enforced local carrier rows, enforced global continuity closure, and "
        "workflow-hard Id/KCL >= 10 acceptance"
    )
    config["output_csv"] = str((run_dir / "curve.csv").resolve())
    config["log_file"] = str((run_dir / "curve.log").resolve())
    for contact in config["contacts"]:
        if contact["name"] == "gate":
            contact["bias"] = bias

    solver = config["solver"]
    solver.update(
        {
            "max_iter": 200,
            # Disable relative-only termination. The combined absolute norm is
            # intentionally below the deep-off carrier-block scale; Poisson
            # ULP-floor acceptance remains possible only after both carrier
            # row and global continuity gates pass.
            "reltol": 0.0,
            "abstol": 1.0e-18,
            "stall_residual_floor": 2.0e-11,
            "poisson_line_search_stall_residual_floor": 2.0e-11,
            "poisson_line_search_stall_carrier_residual_floor": 1.0e-16,
            "carrier_row_qualified_stall_acceptance": True,
            "line_search": True,
            "damping_factor": 1.0,
            "warm_start": True,
            "diagnostics": True,
            "quasi_fermi_reference": "contact_basin",
            "carrier_row_convergence": {
                "mode": "enforce",
                "eps_row": 1.0e-3,
                "scale_floor": 1.0e-300,
                "min_source_scale_fraction": 0.0,
                # eps_row * min_source_scale is an effective 1e-21 absolute
                # per-row residual cap for source-qualified deep-off rows.
                "min_source_scale": 1.0e-18,
                "min_newton_max_iter": 200,
                "diagnostic_csv": str((run_dir / "carrier_row_violations.csv").resolve()),
                "trace_csv": str((run_dir / "carrier_row_trace.csv").resolve()),
                "trace_first_iterations": 20,
                "trace_every_iterations": 10,
                "recovery": {
                    "mode": "gummel_density",
                    "max_attempts": 2,
                    "max_cycles": 2,
                    "density_change_reltol": 1.0e-10,
                },
            },
            "global_continuity_closure": {
                "mode": "enforce",
                "tolerance": 0.1,
                "source_floor": 1.0e-18,
            },
        }
    )
    if variant == "scaled_filter":
        # Preserve the hard convergence gates while allowing carrier-block
        # progress to drive globalization at the Poisson ULP floor.
        solver["line_search_mode"] = "block_filter"
        solver["residual_filter_gamma"] = 1.0e-4
        solver["residual_filter_envelope_factor"] = 2.0
        solver["continuity_row_scaling"] = {
            "enabled": True,
            "flux_fraction": 1.0e-3,
            "scale_floor": 1.0e-30,
            "min_source_scale": 1.0e-18,
            "min_weight": 1.0e-12,
            "max_weight": 1.0e12,
        }

    sweep = config["sweep"]
    sweep.update(
        {
            "start": bias,
            "stop": bias,
            "step": 0.01,
            "bias_points": [bias],
            "initial_state_file": str(state.resolve()),
            "write_state_file": str((run_dir / "final_state.csv").resolve()),
            "write_state_every_point_prefix": str((run_dir / "state").resolve()),
        }
    )
    diagnostics = sweep.setdefault("diagnostics", {})
    diagnostics["terminal_balance"] = {
        "enabled": True,
        "contacts": ["source", "drain", "gate", "substrate"],
        "csv_file": str((run_dir / "terminal_balance.csv").resolve()),
    }
    diagnostics["srh_balance"] = {
        "enabled": True,
        "material": "Si",
        "drain_contact": "drain",
        "substrate_contact": "substrate",
        "kcl_contacts": ["source", "drain", "gate", "substrate"],
        "resolution_margin_ratio": KCL_MARGIN,
        "csv_file": str((run_dir / "srh_balance.csv").resolve()),
    }
    diagnostics["continuity_balance"] = {
        "enabled": True,
        "contacts": ["source", "drain", "gate", "substrate"],
        "csv_file": str((run_dir / "continuity_balance.csv").resolve()),
    }
    diagnostics["newton_history"] = {
        "enabled": True,
        "csv_file": str((run_dir / "newton_history.csv").resolve()),
        "attempts_csv_file": str((run_dir / "newton_attempts.csv").resolve()),
        "iterations_csv_file": str((run_dir / "newton_iterations.csv").resolve()),
    }
    return config


def last_bias_row(rows: list[dict[str, str]], bias: float) -> dict[str, str] | None:
    matches = [
        row for row in rows
        if "bias_V" in row and math.isclose(float(row["bias_V"]), bias, abs_tol=1.0e-10)
    ]
    return matches[-1] if matches else None


def number(row: dict[str, str] | None, key: str, default: float = math.nan) -> float:
    if row is None or not row.get(key):
        return default
    return float(row[key])


def execute_point(
    runner: Path, bias: float, reference: float, variant: str
) -> dict[str, Any]:
    run_dir = OUTPUT / variant / bias_tag(bias)
    run_dir.mkdir(parents=True, exist_ok=True)
    config = strict_config(bias, run_dir, variant)
    config_path = run_dir / "config.json"
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    completed = subprocess.run(
        [str(runner), "--config", str(config_path)],
        cwd=REPO,
        env=runner_environment(),
        text=True,
        capture_output=True,
        check=False,
    )
    (run_dir / "console.log").write_text(
        completed.stdout + completed.stderr, encoding="utf-8"
    )

    curve = last_bias_row(read_csv(run_dir / "curve.csv"), bias)
    srh = last_bias_row(read_csv(run_dir / "srh_balance.csv"), bias)
    vela = abs(number(curve, "current_total_A_per_um"))
    kcl = abs(number(srh, "four_terminal_kcl_residual_A_per_um"))
    ratio = vela / kcl if math.isfinite(vela) and math.isfinite(kcl) and kcl > 0.0 else math.nan
    solver_converged = curve is not None and curve.get("converged") == "1"
    local_rows_pass = curve is not None and int(curve.get("carrier_row_violations") or -1) == 0
    global_pass = curve is not None and curve.get("global_continuity_closure_satisfied") == "1"
    kcl_pass = math.isfinite(ratio) and ratio >= KCL_MARGIN
    hard_acceptance = (
        completed.returncode == 0
        and solver_converged
        and local_rows_pass
        and global_pass
        and kcl_pass
    )
    return {
        "bias_V": bias,
        "variant": variant,
        "runner_returncode": completed.returncode,
        "solver_converged": solver_converged,
        "convergence_reason": curve.get("newton_convergence_reason", "") if curve else "",
        "failure_reason": curve.get("failure_reason", "missing_curve_row") if curve else "missing_curve_row",
        "vela_A_per_um": vela,
        "sentaurus_A_per_um": reference,
        "absolute_relative_error": abs(vela - reference) / reference if math.isfinite(vela) else math.nan,
        "absolute_log_error_dex": abs(math.log10(vela) - math.log10(reference)) if vela > 0.0 else math.nan,
        "final_psi_residual_norm": number(curve, "final_psi_residual_norm"),
        "final_electron_continuity_residual_norm": number(
            curve, "final_electron_continuity_residual_norm"
        ),
        "final_hole_continuity_residual_norm": number(
            curve, "final_hole_continuity_residual_norm"
        ),
        "carrier_row_violations": int(curve.get("carrier_row_violations") or -1) if curve else -1,
        "carrier_row_max_ratio": number(curve, "carrier_row_max_ratio"),
        "global_continuity_closure_satisfied": global_pass,
        "global_electron_continuity_closure_ratio": number(
            curve, "global_electron_continuity_closure_ratio"
        ),
        "global_hole_continuity_closure_ratio": number(
            curve, "global_hole_continuity_closure_ratio"
        ),
        "four_terminal_kcl_residual_A_per_um": kcl,
        "id_to_kcl_residual_ratio": ratio,
        "kcl_hard_gate_pass": kcl_pass,
        "hard_acceptance": hard_acceptance,
        "config": str(config_path.resolve()),
        "run_directory": str(run_dir.resolve()),
    }


def write_summary(rows: list[dict[str, Any]], runner: Path, variant: str) -> None:
    output = OUTPUT / variant
    output.mkdir(parents=True, exist_ok=True)
    csv_path = output / "deep_off_strict_summary.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    summary = {
        "schema": "vela.transportmodels.dg_deep_off_strict.v1",
        "runner": str(runner.resolve()),
        "variant": variant,
        "policy": {
            "relative_termination_disabled": True,
            "combined_abstol": 1.0e-18,
            "effective_source_qualified_carrier_row_abstol": 1.0e-21,
            "carrier_row_relative_tolerance": 1.0e-3,
            "global_carrier_closure_tolerance": 0.1,
            "kcl_rule": "abs(Id) >= 10 * abs(four-terminal KCL residual)",
        },
        "all_points_hard_accepted": all(row["hard_acceptance"] for row in rows),
        "points": rows,
        "artifacts": {"summary_csv": str(csv_path.resolve())},
    }
    (output / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--runner", type=Path, default=DEFAULT_RUNNER)
    parser.add_argument("--bias", type=float, action="append")
    parser.add_argument(
        "--variant", choices=["strict", "scaled_filter"], default="strict"
    )
    args = parser.parse_args()
    runner = args.runner.resolve()
    if not runner.is_file():
        raise FileNotFoundError(runner)
    selected = tuple(args.bias) if args.bias else tuple(BIAS_STATES)
    references = sentaurus_reference()
    rows = [
        execute_point(runner, bias, references[round(bias, 12)], args.variant)
        for bias in selected
    ]
    write_summary(rows, runner, args.variant)
    return 0 if all(row["hard_acceptance"] for row in rows) else 2


if __name__ == "__main__":
    raise SystemExit(main())
