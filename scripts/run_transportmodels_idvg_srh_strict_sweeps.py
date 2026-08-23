#!/usr/bin/env python3
"""Run strict DD forward/reverse sweeps and a forward DG regression.

The workflow is intentionally limited to the TransportModels Id-Vg benchmark.
It enables the silicon SRH/terminal audit added for the deep-off investigation
and preserves one reproducible configuration per sweep direction.
"""

from __future__ import annotations

import argparse
import copy
import json
import os
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
BASELINE = (
    ROOT
    / "build-release/reference_tcad/transportmodels_sentaurus2022/vela_baseline"
)
DD_SOURCE = BASELINE / "dd_phase7_shared_baseline_2026-08-21/idvg/config.json"
DG_SOURCE = BASELINE / "dg_post_p2_regression_v4_2026-08-21/idvg/config.json"
DD_INITIAL = (
    BASELINE
    / "dd_phase7_shared_baseline_2026-08-21/idvg/state_bias_m1p000000.csv"
)
DG_INITIAL = (
    BASELINE
    / "dg_post_p2_regression_v4_2026-08-21/idvg/state_bias_m1p000000.csv"
)
OUTPUT_ROOT = BASELINE / "idvg_srh_strict_2026-08-21"
RUNNER = ROOT / "build-release/vela_example_runner.exe"

FORWARD_POINTS = [-1.0 + 0.16 * i for i in range(21)]
FORWARD_POINTS[-1] = 2.2
REVERSE_POINTS = list(reversed(FORWARD_POINTS[:-1]))
ROUNDTRIP_POINTS = FORWARD_POINTS + REVERSE_POINTS
DG_RESUME_POINTS = [-0.50, -0.48, -0.44, -0.40] + FORWARD_POINTS[4:]


def load_config(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def runner_environment() -> dict[str, str]:
    env = os.environ.copy()
    ucrt = r"D:\msys64\ucrt64\bin"
    usr = r"D:\msys64\usr\bin"
    env["PATH"] = os.pathsep.join([ucrt, usr, env.get("PATH", "")])
    return env


def strict_solver(config: dict, quantum: bool) -> None:
    solver = config["solver"]
    solver["max_iter"] = 200
    solver["reltol"] = 1.0e-9
    solver["abstol"] = 1.0e-14
    # Calibration on DD/DG deep-off checkpoints places the converged numerical
    # plateau near 1.5e-11--1.9e-11.  Keep the floor just above that plateau;
    # carrier-row and global-continuity hard gates must also pass.
    solver["stall_residual_floor"] = 2.0e-11
    # Smaller quasi-Fermi caps did not improve the hard-gated checkpoints, and
    # 0.5 initial damping greatly increased DD/DG iteration counts.  Retain a
    # full Newton proposal with backtracking globalization.
    solver["line_search"] = True
    solver["damping_factor"] = 1.0
    solver["quasi_fermi_update_limit_V"] = 2.5e-2
    solver["warm_start"] = True
    # Store quasi-Fermi unknowns relative to their majority-carrier contact so
    # deep-off current-carrying increments remain representable.  Mixed-
    # material gauge rows and node-source Jacobians are reference invariant.
    solver["quasi_fermi_reference"] = "contact_majority"
    solver["carrier_row_qualified_stall_acceptance"] = True
    solver["carrier_row_convergence"] = {
        "mode": "enforce",
        "eps_row": 1.0e-3,
        "scale_floor": 1.0e-300,
        "min_source_scale_fraction": 0.0,
        "min_source_scale": 1.0e-18,
        "min_newton_max_iter": 200,
        "recovery": {
            "mode": "gummel_density",
            "max_attempts": 2,
            "max_cycles": 2,
            "density_change_reltol": 1.0e-10,
        },
    }
    # A deep-off point is valid only when each carrier's contact flux closes
    # its integrated source within one order of magnitude.  Enforcing this in
    # Newton prevents an unresolved point from entering the comparison curve.
    solver["global_continuity_closure"] = {
        "mode": "enforce",
        "tolerance": 0.1,
        "source_floor": 1.0e-18,
    }
    if not quantum:
        solver.setdefault("electron_quantum_potential", {})["enabled"] = False


def make_config(
    source: Path,
    output_dir: Path,
    points: list[float],
    initial_state: Path,
    quantum: bool,
) -> dict:
    config = copy.deepcopy(load_config(source))
    output_dir.mkdir(parents=True, exist_ok=True)
    config["_comment"] = (
        "TransportModels Id-Vg strict SRH/KCL validation; "
        f"{'DG' if quantum else 'DD'} {points[0]:g} V to {points[-1]:g} V"
    )
    config["output_csv"] = str((output_dir / "curve.csv").resolve())
    config["log_file"] = str((output_dir / "curve.log").resolve())
    strict_solver(config, quantum)

    config["contacts"] = copy.deepcopy(config["contacts"])
    for contact in config["contacts"]:
        if contact["name"] == "gate":
            contact["bias"] = points[0]

    sweep = config["sweep"]
    sweep["start"] = points[0]
    sweep["stop"] = points[-1]
    sweep["step"] = points[1] - points[0]
    sweep["bias_points"] = points
    sweep["initial_state_file"] = str(initial_state.resolve())
    sweep["write_state_file"] = str((output_dir / "final_state.csv").resolve())
    sweep["write_state_every_point_prefix"] = str((output_dir / "state").resolve())
    diagnostics = sweep.setdefault("diagnostics", {})
    diagnostics["terminal_balance"] = {
        "enabled": True,
        "contacts": ["source", "drain", "gate", "substrate"],
        "csv_file": str((output_dir / "terminal_balance.csv").resolve()),
    }
    diagnostics["srh_balance"] = {
        "enabled": True,
        "material": "Si",
        "drain_contact": "drain",
        "substrate_contact": "substrate",
        "kcl_contacts": ["source", "drain", "gate", "substrate"],
        "resolution_margin_ratio": 10.0,
        "csv_file": str((output_dir / "srh_balance.csv").resolve()),
    }
    return config


def execute(name: str, config: dict, output_dir: Path) -> dict:
    config_path = output_dir / "config.json"
    config_path.write_text(json.dumps(config, indent=2), encoding="utf-8")
    completed = subprocess.run(
        [str(RUNNER), "--config", str(config_path)],
        cwd=ROOT,
        env=runner_environment(),
        text=True,
        capture_output=True,
        check=False,
    )
    (output_dir / "console.log").write_text(
        completed.stdout + completed.stderr, encoding="utf-8"
    )
    curve = output_dir / "curve.csv"
    return {
        "name": name,
        "returncode": completed.returncode,
        "config": str(config_path.resolve()),
        "curve": str(curve.resolve()),
        "srh_balance": str((output_dir / "srh_balance.csv").resolve()),
        "curve_exists": curve.exists(),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--phase",
        choices=["all", "dd", "dd_reverse", "dd_roundtrip", "dg", "dg_resume"],
        default="all",
        help="Run DD forward/reverse, DG forward, or all in sequence.",
    )
    args = parser.parse_args()
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    runs: list[dict] = []

    if args.phase in {"all", "dd"}:
        dd_forward_dir = OUTPUT_ROOT / "dd_forward"
        dd_forward = make_config(
            DD_SOURCE, dd_forward_dir, FORWARD_POINTS, DD_INITIAL, quantum=False
        )
        forward_result = execute("dd_forward", dd_forward, dd_forward_dir)
        runs.append(forward_result)
        if forward_result["returncode"] == 0:
            dd_reverse_dir = OUTPUT_ROOT / "dd_reverse"
            dd_reverse = make_config(
                DD_SOURCE,
                dd_reverse_dir,
                REVERSE_POINTS,
                dd_forward_dir / "final_state.csv",
                quantum=False,
            )
            runs.append(execute("dd_reverse", dd_reverse, dd_reverse_dir))

    if args.phase == "dd_reverse":
        dd_forward_dir = OUTPUT_ROOT / "dd_forward"
        dd_reverse_dir = OUTPUT_ROOT / "dd_reverse"
        dd_reverse = make_config(
            DD_SOURCE,
            dd_reverse_dir,
            REVERSE_POINTS,
            dd_forward_dir / "final_state.csv",
            quantum=False,
        )
        runs.append(execute("dd_reverse", dd_reverse, dd_reverse_dir))

    if args.phase == "dd_roundtrip":
        dd_roundtrip_dir = OUTPUT_ROOT / "dd_roundtrip"
        dd_roundtrip = make_config(
            DD_SOURCE,
            dd_roundtrip_dir,
            ROUNDTRIP_POINTS,
            DD_INITIAL,
            quantum=False,
        )
        runs.append(execute("dd_roundtrip", dd_roundtrip, dd_roundtrip_dir))

    if args.phase in {"all", "dg"}:
        dg_forward_dir = OUTPUT_ROOT / "dg_forward"
        dg_forward = make_config(
            DG_SOURCE, dg_forward_dir, FORWARD_POINTS, DG_INITIAL, quantum=True
        )
        runs.append(execute("dg_forward", dg_forward, dg_forward_dir))

    if args.phase == "dg_resume":
        dg_forward_dir = OUTPUT_ROOT / "dg_forward"
        dg_resume_dir = OUTPUT_ROOT / "dg_resume"
        dg_resume = make_config(
            DG_SOURCE,
            dg_resume_dir,
            DG_RESUME_POINTS,
            dg_forward_dir / "state_bias_m0p520000.csv",
            quantum=True,
        )
        runs.append(execute("dg_resume", dg_resume, dg_resume_dir))

    summary = {
        "schema": "vela.transportmodels.idvg_srh_strict.v1",
        "output_root": str(OUTPUT_ROOT.resolve()),
        "runs": runs,
    }
    (OUTPUT_ROOT / "execution_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))
    return 0 if runs and all(run["returncode"] == 0 for run in runs) else 1


if __name__ == "__main__":
    raise SystemExit(main())
