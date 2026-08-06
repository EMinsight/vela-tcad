#!/usr/bin/env python3
"""Run a resumable 0--7 V E2 branch for the BVmethods NMOS benchmark."""

from __future__ import annotations

import argparse
import copy
import csv
import json
import os
import subprocess
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[1]
RUN_ROOT = REPO / "build-release/reference_tcad/bvmethods_sentaurus2018/run01"
DEFAULT_BASE = (
    RUN_ROOT
    / "vela_validation/btbt_e2_semiconductor_cell_20260803/"
      "activate_A_1ep00/simulation.json"
)
DEFAULT_INITIAL = (
    RUN_ROOT
    / "vela_validation/fermi_dirac_20260802/equilibrium_stages/"
      "04_masetti_high_field/state.csv"
)
DEFAULT_OUTPUT = (
    RUN_ROOT / "vela_validation/btbt_e2_adaptive_0_7_20260804"
)


def absolute(path: Path) -> str:
    return str(path.resolve())


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def bias_token(value: float) -> str:
    return f"{value:.1f}".replace("-", "m").replace(".", "p")


def segment_biases(start: float, stop: float, spacing: float) -> list[float]:
    count = round((stop - start) / spacing)
    return [round(start + index * spacing, 12) for index in range(count + 1)]


def configure_outputs(config: dict[str, Any], run_dir: Path) -> None:
    config["output_csv"] = absolute(run_dir / "sweep.csv")
    sweep = config["sweep"]
    sweep["write_state_file"] = absolute(run_dir / "last_state.csv")
    sweep["write_state_every_point_prefix"] = absolute(
        run_dir / "states" / "accepted_state"
    )
    sweep["write_vtk"] = False
    sweep["vtk_prefix"] = absolute(run_dir / "vtk" / "state")
    diagnostics = sweep.setdefault("diagnostics", {})
    diagnostics["qf_bounds"] = {
        "enabled": True,
        "mode": "reject_and_recover",
        "margin_V": 1.0,
        "min_carrier_density_m3": 1.0e6,
        "csv_file": absolute(run_dir / "qf_bounds.csv"),
    }
    diagnostics["newton_history"] = {
        "enabled": True,
        "csv_file": absolute(run_dir / "newton_history.csv"),
    }
    diagnostics["sg_avalanche_edges"] = {
        "enabled": False,
        "csv_file": absolute(run_dir / "sg_avalanche_edges.csv"),
    }
    diagnostics.pop("terminal_current_method_compare", None)
    diagnostics.pop("continuity_balance", None)
    carrier_rows = config["solver"].setdefault("carrier_row_convergence", {})
    carrier_rows["eps_row"] = 2.0e-3
    carrier_rows["diagnostic_csv"] = absolute(run_dir / "carrier_row_convergence.csv")
    carrier_rows["trace_csv"] = absolute(run_dir / "carrier_row_trace.csv")


def build_segment_config(
    base: dict[str, Any],
    run_dir: Path,
    initial_state: Path,
    biases: list[float],
) -> dict[str, Any]:
    config = copy.deepcopy(base)
    solver = config["solver"]
    solver["method"] = "newton"
    solver["max_iter"] = max(160, int(solver.get("max_iter", 0)))
    solver["reltol"] = 0.0
    # A changed 0.1 V contact bias produces an O(1e4) normalized initial
    # residual, whereas a single valid Newton update reaches O(1e-10).  This
    # threshold therefore prevents stale zero-iteration acceptance while
    # avoiding hundreds of iterations against the drift/diffusion noise floor.
    solver["abstol"] = 5.0e-10
    solver["stall_residual_floor"] = 1.0e-5
    solver["carrier_row_qualified_stall_acceptance"] = True
    solver["quasi_fermi_update_limit_V"] = 0.02
    solver["quasi_fermi_update_limit_minority_V"] = 0.02
    solver["band_to_band"] = {
        "model": "e2",
        "A_cm_inv_s_inv_V_inv2": 3.4e21,
        "B_V_per_cm": 22.6e6,
        "jacobian": "frozen_field",
    }
    solver["impact_ionization"]["coupling_mode"] = "postprocess_only"

    sweep = config["sweep"]
    sweep.update({
        "mode": "bv_reverse",
        "contact": "drain",
        "current_contact": "drain",
        "start": biases[0],
        "stop": biases[-1],
        "step": 0.1,
        "bias_points": biases,
        "initial_step": 0.05,
        "min_step": 1.0e-8,
        "max_step": 0.1,
        "growth_factor": 1.35,
        "shrink_factor": 0.5,
        "max_retries": 29,
        "stop_on_failure": True,
        "initial_state_file": absolute(initial_state),
    })
    sweep["continuation"] = {
        "predictor": {
            "mode": "none",
            "fields": ["psi", "phin", "phip"],
            "max_extrapolation_ratio": 2.0,
        }
    }
    configure_outputs(config, run_dir)
    config["_validation_case"] = {
        "purpose": "resumable full E2 adaptive branch",
        "requested_biases_V": biases,
        "sentaurus_reference_BV_V": 6.377494277837012,
    }
    return config


def completed_segment(run_dir: Path, expected_stop: float) -> bool:
    sweep = run_dir / "sweep.csv"
    state = run_dir / "last_state.csv"
    if not sweep.exists() or not state.exists():
        return False
    records = read_rows(sweep)
    if not records:
        return False
    final = records[-1]
    return (
        final.get("converged", "0") not in {"0", "false", "False"}
        and abs(float(final["bias_V"]) - expected_stop) <= 1.0e-10
    )


def run_segment(runner: Path, config_path: Path, log_path: Path) -> None:
    completed = subprocess.run(
        [str(runner), "--config", str(config_path)],
        cwd=config_path.parent,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    log_path.write_text(completed.stdout, encoding="utf-8")
    if completed.returncode != 0:
        tail = "\n".join(completed.stdout.splitlines()[-40:])
        raise RuntimeError(f"segment failed: {config_path}\n{tail}")


def write_branch_summary(output: Path, segment_dirs: list[Path]) -> None:
    records: list[dict[str, str]] = []
    seen: set[float] = set()
    for run_dir in segment_dirs:
        for row in read_rows(run_dir / "sweep.csv"):
            bias = round(float(row["bias_V"]), 12)
            if bias in seen:
                continue
            seen.add(bias)
            records.append(row)
    records.sort(key=lambda row: float(row["bias_V"]))
    if not records:
        return
    path = output / "branch_0_7.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base", type=Path, default=DEFAULT_BASE)
    parser.add_argument("--initial-state", type=Path, default=DEFAULT_INITIAL)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--start", type=float, default=0.0)
    parser.add_argument("--stop", type=float, default=7.0)
    parser.add_argument("--segment-span", type=float, default=0.5)
    parser.add_argument("--record-spacing", type=float, default=0.1)
    parser.add_argument("--runner", type=Path)
    parser.add_argument("--no-resume", action="store_true")
    args = parser.parse_args()

    if args.runner is None:
        executable = "vela_example_runner.exe" if os.name == "nt" else "vela_example_runner"
        args.runner = REPO / "build-release" / executable
    runner = args.runner.resolve()
    base = json.loads(args.base.read_text(encoding="utf-8-sig"))
    args.output.mkdir(parents=True, exist_ok=True)

    current = args.start
    initial_state = args.initial_state.resolve()
    segment_dirs: list[Path] = []
    while current < args.stop - 1.0e-12:
        segment_stop = min(args.stop, current + args.segment_span)
        biases = segment_biases(current, segment_stop, args.record_spacing)
        # The previous segment already solved and persisted ``current``.  Do
        # not force the strict Newton acceptance logic to resolve that exact
        # same bias before advancing; start at the first new requested point.
        if current > args.start + 1.0e-12:
            biases = biases[1:]
        run_dir = args.output / (
            f"segment_{bias_token(current)}_{bias_token(segment_stop)}"
        )
        segment_dirs.append(run_dir)
        if not args.no_resume and completed_segment(run_dir, segment_stop):
            print(f"{run_dir.name}: resume-skip", flush=True)
            initial_state = (run_dir / "last_state.csv").resolve()
            current = segment_stop
            continue

        config = build_segment_config(base, run_dir, initial_state, biases)
        run_dir.mkdir(parents=True, exist_ok=True)
        config_path = run_dir / "simulation.json"
        config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
        print(f"{run_dir.name}: running", flush=True)
        run_segment(runner, config_path.resolve(), run_dir / "run.log")
        if not completed_segment(run_dir, segment_stop):
            raise RuntimeError(f"segment did not reach {segment_stop}: {run_dir}")
        final = read_rows(run_dir / "sweep.csv")[-1]
        print(
            f"{run_dir.name}: converged Id={float(final['current_total_A_per_um']):.6e} A/um",
            flush=True,
        )
        initial_state = (run_dir / "last_state.csv").resolve()
        current = segment_stop
        write_branch_summary(args.output, segment_dirs)

    write_branch_summary(args.output, segment_dirs)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
