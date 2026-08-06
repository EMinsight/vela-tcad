#!/usr/bin/env python3
"""Run staged Sentaurus-compatible E2 activation and the requested bias probes."""

from __future__ import annotations

import argparse
import copy
import csv
import json
import subprocess
from pathlib import Path


DEFAULT_BASE = Path(
    "build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela_validation/"
    "iic_rebuild_fd_gummel_20260803/probe_6p4_full/postprocess_only/simulation.json"
)
DEFAULT_OUTPUT = Path(
    "build-release/reference_tcad/bvmethods_sentaurus2018/run01/vela_validation/"
    "btbt_e2_self_consistent_20260803"
)


def absolute(path: Path) -> str:
    return str(path.resolve())


def configure_outputs(config: dict, run_dir: Path, write_vtk: bool) -> None:
    config["output_csv"] = absolute(run_dir / "sweep.csv")
    sweep = config["sweep"]
    sweep["write_vtk"] = write_vtk
    sweep["vtk_prefix"] = absolute(run_dir / "vtk" / "state")
    sweep["write_state_file"] = absolute(run_dir / "last_state.csv")
    sweep["write_state_every_point_prefix"] = absolute(
        run_dir / "states" / "accepted_state"
    )
    diagnostics = sweep.get("diagnostics", {})
    for name, entry in diagnostics.items():
        if not isinstance(entry, dict):
            continue
        if name == "sg_avalanche_edges":
            entry["enabled"] = False
        if "csv_file" in entry:
            entry["csv_file"] = absolute(run_dir / f"{name}.csv")
    carrier_rows = config["solver"].get("carrier_row_convergence", {})
    if "diagnostic_csv" in carrier_rows:
        carrier_rows["diagnostic_csv"] = absolute(
            run_dir / "carrier_row_convergence.csv"
        )
    if "trace_csv" in carrier_rows:
        carrier_rows["trace_csv"] = absolute(run_dir / "carrier_row_trace.csv")


def write_config(config: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")


def run_case(runner: Path, config_path: Path, log_path: Path) -> None:
    command = [str(runner.resolve()), "--config", str(config_path.resolve())]
    completed = subprocess.run(
        command, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        check=False,
    )
    log_path.write_text(completed.stdout, encoding="utf-8")
    print(f"{config_path.parent.name}: exit={completed.returncode}")
    if completed.returncode != 0:
        tail = "\n".join(completed.stdout.splitlines()[-30:])
        raise RuntimeError(f"simulation failed: {config_path}\n{tail}")


def last_state_from_sweep(path: Path) -> Path:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise RuntimeError(f"empty sweep: {path}")
    row = rows[-1]
    if row.get("converged", "1") in {"0", "false", "False"}:
        raise RuntimeError(f"last sweep point did not converge: {path}")
    return path.parent / "last_state.csv"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", type=Path, default=DEFAULT_BASE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--runner", type=Path,
                        default=Path("build-release/vela_example_runner.exe"))
    parser.add_argument("--scales", nargs="+", type=float,
                        default=[1.0e-3, 1.0e-2, 1.0e-1, 1.0])
    parser.add_argument(
        "--initial-state",
        type=Path,
        help="optional restart state; defaults to the state referenced by --base",
    )
    parser.add_argument("--bias", type=float, default=6.4)
    parser.add_argument(
        "--activation-only",
        action="store_true",
        help="stop after the staged fixed-bias E2 activation",
    )
    args = parser.parse_args()

    base = json.loads(args.base.read_text(encoding="utf-8"))
    initial_state = (
        args.initial_state
        if args.initial_state is not None
        else Path(base["sweep"]["initial_state_file"])
    )
    args.output.mkdir(parents=True, exist_ok=True)

    for scale in args.scales:
        label = f"activate_A_{scale:.0e}".replace("+", "p").replace("-", "m")
        run_dir = args.output / label
        config = copy.deepcopy(base)
        config["solver"]["band_to_band"] = {
            "model": "e2",
            "A_cm_inv_s_inv_V_inv2": 3.4e21 * scale,
            "B_V_per_cm": 22.6e6,
            # The residual is recomputed from the current electric field on every
            # Newton iteration.  Freezing only dG/dpsi in the linear solve keeps
            # this production-size branch affordable; the exact finite-difference
            # option remains available and is covered by a focused unit test.
            "jacobian": "frozen_field",
            "jacobian_relative_step": 1.0e-7,
        }
        config["solver"]["max_iter"] = max(160, config["solver"].get("max_iter", 0))
        config["solver"]["quasi_fermi_update_limit_V"] = 0.02
        config["solver"]["quasi_fermi_update_limit_minority_V"] = 0.02
        config["sweep"]["start"] = args.bias
        config["sweep"]["stop"] = args.bias
        config["sweep"]["bias_points"] = [args.bias]
        config["sweep"]["initial_state_file"] = absolute(initial_state)
        configure_outputs(config, run_dir, write_vtk=(scale == args.scales[-1]))
        config_path = run_dir / "simulation.json"
        write_config(config, config_path)
        run_case(args.runner, config_path, run_dir / "run.log")
        initial_state = last_state_from_sweep(run_dir / "sweep.csv")

    if args.activation_only:
        return 0

    validation_dir = args.output / "bias_validation"
    config = copy.deepcopy(base)
    config["solver"]["band_to_band"] = {
        "model": "e2",
        "A_cm_inv_s_inv_V_inv2": 3.4e21,
        "B_V_per_cm": 22.6e6,
        "jacobian": "frozen_field",
        "jacobian_relative_step": 1.0e-7,
    }
    config["solver"]["max_iter"] = max(160, config["solver"].get("max_iter", 0))
    config["solver"]["quasi_fermi_update_limit_V"] = 0.02
    config["solver"]["quasi_fermi_update_limit_minority_V"] = 0.02
    config["sweep"]["start"] = 6.4
    config["sweep"]["stop"] = 2.0
    config["sweep"]["step"] = -0.1
    config["sweep"]["initial_step"] = 0.1
    config["sweep"]["bias_points"] = [6.4, 6.0, 5.0, 4.0, 2.0]
    config["sweep"]["initial_state_file"] = absolute(initial_state)
    configure_outputs(config, validation_dir, write_vtk=True)
    config_path = validation_dir / "simulation.json"
    write_config(config, config_path)
    run_case(args.runner, config_path, validation_dir / "run.log")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
