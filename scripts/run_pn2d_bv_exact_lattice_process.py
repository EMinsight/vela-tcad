#!/usr/bin/env python3
"""Run Vela PN2D BV branches on a Sentaurus-defined exact bias lattice."""

from __future__ import annotations

import argparse
import copy
import csv
import json
import math
import subprocess
from pathlib import Path
from typing import Any


BRANCHES = ("avalanche_off", "iic_postprocess", "avalanche_on")
EXACT_BIAS_TOLERANCE_V = 1.0e-10


def parse_branch_list(raw: str) -> tuple[str, ...]:
    branches = tuple(item.strip() for item in raw.split(",") if item.strip())
    if not branches:
        raise ValueError("--branches must select at least one branch")
    unknown = sorted(set(branches) - set(BRANCHES))
    if unknown:
        raise ValueError(f"unknown branches: {', '.join(unknown)}")
    if len(set(branches)) != len(branches):
        raise ValueError("--branches contains a duplicate")
    return branches


def exact_bias_lattice(manifest: dict[str, Any]) -> list[float]:
    requested: dict[str, list[float]] = {}
    for branch in manifest.get("branch_records", []):
        name = str(branch.get("branch", ""))
        if name in BRANCHES:
            requested[name] = [float(value) for value in branch["requested_biases_V"]]
    missing = [branch for branch in BRANCHES if branch not in requested]
    if missing:
        raise ValueError(f"Sentaurus manifest is missing branches: {', '.join(missing)}")
    reference = requested[BRANCHES[0]]
    for branch in BRANCHES[1:]:
        if requested[branch] != reference:
            raise ValueError(f"Sentaurus branch lattice differs for {branch}")
    if not reference:
        raise ValueError("Sentaurus exact bias lattice is empty")
    if len(set(reference)) != len(reference):
        raise ValueError("Sentaurus exact bias lattice contains duplicates")
    return reference


def branch_config(
    base: dict[str, Any],
    branch: str,
    biases: list[float],
    case_dir: Path,
    max_iter: int,
    qf_carrier_truncation: float | None = None,
) -> dict[str, Any]:
    config = copy.deepcopy(base)
    solver = config["solver"]
    solver["max_iter"] = max_iter
    solver["verbose"] = False
    solver["diagnostics"] = True
    handoff = solver.setdefault("handoff", {})
    handoff["newton_max_iter"] = max_iter

    impact = copy.deepcopy(solver.get("impact_ionization", {}))
    if branch == "avalanche_off":
        solver["impact_ionization"] = {"model": "none"}
    else:
        impact["model"] = "van_overstraeten"
        impact["coupling_mode"] = (
            "postprocess_only" if branch == "iic_postprocess"
            else "self_consistent"
        )
        if qf_carrier_truncation is not None:
            impact["quasi_fermi_carrier_truncation"] = qf_carrier_truncation
        solver["impact_ionization"] = impact

    case_dir = case_dir.resolve()
    states_dir = case_dir / "states"
    states_dir.mkdir(parents=True, exist_ok=True)
    output_csv = case_dir / "iv.csv"
    attempts_csv = case_dir / "newton_attempts.csv"
    iterations_csv = case_dir / "newton_iterations.csv"
    legacy_newton_csv = case_dir / "newton_history.csv"
    diagnostics: dict[str, Any] = {
        "newton_history": {
            "enabled": True,
            "csv_file": str(legacy_newton_csv),
            "attempts_csv_file": str(attempts_csv),
            "iterations_csv_file": str(iterations_csv),
        },
    }
    if branch != "avalanche_off":
        diagnostics["bv_process_probe"] = {
            "enabled": True,
            "csv_file": str(case_dir / "process_probe.csv"),
        }

    sweep = config["sweep"]
    sweep.update(
        {
            "start": biases[0],
            "stop": biases[-1],
            "step": -0.05 if biases[-1] < biases[0] else 0.05,
            "bias_points": biases,
            "initial_step": 0.05,
            "min_step": 1.0e-10,
            "max_step": 0.05,
            "growth_factor": 1.2,
            "shrink_factor": 0.5,
            "max_retries": 30,
            "stop_on_failure": True,
            "write_vtk": False,
            "write_state_file": str(case_dir / "last_state.csv"),
            "write_state_every_point_prefix": str(states_dir / "state"),
            "diagnostics": diagnostics,
        }
    )
    sweep.pop("initial_state_file", None)
    sweep.setdefault("initialization", {"mode": "poisson_block"})
    config["output_csv"] = str(output_csv)
    config.pop("output_vtk_prefix", None)
    return config


def read_sweep_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def qualify_rows(rows: list[dict[str, str]], biases: list[float]) -> dict[str, Any]:
    observed = [float(row["bias_V"]) for row in rows]
    converged = [row["converged"] == "1" for row in rows]
    exact = (
        len(observed) == len(biases)
        and all(
            abs(actual - requested) <= EXACT_BIAS_TOLERANCE_V
            for actual, requested in zip(observed, biases)
        )
    )
    complete = exact and all(converged)
    first_failure = next(
        (
            {
                "requested_bias_V": biases[index] if index < len(biases) else None,
                "actual_bias_V": observed[index],
                "failure_reason": rows[index].get("failure_reason", ""),
                "newton_failure_class": rows[index].get(
                    "newton_failure_class", ""
                ),
            }
            for index, ok in enumerate(converged)
            if not ok
        ),
        None,
    )
    if first_failure is None and len(observed) < len(biases):
        first_failure = {
            "requested_bias_V": biases[len(observed)],
            "actual_bias_V": None,
            "failure_reason": "missing_bias_row",
            "newton_failure_class": "",
        }
    return {
        "complete_exact_lattice": complete,
        "requested_bias_count": len(biases),
        "observed_bias_count": len(observed),
        "all_converged": all(converged) and len(converged) == len(biases),
        "all_exact": exact,
        "first_failure": first_failure,
        "last_observed_bias_V": observed[-1] if observed else None,
    }


def run_branch(
    runner: Path,
    base: dict[str, Any],
    branch: str,
    biases: list[float],
    output_root: Path,
    max_iter: int,
    qf_carrier_truncation: float | None,
    resume: bool,
) -> dict[str, Any]:
    case_dir = output_root / branch
    case_dir.mkdir(parents=True, exist_ok=True)
    config = branch_config(
        base,
        branch,
        biases,
        case_dir,
        max_iter,
        qf_carrier_truncation,
    )
    config_path = case_dir / "simulation.json"
    config_path.write_text(
        json.dumps(config, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    log_path = case_dir / "runner.log"
    output_csv = case_dir / "iv.csv"
    resumed = False
    returncode = 0
    if resume and output_csv.is_file():
        existing = qualify_rows(read_sweep_rows(output_csv), biases)
        resumed = existing["complete_exact_lattice"]
    if not resumed:
        with log_path.open("w", encoding="utf-8") as log:
            completed = subprocess.run(
                [str(runner), "--config", str(config_path)],
                cwd=case_dir,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=False,
            )
        returncode = completed.returncode
    if not output_csv.is_file():
        qualification = {
            "complete_exact_lattice": False,
            "requested_bias_count": len(biases),
            "observed_bias_count": 0,
            "all_converged": False,
            "all_exact": False,
            "first_failure": {
                "requested_bias_V": biases[0],
                "actual_bias_V": None,
                "failure_reason": "missing_sweep_output",
                "newton_failure_class": "",
            },
            "last_observed_bias_V": None,
        }
    else:
        qualification = qualify_rows(read_sweep_rows(output_csv), biases)
    return {
        "branch": branch,
        "returncode": returncode,
        "resumed": resumed,
        "config": str(config_path),
        "output_csv": str(output_csv),
        "runner_log": str(log_path),
        "max_iter": max_iter,
        **qualification,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--base-config", type=Path, required=True)
    parser.add_argument("--sentaurus-manifest", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--branches", default=",".join(BRANCHES))
    parser.add_argument("--max-iter", type=int, default=80)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument(
        "--qf-carrier-truncation",
        type=float,
        help=(
            "opt-in low-density n,p floor relative to ni used only when "
            "rebuilding the avalanche quasi-Fermi driving field"
        ),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.max_iter <= 0:
        raise ValueError("--max-iter must be positive")
    if (
        args.qf_carrier_truncation is not None
        and (
            not math.isfinite(args.qf_carrier_truncation)
            or args.qf_carrier_truncation < 0.0
        )
    ):
        raise ValueError("--qf-carrier-truncation must be finite and nonnegative")
    runner = args.runner.resolve()
    base_path = args.base_config.resolve()
    sentaurus_path = args.sentaurus_manifest.resolve()
    output_root = args.output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    base = json.loads(base_path.read_text(encoding="utf-8-sig"))
    sentaurus = json.loads(sentaurus_path.read_text(encoding="utf-8"))
    biases = exact_bias_lattice(sentaurus)
    selected = parse_branch_list(args.branches)
    results = [
        run_branch(
            runner,
            base,
            branch,
            biases,
            output_root,
            args.max_iter,
            args.qf_carrier_truncation,
            args.resume,
        )
        for branch in selected
    ]
    complete = all(
        result["returncode"] == 0 and result["complete_exact_lattice"]
        for result in results
    )
    execution = {
        "schema": "vela.pn2d_bv_exact_lattice_execution.v1",
        "status": "passed" if complete else "failed",
        "outcome": (
            "complete_exact_lattice"
            if complete
            else "incomplete_exact_lattice"
        ),
        "runner": str(runner),
        "base_config": str(base_path),
        "sentaurus_manifest": str(sentaurus_path),
        "max_iter": args.max_iter,
        "candidate": {
            "axis": "impact_ionization.quasi_fermi_carrier_truncation",
            "value": args.qf_carrier_truncation,
            "default_unchanged": True,
        },
        "requested_biases_V": biases,
        "branches": results,
    }
    execution_path = output_root / "execution.json"
    execution_path.write_text(
        json.dumps(execution, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(execution, indent=2, sort_keys=True))
    return 0 if complete else 2


if __name__ == "__main__":
    raise SystemExit(main())
