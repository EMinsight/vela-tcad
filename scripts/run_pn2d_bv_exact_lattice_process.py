#!/usr/bin/env python3
"""Run Vela PN2D BV branches on a Sentaurus-defined exact bias lattice."""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import math
import subprocess
from pathlib import Path
from typing import Any


BRANCHES = ("avalanche_off", "iic_postprocess", "avalanche_on")
EXACT_BIAS_TOLERANCE_V = 1.0e-10
CONTINUATION_SCHEDULES = {
    "standard_0p05": {
        "initial_step_V": 0.05,
        "minimum_step_V": 1.0e-10,
        "maximum_step_V": 0.05,
        "growth_factor": 1.2,
        "shrink_factor": 0.5,
    },
    "refined_0p025": {
        "initial_step_V": 0.025,
        "minimum_step_V": 1.0e-10,
        "maximum_step_V": 0.025,
        "growth_factor": 1.2,
        "shrink_factor": 0.5,
    },
}
SCHEDULE_KEYS = {
    "step",
    "initial_step",
    "min_step",
    "max_step",
    "growth_factor",
    "shrink_factor",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def payload_sha256(payload: Any) -> str:
    encoded = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def bias_token(bias: float) -> str:
    sign = "m" if bias < 0.0 else ""
    return f"{sign}{abs(bias):.6f}".replace(".", "p")


def normalized_non_schedule_config(config: dict[str, Any]) -> dict[str, Any]:
    normalized = copy.deepcopy(config)
    normalized.pop("output_csv", None)
    normalized.pop("output_vtk_prefix", None)
    sweep = normalized.get("sweep", {})
    for key in SCHEDULE_KEYS:
        sweep.pop(key, None)
    for key in (
        "write_state_file",
        "write_state_every_point_prefix",
        "vtk_prefix",
        "csv_file",
    ):
        sweep.pop(key, None)
    diagnostics = sweep.get("diagnostics", {})
    for value in diagnostics.values():
        if isinstance(value, dict):
            for key in tuple(value):
                if key.endswith("_file") or key == "csv_file":
                    value.pop(key, None)
    return normalized


def physics_config(config: dict[str, Any]) -> dict[str, Any]:
    normalized = normalized_non_schedule_config(config)
    normalized.pop("sweep", None)
    solver = normalized.get("solver", {})
    for key in ("max_iter", "verbose", "diagnostics", "handoff"):
        solver.pop(key, None)
    return normalized


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


def apply_physical_input_overrides(
    base: dict[str, Any],
    mesh_file: Path | None,
    doping_file: Path | None,
) -> tuple[dict[str, Any], dict[str, dict[str, str]]]:
    config = copy.deepcopy(base)
    records: dict[str, dict[str, str]] = {}
    for key, path in (
        ("mesh_file", mesh_file),
        ("node_doping_file", doping_file),
    ):
        if path is None:
            continue
        resolved = path.resolve()
        if not resolved.is_file():
            raise FileNotFoundError(f"missing physical input override: {resolved}")
        config[key] = str(resolved)
        records[key] = {
            "path": str(resolved),
            "sha256": sha256(resolved),
        }
    return config, records


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
    sg_laux_candidate: bool = False,
    continuation_schedule: str = "standard_0p05",
    terminal_current_method_compare: bool = False,
    newton_reltol: float | None = None,
    newton_abstol: float | None = None,
) -> dict[str, Any]:
    if continuation_schedule not in CONTINUATION_SCHEDULES:
        raise ValueError(
            f"unknown continuation schedule: {continuation_schedule}"
        )
    schedule = CONTINUATION_SCHEDULES[continuation_schedule]
    config = copy.deepcopy(base)
    solver = config["solver"]
    solver["max_iter"] = max_iter
    solver["verbose"] = False
    solver["diagnostics"] = True
    if newton_reltol is not None:
        solver["reltol"] = newton_reltol
    if newton_abstol is not None:
        solver["abstol"] = newton_abstol
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
        if sg_laux_candidate:
            impact["current_approximation"] = "element_edge_sg_gss_laux"
            impact["source_mapping_mode"] = "element_vertex_box_measure"
            # The archived baseline selects gss_logistic only for the retired
            # triangle source.  Once that source is replaced, normalize this
            # inactive compatibility field without introducing another
            # physical candidate axis.
            impact["cell_reconstructed_midpoint_density"] = "bernoulli"
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
    if terminal_current_method_compare:
        diagnostics["terminal_current_method_compare"] = {
            "enabled": True,
            "csv_file": str(case_dir / "terminal_current_method_compare.csv"),
        }

    sweep = config["sweep"]
    sweep.update(
        {
            "start": biases[0],
            "stop": biases[-1],
            "step": (
                -schedule["maximum_step_V"]
                if biases[-1] < biases[0]
                else schedule["maximum_step_V"]
            ),
            "bias_points": biases,
            "initial_step": schedule["initial_step_V"],
            "min_step": schedule["minimum_step_V"],
            "max_step": schedule["maximum_step_V"],
            "growth_factor": schedule["growth_factor"],
            "shrink_factor": schedule["shrink_factor"],
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
    sg_laux_candidate: bool,
    continuation_schedule: str,
    terminal_current_method_compare: bool,
    newton_reltol: float | None,
    newton_abstol: float | None,
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
        sg_laux_candidate,
        continuation_schedule,
        terminal_current_method_compare,
        newton_reltol,
        newton_abstol,
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
    state_files: dict[str, dict[str, str]] = {}
    states_dir = case_dir / "states"
    for bias in biases:
        state_path = states_dir / f"state_bias_{bias_token(bias)}.csv"
        if state_path.is_file():
            state_files[f"{bias:.17g}"] = {
                "path": str(state_path),
                "sha256": sha256(state_path),
            }
    return {
        "branch": branch,
        "returncode": returncode,
        "resumed": resumed,
        "config": str(config_path),
        "output_csv": str(output_csv),
        "runner_log": str(log_path),
        "max_iter": max_iter,
        "continuation_schedule": continuation_schedule,
        "continuation_schedule_parameters": CONTINUATION_SCHEDULES[
            continuation_schedule
        ],
        "physics_config_sha256": payload_sha256(physics_config(config)),
        "non_schedule_config_sha256": payload_sha256(
            normalized_non_schedule_config(config)
        ),
        "config_sha256": sha256(config_path),
        "output_csv_sha256": (
            sha256(output_csv) if output_csv.is_file() else None
        ),
        "state_files": state_files,
        **qualification,
    }


def build_state_manifest(
    results: list[dict[str, Any]],
    biases: list[float],
    output_root: Path,
) -> dict[str, Any]:
    branches: list[dict[str, Any]] = []
    for result in results:
        states = result["state_files"]
        records: list[dict[str, Any]] = []
        for bias in biases:
            matches = [
                record
                for raw_bias, record in states.items()
                if abs(float(raw_bias) - bias) <= EXACT_BIAS_TOLERANCE_V
            ]
            if len(matches) != 1:
                raise ValueError(
                    f"{result['branch']} {bias:g} V: missing exact state file"
                )
            state = matches[0]
            state_path = Path(state["path"]).resolve()
            records.append(
                {
                    "requested_bias_V": bias,
                    "actual_bias_V": bias,
                    "snapshot_tdr": {
                        "path": str(state_path.relative_to(output_root)),
                        "sha256": state["sha256"],
                    },
                }
            )
        branches.append(
            {
                "branch": result["branch"],
                "requested_biases_V": biases,
                "bias_records": records,
            }
        )
    return {
        "schema": "vela.pn2d_bv_exact_state_manifest.v1",
        "status": "passed",
        "outcome": "complete_exact_state_manifest",
        "branch_records": branches,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--base-config", type=Path, required=True)
    parser.add_argument("--sentaurus-manifest", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument(
        "--mesh-file",
        type=Path,
        help="override base-config mesh_file with an exact paired physical mesh",
    )
    parser.add_argument(
        "--doping-file",
        type=Path,
        help="override base-config node_doping_file with paired nodal doping",
    )
    parser.add_argument("--branches", default=",".join(BRANCHES))
    parser.add_argument("--max-iter", type=int, default=80)
    parser.add_argument(
        "--continuation-schedule",
        choices=tuple(CONTINUATION_SCHEDULES),
        default="standard_0p05",
    )
    parser.add_argument("--resume", action="store_true")
    parser.add_argument(
        "--qf-carrier-truncation",
        type=float,
        help=(
            "opt-in low-density n,p floor relative to ni used only when "
            "rebuilding the avalanche quasi-Fermi driving field"
        ),
    )
    parser.add_argument(
        "--sg-laux-candidate",
        action="store_true",
        help=(
            "opt in to the complete element-edge SG/GSS-Laux current vector "
            "with matching element-vertex box source mapping"
        ),
    )
    parser.add_argument(
        "--terminal-current-method-compare",
        action="store_true",
        help=(
            "emit an observation-only comparison of SG-flux and continuity-"
            "residual terminal currents"
        ),
    )
    parser.add_argument(
        "--newton-reltol",
        type=float,
        help="observation-only override of solver.reltol",
    )
    parser.add_argument(
        "--newton-abstol",
        type=float,
        help="observation-only override of solver.abstol",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.max_iter <= 0:
        raise ValueError("--max-iter must be positive")
    for name, value in (
        ("--newton-reltol", args.newton_reltol),
        ("--newton-abstol", args.newton_abstol),
    ):
        if value is not None and (not math.isfinite(value) or value <= 0.0):
            raise ValueError(f"{name} must be positive and finite")
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
    base, physical_input_overrides = apply_physical_input_overrides(
        json.loads(base_path.read_text(encoding="utf-8-sig")),
        args.mesh_file,
        args.doping_file,
    )
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
            args.sg_laux_candidate,
            args.continuation_schedule,
            args.terminal_current_method_compare,
            args.newton_reltol,
            args.newton_abstol,
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
        "runner_sha256": sha256(runner),
        "base_config": str(base_path),
        "base_config_sha256": sha256(base_path),
        "physical_input_overrides": physical_input_overrides,
        "sentaurus_manifest": str(sentaurus_path),
        "sentaurus_manifest_sha256": sha256(sentaurus_path),
        "max_iter": args.max_iter,
        "continuation_schedule": {
            "id": args.continuation_schedule,
            **CONTINUATION_SCHEDULES[args.continuation_schedule],
        },
        "terminal_current_method_compare": {
            "enabled": args.terminal_current_method_compare,
            "observation_only": True,
        },
        "newton_tolerance_override": {
            "reltol": args.newton_reltol,
            "abstol": args.newton_abstol,
            "observation_only": True,
        },
        "candidate": (
            {
                "axis": "impact_ionization.element_edge_sg_gss_laux",
                "current_approximation": "element_edge_sg_gss_laux",
                "source_mapping_mode": "element_vertex_box_measure",
                "inactive_compatibility_normalization": {
                    "cell_reconstructed_midpoint_density": "bernoulli",
                },
                "default_unchanged": True,
            }
            if args.sg_laux_candidate
            else {
                "axis": "impact_ionization.quasi_fermi_carrier_truncation",
                "value": args.qf_carrier_truncation,
                "default_unchanged": True,
            }
        ),
        "requested_biases_V": biases,
        "branches": results,
    }
    execution_path = output_root / "execution.json"
    execution_path.write_text(
        json.dumps(execution, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    if complete:
        state_manifest = build_state_manifest(results, biases, output_root)
        state_manifest_path = output_root / "state_manifest.json"
        state_manifest_path.write_text(
            json.dumps(state_manifest, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    print(json.dumps(execution, indent=2, sort_keys=True))
    return 0 if complete else 2


if __name__ == "__main__":
    raise SystemExit(main())
