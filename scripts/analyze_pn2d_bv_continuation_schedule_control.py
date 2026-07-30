#!/usr/bin/env python3
"""Score the predeclared PN2D continuation schedules without changing physics."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any


SCHEDULE_A = "standard_0p05"
SCHEDULE_B = "refined_0p025"
TARGET_BIASES = (-19.7, -19.8)
QFP_STATE_MAX_TOLERANCE_V = 1.0e-8
PSI_STATE_MAX_TOLERANCE_V = 1.0e-8
DENSITY_LOG_MAX_TOLERANCE_DEX = 1.0e-6
GLOBAL_CURRENT_LOG_TOLERANCE_DEX = 0.02
MATERIAL_IMPROVEMENT_FRACTION = 0.05
DIRECTION_COSINE_FLOOR = 0.05


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    for schedule in ("a", "b"):
        for duplicate in ("a", "b"):
            parser.add_argument(
                f"--schedule-{schedule}-execution-{duplicate}",
                type=Path,
                required=True,
            )
        parser.add_argument(
            f"--schedule-{schedule}-feedback",
            type=Path,
            required=True,
        )
    parser.add_argument("--output-root", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def avalanche_on(execution: dict[str, Any]) -> dict[str, Any]:
    matches = [
        branch
        for branch in execution["branches"]
        if branch["branch"] == "avalanche_on"
    ]
    if len(matches) != 1:
        raise ValueError("execution must contain exactly one avalanche_on branch")
    return matches[0]


def artifact_signature(execution: dict[str, Any]) -> dict[str, Any]:
    branch = avalanche_on(execution)
    return {
        "output_csv_sha256": branch["output_csv_sha256"],
        "state_sha256": {
            bias: record["sha256"]
            for bias, record in branch["state_files"].items()
        },
    }


def schedule_id(execution: dict[str, Any]) -> str:
    return str(execution["continuation_schedule"]["id"])


def exact_execution_valid(execution: dict[str, Any]) -> bool:
    branch = avalanche_on(execution)
    return (
        execution["status"] == "passed"
        and branch["returncode"] == 0
        and branch["complete_exact_lattice"]
        and branch["all_converged"]
        and branch["all_exact"]
        and len(branch["state_files"]) == branch["requested_bias_count"]
    )


def state_path(execution: dict[str, Any], bias: float) -> Path:
    states = avalanche_on(execution)["state_files"]
    matches = [
        Path(record["path"])
        for raw_bias, record in states.items()
        if abs(float(raw_bias) - bias) <= 1.0e-10
    ]
    if len(matches) != 1:
        raise ValueError(f"missing exact state at {bias:g} V")
    return matches[0]


def read_state(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if [int(row["node_id"]) for row in rows] != list(range(len(rows))):
        raise ValueError(f"{path}: noncanonical node order")
    return rows


def rmse(values: list[float]) -> float:
    return math.sqrt(sum(value * value for value in values) / len(values))


def state_metrics(left: Path, right: Path) -> dict[str, Any]:
    left_rows = read_state(left)
    right_rows = read_state(right)
    if len(left_rows) != len(right_rows):
        raise ValueError("state node counts differ")
    psi: list[float] = []
    qfp: list[float] = []
    density_log: list[float] = []
    for a, b in zip(left_rows, right_rows):
        psi.append(float(b["psi"]) - float(a["psi"]))
        for column in ("phin", "phip"):
            qfp.append(float(b[column]) - float(a[column]))
        for column in ("electrons_m3", "holes_m3"):
            av = float(a[column])
            bv = float(b[column])
            if av <= 0.0 or bv <= 0.0:
                raise ValueError(f"{column}: nonpositive density")
            density_log.append(math.log10(bv) - math.log10(av))
    return {
        "node_count": len(left_rows),
        "psi_rmse_V": rmse(psi),
        "psi_max_abs_V": max(abs(value) for value in psi),
        "qfp_rmse_V": rmse(qfp),
        "qfp_max_abs_V": max(abs(value) for value in qfp),
        "density_log_rmse_dex": rmse(density_log),
        "density_log_max_abs_dex": max(abs(value) for value in density_log),
    }


def curve_metrics(left: Path, right: Path) -> dict[str, Any]:
    with left.open(newline="", encoding="utf-8") as handle:
        left_rows = list(csv.DictReader(handle))
    with right.open(newline="", encoding="utf-8") as handle:
        right_rows = list(csv.DictReader(handle))
    if len(left_rows) != len(right_rows):
        raise ValueError("IV row counts differ")
    differences: list[float] = []
    rows: list[dict[str, float]] = []
    for a, b in zip(left_rows, right_rows):
        left_bias = float(a["bias_V"])
        right_bias = float(b["bias_V"])
        if abs(left_bias - right_bias) > 1.0e-10:
            raise ValueError("IV exact bias lattices differ")
        left_current = abs(float(a["current_total_A_per_um"]))
        right_current = abs(float(b["current_total_A_per_um"]))
        difference = math.log10(max(right_current, 1.0e-300)) - math.log10(
            max(left_current, 1.0e-300)
        )
        differences.append(difference)
        rows.append(
            {
                "bias_V": left_bias,
                "absolute_log_current_difference_dex": abs(difference),
            }
        )
    worst = max(
        rows,
        key=lambda row: row["absolute_log_current_difference_dex"],
    )
    return {
        "point_count": len(rows),
        "log_current_rmse_dex": rmse(differences),
        "log_current_max_abs_dex": worst[
            "absolute_log_current_difference_dex"
        ],
        "max_difference_bias_V": worst["bias_V"],
    }


def feedback_bias(
    acceptance: dict[str, Any],
    bias: float,
) -> dict[str, Any]:
    matches = [
        result
        for result in acceptance["bias_results"]
        if abs(float(result["bias_V"]) - bias) <= 1.0e-10
    ]
    if len(matches) != 1:
        raise ValueError(f"feedback scorecard is missing {bias:g} V")
    return matches[0]["variants"]["qfp_only"]


def reversal_metric(variant: dict[str, Any]) -> dict[str, Any]:
    carrier_improvement = float(
        variant["carrier_only_qfp_error_improvement_fraction"]
    )
    carrier_cosine = float(variant["carrier_only_update_direction_cosine"])
    full_improvement = float(variant["qfp_error_improvement_fraction"])
    full_cosine = float(variant["update_direction_cosine"])
    reversal = (
        carrier_improvement >= MATERIAL_IMPROVEMENT_FRACTION
        and carrier_cosine >= DIRECTION_COSINE_FLOOR
        and bool(variant["carrier_only_no_carrier_worsening"])
        and (full_improvement < 0.0 or full_cosine < 0.0)
    )
    resolved = (
        full_improvement >= MATERIAL_IMPROVEMENT_FRACTION
        and full_cosine >= DIRECTION_COSINE_FLOOR
        and bool(variant["no_carrier_worsening"])
    )
    return {
        "carrier_only_qfp_error_improvement_fraction": carrier_improvement,
        "carrier_only_update_direction_cosine": carrier_cosine,
        "full_coupled_qfp_error_improvement_fraction": full_improvement,
        "full_coupled_update_direction_cosine": full_cosine,
        "carrier_only_no_carrier_worsening": bool(
            variant["carrier_only_no_carrier_worsening"]
        ),
        "full_coupled_no_carrier_worsening": bool(
            variant["no_carrier_worsening"]
        ),
        "cross_block_reversal": reversal,
        "cross_block_reversal_resolved": resolved,
    }


def classify(
    *,
    controls_valid: bool,
    state_invariant: bool,
    baseline_reversal: bool,
    refined_reversal: bool,
    refined_resolved: bool,
) -> str:
    if not controls_valid:
        return "insufficient_or_nondeterministic_control"
    if baseline_reversal and refined_resolved and not refined_reversal:
        return "continuation_schedule_resolves_cross_block_reversal"
    if state_invariant and baseline_reversal and refined_reversal:
        return "continuation_invariant_cross_block_reversal"
    if not state_invariant and refined_reversal:
        return "continuation_branch_difference_without_metric_resolution"
    return "no_authorized_continuation_candidate"


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    execution_paths = {
        SCHEDULE_A: (
            args.schedule_a_execution_a.resolve(),
            args.schedule_a_execution_b.resolve(),
        ),
        SCHEDULE_B: (
            args.schedule_b_execution_a.resolve(),
            args.schedule_b_execution_b.resolve(),
        ),
    }
    feedback_paths = {
        SCHEDULE_A: args.schedule_a_feedback.resolve(),
        SCHEDULE_B: args.schedule_b_feedback.resolve(),
    }
    executions = {
        schedule: tuple(load_json(path) for path in paths)
        for schedule, paths in execution_paths.items()
    }
    feedback = {
        schedule: load_json(path)
        for schedule, path in feedback_paths.items()
    }

    schedule_identity = all(
        schedule_id(execution) == expected
        for expected, pair in executions.items()
        for execution in pair
    )
    exact_lattice = all(
        exact_execution_valid(execution)
        for pair in executions.values()
        for execution in pair
    )
    duplicate_determinism = {
        schedule: artifact_signature(pair[0]) == artifact_signature(pair[1])
        for schedule, pair in executions.items()
    }
    non_schedule_hashes = {
        avalanche_on(execution)["non_schedule_config_sha256"]
        for pair in executions.values()
        for execution in pair
    }
    physics_hashes = {
        avalanche_on(execution)["physics_config_sha256"]
        for pair in executions.values()
        for execution in pair
    }
    feedback_determinism = {
        schedule: bool(payload["determinism"]["passed"])
        for schedule, payload in feedback.items()
    }

    state_rows: list[dict[str, Any]] = []
    state_invariant = True
    for bias in TARGET_BIASES:
        metrics = state_metrics(
            state_path(executions[SCHEDULE_A][0], bias),
            state_path(executions[SCHEDULE_B][0], bias),
        )
        passed = (
            metrics["psi_max_abs_V"] <= PSI_STATE_MAX_TOLERANCE_V
            and metrics["qfp_max_abs_V"] <= QFP_STATE_MAX_TOLERANCE_V
            and metrics["density_log_max_abs_dex"]
            <= DENSITY_LOG_MAX_TOLERANCE_DEX
        )
        state_invariant = state_invariant and passed
        state_rows.append({"bias_V": bias, "passed": passed, **metrics})
    curve = curve_metrics(
        Path(avalanche_on(executions[SCHEDULE_A][0])["output_csv"]),
        Path(avalanche_on(executions[SCHEDULE_B][0])["output_csv"]),
    )
    curve_invariant = (
        curve["log_current_max_abs_dex"]
        <= GLOBAL_CURRENT_LOG_TOLERANCE_DEX
    )
    branch_invariant = state_invariant and curve_invariant

    internal_rows: list[dict[str, Any]] = []
    internal_by_schedule: dict[str, list[dict[str, Any]]] = {}
    for schedule, payload in feedback.items():
        rows: list[dict[str, Any]] = []
        for bias in TARGET_BIASES:
            metric = reversal_metric(feedback_bias(payload, bias))
            row = {"schedule": schedule, "bias_V": bias, **metric}
            rows.append(row)
            internal_rows.append(row)
        internal_by_schedule[schedule] = rows

    baseline_reversal = all(
        row["cross_block_reversal"]
        for row in internal_by_schedule[SCHEDULE_A]
    )
    refined_reversal = all(
        row["cross_block_reversal"]
        for row in internal_by_schedule[SCHEDULE_B]
    )
    refined_resolved = all(
        row["cross_block_reversal_resolved"]
        for row in internal_by_schedule[SCHEDULE_B]
    )
    controls_valid = (
        schedule_identity
        and exact_lattice
        and all(duplicate_determinism.values())
        and all(feedback_determinism.values())
        and len(non_schedule_hashes) == 1
        and len(physics_hashes) == 1
    )
    outcome = classify(
        controls_valid=controls_valid,
        state_invariant=branch_invariant,
        baseline_reversal=baseline_reversal,
        refined_reversal=refined_reversal,
        refined_resolved=refined_resolved,
    )
    schedule_resolves = outcome == (
        "continuation_schedule_resolves_cross_block_reversal"
    )
    acceptance = {
        "schema": "vela.pn2d_bv_continuation_schedule_control.v1",
        "status": "passed" if controls_valid else "failed",
        "outcome": outcome,
        "task7_outcome": (
            "solver_path_only_prequalification"
            if schedule_resolves
            else "no_authorized_candidate"
        ),
        "task8_authorized": False,
        "complete_curve_campaign_authorized": schedule_resolves,
        "production_defaults_changed": False,
        "schedules": {
            SCHEDULE_A: executions[SCHEDULE_A][0][
                "continuation_schedule"
            ],
            SCHEDULE_B: executions[SCHEDULE_B][0][
                "continuation_schedule"
            ],
        },
        "gates": {
            "schedule_identity": schedule_identity,
            "exact_lattice_complete": exact_lattice,
            "duplicate_determinism": duplicate_determinism,
            "feedback_determinism": feedback_determinism,
            "single_non_schedule_config_hash": len(non_schedule_hashes) == 1,
            "single_physics_config_hash": len(physics_hashes) == 1,
            "state_invariant_at_adjacent_biases": state_invariant,
            "curve_invariant_on_full_exact_lattice": curve_invariant,
            "baseline_cross_block_reversal_at_adjacent_biases":
                baseline_reversal,
            "refined_cross_block_reversal_at_adjacent_biases":
                refined_reversal,
            "refined_resolves_cross_block_reversal_at_adjacent_biases":
                refined_resolved,
        },
        "hashes": {
            "non_schedule_config_sha256": sorted(non_schedule_hashes),
            "physics_config_sha256": sorted(physics_hashes),
        },
        "thresholds": {
            "psi_state_max_abs_V": PSI_STATE_MAX_TOLERANCE_V,
            "qfp_state_max_abs_V": QFP_STATE_MAX_TOLERANCE_V,
            "density_log_state_max_abs_dex":
                DENSITY_LOG_MAX_TOLERANCE_DEX,
            "global_current_log_max_abs_dex":
                GLOBAL_CURRENT_LOG_TOLERANCE_DEX,
            "material_qfp_error_improvement_fraction":
                MATERIAL_IMPROVEMENT_FRACTION,
            "direction_cosine_floor": DIRECTION_COSINE_FLOOR,
        },
        "state_comparison": state_rows,
        "curve_comparison": curve,
        "internal_metric": internal_rows,
        "artifacts": {
            "executions": {
                schedule: [
                    {"path": str(path), "sha256": sha256(path)}
                    for path in paths
                ]
                for schedule, paths in execution_paths.items()
            },
            "feedback_scorecards": {
                schedule: {
                    "path": str(path),
                    "sha256": sha256(path),
                }
                for schedule, path in feedback_paths.items()
            },
        },
        "next_gate": (
            "task7_complete_curve_campaign"
            if schedule_resolves
            else "retain_task8_stop_no_continuation_candidate"
        ),
    }
    output = args.output_root.resolve()
    output.mkdir(parents=True, exist_ok=True)
    (output / "acceptance.json").write_text(
        json.dumps(acceptance, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_csv(output / "schedule_state_comparison.csv", state_rows)
    write_csv(output / "internal_metric_comparison.csv", internal_rows)
    print(json.dumps(acceptance, indent=2, sort_keys=True))
    return 0 if controls_valid else 2


if __name__ == "__main__":
    raise SystemExit(main())
