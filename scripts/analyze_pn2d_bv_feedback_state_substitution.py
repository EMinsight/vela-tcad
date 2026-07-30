#!/usr/bin/env python3
"""Score cross-bias density/QFP feedback-state substitutions for PN2D Task 6."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable


VARIANTS = (
    "baseline",
    "electron_density_only",
    "hole_density_only",
    "density_only",
    "electron_qfp_only",
    "hole_qfp_only",
    "qfp_only",
    "density_qfp",
)
CARRIERS = ("electron", "hole")
CLOSURE_TOLERANCE = 1.0e-12
MATERIAL_IMPROVEMENT_FRACTION = 0.05
DIRECTION_COSINE_FLOOR = 0.05
NO_WORSENING_RELATIVE_TOLERANCE = 1.0e-12


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--execution", type=Path, required=True)
    parser.add_argument("--duplicate-execution", type=Path)
    parser.add_argument("--output-root", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def vector(rows: Iterable[dict[str, str]], column: str) -> list[float]:
    values = [float(row[column]) for row in rows]
    if not all(math.isfinite(value) for value in values):
        raise ValueError(f"{column}: nonfinite value")
    return values


def l2(values: Iterable[float]) -> float:
    return math.sqrt(sum(value * value for value in values))


def rmse(values: Iterable[float]) -> float:
    items = list(values)
    return l2(items) / math.sqrt(len(items)) if items else math.nan


def cosine(left: list[float], right: list[float]) -> float | None:
    left_norm = l2(left)
    right_norm = l2(right)
    if left_norm == 0.0 or right_norm == 0.0:
        return None
    return sum(a * b for a, b in zip(left, right)) / (left_norm * right_norm)


def improvement(baseline: float, candidate: float) -> float | None:
    if baseline <= 0.0 or not math.isfinite(baseline) or not math.isfinite(candidate):
        return None
    return (baseline - candidate) / baseline


def carrier_columns(carrier: str) -> dict[str, str]:
    prefix = "electron" if carrier == "electron" else "hole"
    return {
        "baseline": f"baseline_phi{'n' if carrier == 'electron' else 'p'}_V",
        "replacement": f"replacement_phi{'n' if carrier == 'electron' else 'p'}_V",
        "residual": f"{prefix}_residual",
        "desired": f"desired_{prefix}_residual",
        "delta": f"delta_phi{'n' if carrier == 'electron' else 'p'}_V",
        "trial": f"trial_phi{'n' if carrier == 'electron' else 'p'}_V",
        "carrier_only_delta":
            f"carrier_only_delta_phi{'n' if carrier == 'electron' else 'p'}_V",
        "carrier_only_trial":
            f"carrier_only_trial_phi{'n' if carrier == 'electron' else 'p'}_V",
        "closure": f"{prefix}_closure_error",
    }


def variant_metrics(
    rows: list[dict[str, str]],
    baseline_rows: list[dict[str, str]],
    variant: str,
) -> dict[str, Any]:
    interior = [row for row in rows if int(row["is_contact"]) == 0]
    baseline_interior = [
        row for row in baseline_rows if int(row["is_contact"]) == 0
    ]
    if not interior or len(interior) != len(baseline_interior):
        raise ValueError(f"{variant}: invalid interior support")
    carrier_metrics: dict[str, Any] = {}
    combined_initial_errors: list[float] = []
    combined_trial_errors: list[float] = []
    combined_baseline_trial_errors: list[float] = []
    combined_residual: list[float] = []
    combined_desired: list[float] = []
    combined_delta: list[float] = []
    combined_carrier_only_delta: list[float] = []
    combined_carrier_only_trial_errors: list[float] = []
    combined_target: list[float] = []
    max_closure = 0.0
    for carrier in CARRIERS:
        columns = carrier_columns(carrier)
        initial = vector(interior, columns["baseline"])
        replacement = vector(interior, columns["replacement"])
        trial = vector(interior, columns["trial"])
        baseline_trial = vector(baseline_interior, columns["trial"])
        target = [right - left for left, right in zip(initial, replacement)]
        initial_error = [left - right for left, right in zip(initial, replacement)]
        trial_error = [left - right for left, right in zip(trial, replacement)]
        baseline_trial_error = [
            left - right for left, right in zip(baseline_trial, replacement)
        ]
        residual = vector(interior, columns["residual"])
        desired = vector(interior, columns["desired"])
        delta = vector(interior, columns["delta"])
        carrier_only_delta = vector(interior, columns["carrier_only_delta"])
        carrier_only_trial = vector(interior, columns["carrier_only_trial"])
        carrier_only_trial_error = [
            left - right
            for left, right in zip(carrier_only_trial, replacement)
        ]
        closure = max(abs(value) for value in vector(rows, columns["closure"]))
        max_closure = max(max_closure, closure)
        initial_rmse = rmse(initial_error)
        trial_rmse = rmse(trial_error)
        baseline_trial_rmse = rmse(baseline_trial_error)
        carrier_metrics[carrier] = {
            "initial_qfp_error_rmse_V": initial_rmse,
            "trial_qfp_error_rmse_V": trial_rmse,
            "baseline_trial_qfp_error_rmse_V": baseline_trial_rmse,
            "qfp_error_improvement_fraction": improvement(
                baseline_trial_rmse, trial_rmse
            ),
            "residual_direction_cosine": cosine(residual, desired),
            "update_direction_cosine": cosine(delta, target),
            "carrier_only_trial_qfp_error_rmse_V": rmse(
                carrier_only_trial_error
            ),
            "carrier_only_qfp_error_improvement_fraction": improvement(
                baseline_trial_rmse, rmse(carrier_only_trial_error)
            ),
            "carrier_only_update_direction_cosine": cosine(
                carrier_only_delta, target
            ),
            "max_closure_error": closure,
        }
        combined_initial_errors.extend(initial_error)
        combined_trial_errors.extend(trial_error)
        combined_baseline_trial_errors.extend(baseline_trial_error)
        combined_residual.extend(residual)
        combined_desired.extend(desired)
        combined_delta.extend(delta)
        combined_carrier_only_delta.extend(carrier_only_delta)
        combined_carrier_only_trial_errors.extend(carrier_only_trial_error)
        combined_target.extend(target)

    combined_trial_rmse = rmse(combined_trial_errors)
    combined_baseline_trial_rmse = rmse(combined_baseline_trial_errors)
    combined_improvement = improvement(
        combined_baseline_trial_rmse, combined_trial_rmse
    )
    no_carrier_worsening = all(
        carrier_metrics[carrier]["trial_qfp_error_rmse_V"]
        <= carrier_metrics[carrier]["baseline_trial_qfp_error_rmse_V"]
        * (1.0 + NO_WORSENING_RELATIVE_TOLERANCE)
        for carrier in CARRIERS
    )
    carrier_only_no_carrier_worsening = all(
        carrier_metrics[carrier][
            "carrier_only_trial_qfp_error_rmse_V"
        ]
        <= carrier_metrics[carrier]["baseline_trial_qfp_error_rmse_V"]
        * (1.0 + NO_WORSENING_RELATIVE_TOLERANCE)
        for carrier in CARRIERS
    )
    residual_cosine = cosine(combined_residual, combined_desired)
    update_cosine = cosine(combined_delta, combined_target)
    causal_gate = (
        variant != "baseline"
        and max_closure <= CLOSURE_TOLERANCE
        and combined_improvement is not None
        and combined_improvement >= MATERIAL_IMPROVEMENT_FRACTION
        and residual_cosine is not None
        and residual_cosine >= DIRECTION_COSINE_FLOOR
        and update_cosine is not None
        and update_cosine >= DIRECTION_COSINE_FLOOR
        and no_carrier_worsening
    )
    return {
        "variant": variant,
        "node_count": len(rows),
        "interior_node_count": len(interior),
        "initial_qfp_error_rmse_V": rmse(combined_initial_errors),
        "trial_qfp_error_rmse_V": combined_trial_rmse,
        "baseline_trial_qfp_error_rmse_V": combined_baseline_trial_rmse,
        "qfp_error_improvement_fraction": combined_improvement,
        "residual_direction_cosine": residual_cosine,
        "update_direction_cosine": update_cosine,
        "carrier_only_trial_qfp_error_rmse_V": rmse(
            combined_carrier_only_trial_errors
        ),
        "carrier_only_qfp_error_improvement_fraction": improvement(
            combined_baseline_trial_rmse,
            rmse(combined_carrier_only_trial_errors),
        ),
        "carrier_only_update_direction_cosine": cosine(
            combined_carrier_only_delta, combined_target
        ),
        "carrier_only_no_carrier_worsening":
            carrier_only_no_carrier_worsening,
        "max_closure_error": max_closure,
        "no_carrier_worsening": no_carrier_worsening,
        "causal_gate_passed": causal_gate,
        "carriers": carrier_metrics,
    }


def boundary_identity(
    by_variant: dict[str, list[dict[str, str]]],
) -> dict[str, Any]:
    columns = ("psi_residual", "electron_residual", "hole_residual")
    baseline = {
        int(row["node_id"]): row
        for row in by_variant["baseline"]
        if int(row["is_contact"]) == 1
    }
    mismatches: list[dict[str, Any]] = []
    for variant in VARIANTS[1:]:
        for row in by_variant[variant]:
            if int(row["is_contact"]) == 0:
                continue
            node = int(row["node_id"])
            for column in columns:
                if row[column] != baseline[node][column]:
                    mismatches.append(
                        {
                            "variant": variant,
                            "node_id": node,
                            "column": column,
                            "baseline": baseline[node][column],
                            "candidate": row[column],
                        }
                    )
    return {
        "contact_node_count": len(baseline),
        "mismatch_count": len(mismatches),
        "passed": not mismatches,
        "first_mismatches": mismatches[:10],
    }


def classify(cross_bias: dict[str, bool]) -> str:
    electron_density = cross_bias.get("electron_density_only", False)
    hole_density = cross_bias.get("hole_density_only", False)
    density = cross_bias.get("density_only", False)
    electron_qfp = cross_bias.get("electron_qfp_only", False)
    hole_qfp = cross_bias.get("hole_qfp_only", False)
    qfp = cross_bias.get("qfp_only", False)
    combined = cross_bias.get("density_qfp", False)
    if electron_density and hole_density:
        return "carrier_split_density_feedback_cause"
    if electron_qfp and hole_qfp:
        return "carrier_split_qfp_feedback_cause"
    if electron_density:
        return "electron_density_feedback_cause"
    if hole_density:
        return "hole_density_feedback_cause"
    if electron_qfp:
        return "electron_qfp_feedback_cause"
    if hole_qfp:
        return "hole_qfp_feedback_cause"
    if density and qfp:
        return "density_and_qfp_independent_causes"
    if density:
        return "density_feedback_cause"
    if qfp:
        return "qfp_feedback_cause"
    if combined:
        return "density_qfp_interaction_cause"
    return "no_cross_bias_causal_substitution"


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    fields = (
        "bias_V",
        "variant",
        "causal_gate_passed",
        "initial_qfp_error_rmse_V",
        "baseline_trial_qfp_error_rmse_V",
        "trial_qfp_error_rmse_V",
        "qfp_error_improvement_fraction",
        "residual_direction_cosine",
        "update_direction_cosine",
        "carrier_only_trial_qfp_error_rmse_V",
        "carrier_only_qfp_error_improvement_fraction",
        "carrier_only_update_direction_cosine",
        "max_closure_error",
        "no_carrier_worsening",
        "carrier_only_no_carrier_worsening",
    )
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field) for field in fields})


def execution_output_hashes(execution: dict[str, Any]) -> dict[str, str]:
    return {
        f"{float(case['bias_V']):.17g}": sha256(Path(case["output_csv"]))
        for case in execution["cases"]
    }


def main() -> int:
    args = parse_args()
    execution_path = args.execution.resolve()
    execution = json.loads(execution_path.read_text(encoding="utf-8"))
    output = args.output_root.resolve()
    output.mkdir(parents=True, exist_ok=True)

    bias_results: list[dict[str, Any]] = []
    score_rows: list[dict[str, Any]] = []
    for case in sorted(execution["cases"], key=lambda item: float(item["bias_V"]), reverse=True):
        rows = read_csv(Path(case["output_csv"]))
        by_variant: dict[str, list[dict[str, str]]] = defaultdict(list)
        for row in rows:
            by_variant[row["variant"]].append(row)
        if set(by_variant) != set(VARIANTS):
            raise ValueError(f"{case['bias_V']}: incomplete variant matrix")
        baseline_rows = by_variant["baseline"]
        metrics = {
            variant: variant_metrics(
                by_variant[variant], baseline_rows, variant
            )
            for variant in VARIANTS
        }
        boundary = boundary_identity(by_variant)
        if not boundary["passed"]:
            for variant in VARIANTS[1:]:
                metrics[variant]["causal_gate_passed"] = False
        bias = float(case["bias_V"])
        bias_results.append(
            {
                "bias_V": bias,
                "boundary_identity": boundary,
                "variants": metrics,
            }
        )
        for variant in VARIANTS:
            score_rows.append({"bias_V": bias, **metrics[variant]})

    duplicate_path = (
        args.duplicate_execution.resolve()
        if args.duplicate_execution is not None
        else None
    )
    duplicate = (
        json.loads(duplicate_path.read_text(encoding="utf-8"))
        if duplicate_path is not None
        else None
    )
    primary_hashes = execution_output_hashes(execution)
    duplicate_hashes = (
        execution_output_hashes(duplicate)
        if duplicate is not None
        else {}
    )
    determinism = {
        "evaluated": duplicate is not None,
        "passed": duplicate is not None and primary_hashes == duplicate_hashes,
        "primary_output_sha256": primary_hashes,
        "duplicate_output_sha256": duplicate_hashes,
        "duplicate_execution": (
            {
                "path": str(duplicate_path),
                "sha256": sha256(duplicate_path),
            }
            if duplicate_path is not None
            else None
        ),
    }
    cross_bias = {
        variant: (
            len(bias_results) >= 2
            and determinism["passed"]
            and all(
                result["variants"][variant]["causal_gate_passed"]
                for result in bias_results
            )
        )
        for variant in VARIANTS[1:]
    }
    qfp_cross_block_reversal = (
        len(bias_results) >= 2
        and determinism["passed"]
        and all(
            result["boundary_identity"]["passed"]
            and result["variants"]["qfp_only"]["max_closure_error"]
                <= CLOSURE_TOLERANCE
            and result["variants"]["qfp_only"][
                "carrier_only_qfp_error_improvement_fraction"
            ] >= MATERIAL_IMPROVEMENT_FRACTION
            and result["variants"]["qfp_only"][
                "carrier_only_update_direction_cosine"
            ] >= DIRECTION_COSINE_FLOOR
            and result["variants"]["qfp_only"][
                "carrier_only_no_carrier_worsening"
            ]
            and (
                result["variants"]["qfp_only"][
                    "qfp_error_improvement_fraction"
                ] < 0.0
                or result["variants"]["qfp_only"][
                    "update_direction_cosine"
                ] < 0.0
            )
            for result in bias_results
        )
    )
    coupled_outcome = classify(cross_bias)
    outcome = (
        "continuation_only_cause"
        if coupled_outcome == "no_cross_bias_causal_substitution"
        and qfp_cross_block_reversal
        else coupled_outcome
    )
    causal_evidence_available = coupled_outcome != "no_cross_bias_causal_substitution"
    acceptance = {
        "schema": "vela.pn2d_bv_feedback_state_substitution_analysis.v1",
        "status": "passed",
        "outcome": outcome,
        "biases_V": [result["bias_V"] for result in bias_results],
        "cross_bias_causal_gate": cross_bias,
        "causal_evidence_available": causal_evidence_available,
        "cross_block_path_evidence": {
            "passed": qfp_cross_block_reversal,
            "localization": (
                "coupled_poisson_qfp_cross_block_reversal"
                if qfp_cross_block_reversal
                else None
            ),
            "variant": "qfp_only",
            "comparison":
                "carrier_only_first_update_vs_full_coupled_first_update",
        },
        "determinism": determinism,
        "task8_authorized": False,
        "next_gate": (
            "return_to_task7_with_evidence_authorized_single_axis_candidate"
            if causal_evidence_available
            else (
                "task7_continuation_schedule_branch_invariance"
                if qfp_cross_block_reversal
                else "retain_task8_stop_and_expand_task6_observation"
            )
        ),
        "thresholds": {
            "closure_absolute": CLOSURE_TOLERANCE,
            "material_qfp_error_improvement_fraction": MATERIAL_IMPROVEMENT_FRACTION,
            "direction_cosine_floor": DIRECTION_COSINE_FLOOR,
            "no_worsening_relative_tolerance": NO_WORSENING_RELATIVE_TOLERANCE,
        },
        "execution": {
            "path": str(execution_path),
            "sha256": sha256(execution_path),
        },
        "bias_results": bias_results,
    }
    acceptance_path = output / "acceptance.json"
    acceptance_path.write_text(
        json.dumps(acceptance, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_csv(output / "substitution_scorecard.csv", score_rows)
    print(json.dumps(acceptance, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
