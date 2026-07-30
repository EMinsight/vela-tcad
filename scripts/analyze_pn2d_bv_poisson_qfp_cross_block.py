#!/usr/bin/env python3
"""Localize PN2D QFP update reversal across Poisson-QFP Jacobian blocks."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any


MODES = (
    "independent",
    "no_psi_qfp",
    "no_qfp_psi",
    "schur",
    "full_raw",
    "full_capped",
)
CARRIERS = ("electron", "hole")
CLOSURE_TOLERANCE = 1.0e-10
SCHUR_FULL_MAX_DIFFERENCE_V = 1.0e-10
MATERIAL_IMPROVEMENT_FRACTION = 0.05
DIRECTION_COSINE_FLOOR = 0.05


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--execution", type=Path, required=True)
    parser.add_argument("--duplicate-execution", type=Path, required=True)
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


def l2(values: list[float]) -> float:
    return math.sqrt(sum(value * value for value in values))


def rmse(values: list[float]) -> float:
    return l2(values) / math.sqrt(len(values))


def cosine(left: list[float], right: list[float]) -> float | None:
    denominator = l2(left) * l2(right)
    if denominator == 0.0:
        return None
    return sum(a * b for a, b in zip(left, right)) / denominator


def optimal_projection_scale(
    delta: list[float],
    target: list[float],
) -> float | None:
    denominator = sum(value * value for value in delta)
    if denominator == 0.0:
        return None
    return sum(a * b for a, b in zip(delta, target)) / denominator


def carrier_columns(mode: str, carrier: str) -> str:
    suffix = "phin" if carrier == "electron" else "phip"
    return f"{mode}_delta_{suffix}_V"


def mode_metrics(
    rows: list[dict[str, str]],
    mode: str,
) -> dict[str, Any]:
    interior = [row for row in rows if int(row["is_contact"]) == 0]
    combined_initial: list[float] = []
    combined_trial: list[float] = []
    combined_delta: list[float] = []
    combined_target: list[float] = []
    carriers: dict[str, Any] = {}
    for carrier in CARRIERS:
        suffix = "phin" if carrier == "electron" else "phip"
        baseline = [float(row[f"baseline_{suffix}_V"]) for row in interior]
        replacement = [
            float(row[f"replacement_{suffix}_V"]) for row in interior
        ]
        target = [
            float(row[f"target_delta_{suffix}_V"]) for row in interior
        ]
        delta = [
            float(row[carrier_columns(mode, carrier)]) for row in interior
        ]
        initial_error = [
            left - right for left, right in zip(baseline, replacement)
        ]
        trial_error = [
            left + step - right
            for left, step, right in zip(baseline, delta, replacement)
        ]
        initial_rmse = rmse(initial_error)
        trial_rmse = rmse(trial_error)
        improvement = (
            (initial_rmse - trial_rmse) / initial_rmse
            if initial_rmse > 0.0
            else 0.0
        )
        carriers[carrier] = {
            "initial_qfp_error_rmse_V": initial_rmse,
            "trial_qfp_error_rmse_V": trial_rmse,
            "qfp_error_improvement_fraction": improvement,
            "update_direction_cosine": cosine(delta, target),
            "target_projection_optimal_scale":
                optimal_projection_scale(delta, target),
        }
        combined_initial.extend(initial_error)
        combined_trial.extend(trial_error)
        combined_delta.extend(delta)
        combined_target.extend(target)
    initial_rmse = rmse(combined_initial)
    trial_rmse = rmse(combined_trial)
    return {
        "mode": mode,
        "initial_qfp_error_rmse_V": initial_rmse,
        "trial_qfp_error_rmse_V": trial_rmse,
        "qfp_error_improvement_fraction": (
            (initial_rmse - trial_rmse) / initial_rmse
            if initial_rmse > 0.0
            else 0.0
        ),
        "update_direction_cosine": cosine(
            combined_delta, combined_target
        ),
        "target_projection_optimal_scale": optimal_projection_scale(
            combined_delta, combined_target
        ),
        "carriers": carriers,
    }


def classify_bias(metrics: dict[str, dict[str, Any]]) -> str:
    independent = metrics["independent"]
    no_b = metrics["no_psi_qfp"]
    no_c = metrics["no_qfp_psi"]
    full = metrics["full_raw"]
    direct_direction_good = (
        independent["update_direction_cosine"] is not None
        and independent["update_direction_cosine"] >= DIRECTION_COSINE_FLOOR
    )
    no_b_good = (
        no_b["update_direction_cosine"] is not None
        and no_b["update_direction_cosine"] >= DIRECTION_COSINE_FLOOR
    )
    no_c_good = (
        no_c["update_direction_cosine"] is not None
        and no_c["update_direction_cosine"] >= DIRECTION_COSINE_FLOOR
    )
    full_bad = (
        full["qfp_error_improvement_fraction"] < 0.0
        or (
            full["update_direction_cosine"] is not None
            and full["update_direction_cosine"] < 0.0
        )
    )
    if direct_direction_good and no_b_good and no_c_good and full_bad:
        return "bidirectional_poisson_qfp_closed_loop_cause"
    if direct_direction_good and not no_b_good:
        return "J_qfp_psi_direct_feed_cause"
    if direct_direction_good and not no_c_good:
        return "J_psi_qfp_counterfactual_cause"
    return "insufficient_cross_block_localization"


def boundary_gate(rows: list[dict[str, str]]) -> dict[str, Any]:
    contacts = [row for row in rows if int(row["is_contact"]) == 1]
    maximum = max(
        (
            abs(float(row[column]))
            for row in contacts
            for column in ("target_delta_phin_V", "target_delta_phip_V")
        ),
        default=0.0,
    )
    return {
        "contact_node_count": len(contacts),
        "maximum_target_delta_V": maximum,
        "passed": maximum == 0.0,
    }


def adverse_hotspot(rows: list[dict[str, str]]) -> dict[str, Any]:
    candidates: list[dict[str, Any]] = []
    for row in rows:
        if int(row["is_contact"]) == 1:
            continue
        electron_cross = float(
            row["full_minus_independent_delta_phin_V"]
        )
        hole_cross = float(row["full_minus_independent_delta_phip_V"])
        projection = (
            electron_cross * float(row["target_delta_phin_V"])
            + hole_cross * float(row["target_delta_phip_V"])
        )
        candidates.append(
            {
                "node_id": int(row["node_id"]),
                "x": float(row["x"]),
                "y": float(row["y"]),
                "cross_block_target_projection_V2": projection,
            }
        )
    return min(
        candidates,
        key=lambda row: row["cross_block_target_projection_V2"],
    )


def case_index(execution: dict[str, Any]) -> dict[float, dict[str, Any]]:
    return {
        float(case["bias_V"]): case
        for case in execution["cases"]
    }


def duplicate_hashes(execution: dict[str, Any]) -> dict[str, dict[str, str]]:
    return {
        f"{float(case['bias_V']):.17g}": {
            "node_csv": sha256(Path(case["output_csv"])),
            "jacobian_blocks_csv": sha256(
                Path(case["jacobian_blocks_csv"])
            ),
        }
        for case in execution["cases"]
    }


def write_scorecard(path: Path, results: list[dict[str, Any]]) -> None:
    fields = (
        "bias_V",
        "mode",
        "qfp_error_improvement_fraction",
        "update_direction_cosine",
        "electron_improvement_fraction",
        "electron_update_direction_cosine",
        "hole_improvement_fraction",
        "hole_update_direction_cosine",
    )
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for result in results:
            for mode in MODES:
                metric = result["modes"][mode]
                writer.writerow(
                    {
                        "bias_V": result["bias_V"],
                        "mode": mode,
                        "qfp_error_improvement_fraction":
                            metric["qfp_error_improvement_fraction"],
                        "update_direction_cosine":
                            metric["update_direction_cosine"],
                        "electron_improvement_fraction":
                            metric["carriers"]["electron"][
                                "qfp_error_improvement_fraction"
                            ],
                        "electron_update_direction_cosine":
                            metric["carriers"]["electron"][
                                "update_direction_cosine"
                            ],
                        "hole_improvement_fraction":
                            metric["carriers"]["hole"][
                                "qfp_error_improvement_fraction"
                            ],
                        "hole_update_direction_cosine":
                            metric["carriers"]["hole"][
                                "update_direction_cosine"
                            ],
                    }
                )


def main() -> int:
    args = parse_args()
    execution_path = args.execution.resolve()
    duplicate_path = args.duplicate_execution.resolve()
    execution = json.loads(execution_path.read_text(encoding="utf-8"))
    duplicate = json.loads(duplicate_path.read_text(encoding="utf-8"))
    cases = case_index(execution)
    duplicate_determinism = (
        duplicate_hashes(execution) == duplicate_hashes(duplicate)
    )

    bias_results: list[dict[str, Any]] = []
    closure_passed = True
    schur_full_passed = True
    boundary_passed = True
    for bias, case in sorted(cases.items(), reverse=True):
        rows = read_csv(Path(case["output_csv"]))
        status = json.loads(Path(case["status"]).read_text(encoding="utf-8"))
        metrics = {mode: mode_metrics(rows, mode) for mode in MODES}
        schur_difference = max(
            abs(
                float(row[f"schur_delta_{suffix}_V"])
                - float(row[f"full_raw_delta_{suffix}_V"])
            )
            for row in rows
            for suffix in ("psi", "phin", "phip")
        )
        closure = float(status["schur_relative_closure"])
        boundary = boundary_gate(rows)
        closure_passed = closure_passed and closure <= CLOSURE_TOLERANCE
        schur_full_passed = (
            schur_full_passed
            and schur_difference <= SCHUR_FULL_MAX_DIFFERENCE_V
        )
        boundary_passed = boundary_passed and boundary["passed"]
        block_products = {
            "J_psi_qfp_delta_qfp_l2": l2(
                [float(row["psi_qfp_product"]) for row in rows]
            ),
            "J_qfp_psi_delta_psi_l2": l2(
                [
                    float(row[column])
                    for row in rows
                    for column in (
                        "qfp_psi_electron_product",
                        "qfp_psi_hole_product",
                    )
                ]
            ),
        }
        bias_results.append(
            {
                "bias_V": bias,
                "classification": classify_bias(metrics),
                "boundary": boundary,
                "schur_relative_closure": closure,
                "schur_full_max_abs_difference_V": schur_difference,
                "jacobian_block_norms": status["jacobian_block_norms"],
                "block_products": block_products,
                "adverse_cross_block_hotspot": adverse_hotspot(rows),
                "modes": metrics,
            }
        )

    classifications = {
        result["classification"] for result in bias_results
    }
    cross_bias_cause = (
        len(bias_results) >= 2
        and len(classifications) == 1
        and next(iter(classifications))
        != "insufficient_cross_block_localization"
    )
    controls_valid = (
        duplicate_determinism
        and closure_passed
        and schur_full_passed
        and boundary_passed
    )
    outcome = (
        next(iter(classifications))
        if controls_valid and cross_bias_cause
        else "insufficient_cross_block_localization"
    )
    acceptance = {
        "schema": "vela.pn2d_bv_poisson_qfp_cross_block_analysis.v1",
        "status": "passed" if controls_valid else "failed",
        "outcome": outcome,
        "cross_bias_causal_evidence": cross_bias_cause,
        "task8_authorized": False,
        "production_defaults_changed": False,
        "gates": {
            "duplicate_determinism": duplicate_determinism,
            "schur_relative_closure": closure_passed,
            "schur_matches_full_raw_step": schur_full_passed,
            "boundary_target_preserved": boundary_passed,
            "same_localization_at_adjacent_biases": cross_bias_cause,
        },
        "thresholds": {
            "closure_relative": CLOSURE_TOLERANCE,
            "schur_full_max_abs_difference_V":
                SCHUR_FULL_MAX_DIFFERENCE_V,
            "material_qfp_error_improvement_fraction":
                MATERIAL_IMPROVEMENT_FRACTION,
            "direction_cosine_floor": DIRECTION_COSINE_FLOOR,
        },
        "bias_results": bias_results,
        "determinism": {
            "primary_hashes": duplicate_hashes(execution),
            "duplicate_hashes": duplicate_hashes(duplicate),
        },
        "artifacts": {
            "execution": {
                "path": str(execution_path),
                "sha256": sha256(execution_path),
            },
            "duplicate_execution": {
                "path": str(duplicate_path),
                "sha256": sha256(duplicate_path),
            },
        },
        "next_gate": (
            "review_bidirectional_cross_block_model_ownership"
            if outcome == "bidirectional_poisson_qfp_closed_loop_cause"
            else "retain_task8_stop_expand_observation"
        ),
    }
    output = args.output_root.resolve()
    output.mkdir(parents=True, exist_ok=True)
    acceptance_path = output / "acceptance.json"
    acceptance_path.write_text(
        json.dumps(acceptance, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_scorecard(output / "cross_block_scorecard.csv", bias_results)
    print(json.dumps(acceptance, indent=2, sort_keys=True))
    return 0 if controls_valid else 2


if __name__ == "__main__":
    raise SystemExit(main())
