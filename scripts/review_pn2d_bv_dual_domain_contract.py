#!/usr/bin/env python3
"""Review and replay the prospective PN2D BV dual-domain acceptance contract."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from pathlib import Path
from typing import Any, Iterable, Mapping


EXACT_TOLERANCE_V = 1.0e-10


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def percentile(values: Iterable[float], fraction: float) -> float:
    ordered = sorted(float(value) for value in values)
    if not ordered:
        raise ValueError("cannot compute a percentile over an empty sequence")
    rank = fraction * (len(ordered) - 1)
    lower = math.floor(rank)
    upper = math.ceil(rank)
    if lower == upper:
        return ordered[lower]
    weight = rank - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def exact_curve_rows(
    rows: Iterable[Mapping[str, Any]],
    biases: Iterable[float],
) -> dict[float, Mapping[str, Any]]:
    source = list(rows)
    result: dict[float, Mapping[str, Any]] = {}
    for target in biases:
        matches = [
            row
            for row in source
            if abs(float(row["bias_V"]) - float(target)) <= EXACT_TOLERANCE_V
        ]
        if len(matches) != 1:
            raise ValueError(
                f"expected one exact curve row at {target:g} V, found {len(matches)}"
            )
        result[float(target)] = matches[0]
    return result


def abs_log_error(left: float, right: float) -> float:
    if left == 0.0 or right == 0.0:
        raise ValueError("log-error inputs must be non-zero")
    return abs(math.log10(abs(left)) - math.log10(abs(right)))


def required_gates(
    observed: Mapping[str, Any],
    names: Iterable[str],
) -> dict[str, bool]:
    return {name: bool(observed.get(name, False)) for name in names}


def evaluate(
    contract: Mapping[str, Any],
    curve: Mapping[str, Any],
    candidate_scorecard: Mapping[str, Any],
    frozen_state: Mapping[str, Any],
    low_current_audit: Mapping[str, Any],
) -> dict[str, Any]:
    if contract.get("schema") != "vela.pn2d_bv_dual_domain_acceptance_contract.v1":
        raise ValueError("unsupported dual-domain contract schema")
    if not bool(contract.get("prospective_only")):
        raise ValueError("contract must be prospective-only")
    if not bool(contract.get("retroactive_score_mutation_forbidden")):
        raise ValueError("contract must preserve the historical Task 7 score")

    domains = contract["domains"]
    bv_contract = domains["bv_model_consistency"]
    low_contract = domains["low_current_solver_precision"]
    bv_biases = tuple(float(value) for value in bv_contract["exact_biases_V"])
    low_biases = tuple(float(value) for value in low_contract["exact_biases_V"])
    if set(bv_biases) & set(low_biases):
        raise ValueError("BV-active and low-current bias domains overlap")

    bv_rows = exact_curve_rows(curve["curve_rows"], bv_biases)
    curve_errors = [
        abs_log_error(
            float(row["vela_on_A_per_um"]),
            float(row["sentaurus_on_A_per_um"]),
        )
        for row in bv_rows.values()
    ]
    gain_errors = [
        abs_log_error(float(row["vela_gain"]), float(row["sentaurus_gain"]))
        for row in bv_rows.values()
    ]
    effective_metrics = {
        "curve_median_abs_log_error_dex": statistics.median(curve_errors),
        "curve_p95_abs_log_error_dex": percentile(curve_errors, 0.95),
        "curve_max_abs_log_error_dex": max(curve_errors),
        "gain_median_abs_log_error_dex": statistics.median(gain_errors),
        "gain_max_abs_log_error_dex": max(gain_errors),
    }

    knee = curve["knee_metrics"]
    estimators = curve["knee_estimators"]
    vela_slope = estimators["vela"].get("V_slope")
    sentaurus_slope = estimators["sentaurus"].get("V_slope")
    if vela_slope is None or sentaurus_slope is None:
        slope_error = math.inf
    else:
        slope_error = abs(float(vela_slope) - float(sentaurus_slope))
    break_error = abs(
        float(estimators["vela"]["V_break"])
        - float(estimators["sentaurus"]["V_break"])
    )
    thresholds = bv_contract["thresholds"]
    metric_gates = {
        "effective_curve_median": (
            effective_metrics["curve_median_abs_log_error_dex"]
            <= float(thresholds["effective_curve_median_abs_log_error_dex"])
        ),
        "effective_curve_p95": (
            effective_metrics["curve_p95_abs_log_error_dex"]
            <= float(thresholds["effective_curve_p95_abs_log_error_dex"])
        ),
        "effective_curve_maximum": (
            effective_metrics["curve_max_abs_log_error_dex"]
            <= float(thresholds["effective_curve_max_abs_log_error_dex"])
        ),
        "effective_gain_median": (
            effective_metrics["gain_median_abs_log_error_dex"]
            <= float(thresholds["effective_gain_median_abs_log_error_dex"])
        ),
        "effective_gain_maximum": (
            effective_metrics["gain_max_abs_log_error_dex"]
            <= float(thresholds["effective_gain_max_abs_log_error_dex"])
        ),
        "knee_median": (
            float(knee["median_absolute_log_error_dex"])
            <= float(thresholds["knee_median_abs_log_error_dex"])
        ),
        "knee_maximum": (
            float(knee["maximum_absolute_log_error_dex"])
            <= float(thresholds["knee_max_abs_log_error_dex"])
        ),
        "V_break": break_error <= float(thresholds["V_break_abs_error_V"]),
        "V_slope": slope_error <= float(thresholds["V_slope_abs_error_V"]),
        "adjacent_slope_rmse": (
            float(curve["adjacent_slope_rmse_dex_per_V"])
            <= float(thresholds["adjacent_slope_rmse_dex_per_V"])
        ),
    }
    fixed_state_gates = required_gates(
        frozen_state["fixed_state_gates"],
        bv_contract["required_fixed_state_gates"],
    )
    self_consistent_gates = required_gates(
        candidate_scorecard["gates"],
        bv_contract["required_self_consistent_gates"],
    )
    bv_gates = {
        **metric_gates,
        **{f"fixed_state:{key}": value for key, value in fixed_state_gates.items()},
        **{
            f"self_consistent:{key}": value
            for key, value in self_consistent_gates.items()
        },
    }
    bv_passed = all(bv_gates.values())

    low_rows = exact_curve_rows(curve["curve_rows"], low_biases)
    current_limit = float(low_contract["eligibility_max_abs_current_A_per_um"])
    low_current_values = [
        abs(float(row[column]))
        for row in low_rows.values()
        for column in (
            "vela_on_A_per_um",
            "vela_off_A_per_um",
            "sentaurus_on_A_per_um",
            "sentaurus_off_A_per_um",
        )
    ]
    low_evidence_gates = required_gates(
        low_current_audit["gates"],
        low_contract["required_evidence_gates"],
    )
    low_gates = {
        "all_exact_currents_below_eligibility_limit": (
            max(low_current_values) <= current_limit
        ),
        "observation_only": bool(low_current_audit.get("observation_only")),
        "production_default_unchanged": (
            not bool(low_current_audit.get("production_default_changed"))
        ),
        "typed_outcome_matches": (
            low_current_audit.get("outcome")
            == low_contract["required_typed_outcome"]
        ),
        **{
            f"evidence:{key}": value
            for key, value in low_evidence_gates.items()
        },
    }
    low_classified = all(low_gates.values())

    policy = contract["decision_policy"]
    if not bv_passed:
        decision = policy["bv_fail"]
    elif low_classified:
        decision = policy["bv_pass_low_current_classified"]
    else:
        decision = policy["bv_pass_low_current_unclassified"]

    independent_review = {
        "bias_domains_disjoint": not bool(set(bv_biases) & set(low_biases)),
        "historical_score_preserved": (
            candidate_scorecard.get("outcome") == "tradeoff_without_parity"
        ),
        "low_current_raw_monotonicity_excluded_from_bv_gate": (
            "low_current_nonmonotonic_interval_removed"
            in bv_contract["excluded_gates"]
            and bool(low_contract["raw_monotonicity_is_not_a_model_gate"])
        ),
        "high_current_runaway_cannot_use_precision_waiver": (
            current_limit > 0.0 and math.isfinite(current_limit)
        ),
        "production_default_not_authorized": (
            not bool(contract.get("production_default_change_authorized"))
            and bool(
                policy["production_default_change_requires_separate_contract"]
            )
        ),
    }
    review_passed = all(independent_review.values())

    return {
        "schema": "vela.pn2d_bv_dual_domain_acceptance_review.v1",
        "status": "passed" if review_passed else "failed",
        "decision": decision,
        "next_stage_authorization": (
            policy["next_stage_authorized_by_pass"]
            if bv_passed and low_classified and review_passed
            else "none"
        ),
        "production_default_change_authorized": False,
        "historical_task7_outcome": candidate_scorecard.get("outcome"),
        "bv_model_consistency": {
            "passed": bv_passed,
            "biases_V": list(bv_biases),
            "metrics": {
                **effective_metrics,
                "knee_median_abs_log_error_dex": float(
                    knee["median_absolute_log_error_dex"]
                ),
                "knee_max_abs_log_error_dex": float(
                    knee["maximum_absolute_log_error_dex"]
                ),
                "V_break_abs_error_V": break_error,
                "V_slope_abs_error_V": slope_error,
                "adjacent_slope_rmse_dex_per_V": float(
                    curve["adjacent_slope_rmse_dex_per_V"]
                ),
            },
            "gates": bv_gates,
        },
        "low_current_solver_precision": {
            "classified_as_precision_floor": low_classified,
            "biases_V": list(low_biases),
            "maximum_abs_current_A_per_um": max(low_current_values),
            "eligibility_limit_A_per_um": current_limit,
            "gates": low_gates,
        },
        "independent_review": {
            "passed": review_passed,
            "checks": independent_review,
        },
    }


def write_gate_rows(path: Path, result: Mapping[str, Any]) -> None:
    rows: list[dict[str, Any]] = []
    for domain in ("bv_model_consistency", "low_current_solver_precision"):
        for gate, passed in result[domain]["gates"].items():
            rows.append({"domain": domain, "gate": gate, "passed": int(passed)})
    for gate, passed in result["independent_review"]["checks"].items():
        rows.append(
            {"domain": "independent_review", "gate": gate, "passed": int(passed)}
        )
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("domain", "gate", "passed"),
        )
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--curve-acceptance", type=Path, required=True)
    parser.add_argument("--candidate-scorecard", type=Path, required=True)
    parser.add_argument("--frozen-state-score", type=Path, required=True)
    parser.add_argument("--low-current-audit", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    inputs = {
        "contract": args.contract,
        "curve_acceptance": args.curve_acceptance,
        "candidate_scorecard": args.candidate_scorecard,
        "frozen_state_score": args.frozen_state_score,
        "low_current_audit": args.low_current_audit,
    }
    result = evaluate(
        load_json(args.contract),
        load_json(args.curve_acceptance),
        load_json(args.candidate_scorecard),
        load_json(args.frozen_state_score),
        load_json(args.low_current_audit),
    )
    result["artifacts"] = {
        name: {
            "path": str(path.resolve()),
            "sha256": sha256(path),
        }
        for name, path in inputs.items()
    }
    output = args.output_root.resolve()
    output.mkdir(parents=True, exist_ok=True)
    (output / "acceptance.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    write_gate_rows(output / "gates.csv", result)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
