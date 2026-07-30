#!/usr/bin/env python3
"""Score one PN2D density/QFP candidate against the frozen Task 7 gates."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Iterable


KNEE_BIASES = (-19.7, -19.8, -19.85, -19.9, -19.95, -20.0)
EXACT_TOLERANCE_V = 1.0e-10


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rmse(values: Iterable[float]) -> float:
    data = list(values)
    if not data:
        raise ValueError("cannot compute RMSE over an empty sequence")
    return math.sqrt(sum(value * value for value in data) / len(data))


def knee_rmse(acceptance: dict[str, Any]) -> float:
    return rmse(
        float(row["absolute_log_error_dex"])
        for row in acceptance["knee_error_rows"]
    )


def estimator_error(
    acceptance: dict[str, Any],
    estimator: str,
) -> float | None:
    estimators = acceptance["knee_estimators"]
    sentaurus = estimators["sentaurus"].get(estimator)
    vela = estimators["vela"].get(estimator)
    if sentaurus is None or vela is None:
        return None
    return abs(float(vela) - float(sentaurus))


def estimator_improved(
    baseline_error: float | None,
    candidate_error: float | None,
) -> bool:
    if candidate_error is None:
        return False
    return baseline_error is None or candidate_error < baseline_error


def feedback_initial_qfp_rmse(payload: dict[str, Any]) -> float:
    values = [
        float(row["variants"]["baseline"]["initial_qfp_error_rmse_V"])
        for row in payload["bias_results"]
    ]
    return rmse(values)


def vela_on_reverse_intervals(
    acceptance: dict[str, Any],
) -> list[dict[str, float]]:
    rows = sorted(
        acceptance["curve_rows"],
        key=lambda row: abs(float(row["bias_V"])),
    )
    result: list[dict[str, float]] = []
    for left, right in zip(rows, rows[1:]):
        left_current = abs(float(left["vela_on_A_per_um"]))
        right_current = abs(float(right["vela_on_A_per_um"]))
        if right_current < left_current:
            result.append(
                {
                    "left_bias_V": float(left["bias_V"]),
                    "right_bias_V": float(right["bias_V"]),
                    "left_abs_current_A_per_um": left_current,
                    "right_abs_current_A_per_um": right_current,
                }
            )
    return result


def chain_index(
    payload: dict[str, Any],
    *,
    branch: str,
    stage: str,
    quantity: str,
) -> dict[tuple[float, str, str], float]:
    result: dict[tuple[float, str, str], float] = {}
    for record in payload["records"]:
        bias = float(record["bias_V"])
        if (
            record["branch"] != branch
            or record["stage"] != stage
            or record["quantity"] != quantity
            or not any(abs(bias - target) <= EXACT_TOLERANCE_V for target in KNEE_BIASES)
        ):
            continue
        values = record.get("values", [])
        if len(values) != 1:
            continue
        key = (bias, str(record["carrier"]), str(record["support_key"]))
        if key in result:
            raise ValueError(f"duplicate chain observation: {key}")
        result[key] = float(values[0])
    return result


def paired_metric(
    sentaurus: dict[str, Any],
    vela: dict[str, Any],
    *,
    stage: str,
    quantity: str,
    logarithmic: bool,
) -> float:
    left = chain_index(
        sentaurus,
        branch="avalanche_on",
        stage=stage,
        quantity=quantity,
    )
    right = chain_index(
        vela,
        branch="avalanche_on",
        stage=stage,
        quantity=quantity,
    )
    common = left.keys() & right.keys()
    if not common:
        raise ValueError(f"{stage}/{quantity} chain supports do not overlap")
    differences: list[float] = []
    for key in common:
        a = left[key]
        b = right[key]
        if logarithmic:
            if a <= 0.0 or b <= 0.0:
                continue
            differences.append(math.log10(b) - math.log10(a))
        else:
            differences.append(b - a)
    return rmse(differences)


def max_global_worsening(
    baseline: dict[str, Any],
    candidate: dict[str, Any],
) -> float:
    base = {
        float(row["bias_V"]): float(row["absolute_log_error_dex"])
        for row in baseline["global_error_rows"]
    }
    cand = {
        float(row["bias_V"]): float(row["absolute_log_error_dex"])
        for row in candidate["global_error_rows"]
    }
    if base.keys() != cand.keys():
        raise ValueError("baseline and candidate global lattices differ")
    return max(cand[bias] - base[bias] for bias in base)


def score(args: argparse.Namespace) -> dict[str, Any]:
    baseline_curve = load_json(args.baseline_curve)
    candidate_curve = load_json(args.candidate_curve)
    sentaurus_chain = load_json(args.sentaurus_chain)
    baseline_chain = load_json(args.baseline_vela_chain)
    candidate_chain = load_json(args.candidate_vela_chain)
    baseline_wp7 = load_json(args.baseline_wp7)
    candidate_wp7 = load_json(args.candidate_wp7)
    baseline_feedback = (
        load_json(args.baseline_feedback)
        if args.baseline_feedback is not None
        else None
    )
    candidate_feedback = (
        load_json(args.candidate_feedback)
        if args.candidate_feedback is not None
        else None
    )

    baseline_knee = knee_rmse(baseline_curve)
    candidate_knee = knee_rmse(candidate_curve)
    improvement = (
        (baseline_knee - candidate_knee) / baseline_knee
        if baseline_knee > 0.0
        else 0.0
    )
    baseline_break = estimator_error(baseline_curve, "V_break")
    candidate_break = estimator_error(candidate_curve, "V_break")
    baseline_slope = estimator_error(baseline_curve, "V_slope")
    candidate_slope = estimator_error(candidate_curve, "V_slope")
    baseline_qfp = paired_metric(
        sentaurus_chain,
        baseline_chain,
        stage="state",
        quantity="quasi_fermi",
        logarithmic=False,
    )
    candidate_qfp = paired_metric(
        sentaurus_chain,
        candidate_chain,
        stage="state",
        quantity="quasi_fermi",
        logarithmic=False,
    )
    baseline_density = paired_metric(
        sentaurus_chain,
        baseline_chain,
        stage="density",
        quantity="density",
        logarithmic=True,
    )
    candidate_density = paired_metric(
        sentaurus_chain,
        candidate_chain,
        stage="density",
        quantity="density",
        logarithmic=True,
    )
    duplicate_deterministic = sha256(args.candidate_iv_a) == sha256(
        args.candidate_iv_b
    )
    chain_records_identical = baseline_chain["records"] == candidate_chain["records"]
    global_worsening = max_global_worsening(baseline_curve, candidate_curve)
    baseline_feedback_qfp = (
        feedback_initial_qfp_rmse(baseline_feedback)
        if baseline_feedback is not None
        else None
    )
    candidate_feedback_qfp = (
        feedback_initial_qfp_rmse(candidate_feedback)
        if candidate_feedback is not None
        else None
    )
    feedback_inputs_complete = (
        baseline_feedback is not None and candidate_feedback is not None
    )
    baseline_reverse_intervals = vela_on_reverse_intervals(baseline_curve)
    candidate_reverse_intervals = vela_on_reverse_intervals(candidate_curve)

    gates = {
        "knee_rmse_improves_at_least_50_percent": improvement >= 0.5,
        "V_break_error_reduced": estimator_improved(
            baseline_break, candidate_break
        ),
        "V_slope_error_reduced": estimator_improved(
            baseline_slope, candidate_slope
        ),
        "nonmonotonic_interval_removed": (
            not candidate_reverse_intervals and bool(baseline_reverse_intervals)
        ),
        "qfp_state_metric_improved": candidate_qfp < baseline_qfp,
        "density_state_metric_improved": candidate_density < baseline_density,
        "wp7_closure_preserved": (
            candidate_wp7["status"] == "passed"
            and int(candidate_wp7["failed_closure_rows"]) == 0
            and not candidate_wp7["missing_observation"]
        ),
        "duplicate_determinism": duplicate_deterministic,
        "global_error_worsening_within_0p02_dex": global_worsening <= 0.02,
        "task6_feedback_state_metric_improved": (
            feedback_inputs_complete
            and candidate_feedback_qfp is not None
            and baseline_feedback_qfp is not None
            and candidate_feedback_qfp < baseline_feedback_qfp
        ),
        "task6_feedback_determinism_preserved": (
            feedback_inputs_complete
            and candidate_feedback["status"] == "passed"
            and bool(candidate_feedback["determinism"]["passed"])
        ),
    }
    authorized = all(gates.values())
    if authorized:
        outcome = "single_causal_candidate_authorized"
    elif improvement >= 0.5 and not (
        gates["qfp_state_metric_improved"]
        and gates["density_state_metric_improved"]
    ):
        outcome = "improves_curve_without_internal_causality"
    elif chain_records_identical and candidate_knee == baseline_knee:
        outcome = "no_authorized_candidate"
    else:
        outcome = "tradeoff_without_parity"

    return {
        "schema": "vela.pn2d_bv_task7_candidate_scorecard.v1",
        "status": "passed",
        "outcome": outcome,
        "candidate_axis": args.candidate_axis,
        "candidate_value": args.candidate_value,
        "production_default_changed": False,
        "metrics": {
            "baseline_knee_log_current_rmse_dex": baseline_knee,
            "candidate_knee_log_current_rmse_dex": candidate_knee,
            "knee_rmse_improvement_fraction": improvement,
            "baseline_V_break_error_V": baseline_break,
            "candidate_V_break_error_V": candidate_break,
            "baseline_V_slope_error_V": baseline_slope,
            "candidate_V_slope_error_V": candidate_slope,
            "baseline_qfp_rmse_V": baseline_qfp,
            "candidate_qfp_rmse_V": candidate_qfp,
            "baseline_density_log_rmse_dex": baseline_density,
            "candidate_density_log_rmse_dex": candidate_density,
            "baseline_task6_initial_qfp_rmse_V": baseline_feedback_qfp,
            "candidate_task6_initial_qfp_rmse_V": candidate_feedback_qfp,
            "baseline_task6_outcome": (
                baseline_feedback["outcome"]
                if baseline_feedback is not None
                else None
            ),
            "candidate_task6_outcome": (
                candidate_feedback["outcome"]
                if candidate_feedback is not None
                else None
            ),
            "baseline_vela_on_reverse_intervals": baseline_reverse_intervals,
            "candidate_vela_on_reverse_intervals": candidate_reverse_intervals,
            "maximum_global_error_worsening_dex": global_worsening,
            "chain_records_identical": chain_records_identical,
            "baseline_wp7_outcome": baseline_wp7["outcome"],
            "candidate_wp7_outcome": candidate_wp7["outcome"],
        },
        "gates": gates,
        "artifacts": {
            "baseline_curve": str(args.baseline_curve.resolve()),
            "candidate_curve": str(args.candidate_curve.resolve()),
            "baseline_vela_chain": str(args.baseline_vela_chain.resolve()),
            "candidate_vela_chain": str(args.candidate_vela_chain.resolve()),
            "candidate_iv_a_sha256": sha256(args.candidate_iv_a),
            "candidate_iv_b_sha256": sha256(args.candidate_iv_b),
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-curve", type=Path, required=True)
    parser.add_argument("--candidate-curve", type=Path, required=True)
    parser.add_argument("--sentaurus-chain", type=Path, required=True)
    parser.add_argument("--baseline-vela-chain", type=Path, required=True)
    parser.add_argument("--candidate-vela-chain", type=Path, required=True)
    parser.add_argument("--baseline-wp7", type=Path, required=True)
    parser.add_argument("--candidate-wp7", type=Path, required=True)
    parser.add_argument("--baseline-feedback", type=Path)
    parser.add_argument("--candidate-feedback", type=Path)
    parser.add_argument("--candidate-iv-a", type=Path, required=True)
    parser.add_argument("--candidate-iv-b", type=Path, required=True)
    parser.add_argument(
        "--candidate-axis",
        default="impact_ionization.quasi_fermi_carrier_truncation",
    )
    parser.add_argument("--candidate-value", default="1.0e-2")
    parser.add_argument("--output-root", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    result = score(args)
    output = args.output_root.resolve()
    output.mkdir(parents=True, exist_ok=True)
    (output / "acceptance.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    with (output / "candidate_scorecard.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=("gate", "passed"))
        writer.writeheader()
        for gate, passed in result["gates"].items():
            writer.writerow({"gate": gate, "passed": int(passed)})
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
