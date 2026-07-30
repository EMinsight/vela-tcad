#!/usr/bin/env python3
"""Audit Task 7 low-current reverse intervals without changing solver physics."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any


BIAS_TOLERANCE_V = 1.0e-9
STATE_COLUMNS = ("psi", "phin", "phip", "electrons_m3", "holes_m3")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def exact_row(
    rows: list[dict[str, str]],
    bias: float,
    *,
    column: str = "bias_V",
) -> dict[str, str]:
    matches = [
        row for row in rows
        if abs(float(row[column]) - bias) <= BIAS_TOLERANCE_V
    ]
    if len(matches) != 1:
        raise ValueError(f"{column}={bias:g}: expected one row, found {len(matches)}")
    return matches[0]


def bias_token(bias: float) -> str:
    sign = "m" if bias < 0.0 else ""
    return f"{sign}{abs(bias):.6f}".replace(".", "p")


def state_path(root: Path, branch: str, bias: float) -> Path:
    return root / branch / "states" / f"state_bias_{bias_token(bias)}.csv"


def state_difference(left: Path, right: Path) -> dict[str, Any]:
    left_rows = read_csv(left)
    right_rows = read_csv(right)
    if len(left_rows) != len(right_rows):
        raise ValueError("state node counts differ")
    result: dict[str, Any] = {
        "left_sha256": sha256(left),
        "right_sha256": sha256(right),
        "byte_identical": sha256(left) == sha256(right),
    }
    for column in STATE_COLUMNS:
        differences: list[float] = []
        for a, b in zip(left_rows, right_rows):
            av = float(a[column])
            bv = float(b[column])
            if column.endswith("_m3"):
                differences.append(
                    abs(math.log10(max(abs(av), 1.0)) -
                        math.log10(max(abs(bv), 1.0)))
                )
            else:
                differences.append(abs(av - bv))
        suffix = "log_dex" if column.endswith("_m3") else "abs"
        result[f"{column}_max_{suffix}"] = max(differences, default=0.0)
        result[f"{column}_rmse_{suffix}"] = (
            math.sqrt(sum(value * value for value in differences) / len(differences))
            if differences else 0.0
        )
    return result


def reverse_intervals(
    rows: list[dict[str, Any]],
    current_key: str,
) -> list[dict[str, float]]:
    ordered = sorted(rows, key=lambda row: abs(float(row["bias_V"])))
    result: list[dict[str, float]] = []
    for left, right in zip(ordered, ordered[1:]):
        left_current = abs(float(left[current_key]))
        right_current = abs(float(right[current_key]))
        if right_current < left_current:
            result.append({
                "left_bias_V": float(left["bias_V"]),
                "right_bias_V": float(right["bias_V"]),
                "left_abs_current_A_per_um": left_current,
                "right_abs_current_A_per_um": right_current,
            })
    return result


def exact_attempt(rows: list[dict[str, str]], bias: float) -> dict[str, str]:
    matches = [
        row for row in rows
        if (
            abs(float(row["requested_target_bias_V"]) - bias) <= BIAS_TOLERANCE_V
            and abs(float(row["actual_target_bias_V"]) - bias) <= BIAS_TOLERANCE_V
            and row["status"] == "accepted"
        )
    ]
    if len(matches) != 1:
        raise ValueError(f"{bias:g} V: expected one exact accepted attempt")
    return matches[0]


def process_summary(
    rows: list[dict[str, str]],
    point_index: int,
) -> dict[str, float | int]:
    selected = [
        row for row in rows
        if (
            int(row["point_index"]) == point_index
            and row["support_kind"] == "element_vertex_gss_laux"
        )
    ]
    if not selected:
        raise ValueError(f"point {point_index}: no element-vertex process rows")
    return {
        "process_row_count": len(selected),
        "qG_contribution_sum": sum(
            float(row["qG_contribution"] or 0.0) for row in selected
        ),
        "source_integral_sum": sum(
            float(row["source_integral"] or 0.0) for row in selected
        ),
        "alpha_max_per_m": max(float(row["alpha"] or 0.0) for row in selected),
        "impact_field_max_V_per_m": max(
            float(row["impact_field"] or 0.0) for row in selected
        ),
        "generation_rate_max_per_m3_s": max(
            float(row["generation_rate"] or 0.0) for row in selected
        ),
    }


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-a", type=Path, required=True)
    parser.add_argument("--run-b", type=Path, required=True)
    parser.add_argument(
        "--method-compare-run",
        type=Path,
        help=(
            "optional observation-only run containing per-branch "
            "terminal_current_method_compare.csv files"
        ),
    )
    parser.add_argument(
        "--strict-tolerance-run",
        type=Path,
        help=(
            "optional observation-only run using stricter Newton tolerances"
        ),
    )
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument(
        "--biases",
        nargs="+",
        type=float,
        default=[-3.0, -4.0, -5.0, -6.0, -7.0],
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    run_a = args.run_a.resolve()
    run_b = args.run_b.resolve()
    output = args.output_root.resolve()
    output.mkdir(parents=True, exist_ok=True)
    biases = [float(value) for value in args.biases]
    method_compare_root = (
        args.method_compare_run.resolve()
        if args.method_compare_run is not None else None
    )
    strict_tolerance_root = (
        args.strict_tolerance_run.resolve()
        if args.strict_tolerance_run is not None else None
    )

    branch_iv = {
        branch: read_csv(run_a / branch / "iv.csv")
        for branch in ("avalanche_off", "iic_postprocess", "avalanche_on")
    }
    attempts = {
        branch: read_csv(run_a / branch / "newton_attempts.csv")
        for branch in branch_iv
    }
    process_rows = read_csv(run_a / "avalanche_on" / "process_probe.csv")
    duplicate_iv = read_csv(run_b / "avalanche_on" / "iv.csv")

    bias_rows: list[dict[str, Any]] = []
    state_rows: list[dict[str, Any]] = []
    continuation_rows: list[dict[str, Any]] = []
    terminal_method_rows: list[dict[str, Any]] = []
    tolerance_rows: list[dict[str, Any]] = []
    selected_by_branch: dict[str, list[dict[str, Any]]] = {
        branch: [] for branch in branch_iv
    }

    for bias in biases:
        rows = {
            branch: exact_row(values, bias)
            for branch, values in branch_iv.items()
        }
        duplicate = exact_row(duplicate_iv, bias)
        point_index = branch_iv["avalanche_on"].index(rows["avalanche_on"])
        process = process_summary(process_rows, point_index)
        off_current = float(rows["avalanche_off"]["current_total_A_per_um"])
        on_current = float(rows["avalanche_on"]["current_total_A_per_um"])
        iic_current = float(rows["iic_postprocess"]["current_total_A_per_um"])
        on_hole = float(rows["avalanche_on"]["current_hole_A_per_um"])
        on_hole_drift = float(
            rows["avalanche_on"]["current_hole_drift_A_per_um"]
        )
        on_hole_diffusion = float(
            rows["avalanche_on"]["current_hole_diffusion_A_per_um"]
        )
        cancellation_ratio = (
            (abs(on_hole_drift) + abs(on_hole_diffusion)) / abs(on_hole)
            if on_hole != 0.0 else math.inf
        )
        record = {
            "bias_V": bias,
            "off_current_A_per_um": off_current,
            "iic_current_A_per_um": iic_current,
            "on_current_A_per_um": on_current,
            "on_minus_off_A_per_um": on_current - off_current,
            "on_over_off_abs": (
                abs(on_current) / abs(off_current)
                if off_current != 0.0 else math.inf
            ),
            "on_electron_A_per_um": float(
                rows["avalanche_on"]["current_electron_A_per_um"]
            ),
            "on_hole_A_per_um": on_hole,
            "on_hole_drift_A_per_um": on_hole_drift,
            "on_hole_diffusion_A_per_um": on_hole_diffusion,
            "hole_drift_diffusion_cancellation_ratio": cancellation_ratio,
            "on_newton_iterations": int(
                rows["avalanche_on"]["newton_iterations"]
            ),
            "on_convergence_reason": rows["avalanche_on"][
                "newton_convergence_reason"
            ],
            "duplicate_current_byte_value_equal": (
                rows["avalanche_on"]["current_total_A_per_um"]
                == duplicate["current_total_A_per_um"]
            ),
            **process,
        }
        bias_rows.append(record)
        for branch, row in rows.items():
            selected_by_branch[branch].append({
                "bias_V": bias,
                "current_A_per_um": float(row["current_total_A_per_um"]),
            })

        comparisons = {
            "off_vs_iic": (
                state_path(run_a, "avalanche_off", bias),
                state_path(run_a, "iic_postprocess", bias),
            ),
            "off_vs_on": (
                state_path(run_a, "avalanche_off", bias),
                state_path(run_a, "avalanche_on", bias),
            ),
            "on_a_vs_on_b": (
                state_path(run_a, "avalanche_on", bias),
                state_path(run_b, "avalanche_on", bias),
            ),
        }
        for comparison, (left, right) in comparisons.items():
            state_rows.append({
                "bias_V": bias,
                "comparison": comparison,
                **state_difference(left, right),
            })

        for branch, rows_for_branch in attempts.items():
            attempt = exact_attempt(rows_for_branch, bias)
            same_target = [
                row for row in rows_for_branch
                if abs(float(row["requested_target_bias_V"]) - bias)
                <= BIAS_TOLERANCE_V
            ]
            continuation_rows.append({
                "bias_V": bias,
                "branch": branch,
                "attempt_count_to_exact_target": len(same_target),
                "retry_count": sum(
                    1 for row in same_target if int(row["retry_number"]) > 0
                ),
                "exact_newton_iterations": int(attempt["newton_iterations"]),
                "exact_final_residual_norm": float(
                    attempt["final_residual_norm"]
                ),
                "exact_reason": attempt["reason"],
                "exact_final_state_hash": attempt["final_state_hash"],
            })

    if method_compare_root is not None:
        for branch in ("avalanche_off", "avalanche_on"):
            compare_rows = read_csv(
                method_compare_root
                / branch
                / "terminal_current_method_compare.csv"
            )
            for bias in biases:
                row = exact_row(compare_rows, bias)
                sg = float(row["I_sgflux_A_per_um"])
                residual = float(row["I_residual_A_per_um"])
                qf_floor = float(row["I_sgflux_with_qf_floor_A_per_um"])
                terminal_method_rows.append({
                    "bias_V": bias,
                    "branch": branch,
                    "I_sgflux_A_per_um": sg,
                    "I_residual_A_per_um": residual,
                    "I_sgflux_with_qf_floor_A_per_um": qf_floor,
                    "sg_minus_residual_A_per_um": sg - residual,
                    "sg_vs_residual_relative": (
                        abs(sg - residual) / abs(sg)
                        if sg != 0.0 else abs(sg - residual)
                    ),
                    "sg_minus_qf_floor_A_per_um": sg - qf_floor,
                    "anode_hole_qf_drop_V": float(
                        row["anode_hole_qf_drop_V"]
                    ),
                })

    strict_reverse: dict[str, list[dict[str, float]]] = {}
    if strict_tolerance_root is not None:
        for branch in ("avalanche_off", "avalanche_on"):
            strict_iv = read_csv(strict_tolerance_root / branch / "iv.csv")
            strict_attempts = read_csv(
                strict_tolerance_root / branch / "newton_attempts.csv"
            )
            strict_selected: list[dict[str, Any]] = []
            for bias in biases:
                standard = exact_row(branch_iv[branch], bias)
                strict = exact_row(strict_iv, bias)
                attempt = exact_attempt(strict_attempts, bias)
                standard_current = float(
                    standard["current_total_A_per_um"]
                )
                strict_current = float(strict["current_total_A_per_um"])
                tolerance_rows.append({
                    "bias_V": bias,
                    "branch": branch,
                    "standard_current_A_per_um": standard_current,
                    "strict_current_A_per_um": strict_current,
                    "strict_over_standard_abs": (
                        abs(strict_current) / abs(standard_current)
                        if standard_current != 0.0 else math.inf
                    ),
                    "strict_newton_iterations": int(
                        attempt["newton_iterations"]
                    ),
                    "strict_final_residual_norm": float(
                        attempt["final_residual_norm"]
                    ),
                    "strict_reason": attempt["reason"],
                })
                strict_selected.append({
                    "bias_V": bias,
                    "current_A_per_um": strict_current,
                })
            strict_reverse[branch] = reverse_intervals(
                strict_selected,
                "current_A_per_um",
            )

    reverse = {
        branch: reverse_intervals(rows, "current_A_per_um")
        for branch, rows in selected_by_branch.items()
    }
    off_pairs = {
        (row["left_bias_V"], row["right_bias_V"])
        for row in reverse["avalanche_off"]
    }
    on_pairs = {
        (row["left_bias_V"], row["right_bias_V"])
        for row in reverse["avalanche_on"]
    }
    shared_reverse = sorted(off_pairs & on_pairs)
    duplicate_passed = all(
        bool(row["byte_identical"])
        for row in state_rows if row["comparison"] == "on_a_vs_on_b"
    )
    off_iic_passed = all(
        bool(row["byte_identical"])
        for row in state_rows if row["comparison"] == "off_vs_iic"
    )
    no_retries = all(int(row["retry_count"]) == 0 for row in continuation_rows)
    high_cancellation = min(
        float(row["hole_drift_diffusion_cancellation_ratio"])
        for row in bias_rows
    ) >= 1.0e5
    source_monotonic = all(
        float(right["qG_contribution_sum"])
        > float(left["qG_contribution_sum"])
        for left, right in zip(bias_rows, bias_rows[1:])
    )
    terminal_methods_agree = (
        bool(terminal_method_rows)
        and max(
            float(row["sg_vs_residual_relative"])
            for row in terminal_method_rows
        ) <= 1.0e-10
        and all(
            float(row["sg_minus_qf_floor_A_per_um"]) == 0.0
            for row in terminal_method_rows
        )
    )
    qf_drop_at_precision_floor = (
        bool(terminal_method_rows)
        and max(
            abs(float(row["anode_hole_qf_drop_V"]))
            for row in terminal_method_rows
        ) <= 1.0e-12
    )
    standard_reverse_pairs = {
        branch: {
            (row["left_bias_V"], row["right_bias_V"])
            for row in reverse[branch]
        }
        for branch in ("avalanche_off", "avalanche_on")
    }
    strict_reverse_pairs = {
        branch: {
            (row["left_bias_V"], row["right_bias_V"])
            for row in strict_reverse.get(branch, [])
        }
        for branch in ("avalanche_off", "avalanche_on")
    }
    strict_tolerance_changes_pattern = (
        bool(tolerance_rows)
        and any(
            strict_reverse_pairs[branch] != standard_reverse_pairs[branch]
            for branch in standard_reverse_pairs
        )
    )
    strict_stalls_at_residual_floor = (
        bool(tolerance_rows)
        and all(
            row["strict_reason"] == "stall_residual_floor"
            for row in tolerance_rows
        )
    )
    if (
        duplicate_passed
        and off_iic_passed
        and no_retries
        and high_cancellation
        and source_monotonic
        and shared_reverse
        and terminal_methods_agree
        and qf_drop_at_precision_floor
        and strict_tolerance_changes_pattern
        and strict_stalls_at_residual_floor
    ):
        outcome = (
            "low_current_state_precision_floor_not_avalanche_operator_"
            "or_terminal_extractor"
        )
    else:
        outcome = "low_current_nonmonotonicity_requires_further_localization"

    result = {
        "schema": "vela.pn2d_task7_low_current_nonmonotonicity_audit.v1",
        "status": "passed",
        "outcome": outcome,
        "observation_only": True,
        "production_default_changed": False,
        "biases_V": biases,
        "reverse_intervals": reverse,
        "strict_tolerance_reverse_intervals": strict_reverse,
        "shared_off_on_reverse_intervals": [
            {"left_bias_V": left, "right_bias_V": right}
            for left, right in shared_reverse
        ],
        "gates": {
            "duplicate_on_state_hashes": duplicate_passed,
            "off_iic_state_hashes_identical": off_iic_passed,
            "no_continuation_retries": no_retries,
            "hole_drift_diffusion_cancellation_above_1e5": high_cancellation,
            "raw_avalanche_source_increases_with_reverse_bias": source_monotonic,
            "sg_and_residual_terminal_currents_agree": terminal_methods_agree,
            "contact_qf_drop_at_precision_floor": qf_drop_at_precision_floor,
            "strict_tolerance_changes_reverse_pattern":
                strict_tolerance_changes_pattern,
            "strict_tolerance_stalls_at_residual_floor":
                strict_stalls_at_residual_floor,
            "at_least_one_reverse_interval_shared_by_off_and_on": bool(
                shared_reverse
            ),
        },
        "artifacts": {
            "run_a": str(run_a),
            "run_b": str(run_b),
            "method_compare_run": (
                str(method_compare_root)
                if method_compare_root is not None else None
            ),
            "strict_tolerance_run": (
                str(strict_tolerance_root)
                if strict_tolerance_root is not None else None
            ),
        },
    }
    write_csv(output / "bias_diagnostics.csv", bias_rows)
    write_csv(output / "state_differences.csv", state_rows)
    write_csv(output / "continuation_diagnostics.csv", continuation_rows)
    write_csv(
        output / "terminal_current_method_compare.csv",
        terminal_method_rows,
    )
    write_csv(output / "tolerance_sensitivity.csv", tolerance_rows)
    (output / "acceptance.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
