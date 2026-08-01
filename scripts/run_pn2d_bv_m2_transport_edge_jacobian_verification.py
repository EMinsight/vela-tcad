#!/usr/bin/env python3
"""Decompose M2 hotspot-edge transport Jacobians without changing physics."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
import sys
from pathlib import Path
from typing import Any

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.diagnose_pn2d_bv_predictor_first_step_audit import make_probe_config
from scripts.run_pn2d_bv_m2_single_family_state_substitution import bias_tag


DEFAULT_BIASES = (-18.0, -19.5, -19.7, -20.0)
VARIANTS = ("vela_baseline", "sent_qfp_only")
CARRIERS = ("electron", "hole")
FD_STEP_V = 1.0e-7


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise RuntimeError(f"refusing to write empty CSV: {path}")
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def relative_error(left: float, right: float) -> float:
    return abs(left - right) / max(abs(left), abs(right), 1.0e-300)


def run_probe(runner: Path, config_path: Path) -> dict[str, Any]:
    completed = subprocess.run(
        [str(runner.resolve()), "--config", str(config_path.resolve())],
        cwd=config_path.parent,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    config_path.with_suffix(".stdout.log").write_text(completed.stdout, encoding="utf-8")
    config_path.with_suffix(".stderr.log").write_text(completed.stderr, encoding="utf-8")
    if completed.returncode:
        raise RuntimeError(
            f"runner failed for {config_path}: "
            f"{completed.stderr.strip() or completed.stdout.strip()}"
        )
    return json.loads(completed.stdout)


def representative_flux(
    rows: list[dict[str, str]],
) -> dict[tuple[int, str], dict[str, float | int]]:
    result: dict[tuple[int, str], dict[str, float | int]] = {}
    for row in rows:
        if int(row["row_endpoint"]) == 0 and int(row["column_endpoint"]) == 0:
            result[(int(row["edge_id"]), row["carrier"])] = {
                "flux": float(row["flux_physical"]),
                "node0": int(row["node0"]),
                "node1": int(row["node1"]),
            }
    return result


def select_hotspots(
    bias: float,
    baseline: list[dict[str, str]],
    replacement: list[dict[str, str]],
    residual_hotspot_nodes: dict[str, int],
) -> list[dict[str, object]]:
    base_flux = representative_flux(baseline)
    replacement_flux = representative_flux(replacement)
    if set(base_flux) != set(replacement_flux):
        raise RuntimeError(f"{bias:g} V: transport edge sets differ")
    output: list[dict[str, object]] = []
    for carrier in CARRIERS:
        candidates = []
        residual_hotspot_node = residual_hotspot_nodes[carrier]
        for (edge_id, item_carrier), left_record in base_flux.items():
            if item_carrier != carrier:
                continue
            if residual_hotspot_node not in {
                int(left_record["node0"]), int(left_record["node1"])
            }:
                continue
            left = float(left_record["flux"])
            right = float(replacement_flux[(edge_id, carrier)]["flux"])
            candidates.append((abs(right - left), edge_id, left, right))
        if not candidates:
            raise RuntimeError(f"{bias:g} V: no {carrier} transport edges")
        absolute_delta, edge_id, left, right = max(
            candidates, key=lambda item: (item[0], -item[1])
        )
        output.append(
            {
                "bias_V": bias,
                "carrier": carrier,
                "edge_id": edge_id,
                "residual_hotspot_node_id": residual_hotspot_node,
                "selection_contract": "incident_to_prior_interior_flux_residual_hotspot__max_abs_baseline_to_joint_qfp_flux_delta__tie_lowest_edge_id",
                "baseline_flux_physical": left,
                "joint_qfp_flux_physical": right,
                "signed_flux_delta_physical": right - left,
                "absolute_flux_delta_physical": absolute_delta,
            }
        )
    return output


def decompose_selected(
    bias: float,
    variant: str,
    rows: list[dict[str, str]],
    selected: dict[str, int],
) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    for row in rows:
        carrier = row["carrier"]
        if int(row["edge_id"]) != selected[carrier]:
            continue
        mobility = float(row["mobility_m2_V_s"])
        d_mobility = float(row["d_mobility_d_column_m2_V2_s"])
        flux = float(row["flux_physical"])
        mobility_expected = flux / mobility * d_mobility if mobility else 0.0
        mobility_observed = float(row["mobility_response_derivative_physical"])
        production = float(row["production_frozen_mobility_derivative_physical"])
        frozen_fd = float(row["frozen_mobility_fd_derivative_physical"])
        live_fd = float(row["live_mobility_fd_derivative_physical"])
        live_expected = production + mobility_expected
        contact_eliminated = float(row["contact_eliminated_production_edge_derivative"])
        weight = float(row["continuity_row_weight"])
        solver_derivative = float(row["solver_production_edge_derivative"])
        item: dict[str, object] = {
            "bias_V": bias,
            "variant": variant,
            **row,
            "analytic_to_frozen_fd_relative_error": relative_error(production, frozen_fd),
            "mobility_product_expected_derivative_physical": mobility_expected,
            "mobility_product_relative_error": relative_error(
                mobility_observed, mobility_expected
            ),
            "live_derivative_expected_from_production_plus_mobility": live_expected,
            "live_derivative_elementwise_relative_error": relative_error(
                live_fd, live_expected
            ),
            "mobility_response_fraction_of_live_derivative": abs(mobility_observed)
            / max(abs(live_fd), 1.0e-300),
            "row_scaling_closure_relative_error": relative_error(
                solver_derivative, contact_eliminated * weight
            ),
        }
        output.append(item)
    if len(output) != 8:
        raise RuntimeError(
            f"{bias:g} V {variant}: expected 8 selected records, got {len(output)}"
        )
    production_scale = max(
        abs(float(row["production_frozen_mobility_derivative_physical"]))
        for row in output
    )
    live_scale = max(
        abs(float(row["live_mobility_fd_derivative_physical"])) for row in output
    )
    for row in output:
        row["analytic_to_frozen_fd_edge_scaled_error"] = abs(
            float(row["production_frozen_mobility_derivative_physical"])
            - float(row["frozen_mobility_fd_derivative_physical"])
        ) / max(production_scale, 1.0e-300)
        row["live_total_edge_scaled_error"] = abs(
            float(row["live_mobility_fd_derivative_physical"])
            - float(row["live_derivative_expected_from_production_plus_mobility"])
        ) / max(live_scale, 1.0e-300)
    return output


def contact_audit(
    bias: float, variant: str, rows: list[dict[str, str]]
) -> dict[str, object]:
    constrained = [row for row in rows if int(row["row_constrained"])]
    eliminated_errors = [
        abs(float(row["contact_eliminated_production_edge_derivative"]))
        for row in constrained
    ]
    identity_errors = []
    for row in constrained:
        expected = int(row["row_node"] == row["column_node"])
        identity_errors.append(abs(float(row["contact_identity_entry"]) - expected))
    preserved_columns = [
        row
        for row in rows
        if not int(row["row_constrained"])
        and int(row["column_constrained"])
        and float(row["production_row_derivative_scaled"]) != 0.0
        and float(row["contact_eliminated_production_edge_derivative"]) != 0.0
    ]
    return {
        "bias_V": bias,
        "variant": variant,
        "constrained_edge_row_record_count": len(constrained),
        "maximum_contact_eliminated_edge_derivative_abs": max(eliminated_errors, default=0.0),
        "maximum_contact_identity_error_abs": max(identity_errors, default=0.0),
        "unconstrained_rows_with_preserved_constrained_column_count": len(preserved_columns),
        "contact_semantics": "replace_constrained_rows_only__preserve_constrained_columns",
    }


def summarize_hotspot_states(
    decomposition: list[dict[str, object]],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    grouped: dict[tuple[float, str, str, int], list[dict[str, object]]] = {}
    for row in decomposition:
        key = (
            float(row["bias_V"]), str(row["variant"]), str(row["carrier"]),
            int(row["edge_id"]),
        )
        grouped.setdefault(key, []).append(row)
    state_rows: list[dict[str, object]] = []
    for (bias, variant, carrier, edge_id), rows in sorted(grouped.items()):
        columns = [row for row in rows if int(row["row_endpoint"]) == 0]
        dominant = max(
            columns,
            key=lambda row: abs(
                float(row["production_frozen_mobility_derivative_physical"])
            ),
        )
        production = abs(
            float(dominant["production_frozen_mobility_derivative_physical"])
        )
        mobility_response = max(
            abs(float(row["mobility_response_derivative_physical"]))
            for row in columns
        )
        state_rows.append({
            "bias_V": bias,
            "variant": variant,
            "carrier": carrier,
            "edge_id": edge_id,
            "dominant_qfp_column_node": int(dominant["column_node"]),
            "qfp_drive_V_m": float(dominant["qfp_drive_V_m"]),
            "mobility_m2_V_s": float(dominant["mobility_m2_V_s"]),
            "bernoulli_node0": float(dominant["bernoulli_node0"]),
            "bernoulli_node1": float(dominant["bernoulli_node1"]),
            "carrier_density_node0_m3": float(dominant["carrier_density_node0_m3"]),
            "carrier_density_node1_m3": float(dominant["carrier_density_node1_m3"]),
            "dominant_production_derivative_abs": production,
            "maximum_mobility_response_derivative_abs": mobility_response,
            "mobility_response_fraction_of_dominant_transport": (
                mobility_response / production if production else 0.0
            ),
            "minimum_continuity_row_weight": min(
                float(row["continuity_row_weight"]) for row in rows
            ),
            "maximum_continuity_row_weight": max(
                float(row["continuity_row_weight"]) for row in rows
            ),
            "constrained_row_record_count": sum(
                int(row["row_constrained"]) for row in rows
            ),
        })
    by_state = {
        (float(row["bias_V"]), str(row["variant"]), str(row["carrier"])): row
        for row in state_rows
    }
    changes: list[dict[str, object]] = []
    for bias in DEFAULT_BIASES:
        for carrier in CARRIERS:
            baseline = by_state[(bias, "vela_baseline", carrier)]
            joint = by_state[(bias, "sent_qfp_only", carrier)]
            b_weights = (
                float(baseline["minimum_continuity_row_weight"]),
                float(baseline["maximum_continuity_row_weight"]),
            )
            j_weights = (
                float(joint["minimum_continuity_row_weight"]),
                float(joint["maximum_continuity_row_weight"]),
            )
            changes.append({
                "bias_V": bias,
                "carrier": carrier,
                "edge_id": int(baseline["edge_id"]),
                "qfp_drive_joint_to_baseline_ratio": float(joint["qfp_drive_V_m"]) / float(baseline["qfp_drive_V_m"]),
                "mobility_joint_to_baseline_ratio": float(joint["mobility_m2_V_s"]) / float(baseline["mobility_m2_V_s"]),
                "dominant_transport_derivative_joint_to_baseline_ratio": float(joint["dominant_production_derivative_abs"]) / float(baseline["dominant_production_derivative_abs"]),
                "bernoulli_node0_absolute_change": abs(float(joint["bernoulli_node0"]) - float(baseline["bernoulli_node0"])),
                "bernoulli_node1_absolute_change": abs(float(joint["bernoulli_node1"]) - float(baseline["bernoulli_node1"])),
                "minimum_row_weight_joint_to_baseline_ratio": j_weights[0] / b_weights[0],
                "maximum_row_weight_joint_to_baseline_ratio": j_weights[1] / b_weights[1],
            })
    return state_rows, changes


def classify(
    decomposition: list[dict[str, object]],
    contact_rows: list[dict[str, object]],
    deterministic: bool,
) -> dict[str, object]:
    analytic_max = max(
        float(row["analytic_to_frozen_fd_edge_scaled_error"]) for row in decomposition
    )
    mobility_product_max = max(
        float(row["mobility_product_relative_error"]) for row in decomposition
    )
    row_scaling_max = max(
        float(row["row_scaling_closure_relative_error"]) for row in decomposition
    )
    bernoulli_max = max(
        abs(float(row["bernoulli_qfp_derivative_physical"])) for row in decomposition
    )
    contact_elimination_max = max(
        float(row["maximum_contact_eliminated_edge_derivative_abs"])
        for row in contact_rows
    )
    contact_identity_max = max(
        float(row["maximum_contact_identity_error_abs"]) for row in contact_rows
    )
    mobility_shares = [
        float(row["mobility_response_fraction_of_live_derivative"])
        for row in decomposition
    ]
    live_total_max = max(
        float(row["live_total_edge_scaled_error"]) for row in decomposition
    )
    passed = (
        analytic_max <= 5.0e-5
        and mobility_product_max <= 5.0e-5
        and live_total_max <= 5.0e-5
        and row_scaling_max <= 1.0e-12
        and bernoulli_max == 0.0
        and contact_elimination_max == 0.0
        and contact_identity_max == 0.0
        and deterministic
    )
    return {
        "passed": passed,
        "analytic_to_frozen_fd_max_edge_scaled_error": analytic_max,
        "mobility_product_max_relative_error": mobility_product_max,
        "live_total_max_edge_scaled_error": live_total_max,
        "row_scaling_max_relative_closure_error": row_scaling_max,
        "bernoulli_qfp_derivative_max_abs": bernoulli_max,
        "contact_eliminated_edge_derivative_max_abs": contact_elimination_max,
        "contact_identity_max_abs_error": contact_identity_max,
        "mobility_response_fraction_min": min(mobility_shares),
        "mobility_response_fraction_max": max(mobility_shares),
        "typed_outcome": (
            "transport_edge_decomposition_verified"
            if passed
            else "transport_edge_decomposition_contract_failed"
        ),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--base-config", type=Path, required=True)
    parser.add_argument("--prior-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--repeats", type=int, default=2)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.repeats < 2:
        raise RuntimeError("at least two repeats are required")
    args.output_root.mkdir(parents=True, exist_ok=True)
    hashes: dict[str, dict[str, str]] = {}
    raw: dict[str, dict[tuple[float, str], list[dict[str, str]]]] = {}
    for repeat in range(args.repeats):
        label = chr(ord("a") + repeat)
        hashes[label] = {}
        raw[label] = {}
        for bias in DEFAULT_BIASES:
            tag = bias_tag(bias)
            for variant in VARIANTS:
                case = args.output_root / f"run-{label}" / tag / variant
                csv_path = case / "transport_edge_jacobian.csv"
                config_path = case / "transport_edge_jacobian.json"
                config = make_probe_config(
                    args.base_config,
                    csv_path,
                    args.prior_root / "inputs" / tag / f"{variant}_fields",
                    "transport_edge_jacobian_probe",
                    bias,
                    "Anode",
                    "Cathode",
                )
                config["physical_finite_difference_step_V"] = FD_STEP_V
                write_json(config_path, config)
                write_json(case / "status.json", run_probe(args.runner, config_path))
                raw[label][(bias, variant)] = read_rows(csv_path)
                hashes[label][f"{tag}/{variant}"] = sha256(csv_path)

    hotspot_rows: list[dict[str, object]] = []
    decomposition_rows: list[dict[str, object]] = []
    contact_rows: list[dict[str, object]] = []
    prior_terms = read_rows(args.prior_root / "carrier_term_decomposition.csv")
    residual_hotspots: dict[tuple[float, str], int] = {}
    for row in prior_terms:
        if (
            row["variant"] == "sent_qfp_only"
            and row["scope"] == "interior"
            and row["term"] == "flux"
        ):
            residual_hotspots[(float(row["bias_V"]), row["carrier"])] = int(
                row["hotspot_node_id"]
            )
    for bias in DEFAULT_BIASES:
        selected_rows = select_hotspots(
            bias,
            raw["a"][(bias, "vela_baseline")],
            raw["a"][(bias, "sent_qfp_only")],
            {carrier: residual_hotspots[(bias, carrier)] for carrier in CARRIERS},
        )
        hotspot_rows.extend(selected_rows)
        selected = {row["carrier"]: int(row["edge_id"]) for row in selected_rows}
        for variant in VARIANTS:
            rows = raw["a"][(bias, variant)]
            decomposition_rows.extend(decompose_selected(bias, variant, rows, selected))
            contact_rows.append(contact_audit(bias, variant, rows))

    deterministic = all(
        len({hashes[label][key] for label in hashes}) == 1
        for key in hashes["a"]
    )
    determinism_rows = [
        {
            "artifact": key,
            "repeat_count": args.repeats,
            "unique_hash_count": len({hashes[label][key] for label in hashes}),
            "byte_identical": int(
                len({hashes[label][key] for label in hashes}) == 1
            ),
            "sha256": hashes["a"][key],
        }
        for key in sorted(hashes["a"])
    ]
    verdict = classify(decomposition_rows, contact_rows, deterministic)
    state_summary_rows, state_change_rows = summarize_hotspot_states(
        decomposition_rows
    )
    result = {
        "schema": "vela.pn2d_bv_m2_transport_edge_jacobian_verification.v1",
        "status": "passed" if verdict["passed"] else "failed",
        "biases_V": list(DEFAULT_BIASES),
        "variants": list(VARIANTS),
        "physical_finite_difference_step_V": FD_STEP_V,
        "hotspot_contract": "per_bias_per_carrier_incident_to_prior_interior_flux_residual_hotspot__max_abs_baseline_to_joint_qfp_edge_flux_delta__tie_lowest_edge_id",
        "production_contract": {
            "mobility_jacobian_field_derivatives": False,
            "meaning": "production transport Jacobian freezes high-field mobility",
            "contact_handling": "replace constrained rows only; constrained columns remain",
        },
        "deterministic": deterministic,
        "verdict": verdict,
        "outputs": {
            "hotspot_selection": "hotspot_selection.csv",
            "hotspot_decomposition": "hotspot_decomposition.csv",
            "contact_row_audit": "contact_row_audit.csv",
            "hotspot_state_summary": "hotspot_state_summary.csv",
            "hotspot_state_change": "hotspot_state_change.csv",
            "determinism": "determinism.csv",
        },
    }
    write_rows(args.output_root / "hotspot_selection.csv", hotspot_rows)
    write_rows(args.output_root / "hotspot_decomposition.csv", decomposition_rows)
    write_rows(args.output_root / "contact_row_audit.csv", contact_rows)
    write_rows(args.output_root / "hotspot_state_summary.csv", state_summary_rows)
    write_rows(args.output_root / "hotspot_state_change.csv", state_change_rows)
    write_rows(args.output_root / "determinism.csv", determinism_rows)
    write_json(args.output_root / "result.json", result)
    print(json.dumps(result, sort_keys=True))
    return 0 if verdict["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
