#!/usr/bin/env python3
"""Validate two independent normalized PN2D Sentaurus process matrices."""

from __future__ import annotations

import argparse
import json
import sys
from collections.abc import Mapping
from pathlib import Path
from typing import Any

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.pn2d_bv_process_contract import (
    EXACT_BIAS_TOLERANCE_V,
    validate_process_run,
)

REQUIRED_BRANCHES = frozenset(
    {
        "avalanche_off",
        "iic_postprocess",
        "avalanche_on",
        "avalanche_on_aval_derivatives",
    }
)
SOURCE_CLOSURE_RELATIVE_TOLERANCE = 1.0e-10

STATE_QUANTITIES = frozenset(
    {
        "potential",
        "density",
        "quasi_fermi",
        "electric_field",
        "quasi_fermi_gradient",
        "mobility",
        "velocity",
        "current_density",
        "doping",
        "charge_density",
        "srh_recombination",
    }
)


def process_matrix_shape(manifest: Mapping[str, Any]) -> dict[str, Any]:
    biases_by_branch = {
        record["branch"]: tuple(float(value) for value in record["requested_biases_V"])
        for record in manifest["branch_records"]
    }
    branch_names = frozenset(biases_by_branch)
    reference_biases = (
        biases_by_branch.get("avalanche_off", ())
        if branch_names == REQUIRED_BRANCHES
        else ()
    )
    common_bias_lattice = bool(reference_biases) and all(
        biases == reference_biases for biases in biases_by_branch.values()
    )
    snapshot_count = sum(
        len(record["bias_records"]) for record in manifest["branch_records"]
    )
    return {
        "required_branches_present": branch_names == REQUIRED_BRANCHES,
        "common_bias_lattice": common_bias_lattice,
        "biases_V": list(reference_biases),
        "snapshot_count": snapshot_count,
        "expected_snapshot_count": len(REQUIRED_BRANCHES) * len(reference_biases),
    }


def load_run(root: Path) -> dict[str, Any]:
    manifest = json.loads((root / "manifest.json").read_text(encoding="ascii"))
    validate_process_run(manifest, base_dir=root)
    return manifest


def state_index(
    manifest: Mapping[str, Any],
    branch: str,
) -> dict[tuple[Any, ...], tuple[float, ...]]:
    result: dict[tuple[Any, ...], tuple[float, ...]] = {}
    for record in manifest["field_records"]:
        if record["branch"] != branch or record["quantity"] not in STATE_QUANTITIES:
            continue
        key = (
            float(record["actual_bias_V"]),
            record["support_kind"],
            record["support_key"],
            record["provenance"],
            record["carrier"],
            record["quantity"],
            tuple(record["components"]),
        )
        result[key] = tuple(float(value) for value in record["values"])
    return result


def compare_states(
    manifest: Mapping[str, Any],
    left: str,
    right: str,
) -> dict[str, Any]:
    left_rows = state_index(manifest, left)
    right_rows = state_index(manifest, right)
    if left_rows.keys() != right_rows.keys():
        raise ValueError(
            f"{left}/{right}: state support mismatch "
            f"left_only={len(left_rows.keys() - right_rows.keys())} "
            f"right_only={len(right_rows.keys() - left_rows.keys())}"
        )
    maximum_absolute = 0.0
    maximum_relative = 0.0
    worst: dict[str, Any] | None = None
    for key in left_rows:
        for component, (left_value, right_value) in enumerate(
            zip(left_rows[key], right_rows[key], strict=True)
        ):
            absolute = abs(left_value - right_value)
            relative = absolute / max(
                abs(left_value),
                abs(right_value),
                1.0e-300,
            )
            maximum_absolute = max(maximum_absolute, absolute)
            if relative > maximum_relative:
                maximum_relative = relative
                worst = {
                    "key": list(key[:-1]) + [list(key[-1])],
                    "component": component,
                    "left": left_value,
                    "right": right_value,
                    "absolute": absolute,
                    "relative": relative,
                }
    return {
        "left": left,
        "right": right,
        "record_count": len(left_rows),
        "maximum_absolute": maximum_absolute,
        "maximum_relative": maximum_relative,
        "worst": worst,
    }


def source_closure(
    manifest: Mapping[str, Any],
) -> dict[str, Any]:
    grouped: dict[tuple[str, float, str], dict[str, float]] = {}
    for record in manifest["aggregate_records"]:
        if record["quantity"] != "integrated_source":
            continue
        key = (
            record["branch"],
            float(record["actual_bias_V"]),
            record["carrier"],
        )
        grouped.setdefault(key, {})[record["provenance"]] = float(record["value"])
    maximum_relative = 0.0
    worst: dict[str, Any] | None = None
    compared = 0
    incomplete = 0
    for key, values in grouped.items():
        if "native" not in values or "operator_replay" not in values:
            incomplete += 1
            continue
        compared += 1
        native = values["native"]
        replay = values["operator_replay"]
        relative = abs(native - replay) / max(abs(native), abs(replay), 1.0e-300)
        if relative > maximum_relative:
            maximum_relative = relative
            worst = {
                "branch": key[0],
                "actual_bias_V": key[1],
                "carrier": key[2],
                "native_A_per_um": native,
                "operator_replay_A_per_um": replay,
                "relative": relative,
            }
    return {
        "grouped_records": len(grouped),
        "compared_records": compared,
        "incomplete_records": incomplete,
        "maximum_relative": maximum_relative,
        "worst": worst,
    }


def max_iic_generation(manifest: Mapping[str, Any]) -> float:
    return max(
        (
            abs(float(value))
            for record in manifest["field_records"]
            if record["branch"] == "iic_postprocess"
            and record["quantity"] == "avalanche_generation"
            for value in record["values"]
        ),
        default=0.0,
    )


def analyze(root_a: Path, root_b: Path) -> dict[str, Any]:
    manifest_a = load_run(root_a)
    manifest_b = load_run(root_b)
    shape_a = process_matrix_shape(manifest_a)
    shape_b = process_matrix_shape(manifest_b)
    shape_equal = shape_a == shape_b
    normalized_equal = (
        manifest_a["normalized_output_hashes"]
        == manifest_b["normalized_output_hashes"]
    )
    input_equal = manifest_a["input_hashes"] == manifest_b["input_hashes"]
    exact_bias_error = max(
        abs(
            float(record["requested_bias_V"])
            - float(record["actual_bias_V"])
        )
        for branch in manifest_a["branch_records"]
        for record in branch["bias_records"]
    )
    snapshot_count = int(shape_a["snapshot_count"])
    iic_state = compare_states(
        manifest_a,
        "avalanche_off",
        "iic_postprocess",
    )
    derivatives_state = compare_states(
        manifest_a,
        "avalanche_on",
        "avalanche_on_aval_derivatives",
    )
    iic_generation = max_iic_generation(manifest_a)
    closure = source_closure(manifest_a)
    expected_source_comparisons = snapshot_count * 3

    outcome = "sentaurus_process_matrix_available"
    if (
        not normalized_equal
        or not input_equal
        or not shape_equal
        or not shape_a["required_branches_present"]
        or not shape_a["common_bias_lattice"]
        or snapshot_count != shape_a["expected_snapshot_count"]
        or exact_bias_error > EXACT_BIAS_TOLERANCE_V
        or closure["compared_records"] != expected_source_comparisons
        or closure["incomplete_records"] != 0
        or closure["maximum_relative"] > SOURCE_CLOSURE_RELATIVE_TOLERANCE
    ):
        outcome = "exact_snapshot_mismatch"
    elif (
        iic_state["maximum_relative"] > 1.0e-10
        and iic_state["maximum_absolute"] > 1.0e-12
    ) or iic_generation <= 0.0:
        outcome = "iic_state_not_decoupled"
    elif (
        derivatives_state["maximum_relative"] > 1.0e-10
        and derivatives_state["maximum_absolute"] > 1.0e-12
    ):
        outcome = "sentaurus_solver_path_difference"

    return {
        "schema": "vela.pn2d_bv_process_matrix_pair_acceptance.v1",
        "outcome": outcome,
        "root_a": str(root_a.resolve()),
        "root_b": str(root_b.resolve()),
        "input_hashes_identical": input_equal,
        "normalized_output_hashes_identical": normalized_equal,
        "matrix_shape_identical": shape_equal,
        "matrix_shape": shape_a,
        "exact_bias_maximum_error_V": exact_bias_error,
        "snapshot_count": snapshot_count,
        "expected_source_comparisons": expected_source_comparisons,
        "iic_maximum_generation_cm^-3_s^-1": iic_generation,
        "avalanche_off_to_iic_state": iic_state,
        "avalanche_on_derivative_control_state": derivatives_state,
        "currentplot_to_tcl_source_closure": closure,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root-a", type=Path, required=True)
    parser.add_argument("--root-b", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = analyze(args.root_a.resolve(), args.root_b.resolve())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n",
        encoding="ascii",
        newline="\n",
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["outcome"] == "sentaurus_process_matrix_available" else 1


if __name__ == "__main__":
    raise SystemExit(main())
