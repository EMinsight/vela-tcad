#!/usr/bin/env python3
"""Evaluate prospective M0/M2 evidence for the PN2D BV template default."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from pathlib import Path
from typing import Any, Iterable, Mapping

from generate_pn2d_config import render_named_template


EXACT_TOLERANCE_V = 1.0e-10
BRANCHES = ("avalanche_off", "iic_postprocess", "avalanche_on")


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def artifact(path: Path) -> dict[str, str]:
    resolved = path.resolve()
    return {"path": str(resolved), "sha256": sha256(resolved)}


def percentile(values: Iterable[float], fraction: float) -> float:
    ordered = sorted(float(value) for value in values)
    rank = fraction * (len(ordered) - 1)
    lower = math.floor(rank)
    upper = math.ceil(rank)
    if lower == upper:
        return ordered[lower]
    weight = rank - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def exact_rows(
    rows: Iterable[Mapping[str, Any]], biases: Iterable[float]
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
                f"expected one exact row at {target:g} V, found {len(matches)}"
            )
        result[float(target)] = matches[0]
    return result


def abs_log_error(left: float, right: float) -> float:
    if (
        left == 0.0
        or right == 0.0
        or not math.isfinite(left)
        or not math.isfinite(right)
    ):
        raise ValueError("log-error inputs must be finite and non-zero")
    return abs(math.log10(abs(left)) - math.log10(abs(right)))


def branch_index(execution: Mapping[str, Any]) -> dict[str, Mapping[str, Any]]:
    result = {str(row["branch"]): row for row in execution["branches"]}
    if set(result) != set(BRANCHES):
        raise ValueError(f"execution branches differ: {sorted(result)}")
    return result


def state_evidence(
    manifest: Mapping[str, Any],
    manifest_path: Path,
    expected_biases: tuple[float, ...],
) -> dict[str, Any]:
    result: dict[tuple[str, float], str] = {}
    requested_ok = True
    actual_hashes_ok = True
    for branch in manifest.get("branch_records", []):
        name = str(branch["branch"])
        requested = tuple(float(value) for value in branch["requested_biases_V"])
        requested_ok = requested_ok and requested == expected_biases
        for record in branch.get("bias_records", []):
            key = (name, float(record["requested_bias_V"]))
            if key in result:
                raise ValueError(f"duplicate state record: {key}")
            snapshot = record["snapshot_tdr"]
            expected_hash = str(snapshot["sha256"])
            state_path = manifest_path.parent / str(snapshot["path"])
            actual_hashes_ok = (
                actual_hashes_ok
                and state_path.is_file()
                and sha256(state_path) == expected_hash
            )
            result[key] = expected_hash
    expected_keys = {
        (branch, bias) for branch in BRANCHES for bias in expected_biases
    }
    gates = {
        "manifest_passed": manifest.get("status") == "passed",
        "branch_set_exact": {
            str(branch["branch"])
            for branch in manifest.get("branch_records", [])
        }
        == set(BRANCHES),
        "requested_biases_exact": requested_ok,
        "snapshot_key_set_exact": set(result) == expected_keys,
        "snapshot_file_hashes_valid": actual_hashes_ok and bool(result),
    }
    return {"passed": all(gates.values()), "gates": gates, "hashes": result}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def process_probe_evidence(
    path: Path, biases: tuple[float, ...]
) -> dict[str, Any]:
    required_columns = {
        "bias_V",
        "configuration_fingerprint",
        "source_integral",
        "electron_residual_contributions",
        "hole_residual_contributions",
    }
    observed_biases: set[float] = set()
    columns: set[str] = set()
    if path.is_file():
        with path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            columns = set(reader.fieldnames or ())
            if required_columns <= columns:
                for row in reader:
                    bias = float(row["bias_V"])
                    if any(
                        abs(bias - target) <= EXACT_TOLERANCE_V
                        for target in biases
                    ):
                        observed_biases.add(bias)
    coverage = all(
        any(abs(observed - target) <= EXACT_TOLERANCE_V for observed in observed_biases)
        for target in biases
    )
    gates = {
        "file_present": path.is_file(),
        "required_columns_present": required_columns <= columns,
        "exact_bias_coverage": coverage,
    }
    return {
        "passed": all(gates.values()),
        "gates": gates,
        "observed_contract_bias_count": len(observed_biases),
        "artifact": artifact(path) if path.is_file() else None,
    }


def closure_evidence(
    iv_path: Path,
    process_probe_path: Path,
    biases: tuple[float, ...],
    contract: Mapping[str, Any],
) -> dict[str, Any]:
    rows = exact_rows(read_csv(iv_path), biases)
    required = (
        "global_continuity_closure_satisfied",
        "global_electron_continuity_closure_ratio",
        "global_hole_continuity_closure_ratio",
    )
    missing = sorted(set(required) - set(next(iter(rows.values()))))
    electron = []
    hole = []
    satisfied = []
    if not missing:
        for row in rows.values():
            satisfied.append(row["global_continuity_closure_satisfied"] == "1")
            electron.append(float(row["global_electron_continuity_closure_ratio"]))
            hole.append(float(row["global_hole_continuity_closure_ratio"]))
    gates = {
        "columns_present": not missing,
        "all_satisfied": bool(satisfied) and all(satisfied),
        "electron_ratio": (
            bool(electron)
            and max(electron)
            <= float(contract["maximum_electron_closure_ratio"])
        ),
        "hole_ratio": (
            bool(hole)
            and max(hole) <= float(contract["maximum_hole_closure_ratio"])
        ),
        "process_probe_bound": process_probe_path.is_file(),
    }
    return {
        "passed": all(gates.values()),
        "gates": gates,
        "missing_columns": missing,
        "maximum_electron_closure_ratio": max(electron) if electron else None,
        "maximum_hole_closure_ratio": max(hole) if hole else None,
        "process_probe": (
            artifact(process_probe_path) if process_probe_path.is_file() else None
        ),
    }


def parity_evidence(
    parity: Mapping[str, Any],
    biases: tuple[float, ...],
    thresholds: Mapping[str, Any],
) -> dict[str, Any]:
    rows = exact_rows(parity["curve_rows"], biases)
    curve_errors = [
        abs_log_error(
            float(row["vela_on_A_per_um"]),
            float(row["sentaurus_on_A_per_um"]),
        )
        for row in rows.values()
    ]
    gain_errors = [
        abs_log_error(float(row["vela_gain"]), float(row["sentaurus_gain"]))
        for row in rows.values()
    ]
    knee = parity["knee_metrics"]
    estimators = parity["knee_estimators"]
    vela_slope = estimators["vela"].get("V_slope")
    sentaurus_slope = estimators["sentaurus"].get("V_slope")
    slope_error = (
        None
        if vela_slope is None or sentaurus_slope is None
        else abs(float(vela_slope) - float(sentaurus_slope))
    )
    break_error = abs(
        float(estimators["vela"]["V_break"])
        - float(estimators["sentaurus"]["V_break"])
    )
    metrics = {
        "effective_curve_median_abs_log_error_dex": statistics.median(
            curve_errors
        ),
        "effective_curve_p95_abs_log_error_dex": percentile(curve_errors, 0.95),
        "effective_curve_max_abs_log_error_dex": max(curve_errors),
        "effective_gain_median_abs_log_error_dex": statistics.median(gain_errors),
        "effective_gain_max_abs_log_error_dex": max(gain_errors),
        "knee_median_abs_log_error_dex": float(
            knee["median_absolute_log_error_dex"]
        ),
        "knee_max_abs_log_error_dex": float(
            knee["maximum_absolute_log_error_dex"]
        ),
        "V_break_abs_error_V": break_error,
        "V_slope_abs_error_V": slope_error,
        "adjacent_slope_rmse_dex_per_V": float(
            parity["adjacent_slope_rmse_dex_per_V"]
        ),
    }
    if any(
        value is not None and not math.isfinite(value)
        for value in metrics.values()
    ):
        raise ValueError("contract-domain parity metrics must be finite")
    gates = {
        name: metrics[name] is not None and metrics[name] <= float(limit)
        for name, limit in thresholds.items()
    }
    return {"passed": all(gates.values()), "metrics": metrics, "gates": gates}


def evaluate_level(
    name: str,
    contract: Mapping[str, Any],
    execution_a_path: Path,
    execution_b_path: Path,
    state_a_path: Path,
    state_b_path: Path,
    parity_path: Path,
    render_manifest_path: Path,
    sentaurus_manifest_path: Path,
    sentaurus_on_csv_path: Path,
    sentaurus_off_csv_path: Path,
) -> dict[str, Any]:
    execution_a = load_json(execution_a_path)
    execution_b = load_json(execution_b_path)
    state_a = load_json(state_a_path)
    state_b = load_json(state_b_path)
    parity = load_json(parity_path)
    render_manifest = load_json(render_manifest_path)
    biases = tuple(float(value) for value in contract["bv_domain"]["exact_biases_V"])
    execution_biases_a = tuple(
        float(value) for value in execution_a["requested_biases_V"]
    )
    execution_biases_b = tuple(
        float(value) for value in execution_b["requested_biases_V"]
    )
    branch_a = branch_index(execution_a)
    branch_b = branch_index(execution_b)
    state_result_a = state_evidence(state_a, state_a_path, execution_biases_a)
    state_result_b = state_evidence(state_b, state_b_path, execution_biases_b)

    on_config_path = Path(str(branch_a["avalanche_on"]["config"]))
    on_config = load_json(on_config_path)
    impact = on_config["solver"]["impact_ionization"]
    required_profile = contract["required_default_profile"]
    atomic_fields = required_profile["atomic_fields"]
    config_gates = {
        f"atomic:{key}": impact.get(key) == expected
        for key, expected in atomic_fields.items()
    }
    config_gates.update(
        {
            "mobility_basis": (
                on_config["solver"]["mobility"].get(
                    "doping_concentration_basis"
                )
                == required_profile["mobility_doping_concentration_basis"]
            ),
            "render_profile": (
                render_manifest["parameters"].get(
                    required_profile["profile_parameter"]
                )
                == required_profile["profile_value"]
            ),
            "template_version_2_or_newer": (
                int(render_manifest["template_version"]) >= 2
            ),
            "default_path_not_cli_opt_in": (
                execution_a.get("current_support", {}).get("origin")
                == "base_config"
                and execution_b.get("current_support", {}).get("origin")
                == "base_config"
            ),
            "execution_current_support_matches_profile": all(
                execution.get("current_support", {}).get(key) == expected
                for execution in (execution_a, execution_b)
                for key, expected in atomic_fields.items()
            ),
        }
    )

    complete_gates = {}
    deterministic_gates = {}
    for branch in BRANCHES:
        complete_gates[f"{branch}:run_a"] = (
            int(branch_a[branch]["returncode"]) == 0
            and bool(branch_a[branch]["complete_exact_lattice"])
        )
        complete_gates[f"{branch}:run_b"] = (
            int(branch_b[branch]["returncode"]) == 0
            and bool(branch_b[branch]["complete_exact_lattice"])
        )
        deterministic_gates[f"{branch}:iv"] = (
            branch_a[branch]["output_csv_sha256"]
            == branch_b[branch]["output_csv_sha256"]
        )
        deterministic_gates[f"{branch}:physics_config"] = (
            branch_a[branch]["physics_config_sha256"]
            == branch_b[branch]["physics_config_sha256"]
        )
    deterministic_gates["execution_bias_lattice"] = (
        execution_biases_a == execution_biases_b
    )
    deterministic_gates["all_state_snapshots"] = (
        state_result_a["passed"]
        and state_result_b["passed"]
        and state_result_a["hashes"] == state_result_b["hashes"]
    )
    probe_a = Path(str(branch_a["avalanche_on"]["config"])).parent / "process_probe.csv"
    probe_b = Path(str(branch_b["avalanche_on"]["config"])).parent / "process_probe.csv"
    probe_result_a = process_probe_evidence(probe_a, biases)
    probe_result_b = process_probe_evidence(probe_b, biases)
    deterministic_gates["avalanche_on:process_probe"] = (
        probe_result_a["passed"]
        and probe_result_b["passed"]
        and probe_result_a["artifact"]["sha256"]
        == probe_result_b["artifact"]["sha256"]
    )

    closure = closure_evidence(
        Path(str(branch_a["avalanche_on"]["output_csv"])),
        on_config_path.parent / "process_probe.csv",
        biases,
        contract["closure"],
    )
    closure["gates"]["process_probe_bound"] = probe_result_a["passed"]
    closure["passed"] = all(closure["gates"].values())
    closure["process_probe_evidence"] = probe_result_a
    parity_result = parity_evidence(
        parity, biases, contract["bv_domain"]["thresholds"]
    )

    sentaurus_hash = sha256(sentaurus_manifest_path)
    parity_inputs = parity.get("inputs", {})
    curve_input_bindings = {}
    expected_curve_paths = {
        "vela_on": Path(str(branch_a["avalanche_on"]["output_csv"])).resolve(),
        "vela_off": Path(str(branch_a["avalanche_off"]["output_csv"])).resolve(),
        "sentaurus_on": sentaurus_on_csv_path.resolve(),
        "sentaurus_off": sentaurus_off_csv_path.resolve(),
    }
    curve_input_hashes_valid = set(parity_inputs) == set(expected_curve_paths)
    curve_input_paths_valid = curve_input_hashes_valid
    if curve_input_hashes_valid:
        for input_name, expected_path in expected_curve_paths.items():
            record = parity_inputs[input_name]
            path = Path(str(record["path"])).resolve()
            curve_input_paths_valid = (
                curve_input_paths_valid and path == expected_path
            )
            curve_input_bindings[input_name] = artifact(path)
            curve_input_hashes_valid = (
                curve_input_hashes_valid
                and curve_input_bindings[input_name]["sha256"]
                == record["sha256"]
            )

    branch_artifact_hashes_valid = True
    for branch in BRANCHES:
        for execution_branch in (branch_a[branch], branch_b[branch]):
            branch_artifact_hashes_valid = (
                branch_artifact_hashes_valid
                and sha256(Path(str(execution_branch["config"])))
                == execution_branch["config_sha256"]
                and sha256(Path(str(execution_branch["output_csv"])))
                == execution_branch["output_csv_sha256"]
            )

    render_manifest_hash = sha256(render_manifest_path)
    render_manifest_execution_binding = all(
        execution.get("base_config_manifest", {}).get("sha256")
        == render_manifest_hash
        and Path(
            str(execution.get("base_config_manifest", {}).get("path"))
        ).resolve()
        == render_manifest_path.resolve()
        for execution in (execution_a, execution_b)
    )
    expected_base_config, expected_render_manifest = render_named_template(
        render_manifest["template"],
        render_manifest["overrides"],
        allow_absolute_paths=True,
    )
    rendered_base_config_matches_manifest = (
        expected_render_manifest == render_manifest
        and load_json(Path(str(execution_a["base_config"])))
        == expected_base_config
        and load_json(Path(str(execution_b["base_config"])))
        == expected_base_config
    )
    required_binding_names = set(
        contract["execution_evidence"]["required_hash_bindings"]
    )
    supported_binding_names = {
        "contract",
        "pn2d_bv_template",
        "render_manifest",
        "base_config",
        "branch_configs",
        "mesh",
        "doping",
        "materials",
        "vela_iv",
        "vela_state_manifest",
        "sentaurus_manifest",
        "sentaurus_on_off_aggregate",
        "curve_acceptance",
    }
    binding_gates = {
        "sentaurus_manifest_matches_run_a": (
            sentaurus_hash == execution_a["sentaurus_manifest_sha256"]
        ),
        "sentaurus_manifest_matches_run_b": (
            sentaurus_hash == execution_b["sentaurus_manifest_sha256"]
        ),
        "base_config_matches_render": (
            sha256(Path(str(execution_a["base_config"])))
            == execution_a["base_config_sha256"]
            and sha256(Path(str(execution_b["base_config"])))
            == execution_b["base_config_sha256"]
        ),
        "render_manifest_matches_executions": render_manifest_execution_binding,
        "base_config_matches_render_manifest": rendered_base_config_matches_manifest,
        "branch_config_and_iv_hashes_valid": branch_artifact_hashes_valid,
        "curve_input_paths_valid": curve_input_paths_valid,
        "curve_input_hashes_valid": curve_input_hashes_valid,
        "required_hash_bindings_supported": (
            required_binding_names <= supported_binding_names
        ),
    }
    physical_inputs = {}
    for key in ("mesh_file", "node_doping_file", "materials_file"):
        path = Path(str(on_config[key]))
        physical_inputs[key] = artifact(path)

    gates = {
        "config_profile": all(config_gates.values()),
        "complete_exact_lattice": all(complete_gates.values()),
        "duplicate_determinism": all(deterministic_gates.values()),
        "global_continuity_closure": closure["passed"],
        "same_grid_sentaurus_parity": parity_result["passed"],
        "artifact_binding": all(binding_gates.values()),
    }
    return {
        "schema": "vela.pn2d_bv_template_default_level_acceptance.v1",
        "level": name,
        "status": "passed" if all(gates.values()) else "failed",
        "gates": gates,
        "configuration": {"gates": config_gates, "observed": impact},
        "complete_exact_lattice": complete_gates,
        "duplicate_determinism": deterministic_gates,
        "state_evidence": {
            "run_a": {
                "passed": state_result_a["passed"],
                "gates": state_result_a["gates"],
            },
            "run_b": {
                "passed": state_result_b["passed"],
                "gates": state_result_b["gates"],
            },
        },
        "process_probe_evidence": {
            "run_a": probe_result_a,
            "run_b": probe_result_b,
        },
        "closure": closure,
        "same_grid_sentaurus_parity": parity_result,
        "binding_checks": binding_gates,
        "artifact_bindings": {
            "execution_a": artifact(execution_a_path),
            "execution_b": artifact(execution_b_path),
            "state_manifest_a": artifact(state_a_path),
            "state_manifest_b": artifact(state_b_path),
            "curve_acceptance": artifact(parity_path),
            "render_manifest": artifact(render_manifest_path),
            "sentaurus_manifest": artifact(sentaurus_manifest_path),
            "sentaurus_on_aggregate": artifact(sentaurus_on_csv_path),
            "sentaurus_off_aggregate": artifact(sentaurus_off_csv_path),
            "base_config_a": artifact(Path(str(execution_a["base_config"]))),
            "base_config_b": artifact(Path(str(execution_b["base_config"]))),
            "branch_configs_a": {
                branch: artifact(Path(str(branch_a[branch]["config"])))
                for branch in BRANCHES
            },
            "branch_configs_b": {
                branch: artifact(Path(str(branch_b[branch]["config"])))
                for branch in BRANCHES
            },
            "vela_iv_a": {
                branch: artifact(Path(str(branch_a[branch]["output_csv"])))
                for branch in BRANCHES
            },
            "vela_iv_b": {
                branch: artifact(Path(str(branch_b[branch]["output_csv"])))
                for branch in BRANCHES
            },
            "curve_inputs": curve_input_bindings,
            "process_probe_a": artifact(probe_a),
            "process_probe_b": artifact(probe_b),
            "physical_inputs": physical_inputs,
        },
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--template", type=Path, required=True)
    parser.add_argument("--low-current-audit", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    for level in ("m0", "m2"):
        parser.add_argument(f"--{level}-execution-a", type=Path, required=True)
        parser.add_argument(f"--{level}-execution-b", type=Path, required=True)
        parser.add_argument(f"--{level}-state-a", type=Path, required=True)
        parser.add_argument(f"--{level}-state-b", type=Path, required=True)
        parser.add_argument(f"--{level}-parity", type=Path, required=True)
        parser.add_argument(f"--{level}-render-manifest", type=Path, required=True)
        parser.add_argument(
            f"--{level}-sentaurus-manifest", type=Path, required=True
        )
        parser.add_argument(f"--{level}-sentaurus-on-csv", type=Path, required=True)
        parser.add_argument(f"--{level}-sentaurus-off-csv", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    contract = load_json(args.contract)
    if (
        contract.get("schema")
        != "vela.pn2d_bv_m0_m2_template_default_acceptance_contract.v1"
        or not contract.get("prospective_only")
        or not contract.get("retroactive_threshold_mutation_forbidden")
    ):
        raise ValueError("invalid prospective M0/M2 template-default contract")
    low_current = load_json(args.low_current_audit)
    low_passed = (
        low_current.get("outcome")
        == contract["low_current_domain"]["required_typed_outcome"]
    )
    levels = {
        level.upper(): evaluate_level(
            level.upper(),
            contract,
            getattr(args, f"{level}_execution_a"),
            getattr(args, f"{level}_execution_b"),
            getattr(args, f"{level}_state_a"),
            getattr(args, f"{level}_state_b"),
            getattr(args, f"{level}_parity"),
            getattr(args, f"{level}_render_manifest"),
            getattr(args, f"{level}_sentaurus_manifest"),
            getattr(args, f"{level}_sentaurus_on_csv"),
            getattr(args, f"{level}_sentaurus_off_csv"),
        )
        for level in ("m0", "m2")
    }
    all_levels_passed = all(
        result["status"] == "passed" for result in levels.values()
    )
    accepted = all_levels_passed and low_passed
    policy = contract["decision_policy"]
    result = {
        "schema": "vela.pn2d_bv_m0_m2_template_default_acceptance.v1",
        "status": "passed" if accepted else "failed",
        "decision": (
            policy["all_levels_pass"] if accepted else policy["any_level_fail"]
        ),
        "production_default_change_authorized": accepted,
        "authorized_surface": (
            policy["production_default_surface_if_passed"] if accepted else "none"
        ),
        "global_cpp_default_change_authorized": False,
        "levels": levels,
        "low_current_classification": {
            "passed": low_passed,
            "observed_outcome": low_current.get("outcome"),
            "artifact": artifact(args.low_current_audit),
        },
        "artifact_bindings": {
            "contract": artifact(args.contract),
            "pn2d_bv_template": artifact(args.template),
        },
        "cross_mesh_convergence": {
            "observation_only": True,
            "not_used_as_acceptance_gate": True,
        },
    }
    output = args.output_root.resolve()
    output.mkdir(parents=True, exist_ok=True)
    (output / "acceptance.json").write_text(
        json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    for level, level_result in levels.items():
        (output / f"{level}_acceptance.json").write_text(
            json.dumps(
                level_result, indent=2, sort_keys=True, allow_nan=False
            )
            + "\n",
            encoding="utf-8",
        )
    print(json.dumps(result, indent=2, sort_keys=True, allow_nan=False))
    return 0 if accepted else 1


if __name__ == "__main__":
    raise SystemExit(main())
