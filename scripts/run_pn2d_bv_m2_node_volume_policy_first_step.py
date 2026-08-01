#!/usr/bin/env python3
"""Compare barycentric and mixed-Voronoi M2 frozen first Newton updates."""

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


BIASES = (-18.0, -19.5, -19.7, -20.0)
VARIANTS = ("vela_baseline", "sent_qfp_only")
POLICIES = ("barycentric", "mixed_voronoi")
BLOCK_SUFFIXES = (
    "_summary.json", "_columns.csv", "_singular_modes.csv",
    "_solve_variants.csv", "_solve_nodes.csv",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--base-config", type=Path, required=True)
    parser.add_argument("--prior-root", type=Path, required=True)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--repeats", type=int, default=2)
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise RuntimeError(f"refusing to write empty CSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def sha256(paths: list[Path]) -> str:
    digest = hashlib.sha256()
    for path in paths:
        digest.update(path.name.encode())
        digest.update(path.read_bytes())
    return digest.hexdigest()


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


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


def probe_config(
    base_config: Path,
    fields_dir: Path,
    simulation_type: str,
    output: Path,
    bias: float,
    policy: str,
) -> dict[str, Any]:
    config = make_probe_config(
        base_config, output, fields_dir, simulation_type,
        bias, "Anode", "Cathode",
    )
    config["mesh_geometry"] = {"node_volume_policy": policy}
    return config


def norm(values: list[float]) -> float:
    return math.sqrt(sum(value * value for value in values))


def cosine(left: list[float], right: list[float]) -> float:
    denominator = norm(left) * norm(right)
    return sum(a * b for a, b in zip(left, right)) / denominator if denominator else 1.0


def relative_change(candidate: float, baseline: float) -> float:
    return abs(candidate - baseline) / max(abs(baseline), 1.0e-300)


def summarize_case(case: Path, bias: float, variant: str, policy: str) -> dict[str, Any]:
    step_status = json.loads((case / "newton_step_status.json").read_text(encoding="utf-8"))
    block_summary = json.loads((case / "carrier_block_summary.json").read_text(encoding="utf-8"))
    block_solves = {
        row["variant"]: row
        for row in read_csv(case / "carrier_block_solve_variants.csv")
    }
    carrier_solve_rows = [
        row for row in read_csv(case / "carrier_block_solve_nodes.csv")
        if row["variant"] == "full"
    ]
    carrier_vector = [
        float(row[field])
        for row in carrier_solve_rows
        for field in ("delta_phin_V", "delta_phip_V")
    ]
    modes = sorted(
        read_csv(case / "carrier_block_singular_modes.csv"),
        key=lambda row: float(row["step_energy_fraction"]),
        reverse=True,
    )
    step_rows = read_csv(case / "newton_step.csv")
    step_vector = [
        float(row[field])
        for row in step_rows
        for field in ("delta_psi_V", "delta_phin_V", "delta_phip_V")
    ]
    term_rows = read_csv(case / "carrier_terms.csv")
    electron_source = sum(float(row["impact_electron_source"]) for row in term_rows)
    hole_source = sum(float(row["impact_hole_source"]) for row in term_rows)
    combined_source = sum(float(row["impact_combined_source"]) for row in term_rows)
    initial_carrier = math.hypot(
        float(step_status["block_residuals"]["phin"]),
        float(step_status["block_residuals"]["phip"]),
    )
    trial_carrier = math.hypot(
        float(step_status["trial_block_residuals"]["phin"]),
        float(step_status["trial_block_residuals"]["phip"]),
    )
    row: dict[str, Any] = {
        "bias_V": bias,
        "variant": variant,
        "node_volume_policy": policy,
        "initial_psi_residual": step_status["block_residuals"]["psi"],
        "initial_phin_residual": step_status["block_residuals"]["phin"],
        "initial_phip_residual": step_status["block_residuals"]["phip"],
        "initial_carrier_residual": initial_carrier,
        "initial_combined_residual": step_status["block_residuals"]["combined"],
        "trial_psi_residual": step_status["trial_block_residuals"]["psi"],
        "trial_phin_residual": step_status["trial_block_residuals"]["phin"],
        "trial_phip_residual": step_status["trial_block_residuals"]["phip"],
        "trial_carrier_residual": trial_carrier,
        "trial_combined_residual": step_status["trial_block_residuals"]["combined"],
        "full_first_step_norm": norm(step_vector),
        "carrier_only_step_norm_V": float(block_solves["full"]["physical_step_norm_V"]),
        "l2_equilibrated_condition_number": block_summary["l2_equilibrated_condition"]["resolved_condition_number"],
        "impact_electron_source": electron_source,
        "impact_hole_source": hole_source,
        "impact_combined_source": combined_source,
        "step_vector": step_vector,
        "carrier_vector": carrier_vector,
    }
    for rank, mode in enumerate(modes[:2], start=1):
        sigma = float(mode["singular_value"])
        row[f"mode{rank}_index"] = int(mode["mode_index"])
        row[f"mode{rank}_step_energy_fraction"] = float(mode["step_energy_fraction"])
        row[f"mode{rank}_relative_singular_value"] = float(mode["relative_singular_value"])
        row[f"mode{rank}_transport_over_sigma"] = float(mode["transport_jacobian_projection"]) / sigma
        row[f"mode{rank}_avalanche_diagonal_over_sigma"] = float(mode["avalanche_diagonal_jacobian_projection"]) / sigma
        row[f"mode{rank}_avalanche_cross_over_sigma"] = float(mode["avalanche_cross_jacobian_projection"]) / sigma
        row[f"mode{rank}_recombination_over_sigma"] = float(mode["recombination_jacobian_projection"]) / sigma
    return row


def main() -> int:
    args = parse_args()
    if args.repeats < 2:
        raise RuntimeError("at least two independent repeats are required")
    contract = json.loads(args.contract.read_text(encoding="utf-8"))
    if tuple(float(value) for value in contract["biases_V"]) != BIASES:
        raise RuntimeError("contract bias lattice does not match the runner")
    args.output_root.mkdir(parents=True, exist_ok=True)
    hashes: dict[str, dict[str, str]] = {}
    first_cases: dict[tuple[float, str, str], Path] = {}

    for repeat in range(args.repeats):
        label = chr(ord("a") + repeat)
        hashes[label] = {}
        for policy in POLICIES:
            for bias in BIASES:
                tag = bias_tag(bias)
                for variant in VARIANTS:
                    case = args.output_root / f"run-{label}" / policy / tag / variant
                    fields = args.prior_root / "inputs" / tag / f"{variant}_fields"
                    case.mkdir(parents=True, exist_ok=True)

                    step_config = probe_config(
                        args.base_config, fields, "newton_step_probe",
                        case / "newton_step.csv", bias, policy,
                    )
                    write_json(case / "newton_step.json", step_config)
                    write_json(
                        case / "newton_step_status.json",
                        run_probe(args.runner, case / "newton_step.json"),
                    )

                    block_prefix = case / "carrier_block"
                    block_config = probe_config(
                        args.base_config, fields,
                        "newton_carrier_block_decomposition_probe",
                        case / "unused.csv", bias, policy,
                    )
                    block_config.pop("output_csv", None)
                    block_config["output_prefix"] = str(block_prefix.resolve())
                    write_json(case / "carrier_block.json", block_config)
                    write_json(
                        case / "carrier_block_status.json",
                        run_probe(args.runner, case / "carrier_block.json"),
                    )

                    term_config = probe_config(
                        args.base_config, fields, "newton_carrier_term_probe",
                        case / "carrier_terms.csv", bias, policy,
                    )
                    write_json(case / "carrier_terms.json", term_config)
                    write_json(
                        case / "carrier_terms_status.json",
                        run_probe(args.runner, case / "carrier_terms.json"),
                    )

                    artifacts = [
                        case / "newton_step.csv",
                        case / "carrier_terms.csv",
                        *(Path(str(block_prefix) + suffix) for suffix in BLOCK_SUFFIXES),
                    ]
                    key = f"{policy}/{tag}/{variant}"
                    hashes[label][key] = sha256(artifacts)
                    if repeat == 0:
                        first_cases[(bias, variant, policy)] = case

    deterministic = all(
        len({hashes[label][key] for label in hashes}) == 1
        for key in hashes["a"]
    )
    cases = [
        summarize_case(first_cases[(bias, variant, policy)], bias, variant, policy)
        for bias in BIASES
        for variant in VARIANTS
        for policy in POLICIES
    ]
    by_key = {
        (float(row["bias_V"]), str(row["variant"]), str(row["node_volume_policy"])): row
        for row in cases
    }
    comparisons: list[dict[str, Any]] = []
    decision_changes: list[float] = []
    for bias in BIASES:
        for variant in VARIANTS:
            bary = by_key[(bias, variant, "barycentric")]
            mixed = by_key[(bias, variant, "mixed_voronoi")]
            changes = {
                "full_first_step_relative_change": relative_change(
                    float(mixed["full_first_step_norm"]), float(bary["full_first_step_norm"])
                ),
                "carrier_only_step_relative_change": relative_change(
                    float(mixed["carrier_only_step_norm_V"]), float(bary["carrier_only_step_norm_V"])
                ),
                "initial_carrier_residual_relative_change": relative_change(
                    float(mixed["initial_carrier_residual"]), float(bary["initial_carrier_residual"])
                ),
                "integrated_impact_source_relative_change": relative_change(
                    float(mixed["impact_combined_source"]), float(bary["impact_combined_source"])
                ),
            }
            decision_changes.extend(changes.values())
            comparisons.append({
                "bias_V": bias,
                "variant": variant,
                **changes,
                "full_first_step_direction_cosine": cosine(
                    mixed["step_vector"], bary["step_vector"]
                ),
                "carrier_only_step_direction_cosine": cosine(
                    mixed["carrier_vector"], bary["carrier_vector"]
                ),
                "initial_psi_residual_ratio": float(mixed["initial_psi_residual"]) /
                    max(float(bary["initial_psi_residual"]), 1.0e-300),
                "trial_combined_residual_ratio": float(mixed["trial_combined_residual"]) /
                    max(float(bary["trial_combined_residual"]), 1.0e-300),
                "l2_condition_number_ratio": float(mixed["l2_equilibrated_condition_number"]) /
                    max(float(bary["l2_equilibrated_condition_number"]), 1.0e-300),
                "mode1_transport_over_sigma_barycentric": bary["mode1_transport_over_sigma"],
                "mode1_transport_over_sigma_mixed_voronoi": mixed["mode1_transport_over_sigma"],
                "mode1_avalanche_sum_over_sigma_barycentric":
                    float(bary["mode1_avalanche_diagonal_over_sigma"]) +
                    float(bary["mode1_avalanche_cross_over_sigma"]),
                "mode1_avalanche_sum_over_sigma_mixed_voronoi":
                    float(mixed["mode1_avalanche_diagonal_over_sigma"]) +
                    float(mixed["mode1_avalanche_cross_over_sigma"]),
            })

    maximum_change = max(decision_changes)
    negligible = float(contract["classification"]["negligible_max_relative_change"])
    material = float(contract["classification"]["material_min_relative_change"])
    classification = (
        "material_node_volume_policy_sensitivity"
        if maximum_change >= material
        else "negligible_node_volume_policy_sensitivity"
        if maximum_change <= negligible
        else "intermediate_node_volume_policy_sensitivity"
    )
    passed = deterministic
    determinism_rows = [{
        "artifact": key,
        "repeat_count": args.repeats,
        "unique_hash_count": len({hashes[label][key] for label in hashes}),
        "byte_identical": int(len({hashes[label][key] for label in hashes}) == 1),
        "sha256": hashes["a"][key],
    } for key in sorted(hashes["a"])]

    csv_cases = [{
        key: value for key, value in row.items()
        if key not in {"step_vector", "carrier_vector"}
    } for row in cases]
    write_csv(args.output_root / "case_summary.csv", csv_cases)
    write_csv(args.output_root / "policy_comparison.csv", comparisons)
    write_csv(args.output_root / "determinism.csv", determinism_rows)
    result = {
        "schema": "vela.pn2d_bv_m2_node_volume_policy_first_step.v1",
        "status": "passed" if passed else "failed",
        "typed_outcome": classification if passed else "determinism_failed",
        "biases_V": list(BIASES),
        "variants": list(VARIANTS),
        "policies": list(POLICIES),
        "contract_sha256": file_sha256(args.contract),
        "verdict": {
            "passed": passed,
            "deterministic": deterministic,
            "maximum_decision_metric_relative_change": maximum_change,
            "maximum_full_first_step_relative_change": max(
                float(row["full_first_step_relative_change"]) for row in comparisons
            ),
            "maximum_carrier_only_step_relative_change": max(
                float(row["carrier_only_step_relative_change"]) for row in comparisons
            ),
            "maximum_initial_carrier_residual_relative_change": max(
                float(row["initial_carrier_residual_relative_change"]) for row in comparisons
            ),
            "maximum_integrated_impact_source_relative_change": max(
                float(row["integrated_impact_source_relative_change"]) for row in comparisons
            ),
            "classification": classification,
            "production_defaults_modified": False,
            "state_advanced": False,
            "doping_redistributed": False,
        },
        "outputs": {
            "case_summary": "case_summary.csv",
            "policy_comparison": "policy_comparison.csv",
            "determinism": "determinism.csv",
        },
    }
    write_json(args.output_root / "result.json", result)
    print(json.dumps(result, sort_keys=True))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
