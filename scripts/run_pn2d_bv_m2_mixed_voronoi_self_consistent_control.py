#!/usr/bin/env python3
"""Run the prospective M2 mixed-Voronoi off -> IIC -> on control gates."""

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


BRANCHES = ("avalanche_off", "iic_postprocess", "avalanche_on")
ATOMIC_DEFAULT_PROFILE = {
    "impact_ionization": {
        "current_approximation": "element_edge_sg_gss_laux",
        "source_mapping_mode": "element_vertex_box_measure",
        "cell_reconstructed_midpoint_density": "bernoulli",
    },
    "mesh_geometry": {
        "node_volume_policy": "mixed_voronoi",
        "require_non_obtuse": True,
    },
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def prepare_mixed_config(base_path: Path, output_path: Path) -> dict[str, Any]:
    config = load_json(base_path)
    geometry = config.setdefault("mesh_geometry", {})
    geometry["node_volume_policy"] = "mixed_voronoi"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(
        json.dumps(config, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return config


def validate_actual_default_render(
    config_path: Path, manifest_path: Path
) -> dict[str, Any]:
    config = load_json(config_path)
    manifest = load_json(manifest_path)
    impact = config.get("solver", {}).get("impact_ionization", {})
    observed = {
        "impact_ionization": {
            name: impact.get(name)
            for name in ATOMIC_DEFAULT_PROFILE["impact_ionization"]
        },
        "mesh_geometry": {
            name: config.get("mesh_geometry", {}).get(name)
            for name in ATOMIC_DEFAULT_PROFILE["mesh_geometry"]
        },
    }
    gates = {
        "template_is_pn2d_bv": manifest.get("template") == "pn2d_bv",
        "template_version_at_least_3": int(manifest.get("template_version", 0)) >= 3,
        "default_render_has_no_profile_override": (
            "avalanche_current_support_profile" not in manifest.get("overrides", {})
        ),
        "default_profile_parameter_is_sg_laux": (
            manifest.get("parameters", {}).get("avalanche_current_support_profile")
            == "element_edge_sg_gss_laux"
        ),
        "resolved_profile_matches": (
            manifest.get("resolved_profile") == ATOMIC_DEFAULT_PROFILE
        ),
        "rendered_config_matches": observed == ATOMIC_DEFAULT_PROFILE,
    }
    if not all(gates.values()):
        failed = sorted(name for name, passed in gates.items() if not passed)
        raise ValueError(
            "actual-default render failed atomic binding checks: " + ", ".join(failed)
        )
    return {"gates": gates, "observed_profile": observed}


def terminal_current_rows(path: Path) -> dict[float, float]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return {
            round(float(row["bias_V"]), 9): float(row["current_total_A_per_um"])
            for row in csv.DictReader(stream)
        }


def sentaurus_terminal_current_rows(path: Path, branch: str) -> dict[float, float]:
    manifest = load_json(path)
    return {
        round(float(row["requested_bias_V"]), 9): float(row["value"])
        for row in manifest.get("aggregate_records", [])
        if row.get("branch") == branch
        and row.get("quantity") == "terminal_current"
        and row.get("carrier") == "total"
        and row.get("provenance") == "native"
    }


def state_hashes(execution: dict[str, Any], branch: str) -> dict[str, str]:
    record = next(item for item in execution["branches"] if item["branch"] == branch)
    return {
        key: value["sha256"]
        for key, value in record.get("state_files", {}).items()
    }


def branch_record(execution: dict[str, Any], branch: str) -> dict[str, Any]:
    return next(item for item in execution["branches"] if item["branch"] == branch)


def continuity_max(iv_path: Path) -> float:
    maximum = 0.0
    with iv_path.open(newline="", encoding="utf-8-sig") as stream:
        for row in csv.DictReader(stream):
            for name in (
                "global_electron_continuity_closure_ratio",
                "global_hole_continuity_closure_ratio",
            ):
                raw = row.get(name, "")
                if raw:
                    maximum = max(maximum, abs(float(raw)))
    return maximum


def determinism_metrics(run_a: dict[str, Any], run_b: dict[str, Any], branch: str) -> dict[str, Any]:
    a = branch_record(run_a, branch)
    b = branch_record(run_b, branch)
    states_a = state_hashes(run_a, branch)
    states_b = state_hashes(run_b, branch)
    return {
        "iv_sha256_equal": a["output_csv_sha256"] == b["output_csv_sha256"],
        "state_hashes_equal": states_a == states_b,
        "state_count_a": len(states_a),
        "state_count_b": len(states_b),
    }


def off_golden_metrics(vela_iv: Path, sentaurus_manifest: Path) -> dict[str, Any]:
    vela = terminal_current_rows(vela_iv)
    sentaurus = sentaurus_terminal_current_rows(sentaurus_manifest, "avalanche_off")
    biases = sorted((set(vela) & set(sentaurus)) - {0.0}, reverse=True)
    deltas = []
    for bias in biases:
        lhs = abs(vela[bias])
        rhs = abs(sentaurus[bias])
        if lhs <= 0.0 or rhs <= 0.0:
            raise ValueError(f"nonpositive current magnitude at {bias:g} V")
        deltas.append(math.log10(lhs) - math.log10(rhs))
    if not deltas:
        raise ValueError("no nonzero common off-branch biases")
    return {
        "compared_bias_count": len(biases),
        "compared_biases_V": biases,
        "log10_current_rmse_dex": math.sqrt(sum(value * value for value in deltas) / len(deltas)),
        "maximum_abs_log10_current_error_dex": max(abs(value) for value in deltas),
    }


def run_exact_lattice(args: argparse.Namespace, candidate_config: Path, branch: str, run_name: str) -> dict[str, Any]:
    run_root = args.output_root.resolve() / branch / run_name
    command = [
        sys.executable,
        str(args.exact_lattice_script.resolve()),
        "--runner", str(args.runner.resolve()),
        "--base-config", str(candidate_config.resolve()),
        "--sentaurus-manifest", str(args.sentaurus_manifest.resolve()),
        "--output-root", str(run_root),
        "--branches", branch,
        "--max-iter", str(args.max_iter),
    ]
    if args.base_config_manifest is not None:
        command.extend(
            ["--base-config-manifest", str(args.base_config_manifest.resolve())]
        )
    if not args.actual_default_render:
        command.append("--sg-laux-candidate")
    completed = subprocess.run(command, text=True, capture_output=True)
    (run_root / "gate_driver.stdout.log").write_text(completed.stdout, encoding="utf-8")
    (run_root / "gate_driver.stderr.log").write_text(completed.stderr, encoding="utf-8")
    if completed.returncode != 0:
        raise RuntimeError(
            f"{branch}/{run_name} failed with code {completed.returncode}; "
            f"see {run_root / 'gate_driver.stderr.log'}"
        )
    return load_json(run_root / "execution.json")


def evaluate_common(contract: dict[str, Any], run_a: dict[str, Any], run_b: dict[str, Any], branch: str) -> tuple[dict[str, Any], bool]:
    requirements = contract["common_requirements"]
    a = branch_record(run_a, branch)
    b = branch_record(run_b, branch)
    deterministic = determinism_metrics(run_a, run_b, branch)
    continuity = max(
        continuity_max(Path(a["output_csv"])),
        continuity_max(Path(b["output_csv"])),
    )
    metrics = {
        "run_a_complete_exact_lattice": bool(a["complete_exact_lattice"]),
        "run_b_complete_exact_lattice": bool(b["complete_exact_lattice"]),
        "run_a_observed_bias_count": int(a["observed_bias_count"]),
        "run_b_observed_bias_count": int(b["observed_bias_count"]),
        "maximum_global_continuity_closure_ratio": continuity,
        "determinism": deterministic,
        "current_support_origin": {
            "run_a": run_a.get("current_support", {}).get("origin"),
            "run_b": run_b.get("current_support", {}).get("origin"),
        },
        "base_config_manifest_sha256": {
            "run_a": (run_a.get("base_config_manifest") or {}).get("sha256"),
            "run_b": (run_b.get("base_config_manifest") or {}).get("sha256"),
        },
    }
    passed = (
        metrics["run_a_complete_exact_lattice"]
        and metrics["run_b_complete_exact_lattice"]
        and metrics["run_a_observed_bias_count"] == requirements["requested_bias_count"]
        and metrics["run_b_observed_bias_count"] == requirements["requested_bias_count"]
        and deterministic["iv_sha256_equal"]
        and deterministic["state_hashes_equal"]
        and continuity <= requirements["maximum_global_continuity_closure_ratio"]
    )
    if requirements.get("require_base_config_current_support_origin", False):
        passed = passed and all(
            origin == "base_config"
            for origin in metrics["current_support_origin"].values()
        )
    if requirements.get("require_bound_base_config_manifest", False):
        manifest_hashes = metrics["base_config_manifest_sha256"]
        passed = passed and bool(manifest_hashes["run_a"])
        passed = passed and manifest_hashes["run_a"] == manifest_hashes["run_b"]
    return metrics, passed


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", required=True, type=Path)
    parser.add_argument("--exact-lattice-script", required=True, type=Path)
    parser.add_argument("--base-config", required=True, type=Path)
    parser.add_argument("--base-config-manifest", type=Path)
    parser.add_argument("--sentaurus-manifest", required=True, type=Path)
    parser.add_argument("--contract", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    parser.add_argument("--max-iter", type=int, default=80)
    parser.add_argument(
        "--actual-default-render",
        action="store_true",
        help=(
            "Use the bound base config exactly as rendered, require a version-3 "
            "default manifest with no profile override, and do not add the SG/Laux "
            "CLI candidate flag."
        ),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.output_root = args.output_root.resolve()
    args.output_root.mkdir(parents=True, exist_ok=True)
    contract = load_json(args.contract.resolve())
    binding = None
    if args.actual_default_render:
        if args.base_config_manifest is None:
            raise ValueError("--actual-default-render requires --base-config-manifest")
        candidate_config = args.base_config.resolve()
        binding = validate_actual_default_render(
            candidate_config, args.base_config_manifest.resolve()
        )
        config = load_json(candidate_config)
        candidate_origin = "pn2d_bv_template_default_render"
    else:
        candidate_config = args.output_root / "inputs" / "simulation_mixed_voronoi.json"
        config = prepare_mixed_config(args.base_config.resolve(), candidate_config)
        candidate_origin = "cli_opt_in"
    report: dict[str, Any] = {
        "schema": "vela.pn2d_bv_m2_mixed_voronoi_self_consistent_control.v1",
        "status": "running",
        "outcome": None,
        "contract": {"path": str(args.contract.resolve()), "sha256": sha256(args.contract.resolve())},
        "base_config": {"path": str(args.base_config.resolve()), "sha256": sha256(args.base_config.resolve())},
        "candidate_config": {"path": str(candidate_config), "sha256": sha256(candidate_config)},
        "candidate_origin": candidate_origin,
        "default_render_binding": binding,
        "base_config_manifest": (
            {
                "path": str(args.base_config_manifest.resolve()),
                "sha256": sha256(args.base_config_manifest.resolve()),
            }
            if args.base_config_manifest is not None
            else None
        ),
        "node_volume_policy": config["mesh_geometry"]["node_volume_policy"],
        "require_non_obtuse": config["mesh_geometry"].get("require_non_obtuse", False),
        "sentaurus_manifest": {"path": str(args.sentaurus_manifest.resolve()), "sha256": sha256(args.sentaurus_manifest.resolve())},
        "stages": {},
    }

    executions: dict[str, tuple[dict[str, Any], dict[str, Any]]] = {}
    for branch in BRANCHES:
        run_a = run_exact_lattice(args, candidate_config, branch, "run-a")
        run_b = run_exact_lattice(args, candidate_config, branch, "run-b")
        executions[branch] = (run_a, run_b)
        common, passed = evaluate_common(contract, run_a, run_b, branch)
        stage: dict[str, Any] = {"common": common}
        if branch == "avalanche_off":
            golden = off_golden_metrics(
                Path(branch_record(run_a, branch)["output_csv"]),
                args.sentaurus_manifest.resolve(),
            )
            gate = contract["avalanche_off_gate"]
            stage["sentaurus_golden"] = golden
            passed = passed and golden["log10_current_rmse_dex"] <= gate["maximum_log10_current_rmse_dex"]
            passed = passed and golden["maximum_abs_log10_current_error_dex"] <= gate["maximum_log10_current_error_dex"]
        elif branch == "iic_postprocess":
            off_a, _ = executions["avalanche_off"]
            iic_record = branch_record(run_a, branch)
            off_record = branch_record(off_a, "avalanche_off")
            equality = {
                "iv_sha256_equal_to_off": iic_record["output_csv_sha256"] == off_record["output_csv_sha256"],
                "state_hashes_equal_to_off": state_hashes(run_a, branch) == state_hashes(off_a, "avalanche_off"),
            }
            stage["off_equivalence"] = equality
            passed = passed and all(equality.values())
        stage["passed"] = passed
        report["stages"][branch] = stage
        if not passed:
            report["status"] = "failed"
            report["outcome"] = f"stopped_at_{branch}_gate"
            break
    else:
        report["status"] = "passed"
        report["outcome"] = "completed_all_stages"

    report_path = args.output_root / "gate_report.json"
    report_path.write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report, indent=2, sort_keys=True, allow_nan=False))
    return 0 if report["status"] == "passed" else 2


if __name__ == "__main__":
    raise SystemExit(main())
