#!/usr/bin/env python3
"""Run the PN2D Poisson-QFP Jacobian cross-block diagnostic at exact states."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path
from typing import Any

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.diagnose_pn2d_bv_predictor_first_step_audit import (
    make_probe_config,
    state_csv_to_fields,
)
from scripts.run_pn2d_bv_feedback_state_substitution import (
    bias_token,
    build_replacement_fields,
    exact_state_record,
    mesh_coordinates,
    read_csv,
    sha256,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sentaurus-chain", type=Path, required=True)
    parser.add_argument("--vela-manifest", type=Path, required=True)
    parser.add_argument("--base-config", type=Path, required=True)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--biases", nargs="+", type=float, default=[-19.7, -19.8])
    parser.add_argument("--branch", default="avalanche_on")
    return parser.parse_args()


def run_case(
    *,
    chain: dict[str, Any],
    manifest: dict[str, Any],
    manifest_root: Path,
    branch: str,
    bias: float,
    base_config: Path,
    runner: Path,
    output: Path,
) -> dict[str, Any]:
    record = exact_state_record(manifest, branch, bias)
    state_path = manifest_root / str(record["snapshot_tdr"]["path"])
    baseline_rows = read_csv(state_path)
    coordinates = mesh_coordinates(base_config)
    if len(coordinates) != len(baseline_rows):
        raise ValueError(f"{state_path}: state/mesh node-count mismatch")

    case_root = output / branch / bias_token(bias)
    baseline_fields = case_root / "baseline_state_fields"
    replacement_fields = case_root / "feedback_state_fields"
    output_csv = case_root / "poisson_qfp_cross_block.csv"
    blocks_csv = case_root / "jacobian_blocks.csv"
    config_path = case_root / "poisson_qfp_cross_block.json"
    status_path = case_root / "status.json"
    case_root.mkdir(parents=True, exist_ok=True)
    state_csv_to_fields(state_path, baseline_fields)
    artifact_hashes = build_replacement_fields(
        chain,
        branch,
        bias,
        baseline_rows,
        coordinates,
        replacement_fields,
    )

    config = make_probe_config(
        base_config,
        output_csv,
        baseline_fields,
        "newton_poisson_qfp_cross_block_probe",
        bias,
        "Anode",
        "Cathode",
    )
    config["feedback_state_fields_dir"] = str(replacement_fields.resolve())
    config["jacobian_blocks_csv"] = str(blocks_csv.resolve())
    config_path.write_text(
        json.dumps(config, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    completed = subprocess.run(
        [str(runner), "--config", str(config_path)],
        cwd=case_root,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"{branch} {bias:g} V: "
            f"{completed.stderr.strip() or completed.stdout.strip()}"
        )
    status = json.loads(completed.stdout)
    status_path.write_text(
        json.dumps(status, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    for path in (
        state_path,
        config_path,
        output_csv,
        blocks_csv,
        status_path,
    ):
        artifact_hashes[str(path)] = sha256(path)
    return {
        "branch": branch,
        "bias_V": bias,
        "baseline_state": str(state_path),
        "config": str(config_path),
        "output_csv": str(output_csv),
        "jacobian_blocks_csv": str(blocks_csv),
        "status": str(status_path),
        "artifact_hashes": artifact_hashes,
    }


def main() -> int:
    args = parse_args()
    sentaurus_chain_path = args.sentaurus_chain.resolve()
    vela_manifest_path = args.vela_manifest.resolve()
    base_config = args.base_config.resolve()
    runner = args.runner.resolve()
    output = args.output_root.resolve()
    output.mkdir(parents=True, exist_ok=True)
    chain = json.loads(sentaurus_chain_path.read_text(encoding="utf-8"))
    manifest = json.loads(vela_manifest_path.read_text(encoding="utf-8"))
    cases = [
        run_case(
            chain=chain,
            manifest=manifest,
            manifest_root=vela_manifest_path.parent,
            branch=args.branch,
            bias=float(bias),
            base_config=base_config,
            runner=runner,
            output=output,
        )
        for bias in args.biases
    ]
    execution = {
        "schema": "vela.pn2d_bv_poisson_qfp_cross_block_execution.v1",
        "status": "passed",
        "outcome": "poisson_qfp_cross_block_observations_available",
        "sentaurus_chain": str(sentaurus_chain_path),
        "sentaurus_chain_sha256": sha256(sentaurus_chain_path),
        "vela_manifest": str(vela_manifest_path),
        "vela_manifest_sha256": sha256(vela_manifest_path),
        "base_config": str(base_config),
        "base_config_sha256": sha256(base_config),
        "runner": str(runner),
        "runner_sha256": sha256(runner),
        "contract": {
            "residual": "qfp_only_frozen_substitution",
            "jacobian": "single_production_baseline_jacobian",
            "blocks": [
                "J_psi_psi",
                "J_psi_qfp",
                "J_qfp_psi",
                "J_qfp_qfp",
            ],
            "counterfactuals": [
                "independent_blocks",
                "remove_J_psi_qfp",
                "remove_J_qfp_psi",
                "full_schur",
            ],
            "production_defaults_changed": False,
        },
        "cases": cases,
    }
    execution_path = output / "execution.json"
    execution_path.write_text(
        json.dumps(execution, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(execution, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
