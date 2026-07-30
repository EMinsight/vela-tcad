#!/usr/bin/env python3
"""Evaluate Vela fixed-transition residuals and first Newton corrections."""

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
from scripts.run_pn2d_bv_sentaurus_newton_probe_vm import (
    BRANCHES,
    EXACT_TOLERANCE_V,
    KNEE_BIASES,
    bias_token,
    parse_branches,
    transition_sources,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-manifest", type=Path, required=True)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--branches", default=",".join(BRANCHES))
    parser.add_argument("--biases", nargs="+", type=float, default=list(KNEE_BIASES))
    parser.add_argument("--resume", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    branches = parse_branches(args.branches)
    targets = tuple(float(value) for value in args.biases)
    source_manifest_path = args.source_manifest.resolve()
    source_root = source_manifest_path.parent
    output = args.output_root.resolve()
    output.mkdir(parents=True, exist_ok=True)
    runner = args.runner.resolve()
    manifest = json.loads(source_manifest_path.read_text(encoding="utf-8"))

    cases: list[dict[str, Any]] = []
    for branch in branches:
        base_config = source_root / branch / "simulation.json"
        for target, predecessor in transition_sources(manifest, branch, targets):
            predecessor_bias = float(predecessor["requested_bias_V"])
            if abs(predecessor_bias - float(predecessor["actual_bias_V"])) > EXACT_TOLERANCE_V:
                raise ValueError(f"{branch}: predecessor is not exact")
            source_state = source_root / predecessor["snapshot_tdr"]["path"]
            token = bias_token(target)
            case = output / branch / token
            fields = case / "state_fields"
            config_path = case / "newton_step_probe.json"
            csv_path = case / "newton_step_probe.csv"
            status_path = case / "status.json"
            case.mkdir(parents=True, exist_ok=True)
            state_csv_to_fields(source_state, fields)
            config = make_probe_config(
                base_config,
                csv_path,
                fields,
                "newton_step_probe",
                target,
                "Anode",
                "Cathode",
            )
            config_path.write_text(
                json.dumps(config, indent=2, sort_keys=True) + "\n",
                encoding="utf-8",
            )
            if not (args.resume and csv_path.is_file() and status_path.is_file()):
                completed = subprocess.run(
                    [str(runner), "--config", str(config_path)],
                    cwd=case,
                    text=True,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    check=False,
                )
                if completed.returncode != 0:
                    raise RuntimeError(
                        f"{branch} {target:g} V: "
                        f"{completed.stderr.strip() or completed.stdout.strip()}"
                    )
                status = json.loads(completed.stdout)
                status_path.write_text(
                    json.dumps(status, indent=2, sort_keys=True) + "\n",
                    encoding="utf-8",
                )
            cases.append(
                {
                    "branch": branch,
                    "target_bias_V": target,
                    "predecessor_bias_V": predecessor_bias,
                    "source_state": str(source_state),
                    "config": str(config_path),
                    "csv": str(csv_path),
                    "status": str(status_path),
                }
            )

    execution = {
        "schema": "vela.pn2d_bv_vela_newton_probe_execution.v1",
        "status": "passed",
        "outcome": "fixed_transition_first_newton_observations_available",
        "source_manifest": str(source_manifest_path),
        "runner": str(runner),
        "cases": cases,
    }
    (output / "execution.json").write_text(
        json.dumps(execution, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(execution, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
