#!/usr/bin/env python
"""Run the branch-resolved PN2D avalanche source-Jacobian audit.

The audit reuses a frozen ``newton_jacobian_block_probe`` configuration and
changes only the numerical-reference mode, finite-difference step, and output
path.  The production Jacobian remains the shared local-forward-AD Jacobian.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import subprocess
import sys
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
EXE = ".exe" if sys.platform.startswith("win") else ""
DEFAULT_BASE_CONFIG = (
    REPO
    / "build-release"
    / "pn2d-task12-direct-source-20260727"
    / "coarse7x3_m20V"
    / "config_1em10.json"
)
DEFAULT_RUNNER = REPO / "build-release" / f"vela_example_runner{EXE}"
DEFAULT_OUT_DIR = (
    REPO / "build-release" / "pn2d-wp6-branch-resolved-jacobian-20260729"
)
DEFAULT_STEPS = (1.0e-14, 3.0e-15, 1.0e-15)
DEFAULT_DOUBLE_STEPS = (1.0e-8, 1.0e-10, 3.0e-11)
RELATIVE_GATE = 1.0e-8
ABSOLUTE_GATE = 1.0e-12


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def resolve(path: Path) -> Path:
    return path.resolve() if path.is_absolute() else (REPO / path).resolve()


def step_token(step: float) -> str:
    return f"{step:.0e}".replace("+", "").replace("-", "m")


def input_paths(config: dict[str, object], config_path: Path) -> dict[str, Path]:
    base = config_path.parent
    result: dict[str, Path] = {}
    for key in ("mesh_file", "node_doping_file", "materials_file", "state_file"):
        value = config.get(key)
        if not isinstance(value, str) or not value:
            raise ValueError(f"base config is missing {key}")
        path = Path(value)
        result[key] = path.resolve() if path.is_absolute() else (base / path).resolve()
    return result


def run_probe(
    runner: Path,
    base_config: dict[str, object],
    out_dir: Path,
    step: float,
    reference_mode: str,
) -> dict[str, object]:
    token = step_token(step)
    mode_token = (
        "mp_branch_resolved"
        if reference_mode == "multiprecision_branch_resolved"
        else "double_symmetric"
    )
    output_csv = out_dir / f"jacobian_{mode_token}_{token}.csv"
    config_path = out_dir / f"config_{mode_token}_{token}.json"
    config = json.loads(json.dumps(base_config))
    config["blocks"] = ["sg_avalanche"]
    config["finite_difference_step"] = step
    config["finite_difference_mode"] = reference_mode
    config["output_csv"] = str(output_csv)
    config_path.write_text(json.dumps(config, indent=2), encoding="utf-8")
    environment = os.environ.copy()
    if sys.platform.startswith("win"):
        ucrt_bins = (Path("D:/msys64/ucrt64/bin"), Path("D:/msys64/usr/bin"))
        existing = environment.get("PATH", "")
        environment["PATH"] = os.pathsep.join(
            [str(path) for path in ucrt_bins if path.is_dir()] + [existing]
        )
    completed = subprocess.run(
        [str(runner), "--config", str(config_path)],
        check=True,
        cwd=REPO,
        capture_output=True,
        text=True,
        env=environment,
    )
    runner_payload = json.loads(completed.stdout)
    with output_csv.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    source_rows = [row for row in rows if row["block"] == "sg_avalanche"]
    if len(source_rows) != 1:
        raise RuntimeError(f"{output_csv} does not contain one sg_avalanche row")
    row = source_rows[0]
    analytic_norm = float(row["analytic_norm"])
    reference_norm = float(row["fd_norm"])
    difference_norm = float(row["diff_norm"])
    reference_scale = max(analytic_norm, reference_norm)
    nonzero = reference_scale > ABSOLUTE_GATE
    true_relative = (
        difference_norm / reference_scale if reference_scale > 0.0 else 0.0
    )
    passed = (
        true_relative <= RELATIVE_GATE
        if nonzero
        else difference_norm <= ABSOLUTE_GATE
    )
    return {
        "relative_step": step,
        "reference_mode": reference_mode,
        "analytic_norm": analytic_norm,
        "reference_norm": reference_norm,
        "difference_norm": difference_norm,
        "true_relative_difference": true_relative,
        "classification": "nonzero" if nonzero else "near_zero",
        "passed": passed,
        "config": str(config_path),
        "config_sha256": sha256(config_path),
        "output_csv": str(output_csv),
        "output_csv_sha256": sha256(output_csv),
        "configuration_fingerprint": runner_payload.get(
            "impact_configuration_fingerprint"
        ),
        "active_branch_fingerprint": runner_payload.get(
            "impact_active_branch_fingerprint"
        ),
    }


def write_summary_csv(path: Path, records: list[dict[str, object]]) -> None:
    fields = [
        "relative_step",
        "reference_mode",
        "analytic_norm",
        "reference_norm",
        "difference_norm",
        "true_relative_difference",
        "classification",
        "passed",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for record in records:
            writer.writerow({field: record[field] for field in fields})


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-config", type=Path, default=DEFAULT_BASE_CONFIG)
    parser.add_argument("--runner", type=Path, default=DEFAULT_RUNNER)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--step", action="append", type=float)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    base_config_path = resolve(args.base_config)
    runner = resolve(args.runner)
    out_dir = resolve(args.out_dir)
    steps = tuple(args.step) if args.step else DEFAULT_STEPS
    if len(steps) < 3 or any(step <= 0.0 or not math.isfinite(step) for step in steps):
        raise ValueError("at least three positive finite --step values are required")
    base_config = json.loads(base_config_path.read_text(encoding="utf-8-sig"))
    paths = input_paths(base_config, base_config_path)
    for label, path in {"runner": runner, **paths}.items():
        if not path.is_file():
            raise FileNotFoundError(f"{label} not found: {path}")

    out_dir.mkdir(parents=True, exist_ok=True)
    branch_resolved_records = [
        run_probe(
            runner,
            base_config,
            out_dir,
            step,
            "multiprecision_branch_resolved",
        )
        for step in steps
    ]
    double_records = [
        run_probe(
            runner,
            base_config,
            out_dir,
            step,
            "double_symmetric",
        )
        for step in DEFAULT_DOUBLE_STEPS
    ]
    write_summary_csv(
        out_dir / "derivative_convergence.csv",
        double_records + branch_resolved_records,
    )
    input_hashes = {
        label: {"path": str(path), "sha256": sha256(path)}
        for label, path in paths.items()
    }
    all_records = double_records + branch_resolved_records
    fingerprints_agree = (
        len(
            {
                (
                    record["configuration_fingerprint"],
                    record["active_branch_fingerprint"],
                )
                for record in all_records
            }
        )
        == 1
        and all(
            bool(record["configuration_fingerprint"])
            and bool(record["active_branch_fingerprint"])
            for record in all_records
        )
    )
    summary = {
        "schema": "vela.pn2d_bv_branch_resolved_jacobian_audit.v1",
        "outcome": (
            "source_jacobian_dependency_identified_and_closed"
            if (
                all(bool(record["passed"]) for record in branch_resolved_records)
                and fingerprints_agree
            )
            else "nonsmooth_branch_derivative"
        ),
        "reference_mode": "multiprecision_branch_resolved",
        "zero_branch_policy": "symmetric_semismooth_zero_derivative",
        "nonzero_branch_policy": "step_below_active_branch_margin",
        "relative_gate": RELATIVE_GATE,
        "absolute_gate": ABSOLUTE_GATE,
        "base_config": {
            "path": str(base_config_path),
            "sha256": sha256(base_config_path),
        },
        "runner": {"path": str(runner), "sha256": sha256(runner)},
        "inputs": input_hashes,
        "configuration_and_active_branch_fingerprints_agree": fingerprints_agree,
        "double_symmetric_records": double_records,
        "branch_resolved_records": branch_resolved_records,
    }
    summary_path = out_dir / "summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0 if summary["outcome"] == "source_jacobian_dependency_identified_and_closed" else 2


if __name__ == "__main__":
    raise SystemExit(main())
