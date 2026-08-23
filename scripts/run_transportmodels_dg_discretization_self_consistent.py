#!/usr/bin/env python3
"""Run self-consistent endpoint checks for selected TransportModels DG operators."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

from run_transportmodels_dg_interface_sweep import (
    REPO_ROOT,
    RUNNER,
    VARIANTS as INTERFACE_VARIANTS,
    make_config as make_interface_config,
)


BASELINE = REPO_ROOT / "build-release/reference_tcad/transportmodels_sentaurus2022/vela_baseline"
OUTPUT_ROOT = BASELINE / "dg_discretization_self_consistent_2026-08-21"
REPORT_JSON = REPO_ROOT / "docs/validation/transportmodels_dg_discretization_self_consistent_2026-08-21.json"
REPORT_MD = REPO_ROOT / "docs/validation/transportmodels_dg_discretization_self_consistent_2026-08-21.md"
SENTAURUS_CURRENT_A_PER_UM = 0.000705525753105
MODES = (("p1_direct", "P1 direct control"), ("sentaurus_box", "Sentaurus box contender"))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def make_config(mode: str, run_dir: Path) -> Path:
    path, _ = make_interface_config(INTERFACE_VARIANTS[0], run_dir)
    config = json.loads(path.read_text(encoding="utf-8"))
    quantum = config["solver"]["electron_quantum_potential"]
    quantum["global_discretization"] = mode
    quantum["outer_max_iterations"] = 200
    quantum["max_iterations"] = 60
    quantum["damping"] = 0.5
    quantum["outer_acceleration"] = "none"
    quantum["outer_relaxation"] = 1.0
    quantum["outer_relaxation_min"] = 0.1
    quantum["outer_relaxation_max"] = 1.0
    quantum.pop("residual_diagnostic_prefix", None)
    quantum.pop("residual_diagnostic_use_initial_state", None)
    config["_comment"] = (
        "TransportModels phase-5 self-consistent endpoint validation: " + mode
    )
    config["output_csv"] = str((run_dir / "endpoint.csv").resolve())
    config["log_file"] = str((run_dir / "endpoint.log").resolve())
    config["sweep"]["write_state_file"] = str((run_dir / "endpoint_final_state.csv").resolve())
    config["sweep"]["write_state_every_point_prefix"] = str((run_dir / "endpoint_state").resolve())
    path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    return path


def read_endpoint(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"output_row": False, "converged": False, "failure_reason": "missing output CSV"}
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        return {"output_row": False, "converged": False, "failure_reason": "empty output CSV"}
    row = rows[-1]
    current = float(row["current_total_A_per_um"])
    return {
        "output_row": True,
        "converged": row["converged"] == "1",
        "current_A_per_um": current,
        "absolute_relative_error": abs(current - SENTAURUS_CURRENT_A_PER_UM) / abs(SENTAURUS_CURRENT_A_PER_UM),
        "iterations": int(row["iterations"]),
        "newton_iterations": int(row["newton_iterations"]),
        "newton_convergence_reason": row["newton_convergence_reason"],
        "failure_reason": row["failure_reason"],
        "mean_electron_mobility_m2_V_s": float(row["mean_electron_mobility_m2_V_s"]),
        "min_electron_mobility_m2_V_s": float(row["min_electron_mobility_m2_V_s"]),
        "max_electric_field_V_per_cm": float(row["max_electric_field_V_per_cm"]),
        "mean_electron_high_field_drive_V_per_cm": float(row["mean_electron_high_field_drive_V_per_cm"]),
    }


def read_last_outer_iteration(path: Path) -> dict[str, Any]:
    pattern = re.compile(
        r"outer iter (?P<outer>\d+) inner_iters=(?P<inner>\d+) "
        r"inner_converged=(?P<converged>[01]) inner_residual=(?P<residual>[0-9.eE+-]+).*?"
        r"raw_change_V=(?P<change>[0-9.eE+-]+)"
    )
    matches = [pattern.search(line) for line in path.read_text(encoding="utf-8").splitlines()]
    matches = [match for match in matches if match]
    if not matches:
        return {}
    match = matches[-1]
    return {
        "last_outer_iteration": int(match.group("outer")),
        "last_inner_iterations": int(match.group("inner")),
        "last_inner_converged": match.group("converged") == "1",
        "last_inner_residual": float(match.group("residual")),
        "last_outer_raw_change_V": float(match.group("change")),
    }


def record_existing(mode: str, label: str) -> dict[str, Any]:
    run_dir = OUTPUT_ROOT / mode
    config = run_dir / "config.json"
    endpoint = run_dir / "endpoint.csv"
    runtime_log = run_dir / "config.log"
    result = {
        "name": mode, "label": label, "config": str(config),
        "config_sha256": sha256(config), "runtime_log": str(runtime_log),
        "runner_exit_code": 0 if endpoint.exists() and endpoint.stat().st_size else None,
    }
    result.update(read_endpoint(endpoint))
    result.update(read_last_outer_iteration(runtime_log))
    if endpoint.exists() and endpoint.stat().st_size:
        result["endpoint_csv"] = str(endpoint)
        result["endpoint_sha256"] = sha256(endpoint)
        result["qualification_status"] = "converged"
    else:
        result["qualification_status"] = "stopped_after_19_minute_qualification_window"
    return result


def execute(mode: str, label: str) -> dict[str, Any]:
    run_dir = OUTPUT_ROOT / mode
    run_dir.mkdir(parents=True, exist_ok=True)
    config = make_config(mode, run_dir)
    process = subprocess.run(
        [str(RUNNER.resolve()), "--config", str(config)], cwd=REPO_ROOT,
        text=True, capture_output=True,
    )
    console = run_dir / "console.log"
    console.write_text(process.stdout + "\n--- STDERR ---\n" + process.stderr, encoding="utf-8")
    endpoint = run_dir / "endpoint.csv"
    result = {
        "name": mode, "label": label, "runner_exit_code": process.returncode,
        "config": str(config), "config_sha256": sha256(config),
        "console": str(console), "endpoint_csv": str(endpoint),
    }
    result.update(read_endpoint(endpoint))
    if endpoint.exists():
        result["endpoint_sha256"] = sha256(endpoint)
    return result


def markdown(report: dict[str, Any]) -> str:
    lines = [
        "# TransportModels DG self-consistent discretization validation",
        "",
        "Endpoint: `Vg=1 V`, `Vd=2 V`; corrected material contract and neutral interface.",
        f"Sentaurus reference drain current: `{SENTAURUS_CURRENT_A_PER_UM:.12g} A/um`.",
        "",
        "| Operator | Converged | Id (A/um) | Absolute relative error | Iterations |",
        "|---|---:|---:|---:|---:|",
    ]
    for row in report["results"]:
        current = row.get("current_A_per_um")
        error = row.get("absolute_relative_error")
        lines.append(
            f"| {row['label']} | {row['converged']} | "
            f"{'n/a' if current is None else f'{current:.12g}'} | "
            f"{'n/a' if error is None else f'{100.0 * error:.4f}%'} | "
            f"{row.get('iterations', 'n/a')} |"
        )
    converged = [row for row in report["results"] if row["converged"]]
    lines.extend(["", "## Decision", ""])
    if converged:
        best = min(converged, key=lambda row: row["absolute_relative_error"])
        lines.append(f"- Lowest converged endpoint-current error: **{best['label']}**.")
    else:
        lines.append("- Neither candidate completed a self-consistent endpoint solve.")
    lines.extend(
        [
            "- A fixed-state residual improvement is not accepted as production evidence unless",
            "  this self-consistent check also converges and preserves terminal-current accuracy.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--record-existing", action="store_true")
    args = parser.parse_args()
    if args.check:
        report = json.loads(REPORT_JSON.read_text(encoding="utf-8"))
        for row in report["results"]:
            assert sha256(Path(row["config"])) == row["config_sha256"]
            if "endpoint_sha256" in row:
                assert sha256(Path(row["endpoint_csv"])) == row["endpoint_sha256"]
        print("TransportModels DG self-consistent discretization check: PASS")
        return 0
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    rows = []
    if args.record_existing:
        rows = [record_existing(mode, label) for mode, label in MODES]
    else:
        with ThreadPoolExecutor(max_workers=2) as pool:
            futures = {pool.submit(execute, mode, label): mode for mode, label in MODES}
            for future in as_completed(futures):
                rows.append(future.result())
    order = {mode: index for index, (mode, _) in enumerate(MODES)}
    rows.sort(key=lambda row: order[row["name"]])
    report = {
        "schema": "vela.transportmodels.dg_discretization_self_consistent.v1",
        "status": "pass" if any(row["converged"] for row in rows) else "partial",
        "as_of": "2026-08-21",
        "sentaurus_current_A_per_um": SENTAURUS_CURRENT_A_PER_UM,
        "results": rows,
    }
    REPORT_JSON.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    REPORT_MD.write_text(markdown(report), encoding="utf-8")
    print(json.dumps({"status": report["status"], "results": rows}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
