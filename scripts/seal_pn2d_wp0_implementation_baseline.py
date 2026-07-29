#!/usr/bin/env python3
"""Seal the current-code PN2D avalanche-off/on implementation baseline."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
from pathlib import Path
from typing import Any


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def normalized_on_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    result = []
    for row in rows:
        normalized = dict(row)
        diagnostic = normalized.get("newton_failure_diagnostics_json", "")
        if diagnostic:
            normalized["newton_failure_diagnostics_json"] = Path(diagnostic).name
        result.append(normalized)
    return result


def normalized_hash(rows: list[dict[str, str]]) -> str:
    payload = json.dumps(
        normalized_on_rows(rows),
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def on_summary(path: Path) -> dict[str, Any]:
    rows = read_csv(path)
    accepted = [float(row["bias_V"]) for row in rows if row["converged"] == "1"]
    sequence_payload = json.dumps(accepted, separators=(",", ":")).encode("ascii")
    global_biases = [float(-value) for value in range(21)]
    knee_biases = [
        -18.0,
        -18.5,
        -19.0,
        -19.25,
        -19.5,
        -19.7,
        -19.8,
        -19.85,
        -19.9,
        -19.95,
        -20.0,
    ]
    coverage = lambda targets: [  # noqa: E731 - compact report helper.
        target
        for target in targets
        if any(abs(value - target) <= 1.0e-10 for value in accepted)
    ]
    failed = next((row for row in rows if row["converged"] != "1"), None)
    return {
        "path": str(path.resolve()),
        "sha256": sha256(path),
        "normalized_sha256": normalized_hash(rows),
        "row_count": len(rows),
        "accepted_bias_count": len(accepted),
        "accepted_bias_sequence_sha256": hashlib.sha256(
            sequence_payload
        ).hexdigest(),
        "global_requested_biases_reached_V": coverage(global_biases),
        "knee_requested_biases_reached_V": coverage(knee_biases),
        "last_accepted_bias_V": accepted[-1],
        "first_failure": (
            {
                "requested_bias_V": float(failed["bias_V"]),
                "last_stable_bias_V": float(failed["last_stable_bias"]),
                "failed_bias_V": float(failed["failed_bias"]),
                "failure_reason": failed["failure_reason"],
                "newton_failure_class": failed["newton_failure_class"],
                "breakdown_failure_reason": failed["breakdown_failure_reason"],
            }
            if failed is not None
            else None
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--off-report", type=Path, required=True)
    parser.add_argument("--off-spatial-summary", type=Path, required=True)
    parser.add_argument("--on-a", type=Path, required=True)
    parser.add_argument("--on-b", type=Path, required=True)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--on-config", type=Path, required=True)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--doping", type=Path, required=True)
    parser.add_argument("--materials", type=Path, required=True)
    parser.add_argument("--sentaurus-curve", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    off = json.loads(args.off_report.read_text(encoding="ascii"))
    off_rows = read_csv(args.off_spatial_summary)
    log_errors = [
        float(row["log10_abs_current_ratio"])
        for row in off_rows
        if float(row["bias_V"]) != 0.0
    ]
    rmse = math.sqrt(sum(value * value for value in log_errors) / len(log_errors))
    maximum = max(abs(value) for value in log_errors)
    on_a = on_summary(args.on_a)
    on_b = on_summary(args.on_b)
    deterministic = (
        on_a["normalized_sha256"] == on_b["normalized_sha256"]
        and on_a["accepted_bias_sequence_sha256"]
        == on_b["accepted_bias_sequence_sha256"]
        and on_a["first_failure"] == on_b["first_failure"]
    )
    off_acceptance = off["acceptance"]
    off_passed = (
        rmse <= 0.001
        and maximum <= 0.002
        and off_acceptance["max_electron_closure_relative"] <= 1.0e-5
        and off_acceptance["max_hole_closure_relative"] <= 1.0e-5
        and off_acceptance["max_total_terminal_closure_A_per_um"] <= 1.0e-20
    )
    outcome = (
        "implementation_baseline_sealed"
        if off_passed and deterministic
        else "baseline_or_determinism_mismatch"
    )
    git_commit = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    inputs = {
        label: {"path": str(path.resolve()), "sha256": sha256(path)}
        for label, path in (
            ("executable", args.executable),
            ("on_config", args.on_config),
            ("mesh", args.mesh),
            ("doping", args.doping),
            ("materials", args.materials),
            ("sentaurus_curve", args.sentaurus_curve),
        )
    }
    result = {
        "schema": "vela.pn2d_wp0_implementation_baseline.v1",
        "outcome": outcome,
        "git_commit": git_commit,
        "inputs": inputs,
        "avalanche_off": {
            "log_current_rmse_dex": rmse,
            "maximum_log_current_error_dex": maximum,
            **off_acceptance,
        },
        "avalanche_on_run_a": on_a,
        "avalanche_on_run_b": on_b,
        "duplicate_runs_deterministic": deterministic,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n",
        encoding="ascii",
        newline="\n",
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if outcome == "implementation_baseline_sealed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
