#!/usr/bin/env python3
"""Run prospective 201-point PN2D forward-IV node-volume acceptance."""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any


ANCHORS = (1, 2, 5, 10, 15, 20)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def sentaurus_current(fields_root: Path, bias: int) -> float:
    rows = csv_rows(
        fields_root / f"{bias}v" / "fields" / "ContactCurrentFlux_region2.csv"
    )
    return abs(float(rows[0]["component0"]))


def prepare_config(base: dict[str, Any], case_dir: Path, policy: str) -> Path:
    config = copy.deepcopy(base)
    config.setdefault("mesh_geometry", {})["node_volume_policy"] = policy
    config["output_csv"] = str((case_dir / "iv.csv").resolve())
    sweep = config["sweep"]
    sweep["write_vtk"] = False
    sweep["vtk_prefix"] = str((case_dir / "state").resolve())
    sweep["write_state_file"] = str((case_dir / "last_state.csv").resolve())
    diagnostics = sweep.setdefault("diagnostics", {})
    if "newton_history" in diagnostics:
        diagnostics["newton_history"]["csv_file"] = str(
            (case_dir / "newton_history.csv").resolve()
        )
    path = case_dir / "simulation.json"
    path.write_text(
        json.dumps(config, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return path


def run_case(runner: Path, base: dict[str, Any], root: Path, name: str, policy: str) -> dict[str, Any]:
    case_dir = root / name
    case_dir.mkdir(parents=True, exist_ok=True)
    config = prepare_config(base, case_dir, policy)
    completed = subprocess.run(
        [str(runner), "--config", str(config)],
        text=True,
        capture_output=True,
    )
    (case_dir / "runner.stdout.log").write_text(completed.stdout, encoding="utf-8")
    (case_dir / "runner.stderr.log").write_text(completed.stderr, encoding="utf-8")
    if completed.returncode != 0:
        raise RuntimeError(f"{name} failed with exit code {completed.returncode}")
    iv_path = case_dir / "iv.csv"
    rows = csv_rows(iv_path)
    return {
        "name": name,
        "node_volume_policy": policy,
        "config": str(config),
        "config_sha256": sha256(config),
        "iv": str(iv_path),
        "iv_sha256": sha256(iv_path),
        "point_count": len(rows),
        "converged_count": sum(row.get("converged", "1") == "1" for row in rows),
    }


def currents(path: Path) -> dict[float, float]:
    return {
        round(float(row["bias_V"]), 9): abs(float(row["current_total_A_per_um"]))
        for row in csv_rows(path)
    }


def anchor_metrics(candidate: dict[float, float], baseline: dict[float, float], sentaurus: dict[int, float]) -> dict[str, Any]:
    rows = []
    for bias in ANCHORS:
        sent = sentaurus[bias]
        candidate_error = abs(candidate[float(bias)] / sent - 1.0)
        baseline_error = abs(baseline[float(bias)] / sent - 1.0)
        rows.append({
            "bias_V": bias,
            "sentaurus_A_per_um": sent,
            "mixed_voronoi_A_per_um": candidate[float(bias)],
            "barycentric_A_per_um": baseline[float(bias)],
            "mixed_voronoi_relative_error": candidate_error,
            "barycentric_relative_error": baseline_error,
            "error_degradation": candidate_error - baseline_error,
        })
    candidate_errors = sorted(row["mixed_voronoi_relative_error"] for row in rows)
    median = (candidate_errors[2] + candidate_errors[3]) / 2.0
    return {
        "rows": rows,
        "candidate_median_relative_error": median,
        "candidate_maximum_relative_error": max(candidate_errors),
        "maximum_error_degradation_over_barycentric": max(
            row["error_degradation"] for row in rows
        ),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--base-config", type=Path, required=True)
    parser.add_argument("--sentaurus-fields", type=Path, required=True)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = args.output_root.resolve()
    root.mkdir(parents=True, exist_ok=True)
    contract = load_json(args.contract.resolve())
    gate = contract["forward_iv_gate"]
    base = load_json(args.base_config.resolve())
    cases = [
        run_case(args.runner.resolve(), base, root, "barycentric", "barycentric"),
        run_case(args.runner.resolve(), base, root, "mixed-run-a", "mixed_voronoi"),
        run_case(args.runner.resolve(), base, root, "mixed-run-b", "mixed_voronoi"),
    ]
    by_name = {case["name"]: case for case in cases}
    baseline = currents(Path(by_name["barycentric"]["iv"]))
    mixed = currents(Path(by_name["mixed-run-a"]["iv"]))
    sentaurus = {
        bias: sentaurus_current(args.sentaurus_fields.resolve(), bias)
        for bias in ANCHORS
    }
    anchors = anchor_metrics(mixed, baseline, sentaurus)
    complete = all(
        case["point_count"] == gate["required_point_count"]
        and case["converged_count"] == gate["required_point_count"]
        for case in cases
    )
    deterministic = (
        by_name["mixed-run-a"]["iv_sha256"]
        == by_name["mixed-run-b"]["iv_sha256"]
    )
    passed = (
        complete
        and deterministic
        and anchors["candidate_median_relative_error"]
        <= gate["maximum_anchor_median_relative_error"]
        and anchors["candidate_maximum_relative_error"]
        <= gate["maximum_anchor_relative_error"]
        and anchors["maximum_error_degradation_over_barycentric"]
        <= gate["maximum_degradation_over_barycentric_at_anchor"]
    )
    report = {
        "schema": "vela.pn2d_forward_node_volume_policy_acceptance.v1",
        "status": "passed" if passed else "failed",
        "outcome": "forward_iv_gate_passed" if passed else "forward_iv_gate_failed",
        "contract": {"path": str(args.contract.resolve()), "sha256": sha256(args.contract.resolve())},
        "base_config": {"path": str(args.base_config.resolve()), "sha256": sha256(args.base_config.resolve())},
        "sentaurus_fields": str(args.sentaurus_fields.resolve()),
        "cases": cases,
        "candidate_iv_deterministic": deterministic,
        "anchors": anchors,
    }
    report_path = root / "acceptance.json"
    report_path.write_text(
        json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(report, indent=2, sort_keys=True, allow_nan=False))
    return 0 if passed else 2


if __name__ == "__main__":
    raise SystemExit(main())
