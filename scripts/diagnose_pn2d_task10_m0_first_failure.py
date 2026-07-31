#!/usr/bin/env python3
"""Audit the first rejected Newton attempt of the Task 10 PN2D M0 run."""

from __future__ import annotations

import argparse
import copy
import csv
import json
import math
import subprocess
from pathlib import Path
from typing import Any


BRANCHES = ("avalanche_off", "iic_postprocess", "avalanche_on")


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_json(path: Path, payload: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def resolve_path(base_path: Path, raw: str) -> str:
    path = Path(raw)
    if path.is_absolute():
        return str(path.resolve())
    return str((base_path.parent / path).resolve())


def write_state_fields(state_csv: Path, fields_dir: Path) -> None:
    state = rows(state_csv)
    fields_dir.mkdir(parents=True, exist_ok=True)
    for state_field, output_name in (
        ("psi", "ElectrostaticPotential_region0.csv"),
        ("phin", "eQuasiFermiPotential_region0.csv"),
        ("phip", "hQuasiFermiPotential_region0.csv"),
    ):
        with (fields_dir / output_name).open(
            "w", newline="", encoding="utf-8"
        ) as handle:
            writer = csv.writer(handle)
            writer.writerow(["node_id", "component0"])
            for row in state:
                writer.writerow([row["node_id"], row[state_field]])


def probe_config(
    base_path: Path,
    base: dict[str, Any],
    branch: str,
    simulation_type: str,
    target_bias: float,
    fields_dir: Path,
    output_csv: Path,
) -> dict[str, Any]:
    config = copy.deepcopy(base)
    config["simulation_type"] = simulation_type
    config["output_csv"] = str(output_csv.resolve())
    config["state_fields_dir"] = str(fields_dir.resolve())
    config.pop("sweep", None)
    for key in ("mesh_file", "node_doping_file", "materials_file"):
        config[key] = resolve_path(base_path, config[key])
    for contact in config["contacts"]:
        contact["bias"] = (
            target_bias if contact["name"].lower() == "anode" else 0.0
        )
    if branch == "avalanche_off":
        config["solver"]["impact_ionization"] = {"model": "none"}
    else:
        impact = copy.deepcopy(config["solver"]["impact_ionization"])
        impact["coupling_mode"] = (
            "postprocess_only"
            if branch == "iic_postprocess"
            else "self_consistent"
        )
        config["solver"]["impact_ionization"] = impact
    return config


def run_probe(
    runner: Path, config_path: Path, config: dict[str, Any]
) -> dict[str, Any]:
    write_json(config_path, config)
    completed = subprocess.run(
        [str(runner), "--config", str(config_path)],
        text=True,
        capture_output=True,
        check=False,
    )
    lines = [line for line in completed.stdout.splitlines() if line.strip()]
    status: dict[str, Any] = {}
    if lines:
        status = json.loads(lines[-1])
    if completed.returncode != 0 or not Path(config["output_csv"]).is_file():
        raise RuntimeError(
            f"{config['simulation_type']} failed ({completed.returncode}): "
            f"{completed.stderr.strip()}"
        )
    status["returncode"] = completed.returncode
    return status


def l2(values: list[float]) -> float:
    return math.sqrt(sum(value * value for value in values))


def carrier_term_summary(path: Path) -> dict[str, Any]:
    data = rows(path)
    result: dict[str, Any] = {}
    for carrier in ("electron", "hole"):
        fields = (
            "flux",
            "recombination",
            "impact",
            "gauge",
            "boundary",
            "residual",
        )
        norms = {
            field: l2(
                [float(row[f"{carrier}_{field}"]) for row in data]
            )
            for field in fields
        }
        top_impact = max(
            data, key=lambda row: abs(float(row[f"{carrier}_impact"]))
        )
        top_residual = max(
            data, key=lambda row: abs(float(row[f"{carrier}_residual"]))
        )
        result[carrier] = {
            "l2_norms": norms,
            "top_impact_node": int(top_impact["node_id"]),
            "top_impact": float(top_impact[f"{carrier}_impact"]),
            "top_residual_node": int(top_residual["node_id"]),
            "top_residual": float(top_residual[f"{carrier}_residual"]),
        }
    return result


def first_rejected_attempt(path: Path) -> dict[str, str]:
    for row in rows(path):
        if row["status"] == "rejected":
            for field in (
                "rejected_parent_state_file",
                "rejected_initial_state_file",
                "rejected_final_state_file",
            ):
                if not row.get(field) or not Path(row[field]).is_file():
                    raise RuntimeError(
                        f"first rejected attempt is missing {field}"
                    )
            return row
    raise RuntimeError("attempt history contains no rejected attempt")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--base-config", type=Path, required=True)
    parser.add_argument("--attempts-csv", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()

    runner = args.runner.resolve()
    base_path = args.base_config.resolve()
    attempts_path = args.attempts_csv.resolve()
    output = args.out_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    base = json.loads(base_path.read_text(encoding="utf-8"))
    attempt = first_rejected_attempt(attempts_path)
    target_bias = float(attempt["actual_target_bias_V"])
    initial_state = Path(attempt["rejected_initial_state_file"]).resolve()
    fields_dir = output / "initial_state_fields"
    write_state_fields(initial_state, fields_dir)

    branch_summaries: dict[str, Any] = {}
    for branch in BRANCHES:
        branch_dir = output / branch
        branch_dir.mkdir(parents=True, exist_ok=True)
        probes: dict[str, Any] = {}
        for simulation_type, filename in (
            ("newton_residual_probe", "residual.csv"),
            ("newton_step_probe", "full_step.csv"),
            ("newton_block_step_probe", "block_steps.csv"),
            ("newton_carrier_term_probe", "carrier_terms.csv"),
            ("sg_edge_flux_probe", "sg_edges.csv"),
        ):
            csv_path = branch_dir / filename
            config = probe_config(
                base_path,
                base,
                branch,
                simulation_type,
                target_bias,
                fields_dir,
                csv_path,
            )
            if simulation_type == "newton_block_step_probe":
                config["block_modes"] = ["poisson_only", "carrier_only"]
            config_path = branch_dir / f"{simulation_type}.json"
            probes[simulation_type] = run_probe(
                runner, config_path, config
            )
        branch_summaries[branch] = {
            "probes": probes,
            "carrier_terms": carrier_term_summary(
                branch_dir / "carrier_terms.csv"
            ),
        }

    summary = {
        "schema": "vela.pn2d_task10_m0_first_failure.v1",
        "base_config": str(base_path),
        "attempts_csv": str(attempts_path),
        "attempt": {
            key: attempt[key]
            for key in (
                "attempt_id",
                "parent_accepted_bias_V",
                "parent_state_hash",
                "requested_target_bias_V",
                "actual_target_bias_V",
                "initial_state_hash",
                "status",
                "reason",
                "attempted_step_V",
                "initial_residual_norm",
                "final_residual_norm",
                "newton_iterations",
                "rejected_parent_state_file",
                "rejected_initial_state_file",
                "rejected_final_state_file",
            )
        },
        "branches": branch_summaries,
    }
    summary_path = output / "summary.json"
    write_json(summary_path, summary)
    print(summary_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
