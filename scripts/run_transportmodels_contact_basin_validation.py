#!/usr/bin/env python3
"""Run the full TransportModels deep-off point with nodal contact-basin QF references."""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import subprocess


def redirect(config: dict, output_dir: Path, trace_nodes: list[int]) -> None:
    solver = config["solver"]
    solver.update({
        "quasi_fermi_reference": "contact_basin",
        "max_iter": 200,
        "reltol": 1.0e-9,
        "abstol": 1.0e-14,
        "warm_start": True,
        "line_search": True,
        "damping_factor": 1.0,
        "quasi_fermi_update_limit_V": 0.025,
        "stall_residual_floor": 0.0,
        "poisson_line_search_stall_residual_floor": 0.0,
        "poisson_line_search_stall_carrier_residual_floor": 0.0,
        "carrier_regularization_scale": 0.0,
        "carrier_diagonal_floor": False,
        "local_update_diagnostics": {
            "enabled": True,
            "csv_file": str(output_dir / "local_updates.csv"),
            "nodes": trace_nodes,
            "first_iterations": 20,
            "every_iterations": 10,
        },
        "global_continuity_closure": {"mode": "off"},
    })
    row_check = solver.setdefault("carrier_row_convergence", {})
    row_check.update({
        "mode": "report",
        "eps_row": 1.0e-5,
        "diagnostic_csv": str(output_dir / "carrier_row_violations.csv"),
        "trace_csv": str(output_dir / "carrier_row_trace.csv"),
        "trace_nodes": trace_nodes,
        "trace_first_iterations": 20,
        "trace_every_iterations": 10,
        "recovery": {"mode": "off"},
    })

    sweep = config["sweep"]
    sweep["max_retries"] = 0
    sweep["write_state_file"] = str(output_dir / "final_state.csv")
    sweep["write_state_every_point_prefix"] = str(output_dir / "state")
    diagnostics = sweep.setdefault("diagnostics", {})
    diagnostics["terminal_balance"] = {
        "enabled": True,
        "contacts": ["source", "drain", "gate", "substrate"],
        "csv_file": str(output_dir / "terminal_balance.csv"),
    }
    diagnostics["srh_balance"] = {
        "enabled": True,
        "material": "Si",
        "drain_contact": "drain",
        "substrate_contact": "substrate",
        "kcl_contacts": ["source", "drain", "gate", "substrate"],
        "resolution_margin_ratio": 10.0,
        "csv_file": str(output_dir / "srh_balance.csv"),
    }
    diagnostics["contact_edge"] = {
        "enabled": True,
        "contacts": ["source", "drain", "substrate"],
        "csv_file": str(output_dir / "contact_edges.csv"),
    }
    diagnostics["newton_history"] = {
        "enabled": True,
        "csv_file": str(output_dir / "newton_history.csv"),
        "rejected_state_directory": str(output_dir / "rejected_states"),
    }
    config["output_csv"] = str(output_dir / "curve.csv")


def read_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--base-config", type=Path, required=True)
    parser.add_argument("--state", type=Path, required=True)
    parser.add_argument("--source-nodes", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--reuse-solve", action="store_true")
    args = parser.parse_args()

    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    config = json.loads(args.base_config.resolve().read_text(encoding="utf-8"))
    config["_comment"] = "Full DD -0.68 V validation with contact-basin QF references"
    config["sweep"]["initial_state_file"] = str(args.state.resolve())
    node_sets = json.loads(args.source_nodes.resolve().read_text(encoding="utf-8"))
    redirect(config, output_dir, node_sets["all"])
    config_path = output_dir / "config.json"
    config_path.write_text(
        json.dumps(config, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8")

    environment = os.environ.copy()
    environment["VELA_LINEAR_SOLVER"] = "sparselu"
    failure_state = output_dir / "rejected_states" / "attempt_1_bias_m0p680000_final.csv"
    if args.reuse_solve and failure_state.exists():
        solve_return_code = 1
    else:
        completed = subprocess.run(
            [str(args.runner.resolve()), "--config", str(config_path)],
            env=environment, check=False)
        solve_return_code = completed.returncode

    replay_state = output_dir / "final_state.csv"
    if not replay_state.exists():
        replay_state = failure_state
    frozen_dir = output_dir / "frozen_replay"
    frozen_terminal: list[dict[str, str]] = []
    if replay_state.exists():
        frozen_dir.mkdir(parents=True, exist_ok=True)
        frozen = json.loads(json.dumps(config))
        frozen["solver"]["method"] = "frozen_state"
        frozen["solver"].pop("local_update_diagnostics", None)
        frozen["sweep"]["initial_state_file"] = str(replay_state)
        frozen["sweep"]["frozen_state_compute_current"] = True
        frozen["sweep"]["write_state_file"] = str(frozen_dir / "replayed_state.csv")
        frozen["sweep"].pop("write_state_every_point_prefix", None)
        frozen["output_csv"] = str(frozen_dir / "curve.csv")
        diagnostics = frozen["sweep"]["diagnostics"]
        diagnostics["terminal_balance"]["csv_file"] = str(
            frozen_dir / "terminal_balance.csv")
        diagnostics["srh_balance"]["csv_file"] = str(
            frozen_dir / "srh_balance.csv")
        diagnostics["contact_edge"]["csv_file"] = str(
            frozen_dir / "contact_edges.csv")
        diagnostics["newton_history"]["enabled"] = False
        frozen_path = frozen_dir / "config.json"
        frozen_path.write_text(
            json.dumps(frozen, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8")
        subprocess.run(
            [str(args.runner.resolve()), "--config", str(frozen_path)],
            env=environment, check=False)
        frozen_terminal = read_rows(frozen_dir / "terminal_balance.csv")

    local_updates = read_rows(output_dir / "local_updates.csv")
    source_electron = [
        row for row in local_updates
        if row["carrier"] == "electron"
        and int(row["node_id"]) in set(node_sets["interior"])
    ]
    terminal = read_rows(output_dir / "terminal_balance.csv")
    history = read_rows(output_dir / "newton_history.csv")
    summary = {
        "runner_return_code": solve_return_code,
        "config": str(config_path),
        "newton_trace_rows": len(history),
        "local_update_rows": len(local_updates),
        "source_electron_nonzero_applied_rows": sum(
            float(row["applied_step_V"]) != 0.0 for row in source_electron),
        "source_electron_max_abs_applied_step_V": max(
            (abs(float(row["applied_step_V"])) for row in source_electron),
            default=0.0),
        "terminal_rows": terminal,
        "frozen_replay_state": str(replay_state) if replay_state.exists() else "",
        "frozen_terminal_rows": frozen_terminal,
        "final_state_written": (output_dir / "final_state.csv").exists(),
        "failure_diagnostics_written":
            (output_dir / "curve_newton_failure_diagnostics.json").exists(),
    }
    (output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8")
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
