#!/usr/bin/env python3
"""Audit source-local Newton updates from one frozen TransportModels restart."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
import subprocess


def source_nodes(mesh_path: Path, contact_edges_path: Path) -> dict[str, list[int]]:
    mesh = json.loads(mesh_path.read_text(encoding="utf-8"))
    contacts = {
        entry["name"]: {int(node) for node in entry["node_ids"]}
        for entry in mesh["contacts"]
    }
    source_contact = contacts["source"]
    endpoints: set[int] = set()
    with contact_edges_path.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            if row["current_contact"] == "source":
                endpoints.add(int(row["node0"]))
                endpoints.add(int(row["node1"]))
    return {
        "contact": sorted(endpoints & source_contact),
        "interior": sorted(endpoints - source_contact),
        "all": sorted(endpoints),
    }


def configure_case(base: dict, state: Path, case_dir: Path, nodes: list[int],
                   *, qf_limit: float, damping: float, line_search: bool,
                   floor_enabled: bool, qf_reference: str) -> dict:
    config = json.loads(json.dumps(base))
    config["_comment"] = "Single-step source-local Newton/linear audit"
    solver = config["solver"]
    solver.update({
        "method": "newton",
        "max_iter": 1,
        "reltol": 0.0,
        "abstol": 0.0,
        "warm_start": True,
        "verbose": False,
        "line_search": line_search,
        "damping_factor": damping,
        "max_update": 0.0,
        "quasi_fermi_update_limit_V": qf_limit,
        "quasi_fermi_update_limit_minority_V": 0.0,
        "stall_residual_floor": 0.0,
        "poisson_line_search_stall_residual_floor": 0.0,
        "poisson_line_search_stall_carrier_residual_floor": 0.0,
        "carrier_regularization_scale": 0.0,
        "quasi_fermi_reference": qf_reference,
        "carrier_diagonal_floor": {
            "enabled": floor_enabled,
            "scale": 1.0,
            "minority_density_ratio": 1.0,
        },
        "local_update_diagnostics": {
            "enabled": True,
            "csv_file": str(case_dir / "local_updates.csv"),
            "nodes": nodes,
            "first_iterations": 2,
            "every_iterations": 10,
        },
    })
    row_check = solver.setdefault("carrier_row_convergence", {})
    row_check["mode"] = "report"
    row_check.pop("diagnostic_csv", None)
    row_check.pop("trace_csv", None)
    row_check.pop("trace_nodes", None)
    row_check["recovery"] = {"mode": "off"}
    solver["global_continuity_closure"] = {"mode": "off"}

    sweep = config["sweep"]
    sweep["initial_state_file"] = str(state)
    sweep["max_retries"] = 0
    sweep["write_state_file"] = str(case_dir / "final_state.csv")
    sweep.pop("write_state_every_point_prefix", None)
    diagnostics = sweep.setdefault("diagnostics", {})
    for section in diagnostics.values():
        if isinstance(section, dict):
            section["enabled"] = False
    diagnostics["newton_history"] = {
        "enabled": True,
        "csv_file": str(case_dir / "newton_history.csv"),
        "rejected_state_directory": str(case_dir / "rejected_states"),
    }
    diagnostics["terminal_balance"] = {
        "enabled": True,
        "contacts": ["source", "drain", "gate", "substrate"],
        "csv_file": str(case_dir / "terminal_balance.csv"),
    }
    config["output_csv"] = str(case_dir / "curve.csv")
    return config


def run_case(runner: Path, config_path: Path, backend: str) -> int:
    environment = os.environ.copy()
    environment["VELA_LINEAR_SOLVER"] = backend
    try:
        completed = subprocess.run(
            [str(runner), "--config", str(config_path)],
            env=environment,
            check=False,
            timeout=120,
        )
        return completed.returncode
    except subprocess.TimeoutExpired:
        return 124


def max_abs(rows: list[dict[str, str]], column: str) -> float:
    values = [abs(float(row[column])) for row in rows]
    return max(values, default=0.0)


def summarize_case(case: dict, case_dir: Path, interior: set[int],
                   contact: set[int], return_code: int) -> dict[str, object]:
    trace_path = case_dir / "local_updates.csv"
    if not trace_path.exists():
        return {
            "case": case["name"], "backend": case["backend"],
            "runner_return_code": return_code, "trace_written": False,
        }
    with trace_path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    electron = [row for row in rows if row["carrier"] == "electron"]
    interior_rows = [row for row in electron if int(row["node_id"]) in interior]
    contact_rows = [row for row in electron if int(row["node_id"]) in contact]
    changed_by_cap = sum(
        not math.isclose(float(row["raw_linear_step_V"]),
                         float(row["capped_step_V"]), rel_tol=0.0, abs_tol=1e-30)
        for row in interior_rows
    )
    changed_after_cap = sum(
        not math.isclose(float(row["capped_step_V"]),
                         float(row["applied_step_V"]), rel_tol=0.0, abs_tol=1e-30)
        for row in interior_rows
    )
    nonzero_raw = [
        row for row in interior_rows if float(row["raw_linear_step_V"]) != 0.0
    ]
    raw_step_over_ulp = [
        abs(float(row["raw_linear_step_V"]))
        / math.ulp(float(row["state_increment_V"]))
        for row in nonzero_raw
    ]
    lost_to_state_rounding = sum(
        float(row["raw_linear_step_V"]) != 0.0
        and float(row["capped_step_V"]) != 0.0
        and float(row["applied_step_V"]) == 0.0
        for row in interior_rows
    )
    return {
        "case": case["name"],
        "backend": case["backend"],
        "runner_return_code": return_code,
        "trace_written": True,
        "qf_limit_V": case["qf_limit"],
        "damping": case["damping"],
        "line_search": case["line_search"],
        "floor_enabled": case["floor_enabled"],
        "qf_reference": case["qf_reference"],
        "source_interior_electron_rows": len(interior_rows),
        "source_contact_electron_rows": len(contact_rows),
        "interior_nonzero_rhs_rows": sum(float(row["rhs"]) != 0.0 for row in interior_rows),
        "interior_nonzero_raw_step_rows": sum(
            float(row["raw_linear_step_V"]) != 0.0 for row in interior_rows),
        "interior_nonzero_steps_lost_to_state_rounding": lost_to_state_rounding,
        "min_abs_raw_step_over_state_ulp": min(raw_step_over_ulp, default=0.0),
        "max_abs_raw_step_over_state_ulp": max(raw_step_over_ulp, default=0.0),
        "interior_rows_changed_by_cap": changed_by_cap,
        "interior_rows_changed_after_cap": changed_after_cap,
        "max_abs_interior_rhs": max_abs(interior_rows, "rhs"),
        "max_abs_interior_raw_step_V": max_abs(interior_rows, "raw_linear_step_V"),
        "max_abs_interior_capped_step_V": max_abs(interior_rows, "capped_step_V"),
        "max_abs_interior_applied_step_V": max_abs(interior_rows, "applied_step_V"),
        "max_abs_interior_raw_linear_residual": max_abs(
            interior_rows, "raw_linear_residual"),
        "max_abs_interior_capped_linear_residual": max_abs(
            interior_rows, "capped_linear_residual"),
        "max_abs_interior_selected_trial_residual": max_abs(
            interior_rows, "selected_trial_residual"),
        "max_abs_contact_applied_step_V": max_abs(contact_rows, "applied_step_V"),
        "raw_linear_residual_l2": float(rows[0]["raw_linear_residual_l2"]),
        "raw_linear_residual_inf": float(rows[0]["raw_linear_residual_inf"]),
        "capped_linear_residual_l2": float(rows[0]["capped_linear_residual_l2"]),
        "capped_linear_residual_inf": float(rows[0]["capped_linear_residual_inf"]),
        "line_search_accepted": bool(int(rows[0]["line_search_accepted"])),
        "line_search_attempts": int(rows[0]["line_search_attempts"]),
        "selected_damping": float(rows[0]["selected_damping"]),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--base-config", type=Path, required=True)
    parser.add_argument("--state", type=Path, required=True)
    parser.add_argument("--contact-edges", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    runner = args.runner.resolve()
    base = json.loads(args.base_config.resolve().read_text(encoding="utf-8"))
    state = args.state.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    node_sets = source_nodes(Path(base["mesh_file"]), args.contact_edges.resolve())
    (output_dir / "source_nodes.json").write_text(
        json.dumps(node_sets, indent=2) + "\n", encoding="utf-8")

    cases = [
        {"name": "baseline_cap", "backend": "sparselu", "qf_limit": 0.025,
         "damping": 1.0, "line_search": False, "floor_enabled": False,
         "qf_reference": "contact_majority"},
        {"name": "no_cap", "backend": "sparselu", "qf_limit": 0.0,
         "damping": 1.0, "line_search": False, "floor_enabled": False,
         "qf_reference": "contact_majority"},
        {"name": "floor_on_no_cap", "backend": "sparselu", "qf_limit": 0.0,
         "damping": 1.0, "line_search": False, "floor_enabled": True,
         "qf_reference": "contact_majority"},
        {"name": "damped_no_cap", "backend": "sparselu", "qf_limit": 0.0,
         "damping": 0.35, "line_search": False, "floor_enabled": False,
         "qf_reference": "contact_majority"},
        {"name": "line_search_no_cap", "backend": "sparselu", "qf_limit": 0.0,
         "damping": 1.0, "line_search": True, "floor_enabled": False,
         "qf_reference": "contact_majority"},
        {"name": "bicgstab_ilut_no_cap", "backend": "bicgstab_ilut", "qf_limit": 0.0,
         "damping": 1.0, "line_search": False, "floor_enabled": False,
         "qf_reference": "contact_majority"},
        {"name": "contact_basin_no_cap", "backend": "sparselu", "qf_limit": 0.0,
         "damping": 1.0, "line_search": False, "floor_enabled": False,
         "qf_reference": "contact_basin"},
        {"name": "contact_basin_line_search", "backend": "sparselu", "qf_limit": 0.0,
         "damping": 1.0, "line_search": True, "floor_enabled": False,
         "qf_reference": "contact_basin"},
    ]
    records: list[dict[str, object]] = []
    for case in cases:
        case_dir = output_dir / case["name"]
        case_dir.mkdir(parents=True, exist_ok=True)
        config = configure_case(
            base, state, case_dir, node_sets["all"],
            qf_limit=case["qf_limit"], damping=case["damping"],
            line_search=case["line_search"],
            floor_enabled=case["floor_enabled"],
            qf_reference=case["qf_reference"],
        )
        config_path = case_dir / "config.json"
        config_path.write_text(
            json.dumps(config, indent=2, ensure_ascii=False) + "\n",
            encoding="utf-8")
        trace_path = case_dir / "local_updates.csv"
        return_code = 0 if trace_path.exists() else run_case(
            runner, config_path, case["backend"])
        record = summarize_case(
            case, case_dir, set(node_sets["interior"]),
            set(node_sets["contact"]), return_code)
        records.append(record)
        print(json.dumps(record, ensure_ascii=False))

    all_fields: list[str] = []
    for record in records:
        all_fields.extend(key for key in record if key not in all_fields)
    with (output_dir / "summary.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=all_fields)
        writer.writeheader()
        writer.writerows(records)
    execution = {
        "runner": str(runner), "base_config": str(args.base_config.resolve()),
        "state": str(state), "node_sets": node_sets, "cases": records,
    }
    (output_dir / "execution.json").write_text(
        json.dumps(execution, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
