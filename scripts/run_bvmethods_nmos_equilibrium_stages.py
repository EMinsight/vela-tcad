#!/usr/bin/env python3
"""Run the imported Sentaurus BVmethods NMOS through staged 0 V physics."""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
from copy import deepcopy
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[1]
STAGES: tuple[tuple[str, dict[str, Any]], ...] = (
    (
        "00_constant",
        {
            "recombination": [],
            "bandgap_narrowing": "none",
            "mobility": "constant",
            "impact_ionization": "none",
        },
    ),
    (
        "01_srh",
        {
            "recombination": ["srh"],
            "bandgap_narrowing": "none",
            "mobility": "constant",
            "impact_ionization": "none",
        },
    ),
    (
        "02_old_slotboom",
        {
            "recombination": ["srh"],
            "bandgap_narrowing": "old_slotboom",
            "mobility": "constant",
            "impact_ionization": "none",
        },
    ),
    (
        "03_masetti",
        {
            "recombination": ["srh"],
            "bandgap_narrowing": "old_slotboom",
            "mobility": {"model": "masetti"},
            "impact_ionization": "none",
        },
    ),
    (
        "04_masetti_high_field",
        {
            "recombination": ["srh"],
            "bandgap_narrowing": "old_slotboom",
            "mobility": {
                "model": "masetti_field",
                "high_field_driving_force": "quasi_fermi_gradient",
                "jacobian_field_derivatives": False,
            },
            "impact_ionization": "none",
        },
    ),
)


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, value: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def resolve_input(config_path: Path, value: str) -> str:
    path = Path(value)
    if not path.is_absolute():
        path = config_path.parent / path
    return str(path.resolve())


def read_last_row(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    with path.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    return rows[-1] if rows else {}


def build_stage_config(
    base_path: Path,
    base: dict[str, Any],
    stage_dir: Path,
    stage_name: str,
    physics: dict[str, Any],
    restart: Path | None,
    carrier_statistics: str,
) -> Path:
    cfg = deepcopy(base)
    for key in ("mesh_file", "node_doping_file", "materials_file"):
        if key in cfg:
            cfg[key] = resolve_input(base_path, str(cfg[key]))
    cfg["doping"] = []
    cfg["output_csv"] = str((stage_dir / "sweep.csv").resolve())

    solver = cfg.setdefault("solver", {})
    solver.update(physics)
    solver.update(
        {
            "method": "newton",
            "max_iter": 120,
            "reltol": 1.0e-8,
            "abstol": 1.0e-8,
            # Keep every restart on the same normalization used by the
            # baseline cold solve.  Otherwise a converged restart residual
            # becomes the next solve's reference and reltol asks for another
            # artificial eight orders of reduction.
            "residual_scales": {
                "psi": 1138.7290351540657,
                "phin": 1.0,
                "phip": 1.0,
            },
            "line_search": True,
            "warm_start": True,
            "quasi_fermi_update_limit_V": 0.1,
            "carrier_statistics": carrier_statistics,
        }
    )

    sweep = cfg.setdefault("sweep", {})
    sweep.update(
        {
            "mode": "bv_reverse",
            "contact": "drain",
            "current_contact": "drain",
            "start": 0.0,
            "stop": 0.0,
            "step": 0.05,
            "write_vtk": True,
            "vtk_prefix": str((stage_dir / "vtk" / "state").resolve()),
            "write_state_file": str((stage_dir / "state.csv").resolve()),
        }
    )
    sweep.pop("bias_points", None)
    if restart is None:
        sweep.pop("initial_state_file", None)
        sweep["initialization"] = {
            "mode": "poisson_block",
            "diagnostic_csv": str((stage_dir / "poisson_block.csv").resolve()),
            "write_state_file": str((stage_dir / "poisson_initial_state.csv").resolve()),
        }
    else:
        sweep.pop("initialization", None)
        sweep["initial_state_file"] = str(restart.resolve())

    cfg["_validation_stage"] = {
        "name": stage_name,
        "purpose": "Sequential 0 V physics activation for Sentaurus BVmethods NMOS",
        "restart": str(restart.resolve()) if restart else "poisson_block",
    }
    path = stage_dir / "simulation.json"
    write_json(path, cfg)
    return path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-config", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument(
        "--node-doping-file",
        type=Path,
        help="Override the base deck node_doping_file (for corrected TDR reimports)",
    )
    parser.add_argument("--runner", type=Path)
    parser.add_argument(
        "--carrier-statistics",
        choices=("boltzmann", "fermi_dirac"),
        default="boltzmann",
    )
    args = parser.parse_args()

    base_path = args.base_config.resolve()
    out_dir = args.out_dir.resolve()
    runner = args.runner
    if runner is None:
        name = "vela_example_runner.exe" if os.name == "nt" else "vela_example_runner"
        runner = REPO / "build-release" / name
    runner = runner.resolve()
    base = read_json(base_path)
    if args.node_doping_file is not None:
        base["node_doping_file"] = str(args.node_doping_file.resolve())

    summary: list[dict[str, Any]] = []
    restart: Path | None = None
    for stage_name, physics in STAGES:
        stage_dir = out_dir / stage_name
        config = build_stage_config(
            base_path, base, stage_dir, stage_name, physics, restart,
            args.carrier_statistics,
        )
        print(f":: running {stage_name}", flush=True)
        completed = subprocess.run(
            [str(runner), "--config", str(config)],
            cwd=stage_dir,
            check=False,
        )
        row = read_last_row(stage_dir / "sweep.csv")
        converged = row.get("converged") == "1"
        state = stage_dir / "state.csv"
        summary.append(
            {
                "stage": stage_name,
                "returncode": completed.returncode,
                "converged": int(converged),
                "iterations": row.get("iterations", ""),
                "current_total_A_per_um": row.get("current_total_A_per_um", ""),
                "max_electric_field_V_per_cm": row.get(
                    "max_electric_field_V_per_cm", ""
                ),
                "state_file": str(state),
            }
        )
        if completed.returncode != 0 or not converged or not state.exists():
            break
        restart = state

    out_dir.mkdir(parents=True, exist_ok=True)
    summary_path = out_dir / "summary.csv"
    with summary_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summary[0]))
        writer.writeheader()
        writer.writerows(summary)
    print(f":: summary {summary_path}", flush=True)
    return 0 if len(summary) == len(STAGES) and all(r["converged"] for r in summary) else 1


if __name__ == "__main__":
    raise SystemExit(main())
