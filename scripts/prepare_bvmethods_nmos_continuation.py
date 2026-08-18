#!/usr/bin/env python3
"""Prepare a BVmethods NMOS pseudo-arclength run without changing physics.

The input is an already-qualified voltage or voltage-to-current Vela deck.
This tool retains its mesh, materials, contacts, and solver block verbatim and
replaces only the sweep controller with a restartable pseudo-arclength branch.
"""

from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path
from typing import Any


def write_text_lf(path: Path, text: str) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write(text)


def prepare_config(
    base: dict[str, Any],
    initial_state: Path,
    initial_secant_state: Path,
    output_dir: Path,
    *,
    start_V: float = 6.056459,
    initial_secant_bias_V: float = 6.0,
    stop_V: float = 6.5,
) -> dict[str, Any]:
    if base.get("simulation_type") != "dc_sweep":
        raise ValueError("base config must use simulation_type=dc_sweep")
    solver = base.get("solver", {})
    if solver.get("method") not in {"newton", "gummel_newton"}:
        raise ValueError("pseudo-arclength requires a Newton-capable base config")
    impact = solver.get("impact_ionization", {})
    if impact.get("coupling_mode") != "self_consistent":
        raise ValueError("BVmethods continuation requires self-consistent avalanche")
    if not initial_state.is_file():
        raise FileNotFoundError(initial_state)
    if not initial_secant_state.is_file():
        raise FileNotFoundError(initial_secant_state)
    if initial_secant_bias_V >= start_V:
        raise ValueError("initial secant bias must be below the start voltage")
    if stop_V <= start_V:
        raise ValueError("stop voltage must exceed start voltage")

    output_dir = output_dir.resolve()
    config = copy.deepcopy(base)
    config["output_csv"] = str(output_dir / "sweep.csv")
    config["sweep"] = {
        "mode": "bv_reverse",
        "start": start_V,
        "stop": stop_V,
        "step": 0.01,
        "initial_step": 0.01,
        "min_step": 1.0e-8,
        "max_step": 0.1,
        "growth_factor": 1.25,
        "shrink_factor": 0.5,
        "max_retries": 20,
        "stop_on_failure": True,
        "contact": "drain",
        "current_contact": "drain",
        "initial_state_file": str(initial_state.resolve()),
        "write_state_file": str(output_dir / "last_state.csv"),
        "write_state_every_point_prefix": str(output_dir / "states" / "state"),
        "write_vtk": False,
        "breakdown": {
            "max_electric_field_V_per_m": 1.0e12,
            "current_jump_ratio": 1.0e12,
            "non_convergence": False,
        },
        "diagnostics": {
            "qf_bounds": {
                "enabled": True,
                "mode": "reject_and_recover",
                "margin_V": 1.0,
                "min_carrier_density_m3": 1.0e6,
                "csv_file": str(output_dir / "qf_bounds.csv"),
            },
            "newton_history": {
                "enabled": True,
                "csv_file": str(output_dir / "newton_history.csv"),
                "attempts_csv_file": str(output_dir / "newton_attempts.csv"),
                "iterations_csv_file": str(output_dir / "newton_iterations.csv"),
            },
        },
        "continuation": {
            "arclength": {
                "enabled": True,
                "predictor": "tangent",
                "initial_step": 0.01,
                "min_step": 1.0e-6,
                "max_step": 0.1,
                "growth_factor": 1.2,
                "shrink_factor": 0.5,
                "max_corrector_iterations": 20,
                "corrector_tolerance": 1.0e-8,
                "max_step_retries": 12,
                "parameter_scale": 1.0,
                # The packed NMOS state contains thousands of scaled carrier
                # unknowns.  The mesh-default 1/N weight makes its tangent norm
                # dominate lambda and advances voltage by only ~1e-9 V/step.
                # This explicit numerical weight balances state and voltage;
                # it does not alter F(x, V) or any physical coefficient.
                "state_weight": 1.0e-15,
                "damping_factor": 0.5,
                "max_line_search_steps": 12,
                "max_parameter_update": 0.02,
                "bias_finite_difference_step_V": 1.0e-4,
                "initial_secant_state_file":
                    str(initial_secant_state.resolve()),
                "initial_secant_bias_V": initial_secant_bias_V,
            }
        },
    }
    config["_validation_case"] = {
        "purpose": "BVmethods NMOS non-transient continuation closure",
        "source_config": "physics copied verbatim from supplied base config",
        "sentaurus_reference_BV_V": 6.383727168968036,
        "current_threshold_A_per_um": 1.0e-4,
        "acceptance_relative_error": 0.02,
        "physics_parameter_scaling": "none",
    }
    return config


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-config", type=Path, required=True)
    parser.add_argument("--initial-state", type=Path, required=True)
    parser.add_argument("--initial-secant-state", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--start-V", type=float, default=6.056459)
    parser.add_argument("--initial-secant-bias-V", type=float, default=6.0)
    parser.add_argument("--stop-V", type=float, default=6.5)
    args = parser.parse_args()

    base = json.loads(args.base_config.read_text(encoding="utf-8"))
    config = prepare_config(
        base, args.initial_state, args.initial_secant_state, args.output_dir,
        start_V=args.start_V,
        initial_secant_bias_V=args.initial_secant_bias_V,
        stop_V=args.stop_V)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    output = args.output_dir / "simulation.json"
    write_text_lf(output, json.dumps(config, indent=2) + "\n")
    print(output.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
