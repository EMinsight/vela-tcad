#!/usr/bin/env python3
"""Generate the A/B/C/D Sentaurus full-physics BV ablation decks."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = (
    REPO / "build-release/reference_tcad/bvmethods_sentaurus2018/run01/full_raw"
)
DEFAULT_OUTPUT = (
    REPO / "build-release/reference_tcad/bvmethods_sentaurus2018/run01"
    / "sentaurus_full_physics_ablations_20260812/bundle"
)
DEFAULT_VELA_BASE = (
    REPO / "build-release/reference_tcad/bvmethods_sentaurus2018/run01"
    / "vela_validation/boundary_voltage_to_current_20260806/simulation.json"
)

VARIANTS = {
    "a_constant_no_enormal": (False, False),
    "b_doping_no_enormal": (True, False),
    "c_constant_enormal": (False, True),
    "d_doping_enormal": (True, True),
}


def build_deck(template: str, name: str, doping_srh: bool, enormal: bool) -> str:
    text = template
    text = text.replace('Plot      = "n6_des.tdr"', f'Plot      = "{name}_des.tdr"')
    text = text.replace('Current   = "n6_des.plt"', f'Current   = "{name}_des.plt"')
    text = text.replace('Output    = "n6_des.log"', f'Output    = "{name}_des.log"')
    text = text.replace(
        'Parameter = "pp6_des.par"',
        f'Parameter = "{"pp6_des.par" if doping_srh else "constant_srh.par"}"',
    )
    if not enormal:
        text = re.sub(r"\s*Enormal\s*", "\n", text, count=1)
    if not doping_srh:
        text = re.sub(r"SRH\s*\(\s*DopingDep\s*\)", "SRH", text, count=1)
    text, density_count = re.subn(
        r"eDensity\s+hDensity",
        "eDensity hDensity eQuasiFermi hQuasiFermi\n  eMobility hMobility",
        text,
        count=1,
    )
    if density_count != 1:
        raise ValueError("expected one eDensity/hDensity Plot entry")
    text = text.replace("Current=1.443e-3", "Current=1e-4")
    text = text.replace("current=1.443e-3", "current=1e-4")
    text = text.replace("AbsVal=1.443e-3", "AbsVal=1.2e-4")
    if "current=1e-4" not in text.lower():
        raise ValueError("failed to set the 1e-4 A/um current target")
    return text


def build_vela_config(
    base: dict, output_root: Path, name: str, doping_srh: bool, enormal: bool
) -> dict:
    config = json.loads(json.dumps(base))
    case_dir = output_root / name
    config["output_csv"] = str((case_dir / "sweep.csv").resolve())
    solver = config["solver"]
    mobility = solver["mobility"]
    mobility["model"] = (
        "masetti_field_lombardi" if enormal else "masetti_field"
    )
    mobility["doping_concentration_basis"] = "total_impurity"
    if enormal:
        mobility["surface"] = {
            "surface_region": "R.Substrate",
            "surface_interface": ["R.Substrate", "R.Gateox"],
        }
    else:
        mobility.pop("surface", None)
    if doping_srh:
        solver["srh_doping_dependence"] = {
            "enabled": True,
            "concentration_basis": "total_impurity",
            "electron": {
                "tau_min_s": 0.0,
                "tau_max_s": 1.0e-7,
                "reference_doping_m3": 1.0e16,
                "gamma": 1.0,
            },
            "hole": {
                "tau_min_s": 0.0,
                "tau_max_s": 1.0e-7,
                "reference_doping_m3": 1.0e16,
                "gamma": 1.0,
            },
        }
    else:
        solver.pop("srh_doping_dependence", None)
    carrier_rows = solver.get("carrier_row_convergence", {})
    if carrier_rows:
        carrier_rows["diagnostic_csv"] = str(
            (case_dir / "carrier_row_convergence.csv").resolve()
        )
        carrier_rows["trace_csv"] = str(
            (case_dir / "carrier_row_trace.csv").resolve()
        )
    sweep = config["sweep"]
    sweep["vtk_prefix"] = str((case_dir / "vtk/state").resolve())
    sweep["write_state_file"] = str((case_dir / "last_state.csv").resolve())
    sweep["write_state_every_point_prefix"] = str(
        (case_dir / "states/accepted_state").resolve()
    )
    boundary = sweep.get("boundary_control", {})
    if boundary:
        boundary["evaluation_csv"] = str(
            (case_dir / "boundary_control_evaluations.csv").resolve()
        )
        boundary["checkpoint_directory"] = str(
            (case_dir / "boundary_control_checkpoints").resolve()
        )
        boundary["resume"] = False
    diagnostics = sweep.get("diagnostics", {})
    qf_bounds = diagnostics.get("qf_bounds", {})
    if qf_bounds:
        qf_bounds["csv_file"] = str((case_dir / "qf_bounds.csv").resolve())
    newton_history = diagnostics.get("newton_history", {})
    if newton_history:
        newton_history["csv_file"] = str(
            (case_dir / "newton_history.csv").resolve()
        )
        newton_history["attempts_csv_file"] = str(
            (case_dir / "newton_attempts.csv").resolve()
        )
        newton_history["iterations_csv_file"] = str(
            (case_dir / "newton_iterations.csv").resolve()
        )
    config["_full_physics_ablation"] = {
        "variant": name,
        "srh": "doping_dependent" if doping_srh else "constant_1e-7_s",
        "enormal": enormal,
        "target_current_A_per_um": 1.0e-4,
    }
    return config


def build_vela_external_resistor_config(
    full_config: dict,
    case_dir: Path,
    initial_state_file: Path,
    outer_voltage_V: float = 1006.0,
    resistance_ohm_um: float = 1.0e7,
) -> dict:
    """Convert the full D voltage-to-current case to an independent load line."""
    config = json.loads(json.dumps(full_config))
    config["output_csv"] = str((case_dir / "sweep.csv").resolve())
    solver = config["solver"]
    carrier_rows = solver.get("carrier_row_convergence", {})
    if carrier_rows:
        carrier_rows["diagnostic_csv"] = str(
            (case_dir / "carrier_row_convergence.csv").resolve()
        )
        carrier_rows["trace_csv"] = str(
            (case_dir / "carrier_row_trace.csv").resolve()
        )

    sweep = config["sweep"]
    sweep["start"] = outer_voltage_V
    sweep["stop"] = outer_voltage_V
    sweep["bias_points"] = [outer_voltage_V]
    sweep["initial_state_file"] = str(initial_state_file.resolve())
    sweep["write_state_file"] = str((case_dir / "last_state.csv").resolve())
    sweep["write_state_every_point_prefix"] = str(
        (case_dir / "states/accepted_state").resolve()
    )
    sweep["vtk_prefix"] = str((case_dir / "vtk/state").resolve())
    sweep.pop("voltage_to_current", None)
    sweep.pop("continuation", None)
    sweep["external_circuit"] = {
        "mode": "series_resistor",
        "resistance_ohm_um": resistance_ohm_um,
        "current_direction": 1.0,
        # The boundary solver's first bracket probe is one max step below this
        # value.  Seed that probe at the converged D voltage-to-current state.
        "initial_inner_voltage_V": 6.4069,
        "max_inner_voltage_step_V": 0.005,
        "residual_tolerance_V": 1.0e-1,
        "voltage_tolerance_V": 1.0e-8,
        "max_bracket_steps": 120,
        "max_iterations": 40,
    }
    boundary = sweep.get("boundary_control", {})
    if boundary:
        boundary["evaluation_csv"] = str(
            (case_dir / "boundary_control_evaluations.csv").resolve()
        )
        boundary["checkpoint_directory"] = str(
            (case_dir / "boundary_control_checkpoints").resolve()
        )
        boundary["resume"] = False
        boundary["predictor_max_step_factor"] = 1.0
    diagnostics = sweep.get("diagnostics", {})
    qf_bounds = diagnostics.get("qf_bounds", {})
    if qf_bounds:
        qf_bounds["csv_file"] = str((case_dir / "qf_bounds.csv").resolve())
    history = diagnostics.get("newton_history", {})
    if history:
        history["csv_file"] = str((case_dir / "newton_history.csv").resolve())
        history["attempts_csv_file"] = str(
            (case_dir / "newton_attempts.csv").resolve()
        )
        history["iterations_csv_file"] = str(
            (case_dir / "newton_iterations.csv").resolve()
        )
    config["_full_physics_cross_check"] = {
        "method": "external_resistor",
        "outer_voltage_V": outer_voltage_V,
        "resistance_ohm_um": resistance_ohm_um,
        "expected_current_A_per_um": 1.0e-4,
    }
    return config


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--vela-base-config", type=Path, default=DEFAULT_VELA_BASE)
    parser.add_argument("--vela-output-dir", type=Path)
    args = parser.parse_args()

    source = args.source_dir.resolve()
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    vela_output = (args.vela_output_dir or output.parent / "vela").resolve()
    vela_output.mkdir(parents=True, exist_ok=True)
    vela_base = json.loads(
        args.vela_base_config.resolve().read_text(encoding="utf-8-sig")
    )
    template = (source / "pp6_des.cmd").read_text(encoding="utf-8")
    constant = (
        REPO
        / "reference_tcad/bvmethods_sentaurus2018/source/full_physics_constant_srh.par"
    ).read_text(encoding="utf-8")
    (output / "constant_srh.par").write_text(constant, encoding="utf-8")
    for name in ("n1_msh.tdr", "pp6_des.par"):
        (output / name).write_bytes((source / name).read_bytes())

    manifest = {"schema": "vela.bvmethods_full_physics_ablations.v1", "variants": {}}
    for name, (doping_srh, enormal) in VARIANTS.items():
        deck = build_deck(template, name, doping_srh, enormal)
        path = output / f"{name}.cmd"
        path.write_text(deck, encoding="utf-8")
        manifest["variants"][name] = {
            "deck": path.name,
            "srh": "doping_dependent" if doping_srh else "constant_1e-7_s",
            "enormal": enormal,
            "target_current_A_per_um": 1.0e-4,
        }
        vela_config = build_vela_config(
            vela_base, vela_output, name, doping_srh, enormal
        )
        vela_path = vela_output / name / "simulation.json"
        vela_path.parent.mkdir(parents=True, exist_ok=True)
        vela_path.write_text(
            json.dumps(vela_config, indent=2) + "\n", encoding="utf-8"
        )
        manifest["variants"][name]["vela_config"] = str(vela_path)
        if name == "d_doping_enormal":
            external_dir = vela_output / name / "external_resistor_final"
            external_config = build_vela_external_resistor_config(
                vela_config,
                external_dir,
                vela_output / name / "last_state.csv",
            )
            external_path = external_dir / "simulation.json"
            external_path.parent.mkdir(parents=True, exist_ok=True)
            external_path.write_text(
                json.dumps(external_config, indent=2) + "\n", encoding="utf-8"
            )
            manifest["variants"][name]["vela_external_resistor_config"] = str(
                external_path
            )
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
