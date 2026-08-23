#!/usr/bin/env python3
"""Run the 42-point DG regression after the Fermi-BGN/P2 changes."""

from __future__ import annotations

import argparse
import importlib.util
import json
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "scripts/run_transportmodels_dg_phase7_regression.py"


def load_module():
    spec = importlib.util.spec_from_file_location("transportmodels_phase7_post_p2", SOURCE)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load {SOURCE}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--report-only", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    phase7 = load_module()
    phase7.OUTPUT_ROOT = (
        phase7.BASELINE / "dg_post_p2_regression_v4_2026-08-21"
    )
    phase7.REPORT_JSON = (
        ROOT / "docs/validation/transportmodels_dg_post_p2_regression_v4_2026-08-21.json"
    )
    phase7.REPORT_MD = (
        ROOT / "docs/validation/transportmodels_dg_post_p2_regression_v4_2026-08-21.md"
    )
    phase7.PLOT_SUBTITLE = (
        "Fermi-BGN on + direct-band-edge DG + sentaurus_box + neutral interface + Lombardi"
    )
    phase7.REPORT_CONFIGURATION_EXTRA = {
        "fermi_statistics_correction": True,
        "transport_coupling": "direct_band_edge",
        "transport_coupling_weight": 1.0,
    }
    original_make_config = phase7.make_config

    def make_config(curve):
        path = original_make_config(curve)
        config = json.loads(path.read_text(encoding="utf-8"))
        config["_comment"] = (
            "TransportModels post-P2 42-point regression: Fermi-BGN on, "
            "direct-band-edge DG compatibility mode"
        )
        config["solver"]["bandgap_narrowing"] = {
            "model": "old_slotboom",
            "fermi_statistics_correction": True,
        }
        quantum = config["solver"]["electron_quantum_potential"]
        quantum["transport_coupling"] = "direct_band_edge"
        quantum["transport_coupling_weight"] = 1.0
        path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
        return path

    phase7.make_config = make_config
    if args.check:
        sys.argv = [sys.argv[0], "--check"]
        return phase7.main()
    if not args.report_only:
        phase7.OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
        phase7.prepare_idvg_restart()
        idvg_warm = phase7.run_idvg_frozen_warmup()
        idvg_terminal = phase7.OUTPUT_ROOT / "idvg/terminal_balance.csv"
        if idvg_terminal.exists():
            shutil.move(
                idvg_terminal,
                phase7.OUTPUT_ROOT / "idvg/warmup_terminal_balance.csv",
            )
        idvg_state = (
            phase7.OUTPUT_ROOT / "idvg" /
            f"state_bias_{phase7.bias_slug(phase7.CURVES[0]['points'][0])}.csv"
        )
        shutil.copyfile(idvg_warm, idvg_state)
        phase7.CURVES[0]["initial_state"] = idvg_warm

        idvd = phase7.CURVES[1]
        idvd_dir = phase7.OUTPUT_ROOT / "idvd"
        base_path = make_config(idvd)
        config = json.loads(base_path.read_text(encoding="utf-8"))
        config["_comment"] = "Post-P2 Id-Vd Frozen-Q Fermi-BGN endpoint warmup"
        config["solver"]["electron_quantum_potential"]["coupling_mode"] = "frozen"
        config["output_csv"] = str((idvd_dir / "warmup.csv").resolve())
        fermi_idvg_seed = (
            phase7.BASELINE /
            "vela_fermi_bgn_ab_2026-08-21/dg_on/state_bias_1p000000.csv"
        )
        drain_warm_points = [1.1 + 0.1 * index for index in range(10)]
        config["sweep"].update({
            "start": 1.1,
            "stop": 2.0,
            "step": 0.1,
            "bias_points": drain_warm_points,
            "initial_state_file": str(fermi_idvg_seed.resolve()),
            "write_state_file": str((idvd_dir / "warmup_final_state.csv").resolve()),
            "write_state_every_point_prefix": str((idvd_dir / "warmup_state").resolve()),
        })
        config["sweep"]["diagnostics"]["terminal_balance"]["csv_file"] = str(
            (idvd_dir / "warmup_terminal_balance.csv").resolve()
        )
        warm_path = idvd_dir / "warmup_config.json"
        warm_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
        process = subprocess.run(
            [str(phase7.RUNNER), "--config", str(warm_path)],
            cwd=ROOT, text=True, capture_output=True,
            env=phase7.runner_environment(),
        )
        (idvd_dir / "warmup.console.log").write_text(
            process.stdout + "\n--- STDERR ---\n" + process.stderr,
            encoding="utf-8",
        )
        idvd_warm = idvd_dir / "warmup_final_state.csv"
        if process.returncode != 0 or not idvd_warm.exists():
            raise RuntimeError("Id-Vd Frozen-Q warmup failed: " + process.stderr)
        idvd_state = idvd_dir / f"state_bias_{phase7.bias_slug(2.0)}.csv"
        shutil.copyfile(idvd_warm, idvd_state)
        idvd["initial_state"] = idvd_warm

        with ThreadPoolExecutor(max_workers=2) as pool:
            futures = [pool.submit(phase7.execute_curve, curve) for curve in phase7.CURVES]
            for future in as_completed(futures):
                future.result()

    sys.argv = [sys.argv[0], "--report-only"]
    return phase7.main()


if __name__ == "__main__":
    raise SystemExit(main())
