#!/usr/bin/env python3
"""Write matched DD/DG Vela VTK diagnostics from the latest baseline states."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
BASELINE = (
    REPO_ROOT
    / "build-release/reference_tcad/transportmodels_sentaurus2022/vela_baseline"
    / "dd_dg_continuous_contact_basin_kcl_v5_2026-08-24"
)
OUTPUT_ROOT = (
    REPO_ROOT
    / "build-release/reference_tcad/transportmodels_sentaurus2022/vela_baseline"
    / "three_regime_spatial_vtk_2026-08-24"
)
RUNNER = REPO_ROOT / "build-release/vela_example_runner.exe"
MANIFEST = OUTPUT_ROOT / "three_regime_spatial_manifest.json"
STATES = (
    ("deep_off", -1.00),
    ("threshold", 0.12),
    ("on", 0.92),
)
MODE_CONFIG = {"dd": "03_dd_idvg_curve.json", "dg": "09_dg_idvg_curve.json"}


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def runner_environment() -> dict[str, str]:
    environment = os.environ.copy()
    environment["PATH"] = r"D:\msys64\ucrt64\bin" + os.pathsep + environment.get("PATH", "")
    return environment


def bias_slug(value: float) -> str:
    return f"{value:.6f}".replace("-", "m").replace(".", "p")


def run_case(mode: str, regime: str, bias: float, report_only: bool) -> dict[str, object]:
    case_dir = OUTPUT_ROOT / mode / regime
    case_dir.mkdir(parents=True, exist_ok=True)
    source_state = BASELINE / f"{mode}_idvg_curve_state_bias_{bias_slug(bias)}.csv"
    if not source_state.is_file():
        raise RuntimeError(f"Missing latest baseline state: {source_state}")
    config = json.loads((BASELINE / MODE_CONFIG[mode]).read_text(encoding="utf-8"))
    config["_comment"] = (
        "Frozen-state diagnostic export from the accepted continuous DD/DG baseline; "
        "this is not a replacement IV sweep.")
    config["output_csv"] = str((case_dir / "point.csv").resolve())
    config["log_file"] = str((case_dir / "point.log").resolve())
    sweep = config["sweep"]
    sweep.update({
        "start": bias,
        "stop": bias,
        "step": 0.01,
        "bias_points": [bias],
        "initial_state_file": str(source_state.resolve()),
        "write_vtk": True,
        "vtk_prefix": str((case_dir / "state").resolve()),
        "write_state_file": str((case_dir / "final_state.csv").resolve()),
        "write_state_every_point_prefix": str((case_dir / "state").resolve()),
    })
    sweep.setdefault("diagnostics", {})["transport"] = {"enabled": True}
    config_path = case_dir / "config.json"
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    returncode = None
    if not report_only:
        completed = subprocess.run(
            [str(RUNNER), "--config", str(config_path)],
            cwd=REPO_ROOT, text=True, capture_output=True, env=runner_environment())
        (case_dir / "console.log").write_text(
            completed.stdout + "\n--- STDERR ---\n" + completed.stderr,
            encoding="utf-8")
        returncode = completed.returncode
        if completed.returncode != 0:
            raise RuntimeError(
                f"Vela {mode}/{regime} diagnostic export failed; see {case_dir / 'console.log'}")
    vtk_files = sorted(case_dir.glob("state_*.vtk"))
    if len(vtk_files) != 1:
        raise RuntimeError(f"Expected one VTK in {case_dir}, found {len(vtk_files)}")
    vtk = vtk_files[0]
    return {
        "mode": mode,
        "regime": regime,
        "gate_bias_V": bias,
        "drain_bias_V": 1.1,
        "source_continuous_state": str(source_state.resolve()),
        "source_continuous_state_sha256": sha256(source_state),
        "config": str(config_path.resolve()),
        "execution_returncode": returncode,
        "vtk": str(vtk.resolve()),
        "vtk_sha256": sha256(vtk),
        "state_csv": str((case_dir / "final_state.csv").resolve()),
        "state_sha256": sha256(case_dir / "final_state.csv"),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report-only", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.check:
        manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
        for state in manifest["states"]:
            if sha256(Path(state["source_continuous_state"])) != state["source_continuous_state_sha256"]:
                raise RuntimeError(f"Source state changed: {state['source_continuous_state']}")
            if sha256(Path(state["vtk"])) != state["vtk_sha256"]:
                raise RuntimeError(f"VTK changed: {state['vtk']}")
        print("TransportModels Vela three-regime spatial diagnostics: PASS")
        return 0

    artifacts = [
        run_case(mode, regime, bias, args.report_only)
        for mode in ("dd", "dg")
        for regime, bias in STATES
    ]
    manifest = {
        "schema": "vela.transportmodels.vela_three_regime_spatial.v1",
        "as_of": "2026-08-24",
        "status": "complete",
        "source_baseline": str(BASELINE.resolve()),
        "states": artifacts,
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "complete", "states": len(artifacts)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
