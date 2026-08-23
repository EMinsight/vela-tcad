#!/usr/bin/env python3
"""Write five-bias Vela DG-on VTK states for the spatial oracle comparison."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
BASELINE = REPO_ROOT / "build-release/reference_tcad/transportmodels_sentaurus2022/vela_baseline"
SOURCE_CONFIG = BASELINE / "vela_fermi_bgn_ab_2026-08-21/dg_on/config.json"
SOURCE_STATE = BASELINE / "vela_fermi_bgn_ab_2026-08-21/dg_on/state_bias_m0p200000.csv"
OUTPUT_ROOT = BASELINE / "idvg_spatial_oracle_2026-08-21/vela"
RUNNER = REPO_ROOT / "build-release/vela_example_runner.exe"
MANIFEST = OUTPUT_ROOT / "spatial_oracle_manifest.json"
POINTS = (-0.20, -0.04, 0.12, 0.28, 1.00)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def runner_environment() -> dict[str, str]:
    environment = os.environ.copy()
    environment["PATH"] = r"D:\msys64\ucrt64\bin" + os.pathsep + environment.get("PATH", "")
    return environment


def bias_slug(value: float) -> str:
    return f"{value:.6f}".replace("-", "m").replace(".", "p")


def prepare_config() -> Path:
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    config = json.loads(SOURCE_CONFIG.read_text(encoding="utf-8"))
    config["_comment"] = "TransportModels DG-on five-point spatial VTK oracle"
    config["output_csv"] = str((OUTPUT_ROOT / "curve.csv").resolve())
    config["log_file"] = str((OUTPUT_ROOT / "curve.log").resolve())
    sweep = config["sweep"]
    sweep.update(
        {
            "start": POINTS[0],
            "stop": POINTS[-1],
            "step": POINTS[1] - POINTS[0],
            "bias_points": list(POINTS),
            "initial_state_file": str(SOURCE_STATE.resolve()),
            "write_vtk": True,
            "vtk_prefix": str((OUTPUT_ROOT / "state").resolve()),
            "write_state_file": str((OUTPUT_ROOT / "final_state.csv").resolve()),
            "write_state_every_point_prefix": str((OUTPUT_ROOT / "state").resolve()),
        }
    )
    sweep.setdefault("diagnostics", {})["transport"] = {"enabled": True}
    config_path = OUTPUT_ROOT / "config.json"
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    return config_path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--report-only", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.check:
        manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
        for artifact in manifest["states"]:
            assert sha256(Path(artifact["vtk"])) == artifact["vtk_sha256"]
            assert sha256(Path(artifact["state_csv"])) == artifact["state_sha256"]
        print("TransportModels Vela spatial oracle check: PASS")
        return 0

    config = prepare_config()
    execution = None
    if not args.report_only:
        completed = subprocess.run(
            [str(RUNNER), "--config", str(config)],
            cwd=REPO_ROOT,
            text=True,
            capture_output=True,
            env=runner_environment(),
        )
        (OUTPUT_ROOT / "console.log").write_text(
            completed.stdout + "\n--- STDERR ---\n" + completed.stderr,
            encoding="utf-8",
        )
        execution = {"returncode": completed.returncode}
        if completed.returncode != 0:
            raise RuntimeError(f"Vela spatial-oracle run failed: {completed.stderr}")

    vtk_files = sorted(OUTPUT_ROOT.glob("state_*.vtk"))
    state_files = [OUTPUT_ROOT / f"state_bias_{bias_slug(bias)}.csv" for bias in POINTS]
    if len(vtk_files) != len(POINTS) or len(state_files) != len(POINTS):
        raise RuntimeError(
            f"Expected {len(POINTS)} VTK/state pairs, got {len(vtk_files)}/{len(state_files)}"
        )
    states = []
    for bias, vtk, state in zip(POINTS, vtk_files, state_files):
        states.append(
            {
                "gate_bias_V": bias,
                "vtk": str(vtk.resolve()),
                "vtk_sha256": sha256(vtk),
                "state_csv": str(state.resolve()),
                "state_sha256": sha256(state),
            }
        )
    manifest = {
        "schema": "vela.transportmodels.vela_idvg_spatial_oracle.v1",
        "as_of": "2026-08-21",
        "status": "complete",
        "fixed_drain_bias_V": 1.1,
        "gate_biases_V": list(POINTS),
        "execution": execution,
        "config": str(config.resolve()),
        "states": states,
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": manifest["status"], "states": states}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
