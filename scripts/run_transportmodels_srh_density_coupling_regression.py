#!/usr/bin/env python3
"""Run TransportModels DG regressions for Sentaurus-default SRH density coupling."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
import os
import subprocess
from pathlib import Path
from typing import Any

import transportmodels_fixed_contract as fixed


REPO = Path(__file__).resolve().parents[1]
REF = REPO / "build-release/reference_tcad/transportmodels_sentaurus2022"
OUTPUT = REF / "vela_baseline/dg_srh_density_coupling_sentaurus_default_2026-08-23"
RUNNER = REPO / "build-release/vela_example_runner.exe"
STRICT_SCRIPT = REPO / "scripts/run_transportmodels_dg_deep_off_strict.py"
FULL_SCRIPT = REPO / "scripts/run_transportmodels_dg_quantum_contract_regression.py"
BASE_RUN = REF / "vela_baseline/dg_quantum_contract_regression_2026-08-23/runs/dg"
BASE_CONFIG = BASE_RUN / "03_dg_idvg_curve.json"
SENT_IDVG = REF / "run02/normalized/dg_idvg.csv"
KEY_BIASES = (-0.20, -0.04, 0.12, 0.28, 1.00)


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def environment() -> dict[str, str]:
    result = os.environ.copy()
    result["PATH"] = os.pathsep.join(
        [r"D:\msys64\ucrt64\bin", r"D:\msys64\usr\bin", result.get("PATH", "")]
    )
    return result


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def reference_idvg() -> dict[float, float]:
    return {
        round(float(row["bias_V"]), 12): abs(float(row["current_total"]))
        for row in read_csv(SENT_IDVG)
    }


def run_three_points(runner: Path) -> dict[str, Any]:
    strict = load_module(STRICT_SCRIPT, "transportmodels_strict_srh_density")
    strict.OUTPUT = OUTPUT / "three_deep_off_points"
    original_config = strict.strict_config

    def patched_config(bias: float, run_dir: Path, variant: str):
        config = fixed.apply_contract(original_config(bias, run_dir, variant), "dg")
        config["_comment"] += "; SRH uses Sentaurus-default classical density"
        return config

    strict.strict_config = patched_config
    references = strict.sentaurus_reference()
    rows = [
        strict.execute_point(runner, bias, references[round(bias, 12)], "scaled_filter")
        for bias in strict.BIAS_STATES
    ]
    strict.write_summary(rows, runner, "scaled_filter")
    return {
        "phase": "three_points",
        "all_hard_accepted": all(row["hard_acceptance"] for row in rows),
        "points": rows,
    }


def run_five_points(runner: Path) -> dict[str, Any]:
    run_dir = OUTPUT / "five_key_vg_points"
    run_dir.mkdir(parents=True, exist_ok=True)
    config = fixed.apply_contract(
        json.loads(BASE_CONFIG.read_text(encoding="utf-8")), "dg"
    )
    config["_comment"] = (
        "Five key DG Id-Vg points with explicit srh_density_coupling="
        "sentaurus_default and analytic generalized Fermi-SRH Jacobian"
    )
    config["output_csv"] = str((run_dir / "curve.csv").resolve())
    config["log_file"] = str((run_dir / "curve.log").resolve())
    config["solver"]["srh_density_coupling"] = "sentaurus_default"
    config["solver"]["verbose"] = False
    for contact in config["contacts"]:
        if contact["name"] == "gate":
            contact["bias"] = KEY_BIASES[0]
    sweep = config["sweep"]
    sweep.update(
        {
            "start": KEY_BIASES[0],
            "stop": KEY_BIASES[-1],
            "step": KEY_BIASES[1] - KEY_BIASES[0],
            "bias_points": list(KEY_BIASES),
            "initial_state_file": str(
                (BASE_RUN / "dg_idvg_curve_state_bias_m0p200000.csv").resolve()
            ),
            "write_state_file": str((run_dir / "final_state.csv").resolve()),
            "write_state_every_point_prefix": str((run_dir / "state").resolve()),
        }
    )
    diagnostics = sweep.setdefault("diagnostics", {})
    diagnostics["srh_balance"] = {
        "enabled": True,
        "material": "Si",
        "drain_contact": "drain",
        "substrate_contact": "substrate",
        "kcl_contacts": ["source", "drain", "gate", "substrate"],
        "resolution_margin_ratio": 10.0,
        "csv_file": str((run_dir / "srh_balance.csv").resolve()),
    }
    config_path = run_dir / "config.json"
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    existing_rows = read_csv(run_dir / "curve.csv") if (run_dir / "curve.csv").is_file() else []
    existing_biases = {
        round(float(row["bias_V"]), 12)
        for row in existing_rows
        if row.get("converged") == "1"
    }
    expected_biases = {round(value, 12) for value in KEY_BIASES}
    if expected_biases.issubset(existing_biases):
        returncode = 0
    else:
        completed = subprocess.run(
            [str(runner), "--config", str(config_path)],
            cwd=REPO,
            env=environment(),
            text=True,
            capture_output=True,
            check=False,
        )
        returncode = completed.returncode
        (run_dir / "console.log").write_text(
            completed.stdout + completed.stderr, encoding="utf-8"
        )
    reference = reference_idvg()
    aligned: list[dict[str, Any]] = []
    if (run_dir / "curve.csv").is_file():
        for row in read_csv(run_dir / "curve.csv"):
            if row.get("converged") != "1":
                continue
            bias = round(float(row["bias_V"]), 12)
            if bias not in {round(value, 12) for value in KEY_BIASES}:
                continue
            vela = abs(float(row["current_total_A_per_um"]))
            if bias in reference:
                sentaurus = reference[bias]
                reference_method = "exact_21_point_lattice"
            else:
                lower = max(value for value in reference if value < bias)
                upper = min(value for value in reference if value > bias)
                fraction = (bias - lower) / (upper - lower)
                sentaurus = reference[lower] + fraction * (
                    reference[upper] - reference[lower]
                )
                reference_method = f"linear_interpolation_{lower:g}_{upper:g}_V"
            aligned.append(
                {
                    "bias_V": bias,
                    "vela_A_per_um": vela,
                    "sentaurus_A_per_um": sentaurus,
                    "sentaurus_reference_method": reference_method,
                    "absolute_relative_error": abs(vela - sentaurus) / sentaurus,
                    "absolute_log_error_dex": abs(
                        math.log10(max(vela, 1.0e-300)) - math.log10(sentaurus)
                    ),
                }
            )
    aligned.sort(key=lambda row: row["bias_V"])
    summary = {
        "schema": "vela.transportmodels.srh_density_coupling.five_vg.v1",
        "runner_returncode": returncode,
        "completed_points": len(aligned),
        "expected_points": len(KEY_BIASES),
        "points": aligned,
        "config": str(config_path.resolve()),
    }
    (run_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    return {"phase": "five_points", **summary}


def run_full_curves(runner: Path) -> dict[str, Any]:
    full = load_module(FULL_SCRIPT, "transportmodels_full_srh_density")
    full.OUTPUT = OUTPUT / "full_21_point_curves"
    full.RUN_DIR = full.OUTPUT / "runs/dg"
    full.DEFAULT_RUNNER = runner
    full.REPORT_JSON = (
        REPO / "docs/validation/transportmodels_dg_srh_density_coupling_2026-08-23.json"
    )
    full.REPORT_MD = (
        REPO / "docs/validation/transportmodels_dg_srh_density_coupling_2026-08-23.md"
    )
    full.RUN_DIR.mkdir(parents=True, exist_ok=True)
    curve_specs = (
        {
            "name": "idvg",
            "source_config": BASE_CONFIG,
            "biases": [round(-1.0 + 0.16 * index, 12) for index in range(21)],
            "initial_state": OUTPUT
                / "three_deep_off_points/scaled_filter/m1p00/final_state.csv",
        },
        {
            "name": "idvd",
            "source_config": BASE_RUN / "05_dg_idvd_curve.json",
            "biases": [round(0.1 * index, 12) for index in range(21)],
            "initial_state": BASE_RUN / "dg_idvd_equilibrium_final_state.csv",
        },
    )
    stages: list[dict[str, Any]] = []
    for spec in curve_specs:
        name = spec["name"]
        config = fixed.apply_contract(
            json.loads(Path(spec["source_config"]).read_text(encoding="utf-8")), "dg"
        )
        config["_comment"] = (
            "Complete 21-point DG curve with analytic Sentaurus-default "
            "classical-density SRH coupling"
        )
        config["output_csv"] = str(
            (full.RUN_DIR / f"dg_{name}_curve_comparison_candidate.csv").resolve()
        )
        config["log_file"] = str((full.RUN_DIR / f"dg_{name}_curve.log").resolve())
        config["solver"]["srh_density_coupling"] = "sentaurus_default"
        config["solver"]["warm_start"] = True
        config["solver"]["verbose"] = False
        sweep = config["sweep"]
        sweep.update(
            {
                "start": spec["biases"][0],
                "stop": spec["biases"][-1],
                "step": spec["biases"][1] - spec["biases"][0],
                "bias_points": spec["biases"],
                "initial_state_file": str(Path(spec["initial_state"]).resolve()),
                "write_state_file": str(
                    (full.RUN_DIR / f"dg_{name}_curve_final_state.csv").resolve()
                ),
                "write_state_every_point_prefix": str(
                    (full.RUN_DIR / f"dg_{name}_curve_state").resolve()
                ),
            }
        )
        diagnostics = sweep.setdefault("diagnostics", {})
        diagnostics["srh_balance"] = {
            "enabled": True,
            "material": "Si",
            "drain_contact": "drain",
            "substrate_contact": "substrate",
            "kcl_contacts": ["source", "drain", "gate", "substrate"],
            "resolution_margin_ratio": 10.0,
            "csv_file": str(
                (full.RUN_DIR / f"dg_{name}_curve_srh_balance.csv").resolve()
            ),
        }
        config_path = full.RUN_DIR / f"dg_{name}_curve.json"
        config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
        candidate_path = full.RUN_DIR / f"dg_{name}_curve_comparison_candidate.csv"
        existing = read_csv(candidate_path) if candidate_path.is_file() else []
        expected = {round(float(value), 12) for value in spec["biases"]}
        converged = {
            round(float(row["bias_V"]), 12)
            for row in existing
            if row.get("converged") == "1"
        }
        if expected.issubset(converged):
            returncode = 0
        else:
            completed = subprocess.run(
                [str(runner), "--config", str(config_path)],
                cwd=REPO,
                env=environment(),
                text=True,
                capture_output=True,
                check=False,
            )
            returncode = completed.returncode
            (full.RUN_DIR / f"dg_{name}_curve.console.log").write_text(
                completed.stdout + completed.stderr, encoding="utf-8"
            )
        stages.append(
            {
                "name": f"dg_{name}_curve",
                "config": str(config_path.resolve()),
                "returncode": returncode,
                "status": "complete" if returncode == 0 else "fail",
                "expected_points": 21,
            }
        )
        if returncode != 0:
            break

    if len(stages) == 2 and all(stage["returncode"] == 0 for stage in stages):
        prefix_biases = (-0.68, -0.52, -0.36, -0.20, -0.04, 0.12, 0.28)
        prefix_config = fixed.apply_contract(
            json.loads(BASE_CONFIG.read_text(encoding="utf-8")), "dg"
        )
        prefix_config["_comment"] = (
            "Strict-seed Id-Vg prefix proving convergence from the hard-gated "
            "-0.68 V state into the transition branch"
        )
        prefix_config["output_csv"] = str(
            (full.RUN_DIR / "dg_idvg_strict_seed_prefix.csv").resolve()
        )
        prefix_config["log_file"] = str(
            (full.RUN_DIR / "dg_idvg_strict_seed_prefix.log").resolve()
        )
        prefix_config["solver"]["srh_density_coupling"] = "sentaurus_default"
        prefix_config["solver"]["warm_start"] = True
        prefix_config["solver"]["verbose"] = False
        prefix_sweep = prefix_config["sweep"]
        prefix_sweep.update(
            {
                "start": prefix_biases[0],
                "stop": prefix_biases[-1],
                "step": 0.16,
                "bias_points": list(prefix_biases),
                "initial_state_file": str(
                    (OUTPUT / "three_deep_off_points/scaled_filter/m0p68/final_state.csv").resolve()
                ),
                "write_state_file": str(
                    (full.RUN_DIR / "dg_idvg_strict_seed_prefix_final_state.csv").resolve()
                ),
                "write_state_every_point_prefix": str(
                    (full.RUN_DIR / "dg_idvg_strict_seed_prefix_state").resolve()
                ),
            }
        )
        prefix_sweep.setdefault("diagnostics", {})["srh_balance"] = {
            "enabled": True,
            "material": "Si",
            "drain_contact": "drain",
            "substrate_contact": "substrate",
            "kcl_contacts": ["source", "drain", "gate", "substrate"],
            "resolution_margin_ratio": 10.0,
            "csv_file": str(
                (full.RUN_DIR / "dg_idvg_strict_seed_prefix_srh_balance.csv").resolve()
            ),
        }
        prefix_config_path = full.RUN_DIR / "dg_idvg_strict_seed_prefix.json"
        prefix_config_path.write_text(
            json.dumps(prefix_config, indent=2) + "\n", encoding="utf-8"
        )
        prefix_rows = read_csv(full.RUN_DIR / "dg_idvg_strict_seed_prefix.csv") \
            if (full.RUN_DIR / "dg_idvg_strict_seed_prefix.csv").is_file() else []
        if len([row for row in prefix_rows if row.get("converged") == "1"]) != 7:
            completed = subprocess.run(
                [str(runner), "--config", str(prefix_config_path)],
                cwd=REPO,
                env=environment(),
                text=True,
                capture_output=True,
                check=False,
            )
            (full.RUN_DIR / "dg_idvg_strict_seed_prefix.console.log").write_text(
                completed.stdout + completed.stderr, encoding="utf-8"
            )
            if completed.returncode != 0:
                raise RuntimeError("Strict-seed Id-Vg prefix failed")
            prefix_rows = read_csv(full.RUN_DIR / "dg_idvg_strict_seed_prefix.csv")

        candidate = full.RUN_DIR / "dg_idvg_curve_comparison_candidate.csv"
        ordinary_rows = read_csv(candidate)
        ordinary_archive = full.RUN_DIR / "dg_idvg_curve_ordinary_tolerance.csv"
        if not ordinary_archive.is_file():
            ordinary_archive.write_text(candidate.read_text(encoding="utf-8"), encoding="utf-8")
        strict_roots = [
            OUTPUT / "three_deep_off_points/scaled_filter/m1p00",
            OUTPUT / "three_deep_off_points/scaled_filter/m0p84",
            OUTPUT / "three_deep_off_points/scaled_filter/m0p68",
        ]
        strict_rows = [read_csv(root / "curve.csv")[-1] for root in strict_roots]
        merged_rows = strict_rows + [
            row for row in prefix_rows if float(row["bias_V"]) > -0.68 + 1.0e-12
        ] + [
            row for row in ordinary_rows if float(row["bias_V"]) > 0.28 + 1.0e-12
        ]
        fieldnames = list(merged_rows[0])
        with candidate.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(merged_rows)

        srh_path = full.RUN_DIR / "dg_idvg_curve_srh_balance.csv"
        ordinary_srh = read_csv(srh_path)
        strict_srh = [read_csv(root / "srh_balance.csv")[-1] for root in strict_roots]
        prefix_srh = read_csv(full.RUN_DIR / "dg_idvg_strict_seed_prefix_srh_balance.csv")
        merged_srh = strict_srh + [
            row for row in prefix_srh if float(row["bias_V"]) > -0.68 + 1.0e-12
        ] + [
            row for row in ordinary_srh if float(row["bias_V"]) > 0.28 + 1.0e-12
        ]
        with srh_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(merged_srh[0]))
            writer.writeheader()
            writer.writerows(merged_srh)
        stages.append(
            {
                "name": "dg_idvg_strict_seed_prefix",
                "config": str(prefix_config_path.resolve()),
                "returncode": 0,
                "status": "complete",
                "expected_points": 7,
            }
        )
    manifest = {
        "schema": "vela.transportmodels.dg_srh_density_coupling_workflow.v1",
        "status": "complete" if len(stages) == 3 and all(
            stage["returncode"] == 0 for stage in stages
        ) else "fail",
        "stages": stages,
        "quantum_contract": {
            "source": str(BASE_CONFIG.resolve()),
            "srh_density_coupling": "sentaurus_default",
        },
    }
    (full.RUN_DIR / "workflow_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    curves = [full.aligned_curve("idvg"), full.aligned_curve("idvd")]
    completed_points = sum(curve["completed_points"] for curve in curves)
    artifacts: dict[str, str] = {}
    if completed_points:
        png, svg = full.make_plot(curves)
        artifacts = {"png": str(png.resolve()), "svg": str(svg.resolve())}
    report = {
        "schema": "vela.transportmodels.dg_srh_density_coupling_regression.v1",
        "as_of": "2026-08-23",
        "execution_status": "complete" if completed_points == 42 else "partial",
        "completed_points": completed_points,
        "expected_points": 42,
        "runner": {"path": str(runner.resolve()), "sha256": full.sha256(runner)},
        "manifest": str((full.RUN_DIR / "workflow_manifest.json").resolve()),
        "quantum_contract": manifest.get("quantum_contract", {}),
        "srh_density_coupling": "sentaurus_default",
        "curves": curves,
        "continuation_history": [],
        "prior_baseline": full.prior_baseline_metrics(),
        "artifacts": artifacts,
    }
    report["acceptance"] = full.acceptance(curves)
    if completed_points:
        full.write_report(report)
    return {
        "phase": "full_curves",
        "execution_status": report["execution_status"],
        "completed_points": completed_points,
        "main_curve_pass": report["acceptance"]["main_curve_pass"],
        "report": str(full.REPORT_JSON.resolve()),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, default=RUNNER)
    parser.add_argument(
        "--phase", choices=("three", "five", "full", "all"), default="all"
    )
    args = parser.parse_args()
    runner = args.runner.resolve()
    if not runner.is_file():
        raise FileNotFoundError(runner)
    OUTPUT.mkdir(parents=True, exist_ok=True)
    phases = ("three", "five", "full") if args.phase == "all" else (args.phase,)
    results = []
    for phase in phases:
        if phase == "three":
            results.append(run_three_points(runner))
        elif phase == "five":
            results.append(run_five_points(runner))
        else:
            results.append(run_full_curves(runner))
        (OUTPUT / "progress.json").write_text(
            json.dumps({"completed_phases": results}, indent=2) + "\n",
            encoding="utf-8",
        )
    print(json.dumps({"completed_phases": results}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
