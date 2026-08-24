#!/usr/bin/env python3
"""Recompute DD and issue a unified DD/DG fixed-contract acceptance report."""

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
UNCONFIGURED_ARTIFACT_ROOT = REPO / "__transportmodels_artifact_root_required__"
REF = UNCONFIGURED_ARTIFACT_ROOT
GENERATED = REF / "__unconfigured_generated__"
OUTPUT = REF / "__unconfigured_output__"
RUN_DIR = OUTPUT / "runs/dd"
WORKFLOW_SCRIPT = REPO / "scripts/run_transportmodels_dd_dg_workflow.py"
DG_REPORT = (
    REPO / "docs/validation/transportmodels_dg_srh_density_coupling_2026-08-23.json"
)
DG_RUN = REF / "__unconfigured_dg_run__"
REPORT_JSON = (
    REPO / "docs/validation/transportmodels_dd_dg_fixed_contract_v1_2026-08-24.json"
)
REPORT_MD = (
    REPO / "docs/validation/transportmodels_dd_dg_fixed_contract_v1_2026-08-24.md"
)
DEFAULT_RUNNER = REPO / "build-release/vela_example_runner.exe"


def configure_artifact_paths(root: Path) -> None:
    """Bind generated regression inputs to one explicit artifact bundle."""
    global REF, GENERATED, OUTPUT, RUN_DIR, DG_RUN
    REF = root.resolve()
    GENERATED = (
        REF / "vela_baseline/dd_dg_srh_corrected_cold_regression_2026-08-23"
        / "generated_corrected"
    )
    OUTPUT = REF / "vela_baseline/dd_dg_fixed_contract_v1_2026-08-24"
    RUN_DIR = OUTPUT / "runs/dd"
    DG_RUN = (
        REF / "vela_baseline/dg_srh_density_coupling_sentaurus_default_2026-08-23"
        / "full_21_point_curves/runs/dg"
    )


def validate_artifact_bundle() -> None:
    required = [
        GENERATED,
        REF / "run02/normalized/dd_idvg.csv",
        REF / "run02/normalized/dd_idvd.csv",
        DG_RUN / "dg_idvg_curve.json",
        DG_RUN / "dg_idvd_curve.json",
    ]
    missing = [str(path) for path in required if not path.exists()]
    if missing:
        raise FileNotFoundError(
            "TransportModels artifact bundle is incomplete:\n  "
            + "\n  ".join(missing)
        )


def load_workflow():
    spec = importlib.util.spec_from_file_location(
        "transportmodels_fixed_workflow", WORKFLOW_SCRIPT
    )
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load {WORKFLOW_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def environment() -> dict[str, str]:
    result = os.environ.copy()
    result["PATH"] = os.pathsep.join(
        [r"D:\msys64\ucrt64\bin", r"D:\msys64\usr\bin",
         result.get("PATH", "")]
    )
    return result


def patch_dd_configs(workflow, manifest: dict[str, Any]) -> None:
    violations: dict[str, list[str]] = {}
    for stage in manifest["stages"]:
        path = Path(stage["config"])
        config = json.loads(path.read_text(encoding="utf-8"))
        config = fixed.apply_contract(config, "dd")
        config["solver"]["verbose"] = False
        diagnostics = config["sweep"].setdefault("diagnostics", {})
        diagnostics["srh_balance"] = {
            "enabled": True,
            "material": "Si",
            "drain_contact": "drain",
            "substrate_contact": "substrate",
            "kcl_contacts": ["source", "drain", "gate", "substrate"],
            "resolution_margin_ratio": 10.0,
            "csv_file": str((RUN_DIR / f"{stage['name']}_srh_balance.csv").resolve()),
        }
        config["_comment"] = (
            "Frozen TransportModels DD baseline: corrected Sentaurus 2022 material, "
            "Fermi/OldSlotboom, sentaurus_default SRH density coupling, exact bias contract"
        )
        path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
        stage["config_sha256"] = workflow.sha256(path)
        problem = fixed.validate_config(config, "dd")
        if problem:
            violations[stage["name"]] = problem
    if violations:
        raise RuntimeError(f"Fixed-contract configuration violations: {violations}")
    manifest["fixed_contract"] = {
        "id": fixed.load_contract()["contract_id"],
        "path": str(fixed.DEFAULT_CONTRACT.resolve()),
        "sha256": fixed.sha256(fixed.DEFAULT_CONTRACT),
        "materials": str(fixed.materials_path()),
        "materials_sha256": fixed.sha256(fixed.materials_path()),
        "all_stage_configs_valid": True,
    }
    workflow.write_json(RUN_DIR / "workflow_manifest.json", manifest)


def materialize_dd() -> tuple[Any, dict[str, Any]]:
    workflow = load_workflow()
    manifest = workflow.materialize(GENERATED, RUN_DIR, ["dd"])
    patch_dd_configs(workflow, manifest)
    return workflow, manifest


def bias_tag(value: float) -> str:
    prefix = "m" if value < 0 else ""
    return prefix + f"{abs(value):.6f}".replace(".", "p")


def resume_idvg_with_internal_bridge(workflow) -> dict[str, Any]:
    """Resume a failed DD transition while keeping the 21-point report lattice."""
    failed_curve = RUN_DIR / "dd_idvg_curve.csv"
    converged = [row for row in read_csv(failed_curve)
                 if row.get("converged") == "1"]
    prefix = OUTPUT / "dd_idvg_completed_prefix.csv"
    old_prefix = read_csv(prefix) if prefix.is_file() else []
    if converged:
        restart_bias = round(float(converged[-1]["bias_V"]), 12)
    elif old_prefix:
        restart_bias = round(float(old_prefix[-1]["bias_V"]), 12)
    else:
        raise RuntimeError("DD Id-Vg failed before producing a restart state")
    restart_state = RUN_DIR / (
        f"dd_idvg_curve_state_bias_{bias_tag(restart_bias)}.csv"
    )
    if not restart_state.is_file():
        raise FileNotFoundError(restart_state)
    relax = [row for row in read_csv(RUN_DIR / "dd_idvg_final_bias_relax.csv")
             if row.get("converged") == "1"]
    if not relax:
        raise RuntimeError("Missing converged DD Id-Vg -1 V seed")
    prefix_rows = old_prefix or [{
        "bias_V": -1.0,
        "current_total_A_per_um": float(relax[-1]["current_total_A_per_um"]),
    }]
    if converged:
        prefix_rows = prefix_rows[:1]
        prefix_rows.extend({
            "bias_V": round(float(row["bias_V"]), 12),
            "current_total_A_per_um": float(row["current_total_A_per_um"]),
        } for row in converged)
    with prefix.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["bias_V", "current_total_A_per_um"]
        )
        writer.writeheader()
        writer.writerows(prefix_rows)
    archive = RUN_DIR / "workflow_manifest_failed_idvg_initial.json"
    archive.write_text(
        (RUN_DIR / "workflow_manifest.json").read_text(encoding="utf-8"),
        encoding="utf-8",
    )
    manifest = workflow.materialize(
        GENERATED, RUN_DIR, ["dd"],
        idvg_curve_restarts={"dd": {
            "bias_V": restart_bias,
            "state": restart_state,
            "prefix": prefix,
        }},
    )
    patch_dd_configs(workflow, manifest)
    curve_stage = next(
        stage for stage in manifest["stages"] if stage["name"] == "dd_idvg_curve"
    )
    path = Path(curve_stage["config"])
    config = json.loads(path.read_text(encoding="utf-8"))
    # The direct -0.36 -> -0.20 V step reaches a line-search non-decrease in
    # the weak-inversion transition. These 2.5 mV internal continuation
    # points are not comparison outputs; prepare_comparison_candidate keeps
    # only the exact frozen 21-point lattice.
    first_exact = float(curve_stage["bias_points"][0])
    bridge_count = int(round((first_exact - restart_bias) / 0.0025))
    bridges = [round(restart_bias + 0.0025 * index, 12)
               for index in range(1, bridge_count)]
    all_biases = bridges + [float(value) for value in curve_stage["bias_points"]]
    config["sweep"]["bias_points"] = all_biases
    config["sweep"]["start"] = all_biases[0]
    config["sweep"]["stop"] = all_biases[-1]
    config["sweep"]["step"] = 0.0025
    config["solver"].update({
        "line_search_mode": "block_filter",
        "residual_filter_gamma": 1.0e-4,
        "residual_filter_envelope_factor": 2.0,
        "continuity_row_scaling": {
            "enabled": True,
            "flux_fraction": 1.0e-3,
            "scale_floor": 1.0e-30,
            "min_source_scale": 1.0e-18,
            "min_weight": 1.0e-12,
            "max_weight": 1.0e12,
        },
    })
    config["_comment"] += (
        "; internal 2.5 mV transition bridge with block-filter globalization, "
        "exact 21-point report lattice unchanged"
    )
    path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    curve_stage["config_sha256"] = workflow.sha256(path)
    manifest["idvg_internal_continuation"] = {
        "restart_bias_V": restart_bias,
        "first_exact_remaining_bias_V": first_exact,
        "internal_step_V": 0.0025,
        "internal_bridge_points": len(bridges),
        "reported_bias_lattice_unchanged": True,
        "prefix": str(prefix.resolve()),
        "restart_state": str(restart_state.resolve()),
    }
    workflow.write_json(RUN_DIR / "workflow_manifest.json", manifest)
    return manifest


def complete_idvg_from_transition_checkpoint(
    runner: Path, workflow
) -> dict[str, Any]:
    """Reclose selected checkpoint states and finish the on-state tail."""
    base_path = RUN_DIR / "00_dd_idvg_curve.json"
    base = json.loads(base_path.read_text(encoding="utf-8"))
    solver = base["solver"]
    for key in (
        "line_search_mode", "residual_filter_gamma",
        "residual_filter_envelope_factor", "continuity_row_scaling",
    ):
        solver.pop(key, None)
    state_by_bias = {
        -1.0: RUN_DIR / "dd_idvg_final_bias_relax_final_state.csv",
        -0.84: RUN_DIR / "dd_idvg_curve_state_bias_m0p840000.csv",
        -0.68: RUN_DIR / "dd_idvg_curve_state_bias_m0p680000.csv",
        -0.20: RUN_DIR / "dd_idvg_curve_state_bias_m0p200000.csv",
        -0.04: RUN_DIR / "dd_idvg_curve_state_bias_m0p040000.csv",
        0.12: RUN_DIR / "dd_idvg_curve_state_bias_0p120000.csv",
    }
    checkpoint_root = OUTPUT / "idvg_checkpoint_reclosure"
    saved_prefix = read_csv(OUTPUT / "dd_idvg_completed_prefix.csv")
    saved_prefix_by_bias = {
        round(float(row["bias_V"]), 12): row for row in saved_prefix
    }
    exact_rows: dict[float, dict[str, str]] = {}
    exact_srh: dict[float, dict[str, str]] = {}
    reclosure_failures: list[float] = []
    for bias, state in state_by_bias.items():
        if not state.is_file():
            raise FileNotFoundError(state)
        tag = bias_tag(bias)
        run_dir = checkpoint_root / tag
        run_dir.mkdir(parents=True, exist_ok=True)
        config = json.loads(json.dumps(base))
        for contact in config["contacts"]:
            if contact["name"] == "gate":
                contact["bias"] = bias
        config["output_csv"] = str((run_dir / "curve.csv").resolve())
        config["log_file"] = str((run_dir / "curve.log").resolve())
        config["sweep"].update({
            "start": bias, "stop": bias, "step": 0.01,
            "bias_points": [bias],
            "initial_state_file": str(state.resolve()),
            "write_state_file": str((run_dir / "final_state.csv").resolve()),
            "write_state_every_point_prefix": str((run_dir / "state").resolve()),
        })
        config["sweep"].pop("initialization", None)
        config["sweep"].setdefault("diagnostics", {})["srh_balance"] = {
            "enabled": True, "material": "Si", "drain_contact": "drain",
            "substrate_contact": "substrate",
            "kcl_contacts": ["source", "drain", "gate", "substrate"],
            "resolution_margin_ratio": 10.0,
            "csv_file": str((run_dir / "srh_balance.csv").resolve()),
        }
        config["solver"].update({
            "line_search_mode": "block_filter",
            "residual_filter_gamma": 1.0e-4,
            "residual_filter_envelope_factor": 2.0,
            "continuity_row_scaling": {
                "enabled": True, "flux_fraction": 1.0e-3,
                "scale_floor": 1.0e-30, "min_source_scale": 1.0e-18,
                "min_weight": 1.0e-12, "max_weight": 1.0e12,
            },
        })
        config_path = run_dir / "config.json"
        config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
        completed = subprocess.run(
            [str(runner), "--config", str(config_path)], cwd=REPO,
            env=environment(), text=True, capture_output=True, check=False,
        )
        (run_dir / "console.log").write_text(
            completed.stdout + completed.stderr, encoding="utf-8"
        )
        rows = [row for row in read_csv(run_dir / "curve.csv")
                if row.get("converged") == "1"]
        bias_key = round(bias, 12)
        if completed.returncode != 0 or not rows:
            if bias_key not in saved_prefix_by_bias:
                raise RuntimeError(f"DD checkpoint reclosure failed at Vg={bias}")
            # Preserve the current from the original converged sweep state.
            # Missing KCL evidence intentionally makes the deep-off hard gate
            # unresolved in the unified report; it is not silently passed.
            exact_rows[bias_key] = saved_prefix_by_bias[bias_key]
            reclosure_failures.append(bias)
            continue
        exact_rows[bias_key] = rows[-1]
        srh_rows = read_csv(run_dir / "srh_balance.csv")
        if srh_rows:
            exact_srh[bias_key] = srh_rows[-1]

    tail_dir = OUTPUT / "idvg_on_tail"
    tail_dir.mkdir(parents=True, exist_ok=True)
    tail_biases = fixed.load_contract()["bias_contract"]["idvg"]["gate_bias_V"][8:]
    tail = json.loads(json.dumps(base))
    for contact in tail["contacts"]:
        if contact["name"] == "gate":
            contact["bias"] = tail_biases[0]
    tail["output_csv"] = str((tail_dir / "curve.csv").resolve())
    tail["log_file"] = str((tail_dir / "curve.log").resolve())
    tail["sweep"].update({
        "start": tail_biases[0], "stop": tail_biases[-1],
        "step": tail_biases[1] - tail_biases[0],
        "bias_points": tail_biases,
        "initial_state_file": str(state_by_bias[0.12].resolve()),
        "write_state_file": str((tail_dir / "final_state.csv").resolve()),
        "write_state_every_point_prefix": str((tail_dir / "state").resolve()),
    })
    tail["sweep"].pop("initialization", None)
    tail["sweep"].setdefault("diagnostics", {})["srh_balance"] = {
        "enabled": True, "material": "Si", "drain_contact": "drain",
        "substrate_contact": "substrate",
        "kcl_contacts": ["source", "drain", "gate", "substrate"],
        "resolution_margin_ratio": 10.0,
        "csv_file": str((tail_dir / "srh_balance.csv").resolve()),
    }
    tail_path = tail_dir / "config.json"
    tail_path.write_text(json.dumps(tail, indent=2) + "\n", encoding="utf-8")
    completed = subprocess.run(
        [str(runner), "--config", str(tail_path)], cwd=REPO,
        env=environment(), text=True, capture_output=True, check=False,
    )
    (tail_dir / "console.log").write_text(
        completed.stdout + completed.stderr, encoding="utf-8"
    )
    tail_rows = [row for row in read_csv(tail_dir / "curve.csv")
                 if row.get("converged") == "1"]
    if completed.returncode != 0 or len(tail_rows) != len(tail_biases):
        raise RuntimeError("DD Id-Vg on-state tail failed")
    for row in tail_rows:
        exact_rows[round(float(row["bias_V"]), 12)] = row
    tail_srh = read_csv(tail_dir / "srh_balance.csv")
    for row in tail_srh:
        exact_srh[round(float(row["bias_V"]), 12)] = row

    for row in saved_prefix:
        exact_rows.setdefault(round(float(row["bias_V"]), 12), row)
    expected = fixed.load_contract()["bias_contract"]["idvg"]["gate_bias_V"]
    missing = [bias for bias in expected if round(bias, 12) not in exact_rows]
    if missing:
        raise RuntimeError(f"DD Id-Vg checkpoint assembly misses {missing}")
    merged = [exact_rows[round(bias, 12)] for bias in expected]
    candidate = RUN_DIR / "dd_idvg_curve_comparison_candidate.csv"
    fields: list[str] = []
    for row in merged:
        fields.extend(key for key in row if key not in fields)
    with candidate.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(merged)
    srh_path = RUN_DIR / "dd_idvg_curve_srh_balance.csv"
    srh_rows = [exact_srh[key] for key in sorted(exact_srh)]
    srh_fields = list(srh_rows[0])
    with srh_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=srh_fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(srh_rows)
    return {
        "status": "pass", "reported_points": len(merged),
        "checkpoint_reclosures": len(state_by_bias),
        "checkpoint_reclosure_failures": reclosure_failures,
        "on_tail_points": len(tail_rows),
        "candidate": str(candidate.resolve()),
    }


def numerical_diagnostics() -> dict[float, dict[str, Any]]:
    rows: dict[float, dict[str, Any]] = {}
    for stem in ("dd_idvg_final_bias_relax", "dd_idvg_curve"):
        path = RUN_DIR / f"{stem}_srh_balance.csv"
        if not path.is_file():
            continue
        for row in read_csv(path):
            bias = round(float(row["bias_V"]), 12)
            rows[bias] = {
                "numerical_status": row.get("numerical_status", "unknown"),
                "four_terminal_kcl_residual_A_per_um": abs(float(
                    row["four_terminal_kcl_residual_A_per_um"]
                )),
                "id_to_kcl_residual_ratio": float(
                    row["id_to_kcl_residual_ratio"]
                ),
            }
    return rows


def percentile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * fraction
    low, high = math.floor(position), math.ceil(position)
    if low == high:
        return ordered[low]
    return ordered[low] * (high - position) + ordered[high] * (position - low)


def aligned_dd_curve(name: str) -> dict[str, Any]:
    candidate = RUN_DIR / f"dd_{name}_curve_comparison_candidate.csv"
    reference_path = REF / "run02/normalized" / f"dd_{name}.csv"
    reference = {
        round(float(row["bias_V"]), 12): abs(float(row["current_total"]))
        for row in read_csv(reference_path)
    }
    diagnostics = numerical_diagnostics() if name == "idvg" else {}
    aligned: list[dict[str, Any]] = []
    for row in read_csv(candidate):
        bias = round(float(row["bias_V"]), 12)
        vela = abs(float(row["current_total_A_per_um"]))
        sentaurus = reference[bias]
        aligned.append({
            "bias_V": bias,
            "vela_A_per_um": vela,
            "sentaurus_A_per_um": sentaurus,
            "absolute_relative_error": abs(vela - sentaurus) / max(sentaurus, 1e-300),
            "absolute_log_error_dex": abs(
                math.log10(max(vela, 1e-300)) - math.log10(max(sentaurus, 1e-300))
            ),
            **diagnostics.get(bias, {}),
        })
    aligned.sort(key=lambda row: row["bias_V"])
    if len(aligned) != 21:
        raise RuntimeError(f"DD {name} has {len(aligned)} aligned points, expected 21")
    if name == "idvg":
        regions = {"deep_off": aligned[:3], "transition": aligned[3:8], "on": aligned[8:]}
        metrics: dict[str, Any] = {
            key: {
                "max_relative_error": max(row["absolute_relative_error"] for row in values),
                "max_absolute_log_error_dex": max(
                    row["absolute_log_error_dex"] for row in values
                ),
                "median_absolute_log_error_dex": percentile(
                    [row["absolute_log_error_dex"] for row in values], 0.5
                ),
            }
            for key, values in regions.items()
        }
    else:
        nonzero = [row for row in aligned if row["bias_V"] > 0]
        metrics = {
            "max_relative_error": max(row["absolute_relative_error"] for row in nonzero),
            "median_relative_error": percentile(
                [row["absolute_relative_error"] for row in nonzero], 0.5
            ),
            "endpoint_relative_error": aligned[-1]["absolute_relative_error"],
        }
    aligned_path = OUTPUT / f"dd_{name}_aligned.csv"
    with aligned_path.open("w", newline="", encoding="utf-8") as handle:
        fields = list(aligned[0])
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(aligned)
    return {
        "branch": "dd", "name": name, "completed_points": len(aligned),
        "candidate_csv": str(candidate.resolve()),
        "reference_csv": str(reference_path.resolve()),
        "aligned_csv": str(aligned_path.resolve()),
        "aligned": aligned, "metrics": metrics,
    }


def branch_acceptance(curves: list[dict[str, Any]]) -> dict[str, Any]:
    limits = fixed.load_contract()["acceptance"]
    idvg = next(row for row in curves if row["name"] == "idvg")
    idvd = next(row for row in curves if row["name"] == "idvd")
    gates = {
        "idvg_transition": idvg["metrics"]["transition"]["max_absolute_log_error_dex"]
        <= limits["idvg_transition_max_absolute_log_error_dex"],
        "idvg_on": idvg["metrics"]["on"]["max_relative_error"]
        <= limits["idvg_on_max_absolute_relative_error"],
        "idvd_full": idvd["metrics"]["max_relative_error"]
        <= limits["idvd_max_absolute_relative_error"],
        "idvd_endpoint": idvd["metrics"]["endpoint_relative_error"]
        <= limits["idvd_2V_max_absolute_relative_error"],
    }
    deep = []
    for row in idvg["aligned"][:3]:
        resolved = (
            row.get("id_to_kcl_residual_ratio", 0.0)
            >= limits["deep_off_min_id_to_kcl_ratio"]
            and row.get("numerical_status") != "numerically_unresolved"
        )
        log_pass = row["absolute_log_error_dex"] <= limits[
            "idvg_transition_max_absolute_log_error_dex"
        ]
        deep.append({
            "bias_V": row["bias_V"],
            "absolute_relative_error": row["absolute_relative_error"],
            "absolute_log_error_dex": row["absolute_log_error_dex"],
            "id_to_kcl_residual_ratio": row.get("id_to_kcl_residual_ratio"),
            "status": "pass" if resolved and log_pass else
                      "numerically_unresolved" if not resolved else "fail",
        })
    return {
        "main_gates": gates,
        "main_curve_pass": all(gates.values()),
        "deep_off": {
            "policy": "log-current error plus Id/KCL >= 10",
            "points": deep,
            "pass": all(row["status"] == "pass" for row in deep),
        },
    }


def verify_frozen_dg(runner: Path) -> tuple[dict[str, Any], dict[str, Any]]:
    report = json.loads(DG_REPORT.read_text(encoding="utf-8"))
    contract = fixed.load_contract()
    decks = [DG_RUN / "dg_idvg_curve.json", DG_RUN / "dg_idvd_curve.json"]
    violations: list[str] = []
    for path in decks:
        config = json.loads(path.read_text(encoding="utf-8"))
        solver = config["solver"]
        for key, value in contract["common_solver_physics"].items():
            if solver.get(key) != value:
                violations.append(f"{path.name}: solver.{key}")
        if solver.get("electron_quantum_potential") != contract["dg_quantum_contract"]:
            violations.append(f"{path.name}: electron_quantum_potential")
        old_materials = json.loads(
            Path(config["materials_file"]).read_text(encoding="utf-8")
        )["materials"]
        frozen_materials = json.loads(
            fixed.materials_path().read_text(encoding="utf-8")
        )["materials"]
        if old_materials != frozen_materials:
            violations.append(f"{path.name}: semantic material payload")
    expected = contract["bias_contract"]
    idvg = json.loads(decks[0].read_text(encoding="utf-8"))["sweep"]["bias_points"]
    idvd = json.loads(decks[1].read_text(encoding="utf-8"))["sweep"]["bias_points"]
    if idvg != expected["idvg"]["gate_bias_V"]:
        violations.append("dg_idvg_curve.json: bias lattice")
    if idvd != expected["idvd"]["drain_bias_V"]:
        violations.append("dg_idvd_curve.json: bias lattice")
    current_runner_sha256 = fixed.sha256(runner)
    evidence_runner_sha256 = report.get("runner", {}).get("sha256")
    if evidence_runner_sha256 != current_runner_sha256:
        violations.append(
            "frozen DG evidence runner SHA-256 does not match the current runner"
        )
    evidence = {
        "report": str(DG_REPORT.resolve()),
        "report_sha256": fixed.sha256(DG_REPORT),
        "completed_points": report["completed_points"],
        "main_curve_pass": report["acceptance"]["main_curve_pass"],
        "deep_off_pass": all(
            row["status"] == "pass"
            for row in report["acceptance"]["deep_off"]["points"]
        ),
        "configuration_contract_violations": violations,
        "configuration_contract_pass": not violations,
        "runner": report["runner"],
        "current_runner_sha256": current_runner_sha256,
        "runner_matches_current": evidence_runner_sha256 == current_runner_sha256,
    }
    return report, evidence


def make_plot(dd_curves: list[dict[str, Any]], dg_curves: list[dict[str, Any]]) -> dict[str, str]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(12.4, 8.4))
    for branch, curves, color in (("DD", dd_curves, "tab:blue"),
                                  ("DG", dg_curves, "tab:orange")):
        idvg = next(row for row in curves if row["name"] == "idvg")["aligned"]
        idvd = next(row for row in curves if row["name"] == "idvd")["aligned"]
        axes[0, 0].semilogy([row["bias_V"] for row in idvg],
                            [row["sentaurus_A_per_um"] for row in idvg],
                            "-", color=color, label=f"Sentaurus {branch}")
        axes[0, 0].semilogy([row["bias_V"] for row in idvg],
                            [row["vela_A_per_um"] for row in idvg],
                            "o", ms=3.5, color=color, label=f"Vela {branch}")
        axes[0, 1].plot([row["bias_V"] for row in idvg],
                        [row["absolute_log_error_dex"] for row in idvg],
                        "o-", ms=3.5, color=color, label=branch)
        axes[1, 0].plot([row["bias_V"] for row in idvd],
                        [1e3 * row["sentaurus_A_per_um"] for row in idvd],
                        "-", color=color, label=f"Sentaurus {branch}")
        axes[1, 0].plot([row["bias_V"] for row in idvd],
                        [1e3 * row["vela_A_per_um"] for row in idvd],
                        "o", ms=3.5, color=color, label=f"Vela {branch}")
        nonzero = [row for row in idvd if row["bias_V"] > 0]
        axes[1, 1].plot([row["bias_V"] for row in nonzero],
                        [100 * row["absolute_relative_error"] for row in nonzero],
                        "o-", ms=3.5, color=color, label=branch)
    axes[0, 0].set(xlabel="Gate voltage Vg (V)", ylabel="Drain current Id (A/um)",
                   title="Id-Vg fixed-contract comparison")
    axes[0, 1].axhline(0.15, color="tab:red", ls="--", label="0.15 dex limit")
    axes[0, 1].set(xlabel="Gate voltage Vg (V)", ylabel="Absolute log error (dex)",
                   title="Id-Vg log error")
    axes[1, 0].set(xlabel="Drain voltage Vd (V)", ylabel="Drain current Id (mA/um)",
                   title="Id-Vd fixed-contract comparison")
    axes[1, 1].axhline(5, color="tab:red", ls="--", label="5% limit")
    axes[1, 1].set(xlabel="Drain voltage Vd (V)", ylabel="Absolute relative error (%)",
                   title="Id-Vd relative error")
    for axis in axes.flat:
        axis.grid(True, which="both", alpha=0.25)
        axis.legend(fontsize=8)
    fig.suptitle("TransportModels Sentaurus 2022 unified DD/DG regression")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    png = OUTPUT / "transportmodels_dd_dg_fixed_contract_comparison.png"
    svg = OUTPUT / "transportmodels_dd_dg_fixed_contract_comparison.svg"
    fig.savefig(png, dpi=180)
    fig.savefig(svg)
    plt.close(fig)
    return {"png": str(png.resolve()), "svg": str(svg.resolve())}


def write_report(report: dict[str, Any]) -> None:
    payload, nonfinite_paths = fixed.strict_json_payload(report)
    if nonfinite_paths:
        payload["serialization"] = {
            "nonfinite_values_replaced_with_null": nonfinite_paths,
            "reason": "source evidence did not provide a finite value",
        }
    REPORT_JSON.write_text(
        json.dumps(payload, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    dd = report["branches"]["dd"]
    dg = report["branches"]["dg"]
    rows = []
    for name, branch in (("DD", dd), ("DG", dg)):
        idvg = next(row for row in branch["curves"] if row["name"] == "idvg")
        idvd = next(row for row in branch["curves"] if row["name"] == "idvd")
        if branch["acceptance"]["overall_pass"]:
            result_text = "通过"
        elif branch["acceptance"]["main_curve_pass"]:
            result_text = "主曲线通过；深关断未解析"
        else:
            result_text = "主曲线未通过"
        rows.append(
            f"| {name} | {idvg['metrics']['transition']['max_absolute_log_error_dex']:.6f} dex | "
            f"{100 * idvg['metrics']['on']['max_relative_error']:.3f}% | "
            f"{100 * idvd['metrics']['max_relative_error']:.3f}% | "
            f"{100 * idvd['metrics']['endpoint_relative_error']:.3f}% | "
            f"{result_text} |"
        )
    deep_rows = []
    for name, branch in (("DD", dd), ("DG", dg)):
        for row in branch["acceptance"]["deep_off"]["points"]:
            ratio = row.get("id_to_kcl_residual_ratio")
            ratio_text = "—" if ratio is None else f"{ratio:.3f}"
            relative = row.get("absolute_relative_error")
            relative_text = "—" if relative is None else f"{100 * relative:.3f}%"
            deep_rows.append(
                f"| {name} | {row['bias_V']:.2f} | "
                f"{relative_text} | "
                f"{row['absolute_log_error_dex']:.6f} dex | "
                f"{ratio_text} | {row['status']} |"
            )
    markdown = f"""# TransportModels DD/DG 固定契约统一验收（2026-08-24）

## 结论

本次固定了 Sentaurus 2022 材料参数、Fermi/OldSlotboom、`sentaurus_default` SRH 密度耦合、偏压点阵及最终 DG 量子势契约。DD 使用同一材料与偏压契约重新计算 21 点 Id–Vg 和 21 点 Id–Vd；DG 使用已完成的 42 点证据并通过契约逐字段审计。

统一验收：**{'通过' if report['acceptance']['overall_pass'] else '未通过'}**。

## 主曲线指标

| 分支 | Id–Vg 过渡区最大对数误差 | Id–Vg 导通区最大相对误差 | Id–Vd 最大相对误差 | Id–Vd 2 V 误差 | 结果 |
|---|---:|---:|---:|---:|---|
{chr(10).join(rows)}

## 深关断前三点

| 分支 | Vg (V) | 相对误差 | 对数误差 | Id/abs(KCL) | 状态 |
|---|---:|---:|---:|---:|---|
{chr(10).join(deep_rows)}

深关断采用独立门槛：对数误差不超过 0.15 dex，且 `Id/abs(KCL) >= 10`；不与导通区相对误差混用。

## 固定基线

- 契约：`{report['contract']['path']}`
- 契约 SHA-256：`{report['contract']['sha256']}`
- 材料 SHA-256：`{report['contract']['materials_sha256']}`
- DD/DG 受控差分：{'通过，仅 electron_quantum_potential 不同' if report['controlled_delta']['pass'] else '未通过'}
- DG 配置逐字段审计：{'通过' if report['branches']['dg']['evidence']['configuration_contract_pass'] else '未通过'}
- DD 运行器 SHA-256：`{report['branches']['dd']['runner']['sha256']}`
- DG 证据运行器 SHA-256：`{report['branches']['dg']['evidence']['runner']['sha256']}`

## 产物

- JSON 报告：`{REPORT_JSON.resolve()}`
- 对比图：`{report['artifacts']['png']}`
- DD 工作流清单：`{report['branches']['dd']['manifest']}`
"""
    REPORT_MD.write_text(markdown, encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, default=DEFAULT_RUNNER)
    parser.add_argument(
        "--artifact-root",
        type=Path,
        default=None,
        help=(
            "Mounted TransportModels artifact bundle. Alternatively set "
            f"{fixed.ARTIFACT_ROOT_ENV}."
        ),
    )
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--report-only", action="store_true")
    args = parser.parse_args()
    if args.report_only:
        if not REPORT_JSON.is_file():
            raise FileNotFoundError(REPORT_JSON)
        write_report(json.loads(REPORT_JSON.read_text(encoding="utf-8")))
        print(str(REPORT_MD.resolve()))
        return 0
    try:
        configure_artifact_paths(fixed.resolve_artifact_root(args.artifact_root))
        validate_artifact_bundle()
    except (ValueError, FileNotFoundError) as error:
        parser.error(str(error))
    runner = args.runner.resolve()
    if not runner.is_file():
        raise FileNotFoundError(runner)
    OUTPUT.mkdir(parents=True, exist_ok=True)
    transition_checkpoint = RUN_DIR / "dd_idvg_curve_state_bias_0p120000.csv"
    cached_failed_state = RUN_DIR / "dd_idvg_curve_state_bias_m0p360000.csv"
    idvg_checkpoint_completion = None
    if transition_checkpoint.is_file() and not REPORT_JSON.is_file():
        workflow = load_workflow()
        idvg_checkpoint_completion = complete_idvg_from_transition_checkpoint(
            runner, workflow
        )
        manifest = workflow.materialize(GENERATED, RUN_DIR, ["dd"])
        patch_dd_configs(workflow, manifest)
        manifest["stages"] = [
            stage for stage in manifest["stages"]
            if stage["name"].startswith("dd_idvd_")
        ]
        manifest["idvg_checkpoint_completion"] = idvg_checkpoint_completion
        workflow.write_json(RUN_DIR / "workflow_manifest.json", manifest)
    elif cached_failed_state.is_file() and not REPORT_JSON.is_file():
        workflow = load_workflow()
        manifest = resume_idvg_with_internal_bridge(workflow)
    else:
        workflow, manifest = materialize_dd()
    previous = None
    if args.resume and (RUN_DIR / "workflow_manifest.previous.json").is_file():
        previous = json.loads((RUN_DIR / "workflow_manifest.previous.json").read_text())
    manifest = workflow.execute(manifest, runner, RUN_DIR, previous)
    failed = next(
        (stage for stage in manifest["stages"] if stage.get("status") == "fail"),
        None,
    )
    if (failed and failed["name"] == "dd_idvg_curve"
            and not manifest.get("idvg_internal_continuation")):
        manifest = resume_idvg_with_internal_bridge(workflow)
        manifest = workflow.execute(manifest, runner, RUN_DIR)
    if manifest.get("status") != "pass":
        raise RuntimeError(f"DD workflow failed; inspect {RUN_DIR / 'workflow_manifest.json'}")
    if idvg_checkpoint_completion is not None:
        manifest["idvg_checkpoint_completion"] = idvg_checkpoint_completion
        workflow.write_json(RUN_DIR / "workflow_manifest.json", manifest)
    dd_curves = [aligned_dd_curve("idvg"), aligned_dd_curve("idvd")]
    dd_acceptance = branch_acceptance(dd_curves)
    dd_acceptance["overall_pass"] = (
        dd_acceptance["main_curve_pass"] and dd_acceptance["deep_off"]["pass"]
    )
    dg_report, dg_evidence = verify_frozen_dg(runner)
    dg_curves = dg_report["curves"]
    dg_acceptance = dg_report["acceptance"]
    dg_acceptance = json.loads(json.dumps(dg_acceptance))
    dg_acceptance["deep_off"]["pass"] = all(
        row["status"] == "pass" for row in dg_acceptance["deep_off"]["points"]
    )
    dg_acceptance["overall_pass"] = (
        dg_acceptance["main_curve_pass"] and dg_acceptance["deep_off"]["pass"]
        and dg_evidence["configuration_contract_pass"]
    )
    dd_deck = json.loads((RUN_DIR / "03_dd_idvg_curve.json").read_text())
    dg_deck = json.loads((DG_RUN / "dg_idvg_curve.json").read_text())
    dd_solver, dg_solver = fixed.controlled_solver_delta(dd_deck, dg_deck)
    controlled_delta = {
        "pass": dd_solver == dg_solver,
        "only_branch_delta": "solver.electron_quantum_potential",
    }
    artifacts = make_plot(dd_curves, dg_curves)
    report = {
        "schema": "vela.transportmodels.dd_dg.fixed_contract.acceptance.v1",
        "as_of": "2026-08-24",
        "contract": manifest["fixed_contract"],
        "controlled_delta": controlled_delta,
        "branches": {
            "dd": {
                "execution_status": "complete",
                "completed_points": 42,
                "runner": {"path": str(runner), "sha256": fixed.sha256(runner)},
                "manifest": str((RUN_DIR / "workflow_manifest.json").resolve()),
                "curves": dd_curves,
                "acceptance": dd_acceptance,
            },
            "dg": {
                "execution_status": "reused_verified_evidence",
                "completed_points": dg_report["completed_points"],
                "curves": dg_curves,
                "acceptance": dg_acceptance,
                "evidence": dg_evidence,
            },
        },
        "acceptance": {
            "overall_pass": controlled_delta["pass"]
            and dd_acceptance["overall_pass"] and dg_acceptance["overall_pass"]
        },
        "artifacts": artifacts,
    }
    write_report(report)
    print(json.dumps({
        "overall_pass": report["acceptance"]["overall_pass"],
        "dd_acceptance": dd_acceptance,
        "dg_acceptance": dg_acceptance,
        "report": str(REPORT_JSON.resolve()),
    }, indent=2))
    return 0 if report["acceptance"]["overall_pass"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
