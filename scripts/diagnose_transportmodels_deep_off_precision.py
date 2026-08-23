#!/usr/bin/env python3
"""Reproduce TransportModels deep-off current and DG residual diagnostics.

The script creates isolated, restart-based diagnostic runs.  It does not alter
the physical model or overwrite the strict-sweep baseline artifacts.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[1]
STRICT_ROOT = (
    REPO
    / "build-release/reference_tcad/transportmodels_sentaurus2022/vela_baseline"
    / "idvg_srh_strict_2026-08-21"
)
OUTPUT_ROOT = (
    REPO
    / "build-release/reference_tcad/transportmodels_sentaurus2022/reports"
    / "idvg_deep_off_precision_20260822"
)
RUNNER = REPO / "build-release/vela_example_runner.exe"
BIAS_LABELS = {
    -1.0: "m1p000000",
    -0.68: "m0p680000",
    -0.52: "m0p520000",
    -0.4: "m0p400000",
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def state_path(model: str, bias: float) -> Path:
    directory = "dg_forward" if model == "dg" and bias != -0.4 else (
        "dg_resume" if model == "dg" else "dd_forward"
    )
    return STRICT_ROOT / directory / f"state_bias_{BIAS_LABELS[bias]}.csv"


def base_config(model: str, bias: float) -> dict[str, Any]:
    directory = "dg_forward" if model == "dg" and bias != -0.4 else (
        "dg_resume" if model == "dg" else "dd_forward"
    )
    return json.loads((STRICT_ROOT / directory / "config.json").read_text())


def prepare_single_point(
    model: str,
    bias: float,
    *,
    qf_reference: bool = False,
    hard_gate: bool = False,
) -> tuple[Path, Path]:
    if hard_gate:
        qf_reference = True
    variant = "_hard_gate" if hard_gate else ("_qfref" if qf_reference else "")
    out = OUTPUT_ROOT / f"{model}_{BIAS_LABELS[bias]}{variant}"
    out.mkdir(parents=True, exist_ok=True)
    cfg = base_config(model, bias)
    cfg["_comment"] = (
        f"Deep-off precision diagnostic: {model.upper()} at Vg={bias} V; "
        "same physics as strict baseline"
        + ("; contact-majority quasi-Fermi reference enabled" if qf_reference else "")
        + ("; carrier rows and global continuity enforced" if hard_gate else "")
    )
    if qf_reference:
        cfg["solver"]["quasi_fermi_reference"] = "contact_majority"
    if hard_gate:
        cfg["solver"]["carrier_row_qualified_stall_acceptance"] = True
        cfg["solver"]["global_continuity_closure"] = {
            "mode": "enforce",
            "tolerance": 0.1,
            "source_floor": 1.0e-18,
        }
    cfg["output_csv"] = str(out / "curve.csv")
    cfg["sweep"].update(
        {
            "start": bias,
            "stop": bias,
            "step": 0.01,
            "bias_points": [bias],
            "initial_state_file": str(state_path(model, bias)),
            "write_state_file": str(out / "final_state.csv"),
            "write_state_every_point_prefix": str(out / "state"),
            "max_retries": 0,
        }
    )
    diagnostics = cfg["sweep"].setdefault("diagnostics", {})
    diagnostics["terminal_balance"] = {
        "enabled": True,
        "contacts": ["source", "drain", "gate", "substrate"],
        "csv_file": str(out / "terminal_balance.csv"),
    }
    diagnostics["contact_edge"] = {
        "enabled": True,
        "contacts": ["source", "drain", "substrate"],
        "csv_file": str(out / "contact_edges.csv"),
    }
    diagnostics["srh_balance"] = {
        "enabled": True,
        "material": "Si",
        "drain_contact": "drain",
        "substrate_contact": "substrate",
        "kcl_contacts": ["source", "drain", "gate", "substrate"],
        "resolution_margin_ratio": 10.0,
        "csv_file": str(out / "srh_balance.csv"),
    }
    diagnostics["newton_history"] = {
        "enabled": True,
        "csv_file": str(out / "newton_history.csv"),
    }
    row_cfg = cfg["solver"].setdefault("carrier_row_convergence", {})
    row_cfg["mode"] = "enforce" if hard_gate else "report"
    if hard_gate:
        row_cfg["eps_row"] = 1.0e-3
        row_cfg["min_source_scale_fraction"] = 0.0
        row_cfg["min_source_scale"] = 1.0e-18
    row_cfg["diagnostic_csv"] = str(out / "carrier_row_violations.csv")
    config_path = out / "config.json"
    config_path.write_text(json.dumps(cfg, indent=2) + "\n")
    return config_path, out


NEWTON_CALIBRATION_VARIANTS: tuple[tuple[str, dict[str, float]], ...] = (
    (
        "floor2e11_qf1e2",
        {
            "stall_residual_floor": 2.0e-11,
            "quasi_fermi_update_limit_V": 1.0e-2,
            "damping_factor": 1.0,
        },
    ),
    (
        "floor2e11_qf5e3",
        {
            "stall_residual_floor": 2.0e-11,
            "quasi_fermi_update_limit_V": 5.0e-3,
            "damping_factor": 1.0,
        },
    ),
    (
        "floor2e11_qf1e2_damp5e1",
        {
            "stall_residual_floor": 2.0e-11,
            "quasi_fermi_update_limit_V": 1.0e-2,
            "damping_factor": 5.0e-1,
        },
    ),
)


def prepare_newton_calibration(
    model: str, bias: float, variant: str, settings: dict[str, float]
) -> tuple[Path, Path]:
    baseline_config, _ = prepare_single_point(model, bias, hard_gate=True)
    cfg = json.loads(baseline_config.read_text())
    out = OUTPUT_ROOT / "newton_calibration" / variant / (
        f"{model}_{BIAS_LABELS[bias]}"
    )
    out.mkdir(parents=True, exist_ok=True)
    cfg["_comment"] += (
        f"; Newton calibration variant={variant}; hard convergence gates unchanged"
    )
    cfg["solver"].update(settings)
    cfg["solver"]["line_search"] = True
    cfg["output_csv"] = str(out / "curve.csv")
    cfg["solver"]["carrier_row_convergence"]["diagnostic_csv"] = str(
        out / "carrier_row_violations.csv"
    )
    sweep = cfg["sweep"]
    sweep["write_state_file"] = str(out / "final_state.csv")
    sweep["write_state_every_point_prefix"] = str(out / "state")
    diagnostics = sweep["diagnostics"]
    diagnostics["terminal_balance"]["csv_file"] = str(out / "terminal_balance.csv")
    diagnostics["contact_edge"]["csv_file"] = str(out / "contact_edges.csv")
    diagnostics["srh_balance"]["csv_file"] = str(out / "srh_balance.csv")
    diagnostics["newton_history"]["csv_file"] = str(out / "newton_history.csv")
    config_path = out / "config.json"
    config_path.write_text(json.dumps(cfg, indent=2) + "\n")
    return config_path, out


def summarize_newton_calibration(
    model: str,
    bias: float,
    variant: str,
    settings: dict[str, float],
    run: dict[str, Any],
    out: Path,
) -> dict[str, Any]:
    history = read_csv(out / "newton_history.csv") if (
        out / "newton_history.csv"
    ).exists() else []
    failure_path = out / "curve_newton_failure_diagnostics.json"
    failure: dict[str, Any] = {}
    if failure_path.exists():
        records = json.loads(failure_path.read_text())
        if records:
            failure = records[-1]
    terminal = read_csv(out / "terminal_balance.csv") if (
        out / "terminal_balance.csv"
    ).exists() else []
    carrier_rows = read_csv(out / "carrier_row_violations.csv") if (
        out / "carrier_row_violations.csv"
    ).exists() else []
    return {
        "model": model.upper(),
        "bias_V": bias,
        "variant": variant,
        **settings,
        "returncode": run["returncode"],
        "converged": run["returncode"] == 0,
        "failure_reason": failure.get("failure_reason", ""),
        "iterations": max((int(row["iteration"]) for row in history), default=0),
        "final_residual_norm": (
            float(history[-1]["residual_norm"]) if history else math.nan
        ),
        "final_carrier_row_violations": (
            int(history[-1]["carrier_row_violations"]) if history else -1
        ),
        "final_carrier_row_max_ratio": (
            float(history[-1]["carrier_row_max_ratio"]) if history else math.nan
        ),
        "reported_violation_rows": len(carrier_rows),
        "four_terminal_kcl_A_per_um": sum(
            finite(row["current_total_A_per_um"]) for row in terminal
        ) if terminal else math.nan,
    }


def prepare_dg_failure_probe() -> tuple[Path, Path]:
    start = -0.4
    target = -0.3987109375
    out = OUTPUT_ROOT / "dg_residual_probe"
    out.mkdir(parents=True, exist_ok=True)
    cfg = base_config("dg", start)
    cfg["_comment"] = "DG strict residual-platform reproduction from the -0.4 V checkpoint"
    cfg["output_csv"] = str(out / "curve.csv")
    cfg["sweep"].update(
        {
            "start": target,
            "stop": target,
            "step": target - start,
            "bias_points": [target],
            "initial_state_file": str(state_path("dg", start)),
            "write_state_file": str(out / "final_state.csv"),
            "write_state_every_point_prefix": str(out / "state"),
            "max_retries": 0,
        }
    )
    diagnostics = cfg["sweep"].setdefault("diagnostics", {})
    diagnostics["newton_history"] = {
        "enabled": True,
        "csv_file": str(out / "newton_history.csv"),
    }
    row_cfg = cfg["solver"].setdefault("carrier_row_convergence", {})
    row_cfg["mode"] = "report"
    row_cfg["diagnostic_csv"] = str(out / "carrier_row_violations.csv")
    config_path = out / "config.json"
    config_path.write_text(json.dumps(cfg, indent=2) + "\n")
    return config_path, out


def prepare_dg_row_scaling_probe() -> tuple[Path, Path]:
    config_path, _ = prepare_dg_failure_probe()
    cfg = json.loads(config_path.read_text())
    out = OUTPUT_ROOT / "dg_residual_probe_row_scaled"
    out.mkdir(parents=True, exist_ok=True)
    cfg["_comment"] += "; source-aware continuity row scaling enabled"
    cfg["output_csv"] = str(out / "curve.csv")
    cfg["solver"]["continuity_row_scaling"] = {
        "enabled": True,
        "flux_fraction": 0.0,
        "scale_floor": 1.0e-30,
        "min_source_scale": 1.0e-30,
        "min_weight": 1.0e-12,
        "max_weight": 1.0e18,
    }
    cfg["solver"]["carrier_row_convergence"]["diagnostic_csv"] = str(
        out / "carrier_row_violations.csv"
    )
    diagnostics = cfg["sweep"]["diagnostics"]
    diagnostics["newton_history"]["csv_file"] = str(out / "newton_history.csv")
    cfg["sweep"]["write_state_file"] = str(out / "final_state.csv")
    cfg["sweep"]["write_state_every_point_prefix"] = str(out / "state")
    result = out / "config.json"
    result.write_text(json.dumps(cfg, indent=2) + "\n")
    return result, out


def execute(
    config_path: Path, out: Path, env_extra: dict[str, str] | None = None
) -> dict[str, Any]:
    env = os.environ.copy()
    if env_extra:
        env.update(env_extra)
    completed = subprocess.run(
        [str(RUNNER), "--config", str(config_path), "--log", str(out / "runner.log")],
        cwd=REPO,
        text=True,
        capture_output=True,
        check=False,
        env=env,
    )
    (out / "stdout.txt").write_text(completed.stdout)
    (out / "stderr.txt").write_text(completed.stderr)
    return {
        "config": str(config_path),
        "returncode": completed.returncode,
        "stdout": str(out / "stdout.txt"),
        "stderr": str(out / "stderr.txt"),
    }


def prepare_dg_backend_probe(backend: str) -> tuple[Path, Path]:
    source, _ = prepare_dg_failure_probe()
    cfg = json.loads(source.read_text())
    out = OUTPUT_ROOT / f"dg_residual_probe_{backend}"
    out.mkdir(parents=True, exist_ok=True)
    cfg["_comment"] += f"; VELA_LINEAR_SOLVER={backend}"
    cfg["output_csv"] = str(out / "curve.csv")
    cfg["solver"]["carrier_row_convergence"]["diagnostic_csv"] = str(
        out / "carrier_row_violations.csv"
    )
    cfg["sweep"]["diagnostics"]["newton_history"]["csv_file"] = str(
        out / "newton_history.csv"
    )
    cfg["sweep"]["write_state_file"] = str(out / "final_state.csv")
    cfg["sweep"]["write_state_every_point_prefix"] = str(out / "state")
    result = out / "config.json"
    result.write_text(json.dumps(cfg, indent=2) + "\n")
    return result, out


def finite(value: str) -> float:
    result = float(value)
    return result if math.isfinite(result) else 0.0


def summarize_point(
    model: str, bias: float, *, qf_reference: bool = False
) -> dict[str, Any]:
    variant = "_qfref" if qf_reference else ""
    out = OUTPUT_ROOT / f"{model}_{BIAS_LABELS[bias]}{variant}"
    curve = read_csv(out / "curve.csv")[-1]
    edges = [
        row for row in read_csv(out / "contact_edges.csv")
        if row["current_contact"] == "drain"
    ]
    active_edges = [row for row in edges if finite(row["mun"]) > 0.0]
    terminal = read_csv(out / "terminal_balance.csv")
    srh = read_csv(out / "srh_balance.csv")[-1]
    electron = finite(curve["current_electron_A_per_um"])
    drift = finite(curve["current_electron_drift_A_per_um"])
    diffusion = finite(curve["current_electron_diffusion_A_per_um"])
    split_scale = abs(drift) + abs(diffusion)
    edge_qf_drops = [
        abs(finite(row["phin1"]) - finite(row["phin0"])) for row in active_edges
    ]
    edge_electron = [finite(row["current_electron"]) for row in active_edges]
    kcl = sum(finite(row["current_total_A_per_um"]) for row in terminal)
    source_terms = [
        row for row in read_csv(out / "carrier_row_violations.csv")
        if row["carrier"] == "electron"
    ] if (out / "carrier_row_violations.csv").exists() else []
    return {
        "model": model.upper(),
        "qf_reference": "contact_majority" if qf_reference else "none",
        "bias_V": bias,
        "Id_A_per_um": finite(curve["current_total_A_per_um"]),
        "electron_current_A_per_um": electron,
        "electron_drift_A_per_um": drift,
        "electron_diffusion_A_per_um": diffusion,
        "net_to_split_magnitude_ratio": abs(electron) / split_scale if split_scale else 0.0,
        "decimal_digits_cancelled": (
            -math.log10(abs(electron) / split_scale)
            if electron != 0.0 and split_scale > 0.0 else math.inf
        ),
        "drain_edge_count": len(edges),
        "active_drain_edge_count": len(active_edges),
        "zero_electron_edge_count": sum(value == 0.0 for value in edge_electron),
        "zero_qf_drop_edge_count": sum(value == 0.0 for value in edge_qf_drops),
        "min_nonzero_qf_drop_V": min(
            (value for value in edge_qf_drops if value > 0.0), default=0.0
        ),
        "max_qf_drop_V": max(edge_qf_drops, default=0.0),
        "four_terminal_kcl_A_per_um": kcl,
        "srh_generation_A_per_um": finite(srh["srh_generation_current_A_per_um"]),
        "substrate_hole_A_per_um": finite(srh["substrate_hole_current_A_per_um"]),
        "electron_row_violation_count": len(source_terms),
        "electron_row_max_ratio": max(
            (finite(row["ratio"]) for row in source_terms), default=0.0
        ),
    }


def incident_regions(mesh: dict[str, Any]) -> dict[int, set[str]]:
    region_names = {region["id"]: region["name"] for region in mesh["regions"]}
    result: dict[int, set[str]] = defaultdict(set)
    for cell in mesh["triangles"]:
        name = region_names[cell["region_id"]]
        for node in cell["node_ids"]:
            result[node].add(name)
    return result


def summarize_dg_failure() -> dict[str, Any]:
    original = json.loads(
        (STRICT_ROOT / "dg_resume/curve_newton_failure_diagnostics.json").read_text()
    )[-1]
    mesh_path = Path(base_config("dg", -0.4)["mesh_file"])
    mesh = json.loads(mesh_path.read_text())
    regions = incident_regions(mesh)
    hotspots = []
    for row in original["top_electron_residual_nodes"][:20]:
        hotspots.append(
            {
                **row,
                "incident_regions": sorted(regions.get(row["node_id"], set())),
            }
        )
    block = original["block_residuals"]
    return {
        "bias_V": original["bias_V"],
        "failure_reason": original["failure_reason"],
        "failed_iteration": original["failed_iteration"],
        "line_search_attempts": original["line_search_attempts"],
        "combined_residual": block["combined"],
        "electron_residual": block["phin"],
        "hole_residual": block["phip"],
        "poisson_residual": block["psi"],
        "electron_fraction_of_squared_norm": (
            block["phin"] ** 2 / block["combined"] ** 2
        ),
        "step_norm": original["step_norm"],
        "hotspots": hotspots,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--summarize-only", action="store_true")
    parser.add_argument("--qf-reference-ab", action="store_true")
    parser.add_argument("--hard-gate-ab", action="store_true")
    parser.add_argument("--row-scaling-probe", action="store_true")
    parser.add_argument("--linear-backend-ab", action="store_true")
    parser.add_argument("--newton-calibration", action="store_true")
    args = parser.parse_args()
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    runs: list[dict[str, Any]] = []
    selected = [(model, bias) for model in ("dd", "dg") for bias in (-1.0, -0.68, -0.52)]
    if args.newton_calibration:
        calibration_rows: list[dict[str, Any]] = []
        for variant, settings in NEWTON_CALIBRATION_VARIANTS:
            for model, bias in selected:
                config, out = prepare_newton_calibration(
                    model, bias, variant, settings
                )
                if args.prepare_only:
                    continue
                run = execute(config, out)
                calibration_rows.append(
                    summarize_newton_calibration(
                        model, bias, variant, settings, run, out
                    )
                )
        if args.prepare_only:
            return 0
        write_csv(
            OUTPUT_ROOT / "newton_calibration_summary.csv", calibration_rows
        )
        (OUTPUT_ROOT / "newton_calibration_summary.json").write_text(
            json.dumps(calibration_rows, indent=2) + "\n"
        )
        print(OUTPUT_ROOT / "newton_calibration_summary.csv")
        return 0
    if args.linear_backend_ab:
        backend_runs = []
        # SparseQR is intentionally omitted from the automated pass: on this
        # imported mesh it did not reach the first Newton record within five
        # minutes.  GMRES-ILUT remains a practical backend comparison.
        for backend in ("gmres_ilut",):
            config, out = prepare_dg_backend_probe(backend)
            backend_runs.append(
                execute(config, out, {"VELA_LINEAR_SOLVER": backend})
            )
        (OUTPUT_ROOT / "dg_linear_backend_execution.json").write_text(
            json.dumps(backend_runs, indent=2) + "\n"
        )
        print(OUTPUT_ROOT)
        return 0

    if args.row_scaling_probe:
        config, out = prepare_dg_row_scaling_probe()
        run = execute(config, out)
        (OUTPUT_ROOT / "dg_row_scaling_execution.json").write_text(
            json.dumps(run, indent=2) + "\n"
        )
        print(OUTPUT_ROOT)
        return run["returncode"]

    if args.qf_reference_ab:
        reference_selected = [
            (model, bias)
            for model in ("dd", "dg")
            for bias in (-1.0, -0.68, -0.52)
        ]
        for model, bias in reference_selected:
            config, out = prepare_single_point(model, bias, qf_reference=True)
            if not args.prepare_only:
                runs.append(execute(config, out))
        if args.prepare_only:
            return 0
        points = [
            summarize_point(model, bias, qf_reference=True)
            for model, bias in reference_selected
        ]
        write_csv(OUTPUT_ROOT / "qf_reference_ab_summary.csv", points)
        (OUTPUT_ROOT / "qf_reference_ab_execution.json").write_text(
            json.dumps({"runs": runs, "points": points}, indent=2) + "\n"
        )
        print(OUTPUT_ROOT)
        return 0

    if args.hard_gate_ab:
        hard_gate_selected = [
            (model, bias)
            for model in ("dd", "dg")
            for bias in (-1.0, -0.68, -0.52)
        ]
        for model, bias in hard_gate_selected:
            config, out = prepare_single_point(model, bias, hard_gate=True)
            if not args.prepare_only:
                runs.append(execute(config, out))
        if args.prepare_only:
            return 0
        (OUTPUT_ROOT / "hard_gate_execution.json").write_text(
            json.dumps({"runs": runs}, indent=2) + "\n"
        )
        print(OUTPUT_ROOT)
        return 0 if all(run["returncode"] == 0 for run in runs) else 1

    if not args.summarize_only:
        for model, bias in selected:
            config, out = prepare_single_point(model, bias)
            if not args.prepare_only:
                runs.append(execute(config, out))
        config, out = prepare_dg_failure_probe()
        if not args.prepare_only:
            runs.append(execute(config, out))
    if args.prepare_only:
        return 0

    points = [summarize_point(model, bias) for model, bias in selected]
    failure = summarize_dg_failure()
    write_csv(OUTPUT_ROOT / "deep_off_precision_summary.csv", points)
    (OUTPUT_ROOT / "dg_residual_summary.json").write_text(
        json.dumps(failure, indent=2) + "\n"
    )
    (OUTPUT_ROOT / "execution_summary.json").write_text(
        json.dumps({"runs": runs, "points": points, "dg_failure": failure}, indent=2)
        + "\n"
    )
    print(OUTPUT_ROOT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
