#!/usr/bin/env python3
"""Recompute all 21-point TransportModels DD/DG curves after the SRH fixes."""

from __future__ import annotations

import csv
import importlib.util
import json
import math
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[1]
BASELINE = REPO / "build-release/reference_tcad/transportmodels_sentaurus2022/vela_baseline"
OUTPUT = BASELINE / "dd_dg_srh_corrected_regression_2026-08-23"
PHASE7_SCRIPT = REPO / "scripts/run_transportmodels_dg_phase7_regression.py"
REFERENCE = BASELINE / "generated/reference_curves"
REPORT_JSON = REPO / "docs/validation/transportmodels_dd_dg_srh_corrected_regression_2026-08-23.json"
REPORT_MD = REPO / "docs/validation/transportmodels_dd_dg_srh_corrected_regression_2026-08-23.md"


def load_phase7(tag: str):
    spec = importlib.util.spec_from_file_location(tag, PHASE7_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load {PHASE7_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def reference_path(branch: str, curve: str) -> Path:
    return REFERENCE / f"transportmodels_sentaurus2022_{branch}_{curve}_reference.csv"


def branch_curves(module, branch: str) -> tuple[dict[str, Any], ...]:
    if branch == "dg":
        idvg_initial = OUTPUT / "dg/idvg/hybrid_restart.csv"
        idvd_initial = (
            BASELINE
            / "dg_discretization_self_consistent_2026-08-21/sentaurus_box"
            / "endpoint_final_state.csv"
        )
    else:
        idvg_initial = (
            BASELINE
            / "dd_phase7_shared_baseline_2026-08-21/idvg/state_bias_2p200000.csv"
        )
        idvd_initial = (
            BASELINE
            / "workflow_dd_vector_run01/dd_idvd_curve_state_bias_2p000000.csv"
        )
    return (
        {
            "name": "idvg",
            "label": f"{branch.upper()} Id-Vg, SRH corrected",
            "config": module.GENERATED / "simulation_dg_idvg.json",
            "reference": reference_path(branch, "idvg"),
            "contact": "gate",
            "current_contact": "drain",
            "points": [2.2 - 0.16 * index for index in range(21)],
            "initial_state": idvg_initial,
        },
        {
            "name": "idvd",
            "label": f"{branch.upper()} Id-Vd, SRH corrected",
            "config": module.GENERATED / "simulation_dg_idvd.json",
            "reference": reference_path(branch, "idvd"),
            "contact": "drain",
            "current_contact": "drain",
            "points": [2.0 - 0.1 * index for index in range(21)],
            "initial_state": idvd_initial,
        },
    )


def add_srh_diagnostics(config: dict[str, Any], run_dir: Path) -> None:
    diagnostics = config["sweep"].setdefault("diagnostics", {})
    diagnostics["srh_balance"] = {
        "enabled": True,
        "material": "Si",
        "drain_contact": "drain",
        "substrate_contact": "substrate",
        "kcl_contacts": ["source", "drain", "gate", "substrate"],
        "resolution_margin_ratio": 10.0,
        "csv_file": str((run_dir / "srh_balance.csv").resolve()),
    }


def configure_module(branch: str):
    module = load_phase7(f"transportmodels_phase7_{branch}_20260823")
    module.OUTPUT_ROOT = OUTPUT / branch
    curves = branch_curves(module, branch)
    module.CURVES = curves
    original_make_config = module.make_config

    def make_config(curve: dict[str, Any]) -> Path:
        path = original_make_config(curve)
        config = json.loads(path.read_text(encoding="utf-8"))
        if branch == "dd":
            config["solver"]["electron_quantum_potential"]["enabled"] = False
        config["_comment"] = (
            f"TransportModels {branch.upper()} corrected 21-point {curve['name']} "
            "with Sentaurus ni, Fermi-BGN, and SRH Nref contract"
        )
        add_srh_diagnostics(config, path.parent)
        path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
        return path

    module.make_config = make_config
    return module, curves


def execute_all(modules: dict[str, Any], curves: dict[str, tuple[dict[str, Any], ...]]) -> list[dict[str, Any]]:
    dg = modules["dg"]
    dg.prepare_idvg_restart()
    dg.run_idvg_frozen_warmup()
    work = [
        (branch, module, curve)
        for branch, module in modules.items()
        for curve in curves[branch]
    ]
    results = []
    with ThreadPoolExecutor(max_workers=2) as pool:
        futures = {
            pool.submit(module.execute_curve, curve): (branch, curve["name"])
            for branch, module, curve in work
        }
        for future in as_completed(futures):
            branch, name = futures[future]
            result = future.result()
            result["branch"] = branch
            result["curve"] = name
            results.append(result)
    return sorted(results, key=lambda row: (row["branch"], row["curve"]))


def branch_metrics(module, branch: str, curves: tuple[dict[str, Any], ...]) -> dict[str, Any]:
    module.OUTPUT_ROOT = OUTPUT / branch
    values = [module.curve_metrics(curve) for curve in curves]
    for value in values:
        candidate = Path(value["candidate_csv"])
        value["candidate_sha256"] = module.sha256(candidate) if candidate.exists() else None
        value["reference_sha256"] = module.sha256(Path(value["reference_csv"]))
    return {
        "branch": branch,
        "complete": all(
            value["all_converged"] and value["aligned_points"] == 21
            for value in values
        ),
        "curves": values,
    }


def plot(report: dict[str, Any]) -> tuple[Path, Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(12.4, 9.2))
    for row_index, branch in enumerate(("dd", "dg")):
        branch_data = next(row for row in report["branches"] if row["branch"] == branch)
        for column_index, curve_name in enumerate(("idvg", "idvd")):
            axis = axes[row_index][column_index]
            curve = next(row for row in branch_data["curves"] if row["name"] == curve_name)
            aligned = curve["aligned"]
            x = [row["bias_V"] for row in aligned]
            sentaurus = [row["sentaurus_A_per_um"] for row in aligned]
            vela = [row["current_A_per_um"] for row in aligned]
            if curve_name == "idvg":
                axis.semilogy(x, sentaurus, "o-", label="Sentaurus 2022")
                axis.semilogy(x, vela, "s--", label="Vela corrected")
                axis.set_xlabel("Gate voltage Vg (V)")
                axis.set_ylabel("Drain current Id (A/um)")
            else:
                axis.plot(x, [value * 1.0e3 for value in sentaurus], "o-", label="Sentaurus 2022")
                axis.plot(x, [value * 1.0e3 for value in vela], "s--", label="Vela corrected")
                axis.set_xlabel("Drain voltage Vd (V)")
                axis.set_ylabel("Drain current Id (mA/um)")
            axis.set_title(f"{branch.upper()} {curve_name.replace('idv', 'Id-V').replace('g', 'g').replace('d', 'd')}")
            axis.grid(True, which="both", alpha=0.25)
            axis.legend()
    fig.suptitle("TransportModels corrected DD/DG 84-point regression")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    png = OUTPUT / "dd_dg_idvg_idvd_comparison.png"
    svg = OUTPUT / "dd_dg_idvg_idvd_comparison.svg"
    fig.savefig(png, dpi=180)
    fig.savefig(svg)
    plt.close(fig)
    return png, svg


def format_metric(value: float | None, percent: bool = False) -> str:
    if value is None:
        return "n/a"
    return f"{value:.4%}" if percent else f"{value:.6g}"


def write_report(report: dict[str, Any]) -> None:
    lines = [
        "# TransportModels corrected DD/DG 84-point regression",
        "",
        f"Status: **{report['status']}**; completed `{report['completed_points']}/84` points.",
        "",
        "Shared corrections: silicon `ni=1.4638914958767616e10 cm^-3`, "
        "Fermi-corrected OldSlotboom BGN, and electron/hole SRH "
        "`Nref=1e16 cm^-3` in Vela internal units.",
        "",
        "| Branch | Curve | Points | Key maximum | Endpoint/secondary |",
        "|---|---|---:|---:|---:|",
    ]
    for branch in report["branches"]:
        for curve in branch["curves"]:
            if curve["metrics"] is None:
                key = secondary = "n/a"
            elif curve["name"] == "idvd":
                key = format_metric(curve["metrics"]["max_relative_error"], True)
                secondary = format_metric(curve["metrics"]["endpoint_relative_error"], True)
            else:
                key = f"off {curve['metrics']['off']['max_absolute_log_error_dex']:.6g} dex"
                secondary = f"on {format_metric(curve['metrics']['on']['max_relative_error'], True)}"
            lines.append(
                f"| {branch['branch'].upper()} | {curve['name']} | "
                f"{curve['completed_points']}/21 | {key} | {secondary} |"
            )
    lines.extend(
        [
            "",
            f"Figure: `{report['artifacts']['png']}`",
            "",
        ]
    )
    REPORT_JSON.parent.mkdir(parents=True, exist_ok=True)
    REPORT_JSON.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    REPORT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    configured = {branch: configure_module(branch) for branch in ("dd", "dg")}
    modules = {branch: value[0] for branch, value in configured.items()}
    curves = {branch: value[1] for branch, value in configured.items()}
    execution = execute_all(modules, curves)
    branches = [branch_metrics(modules[branch], branch, curves[branch]) for branch in ("dd", "dg")]
    completed = sum(
        curve["completed_points"]
        for branch in branches
        for curve in branch["curves"]
    )
    report = {
        "schema": "vela.transportmodels.dd_dg_srh_corrected_regression.v1",
        "as_of": "2026-08-23",
        "status": "complete" if all(branch["complete"] for branch in branches) else "partial",
        "completed_points": completed,
        "expected_points": 84,
        "configuration": {
            "silicon_intrinsic_density_cm3": 1.4638914958767616e10,
            "bandgap_narrowing": {
                "model": "old_slotboom",
                "fermi_statistics_correction": True,
            },
            "srh_reference_doping_cm3": 1.0e16,
            "bias_direction": {
                "idvg": "2.2 V to -1.0 V",
                "idvd": "2.0 V to 0.0 V",
            },
        },
        "execution": execution,
        "branches": branches,
    }
    png, svg = plot(report)
    report["artifacts"] = {"png": str(png.resolve()), "svg": str(svg.resolve())}
    write_report(report)
    print(
        json.dumps(
            {
                "status": report["status"],
                "completed_points": completed,
                "branches": branches,
            },
            indent=2,
        )
    )
    return 0 if report["status"] == "complete" else 1


if __name__ == "__main__":
    raise SystemExit(main())
