#!/usr/bin/env python3
"""Cold-start recomputation of all corrected TransportModels DD/DG curves."""

from __future__ import annotations

import csv
import importlib.util
import json
import math
import os
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[1]
BASELINE = REPO / "build-release/reference_tcad/transportmodels_sentaurus2022/vela_baseline"
SOURCE_GENERATED = BASELINE / "generated"
OUTPUT = BASELINE / "dd_dg_srh_corrected_cold_regression_2026-08-23"
CORRECTED_GENERATED = OUTPUT / "generated_corrected"
RUN_ROOT = OUTPUT / "runs"
RUNNER = REPO / "build-release/vela_example_runner.exe"
MATERIALS = (
    BASELINE
    / "dg_parameter_fixed_state_sweep_2026-08-21/materials_sentaurus2022_dg_band_drive.json"
)
WORKFLOW_SCRIPT = REPO / "scripts/run_transportmodels_dd_dg_workflow.py"
REPORT_JSON = REPO / "docs/validation/transportmodels_dd_dg_srh_corrected_cold_regression_2026-08-23.json"
REPORT_MD = REPO / "docs/validation/transportmodels_dd_dg_srh_corrected_cold_regression_2026-08-23.md"


def load_workflow():
    spec = importlib.util.spec_from_file_location("transportmodels_corrected_workflow", WORKFLOW_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load {WORKFLOW_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def prepare_generated() -> None:
    vela_output = CORRECTED_GENERATED / "vela"
    reference_output = CORRECTED_GENERATED / "reference_curves"
    vela_output.mkdir(parents=True, exist_ok=True)
    reference_output.mkdir(parents=True, exist_ok=True)
    source_vela = SOURCE_GENERATED / "vela"
    for branch in ("dd", "dg"):
        for curve in ("idvg", "idvd"):
            name = f"simulation_{branch}_{curve}.json"
            config = json.loads((source_vela / name).read_text(encoding="utf-8"))
            for key in ("mesh_file", "node_doping_file"):
                value = Path(config[key])
                config[key] = str((source_vela / value).resolve()) if not value.is_absolute() else str(value)
            config["materials_file"] = str(MATERIALS.resolve())
            config["solver"]["bandgap_narrowing"] = {
                "model": "old_slotboom",
                "fermi_statistics_correction": True,
            }
            srh = config["solver"]["srh_doping_dependence"]
            srh["electron"]["reference_doping_m3"] = 1.0e16
            srh["hole"]["reference_doping_m3"] = 1.0e16
            (vela_output / name).write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
            reference_name = f"transportmodels_sentaurus2022_{branch}_{curve}_reference.csv"
            shutil.copyfile(SOURCE_GENERATED / "reference_curves" / reference_name, reference_output / reference_name)


def patch_stage_configs(workflow, manifest: dict[str, Any], run_dir: Path) -> None:
    for stage in manifest["stages"]:
        path = Path(stage["config"])
        config = json.loads(path.read_text(encoding="utf-8"))
        quantum = config["solver"].get("electron_quantum_potential", {})
        if quantum.get("enabled", False):
            quantum["outer_absolute_tolerance_V"] = 0.01
        diagnostics = config["sweep"].setdefault("diagnostics", {})
        diagnostics["srh_balance"] = {
            "enabled": True,
            "material": "Si",
            "drain_contact": "drain",
            "substrate_contact": "substrate",
            "kcl_contacts": ["source", "drain", "gate", "substrate"],
            "resolution_margin_ratio": 10.0,
            "csv_file": str((run_dir / f"{stage['name']}_srh_balance.csv").resolve()),
        }
        path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
        stage["config_sha256"] = workflow.sha256(path)
    workflow.write_json(run_dir / "workflow_manifest.json", manifest)


def run_branch(branch: str) -> dict[str, Any]:
    workflow = load_workflow()
    run_dir = RUN_ROOT / branch
    previous_path = run_dir / "workflow_manifest.json"
    previous = (
        json.loads(previous_path.read_text(encoding="utf-8"))
        if previous_path.exists()
        else None
    )
    manifest = workflow.materialize(
        CORRECTED_GENERATED,
        run_dir,
        [branch],
        quantum_outer_max_iterations=80,
        quantum_outer_acceleration="aitken",
        quantum_outer_relaxation=1.0,
    )
    patch_stage_configs(workflow, manifest, run_dir)
    result = workflow.execute(manifest, RUNNER, run_dir, previous)
    return {
        "branch": branch,
        "status": result["status"],
        "manifest": str((run_dir / "workflow_manifest.json").resolve()),
        "stages": [
            {
                "name": stage["name"],
                "status": stage.get("status", "pending"),
                "execution": stage.get("execution"),
            }
            for stage in result["stages"]
        ],
    }


def read_curve(path: Path) -> list[dict[str, float]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return [
            {
                "bias_V": float(row["bias_V"]),
                "current_A_per_um": abs(float(row["current_total_A_per_um"])),
            }
            for row in csv.DictReader(handle)
        ]


def read_reference(path: Path) -> dict[float, float]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return {
            round(float(row["bias_V"]), 12): abs(float(row["current_total"]))
            for row in csv.DictReader(handle)
        }


def percentile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * fraction
    low, high = math.floor(position), math.ceil(position)
    if low == high:
        return ordered[low]
    return ordered[low] * (high - position) + ordered[high] * (position - low)


def curve_metrics(branch: str, curve: str) -> dict[str, Any]:
    candidate = RUN_ROOT / branch / f"{branch}_{curve}_curve_comparison_candidate.csv"
    reference_path = (
        CORRECTED_GENERATED
        / "reference_curves"
        / f"transportmodels_sentaurus2022_{branch}_{curve}_reference.csv"
    )
    if not candidate.exists():
        return {"branch": branch, "curve": curve, "completed_points": 0, "metrics": None}
    reference = read_reference(reference_path)
    numerical_by_bias: dict[float, dict[str, Any]] = {}
    if curve == "idvg":
        for path in (
            RUN_ROOT / branch / f"{branch}_idvg_final_bias_relax_srh_balance.csv",
            RUN_ROOT / branch / f"{branch}_idvg_curve_srh_balance.csv",
        ):
            with path.open(newline="", encoding="utf-8-sig") as handle:
                for row in csv.DictReader(handle):
                    numerical_by_bias[round(float(row["bias_V"]), 12)] = {
                        "numerical_status": row["numerical_status"],
                        "four_terminal_kcl_residual_A_per_um": float(
                            row["four_terminal_kcl_residual_A_per_um"]
                        ),
                        "id_to_kcl_residual_ratio": float(row["id_to_kcl_residual_ratio"]),
                    }
    aligned = []
    for row in read_curve(candidate):
        bias = round(row["bias_V"], 12)
        sentaurus = reference[bias]
        vela = row["current_A_per_um"]
        aligned.append(
            {
                **row,
                "sentaurus_A_per_um": sentaurus,
                "absolute_relative_error": abs(vela - sentaurus) / sentaurus,
                "absolute_log_error_dex": abs(
                    math.log10(max(vela, 1.0e-30))
                    - math.log10(max(sentaurus, 1.0e-30))
                ),
                **numerical_by_bias.get(bias, {}),
            }
        )
    aligned.sort(key=lambda row: row["bias_V"])
    if curve == "idvg":
        regimes = {"off": aligned[:3], "transition": aligned[3:8], "on": aligned[8:]}
        metrics = {
            name: {
                "max_relative_error": max(row["absolute_relative_error"] for row in rows),
                "max_absolute_log_error_dex": max(row["absolute_log_error_dex"] for row in rows),
                "median_absolute_log_error_dex": percentile(
                    [row["absolute_log_error_dex"] for row in rows], 0.5
                ),
            }
            for name, rows in regimes.items()
        }
    else:
        nonzero = [row for row in aligned if row["bias_V"] > 0.0]
        metrics = {
            "max_relative_error": max(row["absolute_relative_error"] for row in nonzero),
            "median_relative_error": percentile(
                [row["absolute_relative_error"] for row in nonzero], 0.5
            ),
            "endpoint_relative_error": next(
                row["absolute_relative_error"]
                for row in nonzero
                if abs(row["bias_V"] - 2.0) < 1.0e-10
            ),
        }
    return {
        "branch": branch,
        "curve": curve,
        "completed_points": len(aligned),
        "candidate_csv": str(candidate.resolve()),
        "reference_csv": str(reference_path.resolve()),
        "aligned": aligned,
        "unresolved_biases_V": [
            row["bias_V"]
            for row in aligned
            if row.get("numerical_status") == "numerically_unresolved"
        ],
        "metrics": metrics if len(aligned) == 21 else None,
    }


def make_plot(curves: list[dict[str, Any]]) -> tuple[Path, Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(12.4, 9.2))
    for row_index, branch in enumerate(("dd", "dg")):
        for column_index, curve_name in enumerate(("idvg", "idvd")):
            result = next(
                row for row in curves
                if row["branch"] == branch and row["curve"] == curve_name
            )
            axis = axes[row_index][column_index]
            aligned = result.get("aligned", [])
            x = [row["bias_V"] for row in aligned]
            sentaurus = [row["sentaurus_A_per_um"] for row in aligned]
            vela = [row["current_A_per_um"] for row in aligned]
            if curve_name == "idvg":
                axis.semilogy(x, sentaurus, "o-", label="Sentaurus 2022")
                axis.semilogy(x, vela, "s--", label="Vela corrected")
                axis.set_xlabel("Gate voltage Vg (V)")
                axis.set_ylabel("Drain current Id (A/um)")
            else:
                axis.plot(x, [value * 1e3 for value in sentaurus], "o-", label="Sentaurus 2022")
                axis.plot(x, [value * 1e3 for value in vela], "s--", label="Vela corrected")
                axis.set_xlabel("Drain voltage Vd (V)")
                axis.set_ylabel("Drain current Id (mA/um)")
            axis.set_title(f"{branch.upper()} {'Id-Vg' if curve_name == 'idvg' else 'Id-Vd'}")
            axis.grid(True, which="both", alpha=0.25)
            axis.legend()
    fig.suptitle("TransportModels corrected cold-start DD/DG regression")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    png = OUTPUT / "dd_dg_idvg_idvd_comparison.png"
    svg = OUTPUT / "dd_dg_idvg_idvd_comparison.svg"
    fig.savefig(png, dpi=180)
    fig.savefig(svg)
    plt.close(fig)
    return png, svg


def write_report(report: dict[str, Any]) -> None:
    lines = [
        "# TransportModels corrected cold-start DD/DG regression",
        "",
        f"Execution status: **{report['status']}**; acceptance status: "
        f"**{report['acceptance_status']}**; completed "
        f"`{report['completed_points']}/84` comparison points.",
        "",
        "Sweep semantics match the Sentaurus decks: Id-Vg initializes at "
        "`Vg=-1 V`, ramps the drain to `1.1 V`, then sweeps to `2.2 V`; "
        "Id-Vd initializes separately at `Vg=1 V`, `Vd=0 V`, then sweeps to `2 V`.",
        "",
        "| Branch | Curve | Points | Primary metric | Secondary metric | Tertiary metric | Numerical qualification |",
        "|---|---|---:|---:|---:|---:|---|",
    ]
    for curve in report["curves"]:
        metrics = curve["metrics"]
        if metrics is None:
            primary = secondary = tertiary = "n/a"
        elif curve["curve"] == "idvg":
            primary = f"off {metrics['off']['max_absolute_log_error_dex']:.6g} dex"
            secondary = (
                f"transition {metrics['transition']['max_absolute_log_error_dex']:.6g} dex"
            )
            tertiary = f"on {metrics['on']['max_relative_error']:.4%}"
        else:
            primary = f"max {metrics['max_relative_error']:.4%}"
            secondary = f"median {metrics['median_relative_error']:.4%}"
            tertiary = f"endpoint {metrics['endpoint_relative_error']:.4%}"
        lines.append(
            f"| {curve['branch'].upper()} | {curve['curve']} | "
            f"{curve['completed_points']}/21 | {primary} | {secondary} | {tertiary} | "
            + (
                f"unresolved at {curve['unresolved_biases_V']} |"
                if curve["unresolved_biases_V"]
                else "resolved |"
            )
        )
    lines.extend(["", f"Figure: `{report['artifacts']['png']}`", ""])
    REPORT_JSON.parent.mkdir(parents=True, exist_ok=True)
    REPORT_JSON.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    REPORT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    os.environ["PATH"] = r"D:\msys64\ucrt64\bin" + os.pathsep + os.environ.get("PATH", "")
    OUTPUT.mkdir(parents=True, exist_ok=True)
    prepare_generated()
    # The two sparse DD/DG runners can both terminate with Windows status
    # 0xFFFFFFFF when their factorization peaks overlap on this workstation.
    # Run branch dependency chains sequentially; stage-level resume preserves
    # all previously completed equilibrium/ramp work.
    with ThreadPoolExecutor(max_workers=1) as pool:
        futures = {pool.submit(run_branch, branch): branch for branch in ("dd", "dg")}
        execution = [future.result() for future in as_completed(futures)]
    execution.sort(key=lambda row: row["branch"])
    curves = [curve_metrics(branch, curve) for branch in ("dd", "dg") for curve in ("idvg", "idvd")]
    completed = sum(curve["completed_points"] for curve in curves)
    curve_acceptance = []
    for curve in curves:
        metrics = curve["metrics"]
        accepted = metrics is not None and not curve["unresolved_biases_V"]
        if metrics is not None and curve["curve"] == "idvg":
            accepted = accepted and (
                metrics["transition"]["max_absolute_log_error_dex"] <= 0.15
                and metrics["on"]["max_relative_error"] <= 0.10
            )
        elif metrics is not None:
            accepted = accepted and (
                metrics["max_relative_error"] <= 0.05
                and metrics["endpoint_relative_error"] <= 0.03
            )
        curve_acceptance.append(
            {"branch": curve["branch"], "curve": curve["curve"], "accepted": accepted}
        )
    report = {
        "schema": "vela.transportmodels.dd_dg_srh_corrected_cold_regression.v1",
        "as_of": "2026-08-23",
        "status": "complete" if completed == 84 else "partial",
        "acceptance_status": (
            "pass" if completed == 84 and all(row["accepted"] for row in curve_acceptance)
            else "fail"
        ),
        "completed_points": completed,
        "expected_points": 84,
        "configuration": {
            "silicon_intrinsic_density_cm3": 1.4638914958767616e10,
            "bandgap_narrowing": {"model": "old_slotboom", "fermi_statistics_correction": True},
            "srh_reference_doping_cm3": 1.0e16,
            "sweep_direction": "Sentaurus-matched low-to-high",
        },
        "execution": execution,
        "curve_acceptance": curve_acceptance,
        "curves": curves,
    }
    png, svg = make_plot(curves)
    report["artifacts"] = {"png": str(png.resolve()), "svg": str(svg.resolve())}
    write_report(report)
    print(json.dumps({"status": report["status"], "completed_points": completed, "execution": execution}, indent=2))
    return 0 if completed == 84 else 1


if __name__ == "__main__":
    raise SystemExit(main())
