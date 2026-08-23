#!/usr/bin/env python3
"""Run the 21-point TransportModels DG curves with the validated quantum contract."""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
from typing import Any, Iterable


REPO = Path(__file__).resolve().parents[1]
REF = REPO / "build-release/reference_tcad/transportmodels_sentaurus2022"
BASELINE = REF / "vela_baseline/dd_dg_srh_corrected_cold_regression_2026-08-23"
GENERATED = BASELINE / "generated_corrected"
OLD_CONFIG = REF / "vela_baseline/vela_fermi_bgn_ab_2026-08-21/dg_on/config.json"
OUTPUT = REF / "vela_baseline/dg_quantum_contract_regression_2026-08-23"
RUN_DIR = OUTPUT / "runs/dg"
WORKFLOW_SCRIPT = REPO / "scripts/run_transportmodels_dd_dg_workflow.py"
DEFAULT_RUNNER = REPO / "build-release/vela_example_runner_quantum_ab.exe"
REPORT_JSON = REPO / "docs/validation/transportmodels_dg_quantum_contract_regression_2026-08-23.json"
REPORT_MD = REPO / "docs/validation/transportmodels_dg_quantum_contract_regression_2026-08-23.md"
PRIOR_REPORT = REPO / "docs/validation/transportmodels_dd_dg_srh_corrected_cold_regression_2026-08-23.json"

IDVG_LOG_LIMIT_DEX = 0.15
IDVG_ON_RELATIVE_LIMIT = 0.10
IDVD_MAX_RELATIVE_LIMIT = 0.05
IDVD_ENDPOINT_RELATIVE_LIMIT = 0.03
DEEP_OFF_RESOLUTION_MARGIN = 10.0


def load_workflow():
    spec = importlib.util.spec_from_file_location("transportmodels_quantum_workflow", WORKFLOW_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load {WORKFLOW_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def percentile(values: Iterable[float], fraction: float) -> float:
    ordered = sorted(float(value) for value in values)
    if not ordered:
        return math.nan
    position = (len(ordered) - 1) * fraction
    low, high = math.floor(position), math.ceil(position)
    if low == high:
        return ordered[low]
    return ordered[low] * (high - position) + ordered[high] * (position - low)


def patch_configs(workflow, manifest: dict[str, Any], quantum: dict[str, Any]) -> None:
    for stage in manifest["stages"]:
        path = Path(stage["config"])
        config = json.loads(path.read_text(encoding="utf-8"))
        config["solver"]["electron_quantum_potential"] = json.loads(json.dumps(quantum))
        config["solver"]["verbose"] = False
        diagnostics = config["sweep"].setdefault("diagnostics", {})
        diagnostics["srh_balance"] = {
            "enabled": True,
            "material": "Si",
            "drain_contact": "drain",
            "substrate_contact": "substrate",
            "kcl_contacts": ["source", "drain", "gate", "substrate"],
            "resolution_margin_ratio": DEEP_OFF_RESOLUTION_MARGIN,
            "csv_file": str((RUN_DIR / f"{stage['name']}_srh_balance.csv").resolve()),
        }
        config["_comment"] = (
            "Corrected material, Fermi-BGN and SRH contract with the independently "
            "validated include-insulators Sentaurus-box DG quantum contract"
        )
        path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
        stage["config_sha256"] = workflow.sha256(path)
    manifest["quantum_contract"] = {
        "source": str(OLD_CONFIG.resolve()),
        "source_sha256": sha256(OLD_CONFIG),
        "payload_sha256": hashlib.sha256(
            json.dumps(quantum, sort_keys=True).encode("utf-8")
        ).hexdigest(),
        "include_insulators": quantum.get("include_insulators"),
        "global_discretization": quantum.get("global_discretization"),
        "interface_boundary": quantum.get("interface_boundary"),
        "formulation": quantum.get("formulation"),
        "coupling_mode": quantum.get("coupling_mode"),
    }
    workflow.write_json(RUN_DIR / "workflow_manifest.json", manifest)


def materialize() -> tuple[Any, dict[str, Any]]:
    workflow = load_workflow()
    old_quantum = json.loads(OLD_CONFIG.read_text(encoding="utf-8"))["solver"][
        "electron_quantum_potential"
    ]
    manifest = workflow.materialize(
        GENERATED,
        RUN_DIR,
        ["dg"],
        quantum_outer_max_iterations=80,
        quantum_outer_acceleration="none",
        quantum_outer_relaxation=1.0,
    )
    patch_configs(workflow, manifest, old_quantum)
    return workflow, manifest


def bias_tag(value: float) -> str:
    prefix = "m" if value < 0.0 else ""
    return prefix + f"{abs(value):.6f}".replace(".", "p")


def update_idvd_prefix() -> tuple[Path, float]:
    """Preserve the completed exact Id-Vd prefix before a resumed stage overwrites its CSV."""
    prefix_path = OUTPUT / "dg_idvd_completed_prefix.csv"
    rows_by_bias: dict[float, dict[str, float]] = {}
    if prefix_path.is_file():
        for row in read_csv(prefix_path):
            bias = round(float(row["bias_V"]), 12)
            rows_by_bias[bias] = {
                "bias_V": bias,
                "current_total_A_per_um": float(row["current_total_A_per_um"]),
            }
    else:
        equilibrium = [
            row for row in read_csv(RUN_DIR / "dg_idvd_equilibrium.csv")
            if row.get("converged") == "1"
        ]
        if not equilibrium:
            raise RuntimeError("Missing converged Id-Vd equilibrium seed")
        rows_by_bias[0.0] = {
            "bias_V": 0.0,
            "current_total_A_per_um": float(equilibrium[-1]["current_total_A_per_um"]),
        }
    reference_biases = set(reference_curve("idvd"))
    curve_path = RUN_DIR / "dg_idvd_curve.csv"
    if curve_path.is_file():
        for row in read_csv(curve_path):
            if row.get("converged") != "1":
                continue
            bias = round(float(row["bias_V"]), 12)
            if bias not in reference_biases:
                continue
            rows_by_bias[bias] = {
                "bias_V": bias,
                "current_total_A_per_um": float(row["current_total_A_per_um"]),
            }
    ordered = [rows_by_bias[bias] for bias in sorted(rows_by_bias)]
    expected = sorted(reference_biases)[: len(ordered)]
    if [row["bias_V"] for row in ordered] != expected:
        raise RuntimeError(f"Id-Vd completed prefix is not contiguous: {ordered}")
    with prefix_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["bias_V", "current_total_A_per_um"])
        writer.writeheader()
        writer.writerows(ordered)
    return prefix_path, ordered[-1]["bias_V"]


def materialize_idvd_resume(bridge_step_V: float) -> tuple[Any, dict[str, Any], float]:
    workflow = load_workflow()
    prefix_path, restart_bias = update_idvd_prefix()
    if restart_bias >= 2.0 - 1.0e-10:
        raise RuntimeError("Id-Vd prefix is already complete")
    restart_state = RUN_DIR / f"dg_idvd_curve_state_bias_{bias_tag(restart_bias)}.csv"
    if not restart_state.is_file():
        raise FileNotFoundError(restart_state)
    old_manifest = RUN_DIR / "workflow_manifest.json"
    archive = RUN_DIR / f"workflow_manifest_failed_through_{bias_tag(restart_bias)}.json"
    if old_manifest.is_file() and not archive.is_file():
        archive.write_text(old_manifest.read_text(encoding="utf-8"), encoding="utf-8")
    if not 0.0 < bridge_step_V < 0.1:
        raise ValueError("bridge_step_V must lie between 0 and 0.1 V")
    bridge_bias = round(restart_bias + bridge_step_V, 12)
    manifest = workflow.materialize(
        GENERATED,
        RUN_DIR,
        ["dg"],
        idvd_curve_restarts={
            "dg": {
                "state": restart_state,
                "prefix": prefix_path,
                "bias_V": restart_bias,
            }
        },
        idvd_bridge_biases={"dg": [bridge_bias]},
        quantum_outer_max_iterations=80,
        quantum_outer_acceleration="none",
        quantum_outer_relaxation=1.0,
    )
    old_quantum = json.loads(OLD_CONFIG.read_text(encoding="utf-8"))["solver"][
        "electron_quantum_potential"
    ]
    patch_configs(workflow, manifest, old_quantum)
    manifest["resume"] = {
        "restart_bias_V": restart_bias,
        "bridge_bias_V": bridge_bias,
        "restart_state": str(restart_state.resolve()),
        "completed_prefix": str(prefix_path.resolve()),
    }
    workflow.write_json(RUN_DIR / "workflow_manifest.json", manifest)
    return workflow, manifest, restart_bias


def reference_curve(name: str) -> dict[float, float]:
    path = REF / "run02/normalized" / f"dg_{name}.csv"
    return {
        round(float(row["bias_V"]), 12): abs(float(row["current_total"]))
        for row in read_csv(path)
    }


def numerical_diagnostics() -> dict[float, dict[str, Any]]:
    result: dict[float, dict[str, Any]] = {}
    for stem in ("dg_idvg_final_bias_relax", "dg_idvg_curve"):
        path = RUN_DIR / f"{stem}_srh_balance.csv"
        if not path.is_file():
            continue
        for row in read_csv(path):
            bias = round(float(row["bias_V"]), 12)
            result[bias] = {
                "numerical_status": row.get("numerical_status", "unknown"),
                "four_terminal_kcl_residual_A_per_um": abs(
                    float(row["four_terminal_kcl_residual_A_per_um"])
                ),
                "id_to_kcl_residual_ratio": float(row["id_to_kcl_residual_ratio"]),
                "silicon_srh_integral_A_per_um": float(
                    row.get("silicon_srh_integral_A_per_um", "nan")
                ),
            }
    return result


def aligned_curve(name: str) -> dict[str, Any]:
    candidate_path = RUN_DIR / f"dg_{name}_curve_comparison_candidate.csv"
    if not candidate_path.is_file():
        return {"name": name, "completed_points": 0, "aligned": [], "metrics": None}
    reference = reference_curve(name)
    diagnostics = numerical_diagnostics() if name == "idvg" else {}
    aligned: list[dict[str, Any]] = []
    for row in read_csv(candidate_path):
        bias = round(float(row["bias_V"]), 12)
        vela = abs(float(row["current_total_A_per_um"]))
        sentaurus = reference[bias]
        aligned.append(
            {
                "bias_V": bias,
                "vela_A_per_um": vela,
                "sentaurus_A_per_um": sentaurus,
                "absolute_relative_error": abs(vela - sentaurus) / max(sentaurus, 1.0e-300),
                "absolute_log_error_dex": abs(
                    math.log10(max(vela, 1.0e-300))
                    - math.log10(max(sentaurus, 1.0e-300))
                ),
                **diagnostics.get(bias, {}),
            }
        )
    aligned.sort(key=lambda item: item["bias_V"])
    metrics: dict[str, Any] | None = None
    if len(aligned) == 21 and name == "idvg":
        regions = {"deep_off": aligned[:3], "transition": aligned[3:8], "on": aligned[8:]}
        metrics = {
            region: {
                "max_relative_error": max(row["absolute_relative_error"] for row in rows),
                "max_absolute_log_error_dex": max(
                    row["absolute_log_error_dex"] for row in rows
                ),
                "median_absolute_log_error_dex": percentile(
                    [row["absolute_log_error_dex"] for row in rows], 0.5
                ),
            }
            for region, rows in regions.items()
        }
        for row in regions["deep_off"]:
            row["resolved_by_kcl"] = (
                row.get("id_to_kcl_residual_ratio", 0.0) >= DEEP_OFF_RESOLUTION_MARGIN
                and row.get("numerical_status") != "numerically_unresolved"
            )
            row["deep_off_acceptance"] = (
                "pass"
                if row["resolved_by_kcl"]
                and row["absolute_log_error_dex"] <= IDVG_LOG_LIMIT_DEX
                else "numerically_unresolved"
                if not row["resolved_by_kcl"]
                else "fail"
            )
    elif len(aligned) == 21:
        nonzero = [row for row in aligned if row["bias_V"] > 0.0]
        metrics = {
            "max_relative_error": max(row["absolute_relative_error"] for row in nonzero),
            "median_relative_error": percentile(
                [row["absolute_relative_error"] for row in nonzero], 0.5
            ),
            "endpoint_relative_error": next(
                row["absolute_relative_error"]
                for row in nonzero
                if math.isclose(row["bias_V"], 2.0)
            ),
        }
    comparison_path = OUTPUT / f"dg_{name}_aligned.csv"
    if aligned:
        fieldnames = list(aligned[0])
        for row in aligned[1:]:
            for key in row:
                if key not in fieldnames:
                    fieldnames.append(key)
        with comparison_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(aligned)
    return {
        "name": name,
        "completed_points": len(aligned),
        "candidate_csv": str(candidate_path.resolve()),
        "aligned_csv": str(comparison_path.resolve()),
        "aligned": aligned,
        "metrics": metrics,
    }


def make_plot(curves: list[dict[str, Any]]) -> tuple[Path, Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    idvg = next(row for row in curves if row["name"] == "idvg")["aligned"]
    idvd = next(row for row in curves if row["name"] == "idvd")["aligned"]
    fig, axes = plt.subplots(2, 2, figsize=(12.2, 8.4))
    axes[0, 0].semilogy(
        [row["bias_V"] for row in idvg],
        [row["sentaurus_A_per_um"] for row in idvg],
        "o-", label="Sentaurus 2022",
    )
    axes[0, 0].semilogy(
        [row["bias_V"] for row in idvg],
        [row["vela_A_per_um"] for row in idvg],
        "s--", label="Vela validated DG contract",
    )
    axes[0, 0].set(xlabel="Gate voltage Vg (V)", ylabel="Drain current Id (A/um)", title="DG Id-Vg")
    axes[0, 0].legend()
    axes[0, 1].plot(
        [row["bias_V"] for row in idvg],
        [row["absolute_log_error_dex"] for row in idvg], "o-",
    )
    axes[0, 1].axhline(IDVG_LOG_LIMIT_DEX, color="tab:red", linestyle="--", label="0.15 dex")
    axes[0, 1].set(xlabel="Gate voltage Vg (V)", ylabel="Absolute log error (dex)", title="Id-Vg log error")
    axes[0, 1].legend()
    axes[1, 0].plot(
        [row["bias_V"] for row in idvd],
        [1.0e3 * row["sentaurus_A_per_um"] for row in idvd],
        "o-", label="Sentaurus 2022",
    )
    axes[1, 0].plot(
        [row["bias_V"] for row in idvd],
        [1.0e3 * row["vela_A_per_um"] for row in idvd],
        "s--", label="Vela validated DG contract",
    )
    axes[1, 0].set(xlabel="Drain voltage Vd (V)", ylabel="Drain current Id (mA/um)", title="DG Id-Vd")
    axes[1, 0].legend()
    nonzero = [row for row in idvd if row["bias_V"] > 0.0]
    axes[1, 1].plot(
        [row["bias_V"] for row in nonzero],
        [100.0 * row["absolute_relative_error"] for row in nonzero], "o-",
    )
    axes[1, 1].axhline(5.0, color="tab:red", linestyle="--", label="5%")
    axes[1, 1].set(xlabel="Drain voltage Vd (V)", ylabel="Absolute relative error (%)", title="Id-Vd relative error")
    axes[1, 1].legend()
    for axis in axes.flat:
        axis.grid(True, which="both", alpha=0.25)
    fig.suptitle("TransportModels DG validated quantum-contract regression")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    png = OUTPUT / "dg_quantum_contract_idvg_idvd_comparison.png"
    svg = OUTPUT / "dg_quantum_contract_idvg_idvd_comparison.svg"
    fig.savefig(png, dpi=180)
    fig.savefig(svg)
    plt.close(fig)
    return png, svg


def acceptance(curves: list[dict[str, Any]]) -> dict[str, Any]:
    idvg = next(row for row in curves if row["name"] == "idvg")
    idvd = next(row for row in curves if row["name"] == "idvd")
    idvg_metrics = idvg["metrics"]
    idvd_metrics = idvd["metrics"]
    main = {
        "idvg_transition": bool(
            idvg_metrics
            and idvg_metrics["transition"]["max_absolute_log_error_dex"] <= IDVG_LOG_LIMIT_DEX
        ),
        "idvg_on": bool(
            idvg_metrics
            and idvg_metrics["on"]["max_relative_error"] <= IDVG_ON_RELATIVE_LIMIT
        ),
        "idvd_full": bool(
            idvd_metrics
            and idvd_metrics["max_relative_error"] <= IDVD_MAX_RELATIVE_LIMIT
            and idvd_metrics["endpoint_relative_error"] <= IDVD_ENDPOINT_RELATIVE_LIMIT
        ),
    }
    deep_rows = idvg["aligned"][:3]
    deep = {
        "policy": "separate log-current plus KCL-resolution branch",
        "points": [
            {
                "bias_V": row["bias_V"],
                "absolute_log_error_dex": row["absolute_log_error_dex"],
                "id_to_kcl_residual_ratio": row.get("id_to_kcl_residual_ratio"),
                "status": row.get("deep_off_acceptance", "missing"),
            }
            for row in deep_rows
        ],
    }
    return {
        "main_curve_pass": len(idvg["aligned"]) == 21 and len(idvd["aligned"]) == 21 and all(main.values()),
        "main_gates": main,
        "deep_off": deep,
    }


def continuation_history(current_manifest: dict[str, Any]) -> list[dict[str, Any]]:
    manifests: list[tuple[Path, dict[str, Any]]] = []
    for path in sorted(RUN_DIR.glob("workflow_manifest_failed_through_*.json")):
        manifests.append((path, json.loads(path.read_text(encoding="utf-8"))))
    manifests.append((RUN_DIR / "workflow_manifest.json", current_manifest))
    rows: list[dict[str, Any]] = []
    seen: set[tuple[float, float]] = set()
    for path, manifest in manifests:
        resume = manifest.get("resume")
        if not resume:
            continue
        key = (float(resume["restart_bias_V"]), float(resume["bridge_bias_V"]))
        if key in seen:
            continue
        seen.add(key)
        bridge_state = RUN_DIR / f"dg_idvd_curve_state_bias_{bias_tag(key[1])}.csv"
        rows.append(
            {
                "restart_bias_V": key[0],
                "bridge_bias_V": key[1],
                "bridge_state_converged": bridge_state.is_file(),
                "bridge_state": str(bridge_state.resolve()),
                "manifest_status": manifest.get("status", "unknown"),
                "manifest": str(path.resolve()),
            }
        )
    rows.sort(key=lambda row: (row["restart_bias_V"], row["bridge_bias_V"]))
    return rows


def prior_baseline_metrics() -> dict[str, Any]:
    report = json.loads(PRIOR_REPORT.read_text(encoding="utf-8"))
    idvg = next(
        row for row in report["curves"]
        if row["branch"] == "dg" and row["curve"] == "idvg"
    )
    idvd = next(
        row for row in report["curves"]
        if row["branch"] == "dg" and row["curve"] == "idvd"
    )
    return {
        "report": str(PRIOR_REPORT.resolve()),
        "idvg_transition_max_log_error_dex": idvg["metrics"]["transition"][
            "max_absolute_log_error_dex"
        ],
        "idvg_on_max_relative_error": idvg["metrics"]["on"]["max_relative_error"],
        "idvd_max_relative_error": idvd["metrics"]["max_relative_error"],
        "idvd_endpoint_relative_error": idvd["metrics"]["endpoint_relative_error"],
    }


def write_report(report: dict[str, Any]) -> None:
    idvg = next(row for row in report["curves"] if row["name"] == "idvg")
    idvd = next(row for row in report["curves"] if row["name"] == "idvd")
    lines = [
        "# TransportModels DG validated quantum-contract regression",
        "",
        f"Execution: **{report['execution_status']}**; main-curve acceptance: "
        f"**{'pass' if report['acceptance']['main_curve_pass'] else 'fail'}**; "
        f"completed `{idvg['completed_points'] + idvd['completed_points']}/42` points.",
        "",
        "The regression uses the corrected material/Fermi-BGN/SRH contract and the "
        "independently validated `include_insulators + sentaurus_box` DG contract. "
        "Deep-off Id-Vg is reported separately and does not veto the normal-region result.",
        "",
        "| Curve/region | Metric | Result | Limit | Status |",
        "|---|---|---:|---:|---|",
    ]
    if idvg["metrics"] and idvd["metrics"]:
        rows = [
            ("Id-Vg transition", "max log error", idvg["metrics"]["transition"]["max_absolute_log_error_dex"], IDVG_LOG_LIMIT_DEX, "dex", report["acceptance"]["main_gates"]["idvg_transition"]),
            ("Id-Vg on", "max relative error", idvg["metrics"]["on"]["max_relative_error"], IDVG_ON_RELATIVE_LIMIT, "%", report["acceptance"]["main_gates"]["idvg_on"]),
            ("Id-Vd", "max relative error", idvd["metrics"]["max_relative_error"], IDVD_MAX_RELATIVE_LIMIT, "%", idvd["metrics"]["max_relative_error"] <= IDVD_MAX_RELATIVE_LIMIT),
            ("Id-Vd 2 V", "endpoint relative error", idvd["metrics"]["endpoint_relative_error"], IDVD_ENDPOINT_RELATIVE_LIMIT, "%", idvd["metrics"]["endpoint_relative_error"] <= IDVD_ENDPOINT_RELATIVE_LIMIT),
        ]
        for region, metric, value, limit, unit, passed in rows:
            if unit == "%":
                value_text, limit_text = f"{value:.3%}", f"{limit:.1%}"
            else:
                value_text, limit_text = f"{value:.6g} dex", f"{limit:.2f} dex"
            lines.append(f"| {region} | {metric} | {value_text} | {limit_text} | {'pass' if passed else 'fail'} |")
        prior = report["prior_baseline"]
        lines.extend([
            "",
            "## Improvement from the corrected cold baseline",
            "",
            "| Metric | Prior DG baseline | Validated quantum contract |",
            "|---|---:|---:|",
            f"| Id-Vg transition max log error | {prior['idvg_transition_max_log_error_dex']:.6g} dex | {idvg['metrics']['transition']['max_absolute_log_error_dex']:.6g} dex |",
            f"| Id-Vg on max relative error | {prior['idvg_on_max_relative_error']:.3%} | {idvg['metrics']['on']['max_relative_error']:.3%} |",
            f"| Id-Vd max relative error | {prior['idvd_max_relative_error']:.3%} | {idvd['metrics']['max_relative_error']:.3%} |",
            f"| Id-Vd 2 V endpoint error | {prior['idvd_endpoint_relative_error']:.3%} | {idvd['metrics']['endpoint_relative_error']:.3%} |",
        ])
    lines.extend([
        "",
        "## Separate deep-off branch",
        "",
        "| Vg (V) | Log error (dex) | Id/KCL residual | Classification |",
        "|---:|---:|---:|---|",
    ])
    for row in report["acceptance"]["deep_off"]["points"]:
        ratio = row["id_to_kcl_residual_ratio"]
        lines.append(
            f"| {row['bias_V']:.2f} | {row['absolute_log_error_dex']:.6g} | "
            f"{ratio:.6g} | {row['status']} |" if ratio is not None else
            f"| {row['bias_V']:.2f} | {row['absolute_log_error_dex']:.6g} | n/a | {row['status']} |"
        )
    lines.extend([
        "",
        "## Id-Vd continuation history",
        "",
        "The extra points below are solver-path bridges only; they are excluded from the "
        "strict 21-point comparison lattice.",
        "",
        "| Restart Vd (V) | Bridge Vd (V) | Bridge converged | Attempt final status |",
        "|---:|---:|---|---|",
    ])
    for row in report["continuation_history"]:
        lines.append(
            f"| {row['restart_bias_V']:.6g} | {row['bridge_bias_V']:.6g} | "
            f"{'yes' if row['bridge_state_converged'] else 'no'} | "
            f"{row['manifest_status']} |"
        )
    lines.extend([
        "",
        f"Figure: `{report['artifacts']['png']}`",
        "",
        f"Run manifest: `{report['manifest']}`",
        "",
    ])
    REPORT_JSON.parent.mkdir(parents=True, exist_ok=True)
    REPORT_JSON.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    REPORT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, default=DEFAULT_RUNNER)
    parser.add_argument("--report-only", action="store_true")
    parser.add_argument(
        "--resume-idvd",
        action="store_true",
        help="resume a failed Id-Vd curve from its highest completed exact point using a 50 mV bridge",
    )
    parser.add_argument(
        "--bridge-step-V",
        type=float,
        default=0.05,
        help="first bridge offset above the completed exact Id-Vd prefix (default: 0.05 V)",
    )
    args = parser.parse_args()
    os.environ["PATH"] = r"D:\msys64\ucrt64\bin" + os.pathsep + os.environ.get("PATH", "")
    OUTPUT.mkdir(parents=True, exist_ok=True)
    previous_path = RUN_DIR / "workflow_manifest.json"
    previous = json.loads(previous_path.read_text(encoding="utf-8")) if previous_path.is_file() else None
    if args.report_only:
        if previous is None:
            raise FileNotFoundError(previous_path)
        manifest = previous
    elif args.resume_idvd:
        workflow, manifest, restart_bias = materialize_idvd_resume(args.bridge_step_V)
        manifest = workflow.execute(manifest, args.runner.resolve(), RUN_DIR, previous)
    else:
        workflow, manifest = materialize()
        manifest = workflow.execute(manifest, args.runner.resolve(), RUN_DIR, previous)
    curves = [aligned_curve("idvg"), aligned_curve("idvd")]
    completed = sum(row["completed_points"] for row in curves)
    if completed:
        png, svg = make_plot(curves)
        artifacts = {"png": str(png.resolve()), "svg": str(svg.resolve())}
    else:
        artifacts = {}
    report = {
        "schema": "vela.transportmodels.dg_quantum_contract_regression.v1",
        "as_of": "2026-08-23",
        "execution_status": "complete" if completed == 42 else "partial",
        "completed_points": completed,
        "expected_points": 42,
        "runner": {"path": str(args.runner.resolve()), "sha256": sha256(args.runner.resolve())},
        "manifest": str((RUN_DIR / "workflow_manifest.json").resolve()),
        "quantum_contract": manifest.get("quantum_contract", {}),
        "curves": curves,
    }
    report["continuation_history"] = continuation_history(manifest)
    report["prior_baseline"] = prior_baseline_metrics()
    report["acceptance"] = acceptance(curves)
    report["artifacts"] = artifacts
    write_report(report)
    print(json.dumps({
        "execution_status": report["execution_status"],
        "completed_points": completed,
        "main_curve_pass": report["acceptance"]["main_curve_pass"],
        "report": str(REPORT_JSON.resolve()),
    }, indent=2))
    return 0 if completed == 42 else 1


if __name__ == "__main__":
    raise SystemExit(main())
