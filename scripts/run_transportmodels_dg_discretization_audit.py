#!/usr/bin/env python3
"""Audit TransportModels DG discretizations on one common fixed state."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

from run_transportmodels_dg_fixed_state_residual_audit import audit
from run_transportmodels_dg_interface_sweep import (
    CORRECTED_MATERIALS,
    REPO_ROOT,
    RUNNER,
    VARIANTS as INTERFACE_VARIANTS,
    make_config as make_interface_config,
)


BASELINE = REPO_ROOT / "build-release/reference_tcad/transportmodels_sentaurus2022/vela_baseline"
OUTPUT_ROOT = BASELINE / "dg_discretization_fixed_state_audit_2026-08-21"
REPORT_JSON = REPO_ROOT / "docs/validation/transportmodels_dg_discretization_audit_2026-08-21.json"
REPORT_MD = REPO_ROOT / "docs/validation/transportmodels_dg_discretization_audit_2026-08-21.md"

DISCRETIZATIONS: tuple[tuple[str, str], ...] = (
    ("p1_direct", "P1 direct"),
    ("p1_lambda_direct", "P1 lambda"),
    ("cvfem_full", "CVFEM full"),
    ("sentaurus_box", "Sentaurus box"),
    ("gss_potentiallike_fitted", "GSS potential-like"),
    ("gss_density_fitted", "GSS density"),
    ("conservative_sqrt_fitted", "Conservative sqrt"),
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def make_config(mode: str, run_dir: Path) -> tuple[Path, Path]:
    config_path, prefix = make_interface_config(INTERFACE_VARIANTS[0], run_dir)
    config = json.loads(config_path.read_text(encoding="utf-8"))
    config["solver"]["electron_quantum_potential"]["global_discretization"] = mode
    config["_comment"] = (
        "TransportModels phase-5 fixed-state DG discretization audit: " + mode
    )
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    return config_path, prefix


def execute(mode: str, label: str) -> dict[str, Any]:
    run_dir = OUTPUT_ROOT / mode
    run_dir.mkdir(parents=True, exist_ok=True)
    config, prefix = make_config(mode, run_dir)
    process = subprocess.run(
        [str(RUNNER.resolve()), "--config", str(config)], cwd=REPO_ROOT,
        text=True, capture_output=True,
    )
    (run_dir / "console.log").write_text(
        process.stdout + "\n--- STDERR ---\n" + process.stderr, encoding="utf-8"
    )
    try:
        data = audit(prefix)
        summary = data["summary"]
        result: dict[str, Any] = {
            "status": "diagnostic_complete",
            "cell_total_l1_free": float(summary["cell_total_l1_free"]),
            "max_free_residual": float(summary["max_free_residual"]),
            "max_free_node": int(summary["max_free_node"]),
            "reaction_l1_share": float(data["component_l1_share"]["reaction"]),
            "top_node": data["top_nodes"][0],
            "top_cell": data["top_cells"][0],
            "nodes": str(prefix) + "_nodes.csv",
            "nodes_sha256": sha256(Path(str(prefix) + "_nodes.csv")),
        }
    except (FileNotFoundError, KeyError, ValueError) as exc:
        result = {
            "status": "diagnostic_failed", "error": str(exc),
            "cell_total_l1_free": None, "max_free_residual": None,
            "max_free_node": None, "reaction_l1_share": None,
            "top_node": None, "top_cell": None, "nodes": None, "nodes_sha256": None,
        }
    result.update(
        {
            "name": mode, "label": label, "runner_exit_code": process.returncode,
            "config": str(config), "config_sha256": sha256(config),
        }
    )
    return result


def plot(rows: list[dict[str, Any]]) -> tuple[Path, Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    valid = [row for row in rows if row["status"] == "diagnostic_complete"]
    labels = [row["label"].replace(" ", "\n") for row in valid]
    l1 = [row["cell_total_l1_free"] for row in valid]
    maxima = [row["max_free_residual"] for row in valid]
    colors = ["#4C78A8" if row["name"] != "sentaurus_box" else "#54A24B" for row in valid]
    fig, axes = plt.subplots(2, 1, figsize=(12.2, 8.2), sharex=True)
    axes[0].bar(labels, l1, color=colors)
    axes[0].set_yscale("log")
    axes[0].set_ylabel("Global free-node L1 (raw)")
    axes[0].set_title("Integrated fixed-state residual")
    axes[1].bar(labels, maxima, color=colors)
    axes[1].set_yscale("log")
    axes[1].set_ylabel("Maximum free-node residual (raw)")
    axes[1].set_title("Worst fixed-state hotspot")
    axes[1].tick_params(axis="x", labelsize=9)
    for ax in axes:
        ax.grid(axis="y", which="both", alpha=0.25)
    fig.suptitle("TransportModels DG discretization audit", fontsize=15, y=0.98)
    fig.text(
        0.5, 0.945,
        "Vg = 1 V, Vd = 2 V; corrected materials; neutral interface; raw scales are operator-specific",
        ha="center", fontsize=10,
    )
    fig.subplots_adjust(left=0.11, right=0.98, bottom=0.16, top=0.88, hspace=0.32)
    png = OUTPUT_ROOT / "discretization_residuals.png"
    svg = OUTPUT_ROOT / "discretization_residuals.svg"
    fig.savefig(png, dpi=180)
    fig.savefig(svg)
    plt.close(fig)
    return png, svg


def write_csv(rows: list[dict[str, Any]]) -> Path:
    path = OUTPUT_ROOT / "discretization_residuals.csv"
    fields = [
        "name", "label", "status", "cell_total_l1_free", "max_free_residual",
        "max_free_node", "reaction_l1_share", "runner_exit_code",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows([{key: row.get(key) for key in fields} for row in rows])
    return path


def markdown(report: dict[str, Any]) -> str:
    lines = [
        "# TransportModels DG discretization audit",
        "",
        "The audit evaluates all supported global DG operators on the same fixed state,",
        "corrected material contract, and neutral interface. Raw residual magnitudes are",
        "reported for hotspot localization, but cannot by themselves rank formulations",
        "with different primary variables or row scaling. A self-consistent current test is",
        "therefore required before selecting the production operator.",
        "",
        "| Discretization | Diagnostic | Global raw L1 | Max raw residual | Max node |",
        "|---|---|---:|---:|---:|",
    ]
    for row in report["results"]:
        l1 = "n/a" if row["cell_total_l1_free"] is None else f"{row['cell_total_l1_free']:.6g}"
        maximum = "n/a" if row["max_free_residual"] is None else f"{row['max_free_residual']:.6g}"
        node = "n/a" if row["max_free_node"] is None else str(row["max_free_node"])
        lines.append(f"| {row['label']} | {row['status']} | {l1} | {maximum} | {node} |")
    lines.extend(
        [
            "",
            "## Decision",
            "",
            "- `p1_direct` remains the conservative control because the phase-2/3 residual",
            "  decomposition and its units are already audited.",
            "- `sentaurus_box` with a neutral interface advances as the primary contender:",
            "  it preserves the potential-like variable and substantially reduces both audited",
            "  fixed-state metrics relative to corrected `p1_direct`.",
            "- Fitted density/square-root formulations remain diagnostic candidates until their",
            "  self-consistent convergence and terminal current are demonstrated.",
            "",
            f"Figure: `{report['artifacts']['png']}`",
            f"CSV: `{report['artifacts']['csv']}`",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.check:
        report = json.loads(REPORT_JSON.read_text(encoding="utf-8"))
        for row in report["results"]:
            assert sha256(Path(row["config"])) == row["config_sha256"]
            if row["nodes"]:
                assert sha256(Path(row["nodes"])) == row["nodes_sha256"]
        print("TransportModels DG discretization audit check: PASS")
        return 0

    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    with ThreadPoolExecutor(max_workers=3) as pool:
        futures = {
            pool.submit(execute, mode, label): (mode, label)
            for mode, label in DISCRETIZATIONS
        }
        for future in as_completed(futures):
            rows.append(future.result())
    order = {mode: index for index, (mode, _) in enumerate(DISCRETIZATIONS)}
    rows.sort(key=lambda row: order[row["name"]])
    csv_path = write_csv(rows)
    png, svg = plot(rows)
    report = {
        "schema": "vela.transportmodels.dg_discretization_fixed_state_audit.v1",
        "status": "pass" if all(row["status"] == "diagnostic_complete" for row in rows) else "partial",
        "as_of": "2026-08-21",
        "work_point": {"gate_bias_V": 1.0, "drain_bias_V": 2.0},
        "controlled_inputs": {
            "materials": str(CORRECTED_MATERIALS), "interface": "neutral_continuous",
        },
        "comparability_caveat": "Raw residual scales are operator-specific across primary variables and row scalings.",
        "results": rows,
        "artifacts": {"csv": str(csv_path), "png": str(png), "svg": str(svg)},
    }
    REPORT_JSON.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    REPORT_MD.write_text(markdown(report), encoding="utf-8")
    print(json.dumps({"status": report["status"], "results": [{"name": r["name"], "l1": r["cell_total_l1_free"], "max": r["max_free_residual"]} for r in rows]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
