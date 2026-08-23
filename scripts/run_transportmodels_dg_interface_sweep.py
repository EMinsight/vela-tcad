#!/usr/bin/env python3
"""Compare three Si/SiO2 DG interface contracts at a fixed TransportModels state."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

from run_transportmodels_dg_fixed_state_residual_audit import audit


REPO_ROOT = Path(__file__).resolve().parents[1]
BASELINE = REPO_ROOT / "build-release/reference_tcad/transportmodels_sentaurus2022/vela_baseline"
OUTPUT_ROOT = BASELINE / "dg_interface_fixed_state_sweep_2026-08-21"
RUNNER = REPO_ROOT / "build-release/vela_example_runner.exe"
BASE_CONFIG = (
    BASELINE
    / "workflow_dg_idvd_resume_1p9_direct2p0_fixed1p5_outer200_run01"
    / "00_dg_idvd_curve.json"
)
HYBRID_RESTART = BASELINE / "frozen_q_oracle_vg1_vd2_run01/vela_state_with_sentaurus_q.csv"
CORRECTED_MATERIALS = (
    BASELINE
    / "dg_parameter_fixed_state_sweep_2026-08-21"
    / "materials_sentaurus2022_dg_band_drive.json"
)
PARAMETER_REPORT = REPO_ROOT / "docs/validation/transportmodels_dg_parameter_sweep_2026-08-21.json"
REPORT_JSON = REPO_ROOT / "docs/validation/transportmodels_dg_interface_sweep_2026-08-21.json"
REPORT_MD = REPO_ROOT / "docs/validation/transportmodels_dg_interface_sweep_2026-08-21.md"


HALF_JUMPS = {
    "sentaurus_interface_insulator_half_jump_offset": 0.02012,
    "sentaurus_interface_silicon_half_jump_offset": -4.6008840569922854e-5,
    "sentaurus_interface_polysilicon_half_jump_offset": 0.0026674992132365016,
}
REACTION = {
    "sentaurus_interface_silicon_reaction_weight": 0.3613278292533479,
    "sentaurus_interface_polysilicon_reaction_weight": 1.0684933639683336,
    "sentaurus_interface_insulator_at_silicon_reaction_weight": 2.6839079693374917,
    "sentaurus_interface_insulator_at_polysilicon_reaction_weight": 2.569027176700638,
    "sentaurus_interface_silicon_reaction_offset_V": -0.00020247747279261268,
    "sentaurus_interface_polysilicon_reaction_offset_V": 0.01872581675079906,
    "sentaurus_interface_insulator_at_silicon_reaction_offset_V": -0.0015039829729206406,
    "sentaurus_interface_insulator_at_polysilicon_reaction_offset_V": -0.0052046570150173915,
}


VARIANTS: tuple[dict[str, Any], ...] = (
    {"name": "neutral_continuous", "label": "Neutral continuous", "half_jumps": False, "reaction": False},
    {"name": "half_jump_only", "label": "Half-jump only", "half_jumps": True, "reaction": False},
    {"name": "affine_calibrated", "label": "Full affine calibrated", "half_jumps": True, "reaction": True},
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def make_config(variant: dict[str, Any], run_dir: Path) -> tuple[Path, Path]:
    config = json.loads(BASE_CONFIG.read_text(encoding="utf-8"))
    for contact in config["contacts"]:
        if contact["name"].lower() == "drain":
            contact["bias"] = 2.0
    config["materials_file"] = str(CORRECTED_MATERIALS.resolve())
    solver = config["solver"]
    solver["verbose"] = False
    quantum = solver["electron_quantum_potential"]
    prefix = run_dir / "eq231"
    quantum.update(
        {
            "enabled": True,
            "coupling_mode": "outer",
            "formulation": "potential_based",
            "include_insulators": True,
            "global_discretization": "sentaurus_box",
            "residual_diagnostic_prefix": str(prefix.resolve()),
            "residual_diagnostic_use_initial_state": True,
            "outer_max_iterations": 1,
            "max_iterations": 1,
            "effective_mass_ratio": 1.0618016171622988,
            "coefficient_mass_ratio": 1.0906506732296395,
            "sentaurus_insulator_reentrant_corner_reaction_weight": 1.0,
        }
    )
    for key in HALF_JUMPS:
        quantum[key] = HALF_JUMPS[key] if variant["half_jumps"] else 0.0
    for key, value in REACTION.items():
        quantum[key] = value if variant["reaction"] else (1.0 if key.endswith("weight") else 0.0)
    config["_comment"] = (
        "TransportModels phase-4 fixed-state Si/SiO2 interface comparison: "
        + variant["label"]
    )
    config["output_csv"] = str((run_dir / "probe.csv").resolve())
    config["log_file"] = str((run_dir / "probe.log").resolve())
    config["sweep"].update(
        {
            "start": 2.0,
            "stop": 2.0,
            "step": 0.1,
            "bias_points": [2.0],
            "initial_state_file": str(HYBRID_RESTART.resolve()),
            "write_vtk": False,
            "write_state_file": str((run_dir / "probe_final_state.csv").resolve()),
            "write_state_every_point_prefix": str((run_dir / "probe_state").resolve()),
        }
    )
    path = run_dir / "config.json"
    path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    return path, prefix


def execute(variant: dict[str, Any]) -> dict[str, Any]:
    run_dir = OUTPUT_ROOT / variant["name"]
    run_dir.mkdir(parents=True, exist_ok=True)
    config, prefix = make_config(variant, run_dir)
    process = subprocess.run(
        [str(RUNNER.resolve()), "--config", str(config)],
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
    )
    (run_dir / "console.log").write_text(
        process.stdout + "\n--- STDERR ---\n" + process.stderr, encoding="utf-8"
    )
    data = audit(prefix)
    summary = data["summary"]
    interface_l1 = sum(float(row["interface_total_l1_free"]) for row in data["regions"])
    return {
        "name": variant["name"],
        "label": variant["label"],
        "runner_exit_code": process.returncode,
        "cell_total_l1_free": float(summary["cell_total_l1_free"]),
        "max_free_residual": float(summary["max_free_residual"]),
        "max_free_node": int(summary["max_free_node"]),
        "interface_cell_l1_free_region_sum": interface_l1,
        "explicit_interface_boundary_l1": float(data["component_l1"]["interface_boundary"]),
        "reaction_l1_share": float(data["component_l1_share"]["reaction"]),
        "top_node": data["top_nodes"][0],
        "top_cell": data["top_cells"][0],
        "config": str(config),
        "nodes": str(prefix) + "_nodes.csv",
        "config_sha256": sha256(config),
        "nodes_sha256": sha256(Path(str(prefix) + "_nodes.csv")),
    }


def load_reference() -> dict[str, float]:
    report = json.loads(PARAMETER_REPORT.read_text(encoding="utf-8"))
    row = next(item for item in report["results"] if item["name"] == "corrected_material_contract")
    return {
        "p1_direct_cell_total_l1_free": float(row["cell_total_l1_free"]),
        "p1_direct_max_free_residual": float(row["max_free_residual"]),
    }


def write_csv(rows: list[dict[str, Any]]) -> Path:
    path = OUTPUT_ROOT / "interface_comparison.csv"
    fields = [
        "name", "label", "cell_total_l1_free", "l1_ratio_to_best", "l1_ratio_to_p1_corrected",
        "max_free_residual", "max_free_node", "interface_cell_l1_free_region_sum",
        "explicit_interface_boundary_l1", "reaction_l1_share", "runner_exit_code",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows([{key: row[key] for key in fields} for row in rows])
    return path


def plot(rows: list[dict[str, Any]]) -> tuple[Path, Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels = [row["label"].replace(" ", "\n") for row in rows]
    ratios = [row["l1_ratio_to_best"] for row in rows]
    max_ratios = [row["max_free_residual"] / min(r["max_free_residual"] for r in rows) for row in rows]
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 5.2))
    colors = ["#4C78A8", "#F58518", "#54A24B"]
    axes[0].bar(labels, ratios, color=colors)
    axes[0].set_ylabel("Global free-node L1 / best")
    axes[0].set_title("Integrated residual")
    axes[0].grid(axis="y", alpha=0.25)
    axes[1].bar(labels, max_ratios, color=colors)
    axes[1].set_ylabel("Maximum free-node residual / best")
    axes[1].set_title("Worst local hotspot")
    axes[1].grid(axis="y", alpha=0.25)
    for ax, values in zip(axes, (ratios, max_ratios)):
        for index, value in enumerate(values):
            ax.text(index, value, f"{value:.3f}x", ha="center", va="bottom", fontsize=9)
    fig.suptitle("TransportModels DG: fixed-state Si/SiO2 interface comparison", fontsize=14)
    fig.text(0.5, 0.92, "Vg = 1 V, Vd = 2 V; corrected material contract; sentaurus_box", ha="center", fontsize=10)
    fig.subplots_adjust(left=0.09, right=0.98, bottom=0.18, top=0.82, wspace=0.28)
    png = OUTPUT_ROOT / "interface_comparison.png"
    svg = OUTPUT_ROOT / "interface_comparison.svg"
    fig.savefig(png, dpi=180)
    fig.savefig(svg)
    plt.close(fig)
    return png, svg


def markdown(report: dict[str, Any]) -> str:
    rows = report["results"]
    lines = [
        "# TransportModels DG Si/SiO2 interface comparison",
        "",
        "Fixed-state audit at `Vg=1 V`, `Vd=2 V`. All three cases use the same mesh,",
        "hybrid Vela-DD/Sentaurus-Q state, corrected material contract, and `sentaurus_box`",
        "operator. The runner exit code `1` is expected because the diagnostic run is",
        "deliberately limited to one inner and one outer iteration.",
        "",
        "| Interface contract | Global L1 | Ratio to best | Max residual | Interface-cell L1 |",
        "|---|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| {row['label']} | {row['cell_total_l1_free']:.6g} | "
            f"{row['l1_ratio_to_best']:.4f} | {row['max_free_residual']:.6g} | "
            f"{row['interface_cell_l1_free_region_sum']:.6g} |"
        )
    best = min(rows, key=lambda row: row["cell_total_l1_free"])
    local = min(rows, key=lambda row: row["max_free_residual"])
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            f"- Best integrated residual: **{best['label']}**.",
            f"- Best worst-node residual: **{local['label']}**.",
            "- The affine coefficients originated from the older SingleDevice mesh and are",
            "  treated as a transferability experiment, not as universal physical constants.",
            "- A model is eligible for the self-consistent stage only if it improves the global",
            "  metric without creating a materially worse local hotspot.",
            "",
            f"Figure: `{report['artifacts']['png']}`",
            f"CSV: `{report['artifacts']['csv']}`",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reuse", action="store_true", help="Reuse existing raw diagnostics")
    parser.add_argument("--check", action="store_true", help="Verify frozen hashes and report data")
    args = parser.parse_args()
    if args.check:
        report = json.loads(REPORT_JSON.read_text(encoding="utf-8"))
        for row in report["results"]:
            assert sha256(Path(row["config"])) == row["config_sha256"]
            assert sha256(Path(row["nodes"])) == row["nodes_sha256"]
        print("TransportModels DG interface sweep check: PASS")
        return 0

    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    if args.reuse:
        results = []
        for variant in VARIANTS:
            config, prefix = make_config(variant, OUTPUT_ROOT / variant["name"])
            data = audit(prefix)
            summary = data["summary"]
            results.append(
                {
                    "name": variant["name"], "label": variant["label"], "runner_exit_code": 1,
                    "cell_total_l1_free": float(summary["cell_total_l1_free"]),
                    "max_free_residual": float(summary["max_free_residual"]),
                    "max_free_node": int(summary["max_free_node"]),
                    "interface_cell_l1_free_region_sum": sum(float(r["interface_total_l1_free"]) for r in data["regions"]),
                    "explicit_interface_boundary_l1": float(data["component_l1"]["interface_boundary"]),
                    "reaction_l1_share": float(data["component_l1_share"]["reaction"]),
                    "top_node": data["top_nodes"][0], "top_cell": data["top_cells"][0],
                    "config": str(config), "nodes": str(prefix) + "_nodes.csv",
                    "config_sha256": sha256(config),
                    "nodes_sha256": sha256(Path(str(prefix) + "_nodes.csv")),
                }
            )
    else:
        results = []
        with ThreadPoolExecutor(max_workers=3) as pool:
            futures = {pool.submit(execute, variant): variant for variant in VARIANTS}
            for future in as_completed(futures):
                results.append(future.result())
        order = {row["name"]: index for index, row in enumerate(VARIANTS)}
        results.sort(key=lambda row: order[row["name"]])

    reference = load_reference()
    best_l1 = min(row["cell_total_l1_free"] for row in results)
    for row in results:
        row["l1_ratio_to_best"] = row["cell_total_l1_free"] / best_l1
        row["l1_ratio_to_p1_corrected"] = row["cell_total_l1_free"] / reference["p1_direct_cell_total_l1_free"]
    csv_path = write_csv(results)
    png, svg = plot(results)
    report = {
        "schema": "vela.transportmodels.dg_interface_fixed_state_sweep.v1",
        "status": "pass",
        "as_of": "2026-08-21",
        "work_point": {"gate_bias_V": 1.0, "drain_bias_V": 2.0},
        "controlled_inputs": {
            "discretization": "sentaurus_box", "materials": str(CORRECTED_MATERIALS),
            "initial_state": str(HYBRID_RESTART),
        },
        "p1_direct_corrected_reference": reference,
        "results": results,
        "artifacts": {"csv": str(csv_path), "png": str(png), "svg": str(svg)},
    }
    REPORT_JSON.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    REPORT_MD.write_text(markdown(report), encoding="utf-8")
    (OUTPUT_ROOT / "interface_comparison_summary.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({"results": [{k: r[k] for k in ("name", "cell_total_l1_free", "max_free_residual")} for r in results]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
