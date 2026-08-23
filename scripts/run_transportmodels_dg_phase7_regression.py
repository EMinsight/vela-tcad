#!/usr/bin/env python3
"""Run and report the 42-point corrected TransportModels DG regression."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[1]
BASELINE = REPO_ROOT / "build-release/reference_tcad/transportmodels_sentaurus2022/vela_baseline"
OUTPUT_ROOT = BASELINE / "dg_phase7_regression_2026-08-21"
RUNNER = REPO_ROOT / "build-release/vela_example_runner.exe"
GENERATED = BASELINE / "generated/vela"
REFERENCE = BASELINE / "generated/reference_curves"
CORRECTED_MATERIALS = BASELINE / "dg_parameter_fixed_state_sweep_2026-08-21/materials_sentaurus2022_dg_band_drive.json"
REPORT_JSON = REPO_ROOT / "docs/validation/transportmodels_dg_phase7_regression_2026-08-21.json"
REPORT_MD = REPO_ROOT / "docs/validation/transportmodels_dg_phase7_regression_2026-08-21.md"
PLOT_SUBTITLE = (
    "Corrected materials + sentaurus_box + neutral interface + existing Lombardi mobility"
)
REPORT_CONFIGURATION_EXTRA: dict[str, Any] = {}

CURVES: tuple[dict[str, Any], ...] = (
    {
        "name": "idvg", "label": "DG Id-Vg", "config": GENERATED / "simulation_dg_idvg.json",
        "reference": REFERENCE / "transportmodels_sentaurus2022_dg_idvg_reference.csv",
        "contact": "gate", "current_contact": "drain",
        "points": [2.2 - 0.16 * index for index in range(21)],
        "initial_state": OUTPUT_ROOT / "idvg/hybrid_restart.csv",
    },
    {
        "name": "idvd", "label": "DG Id-Vd", "config": GENERATED / "simulation_dg_idvd.json",
        "reference": REFERENCE / "transportmodels_sentaurus2022_dg_idvd_reference.csv",
        "contact": "drain", "current_contact": "drain",
        "points": [2.0 - 0.1 * index for index in range(21)],
        "initial_state": BASELINE / "dg_discretization_self_consistent_2026-08-21/sentaurus_box/endpoint_final_state.csv",
    },
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def runner_environment() -> dict[str, str]:
    environment = os.environ.copy()
    ucrt_bin = r"D:\msys64\ucrt64\bin"
    environment["PATH"] = ucrt_bin + os.pathsep + environment.get("PATH", "")
    return environment


def percentile(values: list[float], p: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * p
    low, high = int(math.floor(position)), int(math.ceil(position))
    if low == high:
        return ordered[low]
    return ordered[low] * (high - position) + ordered[high] * (position - low)


def make_config(curve: dict[str, Any]) -> Path:
    run_dir = OUTPUT_ROOT / curve["name"]
    run_dir.mkdir(parents=True, exist_ok=True)
    config = json.loads(Path(curve["config"]).read_text(encoding="utf-8"))
    source_dir = Path(curve["config"]).parent
    for key in ("mesh_file", "node_doping_file"):
        source_path = Path(config[key])
        if not source_path.is_absolute():
            config[key] = str((source_dir / source_path).resolve())
    config["materials_file"] = str(CORRECTED_MATERIALS.resolve())
    for contact in config["contacts"]:
        if contact["name"].lower() == curve["contact"]:
            contact["bias"] = curve["points"][0]
    config["solver"]["mobility"]["high_field_gradient_discretization"] = "transport_cell_vector"
    config["solver"]["bandgap_narrowing"] = {
        "model": "old_slotboom",
        "fermi_statistics_correction": True,
    }
    srh_doping = config["solver"].get("srh_doping_dependence", {})
    for carrier in ("electron", "hole"):
        if carrier in srh_doping:
            srh_doping[carrier]["reference_doping_m3"] = 1.0e16
    config["solver"]["quasi_fermi_update_limit_V"] = 0.025
    quantum = config["solver"]["electron_quantum_potential"]
    quantum.update(
        {
            "enabled": True, "coupling_mode": "outer", "formulation": "potential_based",
            "include_insulators": True, "global_discretization": "sentaurus_box",
            "effective_mass_ratio": 1.0618016171622988,
            "coefficient_mass_ratio": 1.0906506732296395,
            "outer_max_iterations": 80, "max_iterations": 60,
            "outer_absolute_tolerance_V": 0.01, "damping": 0.5,
            "outer_acceleration": "none", "outer_relaxation": 1.0,
            "outer_relaxation_min": 0.1, "outer_relaxation_max": 1.0,
            "sentaurus_interface_insulator_half_jump_offset": 0.0,
            "sentaurus_interface_silicon_half_jump_offset": 0.0,
            "sentaurus_interface_polysilicon_half_jump_offset": 0.0,
            "sentaurus_interface_silicon_reaction_weight": 1.0,
            "sentaurus_interface_polysilicon_reaction_weight": 1.0,
            "sentaurus_interface_insulator_at_silicon_reaction_weight": 1.0,
            "sentaurus_interface_insulator_at_polysilicon_reaction_weight": 1.0,
            "sentaurus_interface_silicon_reaction_offset_V": 0.0,
            "sentaurus_interface_polysilicon_reaction_offset_V": 0.0,
            "sentaurus_interface_insulator_at_silicon_reaction_offset_V": 0.0,
            "sentaurus_interface_insulator_at_polysilicon_reaction_offset_V": 0.0,
            "sentaurus_insulator_reentrant_corner_reaction_weight": 1.0,
        }
    )
    config["_comment"] = "TransportModels phase-7 corrected 21-point " + curve["label"]
    config["output_csv"] = str((run_dir / curve.get("output_name", "curve.csv")).resolve())
    config["log_file"] = str((run_dir / "curve.log").resolve())
    sweep = config["sweep"]
    sweep.update(
        {
            "mode": "iv", "contact": curve["contact"], "current_contact": curve["current_contact"],
            "start": curve["points"][0], "stop": curve["points"][-1],
            "step": (curve["points"][1] - curve["points"][0]) if len(curve["points"]) > 1 else 0.1,
            "bias_points": curve["points"],
            "initial_state_file": str(Path(curve["initial_state"]).resolve()),
            "write_vtk": False, "write_state_file": str((run_dir / "final_state.csv").resolve()),
            "write_state_every_point_prefix": str((run_dir / "state").resolve()),
        }
    )
    diagnostics = sweep.setdefault("diagnostics", {})
    diagnostics["transport"] = {"enabled": True}
    diagnostics["terminal_balance"] = {
        "enabled": True, "contacts": ["source", "drain", "gate", "substrate"],
        "csv_file": str((run_dir / "terminal_balance.csv").resolve()),
    }
    path = run_dir / "config.json"
    path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    return path


def execute_curve(curve: dict[str, Any]) -> dict[str, Any]:
    run_dir = OUTPUT_ROOT / curve["name"]
    old_log = run_dir / "config.log"
    combined_path = run_dir / "curve_combined.csv"
    prior_rows = load_candidate(combined_path) if combined_path.exists() else []
    prior_rows.extend(parse_terminal_balance(run_dir / "terminal_balance.csv"))
    if old_log.exists():
        prior_rows.extend(parse_solve_traces(old_log))
    if old_log.exists() and old_log.stat().st_size:
        (run_dir / "resume_previous_config.log").write_text(old_log.read_text(encoding="utf-8"), encoding="utf-8")
    completed = 0
    for bias in curve["points"]:
        candidate = run_dir / f"state_bias_{bias_slug(bias)}.csv"
        if not candidate.exists() or not candidate.stat().st_size:
            break
        completed += 1
    known_biases = {round(float(row["bias_V"]), 12) for row in prior_rows}
    known_count = 0
    for bias in curve["points"]:
        if round(float(bias), 12) not in known_biases:
            break
        known_count += 1
    active = dict(curve)
    if known_count == len(curve["points"]):
        write_combined_curve(run_dir / "curve_combined.csv", prior_rows, [])
        return {"name": curve["name"], "runner_exit_code": 0, "resumed_after_points": known_count, "config": None, "console": None}
    if known_count:
        # Re-solve the most recently saved point. Its state is durable, but a
        # Ctrl-C can leave the buffered curve/terminal row incomplete.
        resume_index = max(0, known_count - 1)
        active["points"] = curve["points"][resume_index:]
        active["initial_state"] = run_dir / f"state_bias_{bias_slug(curve['points'][resume_index])}.csv"
        active["output_name"] = "curve_resume.csv"
    if not active["points"]:
        write_combined_curve(run_dir / "curve_combined.csv", prior_rows, [])
        return {"name": curve["name"], "runner_exit_code": 0, "resumed_after_points": known_count, "config": None, "console": None}
    config = make_config(active)
    process = subprocess.run(
        [str(RUNNER), "--config", str(config)], cwd=REPO_ROOT,
        text=True, capture_output=True, env=runner_environment(),
    )
    console = run_dir / "console.log"
    console.write_text(process.stdout + "\n--- STDERR ---\n" + process.stderr, encoding="utf-8")
    new_rows = load_candidate(run_dir / active.get("output_name", "curve.csv"))
    write_combined_curve(run_dir / "curve_combined.csv", prior_rows, new_rows)
    return {"name": curve["name"], "runner_exit_code": process.returncode, "resumed_after_points": known_count, "config": str(config), "console": str(console)}


def bias_slug(value: float) -> str:
    sign = "m" if value < -5.0e-13 else ""
    return sign + f"{abs(value):.6f}".replace(".", "p")


def parse_solve_traces(path: Path) -> list[dict[str, Any]]:
    pattern = re.compile(
        r"solve_trace: index=\d+ bias_V=(?P<bias>[-+0-9.eE]+) converged=(?P<converged>[01]) "
        r"iterations=(?P<iterations>\d+) current_total=(?P<current>[-+0-9.eE]+)"
    )
    rows = []
    for line in path.read_text(encoding="utf-8").splitlines():
        match = pattern.search(line)
        if match:
            rows.append(
                {
                    "bias_V": float(match.group("bias")),
                    "current_A_per_um": abs(float(match.group("current"))) * 1.0e-6,
                    "converged": match.group("converged") == "1",
                    "iterations": int(match.group("iterations")),
                }
            )
    return rows


def parse_terminal_balance(path: Path) -> list[dict[str, Any]]:
    if not path.exists() or not path.stat().st_size:
        return []
    rows = []
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row.get("contact") != "drain" or row.get("converged") != "1":
                continue
            if not row.get("current_total_A_per_um") or not row.get("newton_iterations"):
                continue
            rows.append(
                {
                    "bias_V": float(row["bias_V"]),
                    "current_A_per_um": abs(float(row["current_total_A_per_um"])),
                    "converged": True,
                    "iterations": int(row["newton_iterations"]),
                }
            )
    return rows


def write_combined_curve(path: Path, prior: list[dict[str, Any]], current: list[dict[str, Any]]) -> None:
    by_bias = {round(float(row["bias_V"]), 12): row for row in prior}
    by_bias.update({round(float(row["bias_V"]), 12): row for row in current})
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["bias_V", "current_total_A_per_um", "converged", "iterations"])
        writer.writeheader()
        for bias in sorted(by_bias):
            row = by_bias[bias]
            writer.writerow(
                {
                    "bias_V": row["bias_V"], "current_total_A_per_um": row["current_A_per_um"],
                    "converged": "1" if row["converged"] else "0", "iterations": row["iterations"],
                }
            )


def prepare_idvg_restart() -> Path:
    sentaurus = OUTPUT_ROOT / "idvg/sentaurus_restart.csv"
    hybrid = OUTPUT_ROOT / "idvg/hybrid_restart.csv"
    sentaurus.parent.mkdir(parents=True, exist_ok=True)
    process = subprocess.run(
        [
            sys.executable, str(REPO_ROOT / "scripts/sentaurus_fields_to_restart.py"),
            "--export-dir", str(BASELINE / "generated/sim_fields/dg_idvg"),
            "--output", str(sentaurus), "--preserve-insulator-quantum-potential",
        ],
        cwd=REPO_ROOT, text=True, capture_output=True,
    )
    if process.returncode != 0:
        raise RuntimeError("Id-Vg Sentaurus restart generation failed: " + process.stderr)
    with sentaurus.open(newline="", encoding="utf-8") as handle:
        sent_rows = {int(row["node_id"]): row for row in csv.DictReader(handle)}
    vela_state = BASELINE / "workflow_dg_outer80_resume_m036_run01/dg_idvg_curve_state_bias_2p200000.csv"
    with vela_state.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fieldnames = list(reader.fieldnames or [])
        rows = list(reader)
    if "electron_quantum_potential_like_V" not in fieldnames:
        fieldnames.append("electron_quantum_potential_like_V")
    with hybrid.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            sent = sent_rows[int(row["node_id"])]
            row["electron_quantum_potential_V"] = sent["electron_quantum_potential_V"]
            row["electron_quantum_potential_like_V"] = sent["electron_quantum_potential_like_V"]
            writer.writerow(row)
    return hybrid


def run_idvg_frozen_warmup() -> Path:
    curve = CURVES[0]
    run_dir = OUTPUT_ROOT / "idvg"
    existing = run_dir / "warmup_final_state.csv"
    if existing.exists() and existing.stat().st_size:
        curve["initial_state"] = existing
        return existing
    base_path = make_config(curve)
    config = json.loads(base_path.read_text(encoding="utf-8"))
    config["solver"]["electron_quantum_potential"]["coupling_mode"] = "frozen"
    config["_comment"] = "TransportModels phase-7 Id-Vg Frozen-Q warmup at Vg=2.2 V"
    config["output_csv"] = str((run_dir / "warmup.csv").resolve())
    config["sweep"].update(
        {
            "start": 2.2, "stop": 2.2, "step": 0.16, "bias_points": [2.2],
            "initial_state_file": str(Path(curve["initial_state"]).resolve()),
            "write_state_file": str((run_dir / "warmup_final_state.csv").resolve()),
            "write_state_every_point_prefix": str((run_dir / "warmup_state").resolve()),
        }
    )
    path = run_dir / "warmup_config.json"
    path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    process = subprocess.run(
        [str(RUNNER), "--config", str(path)], cwd=REPO_ROOT,
        text=True, capture_output=True, env=runner_environment(),
    )
    (run_dir / "warmup.console.log").write_text(process.stdout + "\n--- STDERR ---\n" + process.stderr, encoding="utf-8")
    output = run_dir / "warmup_final_state.csv"
    if process.returncode != 0 or not output.exists():
        raise RuntimeError("Id-Vg Frozen-Q warmup failed: " + process.stderr)
    curve["initial_state"] = output
    return output


def load_reference(path: Path) -> dict[float, float]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return {round(float(row["bias_V"]), 12): abs(float(row["current_total"])) for row in csv.DictReader(handle)}


def load_candidate(path: Path) -> list[dict[str, Any]]:
    if not path.exists() or not path.stat().st_size:
        return []
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    return [
        {
            "bias_V": float(row["bias_V"]), "current_A_per_um": abs(float(row["current_total_A_per_um"])),
            "converged": row["converged"] == "1", "iterations": int(row["iterations"]),
        }
        for row in rows
    ]


def curve_metrics(curve: dict[str, Any]) -> dict[str, Any]:
    combined = OUTPUT_ROOT / curve["name"] / "curve_combined.csv"
    candidate_path = combined if combined.exists() else OUTPUT_ROOT / curve["name"] / "curve.csv"
    candidate = load_candidate(candidate_path)
    reference = load_reference(Path(curve["reference"]))
    aligned = []
    for row in candidate:
        bias = round(row["bias_V"], 12)
        if bias not in reference:
            continue
        sent = reference[bias]
        vela = row["current_A_per_um"]
        aligned.append(
            {
                **row, "sentaurus_A_per_um": sent,
                "absolute_relative_error": abs(vela - sent) / sent if sent > 1.0e-16 else None,
                "absolute_log_error_dex": abs(math.log10(max(vela, 1.0e-30)) - math.log10(max(sent, 1.0e-30))),
            }
        )
    aligned.sort(key=lambda row: row["bias_V"])
    result: dict[str, Any] = {
        "name": curve["name"], "expected_points": 21, "completed_points": len(candidate),
        "aligned_points": len(aligned), "all_converged": len(candidate) == 21 and all(row["converged"] for row in candidate),
        "candidate_csv": str(candidate_path), "reference_csv": str(curve["reference"]), "aligned": aligned,
    }
    if curve["name"] == "idvd":
        valid = [row for row in aligned if row["bias_V"] > 0 and row["absolute_relative_error"] is not None]
        result["metrics"] = {
            "max_relative_error": max(row["absolute_relative_error"] for row in valid),
            "median_relative_error": percentile([row["absolute_relative_error"] for row in valid], 0.5),
            "endpoint_relative_error": next(row["absolute_relative_error"] for row in valid if abs(row["bias_V"] - 2.0) < 1.0e-9),
        } if len(valid) == 20 else None
    else:
        regimes = {"off": aligned[:3], "transition": aligned[3:8], "on": aligned[8:]}
        result["metrics"] = {
            key: {
                "points": len(rows),
                "max_absolute_log_error_dex": max(row["absolute_log_error_dex"] for row in rows),
                "median_absolute_log_error_dex": percentile([row["absolute_log_error_dex"] for row in rows], 0.5),
                "max_relative_error": max((row["absolute_relative_error"] for row in rows if row["absolute_relative_error"] is not None), default=None),
            }
            for key, rows in regimes.items()
        } if len(aligned) == 21 else None
    return result


def surface_node_ids() -> set[int]:
    mesh = json.loads((GENERATED / "mesh.json").read_text(encoding="utf-8"))
    edge_regions: dict[tuple[int, int], set[int]] = {}
    for tri in mesh["triangles"]:
        ids = [int(value) for value in tri["node_ids"]]
        for pair in (tuple(sorted((ids[0], ids[1]))), tuple(sorted((ids[1], ids[2]))), tuple(sorted((ids[2], ids[0])))):
            edge_regions.setdefault(pair, set()).add(int(tri["region_id"]))
    return {node for pair, regions in edge_regions.items() if regions == {0, 3} for node in pair}


def read_scalar(path: Path) -> dict[int, float]:
    with path.open(newline="", encoding="utf-8") as handle:
        return {int(row["node_id"]): float(row["component0"]) for row in csv.DictReader(handle)}


def spatial_metrics() -> dict[str, Any] | None:
    # The phase-7 curve is intentionally qualified at a looser 10 mV outer
    # tolerance for runtime.  Use the phase-5 strict 1 mV Vd=2 V endpoint for
    # the field comparison; ``idvd/final_state.csv`` is the last reverse-sweep
    # state (Vd=0 V), so comparing it with the Sentaurus Vd=2 V fields is invalid.
    state_path = BASELINE / "dg_discretization_self_consistent_2026-08-21/sentaurus_box/endpoint_final_state.csv"
    if not state_path.exists() or not state_path.stat().st_size:
        return None
    with state_path.open(newline="", encoding="utf-8") as handle:
        vela = {int(row["node_id"]): row for row in csv.DictReader(handle)}
    fields = BASELINE / "generated/sim_fields/dg_idvd/fields"
    sent_q = read_scalar(fields / "eQuantumPotential_region3.csv")
    sent_n = read_scalar(fields / "eDensity_region3.csv")
    nodes = sorted(surface_node_ids() & sent_q.keys() & sent_n.keys() & vela.keys())
    q_errors = [abs(float(vela[node]["electron_quantum_potential_V"]) - sent_q[node]) * 1.0e3 for node in nodes]
    n_errors = [abs(math.log10(max(float(vela[node]["electrons_m3"]) / 1.0e6, 1.0)) - math.log10(max(sent_n[node], 1.0))) for node in nodes]
    return {
        "vela_state_csv": str(state_path.resolve()),
        "bias": {"gate_V": 1.0, "drain_V": 2.0},
        "outer_absolute_tolerance_V": 0.001,
        "surface_node_count": len(nodes),
        "quantum_potential_absolute_error_mV": {"median": percentile(q_errors, 0.5), "p95": percentile(q_errors, 0.95), "p99": percentile(q_errors, 0.99), "max": max(q_errors)},
        "electron_density_absolute_log_error_dex": {"median": percentile(n_errors, 0.5), "p95": percentile(n_errors, 0.95), "p99": percentile(n_errors, 0.99), "max": max(n_errors)},
    }


def plot(curves: list[dict[str, Any]]) -> tuple[Path, Path]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    idvg = next(curve for curve in curves if curve["name"] == "idvg")["aligned"]
    idvd = next(curve for curve in curves if curve["name"] == "idvd")["aligned"]
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 5.3))
    axes[0].semilogy([row["bias_V"] for row in idvg], [row["sentaurus_A_per_um"] for row in idvg], "o-", label="Sentaurus 2022")
    axes[0].semilogy([row["bias_V"] for row in idvg], [row["current_A_per_um"] for row in idvg], "s--", label="Vela corrected")
    axes[0].set_xlabel("Gate voltage Vg (V)")
    axes[0].set_ylabel("Drain current Id (A/µm)")
    axes[0].set_title("DG Id-Vg")
    axes[1].plot([row["bias_V"] for row in idvd], [row["sentaurus_A_per_um"] * 1e3 for row in idvd], "o-", label="Sentaurus 2022")
    axes[1].plot([row["bias_V"] for row in idvd], [row["current_A_per_um"] * 1e3 for row in idvd], "s--", label="Vela corrected")
    axes[1].set_xlabel("Drain voltage Vd (V)")
    axes[1].set_ylabel("Drain current Id (mA/µm)")
    axes[1].set_title("DG Id-Vd")
    for ax in axes:
        ax.grid(True, which="both", alpha=0.25)
        ax.legend()
    fig.suptitle(
        "TransportModels corrected DG regression\n"
        + PLOT_SUBTITLE,
        fontsize=13,
        y=0.99,
        linespacing=1.45,
    )
    fig.subplots_adjust(left=0.09, right=0.98, bottom=0.13, top=0.81, wspace=0.25)
    png = OUTPUT_ROOT / "dg_idvg_idvd_comparison.png"
    svg = OUTPUT_ROOT / "dg_idvg_idvd_comparison.svg"
    fig.savefig(png, dpi=180)
    fig.savefig(svg)
    plt.close(fig)
    return png, svg


def markdown(report: dict[str, Any]) -> str:
    idvg = next(row for row in report["curves"] if row["name"] == "idvg")
    idvd = next(row for row in report["curves"] if row["name"] == "idvd")
    lines = [
        "# TransportModels corrected DG 42-point regression",
        "",
        f"Execution status: **{report['status']}**. Acceptance status: **{report['acceptance_status']}**. Completed: `{idvg['completed_points']}/21` Id-Vg and `{idvd['completed_points']}/21` Id-Vd points.",
        "",
    ]
    if report["status"] == "complete":
        lines.extend(
            [
                "| Acceptance metric | Result | Limit | Pass |",
                "|---|---:|---:|---:|",
                f"| DG Id-Vd max relative error | {idvd['metrics']['max_relative_error']:.4%} | 5% | {idvd['metrics']['max_relative_error'] <= 0.05} |",
                f"| DG Id-Vd endpoint relative error | {idvd['metrics']['endpoint_relative_error']:.4%} | 3% | {idvd['metrics']['endpoint_relative_error'] <= 0.03} |",
                f"| DG Id-Vg on max relative error | {idvg['metrics']['on']['max_relative_error']:.4%} | 10% | {idvg['metrics']['on']['max_relative_error'] <= 0.10} |",
                f"| DG Id-Vg transition max log error | {idvg['metrics']['transition']['max_absolute_log_error_dex']:.4f} dex | 0.15 dex | {idvg['metrics']['transition']['max_absolute_log_error_dex'] <= 0.15} |",
                "",
            ]
        )
    if report.get("spatial"):
        spatial = report["spatial"]
        lines.extend(
            [
                "## Endpoint surface fields",
                "",
                f"- Quantum-potential p95 absolute error: `{spatial['quantum_potential_absolute_error_mV']['p95']:.6g} mV`.",
                f"- Electron-density p95 absolute log error: `{spatial['electron_density_absolute_log_error_dex']['p95']:.6g} dex`.",
                "",
            ]
        )
    lines.extend([f"Figure: `{report['artifacts']['png']}`", ""])
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--report-only", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.check:
        report = json.loads(REPORT_JSON.read_text(encoding="utf-8"))
        for row in report["curves"]:
            assert sha256(Path(row["candidate_csv"])) == row["candidate_sha256"]
            assert sha256(Path(row["reference_csv"])) == row["reference_sha256"]
        print("TransportModels DG phase-7 regression check: PASS")
        return 0
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    execution = []
    if not args.report_only:
        prepare_idvg_restart()
        run_idvg_frozen_warmup()
        with ThreadPoolExecutor(max_workers=2) as pool:
            futures = {pool.submit(execute_curve, curve): curve["name"] for curve in CURVES}
            for future in as_completed(futures):
                execution.append(future.result())
    metrics = [curve_metrics(curve) for curve in CURVES]
    complete = all(row["all_converged"] and row["aligned_points"] == 21 for row in metrics)
    accepted = complete and (
        next(row for row in metrics if row["name"] == "idvd")["metrics"]["max_relative_error"] <= 0.05
        and next(row for row in metrics if row["name"] == "idvd")["metrics"]["endpoint_relative_error"] <= 0.03
        and next(row for row in metrics if row["name"] == "idvg")["metrics"]["on"]["max_relative_error"] <= 0.10
        and next(row for row in metrics if row["name"] == "idvg")["metrics"]["transition"]["max_absolute_log_error_dex"] <= 0.15
    )
    for row in metrics:
        candidate = Path(row["candidate_csv"])
        row["candidate_sha256"] = sha256(candidate) if candidate.exists() else None
        row["reference_sha256"] = sha256(Path(row["reference_csv"]))
    png, svg = plot(metrics)
    report = {
        "schema": "vela.transportmodels.dg_phase7_regression.v1", "as_of": "2026-08-21",
        "status": "complete" if complete else "partial",
        "acceptance_status": "pass" if accepted else "fail",
        "execution": execution,
        "configuration": {
            "materials": str(CORRECTED_MATERIALS), "discretization": "sentaurus_box",
            "interface": "neutral_continuous", "mobility": "masetti_field_lombardi",
            "outer_absolute_tolerance_V": 0.01,
            **REPORT_CONFIGURATION_EXTRA,
        },
        "curves": metrics, "spatial": spatial_metrics(),
        "artifacts": {"png": str(png), "svg": str(svg)},
    }
    REPORT_JSON.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    REPORT_MD.write_text(markdown(report), encoding="utf-8")
    print(json.dumps({"status": report["status"], "curves": [{"name": row["name"], "completed": row["completed_points"], "all_converged": row["all_converged"], "metrics": row["metrics"]} for row in metrics], "spatial": report["spatial"]}, indent=2))
    return 0 if complete else 1


if __name__ == "__main__":
    raise SystemExit(main())
