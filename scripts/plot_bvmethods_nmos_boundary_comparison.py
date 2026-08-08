#!/usr/bin/env python3
"""Generate BVmethods NMOS field, I-V, and runtime comparison artifacts."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from datetime import datetime
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
from matplotlib.colors import TwoSlopeNorm
from matplotlib.ticker import LogLocator


REPO = Path(__file__).resolve().parents[1]
RUN_ROOT = REPO / "build-release/reference_tcad/bvmethods_sentaurus2018/run01"
DEFAULT_SENT = RUN_ROOT / "sentaurus_boundary_state_20260808/imported/current_1e4"
DEFAULT_SENT_PLT = (
    RUN_ROOT
    / "sentaurus_boundary_state_20260808/raw/boundary_current_1e4_des.plt"
)
DEFAULT_SENT_LOG = (
    RUN_ROOT
    / "sentaurus_boundary_state_20260808/raw/boundary_current_1e4_des.log"
)
DEFAULT_VELA_STATE = (
    RUN_ROOT
    / "vela_validation/boundary_voltage_to_current_20260806"
    / "boundary_control_checkpoints/current_target_0p000100_eval_8.csv"
)
DEFAULT_MESH = RUN_ROOT / "vela/mesh.json"
DEFAULT_PREBIAS = (
    RUN_ROOT
    / "vela_validation/boundary_external_resistor_20260806/prebias/sweep.csv"
)
DEFAULT_CURRENT_SWEEP = (
    RUN_ROOT
    / "vela_validation/boundary_voltage_to_current_20260806/sweep.csv"
)
DEFAULT_CURRENT_LOG = (
    RUN_ROOT
    / "vela_validation/boundary_voltage_to_current_20260806/simulation.log"
)
DEFAULT_RESISTOR_LOG = (
    RUN_ROOT
    / "vela_validation/boundary_external_resistor_20260806/simulation.log"
)
DEFAULT_OUT = RUN_ROOT / "sentaurus_boundary_state_20260808/report_20260808"

BLUE = "#2F5D8C"
ORANGE = "#D9772B"
GOLD = "#B58A2A"
INK = "#25313C"
GRID = "#D9DEE3"


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_rows(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0])
    known = set(fieldnames)
    for row in rows[1:]:
        for field in row:
            if field not in known:
                fieldnames.append(field)
                known.add(field)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def sent_scalar(root: Path, name: str) -> dict[int, float]:
    return {
        int(row["node_id"]): float(row["component0"])
        for row in read_rows(root / "fields" / f"{name}_region3.csv")
    }


def sent_vector(root: Path, name: str) -> dict[int, tuple[float, float]]:
    return {
        int(row["node_id"]): (float(row["component0"]), float(row["component1"]))
        for row in read_rows(root / "fields" / f"{name}_region3.csv")
    }


def load_mesh(mesh_path: Path, semiconductor_nodes: list[int]) -> tuple[
    np.ndarray, np.ndarray, np.ndarray, dict[int, int]
]:
    mesh = json.loads(mesh_path.read_text(encoding="utf-8"))
    coordinates = {
        int(node["id"]): (float(node["x"]), float(node["y"]))
        for node in mesh["nodes"]
    }
    local = {node: index for index, node in enumerate(semiconductor_nodes)}
    triangles = []
    for triangle in mesh["triangles"]:
        if int(triangle["region_id"]) != 3:
            continue
        ids = [int(value) for value in triangle["node_ids"]]
        if all(node in local for node in ids):
            triangles.append([local[node] for node in ids])
    x = np.asarray([coordinates[node][0] for node in semiconductor_nodes])
    y = np.asarray([coordinates[node][1] for node in semiconductor_nodes])
    return x, y, np.asarray(triangles, dtype=int), local


def reconstruct_node_electric_field(
    x_um: np.ndarray,
    y_um: np.ndarray,
    triangles: np.ndarray,
    potential: np.ndarray,
) -> np.ndarray:
    accumulated = np.zeros((len(x_um), 2), dtype=float)
    weights = np.zeros(len(x_um), dtype=float)
    x_m = x_um * 1.0e-6
    y_m = y_um * 1.0e-6
    for triangle in triangles:
        a, b, c = (int(value) for value in triangle)
        matrix = np.asarray([
            [x_m[b] - x_m[a], y_m[b] - y_m[a]],
            [x_m[c] - x_m[a], y_m[c] - y_m[a]],
        ])
        determinant = float(np.linalg.det(matrix))
        if abs(determinant) <= 1.0e-30:
            continue
        delta = np.asarray([potential[b] - potential[a], potential[c] - potential[a]])
        field = -np.linalg.solve(matrix, delta)
        area = 0.5 * abs(determinant)
        for node in (a, b, c):
            accumulated[node] += area * field
            weights[node] += area
    valid = weights > 0.0
    accumulated[valid] /= weights[valid, None]
    return np.linalg.norm(accumulated, axis=1)


def configure_axes(axis: plt.Axes) -> None:
    axis.set_xlabel("x (µm)")
    axis.set_ylabel("y (µm)")
    axis.set_aspect("equal", adjustable="box")
    axis.tick_params(colors=INK, labelsize=8)


def plot_spatial_fields(
    out_dir: Path,
    x: np.ndarray,
    y: np.ndarray,
    triangles: np.ndarray,
    sent_potential: np.ndarray,
    vela_potential: np.ndarray,
    sent_field: np.ndarray,
    vela_field: np.ndarray,
) -> None:
    triangulation = mtri.Triangulation(x, y, triangles)
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "axes.edgecolor": INK,
        "axes.labelcolor": INK,
        "text.color": INK,
        "figure.facecolor": "white",
    })

    potential_min = float(min(sent_potential.min(), vela_potential.min()))
    potential_max = float(max(sent_potential.max(), vela_potential.max()))
    difference = vela_potential - sent_potential
    diff_limit = float(np.percentile(np.abs(difference), 99.5))
    diff_limit = max(diff_limit, 1.0e-6)
    figure, axes = plt.subplots(1, 3, figsize=(14.2, 4.4), constrained_layout=True)
    for axis, values, title in zip(
        axes[:2],
        (sent_potential, vela_potential),
        ("Sentaurus potential", "Vela potential"),
    ):
        image = axis.tripcolor(
            triangulation,
            values,
            shading="gouraud",
            cmap="cividis",
            vmin=potential_min,
            vmax=potential_max,
        )
        axis.set_title(title, fontsize=11)
        configure_axes(axis)
        figure.colorbar(image, ax=axis, label="Potential (V)", shrink=0.86)
    image = axes[2].tripcolor(
        triangulation,
        difference,
        shading="gouraud",
        cmap="coolwarm",
        norm=TwoSlopeNorm(vmin=-diff_limit, vcenter=0.0, vmax=diff_limit),
    )
    axes[2].set_title("Vela − Sentaurus", fontsize=11)
    configure_axes(axes[2])
    figure.colorbar(image, ax=axes[2], label="Potential difference (V)", shrink=0.86)
    figure.suptitle(
        "Electrostatic potential at Id = 1e-4 A/µm",
        fontsize=14,
        fontweight="bold",
    )
    figure.savefig(out_dir / "potential_distribution_comparison.png", dpi=220)
    plt.close(figure)

    field_floor = max(float(np.percentile(sent_field[sent_field > 0.0], 1)), 1.0e3)
    sent_log = np.log10(np.maximum(sent_field, field_floor))
    vela_log = np.log10(np.maximum(vela_field, field_floor))
    log_min = float(min(sent_log.min(), vela_log.min()))
    log_max = float(max(sent_log.max(), vela_log.max()))
    log_ratio = np.log10(np.maximum(vela_field, field_floor) / np.maximum(sent_field, field_floor))
    ratio_limit = max(float(np.percentile(np.abs(log_ratio), 99.0)), 0.05)
    figure, axes = plt.subplots(1, 3, figsize=(14.2, 4.4), constrained_layout=True)
    for axis, values, title in zip(
        axes[:2],
        (sent_log, vela_log),
        ("Sentaurus electric field", "Vela electric field"),
    ):
        image = axis.tripcolor(
            triangulation,
            values,
            shading="gouraud",
            cmap="magma",
            vmin=log_min,
            vmax=log_max,
        )
        axis.set_title(title, fontsize=11)
        configure_axes(axis)
        figure.colorbar(image, ax=axis, label="log₁₀ |E| (V/m)", shrink=0.86)
    image = axes[2].tripcolor(
        triangulation,
        log_ratio,
        shading="gouraud",
        cmap="coolwarm",
        norm=TwoSlopeNorm(vmin=-ratio_limit, vcenter=0.0, vmax=ratio_limit),
    )
    axes[2].set_title("log₁₀(Vela / Sentaurus)", fontsize=11)
    configure_axes(axes[2])
    figure.colorbar(image, ax=axes[2], label="Electric-field log ratio (dex)", shrink=0.86)
    figure.suptitle(
        "Electric-field magnitude at Id = 1e-4 A/µm",
        fontsize=14,
        fontweight="bold",
    )
    figure.savefig(out_dir / "electric_field_distribution_comparison.png", dpi=220)
    plt.close(figure)


def parse_sentaurus_plt(path: Path) -> list[dict[str, float]]:
    text = path.read_text(encoding="utf-8", errors="replace")
    body = text.split("Data {", 1)[1].rsplit("}", 1)[0]
    values = np.asarray([
        float(value)
        for value in re.findall(r"[-+]?\d*\.\d+(?:[Ee][-+]?\d+)?", body)
    ])
    columns = 33
    if len(values) % columns != 0:
        raise ValueError(f"Sentaurus PLT value count {len(values)} is not divisible by {columns}")
    data = values.reshape((-1, columns))
    rows = []
    for point in data:
        rows.append({
            "voltage_V": float(point[18]),  # drain InnerVoltage
            "total_current_A_per_um": abs(float(point[23])),
            "electron_current_A_per_um": abs(float(point[21])),
            "hole_current_A_per_um": abs(float(point[22])),
            "time_parameter": float(point[0]),
        })
    return rows


def vela_iv_rows(prebias_path: Path, current_path: Path) -> list[dict[str, float]]:
    selected: dict[float, dict[str, float]] = {}
    for source, path in (("prebias", prebias_path), ("current_boundary", current_path)):
        for row in read_rows(path):
            if row.get("converged") != "1":
                continue
            voltage = float(row.get("inner_voltage_V") or row["bias_V"])
            current = abs(float(row["current_total_A_per_um"]))
            selected[round(voltage, 12)] = {
                "voltage_V": voltage,
                "total_current_A_per_um": current,
                "source_segment": source,
                "newton_iterations": int(row.get("newton_iterations") or 0),
            }
    return sorted(selected.values(), key=lambda row: row["voltage_V"])


def plot_iv(out_dir: Path, sent_rows: list[dict[str, float]], vela_rows: list[dict[str, float]]) -> None:
    sent_voltage = np.asarray([row["voltage_V"] for row in sent_rows])
    sent_current = np.asarray([max(row["total_current_A_per_um"], 1.0e-18) for row in sent_rows])
    vela_voltage = np.asarray([row["voltage_V"] for row in vela_rows])
    vela_current = np.asarray([max(row["total_current_A_per_um"], 1.0e-18) for row in vela_rows])

    figure, axes = plt.subplots(1, 2, figsize=(12.8, 4.8), constrained_layout=True)
    for axis in axes:
        axis.plot(sent_voltage, sent_current, color=BLUE, linewidth=2.0, label="Sentaurus")
        axis.plot(
            vela_voltage,
            vela_current,
            color=ORANGE,
            linewidth=1.8,
            linestyle="--",
            marker="o",
            markersize=3.5,
            markerfacecolor="white",
            markeredgewidth=1.0,
            label="Vela",
        )
        axis.axhline(1.0e-4, color=INK, linewidth=1.0, linestyle=":", label="BV criterion")
        axis.set_yscale("log")
        axis.yaxis.set_major_locator(LogLocator(base=10))
        axis.grid(True, which="major", color=GRID, linewidth=0.7)
        axis.set_xlabel("Drain voltage (V)")
        axis.set_ylabel("|Drain total current| (A/µm)")
    axes[0].set_title("Full reverse sweep", fontsize=11)
    axes[0].set_xlim(left=0.0)
    axes[0].set_ylim(1.0e-12, 3.0e-4)
    axes[1].set_title("Breakdown-region detail", fontsize=11)
    axes[1].set_xlim(5.75, 6.45)
    axes[1].set_ylim(1.0e-5, 2.0e-4)
    axes[1].legend(frameon=False, loc="lower right")
    figure.suptitle("NMOS breakdown I–V comparison", fontsize=14, fontweight="bold")
    figure.savefig(out_dir / "iv_curve_comparison.png", dpi=220)
    plt.close(figure)


def parse_vela_runtime(path: Path, label: str) -> dict[str, Any]:
    text = path.read_text(encoding="utf-8", errors="replace")
    elapsed = float(re.search(r"elapsed_seconds:\s*([0-9.]+)", text).group(1))
    iteration_rows = re.findall(
        r"^\[(?P<timestamp>[^]]+)\].*Newton iter (?P<iter>\d+) residual=(?P<residual>[^ ]+)",
        text,
        flags=re.MULTILINE,
    )
    updates = sum(int(iteration) > 0 for _, iteration, _ in iteration_rows)
    solve_groups = sum(int(iteration) == 0 for _, iteration, _ in iteration_rows)
    first_time = datetime.fromisoformat(iteration_rows[0][0])
    elapsed_trace = []
    for timestamp, iteration, residual in iteration_rows:
        elapsed_trace.append({
            "solver": label,
            "elapsed_s": (datetime.fromisoformat(timestamp) - first_time).total_seconds(),
            "iteration": int(iteration),
            "residual": float(residual),
        })
    return {
        "solver": label,
        "wallclock_s": elapsed,
        "newton_updates": updates,
        "solve_groups": solve_groups,
        "wallclock_s_per_update": elapsed / updates,
        "trace": elapsed_trace,
    }


def parse_sentaurus_runtime(path: Path) -> dict[str, Any]:
    text = path.read_text(encoding="utf-8", errors="replace")
    wallclock = float(re.search(r"wallclock:\s*([0-9.]+) s", text).group(1))
    cpu = float(re.search(r"total cpu:\s*([0-9.]+) s", text).group(1))
    segments = text.split("Computing step from")[1:]
    step_rows = []
    for segment in segments:
        iteration_lines = re.findall(r"^\s*(\d+)\s+[0-9.eE+-]+", segment, flags=re.MULTILINE)
        total_match = re.search(r"Total time:\s*([0-9.]+) s", segment)
        rhs_match = re.search(r"Rhs time:\s*([0-9.]+) s", segment)
        jac_match = re.search(r"Jacobian time:\s*([0-9.]+) s", segment)
        solve_match = re.search(r"Solve time:\s*([0-9.]+) s", segment)
        if not iteration_lines or not total_match:
            continue
        current_contact = "current-contact-equation" in segment
        step_rows.append({
            "current_contact": current_contact,
            "newton_updates": max(int(value) for value in iteration_lines),
            "total_s": float(total_match.group(1)),
            "rhs_s": float(rhs_match.group(1)) if rhs_match else 0.0,
            "jacobian_s": float(jac_match.group(1)) if jac_match else 0.0,
            "linear_solve_s": float(solve_match.group(1)) if solve_match else 0.0,
        })
    updates = sum(row["newton_updates"] for row in step_rows)
    current_rows = [row for row in step_rows if row["current_contact"]]
    return {
        "solver": "Sentaurus full current-boundary run",
        "wallclock_s": wallclock,
        "cpu_s": cpu,
        "nonlinear_steps": len(step_rows),
        "newton_updates": updates,
        "wallclock_s_per_update": wallclock / updates,
        "reported_step_total_s": sum(row["total_s"] for row in step_rows),
        "reported_rhs_s": sum(row["rhs_s"] for row in step_rows),
        "reported_jacobian_s": sum(row["jacobian_s"] for row in step_rows),
        "reported_linear_solve_s": sum(row["linear_solve_s"] for row in step_rows),
        "current_boundary_steps": len(current_rows),
        "current_boundary_newton_updates": sum(row["newton_updates"] for row in current_rows),
        "current_boundary_reported_total_s": sum(row["total_s"] for row in current_rows),
        "steps": step_rows,
    }


def plot_runtime(
    out_dir: Path,
    sentaurus: dict[str, Any],
    vela_current: dict[str, Any],
    vela_resistor: dict[str, Any],
) -> None:
    labels = ["Sentaurus\nfull run", "Vela\ncurrent switch", "Vela\nexternal resistor"]
    minutes = np.asarray([
        sentaurus["wallclock_s"],
        vela_current["wallclock_s"],
        vela_resistor["wallclock_s"],
    ]) / 60.0
    colors = [BLUE, ORANGE, GOLD]
    figure, axes = plt.subplots(1, 2, figsize=(12.8, 4.7), constrained_layout=True)
    bars = axes[0].bar(labels, minutes, color=colors, edgecolor=INK, linewidth=0.8)
    axes[0].set_ylabel("Wallclock (minutes)")
    axes[0].set_title("Observed end-to-end runtime", fontsize=11)
    axes[0].grid(axis="y", color=GRID, linewidth=0.7)
    for bar, value in zip(bars, minutes):
        axes[0].text(bar.get_x() + bar.get_width() / 2, value + 1.0, f"{value:.2f}", ha="center", fontsize=9)

    update_labels = ["Sentaurus\nall steps", "Vela\ncurrent switch", "Vela\nexternal resistor"]
    seconds_per_update = [
        sentaurus["wallclock_s_per_update"],
        vela_current["wallclock_s_per_update"],
        vela_resistor["wallclock_s_per_update"],
    ]
    bars = axes[1].bar(update_labels, seconds_per_update, color=colors, edgecolor=INK, linewidth=0.8)
    axes[1].set_ylabel("Wallclock / Newton update (s)")
    axes[1].set_title("Effective nonlinear iteration cost", fontsize=11)
    axes[1].grid(axis="y", color=GRID, linewidth=0.7)
    for bar, value in zip(bars, seconds_per_update):
        axes[1].text(bar.get_x() + bar.get_width() / 2, value + 0.06, f"{value:.2f}", ha="center", fontsize=9)
    figure.suptitle("Sentaurus and Vela runtime comparison", fontsize=14, fontweight="bold")
    figure.savefig(out_dir / "runtime_comparison.png", dpi=220)
    plt.close(figure)

    figure, axis = plt.subplots(figsize=(8.6, 4.8), constrained_layout=True)
    for runtime, color, label in (
        (vela_current, ORANGE, "Vela current switch"),
        (vela_resistor, GOLD, "Vela external resistor"),
    ):
        trace = runtime["trace"]
        x = np.asarray([row["elapsed_s"] / 60.0 for row in trace])
        y = np.arange(len(trace))
        axis.plot(x, y, color=color, linewidth=1.8, label=label)
    axis.set_xlabel("Elapsed wallclock (minutes)")
    axis.set_ylabel("Logged Newton states (cumulative)")
    axis.set_title("Vela runtime accumulates almost entirely inside Newton solves", fontsize=11)
    axis.grid(True, color=GRID, linewidth=0.7)
    axis.legend(frameon=False)
    figure.savefig(out_dir / "vela_newton_runtime_accumulation.png", dpi=220)
    plt.close(figure)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sentaurus-root", type=Path, default=DEFAULT_SENT)
    parser.add_argument("--sentaurus-plt", type=Path, default=DEFAULT_SENT_PLT)
    parser.add_argument("--sentaurus-log", type=Path, default=DEFAULT_SENT_LOG)
    parser.add_argument("--vela-state", type=Path, default=DEFAULT_VELA_STATE)
    parser.add_argument("--mesh", type=Path, default=DEFAULT_MESH)
    parser.add_argument("--vela-prebias-sweep", type=Path, default=DEFAULT_PREBIAS)
    parser.add_argument("--vela-current-sweep", type=Path, default=DEFAULT_CURRENT_SWEEP)
    parser.add_argument("--vela-current-log", type=Path, default=DEFAULT_CURRENT_LOG)
    parser.add_argument("--vela-resistor-log", type=Path, default=DEFAULT_RESISTOR_LOG)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    sent_potential_map = sent_scalar(args.sentaurus_root, "ElectrostaticPotential")
    sent_field_map = sent_vector(args.sentaurus_root, "ElectricField")
    nodes = sorted(sent_potential_map)
    x, y, triangles, _ = load_mesh(args.mesh, nodes)
    vela_state = {int(row["node_id"]): row for row in read_rows(args.vela_state)}
    sent_potential = np.asarray([sent_potential_map[node] for node in nodes])
    vela_potential = np.asarray([float(vela_state[node]["psi"]) for node in nodes])
    sent_field = np.asarray([
        100.0 * math.hypot(*sent_field_map[node]) for node in nodes
    ])
    vela_field = reconstruct_node_electric_field(x, y, triangles, vela_potential)
    plot_spatial_fields(
        args.out_dir, x, y, triangles,
        sent_potential, vela_potential, sent_field, vela_field,
    )

    field_rows = [
        {
            "node_id": node,
            "x_um": x[index],
            "y_um": y[index],
            "sentaurus_potential_V": sent_potential[index],
            "vela_potential_V": vela_potential[index],
            "potential_difference_V": vela_potential[index] - sent_potential[index],
            "sentaurus_electric_field_V_per_m": sent_field[index],
            "vela_electric_field_V_per_m": vela_field[index],
            "electric_field_log10_ratio": math.log10(
                max(vela_field[index], 1.0e3) / max(sent_field[index], 1.0e3)
            ),
        }
        for index, node in enumerate(nodes)
    ]
    write_rows(args.out_dir / "field_node_comparison.csv", field_rows)

    sent_iv = parse_sentaurus_plt(args.sentaurus_plt)
    vela_iv = vela_iv_rows(args.vela_prebias_sweep, args.vela_current_sweep)
    plot_iv(args.out_dir, sent_iv, vela_iv)
    write_rows(
        args.out_dir / "iv_curve_comparison.csv",
        [dict(solver="Sentaurus", **row) for row in sent_iv]
        + [dict(solver="Vela", **row) for row in vela_iv],
    )

    sent_runtime = parse_sentaurus_runtime(args.sentaurus_log)
    vela_current_runtime = parse_vela_runtime(args.vela_current_log, "Vela current switch")
    vela_resistor_runtime = parse_vela_runtime(args.vela_resistor_log, "Vela external resistor")
    plot_runtime(args.out_dir, sent_runtime, vela_current_runtime, vela_resistor_runtime)
    runtime_rows = []
    for runtime, scope in (
        (sent_runtime, "full 0-to-6 V prebias plus current-boundary run"),
        (vela_current_runtime, "checkpoint-resumed current-boundary run"),
        (vela_resistor_runtime, "checkpoint-resumed external-resistor run"),
    ):
        runtime_rows.append({
            "solver": runtime["solver"],
            "scope": scope,
            "wallclock_s": runtime["wallclock_s"],
            "newton_updates": runtime["newton_updates"],
            "solve_groups_or_steps": (
                runtime["nonlinear_steps"]
                if "nonlinear_steps" in runtime
                else runtime["solve_groups"]
            ),
            "wallclock_s_per_update": runtime["wallclock_s_per_update"],
        })
    write_rows(args.out_dir / "runtime_summary.csv", runtime_rows)

    performance = {
        "sentaurus": sent_runtime,
        "vela_current_switch": {key: value for key, value in vela_current_runtime.items() if key != "trace"},
        "vela_external_resistor": {key: value for key, value in vela_resistor_runtime.items() if key != "trace"},
        "comparisons": {
            "vela_current_wallclock_over_sentaurus": (
                vela_current_runtime["wallclock_s"] / sent_runtime["wallclock_s"]
            ),
            "vela_resistor_wallclock_over_sentaurus": (
                vela_resistor_runtime["wallclock_s"] / sent_runtime["wallclock_s"]
            ),
            "vela_current_per_update_over_sentaurus": (
                vela_current_runtime["wallclock_s_per_update"]
                / sent_runtime["wallclock_s_per_update"]
            ),
            "vela_resistor_per_update_over_sentaurus": (
                vela_resistor_runtime["wallclock_s_per_update"]
                / sent_runtime["wallclock_s_per_update"]
            ),
        },
        "field_summary": {
            "sentaurus_node_field_peak_V_per_m": float(sent_field.max()),
            "vela_reconstructed_node_field_peak_V_per_m": float(vela_field.max()),
            "potential_abs_error_p50_V": float(np.percentile(np.abs(vela_potential - sent_potential), 50)),
            "potential_abs_error_p95_V": float(np.percentile(np.abs(vela_potential - sent_potential), 95)),
        },
    }
    (args.out_dir / "performance_summary.json").write_text(
        json.dumps(performance, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(performance, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
