#!/usr/bin/env python3
"""Frozen-Q TransportModels audit of Enormal, mobility, and driving field."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[1]
BASELINE = REPO_ROOT / "build-release/reference_tcad/transportmodels_sentaurus2022/vela_baseline"
OUTPUT_ROOT = BASELINE / "dg_surface_mobility_frozen_q_2026-08-21"
RUNNER = REPO_ROOT / "build-release/vela_example_runner.exe"
BASE_CONFIG = BASELINE / "frozen_q_oracle_vg1_vd2_run01/frozen_q_idvd_2V.json"
HYBRID_RESTART = BASELINE / "frozen_q_oracle_vg1_vd2_run01/vela_state_with_sentaurus_q.csv"
CORRECTED_MATERIALS = BASELINE / "dg_parameter_fixed_state_sweep_2026-08-21/materials_sentaurus2022_dg_band_drive.json"
MESH = BASELINE / "generated/vela/mesh.json"
SENTAURUS_FIELDS = BASELINE / "generated/sim_fields/dg_idvd/fields"
REPORT_JSON = REPO_ROOT / "docs/validation/transportmodels_dg_surface_mobility_audit_2026-08-21.json"
REPORT_MD = REPO_ROOT / "docs/validation/transportmodels_dg_surface_mobility_audit_2026-08-21.md"
SENTAURUS_CURRENT = 0.000705525753105

VARIANTS: tuple[dict[str, Any], ...] = (
    {"name": "implicit_all_interfaces", "label": "Lombardi implicit interfaces", "model": "masetti_field_lombardi", "explicit": False},
    {"name": "explicit_channel_interface", "label": "Lombardi explicit channel", "model": "masetti_field_lombardi", "explicit": True},
    {"name": "no_enormal", "label": "No Enormal", "model": "masetti_field", "explicit": False},
    {"name": "no_high_field", "label": "No high-field saturation", "model": "masetti_lombardi", "explicit": True},
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def percentile(values: list[float], p: float) -> float:
    ordered = sorted(values)
    if not ordered:
        return math.nan
    position = (len(ordered) - 1) * p
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    if lower == upper:
        return ordered[lower]
    return ordered[lower] * (upper - position) + ordered[upper] * (position - lower)


def make_config(variant: dict[str, Any], run_dir: Path) -> Path:
    config = json.loads(BASE_CONFIG.read_text(encoding="utf-8"))
    config["materials_file"] = str(CORRECTED_MATERIALS.resolve())
    mobility = config["solver"]["mobility"]
    mobility["model"] = variant["model"]
    if variant["explicit"]:
        mobility["surface"] = {
            "surface_region": "R.Substrate",
            "surface_interface": ["R.Substrate", "R.Gateox"],
        }
    else:
        mobility.pop("surface", None)
    config["_comment"] = "TransportModels phase-6 Frozen-Q mobility audit: " + variant["label"]
    config["output_csv"] = str((run_dir / "endpoint.csv").resolve())
    config["log_file"] = str((run_dir / "endpoint.log").resolve())
    config["sweep"].update(
        {
            "initial_state_file": str(HYBRID_RESTART.resolve()),
            "write_state_file": str((run_dir / "endpoint_final_state.csv").resolve()),
            "write_state_every_point_prefix": str((run_dir / "endpoint_state").resolve()),
        }
    )
    config["sweep"]["diagnostics"]["terminal_balance"]["csv_file"] = str(
        (run_dir / "terminal_balance.csv").resolve()
    )
    path = run_dir / "config.json"
    path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    return path


def make_probe(config_path: Path, run_dir: Path) -> Path:
    config = json.loads(config_path.read_text(encoding="utf-8"))
    fields_dir = run_dir / "state_fields"
    fields_dir.mkdir(parents=True, exist_ok=True)
    with (run_dir / "endpoint_final_state.csv").open(newline="", encoding="utf-8") as handle:
        state_rows = list(csv.DictReader(handle))
    for output_name, state_name in (
        ("ElectrostaticPotential", "psi"),
        ("eQuasiFermiPotential", "phin"),
        ("hQuasiFermiPotential", "phip"),
    ):
        with (fields_dir / f"{output_name}_region0.csv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(("node_id", "component0"))
            writer.writerows((row["node_id"], row[state_name]) for row in state_rows)
    config["simulation_type"] = "edge_mobility_probe"
    config["state_fields_dir"] = str(fields_dir.resolve())
    config["output_csv"] = str((run_dir / "edge_mobility.csv").resolve())
    path = run_dir / "edge_mobility_probe.json"
    path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    return path


def execute(variant: dict[str, Any]) -> dict[str, Any]:
    run_dir = OUTPUT_ROOT / variant["name"]
    run_dir.mkdir(parents=True, exist_ok=True)
    config = make_config(variant, run_dir)
    endpoint = run_dir / "endpoint.csv"
    final_state = run_dir / "endpoint_final_state.csv"
    if endpoint.exists() and endpoint.stat().st_size and final_state.exists() and final_state.stat().st_size:
        process = subprocess.CompletedProcess([], 0, "reused existing endpoint", "")
    else:
        process = subprocess.run([str(RUNNER), "--config", str(config)], cwd=REPO_ROOT, text=True, capture_output=True)
    (run_dir / "console.log").write_text(process.stdout + "\n--- STDERR ---\n" + process.stderr, encoding="utf-8")
    probe = make_probe(config, run_dir)
    probe_process = subprocess.run([str(RUNNER), "--config", str(probe)], cwd=REPO_ROOT, text=True, capture_output=True)
    (run_dir / "edge_probe.console.log").write_text(probe_process.stdout + "\n--- STDERR ---\n" + probe_process.stderr, encoding="utf-8")
    with endpoint.open(newline="", encoding="utf-8") as handle:
        row = list(csv.DictReader(handle))[-1]
    current = float(row["current_total_A_per_um"])
    edges = run_dir / "edge_mobility.csv"
    return {
        "name": variant["name"], "label": variant["label"],
        "runner_exit_code": process.returncode, "probe_exit_code": probe_process.returncode,
        "converged": row["converged"] == "1", "current_A_per_um": current,
        "current_absolute_relative_error": abs(current - SENTAURUS_CURRENT) / SENTAURUS_CURRENT,
        "aggregate": {
            "mean_electron_mobility_m2_V_s": float(row["mean_electron_mobility_m2_V_s"]),
            "min_electron_mobility_m2_V_s": float(row["min_electron_mobility_m2_V_s"]),
            "mean_electron_high_field_drive_V_per_cm": float(row["mean_electron_high_field_drive_V_per_cm"]),
            "max_electric_field_V_per_cm": float(row["max_electric_field_V_per_cm"]),
        },
        "config": str(config), "config_sha256": sha256(config),
        "endpoint_csv": str(endpoint), "endpoint_sha256": sha256(endpoint),
        "edge_csv": str(edges), "edge_sha256": sha256(edges),
        "state_csv": str(run_dir / "endpoint_final_state.csv"),
    }


def channel_geometry() -> tuple[set[tuple[int, int]], dict[int, tuple[float, float]], list[tuple[int, tuple[int, int, int]]]]:
    mesh = json.loads(MESH.read_text(encoding="utf-8"))
    nodes = {int(row["id"]): (float(row["x"]), float(row["y"])) for row in mesh["nodes"]}
    edge_regions: dict[tuple[int, int], set[int]] = {}
    substrate_triangles: list[tuple[int, tuple[int, int, int]]] = []
    for tri in mesh["triangles"]:
        ids = tuple(int(value) for value in tri["node_ids"])
        region = int(tri["region_id"])
        if region == 3:
            substrate_triangles.append((int(tri["id"]), ids))
        for a, b in ((ids[0], ids[1]), (ids[1], ids[2]), (ids[2], ids[0])):
            edge_regions.setdefault(tuple(sorted((a, b))), set()).add(region)
    channel = {pair for pair, regions in edge_regions.items() if regions == {0, 3}}
    return channel, nodes, substrate_triangles


def state_psi(path: Path) -> dict[int, float]:
    with path.open(newline="", encoding="utf-8") as handle:
        return {int(row["node_id"]): float(row["psi"]) for row in csv.DictReader(handle)}


def channel_enormal(
    channel: set[tuple[int, int]], nodes: dict[int, tuple[float, float]],
    triangles: list[tuple[int, tuple[int, int, int]]], psi: dict[int, float],
) -> dict[tuple[int, int], float]:
    adjacent: dict[tuple[int, int], tuple[int, int, int]] = {}
    for _, ids in triangles:
        for pair in (tuple(sorted((ids[0], ids[1]))), tuple(sorted((ids[1], ids[2]))), tuple(sorted((ids[2], ids[0])))):
            if pair in channel:
                adjacent[pair] = ids
    result = {}
    for pair, ids in adjacent.items():
        (x0, y0), (x1, y1), (x2, y2) = (nodes[node] for node in ids)
        det = (x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0)
        gx = ((psi[ids[1]] - psi[ids[0]]) * (y2 - y0) - (psi[ids[2]] - psi[ids[0]]) * (y1 - y0)) / det
        gy = ((x1 - x0) * (psi[ids[2]] - psi[ids[0]]) - (x2 - x0) * (psi[ids[1]] - psi[ids[0]])) / det
        (xa, ya), (xb, yb) = nodes[pair[0]], nodes[pair[1]]
        length = math.hypot(xb - xa, yb - ya)
        nx, ny = -(yb - ya) / length, (xb - xa) / length
        result[pair] = abs(gx * nx + gy * ny) * 1.0e6  # V/um -> V/m
    return result


def spatial_metrics(row: dict[str, Any], channel: set[tuple[int, int]], nodes: dict[int, tuple[float, float]], triangles: list[tuple[int, tuple[int, int, int]]]) -> None:
    psi = state_psi(Path(row["state_csv"]))
    enormal = channel_enormal(channel, nodes, triangles, psi)
    with Path(row["edge_csv"]).open(newline="", encoding="utf-8") as handle:
        edge_rows = [edge for edge in csv.DictReader(handle) if tuple(sorted((int(edge["node0"]), int(edge["node1"])))) in channel]
    mobilities = [float(edge["electron_final_mobility_m2_V_s"]) * 1.0e4 for edge in edge_rows]
    drives = [float(edge["electron_mobility_field_V_m"]) / 100.0 for edge in edge_rows]
    normal_fields = [value / 100.0 for value in enormal.values()]
    row["channel"] = {
        "edge_count": len(edge_rows),
        "enormal_V_per_cm": {"median": statistics.median(normal_fields), "p95": percentile(normal_fields, 0.95)},
        "electron_mobility_cm2_V_s": {"median": statistics.median(mobilities), "p05": percentile(mobilities, 0.05)},
        "electron_drive_V_per_cm": {"median": statistics.median(drives), "p95": percentile(drives, 0.95)},
    }


def sentaurus_metrics(channel: set[tuple[int, int]]) -> dict[str, Any]:
    node_ids = {node for pair in channel for node in pair}
    def values(name: str) -> list[float]:
        with (SENTAURUS_FIELDS / f"{name}_region3.csv").open(newline="", encoding="utf-8") as handle:
            return [float(row["component0"]) for row in csv.DictReader(handle) if int(row["node_id"]) in node_ids]
    normal = [abs(value) for value in values("eEnormal")]
    mobility = values("eMobility")
    drive = [abs(value) for value in values("eEparallel")]
    return {
        "node_count": len(normal),
        "enormal_V_per_cm": {"median": statistics.median(normal), "p95": percentile(normal, 0.95)},
        "electron_mobility_cm2_V_s": {"median": statistics.median(mobility), "p05": percentile(mobility, 0.05)},
        "electron_drive_V_per_cm": {"median": statistics.median(drive), "p95": percentile(drive, 0.95)},
    }


def plot(rows: list[dict[str, Any]], sentaurus: dict[str, Any]) -> tuple[Path, Path]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    labels = [row["label"].replace(" ", "\n") for row in rows]
    colors = ["#4C78A8", "#54A24B", "#E45756", "#F58518"]
    fig, axes = plt.subplots(2, 2, figsize=(13.2, 8.4))
    current_errors = [100.0 * row["current_absolute_relative_error"] for row in rows]
    axes[0, 0].bar(labels, current_errors, color=colors)
    axes[0, 0].set_ylabel("Absolute Id error (%)")
    axes[0, 0].set_title("Frozen-Q terminal current")
    metrics = (
        ("enormal_V_per_cm", "Channel median Enormal (MV/cm)", 1.0e6),
        ("electron_mobility_cm2_V_s", "Channel median electron mobility (cm²/Vs)", 1.0),
        ("electron_drive_V_per_cm", "Channel median drive (kV/cm)", 1.0e3),
    )
    for ax, (key, title, scale) in zip((axes[0, 1], axes[1, 0], axes[1, 1]), metrics):
        values = [row["channel"][key]["median"] / scale for row in rows]
        ax.bar(labels, values, color=colors)
        reference = sentaurus[key]["median"] / scale
        ax.axhline(reference, color="black", linestyle="--", linewidth=1.4, label="Sentaurus")
        ax.set_title(title)
        ax.legend(fontsize=8)
    for ax in axes.flat:
        ax.grid(axis="y", alpha=0.25)
        ax.tick_params(axis="x", labelsize=8)
    fig.suptitle("TransportModels DG Frozen-Q mobility audit", fontsize=15, y=0.98)
    fig.text(0.5, 0.945, "Vg = 1 V, Vd = 2 V; corrected materials; channel = R.Substrate/R.Gateox", ha="center", fontsize=10)
    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.14, top=0.88, hspace=0.45, wspace=0.25)
    png = OUTPUT_ROOT / "surface_mobility_comparison.png"
    svg = OUTPUT_ROOT / "surface_mobility_comparison.svg"
    fig.savefig(png, dpi=180)
    fig.savefig(svg)
    plt.close(fig)
    return png, svg


def markdown(report: dict[str, Any]) -> str:
    lines = [
        "# TransportModels DG surface-mobility audit",
        "",
        "Frozen Sentaurus quantum potential at `Vg=1 V`, `Vd=2 V`; all cases use the",
        "corrected material contract. This isolates the classical mobility/transport path.",
        "",
        "| Variant | Converged | Id error | Channel median Enormal (V/cm) | Mobility (cm2/Vs) | Drive (V/cm) |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in report["results"]:
        channel = row["channel"]
        lines.append(
            f"| {row['label']} | {row['converged']} | {100 * row['current_absolute_relative_error']:.4f}% | "
            f"{channel['enormal_V_per_cm']['median']:.6g} | {channel['electron_mobility_cm2_V_s']['median']:.6g} | "
            f"{channel['electron_drive_V_per_cm']['median']:.6g} |"
        )
    best = min((row for row in report["results"] if row["converged"]), key=lambda row: row["current_absolute_relative_error"])
    lines.extend(
        [
            "", "## Decision", "",
            f"- Best Frozen-Q terminal-current agreement: **{best['label']}**.",
            "- The explicit channel selector is preferred only if it improves current and spatial",
            "  mobility agreement; otherwise the existing model remains frozen to avoid fitting",
            "  one endpoint at the expense of the curve.",
            "- Vela Enormal is reconstructed from the adjacent substrate triangle and the exact",
            "  R.Substrate/R.Gateox interface normal; Sentaurus values are native `eEnormal`.",
            "", f"Figure: `{report['artifacts']['png']}`", "",
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
            for path_key, hash_key in (("config", "config_sha256"), ("endpoint_csv", "endpoint_sha256"), ("edge_csv", "edge_sha256")):
                assert sha256(Path(row[path_key])) == row[hash_key]
        print("TransportModels DG surface mobility audit check: PASS")
        return 0
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    rows = []
    with ThreadPoolExecutor(max_workers=3) as pool:
        futures = {pool.submit(execute, variant): variant for variant in VARIANTS}
        for future in as_completed(futures):
            rows.append(future.result())
    order = {variant["name"]: index for index, variant in enumerate(VARIANTS)}
    rows.sort(key=lambda row: order[row["name"]])
    channel, nodes, triangles = channel_geometry()
    for row in rows:
        spatial_metrics(row, channel, nodes, triangles)
    sentaurus = sentaurus_metrics(channel)
    png, svg = plot(rows, sentaurus)
    report = {
        "schema": "vela.transportmodels.dg_surface_mobility_frozen_q.v1",
        "status": "pass" if all(row["converged"] for row in rows) else "partial",
        "as_of": "2026-08-21", "sentaurus_current_A_per_um": SENTAURUS_CURRENT,
        "sentaurus_channel": sentaurus, "results": rows,
        "artifacts": {"png": str(png), "svg": str(svg)},
    }
    REPORT_JSON.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    REPORT_MD.write_text(markdown(report), encoding="utf-8")
    print(json.dumps({"status": report["status"], "sentaurus_channel": sentaurus, "results": [{"name": row["name"], "current_error_percent": 100 * row["current_absolute_relative_error"], "channel": row["channel"]} for row in rows]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
