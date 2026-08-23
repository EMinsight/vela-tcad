#!/usr/bin/env python3
"""Audit the TransportModels DG Eq. 231 residual at a fixed 2 V state."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
from collections import defaultdict
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
BASELINE = REPO_ROOT / "build-release/reference_tcad/transportmodels_sentaurus2022/vela_baseline"
DEFAULT_OUTPUT = BASELINE / "fixed_state_residual_audit_vg1_vd2_p1_direct_run01"
DEFAULT_RUNNER = REPO_ROOT / "build-release/vela_example_runner.exe"
BASE_CONFIG = (
    BASELINE
    / "workflow_dg_idvd_resume_1p9_direct2p0_fixed1p5_outer200_run01"
    / "00_dg_idvd_curve.json"
)
HYBRID_RESTART = (
    BASELINE
    / "frozen_q_oracle_vg1_vd2_run01"
    / "vela_state_with_sentaurus_q.csv"
)
DEFAULT_REPORT_JSON = (
    REPO_ROOT / "docs/validation/transportmodels_dg_fixed_state_residual_audit_2026-08-20.json"
)
DEFAULT_REPORT_MD = (
    REPO_ROOT / "docs/validation/transportmodels_dg_fixed_state_residual_audit_2026-08-20.md"
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def make_config(output_dir: Path) -> tuple[Path, Path]:
    config = json.loads(BASE_CONFIG.read_text(encoding="utf-8"))
    prefix = output_dir / "fixed_state_eq231"
    config["_comment"] = (
        "Phase-2 fixed-state Eq. 231 residual audit. The converged Vela DD "
        "variables are combined with Sentaurus Qn/Phi-like fields; the "
        "diagnostic is assembled from the initial state before a DG update."
    )
    for contact in config["contacts"]:
        if contact["name"].lower() == "drain":
            contact["bias"] = 2.0
    solver = config["solver"]
    solver["verbose"] = False
    quantum = solver["electron_quantum_potential"]
    quantum.update(
        {
            "enabled": True,
            "coupling_mode": "outer",
            "formulation": "potential_based",
            "include_insulators": True,
            "global_discretization": "p1_direct",
            "residual_diagnostic_prefix": str(prefix.resolve()),
            "residual_diagnostic_use_initial_state": True,
            "outer_max_iterations": 1,
            "max_iterations": 1,
        }
    )
    config["output_csv"] = str((output_dir / "fixed_state_probe.csv").resolve())
    config["log_file"] = str((output_dir / "fixed_state_probe.log").resolve())
    config["sweep"].update(
        {
            "start": 2.0,
            "stop": 2.0,
            "step": 0.1,
            "bias_points": [2.0],
            "initial_state_file": str(HYBRID_RESTART.resolve()),
            "write_vtk": False,
            "write_state_file": str((output_dir / "fixed_state_probe_final.csv").resolve()),
            "write_state_every_point_prefix": str((output_dir / "fixed_state_probe_state").resolve()),
        }
    )
    config_path = output_dir / "fixed_state_residual_audit.json"
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    return config_path, prefix


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def audit(prefix: Path) -> dict[str, object]:
    nodes_path = Path(str(prefix) + "_nodes.csv")
    cells_path = Path(str(prefix) + "_cells.csv")
    regions_path = Path(str(prefix) + "_regions.csv")
    summary_path = Path(str(prefix) + "_summary.txt")
    nodes = read_rows(nodes_path)
    cells = read_rows(cells_path)
    regions = read_rows(regions_path)
    free = [
        row
        for row in nodes
        if row["is_active"] == "1" and row["is_dirichlet"] == "0"
    ]
    free.sort(key=lambda row: abs(float(row["raw_total"])), reverse=True)
    top_nodes = []
    for row in free[:20]:
        components = {
            key: float(row[key])
            for key in (
                "stiffness",
                "gradient_squared",
                "reaction",
                "interface_boundary",
            )
        }
        dominant = max(components, key=lambda key: abs(components[key]))
        top_nodes.append(
            {
                "node_id": int(row["node_id"]),
                "x_um": float(row["x_internal"]),
                "y_um": float(row["y_internal"]),
                "raw_total": float(row["raw_total"]),
                "absolute_raw_total": abs(float(row["raw_total"])),
                "dominant_component": dominant,
                "components": components,
            }
        )

    cell_groups: dict[int, list[dict[str, str]]] = defaultdict(list)
    for row in cells:
        cell_groups[int(row["cell_id"])].append(row)
    cell_hotspots = []
    for cell_id, rows in cell_groups.items():
        score = max(abs(float(row["total"])) for row in rows)
        cell_hotspots.append(
            {
                "cell_id": cell_id,
                "region_name": rows[0]["region_name"],
                "material": rows[0]["material"],
                "is_interface_cell": rows[0]["is_interface_cell"] == "1",
                "interface_pairs": rows[0]["interface_pairs"],
                "max_absolute_local_total": score,
                "centroid_um": {
                    "x": sum(float(row["x_internal"]) for row in rows) / len(rows),
                    "y": sum(float(row["y_internal"]) for row in rows) / len(rows),
                },
            }
        )
    cell_hotspots.sort(key=lambda row: row["max_absolute_local_total"], reverse=True)

    region_metrics = []
    total_region_l1 = sum(float(row["total_l1_free"]) for row in regions)
    for row in regions:
        l1 = float(row["total_l1_free"])
        interface_l1 = float(row["interface_total_l1_free"])
        region_metrics.append(
            {
                "region_id": int(row["region_id"]),
                "region_name": row["region_name"],
                "material": row["material"],
                "cell_count": int(row["cell_count"]),
                "interface_cell_count": int(row["interface_cell_count"]),
                "total_l1_free": l1,
                "global_l1_share": l1 / total_region_l1 if total_region_l1 else 0.0,
                "interface_total_l1_free": interface_l1,
                "interface_share_within_region": interface_l1 / l1 if l1 else 0.0,
                "max_cell_residual_free": float(row["max_cell_residual_free"]),
                "max_cell_id": int(row["max_cell_id"]),
            }
        )
    region_metrics.sort(key=lambda row: row["total_l1_free"], reverse=True)

    summary_values: dict[str, float | int] = {}
    for line in summary_path.read_text(encoding="utf-8").splitlines():
        key, value = line.split("=", 1)
        summary_values[key] = int(value) if key == "max_free_node" else float(value)
    component_l1 = {
        "stiffness": float(summary_values["stiffness_l1_free"]),
        "gradient_squared": float(summary_values["gradient_squared_l1_free"]),
        "reaction": float(summary_values["reaction_l1_free"]),
        "interface_boundary": float(summary_values["interface_boundary_l1_free"]),
    }
    component_total = sum(component_l1.values())
    component_l1_share = {
        key: value / component_total if component_total else 0.0
        for key, value in component_l1.items()
    }
    return {
        "free_active_node_count": len(free),
        "cell_count": len(cell_groups),
        "summary": summary_values,
        "component_l1": component_l1,
        "component_l1_share": component_l1_share,
        "top_nodes": top_nodes,
        "top_cells": cell_hotspots[:20],
        "regions": region_metrics,
        "raw_files": {
            "nodes": str(nodes_path),
            "cells": str(cells_path),
            "regions": str(regions_path),
            "summary": str(summary_path),
        },
    }


def plot_hotspots(prefix: Path, output: Path, audit_data: dict[str, object]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection

    node_rows = read_rows(Path(str(prefix) + "_nodes.csv"))
    cell_rows = read_rows(Path(str(prefix) + "_cells.csv"))
    free = [
        row
        for row in node_rows
        if row["is_active"] == "1" and row["is_dirichlet"] == "0"
    ]
    node_abs = [abs(float(row["raw_total"])) for row in free]
    node_max = max(node_abs)
    node_values = [math.log10(max(value / node_max, 1.0e-12)) for value in node_abs]

    grouped: dict[int, list[dict[str, str]]] = defaultdict(list)
    for row in cell_rows:
        grouped[int(row["cell_id"])].append(row)
    polygons = []
    scores = []
    interface_segments = []
    for rows in grouped.values():
        if len(rows) != 3:
            continue
        polygon = [(float(row["x_internal"]), float(row["y_internal"])) for row in rows]
        polygons.append(polygon)
        scores.append(max(abs(float(row["total"])) for row in rows))
        if rows[0]["is_interface_cell"] == "1":
            interface_segments.append(polygon + [polygon[0]])
    cell_max = max(scores)
    cell_values = [math.log10(max(value / cell_max, 1.0e-12)) for value in scores]

    plt.rcParams.update({"font.size": 10, "axes.titlesize": 12, "axes.labelsize": 10})
    fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.9), constrained_layout=False)
    scatter = axes[0].scatter(
        [float(row["x_internal"]) for row in free],
        [float(row["y_internal"]) for row in free],
        c=node_values,
        cmap="Blues",
        vmin=-12,
        vmax=0,
        s=7,
        linewidths=0,
    )
    top = audit_data["top_nodes"][0]
    axes[0].scatter(
        [top["x_um"]], [top["y_um"]], marker="x", s=75, linewidths=1.8, color="#C45A20"
    )
    axes[0].annotate(
        f"max node {top['node_id']}",
        (top["x_um"], top["y_um"]),
        xytext=(8, 8),
        textcoords="offset points",
        color="#7A3515",
    )
    axes[0].set_title("Free-node Eq. 231 residual magnitude")
    axes[0].set_xlabel("Lateral position (µm)")
    axes[0].set_ylabel("Surface-normal depth (µm)")
    cbar0 = fig.colorbar(scatter, ax=axes[0], orientation="horizontal", pad=0.17, shrink=0.86)
    cbar0.set_label("log10(|raw residual| / node maximum)")

    collection = PolyCollection(
        polygons, array=cell_values, cmap="Blues", clim=(-12, 0), edgecolors="none"
    )
    axes[1].add_collection(collection)
    for segment in interface_segments:
        axes[1].plot(
            [point[0] for point in segment],
            [point[1] for point in segment],
            color="#C45A20",
            linewidth=0.16,
            alpha=0.18,
        )
    axes[1].autoscale_view()
    axes[1].set_title("Cell-local Eq. 231 residual hotspot")
    axes[1].set_xlabel("Lateral position (µm)")
    axes[1].set_ylabel("Surface-normal depth (µm)")
    cbar1 = fig.colorbar(collection, ax=axes[1], orientation="horizontal", pad=0.17, shrink=0.86)
    cbar1.set_label("log10(max local |residual| / cell maximum)")
    for axis in axes:
        axis.invert_yaxis()
        axis.set_aspect("equal", adjustable="box")
        axis.grid(False)
    fig.subplots_adjust(left=0.07, right=0.985, bottom=0.22, top=0.79, wspace=0.28)
    fig.suptitle(
        "TransportModels DG fixed-state residual audit", fontsize=16, y=0.975
    )
    fig.text(
        0.5,
        0.925,
        "Vg = 1.0 V, Vd = 2.0 V; p1_direct Eq. 231; orange outlines mark interface cells",
        ha="center",
        color="#4B5563",
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220, bbox_inches="tight", facecolor="white")
    fig.savefig(output.with_suffix(".svg"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


def render_markdown(report: dict[str, object]) -> str:
    audit_data = report["audit"]
    top = audit_data["top_nodes"][0]
    regions = audit_data["regions"][:5]
    component_share = audit_data["component_l1_share"]
    region_rows = "\n".join(
        f"| {row['region_name']} | {row['material']} | {row['global_l1_share']:.2%} | "
        f"{row['interface_share_within_region']:.2%} | {row['max_cell_id']} |"
        for row in regions
    )
    return f"""# TransportModels DG fixed-state residual audit

Work point: Vg = 1.0 V, Vd = 2.0 V

Status: **{report['status']}**

## Fixed-state definition

The DD variables come from the converged Vela self-consistent 2 V state. The
electron quantum potential and continuous potential-like field come from the
Sentaurus 2022 DG state. The p1_direct Eq. 231 operator is assembled from this
initial state before applying a DG update.

## Main result

- Maximum free-node raw residual: `{top['absolute_raw_total']:.12g}` at node
  `{top['node_id']}` ({top['x_um']:.6g} µm, {top['y_um']:.6g} µm).
- Dominant component at that node: `{top['dominant_component']}`.
- Global component L1 shares: stiffness {component_share['stiffness']:.2%},
  gradient-squared {component_share['gradient_squared']:.2%}, reaction
  {component_share['reaction']:.2%}, explicit interface boundary
  {component_share['interface_boundary']:.2%}.

## Region ranking

| Region | Material | Share of global free-node L1 | Interface share within region | Max cell |
|---|---|---:|---:|---:|
{region_rows}

The raw residual is an integrated discrete equation value and is best used for
relative localization and component attribution. It is not a voltage error and
must not be compared directly with the Qn field error in mV.

## Provenance

- Config SHA-256: `{report['hashes']['config']}`
- Fixed hybrid state SHA-256: `{report['hashes']['hybrid_restart']}`
- Node residual CSV SHA-256: `{report['hashes']['nodes']}`
- Cell residual CSV SHA-256: `{report['hashes']['cells']}`
"""


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, default=DEFAULT_RUNNER)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--report-json", type=Path, default=DEFAULT_REPORT_JSON)
    parser.add_argument("--report-md", type=Path, default=DEFAULT_REPORT_MD)
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--rerender-only", action="store_true")
    args = parser.parse_args()

    if args.check:
        report_path = args.report_json.resolve()
        report = json.loads(report_path.read_text(encoding="utf-8"))
        paths = report["audit"]["raw_files"]
        expected = report["hashes"]
        actual = {
            "config": sha256(Path(report["paths"]["config"])),
            "hybrid_restart": sha256(Path(report["paths"]["hybrid_restart"])),
            "nodes": sha256(Path(paths["nodes"])),
            "cells": sha256(Path(paths["cells"])),
            "regions": sha256(Path(paths["regions"])),
        }
        mismatches = {
            key: {"expected": expected[key], "actual": value}
            for key, value in actual.items()
            if value != expected[key]
        }
        recomputed = audit(Path(paths["nodes"].removesuffix("_nodes.csv")))
        saved_top = report["audit"]["top_nodes"][0]
        current_top = recomputed["top_nodes"][0]
        if current_top["node_id"] != saved_top["node_id"] or not math.isclose(
            current_top["absolute_raw_total"],
            saved_top["absolute_raw_total"],
            rel_tol=1.0e-12,
            abs_tol=0.0,
        ):
            mismatches["top_node"] = {
                "expected": saved_top,
                "actual": current_top,
            }
        if mismatches:
            print(json.dumps(mismatches, indent=2))
            return 1
        print("TransportModels DG fixed-state residual audit check: PASS")
        return 0

    if args.rerender_only:
        report = json.loads(args.report_json.resolve().read_text(encoding="utf-8"))
        prefix = Path(report["audit"]["raw_files"]["nodes"].removesuffix("_nodes.csv"))
        audit_data = audit(prefix)
        figure = Path(report["paths"]["figure_png"])
        plot_hotspots(prefix, figure, audit_data)
        print(figure)
        return 0

    if not HYBRID_RESTART.is_file():
        raise FileNotFoundError(
            f"Run scripts/run_transportmodels_dg_frozen_q_oracle.py first: {HYBRID_RESTART}"
        )
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    config_path, prefix = make_config(output_dir)
    if not args.execute:
        print(config_path)
        return 0

    process = subprocess.run(
        [str(args.runner.resolve()), "--config", str(config_path)],
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        encoding="utf-8",
        errors="replace",
    )
    console = output_dir / "fixed_state_residual_audit.console.log"
    console.write_text(
        process.stdout + ("\n[stderr]\n" + process.stderr if process.stderr else ""),
        encoding="utf-8",
    )
    required = [Path(str(prefix) + suffix) for suffix in ("_nodes.csv", "_cells.csv", "_regions.csv", "_summary.txt")]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise RuntimeError(
            f"Residual diagnostic was not emitted (runner exit {process.returncode}): {missing}; see {console}"
        )

    audit_data = audit(prefix)
    figure = output_dir / "fixed_state_eq231_residual_hotspots.png"
    plot_hotspots(prefix, figure, audit_data)
    report: dict[str, object] = {
        "schema": "vela.transportmodels.dg_fixed_state_residual_audit.v1",
        "status": "pass",
        "work_point": {"gate_bias_V": 1.0, "drain_bias_V": 2.0},
        "operator": {
            "equation": "Sentaurus Eq. 231 potential-based electron DG form",
            "global_discretization": "p1_direct",
            "include_insulators": True,
            "diagnostic_use_initial_state": True,
        },
        "runner": {
            "exit_code": process.returncode,
            "note": "A nonzero exit is acceptable because the audit intentionally limits the subsequent DG solve to one inner and one outer iteration.",
        },
        "audit": audit_data,
        "paths": {
            "config": str(config_path),
            "hybrid_restart": str(HYBRID_RESTART),
            "console": str(console),
            "figure_png": str(figure),
            "figure_svg": str(figure.with_suffix(".svg")),
        },
        "hashes": {
            "config": sha256(config_path),
            "hybrid_restart": sha256(HYBRID_RESTART),
            "nodes": sha256(required[0]),
            "cells": sha256(required[1]),
            "regions": sha256(required[2]),
        },
    }
    markdown = render_markdown(report)
    local_json = output_dir / "fixed_state_residual_audit_summary.json"
    local_md = output_dir / "fixed_state_residual_audit_summary.md"
    local_json.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    local_md.write_text(markdown, encoding="utf-8")
    report_json = args.report_json.resolve()
    report_md = args.report_md.resolve()
    report_json.parent.mkdir(parents=True, exist_ok=True)
    report_md.parent.mkdir(parents=True, exist_ok=True)
    report_json.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    report_md.write_text(markdown, encoding="utf-8")
    print(json.dumps({"status": report["status"], "top_node": audit_data["top_nodes"][0], "top_region": audit_data["regions"][0]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
