#!/usr/bin/env python3
"""Run and summarize paired Vela/Sentaurus PN2D SRH mesh audits."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import subprocess
from pathlib import Path
from typing import Any, Sequence


REPO = Path(__file__).resolve().parents[1]
ANCHORS = (1, 5, 10, 15, 20)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_checked(command: Sequence[str]) -> None:
    completed = subprocess.run(list(command), cwd=REPO, check=False)
    if completed.returncode:
        raise RuntimeError(
            f"command failed ({completed.returncode}): {' '.join(command)}"
        )


def triangle_area(points: Sequence[tuple[float, float]]) -> float:
    (x0, y0), (x1, y1), (x2, y2) = points
    return abs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)) / 2.0


def mesh_metrics(mesh_path: Path, doping_path: Path) -> dict[str, Any]:
    mesh = json.loads(mesh_path.read_text(encoding="utf-8-sig"))
    nodes = {
        int(row["id"]): (float(row["x"]), float(row["y"]))
        for row in mesh["nodes"]
    }
    areas = {node_id: 0.0 for node_id in nodes}
    total_area = 0.0
    for cell in mesh["triangles"]:
        node_ids = tuple(int(value) for value in cell["node_ids"])
        cell_area = triangle_area(tuple(nodes[node_id] for node_id in node_ids))
        total_area += cell_area
        for node_id in node_ids:
            areas[node_id] += cell_area / 3.0
    doping = {
        int(row["node_id"]): (
            float(row["donors_cm3"]),
            float(row["acceptors_cm3"]),
        )
        for row in read_csv(doping_path)
    }
    dose = sum(
        areas[node_id] * (donors + acceptors)
        for node_id, (donors, acceptors) in doping.items()
    )
    xs = [point[0] for point in nodes.values()]
    ys = [point[1] for point in nodes.values()]
    contacts = {
        contact["name"]: sorted(
            (nodes[int(node_id)] for node_id in contact["node_ids"])
        )
        for contact in mesh["contacts"]
    }
    return {
        "mesh_sha256": sha256(mesh_path),
        "doping_sha256": sha256(doping_path),
        "node_count": len(nodes),
        "triangle_count": len(mesh["triangles"]),
        "bounds_um": [min(xs), min(ys), max(xs), max(ys)],
        "area_um2": total_area,
        "total_impurity_dose_cm3_um2": dose,
        "contacts": contacts,
    }


def summarize_level(report_dir: Path) -> dict[str, Any]:
    spatial = read_csv(report_dir / "spatial_summary.csv")
    profiles = read_csv(report_dir / "source_profile_metrics.csv")
    nonzero = [row for row in spatial if int(round(abs(float(row["bias_V"])))) > 0]
    anchor_rows = {
        int(round(abs(float(row["bias_V"])))): row
        for row in spatial
        if int(round(abs(float(row["bias_V"])))) in ANCHORS
    }
    source = {
        simulator: {
            int(round(float(row["reverse_bias_V"]))): abs(
                float(row["integrated_positive_generation_A_per_um"])
            )
            for row in profiles
            if row["simulator"] == simulator
            and int(round(float(row["reverse_bias_V"]))) in ANCHORS
        }
        for simulator in ("vela", "sentaurus")
    }
    errors = [float(row["log10_abs_current_ratio"]) for row in nonzero]
    return {
        "converged_points": sum(
            row["vela_converged"] == "1" and row["sentaurus_converged"] == "1"
            for row in spatial
        ),
        "log_current_rmse_dex": math.sqrt(
            sum(value * value for value in errors) / len(errors)
        ),
        "max_abs_log_current_error_dex": max(abs(value) for value in errors),
        "max_electron_closure_relative": max(
            abs(float(row["electron_closure_relative"])) for row in nonzero
        ),
        "max_hole_closure_relative": max(
            abs(float(row["hole_closure_relative"])) for row in nonzero
        ),
        "anchors": {
            str(bias): {
                "log_current_error_dex": float(
                    anchor_rows[bias]["log10_abs_current_ratio"]
                ),
                "vela_source_A_per_um": source["vela"][bias],
                "sentaurus_source_A_per_um": source["sentaurus"][bias],
            }
            for bias in ANCHORS
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix-root", type=Path, required=True)
    parser.add_argument("--base-config", type=Path, required=True)
    parser.add_argument("--out-root", type=Path, required=True)
    parser.add_argument("--levels", default="M0,M1,M2")
    parser.add_argument("--skip-runs", action="store_true")
    args = parser.parse_args()

    matrix_root = args.matrix_root.resolve()
    out_root = args.out_root.resolve()
    out_root.mkdir(parents=True, exist_ok=True)
    levels = [value.strip() for value in args.levels.split(",") if value.strip()]
    python = Path(os.sys.executable)
    level_rows: dict[str, Any] = {}

    for level in levels:
        imported = matrix_root / "imports" / f"{level}_vela"
        sent_source = (
            matrix_root
            / "sentaurus_runs"
            / f"pn2d_srh_coarse_mesh_{level}_20260728"
            / "source"
        )
        level_out = out_root / level
        level_out.mkdir(parents=True, exist_ok=True)
        curve = level_out / "sentaurus_curve.csv"
        config = json.loads(args.base_config.read_text(encoding="utf-8-sig"))
        config["mesh_file"] = str((imported / "mesh.json").resolve())
        config["node_doping_file"] = str((imported / "doping.csv").resolve())
        base_config = level_out / "base_config.json"
        base_config.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
        if not args.skip_runs:
            run_checked(
                [
                    str(python),
                    str(REPO / "scripts" / "extract_pn2d_bv_terminal_reference.py"),
                    "--sentaurus-plt",
                    str(sent_source / "pn2d_bv.plt"),
                    "--out-csv",
                    str(curve),
                ]
            )
            run_checked(
                [
                    str(python),
                    str(REPO / "scripts" / "run_pn2d_bv_off_srh_spatial_audit.py"),
                    "--base-config",
                    str(base_config),
                    "--sentaurus-source",
                    str(sent_source),
                    "--sentaurus-curve",
                    str(curve),
                    "--sentaurus-intervals",
                    "400",
                    "--out-dir",
                    str(level_out / "audit"),
                ]
            )
        level_rows[level] = {
            "mesh": mesh_metrics(imported / "mesh.json", imported / "doping.csv"),
            "audit": summarize_level(level_out / "audit" / "report"),
        }

    ordered = [level_rows[level] for level in levels]
    baseline_dose = ordered[0]["mesh"]["total_impurity_dose_cm3_um2"]
    for row in ordered:
        row["mesh"]["dose_relative_to_first"] = (
            row["mesh"]["total_impurity_dose_cm3_um2"] / baseline_dose - 1.0
        )
    finest_change: dict[str, dict[str, float]] = {}
    if "M1" in level_rows and "M2" in level_rows:
        for simulator in ("vela", "sentaurus"):
            finest_change[simulator] = {}
            for bias in ANCHORS:
                first = level_rows["M1"]["audit"]["anchors"][str(bias)][
                    f"{simulator}_source_A_per_um"
                ]
                second = level_rows["M2"]["audit"]["anchors"][str(bias)][
                    f"{simulator}_source_A_per_um"
                ]
                finest_change[simulator][str(bias)] = abs(second / first - 1.0)
    summary = {
        "schema": "vela.pn2d_bv_off_srh_mesh_matrix.v1",
        "levels": level_rows,
        "finest_source_relative_change": finest_change,
    }
    (out_root / "mesh_matrix_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
