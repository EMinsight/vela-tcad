#!/usr/bin/env python3
"""Build same-cell PN2D high-bias process-chain evidence from accepted runs."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from typing import Any

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.analyze_pn2d_high_bias_oracle import records
from scripts.pn2d_high_bias_process_contract import EXACT_HIGH_BIAS_V


ACTIVE_FRACTION = 1.0e-3
MATERIAL_RATIO = 1.10
STAGE_METRICS = (
    ("density", ("mean_n_cm3", "mean_p_cm3")),
    ("drive", ("efield_V_cm", "grad_qf_n_V_cm", "grad_qf_p_V_cm")),
    ("mobility", ("mu_n_cm2_Vs", "mu_p_cm2_Vs")),
    ("current", ("current_n_A_cm2", "current_p_A_cm2")),
    ("alpha", ("max_alpha_n_cm_inv", "max_alpha_p_cm_inv")),
    ("generation", ("max_generation_total_cm3_s",)),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root-a", type=Path, required=True)
    parser.add_argument("--root-b", type=Path, required=True)
    parser.add_argument("--mesh-json", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    return parser.parse_args()


def vector_norm(row: dict[str, Any], x: str, y: str) -> float:
    return math.hypot(float(row[x]), float(row[y]))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty CSV {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def implicit_records(root: Path) -> list[dict[str, Any]]:
    path = root / "implicit_default" / "fetched" / "run_implicit_default.out"
    return records(path)


def group_records(parsed: list[dict[str, Any]]) -> dict[tuple[str, float], list[dict[str, Any]]]:
    grouped: dict[tuple[str, float], list[dict[str, Any]]] = {}
    for row in parsed:
        if "bias_V" not in row:
            continue
        grouped.setdefault((str(row["kind"]), float(row["bias_V"])), []).append(row)
    return grouped


def same_support_rows(mesh: dict[str, Any], parsed: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grouped = group_records(parsed)
    contact_nodes = {
        int(node)
        for contact in mesh["contacts"]
        for node in contact["node_ids"]
    }
    rows: list[dict[str, Any]] = []
    for bias in EXACT_HIGH_BIAS_V:
        vertices = {
            int(row["vertex"]): row
            for row in grouped[("vertex", bias)]
            if int(row["vertex"]) < len(mesh["nodes"])
        }
        processes = {
            int(row["vertex"]): row
            for row in grouped[("process", bias)]
            if int(row["vertex"]) < len(mesh["nodes"])
        }
        elements = {
            int(row["element"]): row
            for row in grouped[("element", bias)]
        }
        if len(vertices) != len(mesh["nodes"]):
            raise ValueError(f"{bias:g}: physical vertex mapping is incomplete")
        if len(elements) != len(mesh["triangles"]):
            raise ValueError(f"{bias:g}: element mapping is incomplete")
        bias_rows: list[dict[str, Any]] = []
        for triangle in mesh["triangles"]:
            cell = int(triangle["id"])
            node_ids = [int(node) for node in triangle["node_ids"]]
            vertex_rows = [vertices[node] for node in node_ids]
            process_rows = [processes[node] for node in node_ids]
            element = elements[cell]
            coords = [mesh["nodes"][node] for node in node_ids]
            bias_rows.append(
                {
                    "bias_V": bias,
                    "cell_id": cell,
                    "node_ids": ";".join(str(node) for node in node_ids),
                    "centroid_x_um": sum(float(node["x"]) for node in coords) / 3.0,
                    "centroid_y_um": sum(float(node["y"]) for node in coords) / 3.0,
                    "contact_class": (
                        "contact_adjacent"
                        if any(node in contact_nodes for node in node_ids)
                        else "interior"
                    ),
                    "mean_n_cm3": sum(float(row["n_cm3"]) for row in vertex_rows) / 3.0,
                    "mean_p_cm3": sum(float(row["p_cm3"]) for row in vertex_rows) / 3.0,
                    "efield_V_cm": vector_norm(element, "efield_x_V_cm", "efield_y_V_cm"),
                    "grad_qf_n_V_cm": vector_norm(
                        element, "grad_qf_n_x_V_cm", "grad_qf_n_y_V_cm"
                    ),
                    "grad_qf_p_V_cm": vector_norm(
                        element, "grad_qf_p_x_V_cm", "grad_qf_p_y_V_cm"
                    ),
                    "mu_n_cm2_Vs": abs(float(element["mu_n_cm2_Vs"])),
                    "mu_p_cm2_Vs": abs(float(element["mu_p_cm2_Vs"])),
                    "current_n_A_cm2": vector_norm(
                        element, "current_n_x_A_cm2", "current_n_y_A_cm2"
                    ),
                    "current_p_A_cm2": vector_norm(
                        element, "current_p_x_A_cm2", "current_p_y_A_cm2"
                    ),
                    "mean_velocity_n_cm_s": sum(
                        abs(float(row["velocity_n_cm_s"])) for row in process_rows
                    ) / 3.0,
                    "mean_velocity_p_cm_s": sum(
                        abs(float(row["velocity_p_cm_s"])) for row in process_rows
                    ) / 3.0,
                    "max_alpha_n_cm_inv": max(
                        abs(float(row["alpha_n_cm_inv"])) for row in vertex_rows
                    ),
                    "max_alpha_p_cm_inv": max(
                        abs(float(row["alpha_p_cm_inv"])) for row in vertex_rows
                    ),
                    "max_generation_total_cm3_s": max(
                        abs(float(row["generation_total_cm3_s"]))
                        for row in vertex_rows
                    ),
                }
            )
        peak = max(row["max_generation_total_cm3_s"] for row in bias_rows)
        threshold = peak * ACTIVE_FRACTION
        for row in bias_rows:
            row["active_region"] = int(
                peak > 0.0 and row["max_generation_total_cm3_s"] >= threshold
            )
            row["active_threshold_cm3_s"] = threshold
        rows.extend(bias_rows)
    return rows


def ratio(after: float, before: float) -> float | None:
    return abs(after) / abs(before) if before != 0.0 else None


def support_class_summary(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    for bias in EXACT_HIGH_BIAS_V:
        bias_rows = [row for row in rows if float(row["bias_V"]) == bias]
        for contact_class in ("contact_adjacent", "interior"):
            for active in (0, 1):
                group = [
                    row for row in bias_rows
                    if row["contact_class"] == contact_class
                    and int(row["active_region"]) == active
                ]
                if not group:
                    continue
                result.append(
                    {
                        "bias_V": bias,
                        "contact_class": contact_class,
                        "active_region": active,
                        "cell_count": len(group),
                        "max_current_n_A_cm2": max(
                            float(row["current_n_A_cm2"]) for row in group
                        ),
                        "max_current_p_A_cm2": max(
                            float(row["current_p_A_cm2"]) for row in group
                        ),
                        "max_alpha_n_cm_inv": max(
                            float(row["max_alpha_n_cm_inv"]) for row in group
                        ),
                        "max_alpha_p_cm_inv": max(
                            float(row["max_alpha_p_cm_inv"]) for row in group
                        ),
                        "max_generation_total_cm3_s": max(
                            float(row["max_generation_total_cm3_s"])
                            for row in group
                        ),
                    }
                )
    return result

def hotspot_evidence(rows: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    at_minus20 = [row for row in rows if float(row["bias_V"]) == -20.0]
    hotspot = max(at_minus20, key=lambda row: float(row["max_generation_total_cm3_s"]))
    cell_id = int(hotspot["cell_id"])
    chain = sorted(
        (row for row in rows if int(row["cell_id"]) == cell_id),
        key=lambda row: float(row["bias_V"]),
        reverse=True,
    )
    by_bias = {float(row["bias_V"]): row for row in chain}
    before = by_bias[-19.0]
    after = by_bias[-20.0]
    metric_ratios: dict[str, float | None] = {}
    stage_material: dict[str, bool] = {}
    first_stage: str | None = None
    for stage, metrics in STAGE_METRICS:
        values = [
            ratio(float(after[metric]), float(before[metric]))
            for metric in metrics
        ]
        material = any(value is not None and value >= MATERIAL_RATIO for value in values)
        stage_material[stage] = material
        if first_stage is None and material:
            first_stage = stage
        for metric, value in zip(metrics, values):
            metric_ratios[metric] = value
    summary = {
        "schema": "vela.pn2d.same_support_process_chain.v1",
        "hotspot_cell_id": cell_id,
        "hotspot_node_ids": hotspot["node_ids"],
        "hotspot_contact_class": hotspot["contact_class"],
        "active_fraction": ACTIVE_FRACTION,
        "material_ratio_threshold": MATERIAL_RATIO,
        "comparison_biases_V": [-19.0, -20.0],
        "first_material_stage": first_stage,
        "stage_material": stage_material,
        "metric_ratios_m20_over_m19": metric_ratios,
        "support_claim": (
            "native element mobility/field/current with native vertex quantities "
            "aggregated only over the same Tri3 nodes"
        ),
    }
    return chain, summary


def main() -> int:
    args = parse_args()
    root_a = args.root_a.resolve()
    root_b = args.root_b.resolve()
    parsed_a = implicit_records(root_a)
    parsed_b = implicit_records(root_b)
    normalized_a = [row for row in parsed_a if row["kind"] != "begin"]
    normalized_b = [row for row in parsed_b if row["kind"] != "begin"]
    if normalized_a != normalized_b:
        raise ValueError("paired implicit-default runtime records differ")
    mesh = json.loads(args.mesh_json.resolve().read_text(encoding="utf-8"))
    rows = same_support_rows(mesh, parsed_a)
    chain, summary = hotspot_evidence(rows)
    output = args.output_root.resolve()
    write_csv(output / "same_cell_process.csv", rows)
    write_csv(output / "support_class_summary.csv", support_class_summary(rows))
    write_csv(output / "fixed_hotspot_process_chain.csv", chain)
    output.mkdir(parents=True, exist_ok=True)
    (output / "same_support_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="ascii"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
