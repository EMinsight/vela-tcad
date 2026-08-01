#!/usr/bin/env python3
"""Audit the complete M2 carrier linear solve on frozen baseline/joint-QFP states."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
import sys
from collections import deque
from pathlib import Path
from typing import Any

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.diagnose_pn2d_bv_predictor_first_step_audit import make_probe_config
from scripts.run_pn2d_bv_m2_single_family_state_substitution import bias_tag


BIASES = (-18.0, -19.5, -19.7, -20.0)
VARIANTS = ("vela_baseline", "sent_qfp_only")
OUTPUT_SUFFIXES = (
    "_summary.json", "_columns.csv", "_singular_modes.csv",
    "_solve_variants.csv", "_solve_nodes.csv",
)
ROW_SCALED_STEP_TOLERANCE = 1.0e-8
LINEAR_CLOSURE_TOLERANCE = 1.0e-8
TRANSPORT_CROSS_RELATIVE_TOLERANCE = 1.0e-12
SVD_ENERGY_CLOSURE_TOLERANCE = 1.0e-10
MODAL_COMPONENT_CLOSURE_TOLERANCE = 1.0e-10


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise RuntimeError(f"refusing to write empty CSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def sha256(paths: list[Path]) -> str:
    digest = hashlib.sha256()
    for path in paths:
        digest.update(path.name.encode())
        digest.update(path.read_bytes())
    return digest.hexdigest()


def run_probe(runner: Path, config_path: Path) -> dict[str, Any]:
    result = subprocess.run(
        [str(runner.resolve()), "--config", str(config_path.resolve())],
        cwd=config_path.parent,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    config_path.with_suffix(".stdout.log").write_text(result.stdout, encoding="utf-8")
    config_path.with_suffix(".stderr.log").write_text(result.stderr, encoding="utf-8")
    if result.returncode:
        raise RuntimeError(
            f"runner failed for {config_path}: "
            f"{result.stderr.strip() or result.stdout.strip()}"
        )
    return json.loads(result.stdout)


def f(row: dict[str, str], field: str) -> float:
    return float(row[field])


def load_junction_support(base_config: Path) -> tuple[float, dict[int, dict[str, object]]]:
    config = json.loads(base_config.read_text(encoding="utf-8-sig"))
    mesh_path = Path(config["mesh_file"])
    doping_path = Path(config["node_doping_file"])
    if not mesh_path.is_absolute():
        mesh_path = (base_config.parent / mesh_path).resolve()
    if not doping_path.is_absolute():
        doping_path = (base_config.parent / doping_path).resolve()
    mesh = json.loads(mesh_path.read_text(encoding="utf-8-sig"))
    doping_rows = read_rows(doping_path)
    net = {
        int(row["node_id"]): float(row["donors_cm3"]) - float(row["acceptors_cm3"])
        for row in doping_rows
    }
    doping_scale = max(abs(value) for value in net.values())
    compensated = {
        node_id for node_id, value in net.items()
        if abs(value) <= max(1.0, doping_scale * 1.0e-12)
    }
    if not compensated:
        raise RuntimeError("M2 junction support has no compensated nodes")
    nodes = {int(node["id"]): node for node in mesh["nodes"]}
    junction_x = sum(float(nodes[node_id]["x"]) for node_id in compensated) / len(compensated)
    graph: dict[int, set[int]] = {node_id: set() for node_id in nodes}
    touches_compensated_triangle: dict[int, bool] = {node_id: False for node_id in nodes}
    for triangle in mesh["triangles"]:
        ids = [int(value) for value in triangle["node_ids"]]
        has_compensated = any(node_id in compensated for node_id in ids)
        for node_id in ids:
            graph[node_id].update(other for other in ids if other != node_id)
            touches_compensated_triangle[node_id] |= has_compensated
    distance = {node_id: -1 for node_id in nodes}
    queue: deque[int] = deque()
    for node_id in compensated:
        distance[node_id] = 0
        queue.append(node_id)
    while queue:
        node_id = queue.popleft()
        for neighbor in graph[node_id]:
            if distance[neighbor] < 0:
                distance[neighbor] = distance[node_id] + 1
                queue.append(neighbor)
    support = {}
    for node_id, node in nodes.items():
        support[node_id] = {
            "node_x_um": float(node["x"]),
            "node_y_um": float(node["y"]),
            "net_doping_cm3": net[node_id],
            "distance_to_junction_um": abs(float(node["x"]) - junction_x),
            "junction_graph_distance": distance[node_id],
            "touches_compensated_junction_triangle": int(
                touches_compensated_triangle[node_id]
            ),
        }
    return junction_x, support


def summarize_case(
    bias: float,
    variant: str,
    prefix: Path,
    junction_support: dict[int, dict[str, object]],
) -> tuple[dict[str, object], list[dict[str, object]], list[dict[str, object]]]:
    summary = json.loads(Path(str(prefix) + "_summary.json").read_text(encoding="utf-8"))
    columns = read_rows(Path(str(prefix) + "_columns.csv"))
    modes = read_rows(Path(str(prefix) + "_singular_modes.csv"))
    solves = read_rows(Path(str(prefix) + "_solve_variants.csv"))
    by_solve = {row["variant"]: row for row in solves}
    full = by_solve["full"]
    row_scaled = by_solve["row_scaled"]
    total_norm = (
        float(summary["electron_electron_norm"]) ** 2
        + float(summary["electron_hole_norm"]) ** 2
        + float(summary["hole_electron_norm"]) ** 2
        + float(summary["hole_hole_norm"]) ** 2
    ) ** 0.5
    total_cross_norm = (
        float(summary["electron_hole_norm"]) ** 2
        + float(summary["hole_electron_norm"]) ** 2
    ) ** 0.5
    case = {
        "bias_V": bias,
        "variant": variant,
        "raw_condition_number": summary["raw_condition"]["resolved_condition_number"],
        "row_scaled_condition_number": summary["row_scaled_condition"]["resolved_condition_number"],
        "l2_equilibrated_condition_number": summary["l2_equilibrated_condition"]["resolved_condition_number"],
        "raw_numerical_rank": summary["raw_condition"]["numerical_rank"],
        "row_scaled_numerical_rank": summary["row_scaled_condition"]["numerical_rank"],
        "l2_equilibrated_numerical_rank": summary["l2_equilibrated_condition"]["numerical_rank"],
        "free_unknown_count": int(summary["free_electron_unknowns"]) + int(summary["free_hole_unknowns"]),
        "phin_residual_norm": summary["block_residuals"]["phin"],
        "phip_residual_norm": summary["block_residuals"]["phip"],
        "electron_electron_norm": summary["electron_electron_norm"],
        "electron_hole_norm": summary["electron_hole_norm"],
        "hole_electron_norm": summary["hole_electron_norm"],
        "hole_hole_norm": summary["hole_hole_norm"],
        "free_column_norm_spread": summary["free_column_norm_spread"],
        "free_row_norm_spread": summary["free_row_norm_spread"],
        "row_weight_spread": summary["row_weight_spread"],
        "cross_carrier_norm_fraction": summary["cross_carrier_norm_fraction"],
        "recombination_cross_norm": summary["recombination_cross_norm"],
        "avalanche_cross_norm": summary["avalanche_cross_norm"],
        "avalanche_fraction_of_cross_norm": float(summary["avalanche_cross_norm"]) / max(total_cross_norm, 1.0e-300),
        "transport_cross_norm": summary["transport_cross_norm"],
        "transport_cross_relative_to_full": float(summary["transport_cross_norm"]) / max(total_norm, 1.0e-300),
        "full_physical_step_norm_V": f(full, "physical_step_norm_V"),
        "full_relative_linear_closure": f(full, "relative_linear_closure"),
        "row_scaled_relative_difference_from_full": f(row_scaled, "relative_difference_from_full"),
        "row_scaled_cosine_with_full": f(row_scaled, "cosine_with_full"),
        "singular_rhs_energy_sum": sum(f(row, "rhs_energy_fraction") for row in modes),
        "singular_step_energy_sum": sum(f(row, "step_energy_fraction") for row in modes),
    }
    for name in ("no_cross_carrier", "no_recombination", "no_avalanche", "transport_only"):
        row = by_solve[name]
        case[f"{name}_relative_difference_from_full"] = f(
            row, "relative_difference_from_full"
        )
        case[f"{name}_cosine_with_full"] = f(row, "cosine_with_full")

    dominant_columns: list[dict[str, object]] = []
    for ranking, field in (("column_norm", "column_l2_norm"), ("update", "full_delta_qfp_V")):
        ranked = sorted(columns, key=lambda row: abs(f(row, field)), reverse=True)
        for rank, row in enumerate(ranked[:10], start=1):
            node_id = int(row["node_id"])
            dominant_columns.append({
                "bias_V": bias,
                "variant": variant,
                "ranking": ranking,
                "rank": rank,
                "carrier": row["carrier"],
                "node_id": node_id,
                "x": f(row, "x"),
                "y": f(row, "y"),
                "column_l2_norm": f(row, "column_l2_norm"),
                "diagonal_fraction": f(row, "diagonal_fraction"),
                "cross_carrier_row_fraction": f(row, "cross_carrier_row_fraction"),
                "continuity_row_weight": f(row, "continuity_row_weight"),
                "residual": f(row, "residual"),
                "full_delta_qfp_V": f(row, "full_delta_qfp_V"),
                **junction_support[node_id],
            })

    dominant_modes: list[dict[str, object]] = []
    ranked_modes = sorted(modes, key=lambda row: f(row, "step_energy_fraction"), reverse=True)
    for rank, row in enumerate(ranked_modes[:12], start=1):
        right_node = int(row["top_right_node"])
        left_node = int(row["top_left_node"])
        dominant_modes.append({
            "bias_V": bias,
            "variant": variant,
            "step_energy_rank": rank,
            "mode_index": int(row["mode_index"]),
            "singular_value": f(row, "singular_value"),
            "relative_singular_value": f(row, "relative_singular_value"),
            "rhs_projection": f(row, "rhs_projection"),
            "rhs_energy_fraction": f(row, "rhs_energy_fraction"),
            "step_amplitude": f(row, "step_amplitude"),
            "step_energy_fraction": f(row, "step_energy_fraction"),
            "transport_jacobian_projection": f(row, "transport_jacobian_projection"),
            "recombination_jacobian_projection": f(row, "recombination_jacobian_projection"),
            "avalanche_diagonal_jacobian_projection": f(row, "avalanche_diagonal_jacobian_projection"),
            "avalanche_cross_jacobian_projection": f(row, "avalanche_cross_jacobian_projection"),
            "jacobian_projection_closure": f(row, "jacobian_projection_closure"),
            "transport_rhs_projection": f(row, "transport_rhs_projection"),
            "recombination_rhs_projection": f(row, "recombination_rhs_projection"),
            "avalanche_rhs_projection": f(row, "avalanche_rhs_projection"),
            "rhs_projection_closure": f(row, "rhs_projection_closure"),
            "no_cross_carrier_step_amplitude": f(row, "no_cross_carrier_step_amplitude"),
            "no_recombination_step_amplitude": f(row, "no_recombination_step_amplitude"),
            "no_avalanche_step_amplitude": f(row, "no_avalanche_step_amplitude"),
            "transport_only_step_amplitude": f(row, "transport_only_step_amplitude"),
            "right_electron_fraction": f(row, "right_electron_fraction"),
            "left_electron_fraction": f(row, "left_electron_fraction"),
            "top_right_carrier": row["top_right_carrier"],
            "top_right_node": right_node,
            "top_left_carrier": row["top_left_carrier"],
            "top_left_node": left_node,
            "top_right_distance_to_junction_um": junction_support[right_node]["distance_to_junction_um"],
            "top_right_junction_graph_distance": junction_support[right_node]["junction_graph_distance"],
            "top_right_touches_compensated_triangle": junction_support[right_node]["touches_compensated_junction_triangle"],
            "top_right_net_doping_cm3": junction_support[right_node]["net_doping_cm3"],
            "top_left_distance_to_junction_um": junction_support[left_node]["distance_to_junction_um"],
            "top_left_junction_graph_distance": junction_support[left_node]["junction_graph_distance"],
            "top_left_touches_compensated_triangle": junction_support[left_node]["touches_compensated_junction_triangle"],
            "top_left_net_doping_cm3": junction_support[left_node]["net_doping_cm3"],
        })
    return case, dominant_columns, dominant_modes


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--base-config", type=Path, required=True)
    parser.add_argument("--prior-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--repeats", type=int, default=2)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.repeats < 2:
        raise RuntimeError("at least two independent repeats are required")
    args.output_root.mkdir(parents=True, exist_ok=True)
    junction_x_um, junction_support = load_junction_support(args.base_config)
    hashes: dict[str, dict[str, str]] = {}
    first_prefixes: dict[tuple[float, str], Path] = {}
    for repeat in range(args.repeats):
        label = chr(ord("a") + repeat)
        hashes[label] = {}
        for bias in BIASES:
            tag = bias_tag(bias)
            for variant in VARIANTS:
                case = args.output_root / f"run-{label}" / tag / variant
                prefix = case / "carrier_block"
                config_path = case / "carrier_block.json"
                config = make_probe_config(
                    args.base_config,
                    case / "unused.csv",
                    args.prior_root / "inputs" / tag / f"{variant}_fields",
                    "newton_carrier_block_decomposition_probe",
                    bias,
                    "Anode",
                    "Cathode",
                )
                config.pop("output_csv", None)
                config["output_prefix"] = str(prefix.resolve())
                write_json(config_path, config)
                write_json(case / "status.json", run_probe(args.runner, config_path))
                outputs = [Path(str(prefix) + suffix) for suffix in OUTPUT_SUFFIXES]
                hashes[label][f"{tag}/{variant}"] = sha256(outputs)
                if repeat == 0:
                    first_prefixes[(bias, variant)] = prefix

    deterministic = all(
        len({hashes[label][key] for label in hashes}) == 1
        for key in hashes["a"]
    )
    cases: list[dict[str, object]] = []
    columns: list[dict[str, object]] = []
    modes: list[dict[str, object]] = []
    for bias in BIASES:
        for variant in VARIANTS:
            case, case_columns, case_modes = summarize_case(
                bias, variant, first_prefixes[(bias, variant)], junction_support
            )
            cases.append(case)
            columns.extend(case_columns)
            modes.extend(case_modes)

    maximum_row_scaled_difference = max(
        float(row["row_scaled_relative_difference_from_full"]) for row in cases
    )
    maximum_linear_closure = max(
        float(row["full_relative_linear_closure"]) for row in cases
    )
    maximum_transport_cross_relative = max(
        float(row["transport_cross_relative_to_full"]) for row in cases
    )
    maximum_svd_energy_closure = max(
        max(abs(float(row["singular_rhs_energy_sum"]) - 1.0),
            abs(float(row["singular_step_energy_sum"]) - 1.0))
        for row in cases
    )
    modal_jacobian_closures: list[float] = []
    modal_rhs_closures: list[float] = []
    for prefix in first_prefixes.values():
        for row in read_rows(Path(str(prefix) + "_singular_modes.csv")):
            jacobian_scale = max(
                abs(f(row, "singular_value")),
                sum(abs(f(row, field)) for field in (
                    "transport_jacobian_projection",
                    "recombination_jacobian_projection",
                    "avalanche_diagonal_jacobian_projection",
                    "avalanche_cross_jacobian_projection",
                )),
                1.0e-300,
            )
            rhs_scale = max(
                abs(f(row, "rhs_projection")),
                sum(abs(f(row, field)) for field in (
                    "transport_rhs_projection",
                    "recombination_rhs_projection",
                    "avalanche_rhs_projection",
                )),
                1.0e-300,
            )
            modal_jacobian_closures.append(
                abs(f(row, "jacobian_projection_closure")) / jacobian_scale
            )
            modal_rhs_closures.append(
                abs(f(row, "rhs_projection_closure")) / rhs_scale
            )
    maximum_modal_jacobian_closure = max(modal_jacobian_closures)
    maximum_modal_rhs_closure = max(modal_rhs_closures)
    passed = (
        deterministic
        and maximum_row_scaled_difference <= ROW_SCALED_STEP_TOLERANCE
        and maximum_linear_closure <= LINEAR_CLOSURE_TOLERANCE
        and maximum_transport_cross_relative <= TRANSPORT_CROSS_RELATIVE_TOLERANCE
        and maximum_svd_energy_closure <= SVD_ENERGY_CLOSURE_TOLERANCE
        and maximum_modal_jacobian_closure <= MODAL_COMPONENT_CLOSURE_TOLERANCE
        and maximum_modal_rhs_closure <= MODAL_COMPONENT_CLOSURE_TOLERANCE
    )
    determinism_rows = [{
        "artifact": key,
        "repeat_count": args.repeats,
        "unique_hash_count": len({hashes[label][key] for label in hashes}),
        "byte_identical": int(len({hashes[label][key] for label in hashes}) == 1),
        "sha256": hashes["a"][key],
    } for key in sorted(hashes["a"])]
    result = {
        "schema": "vela.pn2d_bv_m2_carrier_block_decomposition.v1",
        "status": "passed" if passed else "failed",
        "typed_outcome": "carrier_block_linear_solve_decomposed" if passed else "carrier_block_contract_failed",
        "biases_V": list(BIASES),
        "variants": list(VARIANTS),
        "junction_x_um": junction_x_um,
        "contract": {
            "condition_and_svd_domain": "free_carrier_qfp_rows_and_columns_only",
            "contact_identity_rows_excluded": True,
            "residual_frozen_across_counterfactual_solves": True,
            "production_defaults_modified": False,
            "row_scaled_step_tolerance": ROW_SCALED_STEP_TOLERANCE,
            "linear_closure_tolerance": LINEAR_CLOSURE_TOLERANCE,
            "transport_cross_relative_tolerance": TRANSPORT_CROSS_RELATIVE_TOLERANCE,
            "svd_energy_closure_tolerance": SVD_ENERGY_CLOSURE_TOLERANCE,
            "modal_component_closure_tolerance": MODAL_COMPONENT_CLOSURE_TOLERANCE,
        },
        "verdict": {
            "passed": passed,
            "deterministic": deterministic,
            "maximum_row_scaled_step_relative_difference": maximum_row_scaled_difference,
            "maximum_full_linear_closure": maximum_linear_closure,
            "maximum_transport_cross_relative_to_full": maximum_transport_cross_relative,
            "maximum_svd_energy_closure": maximum_svd_energy_closure,
            "maximum_modal_jacobian_projection_closure": maximum_modal_jacobian_closure,
            "maximum_modal_rhs_projection_closure": maximum_modal_rhs_closure,
        },
        "outputs": {
            "case_summary": "case_summary.csv",
            "dominant_columns": "dominant_columns.csv",
            "dominant_singular_modes": "dominant_singular_modes.csv",
            "determinism": "determinism.csv",
        },
    }
    write_rows(args.output_root / "case_summary.csv", cases)
    write_rows(args.output_root / "dominant_columns.csv", columns)
    write_rows(args.output_root / "dominant_singular_modes.csv", modes)
    write_rows(args.output_root / "determinism.csv", determinism_rows)
    write_json(args.output_root / "result.json", result)
    print(json.dumps(result, sort_keys=True))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
