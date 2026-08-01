#!/usr/bin/env python3
"""Audit M2 nodal doping, topology, edge averages, and control-volume semantics."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path
from typing import Any

if not __package__:
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.diagnose_pn2d_bv_fixed_charge_rhs import (
    compute_mixed_voronoi_volumes,
    triangle_area_m2,
)


RELATIVE_TOLERANCE = 1.0e-12
COORDINATE_TOLERANCE_UM = 1.0e-12


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vela-config", type=Path, required=True)
    parser.add_argument("--sentaurus-export", type=Path, required=True)
    parser.add_argument("--soft-modes", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--junction-x-min", type=float, default=0.75)
    parser.add_argument("--junction-x-max", type=float, default=1.25)
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise RuntimeError(f"refusing to write empty CSV: {path}")
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, value: Any) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


def resolve_config_path(config_path: Path, raw: str) -> Path:
    path = Path(raw)
    return path.resolve() if path.is_absolute() else (config_path.parent / path).resolve()


def load_vela_mesh(path: Path) -> tuple[dict[int, dict[str, float]], list[dict[str, Any]]]:
    payload = json.loads(path.read_text(encoding="utf-8-sig"))
    nodes = {
        int(row["id"]): {"x_um": float(row["x"]), "y_um": float(row["y"])}
        for row in payload["nodes"]
    }
    triangles = [
        {"id": int(row["id"]), "node_ids": [int(value) for value in row["node_ids"]]}
        for row in payload["triangles"]
    ]
    return nodes, triangles


def load_sentaurus_mesh(root: Path) -> tuple[dict[int, dict[str, float]], list[dict[str, Any]]]:
    nodes = {
        int(row["id"]): {"x_um": float(row["x_um"]), "y_um": float(row["y_um"])}
        for row in read_csv(root / "nodes.csv")
    }
    triangles = [
        {
            "id": int(row["id"]),
            "node_ids": [int(row["node0"]), int(row["node1"]), int(row["node2"])],
        }
        for row in read_csv(root / "elements.csv")
    ]
    return nodes, triangles


def load_doping(path: Path) -> dict[int, dict[str, float]]:
    return {
        int(row["node_id"]): {
            "donors_cm3": float(row["donors_cm3"]),
            "acceptors_cm3": float(row["acceptors_cm3"]),
            "net_cm3": float(row["donors_cm3"]) - float(row["acceptors_cm3"]),
        }
        for row in read_csv(path)
    }


def load_scalar(path: Path) -> dict[int, float]:
    return {int(row["node_id"]): float(row["component0"]) for row in read_csv(path)}


def barycentric_volumes(
    nodes: dict[int, dict[str, float]], triangles: list[dict[str, Any]]
) -> dict[int, float]:
    volumes = {node_id: 0.0 for node_id in nodes}
    for triangle in triangles:
        points = [
            (nodes[node_id]["x_um"], nodes[node_id]["y_um"])
            for node_id in triangle["node_ids"]
        ]
        share = triangle_area_m2(points) / 3.0
        for node_id in triangle["node_ids"]:
            volumes[node_id] += share
    return volumes


def edges(triangles: list[dict[str, Any]]) -> list[tuple[int, int]]:
    result = set()
    for triangle in triangles:
        ids = triangle["node_ids"]
        for index in range(3):
            result.add(tuple(sorted((ids[index], ids[(index + 1) % 3]))))
    return sorted(result)


def relative_difference(left: float, right: float, scale: float) -> float:
    return abs(left - right) / max(abs(left), abs(right), scale, 1.0e-300)


def main() -> int:
    args = parse_args()
    args.output_root.mkdir(parents=True, exist_ok=True)
    config = json.loads(args.vela_config.read_text(encoding="utf-8-sig"))
    vela_mesh_path = resolve_config_path(args.vela_config, config["mesh_file"])
    vela_doping_path = resolve_config_path(args.vela_config, config["node_doping_file"])
    vela_nodes, vela_triangles = load_vela_mesh(vela_mesh_path)
    sent_nodes, sent_triangles = load_sentaurus_mesh(args.sentaurus_export)
    vela_doping = load_doping(vela_doping_path)
    sent_doping = load_doping(args.sentaurus_export / "doping.csv")

    ids_exact = set(vela_nodes) == set(sent_nodes) == set(vela_doping) == set(sent_doping)
    common_ids = sorted(set(vela_nodes) & set(sent_nodes) & set(vela_doping) & set(sent_doping))
    max_coordinate_error = max(
        max(
            abs(vela_nodes[node_id]["x_um"] - sent_nodes[node_id]["x_um"]),
            abs(vela_nodes[node_id]["y_um"] - sent_nodes[node_id]["y_um"]),
        )
        for node_id in common_ids
    )
    vela_connectivity = {row["id"]: tuple(sorted(row["node_ids"])) for row in vela_triangles}
    sent_connectivity = {row["id"]: tuple(sorted(row["node_ids"])) for row in sent_triangles}
    topology_exact = vela_connectivity == sent_connectivity

    doping_scale = max(
        max(abs(row["donors_cm3"]), abs(row["acceptors_cm3"]), abs(row["net_cm3"]))
        for row in vela_doping.values()
    )
    max_donor_error = max(
        relative_difference(
            vela_doping[node_id]["donors_cm3"], sent_doping[node_id]["donors_cm3"], doping_scale
        )
        for node_id in common_ids
    )
    max_acceptor_error = max(
        relative_difference(
            vela_doping[node_id]["acceptors_cm3"], sent_doping[node_id]["acceptors_cm3"], doping_scale
        )
        for node_id in common_ids
    )
    max_net_error = max(
        relative_difference(
            vela_doping[node_id]["net_cm3"], sent_doping[node_id]["net_cm3"], doping_scale
        )
        for node_id in common_ids
    )

    fields = args.sentaurus_export / "fields"
    sent_reported_net = load_scalar(fields / "DopingConcentration_region0.csv")
    sent_donor_field = load_scalar(fields / "PhosphorusActiveConcentration_region0.csv")
    sent_acceptor_field = load_scalar(fields / "BoronActiveConcentration_region0.csv")
    max_reported_net_error = max(
        relative_difference(sent_reported_net[node_id], sent_doping[node_id]["net_cm3"], doping_scale)
        for node_id in common_ids
    )
    max_donor_field_error = max(
        relative_difference(sent_donor_field[node_id], sent_doping[node_id]["donors_cm3"], doping_scale)
        for node_id in common_ids
    )
    max_acceptor_field_error = max(
        relative_difference(sent_acceptor_field[node_id], sent_doping[node_id]["acceptors_cm3"], doping_scale)
        for node_id in common_ids
    )

    barycentric = barycentric_volumes(vela_nodes, vela_triangles)
    mixed_values = compute_mixed_voronoi_volumes(vela_nodes, vela_triangles)
    mixed = {node_id: mixed_values[node_id] for node_id in vela_nodes}
    total_area = sum(barycentric.values())
    total_mixed_area = sum(mixed.values())
    area_closure = abs(total_area - total_mixed_area) / max(total_area, 1.0e-300)

    mode_rows = read_csv(args.soft_modes)
    mode_hits: Counter[int] = Counter()
    mode_energy: Counter[int] = Counter()
    for row in mode_rows:
        if row["variant"] != "sent_qfp_only" or int(row["step_energy_rank"]) > 2:
            continue
        for field in ("top_right_node", "top_left_node"):
            node_id = int(row[field])
            mode_hits[node_id] += 1
            mode_energy[node_id] += float(row["step_energy_fraction"])

    node_rows: list[dict[str, Any]] = []
    for node_id in common_ids:
        node = vela_nodes[node_id]
        vd = vela_doping[node_id]
        sd = sent_doping[node_id]
        bary = barycentric[node_id]
        mix = mixed[node_id]
        node_rows.append({
            "node_id": node_id,
            "x_um": node["x_um"],
            "y_um": node["y_um"],
            "vela_donors_cm3": vd["donors_cm3"],
            "sentaurus_donors_cm3": sd["donors_cm3"],
            "vela_acceptors_cm3": vd["acceptors_cm3"],
            "sentaurus_acceptors_cm3": sd["acceptors_cm3"],
            "vela_net_cm3": vd["net_cm3"],
            "sentaurus_reported_net_cm3": sent_reported_net[node_id],
            "barycentric_volume_m2": bary,
            "mixed_voronoi_volume_m2": mix,
            "mixed_to_barycentric_ratio": mix / bary if bary > 0.0 else math.nan,
            "barycentric_net_dopants_per_m": vd["net_cm3"] * 1.0e6 * bary,
            "mixed_net_dopants_per_m": vd["net_cm3"] * 1.0e6 * mix,
            "dominant_soft_mode_hit_count": mode_hits[node_id],
            "dominant_soft_mode_energy_sum": mode_energy[node_id],
        })

    junction_edges: list[dict[str, Any]] = []
    for n0, n1 in edges(vela_triangles):
        x0 = vela_nodes[n0]["x_um"]
        x1 = vela_nodes[n1]["x_um"]
        if max(x0, x1) < args.junction_x_min or min(x0, x1) > args.junction_x_max:
            continue
        vela_average = 0.5 * (vela_doping[n0]["net_cm3"] + vela_doping[n1]["net_cm3"])
        sent_average = 0.5 * (sent_doping[n0]["net_cm3"] + sent_doping[n1]["net_cm3"])
        junction_edges.append({
            "node0": n0,
            "node1": n1,
            "x0_um": x0,
            "x1_um": x1,
            "vela_endpoint_average_net_cm3": vela_average,
            "sentaurus_endpoint_average_net_cm3": sent_average,
            "relative_difference": relative_difference(vela_average, sent_average, doping_scale),
        })

    compensated = {
        node_id for node_id in common_ids
        if abs(vela_doping[node_id]["net_cm3"]) <= doping_scale * RELATIVE_TOLERANCE
    }
    junction_triangles: list[dict[str, Any]] = []
    for triangle in vela_triangles:
        ids = triangle["node_ids"]
        xs = [vela_nodes[node_id]["x_um"] for node_id in ids]
        if max(xs) < args.junction_x_min or min(xs) > args.junction_x_max:
            continue
        junction_triangles.append({
            "triangle_id": triangle["id"],
            "node0": ids[0],
            "node1": ids[1],
            "node2": ids[2],
            "sentaurus_unordered_connectivity_exact": int(
                sent_connectivity.get(triangle["id"]) == tuple(sorted(ids))
            ),
            "compensated_node_count": sum(node_id in compensated for node_id in ids),
            "node0_net_cm3": vela_doping[ids[0]]["net_cm3"],
            "node1_net_cm3": vela_doping[ids[1]]["net_cm3"],
            "node2_net_cm3": vela_doping[ids[2]]["net_cm3"],
        })

    mobility_basis = config["solver"]["mobility"]["doping_concentration_basis"]
    node_volume_policy = config.get("mesh_geometry", {}).get("node_volume_policy", "barycentric")
    maximum_edge_average_error = max(float(row["relative_difference"]) for row in junction_edges)
    soft_nodes = [row for row in node_rows if int(row["dominant_soft_mode_hit_count"]) > 0]
    soft_volume_ratios = [float(row["mixed_to_barycentric_ratio"]) for row in soft_nodes]
    passed = (
        ids_exact
        and topology_exact
        and max_coordinate_error <= COORDINATE_TOLERANCE_UM
        and max(max_donor_error, max_acceptor_error, max_net_error) <= RELATIVE_TOLERANCE
        and max(max_reported_net_error, max_donor_field_error, max_acceptor_field_error) <= RELATIVE_TOLERANCE
        and maximum_edge_average_error <= RELATIVE_TOLERANCE
        and area_closure <= RELATIVE_TOLERANCE
        and mobility_basis == "net_doping"
        and node_volume_policy == "barycentric"
    )
    result = {
        "schema": "vela.pn2d_bv_m2_doping_control_volume_semantics.v1",
        "status": "passed" if passed else "failed",
        "typed_outcome": (
            "nodal_doping_topology_and_edge_average_exact__sentaurus_control_volume_not_exported"
            if passed else "doping_or_geometry_semantics_mismatch"
        ),
        "contract": {
            "production_defaults_modified": False,
            "frozen_state_only": True,
            "sentaurus_is_golden": True,
            "relative_tolerance": RELATIVE_TOLERANCE,
            "coordinate_tolerance_um": COORDINATE_TOLERANCE_UM,
        },
        "configuration": {
            "mobility_doping_concentration_basis": mobility_basis,
            "node_volume_policy": node_volume_policy,
            "triangle_cell_reconstructed_total_impurity_active": False,
            "edge_mobility_doping_semantics": "arithmetic_mean_of_endpoint_nodal_net_doping",
        },
        "comparison": {
            "vela_node_count": len(vela_nodes),
            "sentaurus_node_count": len(sent_nodes),
            "vela_triangle_count": len(vela_triangles),
            "sentaurus_triangle_count": len(sent_triangles),
            "node_id_sets_exact": ids_exact,
            "unordered_triangle_connectivity_exact": topology_exact,
            "maximum_coordinate_error_um": max_coordinate_error,
            "maximum_donor_relative_error": max_donor_error,
            "maximum_acceptor_relative_error": max_acceptor_error,
            "maximum_net_relative_error": max_net_error,
            "maximum_sentaurus_reported_net_relative_error": max_reported_net_error,
            "maximum_sentaurus_donor_field_relative_error": max_donor_field_error,
            "maximum_sentaurus_acceptor_field_relative_error": max_acceptor_field_error,
            "maximum_junction_edge_average_relative_error": maximum_edge_average_error,
            "compensated_node_count": len(compensated),
        },
        "control_volume": {
            "vela_active_policy": "barycentric_area_over_three",
            "sentaurus_control_volume_field_present_in_export": False,
            "sentaurus_control_volume_semantics_directly_observable": False,
            "barycentric_total_area_m2": total_area,
            "mixed_voronoi_total_area_m2": total_mixed_area,
            "total_area_relative_closure": area_closure,
            "soft_mode_node_count": len(soft_nodes),
            "soft_mode_mixed_to_barycentric_ratio_min": min(soft_volume_ratios),
            "soft_mode_mixed_to_barycentric_ratio_max": max(soft_volume_ratios),
        },
        "inputs": {
            "vela_config": str(args.vela_config),
            "vela_config_sha256": sha256(args.vela_config),
            "vela_mesh": str(vela_mesh_path),
            "vela_mesh_sha256": sha256(vela_mesh_path),
            "vela_doping": str(vela_doping_path),
            "vela_doping_sha256": sha256(vela_doping_path),
            "sentaurus_nodes_sha256": sha256(args.sentaurus_export / "nodes.csv"),
            "sentaurus_elements_sha256": sha256(args.sentaurus_export / "elements.csv"),
            "sentaurus_doping_sha256": sha256(args.sentaurus_export / "doping.csv"),
            "soft_modes_sha256": sha256(args.soft_modes),
        },
        "outputs": {
            "node_comparison": "node_comparison.csv",
            "junction_edges": "junction_edges.csv",
            "junction_triangles": "junction_triangles.csv",
        },
    }
    write_csv(args.output_root / "node_comparison.csv", node_rows)
    write_csv(args.output_root / "junction_edges.csv", junction_edges)
    write_csv(args.output_root / "junction_triangles.csv", junction_triangles)
    write_json(args.output_root / "result.json", result)
    print(json.dumps(result, sort_keys=True))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
