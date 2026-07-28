#!/usr/bin/env python3
"""Compare avalanche-off terminal currents with the integrated SRH source."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

import diagnose_pn2d_bv_sg_avalanche_edges as vtk_tools


Q = 1.602176634e-19


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--vtk", type=Path, required=True)
    parser.add_argument("--terminal-balance", type=Path, required=True)
    parser.add_argument("--bias", type=float, required=True)
    parser.add_argument("--out", type=Path, required=True)
    return parser.parse_args()


def nodal_areas_um2(mesh: dict[str, Any]) -> tuple[list[float], float]:
    nodes = {
        int(node["id"]): (float(node["x"]), float(node["y"]))
        for node in mesh["nodes"]
    }
    areas = [0.0] * len(nodes)
    total = 0.0
    for triangle in mesh["triangles"]:
        ids = [int(value) for value in triangle["node_ids"]]
        a, b, c = (nodes[node_id] for node_id in ids)
        area = 0.5 * abs(
            (b[0] - a[0]) * (c[1] - a[1])
            - (b[1] - a[1]) * (c[0] - a[0])
        )
        total += area
        for node_id in ids:
            areas[node_id] += area / 3.0
    return areas, total


def terminal_rows(path: Path, bias: float) -> dict[str, dict[str, str]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    selected: dict[str, dict[str, str]] = {}
    for row in rows:
        row_bias = float(row["bias_V"])
        if math.isclose(row_bias, bias, rel_tol=0.0, abs_tol=1.0e-9):
            selected[row["contact"]] = row
    return selected


def main() -> None:
    args = parse_args()
    mesh = json.loads(args.mesh.read_text())
    areas, device_area = nodal_areas_um2(mesh)
    scalars = vtk_tools.parse_vtk_scalars(args.vtk)
    srh = scalars["SRHRecombination"]
    if len(srh) != len(areas):
        raise RuntimeError("SRHRecombination length does not match mesh nodes")

    # Unit-scaling VTK rates are in cm^-3 s^-1.  Multiplication by a nodal
    # area in um^2 and a 1-um device depth gives a volume factor of 1e-12 cm^3.
    integrated_generation_cm3_s_um2 = sum(
        -rate * area for rate, area in zip(srh, areas)
    )
    contact_nodes = {
        int(node_id)
        for contact in mesh.get("contacts", [])
        for node_id in contact.get("node_ids", [])
    }
    contact_generation_cm3_s_um2 = sum(
        -rate * areas[node_id]
        for node_id, rate in enumerate(srh)
        if node_id in contact_nodes
    )
    free_generation_cm3_s_um2 = (
        integrated_generation_cm3_s_um2 - contact_generation_cm3_s_um2
    )
    source_current_a_per_um = Q * integrated_generation_cm3_s_um2 * 1.0e-12
    contact_source_current_a_per_um = (
        Q * contact_generation_cm3_s_um2 * 1.0e-12
    )
    free_source_current_a_per_um = Q * free_generation_cm3_s_um2 * 1.0e-12

    contacts = terminal_rows(args.terminal_balance, args.bias)
    terminal = {
        name: {
            "electron": float(row["current_electron_A_per_um"]),
            "hole": float(row["current_hole_A_per_um"]),
            "total": float(row["current_total_A_per_um"]),
        }
        for name, row in contacts.items()
    }
    electron_outward_sum = sum(item["electron"] for item in terminal.values())
    hole_outward_sum = sum(item["hole"] for item in terminal.values())
    total_outward_sum = sum(item["total"] for item in terminal.values())
    # ContactCurrent reports both carrier particle-flux contributions with the
    # outward-generation convention here, so each carrier sum must equal the
    # positive integrated generation current.
    electron_closure_error = electron_outward_sum - source_current_a_per_um
    hole_closure_error = hole_outward_sum - source_current_a_per_um
    summary = {
        "bias_V": args.bias,
        "device_area_um2": device_area,
        "srh_rate_min_cm3_s": min(srh),
        "srh_rate_max_cm3_s": max(srh),
        "integrated_srh_generation_cm_minus1_s_minus1": (
            integrated_generation_cm3_s_um2 * 1.0e-8
        ),
        "integrated_srh_source_current_A_per_um": source_current_a_per_um,
        "contact_node_srh_source_current_A_per_um": (
            contact_source_current_a_per_um
        ),
        "free_node_srh_source_current_A_per_um": free_source_current_a_per_um,
        "contact_node_source_fraction": (
            contact_source_current_a_per_um / source_current_a_per_um
            if source_current_a_per_um != 0.0
            else None
        ),
        "terminal_current_A_per_um": terminal,
        "electron_outward_sum_A_per_um": electron_outward_sum,
        "hole_outward_sum_A_per_um": hole_outward_sum,
        "total_outward_sum_A_per_um": total_outward_sum,
        "electron_continuity_closure_error_A_per_um": electron_closure_error,
        "hole_continuity_closure_error_A_per_um": hole_closure_error,
        "electron_closure_error_over_source": (
            electron_closure_error / source_current_a_per_um
            if source_current_a_per_um != 0.0
            else None
        ),
        "hole_closure_error_over_source": (
            hole_closure_error / source_current_a_per_um
            if source_current_a_per_um != 0.0
            else None
        ),
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
