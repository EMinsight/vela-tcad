#!/usr/bin/env python3
"""Compare Sentaurus Eq. 231 quadratic source with recovered vertex fields."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

from audit_singledevice_eq231_box_formula0 import mixed_voronoi_shares


def read_vector_field(export_dir: Path, field_name: str) -> dict[int, tuple[float, float]]:
    manifest = json.loads(
        (export_dir / "field_manifest.json").read_text(encoding="utf-8")
    )
    result: dict[int, tuple[float, float]] = {}
    for entry in manifest["fields"]:
        if entry["name"] != field_name or entry["mapping_status"] != "complete":
            continue
        with (export_dir / "fields" / entry["csv_file"]).open(
            newline="", encoding="utf-8"
        ) as stream:
            for row in csv.DictReader(stream):
                result[int(row["node_id"])] = (
                    float(row["component0"]), float(row["component1"])
                )
    if not result:
        raise ValueError(f"field {field_name!r} is missing from {export_dir}")
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--cells", type=Path, required=True)
    parser.add_argument("--state-export", type=Path, required=True)
    parser.add_argument("--phi-gradient-export", type=Path, required=True)
    parser.add_argument("--thermal-voltage", type=float, default=0.025851999786)
    parser.add_argument("--theta", type=float, default=0.5)
    parser.add_argument("--nodes", type=int, nargs="*", default=[])
    parser.add_argument("--top", type=int, default=30)
    args = parser.parse_args()

    mesh = json.loads(args.mesh.read_text(encoding="utf-8"))
    nodes = {int(item["id"]): item for item in mesh["nodes"]}
    triangles = {int(item["id"]): item for item in mesh["triangles"]}
    box_shares: dict[tuple[int, int], float] = {}
    for cell_id, triangle in triangles.items():
        triangle_nodes = triangle.get("nodes", triangle.get("node_ids"))
        points = [nodes[int(node)] for node in triangle_nodes]
        for local, share in enumerate(mixed_voronoi_shares(points)):
            box_shares[(cell_id, local)] = share

    terms = defaultdict(lambda: {"stiffness": 0.0, "reaction": 0.0, "volume": 0.0})
    metadata: dict[int, dict[str, str]] = {}
    with args.cells.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            if row["is_active"] != "1" or row["is_dirichlet"] == "1":
                continue
            node = int(row["node_id"])
            cell_id = int(row["cell_id"])
            local = int(row["local_node"])
            terms[node]["stiffness"] += float(row["stiffness"])
            terms[node]["reaction"] += float(row["reaction"])
            terms[node]["volume"] += box_shares[(cell_id, local)]
            metadata[node] = {
                "region_name": row["region_name"],
                "material": row["material"],
            }

    psi_field = read_vector_field(args.state_export, "ElectricField")
    phi_field = read_vector_field(args.phi_gradient_export, "ElectricField")
    rows = []
    for node, value in terms.items():
        if node not in psi_field or node not in phi_field:
            continue
        # Sentaurus ElectricField is -grad(field) in V/cm. In an insulator,
        # Eq. 231 uses A=-grad(Phi)-grad(psi), so A=(E_Phi+E_psi)/Vt.
        ax = 1.0e-4 * (phi_field[node][0] + psi_field[node][0]) / args.thermal_voltage
        ay = 1.0e-4 * (phi_field[node][1] + psi_field[node][1]) / args.thermal_voltage
        predicted = args.theta * (ax * ax + ay * ay) * value["volume"]
        required = -(value["stiffness"] + value["reaction"])
        rows.append({
            "node_id": node,
            **metadata[node],
            "stiffness": value["stiffness"],
            "reaction": value["reaction"],
            "required_quadratic": required,
            "recovered_ax_per_um": ax,
            "recovered_ay_per_um": ay,
            "recovered_quadratic": predicted,
            "quadratic_error": predicted - required,
            "absolute_error": abs(predicted - required),
        })
    rows.sort(key=lambda row: row["absolute_error"], reverse=True)
    selected = [row for row in rows if row["node_id"] in set(args.nodes)]
    print(json.dumps({
        "top": rows[: args.top],
        "selected": selected,
    }, indent=2))


if __name__ == "__main__":
    main()
