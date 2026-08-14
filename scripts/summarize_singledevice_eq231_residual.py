#!/usr/bin/env python3
"""Summarize fixed-state SingleDevice Eq. 231 residual diagnostics."""

from __future__ import annotations

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def f(row: dict[str, str], key: str) -> float:
    return float(row[key])


def summarize(prefix: Path) -> dict[str, object]:
    cell_rows = rows(Path(f"{prefix}_cells.csv"))
    node_rows = rows(Path(f"{prefix}_nodes.csv"))
    region_rows = rows(Path(f"{prefix}_regions.csv"))
    free = [
        row for row in cell_rows
        if row["is_active"] == "1" and row["is_dirichlet"] == "0"
    ]
    term_l1 = {
        key: sum(abs(f(row, key)) for row in free)
        for key in ("stiffness", "gradient_squared", "reaction", "total")
    }
    component_sum = sum(term_l1[key] for key in (
        "stiffness", "gradient_squared", "reaction"))
    interface_l1 = sum(
        abs(f(row, "total")) for row in free
        if row["is_interface_cell"] == "1"
    )
    domain_l1 = {
        "transport": sum(
            abs(f(row, "total")) for row in free
            if row["is_transport"] == "1"),
        "insulator": sum(
            abs(f(row, "total")) for row in free
            if row["is_transport"] == "0"),
    }

    domain_node_residual: dict[str, dict[int, float]] = {
        "transport": defaultdict(float),
        "insulator": defaultdict(float),
    }
    for row in free:
        domain = "transport" if row["is_transport"] == "1" else "insulator"
        domain_node_residual[domain][int(row["node_id"])] += f(row, "total")
    domain_max = {}
    for domain, values in domain_node_residual.items():
        node_id, value = max(values.items(), key=lambda item: abs(item[1]))
        domain_max[domain] = {"node_id": node_id, "residual": value}

    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in free:
        grouped[row["cell_id"]].append(row)
    cells = []
    for cell_id, group in grouped.items():
        cells.append({
            "cell_id": int(cell_id),
            "region": group[0]["region_name"],
            "material": group[0]["material"],
            "interface_cell": group[0]["is_interface_cell"] == "1",
            "interface_pairs": group[0]["interface_pairs"],
            "transport": group[0]["is_transport"] == "1",
            "centroid_x_internal": sum(f(row, "x_internal") for row in group) / len(group),
            "centroid_y_internal": sum(f(row, "y_internal") for row in group) / len(group),
            "max_abs_local_total": max(abs(f(row, "total")) for row in group),
            "total_l1": sum(abs(f(row, "total")) for row in group),
            "stiffness_l1": sum(abs(f(row, "stiffness")) for row in group),
            "gradient_squared_l1": sum(abs(f(row, "gradient_squared")) for row in group),
            "reaction_l1": sum(abs(f(row, "reaction")) for row in group),
        })
    cells.sort(key=lambda row: row["max_abs_local_total"], reverse=True)

    free_nodes = [
        row for row in node_rows
        if row["is_active"] == "1" and row["is_dirichlet"] == "0"
    ]
    free_nodes.sort(key=lambda row: abs(f(row, "raw_total")), reverse=True)
    top_nodes = [{
        "node_id": int(row["node_id"]),
        "x_internal": f(row, "x_internal"),
        "y_internal": f(row, "y_internal"),
        "stiffness": f(row, "stiffness"),
        "gradient_squared": f(row, "gradient_squared"),
        "reaction": f(row, "reaction"),
        "raw_total": f(row, "raw_total"),
    } for row in free_nodes[:10]]

    regions = [{
        "region": row["region_name"],
        "material": row["material"],
        "cell_count": int(row["cell_count"]),
        "interface_cell_count": int(row["interface_cell_count"]),
        "stiffness_l1_free": f(row, "stiffness_l1_free"),
        "gradient_squared_l1_free": f(row, "gradient_squared_l1_free"),
        "reaction_l1_free": f(row, "reaction_l1_free"),
        "total_l1_free": f(row, "total_l1_free"),
        "interface_total_l1_free": f(row, "interface_total_l1_free"),
        "max_cell_residual_free": f(row, "max_cell_residual_free"),
        "max_cell_id": int(row["max_cell_id"]),
    } for row in region_rows]
    regions.sort(key=lambda row: row["total_l1_free"], reverse=True)

    return {
        "prefix": str(prefix),
        "free_cell_vertex_contribution_count": len(free),
        "free_node_count": len(free_nodes),
        "max_free_node": top_nodes[0],
        "component_l1": term_l1,
        "component_l1_fraction": {
            key: term_l1[key] / component_sum
            for key in ("stiffness", "gradient_squared", "reaction")
        },
        "interface_total_l1": interface_l1,
        "interface_fraction_of_cell_total_l1": interface_l1 / term_l1["total"],
        "domain_total_l1": domain_l1,
        "domain_fraction_of_cell_total_l1": {
            key: value / term_l1["total"] for key, value in domain_l1.items()
        },
        "domain_max_abs_aggregated_node_residual": domain_max,
        "regions": regions,
        "top_cells": cells[:10],
        "top_nodes": top_nodes,
    }


def number(value: float) -> str:
    return f"{value:.6g}"


def markdown(results: dict[str, dict[str, object]]) -> str:
    out = ["# SingleDevice fixed-state Eq. 231 residual decomposition", ""]
    out += [
        "| endpoint | max free residual | stiffness L1 | gradient-square L1 | reaction L1 | interface share | insulator share |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for tag, result in results.items():
        node = result["max_free_node"]
        l1 = result["component_l1"]
        out.append(
            f"| {tag} | {number(abs(node['raw_total']))} (node {node['node_id']}) | "
            f"{number(l1['stiffness'])} | {number(l1['gradient_squared'])} | "
            f"{number(l1['reaction'])} | "
            f"{100.0 * result['interface_fraction_of_cell_total_l1']:.3f}% | "
            f"{100.0 * result['domain_fraction_of_cell_total_l1']['insulator']:.3f}% |"
        )
    for tag, result in results.items():
        out += ["", f"## {tag}: regions", "", "| region | material | total L1 | gradient-square L1 | interface L1 | max cell |", "|---|---|---:|---:|---:|---:|"]
        for region in result["regions"]:
            out.append(
                f"| {region['region']} | {region['material']} | "
                f"{number(region['total_l1_free'])} | "
                f"{number(region['gradient_squared_l1_free'])} | "
                f"{number(region['interface_total_l1_free'])} | "
                f"{region['max_cell_id']} ({number(region['max_cell_residual_free'])}) |"
            )
        out += ["", f"## {tag}: largest cells", "", "| cell | region | interface | centroid | max local residual | gradient-square L1 |", "|---:|---|---|---|---:|---:|"]
        for cell in result["top_cells"]:
            out.append(
                f"| {cell['cell_id']} | {cell['region']} | "
                f"{'yes' if cell['interface_cell'] else 'no'} | "
                f"({number(cell['centroid_x_internal'])}, {number(cell['centroid_y_internal'])}) | "
                f"{number(cell['max_abs_local_total'])} | "
                f"{number(cell['gradient_squared_l1'])} |"
            )
    out.append("")
    return "\n".join(out)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    args = parser.parse_args()
    results = {tag: summarize(args.root / tag) for tag in ("lin", "sat")}
    (args.root / "summary.json").write_text(
        json.dumps(results, indent=2), encoding="utf-8")
    (args.root / "summary.md").write_text(markdown(results), encoding="utf-8")
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
