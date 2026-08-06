#!/usr/bin/env python3
"""Audit majority-carrier contact/interior QF drops in saved Vela states."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from pathlib import Path


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def bias_from_path(path: Path) -> float:
    match = re.search(r"bias_([mp]?\d+(?:p\d+)?)", path.stem)
    if not match:
        return math.nan
    token = match.group(1)
    sign = -1.0 if token.startswith("m") else 1.0
    token = token.lstrip("mp").replace("p", ".")
    return sign * float(token)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--doping", type=Path, required=True)
    parser.add_argument("--state-dir", type=Path, action="append", required=True)
    parser.add_argument("--out-csv", type=Path, required=True)
    args = parser.parse_args()

    mesh = json.loads(args.mesh.read_text(encoding="utf-8-sig"))
    silicon_regions = {
        int(region["id"])
        for region in mesh["regions"]
        if region.get("material") == "Si"
    }
    transport_edges: set[tuple[int, int]] = set()
    for triangle in mesh["triangles"]:
        if int(triangle["region_id"]) not in silicon_regions:
            continue
        nodes = [int(value) for value in triangle["node_ids"]]
        for a, b in ((nodes[0], nodes[1]), (nodes[1], nodes[2]), (nodes[2], nodes[0])):
            transport_edges.add(tuple(sorted((a, b))))

    net_doping: dict[int, float] = {}
    for row in read_rows(args.doping):
        node = int(row["node_id"])
        net_doping[node] = float(row["donors_cm3"]) - float(row["acceptors_cm3"])

    contacts: list[tuple[str, set[int], bool]] = []
    for contact in mesh["contacts"]:
        name = str(contact["name"])
        if name.lower() == "gate":
            continue
        nodes = {int(value) for value in contact["node_ids"]}
        mean = sum(net_doping.get(node, 0.0) for node in nodes) / max(len(nodes), 1)
        contacts.append((name, nodes, mean >= 0.0))

    paths: set[Path] = set()
    for directory in args.state_dir:
        paths.update(directory.rglob("accepted_state_bias_*.csv"))

    output: list[dict[str, object]] = []
    for path in sorted(paths, key=lambda value: (bias_from_path(value), str(value))):
        state = {int(row["node_id"]): row for row in read_rows(path)}
        best: dict[str, object] = {
            "bias_V": bias_from_path(path),
            "max_contact_majority_qf_drop_V": 0.0,
            "contact": "",
            "contact_node": -1,
            "interior_node": -1,
            "carrier": "",
            "state_file": str(path.resolve()),
        }
        for name, contact_nodes, electron_majority in contacts:
            qf_key = "phin" if electron_majority else "phip"
            density_key = "electrons_m3" if electron_majority else "holes_m3"
            for a, b in transport_edges:
                if (a in contact_nodes) == (b in contact_nodes):
                    continue
                row_a = state.get(a)
                row_b = state.get(b)
                if row_a is None or row_b is None:
                    continue
                density_a = float(row_a[density_key])
                density_b = float(row_b[density_key])
                if not (density_a > 0.0 and density_b > 0.0 and
                        math.isfinite(density_a) and math.isfinite(density_b)):
                    continue
                drop = abs(float(row_a[qf_key]) - float(row_b[qf_key]))
                if drop > float(best["max_contact_majority_qf_drop_V"]):
                    best.update({
                        "max_contact_majority_qf_drop_V": drop,
                        "contact": name,
                        "contact_node": a if a in contact_nodes else b,
                        "interior_node": b if a in contact_nodes else a,
                        "carrier": "electron" if electron_majority else "hole",
                    })
        output.append(best)

    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    fields = list(output[0]) if output else ["bias_V", "max_contact_majority_qf_drop_V"]
    with args.out_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(output)
    for row in output:
        print(f"{row['bias_V']:.6f},{row['max_contact_majority_qf_drop_V']:.17g},"
              f"{row['contact']},{row['contact_node']},{row['interior_node']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
