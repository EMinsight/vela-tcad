#!/usr/bin/env python3
"""Compare matched coarse7x3 Vela/Sentaurus forward states and geometry."""

from __future__ import annotations

import csv
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "build-release" / "pn2d-forward-field-audit-20260727"
MESH_PATH = (
    ROOT
    / "build-release"
    / "pn2d-general-tri3-task7-imported-mesh-20260726"
    / "vela"
    / "mesh.json"
)
BIAS_INDEX = {0.0: 0, 1.0: 10, 2.0: 20, 5.0: 50, 10.0: 100, 15.0: 150, 20.0: 200}
Q_OVER_K_300 = 1.0 / (8.617333262145e-5 * 300.0)
BASE_NI_CM3 = 14638914958.767616


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def read_scalar(path: Path) -> dict[int, float]:
    return {int(float(r["node_id"])): float(r["component0"]) for r in read_csv(path)}


def parse_vtk_points(path: Path) -> dict[str, list[float]]:
    lines = path.read_text(encoding="ascii").splitlines()
    point_start = next(i for i, line in enumerate(lines) if line.startswith("POINT_DATA "))
    count = int(lines[point_start].split()[1])
    result: dict[str, list[float]] = {}
    i = point_start + 1
    while i < len(lines):
        parts = lines[i].split()
        if len(parts) >= 2 and parts[0] == "SCALARS":
            name = parts[1]
            i += 1
            if i < len(lines) and lines[i].startswith("LOOKUP_TABLE"):
                i += 1
            values: list[float] = []
            while i < len(lines) and len(values) < count:
                values.extend(float(v) for v in lines[i].split())
                i += 1
            result[name] = values[:count]
            continue
        if len(parts) >= 2 and parts[0] == "VECTORS":
            i += count + 1
            continue
        i += 1
    return result


def vela_path(bias: float) -> Path:
    index = BIAS_INDEX[bias]
    token = f"{bias:g}V"
    return AUDIT / f"vela_forward_fields_{index:04d}_{token}.vtk"


def sent_field(bias: float, name: str) -> dict[int, float]:
    token = f"{bias:g}v"
    return read_scalar(AUDIT / "sentaurus_fields" / token / "fields" / f"{name}_region0.csv")


def rms(values: list[float]) -> float:
    return math.sqrt(sum(v * v for v in values) / len(values)) if values else 0.0


def field_metrics(
    reference: dict[int, float],
    candidate: dict[int, float],
    log: bool,
    relative_floor: float = 0.0,
) -> dict[str, float]:
    nodes = sorted(set(reference) & set(candidate))
    if log:
        scale = max(abs(reference[n]) for n in nodes)
        floor = scale * relative_floor
        selected = [
            n for n in nodes
            if abs(reference[n]) > floor and abs(candidate[n]) > floor
        ]
        errors = [
            math.log10(abs(candidate[n]) / abs(reference[n]))
            for n in selected
        ]
    else:
        errors = [candidate[n] - reference[n] for n in nodes]
    return {
        "count": len(errors),
        "mean": sum(errors) / len(errors),
        "rms": rms(errors),
        "median_abs": statistics.median(abs(v) for v in errors),
        "max_abs": max(abs(v) for v in errors),
    }


def triangle_area(points: list[tuple[float, float]], ids: list[int]) -> float:
    a, b, c = (points[i] for i in ids)
    return 0.5 * abs((b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1]))


def cotangent(a: tuple[float, float], b: tuple[float, float], o: tuple[float, float]) -> float:
    ux, uy = a[0] - o[0], a[1] - o[1]
    vx, vy = b[0] - o[0], b[1] - o[1]
    cross = ux * vy - uy * vx
    return (ux * vx + uy * vy) / abs(cross)


def geometry_audit(mesh: dict[str, Any]) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    nodes_by_id = {int(n["id"]): (float(n["x"]), float(n["y"])) for n in mesh["nodes"]}
    points = [nodes_by_id[i] for i in range(len(nodes_by_id))]
    node_volume = [0.0] * len(points)
    edge_data: dict[tuple[int, int], dict[str, float]] = {}
    min_angle, max_angle = 180.0, 0.0
    negative_cot = 0
    total_area = 0.0
    for tri in mesh["triangles"]:
        ids = [int(v) for v in tri["node_ids"]]
        area = triangle_area(points, ids)
        total_area += area
        for node in ids:
            node_volume[node] += area / 3.0
        for k in range(3):
            a, b, opp = ids[k], ids[(k + 1) % 3], ids[(k + 2) % 3]
            key = tuple(sorted((a, b)))
            length = math.dist(points[a], points[b])
            cot = cotangent(points[a], points[b], points[opp])
            if cot < 0.0:
                negative_cot += 1
                local_couple = area / (3.0 * length)
            else:
                local_couple = 0.5 * cot * length
            item = edge_data.setdefault(key, {"length": length, "couple": 0.0, "cells": 0.0})
            item["couple"] += max(local_couple, 0.0)
            item["cells"] += 1.0
        a, b, c = (points[i] for i in ids)
        for p0, p1, p2 in ((a, b, c), (b, c, a), (c, a, b)):
            u, v = (p1[0] - p0[0], p1[1] - p0[1]), (p2[0] - p0[0], p2[1] - p0[1])
            cosine = max(-1.0, min(1.0, (u[0] * v[0] + u[1] * v[1]) / (math.hypot(*u) * math.hypot(*v))))
            angle = math.degrees(math.acos(cosine))
            min_angle, max_angle = min(min_angle, angle), max(max_angle, angle)
    edge_rows = []
    for (a, b), item in sorted(edge_data.items()):
        edge_rows.append({
            "node0": a,
            "node1": b,
            "length_um": item["length"],
            "couple_um": item["couple"],
            "couple_over_length": item["couple"] / item["length"],
            "adjacent_cells": int(item["cells"]),
        })
    summary = {
        "nodes": len(points),
        "triangles": len(mesh["triangles"]),
        "edges": len(edge_rows),
        "domain_area_um2": total_area,
        "node_volume_sum_um2": sum(node_volume),
        "node_volume_min_um2": min(node_volume),
        "node_volume_max_um2": max(node_volume),
        "min_angle_deg": min_angle,
        "max_angle_deg": max_angle,
        "negative_cotangent_contributions": negative_cot,
        "zero_or_negative_couples": sum(1 for r in edge_rows if r["couple_um"] <= 0.0),
        "couple_over_length_min": min(r["couple_over_length"] for r in edge_rows),
        "couple_over_length_max": max(r["couple_over_length"] for r in edge_rows),
    }
    return summary, edge_rows


def group_for_x(x: float) -> str:
    if abs(x - 1.0) < 1e-12:
        return "junction"
    return "p_side" if x < 1.0 else "n_side"


def main() -> None:
    mesh = json.loads(MESH_PATH.read_text(encoding="utf-8-sig"))
    nodes = {int(n["id"]): (float(n["x"]), float(n["y"])) for n in mesh["nodes"]}
    geometry, edge_rows = geometry_audit(mesh)
    with (AUDIT / "geometry_edges.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(edge_rows[0]))
        writer.writeheader()
        writer.writerows(edge_rows)

    quantities = {
        "psi": ("ElectrostaticPotential", "Potential", 1.0, False, 0.0),
        "n": ("eDensity", "Electrons", 1.0, True, 0.0),
        "p": ("hDensity", "Holes", 1.0, True, 0.0),
        "phin": ("eQuasiFermiPotential", "ElectronQuasiFermi", 1.0, False, 0.0),
        "phip": ("hQuasiFermiPotential", "HoleQuasiFermi", 1.0, False, 0.0),
        "srh": ("srhRecombination", "SRHRecombination", 1.0, True, 1e-12),
    }
    comparisons: dict[str, Any] = {}
    detail_rows: list[dict[str, Any]] = []
    state_cache: dict[float, dict[str, dict[int, float]]] = {}
    for bias in BIAS_INDEX:
        vtk = parse_vtk_points(vela_path(bias))
        states: dict[str, dict[int, float]] = {}
        comparisons[f"{bias:g}V"] = {}
        for short, (sent_name, vela_name, factor, use_log, relative_floor) in quantities.items():
            sent = sent_field(bias, sent_name)
            vela = {i: factor * value for i, value in enumerate(vtk[vela_name])}
            states[f"sent_{short}"] = sent
            states[f"vela_{short}"] = vela
            comparisons[f"{bias:g}V"][short] = field_metrics(sent, vela, use_log, relative_floor)
        state_cache[bias] = states
        for node in sorted(nodes):
            row: dict[str, Any] = {
                "bias_V": bias,
                "node_id": node,
                "x_um": nodes[node][0],
                "y_um": nodes[node][1],
                "group": group_for_x(nodes[node][0]),
            }
            for short in quantities:
                sent = states[f"sent_{short}"][node]
                vela = states[f"vela_{short}"][node]
                row[f"sent_{short}"] = sent
                row[f"vela_{short}"] = vela
                row[f"vela_minus_sent_{short}"] = vela - sent
                if short in {"n", "p", "srh"}:
                    row[f"vela_over_sent_{short}"] = vela / sent if sent != 0 else None
            detail_rows.append(row)

    with (AUDIT / "node_field_comparison.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(detail_rows[0]))
        writer.writeheader()
        writer.writerows(detail_rows)

    zero = state_cache[0.0]
    ni_groups: dict[str, dict[str, list[float]]] = defaultdict(lambda: {"sent": [], "vela": []})
    for node, (x, _) in nodes.items():
        group = group_for_x(x)
        ni_groups[group]["sent"].append(math.sqrt(max(zero["sent_n"][node] * zero["sent_p"][node], 0.0)))
        ni_groups[group]["vela"].append(math.sqrt(max(zero["vela_n"][node] * zero["vela_p"][node], 0.0)))
    inferred_ni = {
        group: {
            "sentaurus_median_cm3": statistics.median(values["sent"]),
            "vela_median_cm3": statistics.median(values["vela"]),
            "vela_over_sentaurus": statistics.median(values["vela"]) / statistics.median(values["sent"]),
        }
        for group, values in ni_groups.items()
    }

    contact_nodes = {
        contact["name"]: [int(n) for n in contact["node_ids"]]
        for contact in mesh["contacts"]
    }
    contact_built_in: dict[str, Any] = {}
    for name, ids in contact_nodes.items():
        contact_built_in[name] = {
            "sentaurus_psi_minus_phin_V": statistics.mean(
                zero["sent_psi"][n] - zero["sent_phin"][n] for n in ids
            ),
            "vela_psi_minus_phin_V": statistics.mean(
                zero["vela_psi"][n] - zero["vela_phin"][n] for n in ids
            ),
            "sentaurus_psi_minus_phip_V": statistics.mean(
                zero["sent_psi"][n] - zero["sent_phip"][n] for n in ids
            ),
            "vela_psi_minus_phip_V": statistics.mean(
                zero["vela_psi"][n] - zero["vela_phip"][n] for n in ids
            ),
        }

    centerline = []
    high = state_cache[20.0]
    for node, (x, y) in sorted(nodes.items(), key=lambda item: (item[1][0], item[1][1])):
        if abs(y - 0.25) > 1e-12:
            continue
        centerline.append({
            "node": node,
            "x_um": x,
            **{
                f"{side}_{short}": high[f"{side}_{short}"][node]
                for short in ("psi", "n", "p", "phin", "phip")
                for side in ("sent", "vela")
            },
        })
    with (AUDIT / "centerline_20v.csv").open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(centerline[0]))
        writer.writeheader()
        writer.writerows(centerline)

    summary = {
        "schema": "vela.pn2d.forward_field_audit.v1",
        "biases_V": list(BIAS_INDEX),
        "geometry": geometry,
        "old_slotboom": {
            "base_ni_cm3": BASE_NI_CM3,
            "vela_parameters": {"offset_eV": 0.0, "Ebgn_eV": 0.009, "Nref_cm3": 1e17, "C": 0.5},
            "sentaurus_parameter_file": {"dEg0_eV": -0.01595, "Ebgn_eV": 0.009, "Nref_cm3": 1e17, "C": 0.5},
            "inferred_ni_eff": inferred_ni,
        },
        "contact_built_in": contact_built_in,
        "field_comparison": comparisons,
    }
    (AUDIT / "forward_field_audit_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )

    lines = [
        "# PN2D coarse7x3 forward field audit",
        "",
        "Matched physics: DopingDependence + SRH + OldSlotboom, avalanche off, no high-field mobility.",
        "",
        "## Geometry",
        "",
        f"- topology: {geometry['nodes']} nodes, {geometry['triangles']} triangles, {geometry['edges']} edges",
        f"- domain area / node-volume closure: {geometry['domain_area_um2']:.12g} / {geometry['node_volume_sum_um2']:.12g} um^2",
        f"- angle range: {geometry['min_angle_deg']:.6g} to {geometry['max_angle_deg']:.6g} deg",
        f"- negative cotangent fallbacks: {geometry['negative_cotangent_contributions']}",
        f"- zero/negative SG couples: {geometry['zero_or_negative_couples']}",
        f"- couple/length range: {geometry['couple_over_length_min']:.6g} to {geometry['couple_over_length_max']:.6g}",
        "",
        "## Equilibrium effective intrinsic density",
        "",
        "| region | Sentaurus median (cm^-3) | Vela median (cm^-3) | Vela/Sentaurus |",
        "|---|---:|---:|---:|",
    ]
    for group, values in inferred_ni.items():
        lines.append(
            f"| {group} | {values['sentaurus_median_cm3']:.6e} | "
            f"{values['vela_median_cm3']:.6e} | {values['vela_over_sentaurus']:.9f} |"
        )
    lines += [
        "",
        "## Contact built-in potential at 0 V",
        "",
        "| contact | Sent psi-phin | Vela psi-phin | Sent psi-phip | Vela psi-phip |",
        "|---|---:|---:|---:|---:|",
    ]
    for contact, values in contact_built_in.items():
        lines.append(
            f"| {contact} | {values['sentaurus_psi_minus_phin_V']:.9f} | "
            f"{values['vela_psi_minus_phin_V']:.9f} | "
            f"{values['sentaurus_psi_minus_phip_V']:.9f} | "
            f"{values['vela_psi_minus_phip_V']:.9f} |"
        )
    lines += [
        "",
        "## Spatial error summary",
        "",
        "Density entries are log10(Vela/Sentaurus). SRH uses log10(abs(Vela)/abs(Sentaurus)) after filtering numerical-zero nodes; potential/QF entries are Vela-Sentaurus in volts.",
        "",
        "| bias | psi RMS | n RMS log10 | p RMS log10 | phin RMS | phip RMS | SRH RMS log10 |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for bias in BIAS_INDEX:
        entry = comparisons[f"{bias:g}V"]
        lines.append(
            f"| {bias:g} | {entry['psi']['rms']:.6g} | {entry['n']['rms']:.6g} | "
            f"{entry['p']['rms']:.6g} | {entry['phin']['rms']:.6g} | "
            f"{entry['phip']['rms']:.6g} | {entry['srh']['rms']:.6g} |"
        )
    lines += [
        "",
        "## 20 V centerline",
        "",
        "| x um | Sent n | Vela n | Sent p | Vela p | Sent phin | Vela phin | Sent phip | Vela phip |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in centerline:
        lines.append(
            f"| {row['x_um']:.2f} | {row['sent_n']:.4e} | {row['vela_n']:.4e} | "
            f"{row['sent_p']:.4e} | {row['vela_p']:.4e} | "
            f"{row['sent_phin']:.6f} | {row['vela_phin']:.6f} | "
            f"{row['sent_phip']:.6f} | {row['vela_phip']:.6f} |"
        )
    (AUDIT / "pn2d_forward_field_audit.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps({
        "report": str(AUDIT / "pn2d_forward_field_audit.md"),
        "summary": str(AUDIT / "forward_field_audit_summary.json"),
        "detail": str(AUDIT / "node_field_comparison.csv"),
    }, indent=2))


if __name__ == "__main__":
    main()
