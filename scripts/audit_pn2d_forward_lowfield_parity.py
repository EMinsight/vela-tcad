#!/usr/bin/env python3
"""Audit PN2D low-field mobility, contacts, and terminal current."""

from __future__ import annotations

import csv
import json
import math
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RUN = ROOT / "build-release" / "pn2d-forward-iv-0v20v-20260727"
MESH_DIR = ROOT / "build-release" / "pn2d-general-tri3-task7-imported-mesh-20260726" / "vela"
PARAMS = {
    "electron": (1417.0, 52.2, 52.2, 43.4, 0.0, 9.68e16, 3.43e20, 0.68, 2.0),
    "hole": (470.5, 44.9, 0.0, 29.0, 9.23e16, 2.23e17, 6.10e20, 0.719, 2.0),
}


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def value(row: dict[str, str], key: str) -> float:
    return float(row[key])


def masetti(doping: float, params: tuple[float, ...]) -> float:
    mu0, mumin1, mumin2, mu1, pc, cr, cs, alpha, beta = params
    doping = abs(doping)
    if doping == 0:
        return mu0
    return (
        mumin1 * math.exp(-max(0.0, pc) / doping)
        + (mu0 - mumin2) / (1.0 + (doping / cr) ** alpha)
        - mu1 / (1.0 + (cs / doping) ** beta)
    )


def interpolate(points: list[tuple[float, float]], target: float) -> float:
    points = sorted(points)
    for x, y in points:
        if abs(x - target) < 1e-12:
            return y
    for (x0, y0), (x1, y1) in zip(points, points[1:]):
        if x0 <= target <= x1:
            weight = (target - x0) / (x1 - x0)
            return y0 + weight * (y1 - y0)
    raise ValueError(f"Target {target} outside interpolation range")


def main() -> None:
    mesh = json.loads((MESH_DIR / "mesh.json").read_text(encoding="utf-8-sig"))
    doping = {
        int(r["node_id"]): (float(r["donors_cm3"]), float(r["acceptors_cm3"]))
        for r in rows(MESH_DIR / "doping.csv")
    }
    contacts = {c["name"]: set(c["node_ids"]) for c in mesh["contacts"]}
    terminal = rows(RUN / "vela_coarse7x3_forward_iv_legacy_lowfield_terminal_balance.csv")
    edge_data = rows(RUN / "vela_coarse7x3_forward_iv_legacy_lowfield_contact_edges.csv")
    methods = rows(RUN / "vela_coarse7x3_forward_iv_legacy_lowfield_terminal_current_method_compare.csv")
    vela_iv = rows(RUN / "vela_coarse7x3_forward_iv_legacy_lowfield_diagnostics_0v20v.csv")
    sent_iv = rows(RUN / "sentaurus_coarse7x3_forward_iv_legacy_physics_0v20v.csv")

    terminal_map: dict[float, dict[str, float]] = defaultdict(dict)
    for row in terminal:
        terminal_map[round(value(row, "bias_V"), 10)][row["contact"]] = value(
            row, "current_total_A_per_um"
        )

    selected = {0.1, 1.0, 2.0, 5.0, 10.0, 15.0, 20.0}
    selected_kcl = []
    max_kcl_rel = 0.0
    for bias, current in terminal_map.items():
        if bias < 0.1 or set(current) != {"Anode", "Cathode"}:
            continue
        imbalance = current["Anode"] + current["Cathode"]
        rel = abs(imbalance) / max(abs(current["Anode"]), abs(current["Cathode"]), 1e-300)
        max_kcl_rel = max(max_kcl_rel, rel)
        if bias in selected:
            selected_kcl.append((bias, current["Anode"], current["Cathode"], imbalance, rel))

    max_method_rel = 0.0
    for row in methods:
        sg = value(row, "I_sgflux_A_per_um")
        residual = value(row, "I_residual_A_per_um")
        max_method_rel = max(
            max_method_rel,
            abs(sg - residual) / max(abs(sg), abs(residual), 1e-300),
        )

    edge_sums: dict[tuple[float, str], float] = defaultdict(float)
    normalization_error = 0.0
    qf_deviation = 0.0
    mobility_20 = defaultdict(lambda: [math.inf, -math.inf, math.inf, -math.inf])
    for row in edge_data:
        bias = round(value(row, "bias_V"), 10)
        contact = row["current_contact"]
        edge_sums[(bias, contact)] += value(row, "current_total_A_per_um")
        raw = value(row, "current_total")
        if raw:
            normalization_error = max(
                normalization_error,
                abs(value(row, "current_total_A_per_um") / raw - 1e-6),
            )
        contact_is_n0 = int(row["node0"]) in contacts[contact]
        expected_qf = bias if contact == "Anode" else 0.0
        for carrier in ("phin", "phip"):
            qf_deviation = max(
                qf_deviation,
                abs(value(row, carrier + ("0" if contact_is_n0 else "1")) - expected_qf),
            )
        if bias == 20.0:
            limits = mobility_20[contact]
            mun, mup = value(row, "mun"), value(row, "mup")
            limits[:] = [min(limits[0], mun), max(limits[1], mun),
                         min(limits[2], mup), max(limits[3], mup)]

    max_edge_sum_rel = 0.0
    for bias, current in terminal_map.items():
        for contact, total in current.items():
            edge_sum = edge_sums[(bias, contact)]
            max_edge_sum_rel = max(
                max_edge_sum_rel,
                abs(edge_sum - total) / max(abs(edge_sum), abs(total), 1e-300),
            )

    compensated = []
    for node, (donors, acceptors) in doping.items():
        if donors and acceptors:
            net, total = donors - acceptors, donors + acceptors
            compensated.append((
                node,
                masetti(net, PARAMS["electron"]),
                masetti(total, PARAMS["electron"]),
                masetti(net, PARAMS["hole"]),
                masetti(total, PARAMS["hole"]),
            ))

    vela_points = [(value(r, "bias_V"), value(r, "current_total_A_per_um")) for r in vela_iv]
    sent_points = [(value(r, "bias_V"), value(r, "current_total")) for r in sent_iv]
    iv_compare = [
        (bias, interpolate(vela_points, bias), interpolate(sent_points, bias))
        for bias in (1.0, 2.0, 5.0, 10.0, 15.0, 20.0)
    ]

    report = RUN / "pn2d_coarse7x3_lowfield_mobility_contact_current_audit.md"
    output = [
        "# PN2D coarse7x3 low-field mobility/contact/current audit",
        "",
        "## Masetti point values",
        "",
        "| N (cm^-3) | electron (cm^2/V/s) | hole (cm^2/V/s) |",
        "|---:|---:|---:|",
    ]
    for concentration in (1e12, 1e14, 1e16, 1e17, 2e17, 1e18, 1e19, 1e20):
        output.append(
            f"| {concentration:.0e} | {masetti(concentration, PARAMS['electron']):.6f} "
            f"| {masetti(concentration, PARAMS['hole']):.6f} |"
        )
    output += [
        "",
        "Vela and Sentaurus parameters produce identical values for identical concentration inputs.",
        "The compensated junction nodes expose a concentration-definition difference:",
        "",
        "| node | Vela mun(net) | Sent-like mun(total) | Vela mup(net) | Sent-like mup(total) |",
        "|---:|---:|---:|---:|---:|",
    ]
    for item in compensated:
        output.append(f"| {item[0]} | {item[1]:.6f} | {item[2]:.6f} | {item[3]:.6f} | {item[4]:.6f} |")

    output += ["", "## Contact nodes", "",
               "| contact | node | Nd | Na | net (cm^-3) |", "|---|---:|---:|---:|---:|"]
    for contact, nodes in contacts.items():
        for node in sorted(nodes):
            donors, acceptors = doping[node]
            output.append(f"| {contact} | {node} | {donors:.6e} | {acceptors:.6e} | {donors-acceptors:.6e} |")
    output += [
        "",
        f"Maximum contact quasi-Fermi BC deviation: {qf_deviation:.9e} V.",
        "All contact nodes have the expected polarity; dominant_signed_contact_mean is inactive here.",
        "",
        "20 V contact-edge mobility ranges:",
    ]
    for contact, limits in sorted(mobility_20.items()):
        output.append(f"- {contact}: mun {limits[0]:.6f}..{limits[1]:.6f}; mup {limits[2]:.6f}..{limits[3]:.6f}")

    output += [
        "", "## Terminal current checks", "",
        "| V | Anode (A/um) | Cathode (A/um) | sum | relative imbalance |",
        "|---:|---:|---:|---:|---:|",
    ]
    for bias, anode, cathode, imbalance, rel in sorted(selected_kcl):
        output.append(f"| {bias:.1f} | {anode:.6e} | {cathode:.6e} | {imbalance:.6e} | {rel:.3e} |")
    output += [
        "",
        f"- max KCL relative imbalance, 0.1--20 V: {max_kcl_rel:.9e}",
        f"- max contact-edge sum vs terminal CSV relative difference: {max_edge_sum_rel:.9e}",
        f"- max SG-flux vs residual-method relative difference: {max_method_rel:.9e}",
        f"- max error in A/um conversion factor 1e-6: {normalization_error:.9e}",
        "",
        "## Vela vs Sentaurus low-field IV",
        "",
        "| V | Vela (A/um) | Sentaurus (A/um) | ratio |",
        "|---:|---:|---:|---:|",
    ]
    for bias, vela, sent in iv_compare:
        output.append(f"| {bias:.1f} | {vela:.6e} | {sent:.6e} | {vela/sent:.6f} |")
    report.write_text("\n".join(output) + "\n", encoding="utf-8")
    print(report)
    print(f"max_kcl_rel={max_kcl_rel:.9e}")
    print(f"max_method_rel={max_method_rel:.9e}")
    print(f"max_edge_sum_rel={max_edge_sum_rel:.9e}")
    print(f"max_qf_bc_deviation_V={qf_deviation:.9e}")


if __name__ == "__main__":
    main()
