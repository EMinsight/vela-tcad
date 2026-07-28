#!/usr/bin/env python3
"""Replay exact Sentaurus states through the Vela general-Tri3 SG operator."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
import subprocess
from pathlib import Path


Q = 1.602176634e-19
CANDIDATES = (
    "net_doping",
    "total_impurity",
    "cell_reconstructed_total_impurity",
)
BIASES = (1, 2, 5, 10, 15, 20)


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def scalar_map(fields: Path, bias: int, name: str) -> dict[int, float]:
    return {
        int(row["node_id"]): float(row["component0"])
        for row in rows(fields / f"{bias}v" / "fields" / f"{name}_region0.csv")
    }


def write_state(fields: Path, bias: int, path: Path) -> None:
    names = {
        "psi_V": "ElectrostaticPotential",
        "phin_V": "eQuasiFermiPotential",
        "phip_V": "hQuasiFermiPotential",
        "n_m3": "eDensity",
        "p_m3": "hDensity",
    }
    values = {
        column: scalar_map(fields, bias, source)
        for column, source in names.items()
    }
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=["node_id", "psi_V", "phin_V", "phip_V", "n_m3", "p_m3"],
        )
        writer.writeheader()
        for node in sorted(values["psi_V"]):
            writer.writerow(
                {
                    "node_id": node,
                    "psi_V": values["psi_V"][node],
                    "phin_V": values["phin_V"][node],
                    "phip_V": values["phip_V"][node],
                    "n_m3": values["n_m3"][node] * 1.0e6,
                    "p_m3": values["p_m3"][node] * 1.0e6,
                }
            )


def node_volumes_m2(mesh: dict) -> list[float]:
    points = {int(node["id"]): (float(node["x"]), float(node["y"]))
              for node in mesh["nodes"]}
    volumes = [0.0] * len(points)
    for cell in mesh["triangles"]:
        nodes = cell["node_ids"]
        a, b, c = (points[int(node)] for node in nodes)
        area_um2 = abs(
            (b[0] - a[0]) * (c[1] - a[1])
            - (b[1] - a[1]) * (c[0] - a[0])
        ) / 2.0
        for node in nodes:
            volumes[int(node)] += area_um2 / 3.0
    return [value * 1.0e-12 for value in volumes]


def contact_nodes(mesh: dict) -> set[int]:
    return {
        int(node)
        for contact in mesh.get("contacts", [])
        for node in contact["node_ids"]
    }


def normalized_divergence(
    edge_rows: list[dict[str, str]],
    geometry: dict[tuple[int, int], dict[str, str]],
    source_cm3_s: dict[int, float],
    volumes: list[float],
    contacts: set[int],
    carrier: str,
    selected_nodes: set[int] | None = None,
) -> tuple[float, float]:
    divergence = [0.0] * len(volumes)
    incident = [0.0] * len(volumes)
    column = f"{carrier}_raw_signed_flux_per_m2_s"
    for edge in edge_rows:
        n0, n1 = int(edge["node0"]), int(edge["node1"])
        geom = geometry[(n0, n1)]
        face_m = float(geom["couple_um"]) * 1.0e-6
        integrated = float(edge[column]) * face_m
        divergence[n0] += integrated
        divergence[n1] -= integrated
        incident[n0] += abs(integrated)
        incident[n1] += abs(integrated)
    residuals = []
    scales = []
    for node in range(len(volumes)):
        if node in contacts or (
            selected_nodes is not None and node not in selected_nodes
        ):
            continue
        source = source_cm3_s[node] * 1.0e6 * volumes[node]
        residuals.append(abs(divergence[node] - source))
        scales.append(incident[node] + abs(source))
    total_scale = sum(scales)
    normalized_l1 = sum(residuals) / total_scale if total_scale else 0.0
    normalized_max = max(
        (residual / scale for residual, scale in zip(residuals, scales) if scale),
        default=0.0,
    )
    return normalized_l1, normalized_max


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--doping", type=Path, required=True)
    parser.add_argument("--baseline-config", type=Path, required=True)
    parser.add_argument("--sentaurus-fields", type=Path, required=True)
    parser.add_argument("--geometry-edges", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    mesh = json.loads(args.mesh.read_text(encoding="utf-8-sig"))
    baseline = json.loads(args.baseline_config.read_text(encoding="utf-8-sig"))
    geometry = {
        (int(row["node0"]), int(row["node1"])): row
        for row in rows(args.geometry_edges)
    }
    points = {
        int(node["id"]): (float(node["x"]), float(node["y"]))
        for node in mesh["nodes"]
    }
    volumes = node_volumes_m2(mesh)
    contacts = contact_nodes(mesh)
    junction_nodes = set(range(7, 16))
    summaries = []
    edge_comparisons = []

    for bias in BIASES:
        state_path = args.out_dir / f"sentaurus_state_{bias}V.csv"
        write_state(args.sentaurus_fields, bias, state_path)
        sent_e_current = scalar_map(
            args.sentaurus_fields, bias, "eCurrentDensity"
        )
        sent_h_current = scalar_map(
            args.sentaurus_fields, bias, "hCurrentDensity"
        )
        sent_srh = scalar_map(args.sentaurus_fields, bias, "srhRecombination")

        for basis in CANDIDATES:
            case = args.out_dir / basis / f"{bias}V"
            case.mkdir(parents=True, exist_ok=True)
            config = json.loads(json.dumps(baseline))
            config["solver"]["mobility"]["doping_concentration_basis"] = basis
            config_path = case / "audit.json"
            config_path.write_text(json.dumps(config, indent=2), encoding="utf-8")
            edge_path = case / "vela_edge_audit.csv"
            command = [
                str(args.audit.resolve()),
                "--mesh", str(args.mesh.resolve()),
                "--doping", str(args.doping.resolve()),
                "--state", str(state_path.resolve()),
                "--config", str(config_path.resolve()),
                "--node-out", str((case / "vela_node_state.csv").resolve()),
                "--edge-out", str(edge_path.resolve()),
                "--triangle-out", str((case / "vela_triangle_audit.csv").resolve()),
                "--scope", "general_tri3",
            ]
            completed = subprocess.run(
                command, text=True, capture_output=True, check=False
            )
            (case / "command.json").write_text(
                json.dumps(command, indent=2), encoding="utf-8"
            )
            (case / "stderr.log").write_text(
                completed.stderr, encoding="utf-8"
            )
            if completed.returncode:
                raise RuntimeError(
                    f"{basis} {bias} V fixed-state audit failed: "
                    f"{completed.stderr[-1000:]}"
                )
            edge_rows = rows(edge_path)
            electron_ratios = []
            hole_ratios = []
            junction_electron_ratios = []
            junction_hole_ratios = []
            zero_couple_raw_A_m2 = []
            for edge in edge_rows:
                n0, n1 = int(edge["node0"]), int(edge["node1"])
                geom = geometry[(n0, n1)]
                x0, y0 = points[n0]
                x1, y1 = points[n1]
                length = math.hypot(x1 - x0, y1 - y0)
                ux = (x1 - x0) / length
                e_vela = Q * abs(
                    float(edge["electron_raw_signed_flux_per_m2_s"])
                )
                h_vela = Q * abs(
                    float(edge["hole_raw_signed_flux_per_m2_s"])
                )
                e_sent = abs(0.5 * (
                    sent_e_current[n0] + sent_e_current[n1]
                ) * 1.0e4 * ux)
                h_sent = abs(0.5 * (
                    sent_h_current[n0] + sent_h_current[n1]
                ) * 1.0e4 * ux)
                positive_support = float(geom["couple_um"]) > 0.0
                horizontal = abs(ux) > 0.9
                junction = n0 in junction_nodes or n1 in junction_nodes
                if positive_support and horizontal and e_sent > 0.0:
                    electron_ratios.append(e_vela / e_sent)
                    if junction:
                        junction_electron_ratios.append(e_vela / e_sent)
                if positive_support and horizontal and h_sent > 0.0:
                    hole_ratios.append(h_vela / h_sent)
                    if junction:
                        junction_hole_ratios.append(h_vela / h_sent)
                if not positive_support:
                    zero_couple_raw_A_m2.append(e_vela + h_vela)
                edge_comparisons.append(
                    {
                        "basis": basis,
                        "bias_V": bias,
                        "edge_id": edge["edge_id"],
                        "node0": n0,
                        "node1": n1,
                        "couple_um": geom["couple_um"],
                        "junction_edge": int(junction),
                        "vela_e_abs_A_m2": e_vela,
                        "sentaurus_e_projected_abs_A_m2": e_sent,
                        "vela_h_abs_A_m2": h_vela,
                        "sentaurus_h_projected_abs_A_m2": h_sent,
                    }
                )
            e_l1, e_max = normalized_divergence(
                edge_rows, geometry, sent_srh, volumes, contacts, "electron"
            )
            h_l1, h_max = normalized_divergence(
                edge_rows, geometry, sent_srh, volumes, contacts, "hole"
            )
            je_l1, je_max = normalized_divergence(
                edge_rows, geometry, sent_srh, volumes, contacts, "electron",
                junction_nodes
            )
            jh_l1, jh_max = normalized_divergence(
                edge_rows, geometry, sent_srh, volumes, contacts, "hole",
                junction_nodes
            )
            summaries.append(
                {
                    "basis": basis,
                    "bias_V": bias,
                    "supported_horizontal_e_current_ratio_median": statistics.median(
                        electron_ratios
                    ),
                    "supported_horizontal_h_current_ratio_median": statistics.median(
                        hole_ratios
                    ),
                    "junction_supported_horizontal_e_ratio_median": statistics.median(
                        junction_electron_ratios
                    ),
                    "junction_supported_horizontal_h_ratio_median": statistics.median(
                        junction_hole_ratios
                    ),
                    "electron_continuity_normalized_L1": e_l1,
                    "hole_continuity_normalized_L1": h_l1,
                    "junction_electron_continuity_normalized_L1": je_l1,
                    "junction_hole_continuity_normalized_L1": jh_l1,
                    "electron_continuity_max_node_ratio": e_max,
                    "hole_continuity_max_node_ratio": h_max,
                    "junction_electron_continuity_max_node_ratio": je_max,
                    "junction_hole_continuity_max_node_ratio": jh_max,
                    "zero_couple_edge_count": len(zero_couple_raw_A_m2),
                    "zero_couple_raw_total_current_density_median_A_m2": (
                        statistics.median(zero_couple_raw_A_m2)
                    ),
                }
            )

    def write_csv(path: Path, values: list[dict]) -> None:
        with path.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=values[0].keys())
            writer.writeheader()
            writer.writerows(values)

    write_csv(args.out_dir / "fixed_state_summary.csv", summaries)
    write_csv(args.out_dir / "fixed_state_edge_comparison.csv", edge_comparisons)
    (args.out_dir / "fixed_state_summary.json").write_text(
        json.dumps({"cases": summaries}, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
