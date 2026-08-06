#!/usr/bin/env python3
"""Compare millivolt Vela/Sentaurus QF states and SG edge fluxes."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import subprocess
from collections import deque
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[1]
Q = 1.602176634e-19


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, fields: list[str], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def bias_tag(value: float) -> str:
    return f"{value:.6f}".replace("-", "m").replace(".", "p")


def load_scalar(path: Path) -> dict[int, float]:
    result: dict[int, float] = {}
    for row in read_csv(path):
        result[int(row["node_id"])] = float(row["component0"])
    return result


def load_vector(path: Path) -> dict[int, tuple[float, float]]:
    result: dict[int, tuple[float, float]] = {}
    for row in read_csv(path):
        result[int(row["node_id"])] = (
            float(row["component0"]),
            float(row["component1"]),
        )
    return result


def write_probe_fields(path: Path, values: dict[str, dict[int, float]], count: int) -> None:
    path.mkdir(parents=True, exist_ok=True)
    for field, data in values.items():
        rows = [{"node_id": node, "component0": data.get(node, 0.0)} for node in range(count)]
        write_csv(path / f"{field}_region0.csv", ["node_id", "component0"], rows)


def merge_sentaurus_state(imported: Path, count: int) -> dict[str, dict[int, float]]:
    fields = imported / "fields"
    result: dict[str, dict[int, float]] = {}
    for name in ("ElectrostaticPotential", "eQuasiFermiPotential", "hQuasiFermiPotential"):
        merged: dict[int, float] = {}
        # Potential is continuous. For QF, process region 3 (Si substrate)
        # last so shared Si/insulator nodes receive the transport-side value.
        paths = sorted(fields.glob(f"{name}_region*.csv"))
        paths.sort(key=lambda path: path.name.endswith("region3.csv"))
        for path in paths:
            merged.update(load_scalar(path))
        if len(merged) != count:
            missing = count - len(merged)
            raise RuntimeError(f"{name}: merged node coverage incomplete by {missing}")
        result[name] = merged
    return result


def vela_state_fields(path: Path, count: int) -> dict[str, dict[int, float]]:
    rows = read_csv(path)
    if len(rows) != count:
        raise RuntimeError(f"Vela state {path} has {len(rows)} rows, expected {count}")
    columns = {
        "ElectrostaticPotential": "psi",
        "eQuasiFermiPotential": "phin",
        "hQuasiFermiPotential": "phip",
    }
    return {
        field: {int(row["node_id"]): float(row[column]) for row in rows}
        for field, column in columns.items()
    }


def run_probe(
    runner: Path,
    base: dict[str, Any],
    bias: float,
    state_fields: Path,
    output_csv: Path,
    config_path: Path,
) -> None:
    cfg = json.loads(json.dumps(base))
    cfg["simulation_type"] = "sg_edge_flux_probe"
    cfg["state_fields_dir"] = str(state_fields.resolve())
    cfg["output_csv"] = str(output_csv.resolve())
    cfg.pop("sweep", None)
    for contact in cfg["contacts"]:
        if contact["name"].lower() == "drain":
            contact["bias"] = bias
    config_path.parent.mkdir(parents=True, exist_ok=True)
    output_csv.parent.mkdir(parents=True, exist_ok=True)
    config_path.write_text(json.dumps(cfg, indent=2) + "\n", encoding="utf-8")
    completed = subprocess.run(
        [str(runner.resolve()), "--config", str(config_path.resolve())],
        cwd=config_path.parent,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(completed.stderr.strip() or completed.stdout.strip())


def hop_distance(mesh: dict[str, Any], edges: list[dict[str, str]]) -> dict[int, int]:
    adjacency: dict[int, set[int]] = {int(node["id"]): set() for node in mesh["nodes"]}
    for edge in edges:
        n0, n1 = int(edge["node0"]), int(edge["node1"])
        adjacency[n0].add(n1)
        adjacency[n1].add(n0)
    drain = next(contact for contact in mesh["contacts"] if contact["name"].lower() == "drain")
    distance = {int(node): 0 for node in drain["node_ids"]}
    queue = deque(distance)
    while queue:
        node = queue.popleft()
        for other in adjacency[node]:
            if other not in distance:
                distance[other] = distance[node] + 1
                queue.append(other)
    return distance


def semiconductor_edges(mesh: dict[str, Any]) -> set[tuple[int, int]]:
    semiconductor_regions = {
        int(region["id"])
        for region in mesh["regions"]
        if region["material"].strip().lower() in {"si", "silicon"}
    }
    result: set[tuple[int, int]] = set()
    for triangle in mesh["triangles"]:
        if int(triangle["region_id"]) not in semiconductor_regions:
            continue
        a, b, c = (int(node) for node in triangle["node_ids"])
        result.update(tuple(sorted(edge)) for edge in ((a, b), (b, c), (c, a)))
    return result


def signed_dex(numerator: float, denominator: float, floor: float) -> float:
    return math.log10((abs(numerator) + floor) / (abs(denominator) + floor))


def analyze_bias(
    bias: float,
    vela_edge_path: Path,
    sent_edge_path: Path,
    sent_imported: Path,
    mesh: dict[str, Any],
    output: Path,
) -> dict[str, Any]:
    vela_rows = read_csv(vela_edge_path)
    sent_rows = {int(row["edge_id"]): row for row in read_csv(sent_edge_path)}
    distances = hop_distance(mesh, vela_rows)
    silicon_edges = semiconductor_edges(mesh)
    drain_nodes = {
        int(node)
        for contact in mesh["contacts"]
        if contact["name"].lower() == "drain"
        for node in contact["node_ids"]
    }
    sent_current = load_vector(sent_imported / "fields/eCurrentDensity_region3.csv")
    max_native = 0.0
    max_state_line_flux = 0.0
    native_by_edge: dict[int, float] = {}
    for row in vela_rows:
        edge_id = int(row["edge_id"])
        n0, n1 = int(row["node0"]), int(row["node1"])
        if tuple(sorted((n0, n1))) not in silicon_edges or float(row["couple_m"]) <= 0.0:
            continue
        sent_row = sent_rows[edge_id]
        max_state_line_flux = max(
            max_state_line_flux,
            abs(float(row["electron_particle_line_flux_per_m_s"])),
            abs(float(sent_row["electron_particle_line_flux_per_m_s"])),
        )
        if n0 not in sent_current or n1 not in sent_current:
            continue
        dx = float(row["x1"]) - float(row["x0"])
        dy = float(row["y1"]) - float(row["y0"])
        length = math.hypot(dx, dy)
        if length == 0.0:
            continue
        jx = 0.5 * (sent_current[n0][0] + sent_current[n1][0])
        jy = 0.5 * (sent_current[n0][1] + sent_current[n1][1])
        # Sentaurus current is A/cm^2. Convert its edge-parallel magnitude to
        # particle flux m^-2 s^-1 for comparison with Vela SG electronFlux.
        particle_flux = abs((jx * dx + jy * dy) / length) * 1.0e4 / Q
        particle_line_flux = particle_flux * float(row["couple_m"])
        native_by_edge[edge_id] = particle_line_flux
        max_native = max(max_native, particle_line_flux)

    compared: list[dict[str, Any]] = []
    for vela in vela_rows:
        edge_id = int(vela["edge_id"])
        sent = sent_rows[edge_id]
        n0, n1 = int(vela["node0"]), int(vela["node1"])
        is_silicon = tuple(sorted((n0, n1))) in silicon_edges
        active_sg = is_silicon and float(vela["couple_m"]) > 0.0
        drain_endpoint_count = int(n0 in drain_nodes) + int(n1 in drain_nodes)
        vflux = float(vela["electron_particle_line_flux_per_m_s"])
        sflux = float(sent["electron_particle_line_flux_per_m_s"])
        native = native_by_edge.get(edge_id, 0.0)
        vdrop = float(vela["phin1_V"]) - float(vela["phin0_V"])
        sdrop = float(sent["phin1_V"]) - float(sent["phin0_V"])
        mobility = max(
            float(vela["electron_mobility_m2_V_s"]),
            float(sent["electron_mobility_m2_V_s"]),
        )
        native_relevant = native >= max(max_native * 1.0e-8, 1.0)
        state_relevant = max(abs(vflux), abs(sflux)) >= max(
            max_state_line_flux * 1.0e-8, 1.0
        )
        state_dex = signed_dex(vflux, sflux, 1.0e-30)
        native_dex = signed_dex(vflux, native, 1.0e-30)
        qf_dex = signed_dex(vdrop, sdrop, 1.0e-18)
        # A nodal Sentaurus current vector projected onto a Vela dual edge is
        # useful context but is not the same discrete support as an SG box
        # flux.  Rank anomalies using same-operator state replay and QF drops;
        # do not let the reconstructed native projection create false edges.
        abnormal = bias != 0.0 and active_sg and drain_endpoint_count < 2 and (
            mobility > 0.0 and state_relevant and
            (abs(state_dex) >= 1.0 or abs(qf_dex) >= 1.0)
        )
        compared.append(
            {
                "bias_V": bias,
                "edge_id": edge_id,
                "node0": int(vela["node0"]),
                "node1": int(vela["node1"]),
                "is_semiconductor_edge": int(is_silicon),
                "active_sg_edge": int(active_sg),
                "drain_contact_endpoint_count": drain_endpoint_count,
                "drain_hop": min(distances[int(vela["node0"])], distances[int(vela["node1"])]),
                "x_mid_um": 0.5 * (float(vela["x0"]) + float(vela["x1"])) * 1.0e6,
                "y_mid_um": 0.5 * (float(vela["y0"]) + float(vela["y1"])) * 1.0e6,
                "length_m": float(vela["length_m"]),
                "electron_mobility_m2_V_s": mobility,
                "vela_psi0_V": float(vela["psi0_V"]),
                "vela_psi1_V": float(vela["psi1_V"]),
                "sentaurus_psi0_V": float(sent["psi0_V"]),
                "sentaurus_psi1_V": float(sent["psi1_V"]),
                "vela_phin0_V": float(vela["phin0_V"]),
                "vela_phin1_V": float(vela["phin1_V"]),
                "sentaurus_phin0_V": float(sent["phin0_V"]),
                "sentaurus_phin1_V": float(sent["phin1_V"]),
                "vela_phin_drop_V": vdrop,
                "sentaurus_phin_drop_V": sdrop,
                "vela_over_sentaurus_qf_drop_dex": qf_dex,
                "vela_sg_particle_line_flux_per_m_s": vflux,
                "sentaurus_state_vela_sg_particle_line_flux_per_m_s": sflux,
                "sentaurus_native_edge_parallel_particle_line_flux_per_m_s": native,
                "sentaurus_native_projection_relevant": int(native_relevant),
                "state_flux_relevant": int(state_relevant),
                "vela_over_sentaurus_state_sg_flux_dex": state_dex,
                "vela_over_sentaurus_native_flux_dex": native_dex,
                "abnormal": int(abnormal),
            }
        )

    compared.sort(key=lambda row: (row["drain_hop"], row["edge_id"]))
    fields = list(compared[0])
    write_csv(output / f"edge_compare_{bias_tag(bias)}.csv", fields, compared)
    abnormal_rows = [row for row in compared if row["abnormal"]]
    first = min(
        abnormal_rows,
        key=lambda row: (
            row["drain_hop"],
            -abs(row["vela_over_sentaurus_native_flux_dex"]),
            row["edge_id"],
        ),
        default=None,
    )
    zero_couple_exit = next(
        (
            row for row in compared
            if row["is_semiconductor_edge"]
            and not row["active_sg_edge"]
            and row["drain_contact_endpoint_count"] == 1
        ),
        None,
    )
    return {
        "bias_V": bias,
        "edge_count": len(compared),
        "abnormal_edge_count": len(abnormal_rows),
        "first_abnormal_edge_id": first["edge_id"] if first else "",
        "first_abnormal_node0": first["node0"] if first else "",
        "first_abnormal_node1": first["node1"] if first else "",
        "first_abnormal_drain_hop": first["drain_hop"] if first else "",
        "first_abnormal_x_mid_um": first["x_mid_um"] if first else "",
        "first_abnormal_y_mid_um": first["y_mid_um"] if first else "",
        "first_abnormal_qf_drop_dex": first["vela_over_sentaurus_qf_drop_dex"] if first else "",
        "first_abnormal_state_sg_flux_dex": first["vela_over_sentaurus_state_sg_flux_dex"] if first else "",
        "first_abnormal_native_flux_dex": first["vela_over_sentaurus_native_flux_dex"] if first else "",
        "first_zero_couple_drain_exit_edge_id": zero_couple_exit["edge_id"] if zero_couple_exit else "",
        "first_zero_couple_drain_exit_node0": zero_couple_exit["node0"] if zero_couple_exit else "",
        "first_zero_couple_drain_exit_node1": zero_couple_exit["node1"] if zero_couple_exit else "",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-config", type=Path, required=True)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--vela-state-dir", type=Path, required=True)
    parser.add_argument("--sentaurus-imported-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--biases", default="0,0.001,0.002,0.005,0.01")
    parser.add_argument(
        "--runner",
        type=Path,
        default=REPO / "build-release/vela_example_runner.exe",
    )
    args = parser.parse_args()

    biases = [float(item) for item in args.biases.split(",") if item.strip()]
    base = json.loads(args.base_config.read_text(encoding="utf-8-sig"))
    mesh = json.loads(args.mesh.read_text(encoding="utf-8-sig"))
    count = len(mesh["nodes"])
    output = args.out_dir.resolve()
    summaries: list[dict[str, Any]] = []

    for bias in biases:
        tag = bias_tag(bias)
        vela_state = args.vela_state_dir / f"accepted_state_bias_{tag}.csv"
        sent_imported = args.sentaurus_imported_dir / f"lowbias_{tag}"
        vela_fields = output / "state_fields" / tag / "vela"
        sent_fields = output / "state_fields" / tag / "sentaurus"
        write_probe_fields(vela_fields, vela_state_fields(vela_state, count), count)
        write_probe_fields(sent_fields, merge_sentaurus_state(sent_imported, count), count)

        vela_edge = output / "edge_flux" / f"vela_{tag}.csv"
        sent_edge = output / "edge_flux" / f"sentaurus_state_{tag}.csv"
        run_probe(
            args.runner,
            base,
            bias,
            vela_fields,
            vela_edge,
            output / "configs" / f"vela_{tag}.json",
        )
        run_probe(
            args.runner,
            base,
            bias,
            sent_fields,
            sent_edge,
            output / "configs" / f"sentaurus_state_{tag}.json",
        )
        summaries.append(
            analyze_bias(bias, vela_edge, sent_edge, sent_imported, mesh, output)
        )

    write_csv(output / "first_abnormal_edge_summary.csv", list(summaries[0]), summaries)
    (output / "first_abnormal_edge_summary.json").write_text(
        json.dumps(summaries, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summaries, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
