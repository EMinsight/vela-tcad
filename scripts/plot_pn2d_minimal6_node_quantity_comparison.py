#!/usr/bin/env python3
"""Create six-node PN2D fixed-state quantity tables and PNG comparisons."""

from __future__ import annotations
import argparse, csv, hashlib, json, math
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image, ImageStat

TOPOLOGIES = ["sketch", "mirror"]
BIASES = [0.0, -12.0, -19.0]
NODES = list(range(1, 7))
INTENSIVE = [
    "vela_electron_impact_field_V_per_m", "vela_hole_impact_field_V_per_m",
    "vela_electron_alpha_per_m", "vela_hole_alpha_per_m",
    "vela_electron_flux_per_m2_s", "vela_hole_flux_per_m2_s",
]
SOURCES = ["vela_electron_edge_source_per_s", "vela_hole_edge_source_per_s"]


def aggregate_edges(rows: Iterable[dict[str, str]],
                    node_ids: set[int]) -> dict[int, dict[str, float]]:
    counts = {node: 0 for node in node_ids}
    sums = {node: {field: 0.0 for field in INTENSIVE} for node in node_ids}
    maxima = {node: {field: 0.0 for field in INTENSIVE} for node in node_ids}
    source = {node: {field: 0.0 for field in SOURCES} for node in node_ids}
    for row in rows:
        endpoints = [int(row["node0"]), int(row["node1"])]
        if any(node not in node_ids for node in endpoints):
            raise ValueError(f"edge references node outside 1..6: {endpoints}")
        for node in endpoints:
            counts[node] += 1
            for field in INTENSIVE:
                value = float(row[field])
                sums[node][field] += value
                maxima[node][field] = max(maxima[node][field], abs(value))
            for field in SOURCES:
                source[node][field] += 0.5 * float(row[field])
    result = {}
    for node in sorted(node_ids):
        count = counts[node]
        values: dict[str, float] = {"incident_edge_count": count}
        for field in INTENSIVE:
            stem = field.replace("vela_", "")
            if stem.endswith("_V_per_m"):
                stem = stem[:-8]
                mean_name = f"vela_{stem}_mean_V_per_m"
                max_name = f"vela_{stem}_max_abs_V_per_m"
            elif stem.endswith("_per_m2_s"):
                stem = stem[:-9]
                mean_name = f"vela_{stem}_mean_per_m2_s"
                max_name = f"vela_{stem}_max_abs_per_m2_s"
            elif stem.endswith("_per_m"):
                stem = stem[:-6]
                mean_name = f"vela_{stem}_mean_per_m"
                max_name = f"vela_{stem}_max_abs_per_m"
            values[mean_name] = sums[node][field] / count if count else 0.0
            values[max_name] = maxima[node][field]
        values["vela_electron_source_node_per_s"] = source[node][SOURCES[0]]
        values["vela_hole_source_node_per_s"] = source[node][SOURCES[1]]
        values["vela_total_source_node_per_s"] = sum(source[node].values())
        result[node] = values
    return result


def validate_node_groups(rows: list[dict[str, Any]]) -> None:
    groups: dict[tuple[str, float], set[int]] = defaultdict(set)
    for row in rows:
        groups[(str(row["topology_id"]), float(row["bias_V"]))].add(int(row["node_id"]))
    expected_nodes = set(NODES)
    for key, nodes in groups.items():
        if nodes != expected_nodes:
            raise ValueError(f"{key} must contain exactly nodes 1..6, got {sorted(nodes)}")
    expected_keys = {(topology, bias) for topology in TOPOLOGIES for bias in BIASES}
    if set(groups) != expected_keys:
        raise ValueError(f"topology/bias groups differ: {sorted(groups)}")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def build_rows(node_rows: list[dict[str, str]],
               edge_rows: list[dict[str, str]]) -> list[dict[str, Any]]:
    validate_node_groups(node_rows)
    edges: dict[tuple[str, float], list[dict[str, str]]] = defaultdict(list)
    for row in edge_rows:
        edges[(row["topology_id"], float(row["bias_V"]))].append(row)
    output = []
    for topology in TOPOLOGIES:
        for bias in BIASES:
            agg = aggregate_edges(edges[(topology, bias)], set(NODES))
            selected = [row for row in node_rows
                        if row["topology_id"] == topology and float(row["bias_V"]) == bias]
            selected.sort(key=lambda row: int(row["node_id"]))
            for row in selected:
                node = int(row["node_id"])
                ex = float(row["sentaurus_electric_field_x_V_per_m"])
                ey = float(row["sentaurus_electric_field_y_V_per_m"])
                item = {
                    "topology_id": topology, "bias_V": bias, "node_id": node,
                    "x_um": float(row["x_um"]), "y_um": float(row["y_um"]),
                    "sentaurus_psi_V": float(row["sentaurus_psi_V"]),
                    "sentaurus_phin_V": float(row["sentaurus_phin_V"]),
                    "sentaurus_phip_V": float(row["sentaurus_phip_V"]),
                    "sentaurus_n_m3": float(row["sentaurus_n_m3"]),
                    "sentaurus_p_m3": float(row["sentaurus_p_m3"]),
                    "vela_psi_V": float(row["vela_psi_V"]),
                    "vela_phin_V": float(row["vela_phin_V"]),
                    "vela_phip_V": float(row["vela_phip_V"]),
                    "vela_n_m3": float(row["vela_n_m3"]),
                    "vela_p_m3": float(row["vela_p_m3"]),
                    "sentaurus_electric_field_magnitude_V_per_m": math.hypot(ex, ey),
                    "sentaurus_electron_mobility_m2_per_V_s": float(row["sentaurus_electron_mobility_m2_per_V_s"]),
                    "sentaurus_hole_mobility_m2_per_V_s": float(row["sentaurus_hole_mobility_m2_per_V_s"]),
                    "sentaurus_electron_alpha_per_m": float(row["sentaurus_electron_alpha_per_m"]),
                    "sentaurus_hole_alpha_per_m": float(row["sentaurus_hole_alpha_per_m"]),
                }
                item.update(agg[node]); output.append(item)
    if len(output) != 36:
        raise ValueError(f"expected 36 rows, got {len(output)}")
    return output


def transform(quantity: str, values: list[float]) -> list[float]:
    if quantity.endswith(("_V", "_V_per_m")) and "alpha" not in quantity:
        return values
    if "flux" in quantity:
        return [math.copysign(math.log10(1 + abs(v)), v) if v else 0.0 for v in values]
    return [math.log10(max(abs(v), 1e-300)) for v in values]


def plot_grid(rows: list[dict[str, Any]], quantities: list[str],
              path: Path, title: str) -> None:
    indexed = {(r["topology_id"], float(r["bias_V"]), int(r["node_id"])): r for r in rows}
    colors = {0.0: "#374151", -12.0: "#2563eb", -19.0: "#d97706"}
    styles = {0.0: ":", -12.0: "-", -19.0: "--"}
    fig, axes = plt.subplots(len(quantities), 2, figsize=(16, max(10, 3*len(quantities))),
                             squeeze=False, constrained_layout=True)
    for qi, quantity in enumerate(quantities):
        for ti, topology in enumerate(TOPOLOGIES):
            ax = axes[qi, ti]
            for bias in BIASES:
                values = [float(indexed[(topology,bias,node)][quantity]) for node in NODES]
                ax.plot(NODES, transform(quantity, values), marker="o",
                        color=colors[bias], linestyle=styles[bias], label=f"{bias:g} V")
            ax.set_title(f"{topology}: {quantity}", fontsize=10)
            ax.set_xlabel("Canonical node ID"); ax.set_xticks(NODES)
            ax.grid(True, color="#d1d5db", linewidth=.6)
            if qi == 0: ax.legend(ncol=3)
    fig.suptitle(title, fontsize=16)
    fig.savefig(path, dpi=180, facecolor="white"); plt.close(fig)


def sha256(path: Path) -> str:
    d=hashlib.sha256()
    with path.open("rb") as h:
        for chunk in iter(lambda:h.read(1024*1024),b""): d.update(chunk)
    return d.hexdigest()


def validate_png(path: Path) -> dict[str, Any]:
    im=Image.open(path).convert("RGB"); stat=ImageStat.Stat(im)
    if im.width<1600 or im.height<900 or max(stat.var)<=100:
        raise ValueError(f"PNG QA failed: {path}")
    return {"path":str(path.resolve()),"width_px":im.width,"height_px":im.height,
            "max_pixel_variance":max(stat.var),"size_bytes":path.stat().st_size}


def parse_args():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--audit-root",type=Path,required=True)
    p.add_argument("--out-dir",type=Path,required=True)
    return p.parse_args()


def main():
    a=parse_args(); node_path=a.audit_root/"node_state.csv"; edge_path=a.audit_root/"edge_audit.csv"
    rows=build_rows(read_csv(node_path),read_csv(edge_path))
    a.out_dir.mkdir(parents=True,exist_ok=True)
    write_csv(a.out_dir/"minimal6_all_nodes.csv",rows)
    table_names=["minimal6_all_nodes.csv"]
    for bias,token in [(0.0,"0V"),(-12.0,"minus12V"),(-19.0,"minus19V")]:
        name=f"minimal6_nodes_{token}.csv"
        write_csv(a.out_dir/name,[row for row in rows if row["bias_V"]==bias])
        table_names.append(name)
    state=["sentaurus_psi_V","sentaurus_phin_V","sentaurus_phip_V",
           "sentaurus_n_m3","sentaurus_p_m3"]
    ion=["sentaurus_electric_field_magnitude_V_per_m",
         "sentaurus_electron_alpha_per_m","sentaurus_hole_alpha_per_m",
         "vela_electron_impact_field_mean_V_per_m","vela_hole_impact_field_mean_V_per_m",
         "vela_electron_alpha_mean_per_m","vela_hole_alpha_mean_per_m"]
    transport=["vela_electron_flux_mean_per_m2_s","vela_hole_flux_mean_per_m2_s",
               "vela_electron_source_node_per_s","vela_hole_source_node_per_s",
               "vela_total_source_node_per_s"]
    pngs=["minimal6_state_by_node.png","minimal6_ionization_by_node.png",
          "minimal6_source_flux_by_node.png"]
    plot_grid(rows,state,a.out_dir/pngs[0],"PN2D minimal6 state quantities: 6 nodes, 4 triangles")
    plot_grid(rows,ion,a.out_dir/pngs[1],"PN2D minimal6 collision-ionization quantities")
    plot_grid(rows,transport,a.out_dir/pngs[2],"PN2D minimal6 flux and avalanche source")
    qa=[validate_png(a.out_dir/name) for name in pngs]
    manifest={"schema":"vela.pn2d_minimal6_node_quantity_comparison.v1",
              "audit_root":str(a.audit_root.resolve()),"topologies":TOPOLOGIES,
              "biases_V":BIASES,"node_ids":NODES,"node_count_per_state":6,
              "triangle_count_per_topology":4,"rows":len(rows),
              "source_sha256":{"node_state.csv":sha256(node_path),"edge_audit.csv":sha256(edge_path)},
              "tables":table_names,"png_qa":qa,
              "edge_to_node":{"source":"half to each endpoint then sum",
                              "intensive":"incident-edge mean and max-abs",
                              "interpolation":False}}
    (a.out_dir/"minimal6_node_quantity_manifest.json").write_text(
        json.dumps(manifest,indent=2,sort_keys=True)+"\n",encoding="utf-8")
    print(json.dumps(manifest,sort_keys=True)); return 0

if __name__=="__main__": raise SystemExit(main())
