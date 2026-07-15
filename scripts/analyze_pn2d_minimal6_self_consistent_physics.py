#!/usr/bin/env python3
"""Compare self-consistent Minimal6 Vela states with Sentaurus exports."""
from __future__ import annotations
import argparse, csv, json, math
from pathlib import Path
from typing import Any

BIASES = (-1, -12, -19, -20)
TOPOLOGIES = ("sketch", "mirror")
VT_300_K = 0.025851999786

def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))

def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

def scalar_field(export: Path, name: str) -> dict[int, float]:
    return {int(r["node_id"]): float(r["component0"])
            for r in read_csv(export / "fields" / f"{name}_region0.csv")}

def vector_magnitude_field(export: Path, name: str) -> dict[int, float]:
    return {int(r["node_id"]): math.hypot(float(r["component0"]), float(r["component1"]))
            for r in read_csv(export / "fields" / f"{name}_region0.csv")}

def state_rows(path: Path, *, sentaurus: bool) -> dict[int, dict[str, float]]:
    result = {}
    for row in read_csv(path):
        node = int(row["node_id"])
        if sentaurus:
            result[node] = {
                "psi": float(row["psi_V"]), "phin": float(row["phin_V"]),
                "phip": float(row["phip_V"]), "n_m3": float(row["n_m3"]),
                "p_m3": float(row["p_m3"])}
        else:
            # Unit-scaling serialization currently puts internal cm^-3 values
            # under *_m3 headers, so convert them numerically before comparison.
            result[node] = {
                "psi": float(row["psi"]), "phin": float(row["phin"]),
                "phip": float(row["phip"]),
                "n_m3": float(row["electrons_m3"]) * 1.0e6,
                "p_m3": float(row["holes_m3"]) * 1.0e6}
    return result

def sentaurus_state_from_export(export: Path) -> dict[int, dict[str, float]]:
    psi = scalar_field(export, "ElectrostaticPotential")
    phin = scalar_field(export, "eQuasiFermiPotential")
    phip = scalar_field(export, "hQuasiFermiPotential")
    electrons = scalar_field(export, "eDensity")
    holes = scalar_field(export, "hDensity")
    return {node: {"psi": psi[node], "phin": phin[node], "phip": phip[node],
                   "n_m3": electrons[node] * 1.0e6, "p_m3": holes[node] * 1.0e6}
            for node in psi}

def mesh_edges(mesh: dict[str, Any]) -> list[tuple[int, int, float]]:
    xy = {int(n["id"]): (float(n["x"]), float(n["y"])) for n in mesh["nodes"]}
    unique = set()
    for tri in mesh["triangles"]:
        ids = [int(v) for v in tri["node_ids"]]
        unique.update(tuple(sorted(pair)) for pair in
                      ((ids[0], ids[1]), (ids[1], ids[2]), (ids[2], ids[0])))
    return [(a, b, math.hypot(xy[b][0]-xy[a][0], xy[b][1]-xy[a][1]) * 1.0e-6)
            for a, b in sorted(unique)]

def max_edge_gradient(state, edges, key):
    return max(abs(state[b][key] - state[a][key]) / length for a, b, length in edges)

def geomean(values):
    return math.exp(sum(math.log(value) for value in values) / len(values))

def endpoint_file(root, topology, bias):
    return root / "vela" / topology / "states" / f"segment_{abs(bias)-1:02d}_bias_m{abs(bias)}p000000.csv"

def diagnostic_file(root, topology, bias):
    return root / "vela" / topology / "diagnostics" / f"segment_{abs(bias)-1:02d}_sg_avalanche_edges.csv"

def sentaurus_export(root, topology, bias):
    return root / "sentaurus" / topology / "exports" / f"{topology}_m{abs(bias)}p000000"

def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--vela-root", type=Path, required=True)
    p.add_argument("--sentaurus-root", type=Path, required=True)
    p.add_argument("--comparison-csv", type=Path, required=True)
    p.add_argument("--out-dir", type=Path, required=True)
    a = p.parse_args()
    a.out_dir.mkdir(parents=True, exist_ok=True)
    ratios = {(r["topology"], int(float(r["bias_V"]))): r for r in read_csv(a.comparison_csv)
              if r["classification"] == "common_exact"}
    summaries, nodes_out = [], []
    for topology in TOPOLOGIES:
        mesh = json.loads((a.vela_root / "inputs" / topology / "mesh.json").read_text(encoding="utf-8"))
        edges = mesh_edges(mesh)
        contacts = {int(n) for c in mesh["contacts"] for n in c["node_ids"]}
        interior = sorted(set(range(len(mesh["nodes"]))) - contacts)
        for bias in BIASES:
            export = sentaurus_export(a.sentaurus_root, topology, bias)
            vela = state_rows(endpoint_file(a.vela_root, topology, bias), sentaurus=False)
            sent = sentaurus_state_from_export(export)
            sent_field = vector_magnitude_field(export, "ElectricField")
            sent_e_alpha = scalar_field(export, "eAlphaAvalanche")
            sent_h_alpha = scalar_field(export, "hAlphaAvalanche")
            edge_rows = [r for r in read_csv(diagnostic_file(a.vela_root, topology, bias))
                         if abs(float(r["bias_V"]) - bias) <= 1.0e-10]
            vela_field = max(float(r["electric_field_V_per_m"]) for r in edge_rows)
            vela_e_alpha = max(float(r["electron_alpha_m_inv"]) for r in edge_rows) / 100.0
            vela_h_alpha = max(float(r["hole_alpha_m_inv"]) for r in edge_rows) / 100.0
            for node in sorted(vela):
                v, s = vela[node], sent[node]
                nodes_out.append({
                    "topology": topology, "bias_V": bias, "node_id": node,
                    "is_contact": node in contacts,
                    "vela_psi_V": v["psi"], "sentaurus_psi_V": s["psi"],
                    "psi_error_V": v["psi"]-s["psi"],
                    "vela_phin_V": v["phin"], "sentaurus_phin_V": s["phin"],
                    "phin_error_V": v["phin"]-s["phin"],
                    "vela_phip_V": v["phip"], "sentaurus_phip_V": s["phip"],
                    "phip_error_V": v["phip"]-s["phip"],
                    "vela_n_m3": v["n_m3"], "sentaurus_n_m3": s["n_m3"],
                    "n_ratio": v["n_m3"]/s["n_m3"],
                    "vela_p_m3": v["p_m3"], "sentaurus_p_m3": s["p_m3"],
                    "p_ratio": v["p_m3"]/s["p_m3"]})
            current = ratios[(topology, bias)]
            phin_err = sum(abs(vela[n]["phin"]-sent[n]["phin"]) for n in interior)/len(interior)
            phip_err = sum(abs(vela[n]["phip"]-sent[n]["phip"]) for n in interior)/len(interior)
            sent_max_field = max(sent_field.values()) * 100.0
            summaries.append({
                "topology": topology, "bias_V": bias,
                "max_abs_psi_error_all_V": max(abs(vela[n]["psi"]-sent[n]["psi"]) for n in vela),
                "max_abs_psi_error_interior_V": max(abs(vela[n]["psi"]-sent[n]["psi"]) for n in interior),
                "mean_abs_phin_error_interior_V": phin_err,
                "mean_abs_phip_error_interior_V": phip_err,
                "vela_mean_qf_split_interior_V": sum(abs(vela[n]["phip"]-vela[n]["phin"]) for n in interior)/len(interior),
                "sentaurus_mean_qf_split_interior_V": sum(abs(sent[n]["phip"]-sent[n]["phin"]) for n in interior)/len(interior),
                "n_ratio_geomean_interior": geomean([vela[n]["n_m3"]/sent[n]["n_m3"] for n in interior]),
                "p_ratio_geomean_interior": geomean([vela[n]["p_m3"]/sent[n]["p_m3"] for n in interior]),
                "boltzmann_factor_from_mean_phin_error": math.exp(phin_err/VT_300_K),
                "vela_max_electric_field_V_per_m": vela_field,
                "sentaurus_max_electric_field_V_per_m": sent_max_field,
                "electric_field_ratio": vela_field/sent_max_field,
                "vela_max_abs_phin_gradient_V_per_m": max_edge_gradient(vela, edges, "phin"),
                "sentaurus_max_abs_phin_gradient_V_per_m": max_edge_gradient(sent, edges, "phin"),
                "vela_max_abs_phip_gradient_V_per_m": max_edge_gradient(vela, edges, "phip"),
                "sentaurus_max_abs_phip_gradient_V_per_m": max_edge_gradient(sent, edges, "phip"),
                "vela_max_electron_alpha_cm_inv": vela_e_alpha,
                "sentaurus_max_electron_alpha_cm_inv": max(sent_e_alpha.values()),
                "electron_alpha_ratio": vela_e_alpha/max(sent_e_alpha.values()),
                "vela_max_hole_alpha_cm_inv": vela_h_alpha,
                "sentaurus_max_hole_alpha_cm_inv": max(sent_h_alpha.values()),
                "hole_alpha_ratio": vela_h_alpha/max(sent_h_alpha.values()),
                "terminal_current_ratio": float(current["terminal_current_ratio"]),
                "source_integral_ratio": float(current["native_source_ratio"])})
    write_csv(a.out_dir / "self_consistent_node_comparison.csv", nodes_out)
    write_csv(a.out_dir / "self_consistent_summary.csv", summaries)
    payload = {
        "schema": "vela.pn2d_minimal6_self_consistent_physics_error.v1",
        "biases_V": list(BIASES), "topologies": list(TOPOLOGIES),
        "inputs": {
            "vela_state_root": str(a.vela_root),
            "sentaurus_full_field_export_root": str(a.sentaurus_root),
            "corrected_sweep_comparison_csv": str(a.comparison_csv),
        },
        "density_unit_note": "Vela unit_scaling state headers say m3 but values are cm^-3; multiplied by 1e6 before comparison.",
        "formula_audit": {
            "carrier_statistics": "n=ni*exp((psi-phin)/Vt); p=ni*exp((phip-psi)/Vt)",
            "electric_field": "E=-grad(psi)",
            "impact_coefficient": "alpha=gamma*A*exp(-gamma*B/abs(F)) with configured low/high branch",
            "avalanche_source": "G=(alpha_n*abs(Jn)+alpha_p*abs(Jp))/q",
            "tcad_internal_transport_geometry_cm_per_um": 1.0e-4,
            "tcad_internal_source_geometry_cm2_per_um2": 1.0e-8,
            "raw_source_relative_overweight_in_continuity": 1.0e4,
            "task8_source_diagnostic_area_conversion": 1.0e-8},
        "summary_rows": summaries}
    (a.out_dir / "self_consistent_analysis.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
