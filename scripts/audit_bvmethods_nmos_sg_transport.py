#!/usr/bin/env python3
"""Audit the 6.4 V BVmethods NMOS edge-current support term by term."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any


Q_C = 1.602176634e-19
KB_J_K = 1.380649e-23
TEMPERATURE_K = 300.0
VT_V = KB_J_K * TEMPERATURE_K / Q_C

REPO = Path(__file__).resolve().parents[1]
RUN_ROOT = REPO / "build-release/reference_tcad/bvmethods_sentaurus2018/run01"
DEFAULT_SENT = RUN_ROOT / "sentaurus_iic_multibias_exact_extended_20260803/imported/iic_v6p400000"
DEFAULT_VELA = RUN_ROOT / "vela_validation/iic_rebuild_fd_gummel_20260803/probe_6p4_full/postprocess_only"
DEFAULT_OUT = RUN_ROOT / "vela_validation/iic_rebuild_fd_gummel_20260803/transport_audit_6p4"

DEFAULT_BRANCH_STATES = {
    1.0: RUN_ROOT / "vela_validation/iic_rebuild_20260803/trunk_0p95_1p0/postprocess_only/states/accepted_state_bias_1p000000.csv",
    2.0: RUN_ROOT / "vela_validation/iic_rebuild_20260803/trunk_1p9_2p0_earlyfloor/postprocess_only/states/accepted_state_bias_2p000000.csv",
    4.0: RUN_ROOT / "vela_validation/iic_rebuild_20260803/trunk_3p1_4p0_earlyfloor/postprocess_only/states/accepted_state_bias_4p000000.csv",
    5.0: RUN_ROOT / "vela_validation/iic_rebuild_20260803/trunk_4p55_5p0_earlyfloor/postprocess_only/states/accepted_state_bias_5p000000.csv",
    6.0: RUN_ROOT / "vela_validation/iic_rebuild_fd_gummel_20260803/trunk_5p8_6p4_newton/postprocess_only/states/accepted_state_bias_6p000000.csv",
    6.4: RUN_ROOT / "vela_validation/iic_rebuild_fd_gummel_20260803/probe_6p4_full/postprocess_only/states/accepted_state_bias_6p400000.csv",
}


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def scalar_field(root: Path, name: str) -> dict[int, float]:
    return {
        int(row["node_id"]): float(row["component0"])
        for row in read_rows(root / "fields" / f"{name}_region3.csv")
    }


def sent_tag(bias: float) -> str:
    return f"iic_v{bias:.6f}".replace(".", "p")


def vector_field(root: Path, name: str) -> dict[int, tuple[float, float]]:
    return {
        int(row["node_id"]): (float(row["component0"]), float(row["component1"]))
        for row in read_rows(root / "fields" / f"{name}_region3.csv")
    }


def log_mean(left: float, right: float) -> float:
    if left <= 0.0 or right <= 0.0:
        return 0.0
    delta = math.log(right) - math.log(left)
    if abs(delta) < 1.0e-12:
        return 0.5 * (left + right)
    return (right - left) / delta


def generalized_factor(density0: float, density1: float, eta_delta: float) -> float:
    if density0 <= 0.0 or density1 <= 0.0:
        return 1.0
    log_density_delta = math.log(density1 / density0)
    if abs(log_density_delta) <= 1.0e-10:
        return 1.0
    factor = eta_delta / log_density_delta
    return factor if math.isfinite(factor) and factor > 0.0 else 1.0


def projected_edge_value(
    vectors: dict[int, tuple[float, float]], row: dict[str, str]
) -> float:
    node0 = int(row["node0"])
    node1 = int(row["node1"])
    dx = float(row["x1_um"]) - float(row["x0_um"])
    dy = float(row["y1_um"]) - float(row["y0_um"])
    length = math.hypot(dx, dy)
    if length <= 0.0:
        return 0.0
    vx = 0.5 * (vectors[node0][0] + vectors[node1][0])
    vy = 0.5 * (vectors[node0][1] + vectors[node1][1])
    return abs((vx * dx + vy * dy) / length)


def integrate_semiconductor_generation_A_per_um(
    mesh: dict[str, Any], generation_cm3_s: dict[int, float]
) -> float:
    nodes = {
        int(node["id"]): (float(node["x"]), float(node["y"]))
        for node in mesh["nodes"]
    }
    integral_per_cm_s = 0.0
    for cell in mesh["triangles"]:
        if int(cell["region_id"]) != 3 or len(cell["node_ids"]) != 3:
            continue
        node_ids = [int(node) for node in cell["node_ids"]]
        if any(node not in generation_cm3_s for node in node_ids):
            continue
        (x0, y0), (x1, y1), (x2, y2) = (nodes[node] for node in node_ids)
        area_um2 = 0.5 * abs(
            (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)
        )
        generation_average = sum(generation_cm3_s[node] for node in node_ids) / 3.0
        # cm^-3 s^-1 * um^2 * 1 um depth = cm^-3 s^-1 * 1e-12 cm^3.
        integral_per_cm_s += generation_average * area_um2 * 1.0e-12
    return Q_C * integral_per_cm_s


def ratio(numerator: float, denominator: float) -> float:
    if denominator == 0.0:
        return math.inf if numerator != 0.0 else 1.0
    return numerator / denominator


def percentile(values: list[float], fraction: float) -> float:
    finite = sorted(value for value in values if math.isfinite(value))
    if not finite:
        return math.nan
    position = fraction * (len(finite) - 1)
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return finite[lower]
    weight = position - lower
    return finite[lower] * (1.0 - weight) + finite[upper] * weight


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sentaurus-root", type=Path, default=DEFAULT_SENT)
    parser.add_argument("--vela-root", type=Path, default=DEFAULT_VELA)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--top-n", type=int, default=30)
    args = parser.parse_args()

    sent_psi = scalar_field(args.sentaurus_root, "ElectrostaticPotential")
    sent_phin = scalar_field(args.sentaurus_root, "eQuasiFermiPotential")
    sent_phip = scalar_field(args.sentaurus_root, "hQuasiFermiPotential")
    # Sentaurus TDR stores carrier density in cm^-3; Vela accepted states use m^-3.
    sent_n = {
        node: value * 1.0e6
        for node, value in scalar_field(args.sentaurus_root, "eDensity").items()
    }
    sent_p = {
        node: value * 1.0e6
        for node, value in scalar_field(args.sentaurus_root, "hDensity").items()
    }
    sent_mun = scalar_field(args.sentaurus_root, "eMobility")
    sent_mup = scalar_field(args.sentaurus_root, "hMobility")
    sent_jn = vector_field(args.sentaurus_root, "eCurrentDensity")
    sent_jp = vector_field(args.sentaurus_root, "hCurrentDensity")
    vela_state = {
        int(row["node_id"]): row
        for row in read_rows(args.vela_root / "states" / "accepted_state_bias_6p400000.csv")
    }

    audit_rows: list[dict[str, Any]] = []
    for edge in read_rows(args.vela_root / "sg_avalanche_edges.csv"):
        node0 = int(edge["node0"])
        node1 = int(edge["node1"])
        if node0 not in sent_psi or node1 not in sent_psi:
            continue
        state0 = vela_state[node0]
        state1 = vela_state[node1]
        length_m = float(edge["edge_length_m"])
        if length_m <= 0.0:
            continue

        ni0 = float(edge["electron_sg_ni0"])
        ni1 = float(edge["electron_sg_ni1"])
        sent_dphin = sent_phin[node1] - sent_phin[node0]
        vela_dphin = float(state1["phin"]) - float(state0["phin"])
        sent_dphip = sent_phip[node1] - sent_phip[node0]
        vela_dphip = float(state1["phip"]) - float(state0["phip"])
        sent_dpsi = sent_psi[node1] - sent_psi[node0]
        vela_dpsi = float(state1["psi"]) - float(state0["psi"])

        sent_n_mid = log_mean(sent_n[node0], sent_n[node1])
        vela_n_mid = log_mean(float(state0["electrons_m3"]), float(state1["electrons_m3"]))
        sent_p_mid = log_mean(sent_p[node0], sent_p[node1])
        vela_p_mid = log_mean(float(state0["holes_m3"]), float(state1["holes_m3"]))
        sent_mun_edge = 0.5 * (sent_mun[node0] + sent_mun[node1]) * 1.0e-4
        sent_mup_edge = 0.5 * (sent_mup[node0] + sent_mup[node1]) * 1.0e-4
        vela_mun_edge = float(edge["electron_mobility_m2_V_s"])
        vela_mup_edge = float(edge["hole_mobility_m2_V_s"])

        sent_eta_n_delta = (
            (sent_dpsi - sent_dphin) / VT_V + math.log(ni1 / ni0)
        )
        vela_eta_n_delta = (
            (vela_dpsi - vela_dphin) / VT_V + math.log(ni1 / ni0)
        )
        sent_eta_p_delta = (
            (sent_dphip - sent_dpsi) / VT_V + math.log(ni1 / ni0)
        )
        vela_eta_p_delta = (
            (vela_dphip - vela_dpsi) / VT_V + math.log(ni1 / ni0)
        )
        sent_gn = generalized_factor(sent_n[node0], sent_n[node1], sent_eta_n_delta)
        vela_gn = generalized_factor(
            float(state0["electrons_m3"]), float(state1["electrons_m3"]),
            vela_eta_n_delta,
        )
        sent_gp = generalized_factor(sent_p[node0], sent_p[node1], sent_eta_p_delta)
        vela_gp = generalized_factor(
            float(state0["holes_m3"]), float(state1["holes_m3"]),
            vela_eta_p_delta,
        )
        sent_un = (sent_dpsi + VT_V * math.log(ni1 / ni0)) / (VT_V * sent_gn)
        vela_un = (vela_dpsi + VT_V * math.log(ni1 / ni0)) / (VT_V * vela_gn)
        sent_up = (sent_dpsi + VT_V * math.log(ni0 / ni1)) / (VT_V * sent_gp)
        vela_up = (vela_dpsi + VT_V * math.log(ni0 / ni1)) / (VT_V * vela_gp)

        sent_jn_proxy = Q_C * sent_mun_edge * sent_n_mid * abs(sent_dphin) / length_m
        vela_jn_proxy = Q_C * vela_mun_edge * vela_n_mid * abs(vela_dphin) / length_m
        sent_jp_proxy = Q_C * sent_mup_edge * sent_p_mid * abs(sent_dphip) / length_m
        vela_jp_proxy = Q_C * vela_mup_edge * vela_p_mid * abs(vela_dphip) / length_m
        sent_jn_edge = projected_edge_value(sent_jn, edge) * 1.0e4
        sent_jp_edge = projected_edge_value(sent_jp, edge) * 1.0e4
        vela_jn_edge = abs(float(
            edge["electron_sg_production_signed_conventional_current_density_A_per_m2"]
        ))
        vela_jp_edge = Q_C * 1.0e4 * abs(float(edge["hole_flux_proxy"]))

        audit_rows.append({
            "edge_id": int(edge["edge_id"]),
            "node0": node0,
            "node1": node1,
            "x_mid_um": 0.5 * (float(edge["x0_um"]) + float(edge["x1_um"])),
            "y_mid_um": 0.5 * (float(edge["y0_um"]) + float(edge["y1_um"])),
            "edge_length_m": length_m,
            "edge_couple_m": float(edge["edge_couple_m"]),
            "geometry_ratio_vela_over_sentaurus": 1.0,
            "sentaurus_electron_current_density_A_m2": sent_jn_edge,
            "vela_electron_current_density_A_m2": vela_jn_edge,
            "electron_current_ratio_vela_over_sentaurus": ratio(vela_jn_edge, sent_jn_edge),
            "sentaurus_electron_qf_drop_V": abs(sent_dphin),
            "vela_electron_qf_drop_V": abs(vela_dphin),
            "electron_qf_drop_ratio": ratio(abs(vela_dphin), abs(sent_dphin)),
            "sentaurus_electron_density_logmean_m3": sent_n_mid,
            "vela_electron_density_logmean_m3": vela_n_mid,
            "electron_density_ratio": ratio(vela_n_mid, sent_n_mid),
            "sentaurus_electron_mobility_m2_V_s": sent_mun_edge,
            "vela_electron_mobility_m2_V_s": vela_mun_edge,
            "electron_mobility_ratio": ratio(vela_mun_edge, sent_mun_edge),
            "sentaurus_electron_generalized_einstein": sent_gn,
            "vela_electron_generalized_einstein": vela_gn,
            "electron_generalized_einstein_ratio": ratio(vela_gn, sent_gn),
            "sentaurus_electron_bernoulli_argument": sent_un,
            "vela_electron_bernoulli_argument": vela_un,
            "sentaurus_electron_qf_proxy_A_m2": sent_jn_proxy,
            "vela_electron_qf_proxy_A_m2": vela_jn_proxy,
            "sentaurus_hole_current_density_A_m2": sent_jp_edge,
            "vela_hole_current_density_A_m2": vela_jp_edge,
            "hole_current_ratio_vela_over_sentaurus": ratio(vela_jp_edge, sent_jp_edge),
            "sentaurus_hole_qf_drop_V": abs(sent_dphip),
            "vela_hole_qf_drop_V": abs(vela_dphip),
            "hole_qf_drop_ratio": ratio(abs(vela_dphip), abs(sent_dphip)),
            "sentaurus_hole_density_logmean_m3": sent_p_mid,
            "vela_hole_density_logmean_m3": vela_p_mid,
            "hole_density_ratio": ratio(vela_p_mid, sent_p_mid),
            "sentaurus_hole_mobility_m2_V_s": sent_mup_edge,
            "vela_hole_mobility_m2_V_s": vela_mup_edge,
            "hole_mobility_ratio": ratio(vela_mup_edge, sent_mup_edge),
            "sentaurus_hole_generalized_einstein": sent_gp,
            "vela_hole_generalized_einstein": vela_gp,
            "hole_generalized_einstein_ratio": ratio(vela_gp, sent_gp),
            "sentaurus_hole_bernoulli_argument": sent_up,
            "vela_hole_bernoulli_argument": vela_up,
            "sentaurus_hole_qf_proxy_A_m2": sent_jp_proxy,
            "vela_hole_qf_proxy_A_m2": vela_jp_proxy,
        })

    audit_rows.sort(
        key=lambda row: row["sentaurus_electron_current_density_A_m2"], reverse=True
    )
    write_csv(args.out_dir / "edge_transport_audit.csv", audit_rows)
    top = audit_rows[: args.top_n]
    write_csv(args.out_dir / "top_sentaurus_electron_current_edges.csv", top)

    summary_values = {
        "electron_current_ratio": [row["electron_current_ratio_vela_over_sentaurus"] for row in top],
        "electron_qf_drop_ratio": [row["electron_qf_drop_ratio"] for row in top],
        "electron_density_ratio": [row["electron_density_ratio"] for row in top],
        "electron_mobility_ratio": [row["electron_mobility_ratio"] for row in top],
        "electron_generalized_einstein_ratio": [
            row["electron_generalized_einstein_ratio"] for row in top
        ],
    }
    summary_rows = [
        {
            "population": f"top_{len(top)}_sentaurus_electron_current_edges",
            "quantity": name,
            "p50": percentile(values, 0.5),
            "p05": percentile(values, 0.05),
            "p95": percentile(values, 0.95),
        }
        for name, values in summary_values.items()
    ]
    write_csv(args.out_dir / "component_ratio_summary.csv", summary_rows)

    hotspot = top[0]
    hotspot_edge = next(
        edge for edge in read_rows(args.vela_root / "sg_avalanche_edges.csv")
        if int(edge["edge_id"]) == int(hotspot["edge_id"])
    )
    hotspot_evolution: list[dict[str, Any]] = []
    sentaurus_import_root = args.sentaurus_root.parent
    mesh = json.loads((RUN_ROOT / "vela/mesh.json").read_text(encoding="utf-8"))
    for bias, state_path in DEFAULT_BRANCH_STATES.items():
        sent_root = sentaurus_import_root / sent_tag(bias)
        if not state_path.exists() or not sent_root.exists():
            continue
        node0 = int(hotspot_edge["node0"])
        node1 = int(hotspot_edge["node1"])
        state = {int(row["node_id"]): row for row in read_rows(state_path)}
        psi_ref = scalar_field(sent_root, "ElectrostaticPotential")
        phin_ref = scalar_field(sent_root, "eQuasiFermiPotential")
        n_ref = {
            node: value * 1.0e6
            for node, value in scalar_field(sent_root, "eDensity").items()
        }
        current_ref = vector_field(sent_root, "eCurrentDensity")
        btbt_ref = scalar_field(sent_root, "Band2BandGeneration")
        sent_qf_drop = abs(phin_ref[node1] - phin_ref[node0])
        vela_qf_drop = abs(float(state[node1]["phin"]) - float(state[node0]["phin"]))
        sent_n_mid = log_mean(n_ref[node0], n_ref[node1])
        vela_n_mid = log_mean(
            float(state[node0]["electrons_m3"]),
            float(state[node1]["electrons_m3"]),
        )
        sent_qf_average = 0.5 * (phin_ref[node0] + phin_ref[node1])
        vela_qf_average = 0.5 * (
            float(state[node0]["phin"]) + float(state[node1]["phin"])
        )
        sent_psi_average = 0.5 * (psi_ref[node0] + psi_ref[node1])
        vela_psi_average = 0.5 * (
            float(state[node0]["psi"]) + float(state[node1]["psi"])
        )
        sent_drain_current = float(read_rows(
            sent_root / "fields" / "ContactCurrentFlux_region8.csv"
        )[0]["component0"])
        vela_sweep_rows = read_rows(state_path.parent.parent / "sweep.csv")
        vela_drain_current = float(vela_sweep_rows[-1]["current_total_A_per_um"])
        hotspot_evolution.append({
            "bias_V": bias,
            "edge_id": int(hotspot_edge["edge_id"]),
            "sentaurus_drain_current_A_per_um": sent_drain_current,
            "vela_drain_current_A_per_um": vela_drain_current,
            "drain_current_ratio_vela_over_sentaurus": ratio(
                abs(vela_drain_current), abs(sent_drain_current)
            ),
            "sentaurus_edge_projected_electron_current_A_m2":
                projected_edge_value(current_ref, hotspot_edge) * 1.0e4,
            "sentaurus_integrated_btbt_generation_A_per_um":
                integrate_semiconductor_generation_A_per_um(mesh, btbt_ref),
            "sentaurus_electron_qf_average_V": sent_qf_average,
            "vela_electron_qf_average_V": vela_qf_average,
            "electron_qf_average_offset_vela_minus_sentaurus_V":
                vela_qf_average - sent_qf_average,
            "sentaurus_electrostatic_average_V": sent_psi_average,
            "vela_electrostatic_average_V": vela_psi_average,
            "electrostatic_average_offset_vela_minus_sentaurus_V":
                vela_psi_average - sent_psi_average,
            "electron_psi_minus_qf_offset_error_V":
                (vela_psi_average - vela_qf_average)
                - (sent_psi_average - sent_qf_average),
            "electron_qf_drop_ratio_vela_over_sentaurus": ratio(
                vela_qf_drop, sent_qf_drop
            ),
            "electron_density_ratio_vela_over_sentaurus": ratio(
                vela_n_mid, sent_n_mid
            ),
        })
    write_csv(args.out_dir / "hotspot_bias_evolution.csv", hotspot_evolution)

    with (args.out_dir / "conclusion.md").open("w", encoding="utf-8") as handle:
        handle.write("# BVmethods NMOS 6.4 V SG 逐边审计\n\n")
        handle.write(
            f"Sentaurus 最大电子电流边为 {hotspot['edge_id']} "
            f"({hotspot['node0']}--{hotspot['node1']})。\n\n"
        )
        handle.write("该边 Vela/Sentaurus 分量比：\n\n")
        for key, label in (
            ("electron_current_ratio_vela_over_sentaurus", "电子电流密度"),
            ("electron_qf_drop_ratio", "电子 QF 边差"),
            ("electron_density_ratio", "电子密度对数均值"),
            ("electron_mobility_ratio", "电子迁移率"),
            ("electron_generalized_einstein_ratio", "广义 Einstein 因子"),
            ("geometry_ratio_vela_over_sentaurus", "几何权重"),
        ):
            handle.write(f"- {label}: `{hotspot[key]:.9e}`\n")
        handle.write(
            "\n若 QF 边差、迁移率、Einstein 因子和几何量接近 1，而密度比与电流比"
            "同时极小，则第一处分叉是 QF 绝对位置导致的载流子人口不足，不是 SG "
            "Bernoulli、迁移率或几何系数。\n"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
