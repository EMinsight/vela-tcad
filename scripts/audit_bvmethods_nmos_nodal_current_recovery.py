#!/usr/bin/env python3
"""Audit geometry-based nodal current recovery against Sentaurus BVmethods.

The comparison uses the signed SG edge particle flux exported by Vela and
reconstructs a conventional electron-current vector at every silicon vertex.
No empirical current or avalanche-parameter scale is applied.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path
from typing import Any


Q_C = 1.602176634e-19
REPO = Path(__file__).resolve().parents[1]
RUN = REPO / "build-release/reference_tcad/bvmethods_sentaurus2018/run01"
DEFAULT_MESH = RUN / "vela/mesh.json"
DEFAULT_SENT = RUN / "sentaurus_iic_multibias_exact_extended_20260803/imported"
P1_BRANCH = RUN / (
    "vela_validation/"
    "btbt_e2_iic_qf_vector_nodal_vertex_star_p1_colocated_branch_6p4_7p1_20260805"
)
DEFAULT_VELA = {
    6.4: P1_BRANCH / "fixed_6p4/sg_avalanche_edges.csv",
    7.0: P1_BRANCH / "fixed_7p0/sg_avalanche_edges.csv",
}
DEFAULT_OUT = RUN / "vela_validation/nodal_current_recovery_audit_20260806"


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_rows(path: Path, records: list[dict[str, Any]]) -> None:
    if not records:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def solve_vector(
    projections: list[tuple[float, float, float, float]],
) -> tuple[float, float]:
    """Weighted least-squares vector from (tx, ty, flux, weight)."""
    a00 = a01 = a11 = b0 = b1 = 0.0
    fallback_x = fallback_y = fallback_weight = 0.0
    for tx, ty, flux, weight in projections:
        a00 += weight * tx * tx
        a01 += weight * tx * ty
        a11 += weight * ty * ty
        b0 += weight * tx * flux
        b1 += weight * ty * flux
        fallback_x += weight * flux * tx
        fallback_y += weight * flux * ty
        fallback_weight += weight
    determinant = a00 * a11 - a01 * a01
    scale = max(abs(a00 * a11), abs(a01 * a01), 1.0e-300)
    if len(projections) >= 2 and abs(determinant) > 1.0e-24 * scale:
        return (
            (b0 * a11 - b1 * a01) / determinant,
            (a00 * b1 - a01 * b0) / determinant,
        )
    if fallback_weight > 0.0:
        return fallback_x / fallback_weight, fallback_y / fallback_weight
    return 0.0, 0.0


def percentile(values: list[float], fraction: float) -> float:
    values = sorted(value for value in values if math.isfinite(value))
    if not values:
        return math.nan
    position = fraction * (len(values) - 1)
    lo = int(math.floor(position))
    hi = int(math.ceil(position))
    weight = position - lo
    return values[lo] * (1.0 - weight) + values[hi] * weight


def sentaurus_state(root: Path, bias: float) -> Path:
    return root / f"iic_v{int(bias)}p{round((bias - int(bias)) * 1.0e6):06d}"


def scalar_field(state: Path, name: str) -> dict[int, float]:
    return {
        int(row["node_id"]): float(row["component0"])
        for row in rows(state / "fields" / f"{name}_region3.csv")
    }


def vector_field(state: Path, name: str) -> dict[int, tuple[float, float]]:
    return {
        int(row["node_id"]): (
            float(row["component0"]), float(row["component1"])
        )
        for row in rows(state / "fields" / f"{name}_region3.csv")
    }


def p1_measures(state: Path) -> dict[int, float]:
    points = {
        int(row["id"]): (float(row["x_um"]), float(row["y_um"]))
        for row in rows(state / "nodes.csv")
    }
    measures: dict[int, float] = defaultdict(float)
    for cell in rows(state / "elements.csv"):
        if cell["material"].lower() not in {"si", "silicon"}:
            continue
        ids = [int(cell[f"node{i}"]) for i in range(3)]
        p0, p1, p2 = (points[node] for node in ids)
        area = 0.5e-12 * abs(
            (p1[0] - p0[0]) * (p2[1] - p0[1])
            - (p1[1] - p0[1]) * (p2[0] - p0[0])
        )
        for node in ids:
            measures[node] += area / 3.0
    return dict(measures)


def recovered_electric_field(
    mesh_path: Path, state_path: Path
) -> dict[int, tuple[float, float]]:
    """Area-weighted P1 vertex-star E field used by Vela nodal mode."""
    mesh = json.loads(mesh_path.read_text(encoding="utf-8-sig"))
    points = {
        int(node["id"]): (float(node["x"]), float(node["y"]))
        for node in mesh["nodes"]
    }
    psi = {
        int(row["node_id"]): float(row["psi"])
        for row in rows(state_path)
    }
    weighted: dict[int, list[float]] = defaultdict(lambda: [0.0, 0.0, 0.0])
    for cell in mesh["triangles"]:
        if int(cell["region_id"]) != 3:
            continue
        n0, n1, n2 = (int(value) for value in cell["node_ids"])
        x0, y0 = points[n0]
        x1, y1 = points[n1]
        x2, y2 = points[n2]
        double_area = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)
        if abs(double_area) <= 1.0e-30:
            continue
        # Linear-triangle gradient in V/um; convert to E in V/m.
        gx = (
            psi[n0] * (y1 - y2)
            + psi[n1] * (y2 - y0)
            + psi[n2] * (y0 - y1)
        ) / double_area
        gy = (
            psi[n0] * (x2 - x1)
            + psi[n1] * (x0 - x2)
            + psi[n2] * (x1 - x0)
        ) / double_area
        area = 0.5 * abs(double_area)
        for node in (n0, n1, n2):
            weighted[node][0] += area * gx
            weighted[node][1] += area * gy
            weighted[node][2] += area
    return {
        node: (-1.0e6 * values[0] / values[2], -1.0e6 * values[1] / values[2])
        for node, values in weighted.items() if values[2] > 0.0
    }


def eparallel(
    electric: tuple[float, float], current: tuple[float, float]
) -> float:
    magnitude = math.hypot(*current)
    if magnitude <= 1.0e-300:
        return 0.0
    return max((electric[0] * current[0] + electric[1] * current[1]) / magnitude, 0.0)


def electron_alpha(field_v_m: float) -> float:
    if field_v_m <= 0.0:
        return 0.0
    return 7.03e7 * math.exp(-1.231e8 / field_v_m)


def nodal_edge_ls_with_length(
    edge_rows: list[dict[str, str]],
    bias: float,
    mode: str,
) -> dict[int, tuple[float, float]]:
    projections: dict[int, list[tuple[float, float, float, float]]] = defaultdict(list)
    for row in edge_rows:
        if not math.isclose(float(row["bias_V"]), bias, abs_tol=1.0e-10):
            continue
        length = float(row["edge_length_m"])
        couple = float(row["edge_couple_m"])
        mobility = float(row["electron_mobility_m2_V_s"])
        if length <= 1.0e-30 or couple <= 0.0 or mobility <= 0.0:
            continue
        dx = (float(row["x1_um"]) - float(row["x0_um"])) * 1.0e-6
        dy = (float(row["y1_um"]) - float(row["y0_um"])) * 1.0e-6
        flux = float(row["electron_raw_signed_flux_proxy"])
        if mode == "uniform":
            weight = 1.0
        elif mode == "edge_length":
            weight = length
        elif mode == "dual_face":
            weight = couple
        elif mode == "box_area":
            weight = length * couple
        elif mode == "inverse_length":
            weight = 1.0 / length
        else:
            raise ValueError(mode)
        item = (dx / length, dy / length, flux, weight)
        projections[int(row["node0"])].append(item)
        projections[int(row["node1"])].append(item)
    output = {}
    for node, values in projections.items():
        px, py = solve_vector(values)
        output[node] = (-Q_C * px * 1.0e4, -Q_C * py * 1.0e4)
    return output


def audit_bias(
    bias: float, vela_path: Path, sent_root: Path, mesh_path: Path
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    state = sentaurus_state(sent_root, bias)
    sent_current_cm2 = vector_field(state, "eCurrentDensity")
    sent_current = {
        node: (value[0] * 1.0e4, value[1] * 1.0e4)
        for node, value in sent_current_cm2.items()
    }
    sent_generation = {
        node: abs(value) * 1.0e6
        for node, value in scalar_field(state, "eImpactIonization").items()
    }
    sent_alpha = {
        node: abs(value) * 100.0
        for node, value in scalar_field(state, "eAlphaAvalanche").items()
    }
    measure = p1_measures(state)
    sent_electric_cm = vector_field(state, "ElectricField")
    sent_electric = {
        node: (value[0] * 100.0, value[1] * 100.0)
        for node, value in sent_electric_cm.items()
    }
    edge_rows = rows(vela_path)
    candidates = {
        mode: nodal_edge_ls_with_length(edge_rows, bias, mode)
        for mode in ("uniform", "edge_length", "dual_face", "box_area", "inverse_length")
    }
    vela_electric = recovered_electric_field(mesh_path, vela_path.parent / "last_state.csv")
    peak_generation = max(sent_generation.values())
    details: list[dict[str, Any]] = []
    summaries: list[dict[str, Any]] = []
    sent_source = sum(
        sent_generation.get(node, 0.0) * volume for node, volume in measure.items()
    )
    for mode, recovered in candidates.items():
        for threshold in (0.1, 0.3, 0.5, 0.8):
            selected = [
                node for node, value in sent_generation.items()
                if value >= threshold * peak_generation
                and node in sent_current and node in recovered
            ]
            ratios = []
            cosine = []
            drive_ratios = []
            alpha_ratios = []
            sent_formula_alpha_ratios = []
            for node in selected:
                sx, sy = sent_current[node]
                vx, vy = recovered[node]
                smag, vmag = math.hypot(sx, sy), math.hypot(vx, vy)
                if smag > 0.0:
                    ratios.append(vmag / smag)
                if smag > 0.0 and vmag > 0.0:
                    cosine.append((sx * vx + sy * vy) / (smag * vmag))
                sent_drive = eparallel(
                    sent_electric.get(node, (0.0, 0.0)), sent_current[node]
                )
                vela_drive = eparallel(
                    vela_electric.get(node, (0.0, 0.0)), recovered[node]
                )
                if sent_drive > 0.0:
                    drive_ratios.append(vela_drive / sent_drive)
                exported_alpha = sent_alpha.get(node, 0.0)
                if exported_alpha > 0.0:
                    alpha_ratios.append(electron_alpha(vela_drive) / exported_alpha)
                    sent_formula_alpha_ratios.append(
                        electron_alpha(sent_drive) / exported_alpha
                    )
            reconstructed_source = sum(
                sent_alpha.get(node, 0.0)
                * math.hypot(*recovered.get(node, (0.0, 0.0))) / Q_C
                * volume
                for node, volume in measure.items()
            )
            vela_alpha_source = sum(
                electron_alpha(eparallel(vela_electric[node], recovered[node]))
                * math.hypot(*recovered[node]) / Q_C * volume
                for node, volume in measure.items()
                if node in vela_electric and node in recovered
            )
            summaries.append({
                "bias_V": bias,
                "mode": mode,
                "threshold_fraction_of_peak_generation": threshold,
                "selected_nodes": len(selected),
                "magnitude_ratio_p10": percentile(ratios, 0.10),
                "magnitude_ratio_p50": percentile(ratios, 0.50),
                "magnitude_ratio_p90": percentile(ratios, 0.90),
                "direction_cosine_p50": percentile(cosine, 0.50),
                "eparallel_ratio_p10": percentile(drive_ratios, 0.10),
                "eparallel_ratio_p50": percentile(drive_ratios, 0.50),
                "eparallel_ratio_p90": percentile(drive_ratios, 0.90),
                "vela_formula_alpha_over_sentaurus_export_p50": percentile(
                    alpha_ratios, 0.50
                ),
                "sentaurus_formula_alpha_over_export_p50": percentile(
                    sent_formula_alpha_ratios, 0.50
                ),
                "sentaurus_electron_source_per_m_s": sent_source,
                "sentaurus_alpha_recovered_current_source_per_m_s": reconstructed_source,
                "reconstructed_over_sentaurus_source": reconstructed_source / sent_source,
                "vela_field_alpha_recovered_current_source_per_m_s": vela_alpha_source,
                "vela_field_alpha_source_over_sentaurus": vela_alpha_source / sent_source,
            })
        for node, (vx, vy) in recovered.items():
            if node not in sent_current or node not in sent_generation:
                continue
            sx, sy = sent_current[node]
            sent_drive = eparallel(
                sent_electric.get(node, (0.0, 0.0)), sent_current[node]
            )
            vela_drive = eparallel(
                vela_electric.get(node, (0.0, 0.0)), recovered[node]
            )
            details.append({
                "bias_V": bias,
                "mode": mode,
                "node_id": node,
                "sentaurus_generation_m3_s": sent_generation[node],
                "sentaurus_jx_A_m2": sx,
                "sentaurus_jy_A_m2": sy,
                "vela_jx_A_m2": vx,
                "vela_jy_A_m2": vy,
                "sentaurus_eparallel_from_export_V_m": sent_drive,
                "vela_eparallel_V_m": vela_drive,
                "vela_over_sentaurus_eparallel": (
                    vela_drive / sent_drive if sent_drive > 0.0 else math.nan
                ),
                "sentaurus_exported_alpha_m_inv": sent_alpha.get(node, 0.0),
                "sentaurus_formula_alpha_m_inv": electron_alpha(sent_drive),
                "vela_formula_alpha_m_inv": electron_alpha(vela_drive),
                "p1_measure_m2": measure.get(node, 0.0),
                "sentaurus_source_per_m_s": (
                    sent_generation[node] * measure.get(node, 0.0)
                ),
                "vela_formula_source_per_m_s": (
                    electron_alpha(vela_drive) * math.hypot(vx, vy) / Q_C
                    * measure.get(node, 0.0)
                ),
                "vela_minus_sentaurus_source_per_m_s": (
                    electron_alpha(vela_drive) * math.hypot(vx, vy) / Q_C
                    * measure.get(node, 0.0)
                    - sent_generation[node] * measure.get(node, 0.0)
                ),
                "vela_over_sentaurus_magnitude": (
                    math.hypot(vx, vy) / math.hypot(sx, sy)
                    if math.hypot(sx, sy) > 0.0 else math.nan
                ),
            })
    return summaries, details


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sentaurus-states", type=Path, default=DEFAULT_SENT)
    parser.add_argument("--mesh", type=Path, default=DEFAULT_MESH)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--bias", type=float, action="append")
    args = parser.parse_args()
    biases = args.bias or sorted(DEFAULT_VELA)
    summaries: list[dict[str, Any]] = []
    details: list[dict[str, Any]] = []
    for bias in biases:
        local_summary, local_details = audit_bias(
            bias, DEFAULT_VELA[bias], args.sentaurus_states, args.mesh
        )
        summaries.extend(local_summary)
        details.extend(local_details)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_rows(args.out_dir / "summary.csv", summaries)
    write_rows(args.out_dir / "node_details.csv", details)
    compact = {
        f"{row['bias_V']:.1f}V/{row['mode']}": row
        for row in summaries
        if row["threshold_fraction_of_peak_generation"] == 0.3
    }
    (args.out_dir / "summary.json").write_text(
        json.dumps(compact, indent=2) + "\n", encoding="utf-8"
    )
    print(args.out_dir / "summary.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
