#!/usr/bin/env python3
"""Compare exact Sentaurus BVmethods checkpoints with Vela node fields."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np


REPO = Path(__file__).resolve().parents[1]
RUN_ROOT = REPO / "build-release/reference_tcad/bvmethods_sentaurus2018/run01"
DEFAULT_SENT = RUN_ROOT / "sentaurus_iic_multibias_exact_extended_20260803/imported"
DEFAULT_VELA = RUN_ROOT / "vela_validation/iic_postprocess_20260803/probes"
DEFAULT_OUT = RUN_ROOT / "vela_validation/iic_postprocess_20260803/analysis/multibias_sentaurus"
DEFAULT_BIASES = [1.0, 2.0, 4.0, 5.0, 6.0, 6.32, 6.34, 6.36, 6.37, 6.38, 6.39, 6.4]


EDGE_FIELD_SPECS = {
    "electric_field": ("ElectricField", 100.0),
    "electron_current_density": ("eCurrentDensity", 1.0e4),
    "hole_current_density": ("hCurrentDensity", 1.0e4),
    "electron_alpha": ("eAlphaAvalanche", 100.0),
    "hole_alpha": ("hAlphaAvalanche", 100.0),
    "avalanche_generation": ("ImpactIonization", 1.0e6),
}

Q_C = 1.602176634e-19
# unit_scaling stores alpha in cm^-1, particle flux in cm^-2 s^-1, and the
# 2-D source support in um^2.  Their product therefore needs
#   1e2 * 1e4 * 1e-12 = 1e-6
# to become an SI line source in m^-1 s^-1.
VELA_SOURCE_INTEGRAL_TO_PER_M_S = 1.0e-6


def read_rows(path: Path) -> list[dict[str, str]]:
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


def sent_tag(bias: float) -> str:
    return f"iic_v{bias:.6f}".replace(".", "p")


def vela_tag(bias: float) -> str:
    raw = f"{bias:.6f}".rstrip("0").rstrip(".").replace(".", "p")
    return f"v{raw}"


def discover_vela_cases(root: Path) -> dict[float, Path]:
    cases: dict[float, Path] = {}
    for path in root.iterdir():
        summary = path / "postprocess_only/avalanche_summary.csv"
        if not summary.exists():
            continue
        case_rows = read_rows(summary)
        if case_rows:
            cases[round(float(case_rows[0]["bias_V"]), 12)] = path
    return cases


def parse_vtk(path: Path) -> dict[str, Any]:
    lines = path.read_text(encoding="utf-8").splitlines()
    points: list[tuple[float, float]] = []
    fields: dict[str, np.ndarray] = {}
    index = 0
    while index < len(lines):
        parts = lines[index].split()
        if not parts:
            index += 1
            continue
        if parts[0] == "POINTS":
            count = int(parts[1])
            for raw in lines[index + 1:index + 1 + count]:
                x, y, *_ = raw.split()
                # This imported Sentaurus-derived Vela mesh is stored in um in
                # both the legacy VTK coordinates and nodes.csv.
                points.append((float(x), float(y)))
            index += 1 + count
            continue
        if len(parts) >= 3 and parts[0] == "SCALARS":
            name = parts[1]
            index += 2
            values: list[float] = []
            while index < len(lines):
                tokens = lines[index].split()
                if not tokens or tokens[0] in {"SCALARS", "VECTORS", "FIELD", "CELL_DATA", "POINT_DATA"}:
                    break
                values.extend(float(value) for value in tokens)
                index += 1
            fields[name] = np.asarray(values, dtype=float)
            continue
        if len(parts) >= 3 and parts[0] == "VECTORS":
            name = parts[1]
            index += 1
            values: list[float] = []
            while index < len(lines):
                tokens = lines[index].split()
                if not tokens or tokens[0] in {"SCALARS", "VECTORS", "FIELD", "CELL_DATA", "POINT_DATA"}:
                    break
                vector = [float(value) for value in tokens[:3]]
                values.append(math.sqrt(sum(value * value for value in vector)))
                index += 1
            fields[name] = np.asarray(values, dtype=float)
            continue
        index += 1
    return {"points": np.asarray(points, dtype=float), "fields": fields}


def sent_field(root: Path, name: str) -> dict[int, float]:
    path = root / "fields" / f"{name}_region3.csv"
    values: dict[int, float] = {}
    for row in read_rows(path):
        components = [float(value) for key, value in row.items() if key != "node_id" and value not in (None, "")]
        values[int(row["node_id"])] = (
            abs(components[0]) if len(components) == 1
            else math.sqrt(sum(value * value for value in components))
        )
    return values


def sent_vector_field(root: Path, name: str) -> dict[int, tuple[float, float]]:
    values: dict[int, tuple[float, float]] = {}
    for row in read_rows(root / "fields" / f"{name}_region3.csv"):
        values[int(row["node_id"])] = (float(row["component0"]), float(row["component1"]))
    return values


def sent_edge_reference(
    row: dict[str, str],
    quantity: str,
    scalar_values: dict[int, float] | None,
    vector_values: dict[int, tuple[float, float]] | None,
    scale: float,
) -> float:
    node0 = int(row["node0"])
    node1 = int(row["node1"])
    if vector_values is None:
        assert scalar_values is not None
        return 0.5 * (scalar_values[node0] + scalar_values[node1]) * scale
    dx = float(row["x1_um"]) - float(row["x0_um"])
    dy = float(row["y1_um"]) - float(row["y0_um"])
    length = math.hypot(dx, dy)
    if length <= 0.0:
        return 0.0
    vx = 0.5 * (vector_values[node0][0] + vector_values[node1][0])
    vy = 0.5 * (vector_values[node0][1] + vector_values[node1][1])
    return abs((vx * dx + vy * dy) / length) * scale


def percentile(values: np.ndarray, q: float) -> float:
    return float(np.nanpercentile(values, q)) if len(values) else math.nan


def field_metrics(reference: np.ndarray, candidate: np.ndarray) -> dict[str, float]:
    ref_abs = np.abs(reference)
    cand_abs = np.abs(candidate)
    ref_p99 = percentile(ref_abs, 99)
    threshold = max(ref_p99 * 1.0e-9, 1.0e-300)
    active = np.isfinite(reference) & np.isfinite(candidate) & (ref_abs >= threshold)
    if not np.any(active):
        return {key: math.nan for key in (
            "active_nodes", "relative_error_p50", "relative_error_p95",
            "abs_log10_ratio_p50", "abs_log10_ratio_p95", "normalized_rmse",
            "correlation",
        )}
    rel = np.abs(candidate[active] - reference[active]) / np.maximum(ref_abs[active], threshold)
    logerr = np.abs(
        np.log10(np.maximum(cand_abs[active], 1.0e-300))
        - np.log10(np.maximum(ref_abs[active], 1.0e-300))
    )
    scale = max(ref_p99, 1.0e-300)
    corr = math.nan
    if np.std(reference[active]) > 0.0 and np.std(candidate[active]) > 0.0:
        corr = float(np.corrcoef(reference[active], candidate[active])[0, 1])
    return {
        "active_nodes": int(np.count_nonzero(active)),
        "relative_error_p50": percentile(rel, 50),
        "relative_error_p95": percentile(rel, 95),
        "abs_log10_ratio_p50": percentile(logerr, 50),
        "abs_log10_ratio_p95": percentile(logerr, 95),
        "normalized_rmse": float(np.sqrt(np.mean((candidate - reference) ** 2)) / scale),
        "correlation": corr,
    }


def contact_scalar(root: Path, field: str, region: int) -> float:
    return float(read_rows(root / "fields" / f"{field}_region{region}.csv")[0]["component0"])


def vela_edge_value(row: dict[str, str], quantity: str) -> float:
    if quantity == "electric_field":
        return abs(float(row["electric_field_V_per_m"]))
    if quantity == "electron_alpha":
        return abs(float(row["electron_alpha_m_inv"]))
    if quantity == "hole_alpha":
        return abs(float(row["hole_alpha_m_inv"]))
    if quantity == "electron_current_density":
        return abs(float(row["electron_sg_production_signed_conventional_current_density_A_per_m2"]))
    if quantity == "hole_current_density":
        return Q_C * 1.0e4 * abs(float(row["hole_flux_proxy"]))
    if quantity == "avalanche_generation":
        area = float(row["edge_area_proxy_m2"])
        return (
            float(row["edge_source_integral"])
            * VELA_SOURCE_INTEGRAL_TO_PER_M_S
            / area
            if area > 0.0 else 0.0
        )
    raise KeyError(quantity)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sentaurus-root", type=Path, default=DEFAULT_SENT)
    parser.add_argument("--vela-root", type=Path, default=DEFAULT_VELA)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--biases", default=",".join(str(value) for value in DEFAULT_BIASES))
    parser.add_argument("--top-n", type=int, default=20)
    args = parser.parse_args()
    biases = [float(value) for value in args.biases.split(",") if value.strip()]

    summary: list[dict[str, Any]] = []
    matched_edges: list[dict[str, Any]] = []
    node_potential: list[dict[str, Any]] = []
    potential_summary: list[dict[str, Any]] = []
    node_carrier_state: list[dict[str, Any]] = []
    carrier_state_summary: list[dict[str, Any]] = []
    top_errors: list[dict[str, Any]] = []
    checkpoint: list[dict[str, Any]] = []
    coordinate_errors: list[float] = []
    vela_cases = discover_vela_cases(args.vela_root)

    for bias in biases:
        sent = args.sentaurus_root / sent_tag(bias)
        vela_case = vela_cases.get(round(bias, 12), args.vela_root / vela_tag(bias))
        vtk_paths = sorted((vela_case / "postprocess_only/vtk").glob("*.vtk"))
        if not sent.exists() or len(vtk_paths) != 1:
            raise FileNotFoundError(f"missing unique inputs at {bias:g} V")
        vtk = parse_vtk(vtk_paths[0])
        points = vtk["points"]
        node_rows = read_rows(sent / "nodes.csv")
        coordinates = {int(row["id"]): (float(row["x_um"]), float(row["y_um"])) for row in node_rows}
        semiconductor_nodes = sorted(sent_field(sent, "ElectrostaticPotential"))
        max_coordinate_error = max(
            math.hypot(points[node][0] - coordinates[node][0], points[node][1] - coordinates[node][1])
            for node in semiconductor_nodes
        )
        coordinate_errors.append(max_coordinate_error)

        checkpoint.append({
            "bias_V": bias,
            "sentaurus_drain_current_A_per_um": contact_scalar(sent, "ContactCurrentFlux", 8),
            "sentaurus_source_current_A_per_um": contact_scalar(sent, "ContactCurrentFlux", 9),
            "sentaurus_substrate_current_A_per_um": contact_scalar(sent, "ContactCurrentFlux", 7),
            "sentaurus_max_electron_ion_integral": max(sent_field(sent, "eIonIntegral").values()),
            "sentaurus_max_hole_ion_integral": max(sent_field(sent, "hIonIntegral").values()),
            "sentaurus_max_mean_ion_integral": max(sent_field(sent, "MeanIonIntegral").values()),
            "max_coordinate_error_um": max_coordinate_error,
        })

        potential = sent_field(sent, "ElectrostaticPotential")
        vela_state_rows = read_rows(vela_case / "postprocess_only/last_state.csv")
        vela_state = {int(row["node_id"]): row for row in vela_state_rows}
        vela_potential = np.asarray([
            float(vela_state[node]["psi"]) for node in range(len(points))
        ])
        potential_reference = np.asarray([potential[node] for node in semiconductor_nodes])
        potential_candidate = np.asarray([vela_potential[node] for node in semiconductor_nodes])
        potential_error = np.abs(potential_candidate - potential_reference)
        potential_summary.append({
            "bias_V": bias,
            "nodes": len(semiconductor_nodes),
            "absolute_error_p50_V": percentile(potential_error, 50),
            "absolute_error_p95_V": percentile(potential_error, 95),
            "absolute_error_max_V": float(np.max(potential_error)),
            "rms_error_V": float(np.sqrt(np.mean(potential_error ** 2))),
            "correlation": float(np.corrcoef(potential_reference, potential_candidate)[0, 1]),
        })
        for node in semiconductor_nodes:
            node_potential.append({
                "bias_V": bias,
                "node_id": node,
                "x_um": coordinates[node][0],
                "y_um": coordinates[node][1],
                "sentaurus_potential_V": potential[node],
                "vela_potential_V": float(vela_potential[node]),
                "abs_error_V": abs(float(vela_potential[node]) - potential[node]),
            })

        state_specs = (
            ("electron_quasi_fermi", "eQuasiFermiPotential", "phin", "V", 1.0),
            ("hole_quasi_fermi", "hQuasiFermiPotential", "phip", "V", 1.0),
            # Sentaurus TDR carrier densities are exported in cm^-3, whereas
            # DDSolutionCsv explicitly writes physical m^-3.
            ("electron_density", "eDensity", "electrons_m3", "dex", 1.0e6),
            ("hole_density", "hDensity", "holes_m3", "dex", 1.0e6),
        )
        for quantity, sent_name, vela_name, error_unit, sent_scale in state_specs:
            sent_state = sent_field(sent, sent_name)
            reference = np.asarray([
                sent_state[node] * sent_scale for node in semiconductor_nodes
            ])
            candidate = np.asarray([
                float(vela_state[node][vela_name]) for node in semiconductor_nodes
            ])
            if error_unit == "V":
                errors = np.abs(candidate - reference)
            else:
                errors = np.abs(
                    np.log10(np.maximum(np.abs(candidate), 1.0e-300))
                    - np.log10(np.maximum(np.abs(reference), 1.0e-300))
                )
            corr = math.nan
            if np.std(reference) > 0.0 and np.std(candidate) > 0.0:
                corr = float(np.corrcoef(reference, candidate)[0, 1])
            carrier_state_summary.append({
                "bias_V": bias,
                "quantity": quantity,
                "nodes": len(semiconductor_nodes),
                "error_unit": error_unit,
                "absolute_error_p50": percentile(errors, 50),
                "absolute_error_p95": percentile(errors, 95),
                "absolute_error_max": float(np.max(errors)),
                "correlation": corr,
                "vela_nonpositive_count": int(np.count_nonzero(candidate <= 0.0))
                    if error_unit == "dex" else 0,
            })
            for index, node in enumerate(semiconductor_nodes):
                node_carrier_state.append({
                    "bias_V": bias,
                    "quantity": quantity,
                    "node_id": node,
                    "x_um": coordinates[node][0],
                    "y_um": coordinates[node][1],
                    "sentaurus_value": reference[index],
                    "vela_value": candidate[index],
                    "error_unit": error_unit,
                    "absolute_error": errors[index],
                })

        edge_rows = read_rows(vela_case / "postprocess_only/sg_avalanche_edges.csv")
        for quantity, (sent_name, sent_scale) in EDGE_FIELD_SPECS.items():
            is_vector = quantity in {
                "electric_field", "electron_current_density", "hole_current_density"
            }
            sent_values = None if is_vector else sent_field(sent, sent_name)
            sent_vectors = sent_vector_field(sent, sent_name) if is_vector else None
            supported_nodes = set(sent_vectors if sent_vectors is not None else sent_values or {})
            selected = [
                row for row in edge_rows
                if int(row["node0"]) in supported_nodes and int(row["node1"]) in supported_nodes
                and float(row["edge_area_proxy_m2"]) > 0.0
            ]
            reference = np.asarray([
                sent_edge_reference(
                    row, quantity, sent_values, sent_vectors, sent_scale
                )
                for row in selected
            ])
            candidate = np.asarray([vela_edge_value(row, quantity) for row in selected])
            metrics = field_metrics(reference, candidate)
            s_peak_index = int(np.nanargmax(np.abs(reference)))
            v_peak_index = int(np.nanargmax(np.abs(candidate)))
            row = {
                "bias_V": bias,
                "quantity": quantity,
                "matched_edges": len(selected),
                "sentaurus_peak": float(np.abs(reference[s_peak_index])),
                "vela_peak": float(np.abs(candidate[v_peak_index])),
                "abs_vela_over_sentaurus_peak": (
                    float(np.abs(candidate[v_peak_index]) / np.abs(reference[s_peak_index]))
                    if reference[s_peak_index] != 0.0 else math.nan
                ),
                "sentaurus_peak_edge": int(selected[s_peak_index]["edge_id"]),
                "vela_peak_edge": int(selected[v_peak_index]["edge_id"]),
                **metrics,
            }
            summary.append(row)

            ref_abs = np.abs(reference)
            ref_floor = max(percentile(ref_abs, 99) * 1.0e-9, 1.0e-300)
            logerr = np.abs(
                np.log10(np.maximum(np.abs(candidate), 1.0e-300))
                - np.log10(np.maximum(ref_abs, 1.0e-300))
            )
            active = ref_abs >= ref_floor
            ranking = np.argsort(np.where(active, logerr, -1.0))[::-1][: args.top_n]
            for rank, index in enumerate(ranking, start=1):
                if not active[index]:
                    continue
                edge = selected[index]
                top_errors.append({
                    "bias_V": bias,
                    "quantity": quantity,
                    "rank": rank,
                    "edge_id": edge["edge_id"],
                    "node0": edge["node0"],
                    "node1": edge["node1"],
                    "x_mid_um": 0.5 * (float(edge["x0_um"]) + float(edge["x1_um"])),
                    "y_mid_um": 0.5 * (float(edge["y0_um"]) + float(edge["y1_um"])),
                    "sentaurus_value": reference[index],
                    "vela_value": candidate[index],
                    "abs_log10_ratio": logerr[index],
                })
            for index, edge in enumerate(selected):
                matched_edges.append({
                    "bias_V": bias,
                    "quantity": quantity,
                    "edge_id": edge["edge_id"],
                    "node0": edge["node0"],
                    "node1": edge["node1"],
                    "x_mid_um": 0.5 * (float(edge["x0_um"]) + float(edge["x1_um"])),
                    "y_mid_um": 0.5 * (float(edge["y0_um"]) + float(edge["y1_um"])),
                    "sentaurus_value": reference[index],
                    "vela_value": candidate[index],
                    "abs_log10_ratio": logerr[index] if active[index] else "",
                    "sentaurus_active": int(active[index]),
                })

    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_rows(args.out_dir / "checkpoint_scalars.csv", checkpoint)
    write_rows(args.out_dir / "matched_edge_field_summary.csv", summary)
    write_rows(args.out_dir / "top_matched_edge_errors.csv", top_errors)
    write_rows(args.out_dir / "matched_edge_fields.csv", matched_edges)
    write_rows(args.out_dir / "same_node_potential.csv", node_potential)
    write_rows(args.out_dir / "same_node_potential_summary.csv", potential_summary)
    write_rows(args.out_dir / "same_node_carrier_state.csv", node_carrier_state)
    write_rows(args.out_dir / "same_node_carrier_state_summary.csv", carrier_state_summary)
    result = {
        "biases_V": biases,
        "edge_field_quantities": list(EDGE_FIELD_SPECS),
        "same_node_mapping": True,
        "maximum_coordinate_error_um": max(coordinate_errors),
        "node_field_status": (
            "potential is valid; Vela derived node E/alpha/current fields are excluded because "
            "DCSweep did not pass unit_scaling into writeDDSolutionVTK"
        ),
        "sentaurus_root": str(args.sentaurus_root.resolve()),
        "vela_root": str(args.vela_root.resolve()),
        "summary_csv": str((args.out_dir / "matched_edge_field_summary.csv").resolve()),
        "matched_edge_csv": str((args.out_dir / "matched_edge_fields.csv").resolve()),
        "same_node_potential_csv": str((args.out_dir / "same_node_potential.csv").resolve()),
        "same_node_potential_summary_csv": str(
            (args.out_dir / "same_node_potential_summary.csv").resolve()
        ),
        "same_node_carrier_state_summary_csv": str(
            (args.out_dir / "same_node_carrier_state_summary.csv").resolve()
        ),
    }
    (args.out_dir / "result.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(args.out_dir / "result.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
