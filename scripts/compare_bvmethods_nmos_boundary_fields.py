#!/usr/bin/env python3
"""Compare Vela and imported Sentaurus fields at the NMOS BV current boundary."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from statistics import median
from typing import Any, Iterable


REPO = Path(__file__).resolve().parents[1]
RUN_ROOT = REPO / "build-release/reference_tcad/bvmethods_sentaurus2018/run01"
DEFAULT_SENT = RUN_ROOT / "sentaurus_boundary_state_20260808/imported/current_1e4"
DEFAULT_VELA_STATE = (
    RUN_ROOT
    / "vela_validation/boundary_voltage_to_current_20260806"
    / "boundary_control_checkpoints/current_target_0p000100_eval_8.csv"
)
DEFAULT_MESH = RUN_ROOT / "vela/mesh.json"
DEFAULT_OUTPUT = (
    RUN_ROOT
    / "sentaurus_boundary_state_20260808/analysis/vela_sentaurus_field_comparison.json"
)


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def percentile(values: Iterable[float], q: float) -> float:
    ordered = sorted(value for value in values if math.isfinite(value))
    if not ordered:
        return math.nan
    index = (len(ordered) - 1) * q / 100.0
    lower = math.floor(index)
    upper = math.ceil(index)
    if lower == upper:
        return ordered[lower]
    return ordered[lower] + (index - lower) * (ordered[upper] - ordered[lower])


def correlation(reference: list[float], candidate: list[float]) -> float:
    if len(reference) != len(candidate) or len(reference) < 2:
        return math.nan
    ref_mean = sum(reference) / len(reference)
    cand_mean = sum(candidate) / len(candidate)
    ref_var = sum((value - ref_mean) ** 2 for value in reference)
    cand_var = sum((value - cand_mean) ** 2 for value in candidate)
    if ref_var <= 0.0 or cand_var <= 0.0:
        return math.nan
    covariance = sum(
        (ref - ref_mean) * (cand - cand_mean)
        for ref, cand in zip(reference, candidate)
    )
    return covariance / math.sqrt(ref_var * cand_var)


def absolute_metrics(reference: list[float], candidate: list[float]) -> dict[str, float]:
    errors = [abs(cand - ref) for ref, cand in zip(reference, candidate)]
    return {
        "count": len(errors),
        "absolute_error_p50": percentile(errors, 50),
        "absolute_error_p95": percentile(errors, 95),
        "absolute_error_max": max(errors, default=math.nan),
        "rms_error": (
            math.sqrt(sum(error * error for error in errors) / len(errors))
            if errors else math.nan
        ),
        "correlation": correlation(reference, candidate),
    }


def relative_metrics(reference: list[float], candidate: list[float]) -> dict[str, float]:
    relative = [
        abs(cand - ref) / abs(ref)
        for ref, cand in zip(reference, candidate)
        if ref != 0.0
    ]
    return {
        "relative_error_p50": percentile(relative, 50),
        "relative_error_p95": percentile(relative, 95),
        "relative_error_max": max(relative, default=math.nan),
    }


def sent_scalar(root: Path, name: str) -> dict[int, float]:
    return {
        int(row["node_id"]): float(row["component0"])
        for row in read_rows(root / "fields" / f"{name}_region3.csv")
    }


def sent_vector(root: Path, name: str) -> dict[int, tuple[float, float]]:
    return {
        int(row["node_id"]): (float(row["component0"]), float(row["component1"]))
        for row in read_rows(root / "fields" / f"{name}_region3.csv")
    }


def contact_scalar(root: Path, name: str, region: int = 8) -> float:
    rows = read_rows(root / "fields" / f"{name}_region{region}.csv")
    return float(rows[0]["component0"])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sentaurus-root", type=Path, default=DEFAULT_SENT)
    parser.add_argument("--vela-state", type=Path, default=DEFAULT_VELA_STATE)
    parser.add_argument("--mesh", type=Path, default=DEFAULT_MESH)
    parser.add_argument("--output-json", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    vela = {int(row["node_id"]): row for row in read_rows(args.vela_state)}
    mesh = json.loads(args.mesh.read_text(encoding="utf-8"))
    coordinates = {
        int(node["id"]): (float(node["x"]), float(node["y"]))
        for node in mesh["nodes"]
    }

    potential = sent_scalar(args.sentaurus_root, "ElectrostaticPotential")
    electron_qf = sent_scalar(args.sentaurus_root, "eQuasiFermiPotential")
    hole_qf = sent_scalar(args.sentaurus_root, "hQuasiFermiPotential")
    # Sentaurus exports density in cm^-3; Vela checkpoints store m^-3.
    electron_density = {
        node: value * 1.0e6
        for node, value in sent_scalar(args.sentaurus_root, "eDensity").items()
    }
    hole_density = {
        node: value * 1.0e6
        for node, value in sent_scalar(args.sentaurus_root, "hDensity").items()
    }
    semiconductor_nodes = sorted(potential)

    sent_potential = [potential[node] for node in semiconductor_nodes]
    vela_potential = [float(vela[node]["psi"]) for node in semiconductor_nodes]
    potential_offset = median(
        candidate - reference
        for reference, candidate in zip(sent_potential, vela_potential)
    )
    aligned_potential = [value - potential_offset for value in vela_potential]

    scalar_results: dict[str, Any] = {
        "potential_V": {
            "raw": absolute_metrics(sent_potential, vela_potential),
            "median_candidate_minus_reference_V": potential_offset,
            "median_offset_aligned": absolute_metrics(sent_potential, aligned_potential),
        }
    }

    density_specs = {
        "electron": (electron_density, electron_qf, "electrons_m3", "phin"),
        "hole": (hole_density, hole_qf, "holes_m3", "phip"),
    }
    for carrier, (sent_density, sent_qf, vela_density_name, vela_qf_name) in density_specs.items():
        nodes = sorted(set(semiconductor_nodes) & set(sent_density) & set(sent_qf))
        reference_density = [sent_density[node] for node in nodes]
        candidate_density = [float(vela[node][vela_density_name]) for node in nodes]
        positive = [
            index for index, (ref, cand) in enumerate(zip(reference_density, candidate_density))
            if ref > 0.0 and cand > 0.0
        ]
        density_log_reference = [math.log10(reference_density[index]) for index in positive]
        density_log_candidate = [math.log10(candidate_density[index]) for index in positive]
        peak_density = max(reference_density)
        populated_threshold = peak_density * 1.0e-6
        populated = [
            index for index in positive if reference_density[index] >= populated_threshold
        ]
        qf_reference = [sent_qf[node] for node in nodes]
        qf_candidate = [float(vela[node][vela_qf_name]) for node in nodes]
        scalar_results[f"{carrier}_density_dex"] = {
            "all_positive": absolute_metrics(density_log_reference, density_log_candidate),
            "sentaurus_peak_m3": peak_density,
            "populated_threshold_m3": populated_threshold,
            "populated": absolute_metrics(
                [math.log10(reference_density[index]) for index in populated],
                [math.log10(candidate_density[index]) for index in populated],
            ),
        }
        scalar_results[f"{carrier}_quasi_fermi_V"] = {
            "all_semiconductor_nodes": absolute_metrics(qf_reference, qf_candidate),
            "populated_carrier_nodes": absolute_metrics(
                [qf_reference[index] for index in populated],
                [qf_candidate[index] for index in populated],
            ),
        }

    sent_electric_field = sent_vector(args.sentaurus_root, "ElectricField")
    seen_edges: set[tuple[int, int]] = set()
    sent_edge_field: list[float] = []
    vela_edge_field: list[float] = []
    for triangle in mesh["triangles"]:
        if int(triangle["region_id"]) != 3:
            continue
        node_ids = [int(value) for value in triangle["node_ids"]]
        for first, second in (
            (node_ids[0], node_ids[1]),
            (node_ids[1], node_ids[2]),
            (node_ids[2], node_ids[0]),
        ):
            edge = (min(first, second), max(first, second))
            if edge in seen_edges or first not in sent_electric_field or second not in sent_electric_field:
                continue
            seen_edges.add(edge)
            x0, y0 = coordinates[first]
            x1, y1 = coordinates[second]
            dx_um = x1 - x0
            dy_um = y1 - y0
            length_um = math.hypot(dx_um, dy_um)
            if length_um <= 0.0:
                continue
            ex0, ey0 = sent_electric_field[first]
            ex1, ey1 = sent_electric_field[second]
            ex = 0.5 * (ex0 + ex1)
            ey = 0.5 * (ey0 + ey1)
            # Sentaurus vector is V/cm. Project it along the edge and convert to V/m.
            sent_projection = abs((ex * dx_um + ey * dy_um) / length_um) * 100.0
            vela_projection = abs(
                (float(vela[second]["psi"]) - float(vela[first]["psi"]))
                / (length_um * 1.0e-6)
            )
            sent_edge_field.append(sent_projection)
            vela_edge_field.append(vela_projection)

    field_peak = max(sent_edge_field)
    field_results: dict[str, Any] = {
        "comparison": "absolute electric-field projection along matching mesh edges",
        "sentaurus_peak_V_per_m": field_peak,
        "vela_peak_V_per_m": max(vela_edge_field),
        "all_edges": {
            **absolute_metrics(sent_edge_field, vela_edge_field),
            **relative_metrics(sent_edge_field, vela_edge_field),
        },
    }
    for fraction in (0.01, 0.1):
        active = [
            index for index, value in enumerate(sent_edge_field)
            if value >= field_peak * fraction
        ]
        active_reference = [sent_edge_field[index] for index in active]
        active_candidate = [vela_edge_field[index] for index in active]
        field_results[f"above_{fraction:g}_of_sentaurus_peak"] = {
            **absolute_metrics(active_reference, active_candidate),
            **relative_metrics(active_reference, active_candidate),
        }

    sentaurus_voltage = contact_scalar(args.sentaurus_root, "ContactExternalVoltage")
    sentaurus_current = contact_scalar(args.sentaurus_root, "ContactCurrentFlux")
    vela_voltage = max(float(row["phin"]) for row in vela.values())
    result = {
        "operating_point": {
            "sentaurus_voltage_V": sentaurus_voltage,
            "sentaurus_current_A_per_um": sentaurus_current,
            "vela_voltage_V": vela_voltage,
            "vela_current_A_per_um": 0.00010000000616974712,
            "voltage_delta_V": vela_voltage - sentaurus_voltage,
            "voltage_relative_error": (vela_voltage - sentaurus_voltage) / sentaurus_voltage,
        },
        "node_mapping": {
            "semiconductor_nodes": len(semiconductor_nodes),
            "same_node_ids": all(node in vela and node in coordinates for node in semiconductor_nodes),
        },
        "scalar_fields": scalar_results,
        "electric_field": field_results,
        "inputs": {
            "sentaurus_root": str(args.sentaurus_root.resolve()),
            "vela_state": str(args.vela_state.resolve()),
            "mesh": str(args.mesh.resolve()),
        },
    }
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
