#!/usr/bin/env python3
"""Diagnose the SingleDevice Vds=0.1 V, Vg=-0.5 V current mismatch."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Iterable

from sentaurus_import import parse_quoted_list, parse_values_block


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def scalar_field(path: Path) -> dict[int, float]:
    return {
        int(row["node_id"]): float(row["component0"])
        for row in read_rows(path)
    }


def vector_magnitude_field(path: Path) -> dict[int, float]:
    result: dict[int, float] = {}
    for row in read_rows(path):
        components = [
            float(value) for name, value in row.items()
            if name.startswith("component") and value not in (None, "")
        ]
        result[int(row["node_id"])] = math.sqrt(sum(value * value for value in components))
    return result


def percentile(values: list[float], fraction: float) -> float:
    if not values:
        return 0.0
    ordered = sorted(values)
    position = fraction * (len(ordered) - 1)
    lo = int(math.floor(position))
    hi = int(math.ceil(position))
    if lo == hi:
        return ordered[lo]
    weight = position - lo
    return ordered[lo] * (1.0 - weight) + ordered[hi] * weight


def log_comparison(pairs: Iterable[tuple[float, float]]) -> dict[str, float | int]:
    logs = [math.log10(vela / sentaurus) for sentaurus, vela in pairs
            if sentaurus > 0.0 and vela > 0.0]
    absolute = [abs(value) for value in logs]
    return {
        "count": len(logs),
        "median_log10_vela_over_sentaurus": percentile(logs, 0.5),
        "geometric_mean_ratio_vela_over_sentaurus": (
            10.0 ** (sum(logs) / len(logs)) if logs else 0.0),
        "median_abs_log10_error": percentile(absolute, 0.5),
        "p95_abs_log10_error": percentile(absolute, 0.95),
        "max_abs_log10_error": max(absolute, default=0.0),
    }


def weighted_log_ratio(
    node_ids: Iterable[int],
    sentaurus: dict[int, float],
    vela: dict[int, float],
    weights: dict[int, float],
) -> dict[str, float | int]:
    samples: list[tuple[float, float]] = []
    for node in node_ids:
        sent = sentaurus.get(node, 0.0)
        candidate = vela.get(node, 0.0)
        weight = weights.get(node, 0.0)
        if sent > 0.0 and candidate > 0.0 and weight > 0.0:
            samples.append((math.log10(candidate / sent), weight))
    total_weight = sum(weight for _, weight in samples)
    mean = (sum(value * weight for value, weight in samples) / total_weight
            if total_weight > 0.0 else 0.0)
    return {
        "count": len(samples),
        "current_weighted_mean_log10_vela_over_sentaurus": mean,
        "current_weighted_geometric_mean_ratio_vela_over_sentaurus": 10.0 ** mean,
    }


def parse_plt_final(path: Path) -> dict[str, dict[str, float]]:
    text = path.read_text(encoding="utf-8", errors="ignore")
    datasets = parse_quoted_list(text, "datasets")
    rows = parse_values_block(text, len(datasets))
    if not rows:
        raise ValueError(f"no Data rows in {path}")
    row = rows[-1]
    values = dict(zip(datasets, row))
    result: dict[str, dict[str, float]] = {}
    for contact in ("gate", "substrate", "drain", "source"):
        result[contact] = {
            "voltage_V": values[f"{contact} OuterVoltage"],
            "electron_A_per_um": values[f"{contact} eCurrent"],
            "hole_A_per_um": values[f"{contact} hCurrent"],
            "displacement_A_per_um": values[f"{contact} DisplacementCurrent"],
            "total_A_per_um": values[f"{contact} TotalCurrent"],
        }
    return result


def sentaurus_contact_fluxes(export_dir: Path) -> dict[str, float]:
    metadata = json.loads((export_dir / "metadata.json").read_text(encoding="utf-8"))
    regions = {
        str(region["name"]): int(region["index"])
        for region in metadata["regions"] if int(region.get("type", 0)) == 1
    }
    result: dict[str, float] = {}
    for contact in ("gate", "substrate", "drain", "source"):
        field = export_dir / "fields" / f"ContactCurrentFlux_region{regions[contact]}.csv"
        values = scalar_field(field)
        if len(values) != 1:
            raise ValueError(f"expected one ContactCurrentFlux value in {field}")
        result[contact] = next(iter(values.values()))
    return result


def vela_terminal_currents(path: Path) -> dict[str, dict[str, float]]:
    result: dict[str, dict[str, float]] = {}
    for row in read_rows(path):
        result[row["contact"]] = {
            "electron_A_per_um": float(row["current_electron_A_per_um"]),
            "hole_A_per_um": float(row["current_hole_A_per_um"]),
            "total_A_per_um": float(row["current_total_A_per_um"]),
        }
    return result


def point_row(path: Path) -> dict[str, str]:
    rows = read_rows(path)
    if len(rows) != 1:
        raise ValueError(f"expected one sweep point in {path}, got {len(rows)}")
    return rows[0]


def relative_error(reference: float, candidate: float) -> float:
    return abs(candidate - reference) / max(abs(reference), 1.0e-300)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sentaurus-export", required=True, type=Path)
    parser.add_argument("--sentaurus-plt", required=True, type=Path)
    parser.add_argument("--vela-state", required=True, type=Path)
    parser.add_argument("--vela-doping", required=True, type=Path)
    parser.add_argument("--vela-point", required=True, type=Path)
    parser.add_argument("--vela-terminal-balance", required=True, type=Path)
    parser.add_argument("--vela-no-srh-point", required=True, type=Path)
    parser.add_argument("--vela-zero-vds-point", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--aligned-csv", required=True, type=Path)
    args = parser.parse_args()

    fields = args.sentaurus_export / "fields"
    sent_e = scalar_field(fields / "eDensity_region3.csv")
    sent_h = scalar_field(fields / "hDensity_region3.csv")
    sent_srh_si = scalar_field(fields / "srhRecombination_region3.csv")
    sent_srh_poly = scalar_field(fields / "srhRecombination_region4.csv")
    sent_e_current = vector_magnitude_field(fields / "eCurrentDensity_region3.csv")

    state_rows = read_rows(args.vela_state)
    vela_e = {int(row["node_id"]): float(row["electrons_m3"]) / 1.0e6
              for row in state_rows}
    vela_h = {int(row["node_id"]): float(row["holes_m3"]) / 1.0e6
              for row in state_rows}
    doping_rows = read_rows(args.vela_doping)
    donors = {int(row["node_id"]): float(row["donors_cm3"]) for row in doping_rows}
    acceptors = {int(row["node_id"]): float(row["acceptors_cm3"]) for row in doping_rows}
    coordinates = {
        int(row["id"]): (float(row["x_um"]), float(row["y_um"]))
        for row in read_rows(args.sentaurus_export / "nodes.csv")
    }

    silicon_nodes = sorted(sent_e)
    p_type = [node for node in silicon_nodes if acceptors[node] > donors[node]]
    n_type = [node for node in silicon_nodes if donors[node] >= acceptors[node]]
    minority_pairs = [
        (sent_e[node], vela_e[node]) if node in p_type else (sent_h[node], vela_h[node])
        for node in silicon_nodes
    ]

    sentaurus_contacts = parse_plt_final(args.sentaurus_plt)
    sentaurus_tdr_fluxes = sentaurus_contact_fluxes(args.sentaurus_export)
    vela_contacts = vela_terminal_currents(args.vela_terminal_balance)
    vela_point = point_row(args.vela_point)
    no_srh = point_row(args.vela_no_srh_point)
    zero_vds = point_row(args.vela_zero_vds_point)

    sentaurus_drain = sentaurus_contacts["drain"]
    vela_drain = vela_contacts["drain"]
    sentaurus_kcl = sum(item["total_A_per_um"] for item in sentaurus_contacts.values())
    vela_kcl = sum(item["total_A_per_um"] for item in vela_contacts.values())
    vela_current = float(vela_point["current_total_A_per_um"])
    no_srh_current = float(no_srh["current_total_A_per_um"])
    zero_vds_current = float(zero_vds["current_total_A_per_um"])

    sent_srh_values = list(sent_srh_si.values()) + list(sent_srh_poly.values())
    sent_srh_max_m3 = max((abs(value) for value in sent_srh_values), default=0.0) * 1.0e6
    sent_srh_mean_all_nodes_m3 = (
        sum(abs(value) for value in sent_srh_values) / len(state_rows) * 1.0e6)
    vela_srh_max_m3 = float(vela_point["recombination_max_abs_rate_m3_per_s"])
    vela_srh_mean_m3 = float(vela_point["recombination_mean_abs_rate_m3_per_s"])

    electron_gap = relative_error(
        sentaurus_drain["electron_A_per_um"], vela_drain["electron_A_per_um"])
    total_gap = relative_error(
        sentaurus_drain["total_A_per_um"], vela_drain["total_A_per_um"])
    tdr_flux_gap = relative_error(
        sentaurus_tdr_fluxes["drain"], vela_drain["total_A_per_um"])
    sentaurus_output_mismatch = relative_error(
        sentaurus_drain["total_A_per_um"], sentaurus_tdr_fluxes["drain"])
    srh_current_delta = relative_error(vela_current, no_srh_current)
    floor_ratio = abs(zero_vds_current) / max(abs(vela_current), 1.0e-300)
    vela_kcl_ratio = abs(vela_kcl) / max(abs(vela_current), 1.0e-300)

    report = {
        "schema": "vela.singledevice.deep_off_diagnostic.v1",
        "bias": {"Vds_V": 0.1, "Vg_V": -0.5},
        "terminal_currents": {
            "sentaurus": sentaurus_contacts,
            "sentaurus_tdr_contact_flux_A_per_um": sentaurus_tdr_fluxes,
            "vela": vela_contacts,
            "drain_total_relative_error": total_gap,
            "drain_tdr_flux_relative_error": tdr_flux_gap,
            "drain_electron_relative_error": electron_gap,
            "sentaurus_plt_vs_tdr_drain_relative_difference": sentaurus_output_mismatch,
            "sentaurus_kcl_sum_A_per_um": sentaurus_kcl,
            "vela_kcl_sum_A_per_um": vela_kcl,
            "vela_kcl_residual_over_drain_current": vela_kcl_ratio,
        },
        "srh": {
            "sentaurus_max_abs_rate_m3_per_s": sent_srh_max_m3,
            "sentaurus_mean_abs_rate_all_nodes_m3_per_s": sent_srh_mean_all_nodes_m3,
            "vela_max_abs_rate_m3_per_s": vela_srh_max_m3,
            "vela_mean_abs_rate_m3_per_s": vela_srh_mean_m3,
            "vela_no_srh_current_A_per_um": no_srh_current,
            "vela_srh_enabled_current_A_per_um": vela_current,
            "no_srh_current_relative_change": srh_current_delta,
        },
        "minority_carriers": {
            "all_silicon_nodes": log_comparison(minority_pairs),
            "p_type_silicon_electron_density": log_comparison(
                (sent_e[node], vela_e[node]) for node in p_type),
            "n_type_silicon_hole_density": log_comparison(
                (sent_h[node], vela_h[node]) for node in n_type),
            "p_type_electron_density_current_weighted": weighted_log_ratio(
                p_type, sent_e, vela_e, sent_e_current),
        },
        "numerical_floor": {
            "vela_zero_vds_drain_current_A_per_um": zero_vds_current,
            "zero_vds_floor_over_deep_off_current": floor_ratio,
        },
        "classification": {
            "srh_is_dominant": srh_current_delta >= 0.01,
            "contact_balance_is_dominant": vela_kcl_ratio >= 0.5 * total_gap,
            "numerical_floor_is_dominant": floor_ratio >= 0.1,
            "gap_is_electron_transport": abs(electron_gap - total_gap) <= 1.0e-4,
            "relative_error_policy": (
                "retain logarithmic/absolute-current criteria below 1e-13 A/um; "
                "do not use a strict 10% pointwise relative-current gate alone"),
        },
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")

    args.aligned_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.aligned_csv.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow([
            "node_id", "x_um", "y_um", "conductivity_type",
            "donors_cm3", "acceptors_cm3",
            "sentaurus_eDensity_cm3", "vela_eDensity_cm3",
            "sentaurus_hDensity_cm3", "vela_hDensity_cm3",
            "sentaurus_srh_cm3_per_s", "sentaurus_eCurrentDensity_A_per_cm2",
        ])
        for node in silicon_nodes:
            x_um, y_um = coordinates[node]
            writer.writerow([
                node, x_um, y_um, "p" if node in p_type else "n",
                donors[node], acceptors[node], sent_e[node], vela_e[node],
                sent_h[node], vela_h[node], sent_srh_si[node], sent_e_current[node],
            ])

    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
