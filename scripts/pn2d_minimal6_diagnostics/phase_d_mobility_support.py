"""Phase D mobility parameter and support alignment for PN2D Minimal6.

This diagnostic is intentionally read-only with respect to production code.
It compares the sealed Sentaurus element mobility with the documented Vela
Masetti plus Caughey-Thomas law, then propagates both onto the same box-edge
support using the already verified cotangent coefficients.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import statistics
from pathlib import Path
from typing import Iterable

from .mobility_diagnosis import FIELD, MASETTI, field_limited_mobility
from .mobility_diagnosis import masetti_low_field_mobility


CARRIERS = ("electron", "hole")
TOPOLOGIES = ("mirror", "sketch")
BIASES = tuple(-float(value) for value in range(1, 21))


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty CSV {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _quantile(values: Iterable[float], fraction: float) -> float:
    ordered = sorted(float(value) for value in values)
    if not ordered:
        raise ValueError("quantile requires non-empty values")
    position = (len(ordered) - 1) * fraction
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    return ordered[lower] + (ordered[upper] - ordered[lower]) * (
        position - lower
    )


def _abs_dex(candidate: float, reference: float) -> float:
    if candidate <= 0.0 or reference <= 0.0:
        raise ValueError("dex comparison requires positive values")
    return abs(math.log10(candidate / reference))


def _signed_dex(candidate: float, reference: float) -> float:
    if candidate <= 0.0 or reference <= 0.0:
        raise ValueError("dex comparison requires positive values")
    return math.log10(candidate / reference)


def _metric(values: list[float]) -> dict[str, float | int]:
    return {
        "count": len(values),
        "median": statistics.median(values),
        "p95": _quantile(values, 0.95),
        "maximum": max(values),
    }


def _weighted_mean(values: list[float], weights: list[float]) -> float:
    if len(values) != len(weights) or not values:
        raise ValueError("weighted mean requires paired non-empty values")
    denominator = sum(weights)
    if denominator <= 0.0:
        raise ValueError("weighted mean requires positive total weight")
    return sum(value * weight for value, weight in zip(values, weights)) / denominator


def _pair_values(section: str, name: str) -> tuple[float, float]:
    pattern = rf"(?m)^\s*{re.escape(name)}\s*=\s*([^,#\s]+)\s*,\s*([^#\s]+)"
    match = re.search(pattern, section)
    if match is None:
        raise ValueError(f"missing pair parameter {name}")
    return float(match.group(1)), float(match.group(2))


def _section(text: str, name: str) -> str:
    match = re.search(
        rf"(?ms)^{re.escape(name)}:\s*\{{(.*?)^\}}", text
    )
    if match is None:
        raise ValueError(f"missing parameter section {name}")
    return match.group(1)


def _parameter_rows(models_par: Path) -> list[dict[str, object]]:
    text = models_par.read_text(encoding="utf-8", errors="strict")
    constant = _section(text, "ConstantMobility")
    doping = _section(text, "DopingDependence")
    high = _section(text, "HighFieldDependence")

    sent_pairs = {
        "mu_const": tuple(value * 1.0e-4 for value in _pair_values(constant, "mumax")),
        "mu_min1": tuple(value * 1.0e-4 for value in _pair_values(doping, "mumin1")),
        "mu_min2": tuple(value * 1.0e-4 for value in _pair_values(doping, "mumin2")),
        "mu1": tuple(value * 1.0e-4 for value in _pair_values(doping, "mu1")),
        "pc": tuple(value * 1.0e6 for value in _pair_values(doping, "Pc")),
        "cr": tuple(value * 1.0e6 for value in _pair_values(doping, "Cr")),
        "cs": tuple(value * 1.0e6 for value in _pair_values(doping, "Cs")),
        "masetti_alpha": _pair_values(doping, "alpha"),
        "masetti_beta": _pair_values(doping, "beta"),
        "field_beta": _pair_values(high, "beta0"),
        "saturation_velocity": tuple(
            value * 1.0e-2 for value in _pair_values(high, "vsat0")
        ),
    }
    vela_names = {
        "mu_const": "mu_const",
        "mu_min1": "mu_min1",
        "mu_min2": "mu_min2",
        "mu1": "mu1",
        "pc": "pc",
        "cr": "cr",
        "cs": "cs",
        "masetti_alpha": "alpha",
        "masetti_beta": "beta",
    }
    units = {
        "mu_const": "m2_per_Vs",
        "mu_min1": "m2_per_Vs",
        "mu_min2": "m2_per_Vs",
        "mu1": "m2_per_Vs",
        "pc": "m-3",
        "cr": "m-3",
        "cs": "m-3",
        "masetti_alpha": "1",
        "masetti_beta": "1",
        "field_beta": "1",
        "saturation_velocity": "m_per_s",
    }
    result: list[dict[str, object]] = []
    for index, carrier in enumerate(CARRIERS):
        for parameter, sent_values in sent_pairs.items():
            if parameter == "field_beta":
                vela = FIELD[carrier]["beta"]
            elif parameter == "saturation_velocity":
                vela = FIELD[carrier]["saturation_velocity"]
            else:
                vela = MASETTI[carrier][vela_names[parameter]]
            sent = sent_values[index]
            result.append(
                {
                    "carrier": carrier,
                    "parameter": parameter,
                    "unit": units[parameter],
                    "sentaurus_value_si": sent,
                    "vela_value_si": vela,
                    "relative_difference": 0.0
                    if sent == vela == 0.0
                    else abs(vela - sent) / max(abs(sent), 1.0e-300),
                    "status": "exact_match"
                    if math.isclose(vela, sent, rel_tol=1.0e-14, abs_tol=0.0)
                    or sent == vela
                    else "mismatch",
                }
            )
    result.extend(
        [
            {
                "carrier": "both",
                "parameter": "temperature",
                "unit": "K",
                "sentaurus_value_si": 300.0,
                "vela_value_si": 300.0,
                "relative_difference": 0.0,
                "status": "matched_default",
            },
            {
                "carrier": "both",
                "parameter": "high_field_driving_force",
                "unit": "enum",
                "sentaurus_value_si": "GradQuasiFermi_default",
                "vela_value_si": "quasi_fermi_gradient",
                "relative_difference": "",
                "status": "matched_semantics",
            },
            {
                "carrier": "both",
                "parameter": "element_interpolation",
                "unit": "enum",
                "sentaurus_value_si": "not_exposed",
                "vela_value_si": "explicit_controls_below",
                "relative_difference": "",
                "status": "proprietary_default",
            },
        ]
    )
    return result


def _mesh_contract(
    inverse_root: Path,
) -> tuple[
    dict[str, dict[int, tuple[int, int, int]]],
    dict[str, dict[int, float]],
]:
    triangles: dict[str, dict[int, tuple[int, int, int]]] = {}
    doping: dict[str, dict[int, float]] = {}
    for topology in TOPOLOGIES:
        root = inverse_root / "vela" / "source" / "topologies" / topology
        mesh = json.loads((root / "mesh.json").read_text(encoding="utf-8"))
        triangles[topology] = {
            int(cell["id"]): tuple(int(node) for node in cell["node_ids"])
            for cell in mesh["triangles"]
        }
        doping[topology] = {
            int(row["node_id"]): (
                float(row["donors_cm3"]) - float(row["acceptors_cm3"])
            )
            * 1.0e6
            for row in _rows(root / "doping.csv")
        }
        if set(triangles[topology]) != set(range(4)):
            raise ValueError(f"{topology} mesh does not contain cells 0..3")
        if set(doping[topology]) != set(range(6)):
            raise ValueError(f"{topology} doping does not contain nodes 0..5")
    return triangles, doping


def _invert_high_field(
    carrier: str, final_mobility: float, field: float
) -> float | None:
    beta = FIELD[carrier]["beta"]
    velocity = FIELD[carrier]["saturation_velocity"]
    inverse_power = final_mobility ** (-beta) - (field / velocity) ** beta
    if inverse_power <= 0.0:
        return None
    return inverse_power ** (-1.0 / beta)


def _native_rows(
    transport_csv: Path,
    triangles: dict[str, dict[int, tuple[int, int, int]]],
    doping: dict[str, dict[int, float]],
) -> list[dict[str, object]]:
    source = _rows(transport_csv)
    if len(source) != 160:
        raise ValueError("mapped transport input must contain 160 elements")
    output: list[dict[str, object]] = []
    seen: set[tuple[str, float, int]] = set()
    for row in source:
        topology = row["topology"]
        bias = float(row["bias_V"])
        cell = int(row["cell_id"])
        key = (topology, bias, cell)
        if key in seen:
            raise ValueError(f"duplicate mapped element {key}")
        seen.add(key)
        nodes = triangles[topology][cell]
        nodal_doping = [doping[topology][node] for node in nodes]
        mean_doping = sum(nodal_doping) / 3.0
        for carrier in CARRIERS:
            gx = float(row[f"{carrier}_grad_qf_x_V_per_m"])
            gy = float(row[f"{carrier}_grad_qf_y_V_per_m"])
            field = math.hypot(gx, gy)
            sent = float(row[f"{carrier}_mobility_m2_per_Vs"])
            cell_low = masetti_low_field_mobility(carrier, mean_doping)
            cell_final = field_limited_mobility(carrier, cell_low, field)
            node_low_values = [
                masetti_low_field_mobility(carrier, value)
                for value in nodal_doping
            ]
            node_final_values = [
                field_limited_mobility(carrier, value, field)
                for value in node_low_values
            ]
            node_low = sum(node_low_values) / 3.0
            node_final = sum(node_final_values) / 3.0
            inferred_low = _invert_high_field(carrier, sent, field)
            replay = (
                None
                if inferred_low is None
                else field_limited_mobility(carrier, inferred_low, field)
            )
            output.append(
                {
                    "topology": topology,
                    "bias_V": bias,
                    "cell_id": cell,
                    "carrier": carrier,
                    "node_ids": ";".join(str(node) for node in nodes),
                    "cell_average_net_doping_m3": mean_doping,
                    "sentaurus_native_qf_field_V_per_m": field,
                    "vela_cell_average_low_field_m2_per_Vs": cell_low,
                    "vela_cell_average_saturation_factor": cell_final / cell_low,
                    "vela_cell_average_final_m2_per_Vs": cell_final,
                    "vela_node_average_low_field_m2_per_Vs": node_low,
                    "vela_node_average_saturation_factor": node_final / node_low,
                    "vela_node_average_final_m2_per_Vs": node_final,
                    "sentaurus_native_final_m2_per_Vs": sent,
                    "sentaurus_inferred_low_field_m2_per_Vs": ""
                    if inferred_low is None
                    else inferred_low,
                    "sentaurus_inferred_saturation_factor": ""
                    if inferred_low is None
                    else sent / inferred_low,
                    "inferred_low_field_replay_m2_per_Vs": ""
                    if replay is None
                    else replay,
                    "inferred_replay_relative_error": ""
                    if replay is None
                    else abs(replay - sent) / sent,
                    "cell_average_signed_log10_ratio_dex": _signed_dex(
                        cell_final, sent
                    ),
                    "cell_average_abs_log10_error_dex": _abs_dex(
                        cell_final, sent
                    ),
                    "node_average_abs_log10_error_dex": _abs_dex(
                        node_final, sent
                    ),
                    "status": "valid" if inferred_low is not None else "noninvertible",
                }
            )
    if len(output) != 320:
        raise ValueError("native carrier-element output must contain 320 rows")
    return output


def _local_index(path: Path) -> dict[
    tuple[str, float, int, str, tuple[int, int]], dict[str, str]
]:
    output: dict[
        tuple[str, float, int, str, tuple[int, int]], dict[str, str]
    ] = {}
    for row in _rows(path):
        pair = tuple(sorted((int(row["node0"]), int(row["node1"]))))
        key = (
            row["topology"],
            float(row["bias_V"]),
            int(row["cell_id"]),
            row["carrier"],
            pair,
        )
        output[key] = row
    if len(output) != 960:
        raise ValueError("local mobility input must contain 960 rows")
    return output


def _global_index(path: Path) -> dict[
    tuple[str, float, str, tuple[int, int]], dict[str, str]
]:
    output: dict[
        tuple[str, float, str, tuple[int, int]], dict[str, str]
    ] = {}
    for row in _rows(path):
        pair = tuple(sorted((int(row["node0"]), int(row["node1"]))))
        key = (row["topology"], float(row["bias_V"]), row["carrier"], pair)
        output[key] = row
    if len(output) != 720:
        raise ValueError("global mobility input must contain 720 rows")
    return output


def _current_index(path: Path) -> dict[
    tuple[str, float, str, str, tuple[int, int]], dict[str, str]
]:
    output = {}
    for row in _rows(path):
        pair = tuple(sorted((int(row["node0"]), int(row["node1"]))))
        output[
            (
                row["topology"],
                float(row["bias_V"]),
                row["stage"],
                row["carrier"],
                pair,
            )
        ] = row
    return output


def _edge_rows(
    native_rows: list[dict[str, object]],
    geometry_csv: Path,
    local_csv: Path,
    global_csv: Path,
    stage_csv: Path,
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    native = {
        (
            str(row["topology"]),
            float(row["bias_V"]),
            int(row["cell_id"]),
            str(row["carrier"]),
        ): row
        for row in native_rows
    }
    geometry: dict[str, dict[tuple[int, int], list[dict[str, str]]]] = {
        topology: {} for topology in TOPOLOGIES
    }
    edge_ids: dict[tuple[str, tuple[int, int]], int] = {}
    for row in _rows(geometry_csv):
        topology = row["topology"]
        pair = tuple(sorted((int(row["node0"]), int(row["node1"]))))
        geometry[topology].setdefault(pair, []).append(row)
    for topology in TOPOLOGIES:
        for edge_id, pair in enumerate(sorted(geometry[topology])):
            edge_ids[(topology, pair)] = edge_id
    local = _local_index(local_csv)
    global_values = _global_index(global_csv)
    currents = _current_index(stage_csv)

    output: list[dict[str, object]] = []
    adjacent_output: list[dict[str, object]] = []
    for topology in TOPOLOGIES:
        for bias in BIASES:
            for pair, adjacent in sorted(geometry[topology].items()):
                kappa_sum = sum(float(row["kappa"]) for row in adjacent)
                for carrier in CARRIERS:
                    global_row = global_values[(topology, bias, carrier, pair)]
                    global_raw = global_row["baseline_mobility_m2_per_Vs"].strip()
                    global_mu = None if not global_raw else float(global_raw)
                    sent_weighted = None
                    cell_weighted = None
                    local_weighted = None
                    if kappa_sum != 0.0:
                        sent_weighted = sum(
                            float(item["kappa"])
                            * float(
                                native[
                                    (
                                        topology,
                                        bias,
                                        int(item["element"]),
                                        carrier,
                                    )
                                ]["sentaurus_native_final_m2_per_Vs"]
                            )
                            for item in adjacent
                        ) / kappa_sum
                        cell_weighted = sum(
                            float(item["kappa"])
                            * float(
                                native[
                                    (
                                        topology,
                                        bias,
                                        int(item["element"]),
                                        carrier,
                                    )
                                ]["vela_cell_average_final_m2_per_Vs"]
                            )
                            for item in adjacent
                        ) / kappa_sum
                        local_weighted = sum(
                            float(item["kappa"])
                            * float(
                                local[
                                    (
                                        topology,
                                        bias,
                                        int(item["element"]),
                                        carrier,
                                        pair,
                                    )
                                ]["direct_cpp_mobility_m2_per_Vs"]
                            )
                            for item in adjacent
                        ) / kappa_sum
                    current = currents.get(
                        (
                            topology,
                            bias,
                            "sentaurus_qfp_density",
                            carrier,
                            pair,
                        )
                    )
                    reference_current = (
                        None
                        if current is None
                        else abs(float(current["reference_A_per_um"]))
                    )
                    candidate_current = (
                        None
                        if current is None
                        else abs(float(current["candidate_A_per_um"]))
                    )
                    current_error = (
                        None
                        if reference_current is None
                        or candidate_current is None
                        or reference_current <= 0.0
                        or candidate_current <= 0.0
                        else _abs_dex(candidate_current, reference_current)
                    )
                    status = (
                        "geometric_zero"
                        if kappa_sum == 0.0
                        else "missing_global_mobility"
                        if global_mu is None
                        else "valid"
                    )
                    model_signed = (
                        None
                        if sent_weighted is None or cell_weighted is None
                        else _signed_dex(cell_weighted, sent_weighted)
                    )
                    support_signed = (
                        None
                        if cell_weighted is None or global_mu is None
                        else _signed_dex(global_mu, cell_weighted)
                    )
                    total_signed = (
                        None
                        if sent_weighted is None or global_mu is None
                        else _signed_dex(global_mu, sent_weighted)
                    )
                    output.append(
                        {
                            "topology": topology,
                            "bias_V": bias,
                            "carrier": carrier,
                            "edge_id": edge_ids[(topology, pair)],
                            "node0": pair[0],
                            "node1": pair[1],
                            "adjacent_element_ids": ";".join(
                                str(item["element"]) for item in adjacent
                            ),
                            "kappa_values": ";".join(
                                item["kappa"] for item in adjacent
                            ),
                            "kappa_sum": kappa_sum,
                            "sentaurus_box_edge_mobility_m2_per_Vs": ""
                            if sent_weighted is None
                            else sent_weighted,
                            "vela_native_cell_box_edge_mobility_m2_per_Vs": ""
                            if cell_weighted is None
                            else cell_weighted,
                            "vela_triangle_local_box_edge_mobility_m2_per_Vs": ""
                            if local_weighted is None
                            else local_weighted,
                            "vela_production_global_edge_mobility_m2_per_Vs": ""
                            if global_mu is None
                            else global_mu,
                            "model_signed_log10_ratio_dex": ""
                            if model_signed is None
                            else model_signed,
                            "support_signed_log10_ratio_dex": ""
                            if support_signed is None
                            else support_signed,
                            "total_signed_log10_ratio_dex": ""
                            if total_signed is None
                            else total_signed,
                            "decomposition_closure_dex": ""
                            if None in (model_signed, support_signed, total_signed)
                            else total_signed - model_signed - support_signed,
                            "reference_abs_current_A_per_um": ""
                            if reference_current is None
                            else reference_current,
                            "vela_global_candidate_abs_current_A_per_um": ""
                            if candidate_current is None
                            else candidate_current,
                            "current_abs_log10_error_dex": ""
                            if current_error is None
                            else current_error,
                            "status": status,
                        }
                    )
                    for item in adjacent:
                        cell = int(item["element"])
                        native_row = native[(topology, bias, cell, carrier)]
                        local_row = local[(topology, bias, cell, carrier, pair)]
                        adjacent_output.append(
                            {
                                "topology": topology,
                                "bias_V": bias,
                                "carrier": carrier,
                                "edge_id": edge_ids[(topology, pair)],
                                "node0": pair[0],
                                "node1": pair[1],
                                "cell_id": cell,
                                "local_edge": int(item["local_edge"]),
                                "kappa": float(item["kappa"]),
                                "sentaurus_native_qf_field_V_per_m": native_row[
                                    "sentaurus_native_qf_field_V_per_m"
                                ],
                                "sentaurus_native_mobility_m2_per_Vs": native_row[
                                    "sentaurus_native_final_m2_per_Vs"
                                ],
                                "vela_native_cell_mobility_m2_per_Vs": native_row[
                                    "vela_cell_average_final_m2_per_Vs"
                                ],
                                "vela_triangle_local_edge_qf_field_V_per_m": float(
                                    local_row["edge_qf_field_V_per_m"]
                                ),
                                "vela_triangle_local_mobility_m2_per_Vs": float(
                                    local_row["direct_cpp_mobility_m2_per_Vs"]
                                ),
                                "vela_production_global_edge_mobility_m2_per_Vs": ""
                                if global_mu is None
                                else global_mu,
                                "reference_abs_current_A_per_um": ""
                                if reference_current is None
                                else reference_current,
                                "vela_global_candidate_abs_current_A_per_um": ""
                                if candidate_current is None
                                else candidate_current,
                                "current_abs_log10_error_dex": ""
                                if current_error is None
                                else current_error,
                                "status": "geometric_zero"
                                if float(item["kappa"]) == 0.0
                                else "valid",
                            }
                        )
    if len(output) != 720 or len(adjacent_output) != 960:
        raise ValueError("unexpected edge or adjacent-element row count")
    return output, adjacent_output


def _summary_rows(
    native: list[dict[str, object]],
    edges: list[dict[str, object]],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for carrier in CARRIERS:
        carrier_native = [
            row for row in native if row["carrier"] == carrier
        ]
        for branch, column in (
            ("cell_average", "cell_average_abs_log10_error_dex"),
            ("node_average", "node_average_abs_log10_error_dex"),
        ):
            metric = _metric([float(row[column]) for row in carrier_native])
            rows.append(
                {
                    "support": "sentaurus_native_element",
                    "carrier": carrier,
                    "branch": branch,
                    "sample_count": metric["count"],
                    "median_abs_dex": metric["median"],
                    "p95_abs_dex": metric["p95"],
                    "maximum_abs_dex": metric["maximum"],
                    "absolute_current_weighted_mean_abs_dex": "",
                }
            )
        carrier_edges = [
            row
            for row in edges
            if row["carrier"] == carrier and row["status"] == "valid"
        ]
        for branch, column in (
            ("native_cell_model", "model_signed_log10_ratio_dex"),
            ("global_support", "support_signed_log10_ratio_dex"),
            ("production_total", "total_signed_log10_ratio_dex"),
        ):
            values = [abs(float(row[column])) for row in carrier_edges]
            weights = [
                float(row["reference_abs_current_A_per_um"])
                for row in carrier_edges
            ]
            metric = _metric(values)
            rows.append(
                {
                    "support": "active_box_edge",
                    "carrier": carrier,
                    "branch": branch,
                    "sample_count": metric["count"],
                    "median_abs_dex": metric["median"],
                    "p95_abs_dex": metric["p95"],
                    "maximum_abs_dex": metric["maximum"],
                    "absolute_current_weighted_mean_abs_dex": _weighted_mean(
                        values, weights
                    ),
                }
            )
    return rows


def _central_rows(
    adjacent: list[dict[str, object]],
    edges: list[dict[str, object]],
) -> list[dict[str, object]]:
    edge_index = {
        (
            str(row["topology"]),
            float(row["bias_V"]),
            str(row["carrier"]),
            int(row["node0"]),
            int(row["node1"]),
        ): row
        for row in edges
    }
    output = []
    for row in adjacent:
        if (int(row["node0"]), int(row["node1"])) != (1, 5):
            continue
        aggregate = edge_index[
            (
                str(row["topology"]),
                float(row["bias_V"]),
                str(row["carrier"]),
                1,
                5,
            )
        ]
        output.append(
            {
                **row,
                "box_model_signed_log10_ratio_dex": aggregate[
                    "model_signed_log10_ratio_dex"
                ],
                "box_support_signed_log10_ratio_dex": aggregate[
                    "support_signed_log10_ratio_dex"
                ],
                "box_total_signed_log10_ratio_dex": aggregate[
                    "total_signed_log10_ratio_dex"
                ],
            }
        )
    if len(output) != 160:
        raise ValueError("central edge decomposition must contain 160 rows")
    return output


def _control_rows(
    native: list[dict[str, object]],
) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    for carrier in CARRIERS:
        carrier_rows = [row for row in native if row["carrier"] == carrier]
        documented_errors = [
            float(row["cell_average_abs_log10_error_dex"])
            for row in carrier_rows
        ]
        documented = _metric(documented_errors)
        replay_errors = [
            _abs_dex(
                float(row["inferred_low_field_replay_m2_per_Vs"]),
                float(row["sentaurus_native_final_m2_per_Vs"]),
            )
            for row in carrier_rows
            if row["inferred_low_field_replay_m2_per_Vs"] != ""
        ]
        replay = _metric(replay_errors)
        output.extend(
            [
                {
                    "carrier": carrier,
                    "control": "sentaurus_documented_parameter_substitution",
                    "sample_count": documented["count"],
                    "median_abs_dex": documented["median"],
                    "p95_abs_dex": documented["p95"],
                    "maximum_abs_dex": documented["maximum"],
                    "maximum_change_from_vela_current_parameters": 0.0,
                    "production_candidate": "yes",
                    "interpretation": "no_effect_parameters_already_identical",
                },
                {
                    "carrier": carrier,
                    "control": "inferred_low_field_documented_high_field_replay",
                    "sample_count": replay["count"],
                    "median_abs_dex": replay["median"],
                    "p95_abs_dex": replay["p95"],
                    "maximum_abs_dex": replay["maximum"],
                    "maximum_change_from_vela_current_parameters": "",
                    "production_candidate": "no",
                    "interpretation": "algebraic_closure_only",
                },
            ]
        )
    return output


def _markdown(
    parameters: list[dict[str, object]],
    summary: list[dict[str, object]],
    central: list[dict[str, object]],
    controls: list[dict[str, object]],
    outcome: dict[str, object],
) -> str:
    matched = sum(row["status"] in ("exact_match", "matched_default", "matched_semantics") for row in parameters)
    numeric = sum(row["status"] == "exact_match" for row in parameters)
    lines = [
        "# PN2D Minimal6 Phase D mobility and support alignment",
        "",
        "## Result",
        "",
        f"- Status: `{outcome['status']}`.",
        f"- Primary outcome: `{outcome['primary_outcome']}`.",
        f"- Additional outcome: `{outcome['secondary_outcome']}`.",
        f"- Sealed parameter entries: {len(parameters)}; matched entries: {matched}; exact numeric matches: {numeric}.",
        "- No production formula or fitted parameter was used.",
        "",
        "## Common-support metrics",
        "",
        "| Support | Carrier | Branch | N | Median dex | P95 dex | Max dex | Current-weighted mean dex |",
        "|---|---|---|---:|---:|---:|---:|---:|",
    ]
    for row in summary:
        weighted = row["absolute_current_weighted_mean_abs_dex"]
        lines.append(
            f"| {row['support']} | {row['carrier']} | {row['branch']} | "
            f"{row['sample_count']} | {float(row['median_abs_dex']):.6f} | "
            f"{float(row['p95_abs_dex']):.6f} | "
            f"{float(row['maximum_abs_dex']):.6f} | "
            f"{'' if weighted == '' else f'{float(weighted):.6f}'} |"
        )
    high = [
        row
        for row in central
        if float(row["bias_V"]) == -20.0 and int(row["cell_id"]) in (1, 2)
    ]
    lines.extend(
        [
            "",
            "## Central edge 1-5 at -20 V",
            "",
            "| Topology | Carrier | Cell | Kappa | Sent field V/m | Sent mobility | Vela cell mobility | Vela local mobility | Vela global mobility | Reference abs A/um | Candidate abs A/um | Current error dex |",
            "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
    )
    for row in high:
        lines.append(
            f"| {row['topology']} | {row['carrier']} | {row['cell_id']} | "
            f"{float(row['kappa']):.3f} | "
            f"{float(row['sentaurus_native_qf_field_V_per_m']):.6e} | "
            f"{float(row['sentaurus_native_mobility_m2_per_Vs']):.6e} | "
            f"{float(row['vela_native_cell_mobility_m2_per_Vs']):.6e} | "
            f"{float(row['vela_triangle_local_mobility_m2_per_Vs']):.6e} | "
            f"{float(row['vela_production_global_edge_mobility_m2_per_Vs']):.6e} | "
            f"{float(row['reference_abs_current_A_per_um']):.6e} | "
            f"{float(row['vela_global_candidate_abs_current_A_per_um']):.6e} | "
            f"{float(row['current_abs_log10_error_dex']):.6f} |"
        )
    lines.extend(
        [
            "",
            "## No-fit parameter controls",
            "",
            "| Carrier | Control | N | Median dex | P95 dex | Max dex | Production candidate | Interpretation |",
            "|---|---|---:|---:|---:|---:|---|---|",
        ]
    )
    for row in controls:
        lines.append(
            f"| {row['carrier']} | {row['control']} | {row['sample_count']} | "
            f"{float(row['median_abs_dex']):.6e} | "
            f"{float(row['p95_abs_dex']):.6e} | "
            f"{float(row['maximum_abs_dex']):.6e} | "
            f"{row['production_candidate']} | {row['interpretation']} |"
        )
    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "- The documented numeric mobility parameters are already identical after unit conversion, so the no-fit Sentaurus parameter substitution changes no value.",
            "- The native-element parity target is evaluated on all 320 carrier-element samples. Any residual beyond the target is therefore not attributable to the sealed numeric parameter table.",
            "- `native_cell_model` is the residual after applying the Vela law directly to the native Sentaurus element field. `global_support` is the additional change caused by replacing coefficient-weighted cell mobility with one Vela global-edge mobility. Signed terms close algebraically to `production_total`.",
            "- The inferred-low-field replay is an algebraic inversion control, not a fitted production model. It proves that the exported final mobility and documented high-field form can be replayed, but it does not reveal Sentaurus's proprietary element interpolation.",
            "",
            "## Decision",
            "",
            f"`{outcome['primary_outcome']}` is the primary typed result because the native-element target is not met despite exact documented parameters. `{outcome['secondary_outcome']}` is retained because element-to-global-edge aggregation is separately nonzero.",
            "",
        ]
    )
    return "\n".join(lines)


def run_phase_d(
    *,
    models_par: str | Path,
    sdevice_cmd: str | Path,
    vela_deck: str | Path,
    mobility_source: str | Path,
    inverse_inputs_root: str | Path,
    mapped_transport_csv: str | Path,
    cell_mapping_csv: str | Path,
    geometry_csv: str | Path,
    local_mobility_csv: str | Path,
    global_mobility_csv: str | Path,
    stage_edge_csv: str | Path,
    output_root: str | Path,
) -> dict[str, object]:
    paths = {
        "models_par": Path(models_par).resolve(),
        "sdevice_cmd": Path(sdevice_cmd).resolve(),
        "vela_deck": Path(vela_deck).resolve(),
        "mobility_source": Path(mobility_source).resolve(),
        "mapped_transport_csv": Path(mapped_transport_csv).resolve(),
        "cell_mapping_csv": Path(cell_mapping_csv).resolve(),
        "geometry_csv": Path(geometry_csv).resolve(),
        "local_mobility_csv": Path(local_mobility_csv).resolve(),
        "global_mobility_csv": Path(global_mobility_csv).resolve(),
        "stage_edge_csv": Path(stage_edge_csv).resolve(),
    }
    inverse_root = Path(inverse_inputs_root).resolve()
    output = Path(output_root).resolve()
    output.mkdir(parents=True, exist_ok=True)
    for name, path in paths.items():
        if not path.is_file():
            raise FileNotFoundError(f"missing {name}: {path}")

    command = paths["sdevice_cmd"].read_text(encoding="utf-8")
    required_tokens = (
        "DopingDependence",
        "HighFieldSaturation",
        "eGradQuasiFermi/Element/Vector",
        "hGradQuasiFermi/Element/Vector",
        "eMobility/Element",
        "hMobility/Element",
    )
    if any(token not in command for token in required_tokens):
        raise ValueError("Sentaurus command does not seal the Phase D fields/models")
    deck = json.loads(paths["vela_deck"].read_text(encoding="utf-8"))
    mobility_config = deck["solver"]["mobility"]
    if mobility_config != {
        "high_field_driving_force": "quasi_fermi_gradient",
        "model": "masetti_field",
    }:
        raise ValueError("Vela mobility configuration is not the sealed Phase D config")

    mapping = _rows(paths["cell_mapping_csv"])
    if len(mapping) != 8:
        raise ValueError("cell mapping must contain eight topology-cell rows")
    mapping_max = max(float(row["electric_field_relative_residual"]) for row in mapping)
    if mapping_max > 1.0e-12:
        raise ValueError("cell mapping electric-field gate failed")

    parameters = _parameter_rows(paths["models_par"])
    numeric_mismatch = [row for row in parameters if row["status"] == "mismatch"]
    if numeric_mismatch:
        raise ValueError(f"documented mobility parameter mismatch: {numeric_mismatch}")
    triangles, doping = _mesh_contract(inverse_root)
    native = _native_rows(paths["mapped_transport_csv"], triangles, doping)
    edges, adjacent = _edge_rows(
        native,
        paths["geometry_csv"],
        paths["local_mobility_csv"],
        paths["global_mobility_csv"],
        paths["stage_edge_csv"],
    )
    central = _central_rows(adjacent, edges)
    summary = _summary_rows(native, edges)
    controls = _control_rows(native)

    native_cell = [
        row for row in summary
        if row["support"] == "sentaurus_native_element"
        and row["branch"] == "cell_average"
    ]
    parity = all(
        float(row["median_abs_dex"]) <= 0.03
        and float(row["p95_abs_dex"]) <= 0.10
        for row in native_cell
    )
    support_nonzero = any(
        abs(float(row["support_signed_log10_ratio_dex"])) > 1.0e-12
        for row in edges if row["status"] == "valid"
    )
    max_closure = max(
        abs(float(row["decomposition_closure_dex"]))
        for row in edges if row["status"] == "valid"
    )
    max_replay = max(
        float(row["inferred_replay_relative_error"])
        for row in native if row["inferred_replay_relative_error"] != ""
    )
    outcome = {
        "status": "valid",
        "primary_outcome": "parity_candidate"
        if parity
        else "proprietary_model_difference",
        "secondary_outcome": "support_mismatch"
        if support_nonzero
        else "support_matched",
        "parameter_mismatch": False,
        "native_element_parity_target_passed": parity,
    }

    outputs = {
        "parameter_comparison.csv": parameters,
        "native_element_decomposition.csv": native,
        "box_edge_mobility_decomposition.csv": edges,
        "box_edge_adjacent_elements.csv": adjacent,
        "central_edge_1_5_decomposition.csv": central,
        "parameter_substitution_controls.csv": controls,
        "summary.csv": summary,
    }
    for name, rows in outputs.items():
        _write_csv(output / name, rows)
    report = _markdown(parameters, summary, central, controls, outcome)
    (output / "report.md").write_text(report, encoding="utf-8", newline="\n")

    manifest: dict[str, object] = {
        "schema_version": 1,
        "status": "valid",
        "experiment": "pn2d_minimal6_phase_d_mobility_support_alignment",
        "outcome": outcome,
        "contracts": {
            "state_count": 40,
            "native_element_count": 160,
            "native_carrier_element_count": 320,
            "global_carrier_edge_count": 720,
            "active_box_edge_count": sum(
                row["status"] == "valid" for row in edges
            ),
            "geometric_zero_box_edge_count": sum(
                row["status"] == "geometric_zero" for row in edges
            ),
            "adjacent_carrier_element_edge_count": 960,
            "central_edge_adjacent_count": 160,
            "parameter_substitution_control_count": 4,
        },
        "gates": {
            "cell_mapping_max_relative_residual": mapping_max,
            "parameter_numeric_mismatch_count": len(numeric_mismatch),
            "native_inferred_replay_max_relative_error": max_replay,
            "edge_decomposition_max_abs_closure_dex": max_closure,
            "native_element_parity_target_passed": parity,
        },
        "limitations": [
            "Sentaurus native directed-edge current is unavailable",
            "Sentaurus element interpolation details are not exposed by the parameter file",
            "the inferred low-field branch is an algebraic diagnostic and not a fitted production candidate",
        ],
        "inputs": {
            name: {"path": str(path), "sha256": _sha256(path)}
            for name, path in paths.items()
        },
        "inverse_inputs_root": str(inverse_root),
        "outputs": {},
    }
    for name in outputs:
        manifest["outputs"][name] = _sha256(output / name)
    manifest["outputs"]["report.md"] = _sha256(output / "report.md")
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8", newline="\n"
    )
    return manifest
