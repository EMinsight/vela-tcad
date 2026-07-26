#!/usr/bin/env python3
"""Replay Vela upstream physics on imported general-Tri3 Sentaurus states."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.diagnose_pn2d_general_tri3_element_edge_avalanche import (
    geometry_rows,
    parse_log,
)
from scripts.pn2d_general_tri3_contract import (
    EXACT_BIASES_V,
    SCHEMA_ID,
    SENTAURUS_RELEASE,
)


NI_300K_CM3 = 14638914958.767616
VT_300K_V = 0.025851999786435
IMPURITY_CM3 = 1.0e17
SLOTBOOM_REFERENCE_CM3 = 1.0e17
SLOTBOOM_COEFFICIENT_EV = 9.0e-3
SLOTBOOM_SMOOTHING = 0.5

MASETTI = {
    "electron": {
        "mu_const": 1417.0,
        "mu_min1": 52.2,
        "mu_min2": 52.2,
        "mu1": 43.4,
        "pc": 0.0,
        "cr": 9.68e16,
        "cs": 3.43e20,
        "alpha": 0.68,
        "beta": 2.0,
        "vsat": 1.07e7,
        "field_beta": 1.109,
    },
    "hole": {
        "mu_const": 470.5,
        "mu_min1": 44.9,
        "mu_min2": 0.0,
        "mu1": 29.0,
        "pc": 9.23e16,
        "cr": 2.23e17,
        "cs": 6.10e20,
        "alpha": 0.719,
        "beta": 2.0,
        "vsat": 8.37e6,
        "field_beta": 1.213,
    },
}

VAN_OVERSTRAETEN = {
    "electron": {
        "a_low": 7.03e5,
        "a_high": 7.03e5,
        "b_low": 1.231e6,
        "b_high": 1.231e6,
    },
    "hole": {
        "a_low": 1.582e6,
        "a_high": 6.71e5,
        "b_low": 2.036e6,
        "b_high": 1.693e6,
    },
}
ALPHA_SWITCH_FIELD_V_CM = 4.0e5


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    return parser.parse_args()


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty CSV: {path}")
    columns: list[str] = []
    for row in rows:
        for key in row:
            if key not in columns:
                columns.append(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def abs_dex(candidate: float, reference: float) -> float | None:
    if candidate <= 0.0 or reference <= 0.0:
        return None
    return abs(math.log10(candidate / reference))


def percentile(values: list[float], fraction: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    position = fraction * (len(ordered) - 1)
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def error_summary(values: list[float | None]) -> dict[str, Any]:
    finite = [
        value for value in values
        if value is not None and math.isfinite(value)
    ]
    return {
        "finite_count": len(finite),
        "median": statistics.median(finite) if finite else None,
        "p95": percentile(finite, 0.95),
        "maximum": max(finite) if finite else None,
    }


def node_doping(x_um: float, junction_um: float) -> tuple[float, float]:
    if abs(x_um - junction_um) <= 1.0e-12:
        return 0.0, 2.0 * IMPURITY_CM3
    return (
        (-IMPURITY_CM3 if x_um < junction_um else IMPURITY_CM3),
        IMPURITY_CM3,
    )


def effective_ni(total_impurity_cm3: float) -> tuple[float, float]:
    x = math.log(total_impurity_cm3 / SLOTBOOM_REFERENCE_CM3)
    delta_ev = SLOTBOOM_COEFFICIENT_EV * (
        x + math.sqrt(x * x + SLOTBOOM_SMOOTHING)
    )
    return (
        NI_300K_CM3 * math.exp(delta_ev / (2.0 * VT_300K_V)),
        delta_ev,
    )


def carrier_density(
    ni_cm3: float,
    psi_v: float,
    qfp_v: float,
    carrier: str,
) -> float:
    argument = (
        (psi_v - qfp_v) / VT_300K_V
        if carrier == "electron"
        else (qfp_v - psi_v) / VT_300K_V
    )
    return ni_cm3 * math.exp(max(-500.0, min(500.0, argument)))


def p1_gradient(
    points: list[tuple[float, float]],
    values: list[float],
) -> tuple[float, float]:
    (x0, y0), (x1, y1), (x2, y2) = points
    denominator = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)
    if denominator == 0.0:
        raise ValueError("degenerate Tri3")
    gx = (
        (values[1] - values[0]) * (y2 - y0)
        - (values[2] - values[0]) * (y1 - y0)
    ) / denominator
    gy = (
        (x1 - x0) * (values[2] - values[0])
        - (x2 - x0) * (values[1] - values[0])
    ) / denominator
    return gx, gy


def vector_angle_deg(
    candidate: tuple[float, float],
    reference: tuple[float, float],
) -> float | None:
    candidate_mag = math.hypot(*candidate)
    reference_mag = math.hypot(*reference)
    if candidate_mag == 0.0 or reference_mag == 0.0:
        return None
    cosine = (
        candidate[0] * reference[0] + candidate[1] * reference[1]
    ) / (candidate_mag * reference_mag)
    return math.degrees(math.acos(max(-1.0, min(1.0, cosine))))


def masetti(net_doping_cm3: float, carrier: str) -> float:
    params = MASETTI[carrier]
    doping = abs(net_doping_cm3)
    if doping <= 0.0:
        return params["mu_const"]
    exponential = params["mu_min1"] * math.exp(
        -max(0.0, params["pc"]) / doping
    )
    rolloff = (params["mu_const"] - params["mu_min2"]) / (
        1.0 + (doping / params["cr"]) ** params["alpha"]
    )
    correction = params["mu1"] / (
        1.0 + (params["cs"] / doping) ** params["beta"]
    )
    return max(0.0, exponential + rolloff - correction)


def field_limit(
    low_field_mobility: float,
    field_v_cm: float,
    carrier: str,
) -> float:
    params = MASETTI[carrier]
    ratio = low_field_mobility * abs(field_v_cm) / params["vsat"]
    return low_field_mobility / (
        1.0 + ratio ** params["field_beta"]
    ) ** (1.0 / params["field_beta"])


def alpha(field_v_cm: float, carrier: str) -> float:
    field = abs(field_v_cm)
    if field <= 0.0:
        return 0.0
    branch = "low" if field < ALPHA_SWITCH_FIELD_V_CM else "high"
    params = VAN_OVERSTRAETEN[carrier]
    return params[f"a_{branch}"] * math.exp(
        -params[f"b_{branch}"] / field
    )


def main() -> int:
    args = parse_args()
    raw_root = args.raw_root.resolve()
    output_root = args.output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    raw_manifest_path = raw_root / "manifest.json"
    raw_manifest = json.loads(raw_manifest_path.read_text(encoding="ascii"))
    if raw_manifest.get("schema") != SCHEMA_ID:
        raise ValueError("raw schema mismatch")
    if raw_manifest.get("status") != "passed":
        raise ValueError("raw root is not complete")
    if raw_manifest.get("sentaurus_release") != SENTAURUS_RELEASE:
        raise ValueError("Sentaurus release mismatch")
    if tuple(raw_manifest.get("exact_biases_V", ())) != EXACT_BIASES_V:
        raise ValueError("exact bias mismatch")
    case_names = tuple(raw_manifest["cases"])
    if len(case_names) != 1:
        raise ValueError(f"expected one case, got {case_names}")
    case_name = case_names[0]
    log_path = (
        raw_root
        / case_name
        / "implicit_default"
        / "fetched"
        / "run_implicit_default.out"
    )
    groups = parse_log(log_path)
    geometry = geometry_rows(groups, EXACT_BIASES_V[0])
    geometry_by_element = {
        int(row["element"]): row for row in geometry
    }

    vertices = {
        (float(row["bias_V"]), int(row["vertex"])): row
        for row in groups["vertices"]
    }
    local_vertices: dict[tuple[float, int], list[tuple[int, int]]] = (
        defaultdict(list)
    )
    used_vertices: set[tuple[float, int]] = set()
    for row in groups["measures"]:
        key = (float(row["bias_V"]), int(row["element"]))
        vertex = int(row["vertex"])
        local_vertices[key].append((int(row["local_vertex"]), vertex))
        used_vertices.add((key[0], vertex))
    x_coordinates = [
        float(vertices[key]["x_um"]) for key in used_vertices
    ]
    junction_um = 0.5 * (min(x_coordinates) + max(x_coordinates))

    density_rows: list[dict[str, Any]] = []
    recomputed_density: dict[
        tuple[float, int], tuple[float, float]
    ] = {}
    for key in sorted(used_vertices):
        row = vertices[key]
        x_um = float(row["x_um"])
        net_doping, total_impurity = node_doping(x_um, junction_um)
        ni_eff, delta_ev = effective_ni(total_impurity)
        n_cm3 = carrier_density(
            ni_eff,
            float(row["psi_V"]),
            float(row["eQFP_V"]),
            "electron",
        )
        p_cm3 = carrier_density(
            ni_eff,
            float(row["psi_V"]),
            float(row["hQFP_V"]),
            "hole",
        )
        recomputed_density[key] = (n_cm3, p_cm3)
        density_rows.append(
            {
                "case": case_name,
                "bias_V": key[0],
                "vertex": key[1],
                "x_um": x_um,
                "y_um": float(row["y_um"]),
                "net_doping_cm3": net_doping,
                "total_impurity_cm3": total_impurity,
                "slotboom_delta_eV": delta_ev,
                "vela_ni_eff_cm3": ni_eff,
                "sentaurus_n_cm3": float(row["n_cm3"]),
                "vela_n_cm3": n_cm3,
                "electron_absolute_error_dex": abs_dex(
                    n_cm3, float(row["n_cm3"])
                ),
                "sentaurus_p_cm3": float(row["p_cm3"]),
                "vela_p_cm3": p_cm3,
                "hole_absolute_error_dex": abs_dex(
                    p_cm3, float(row["p_cm3"])
                ),
                "observation_label": "native_node_vs_vela_recomputation",
            }
        )

    vector_rows: list[dict[str, Any]] = []
    mobility_rows: list[dict[str, Any]] = []
    alpha_rows: list[dict[str, Any]] = []
    element_rows = {
        (float(row["bias_V"]), int(row["element"])): row
        for row in groups["elements"]
    }
    edge_rows: dict[tuple[float, int], list[dict[str, Any]]] = defaultdict(list)
    for row in groups["edges"]:
        edge_rows[(float(row["bias_V"]), int(row["element"]))].append(row)

    vector_specs = (
        (
            "electric_field",
            "psi_V",
            "efield_x_V_cm",
            "efield_y_V_cm",
        ),
        (
            "electron_qfp_gradient",
            "eQFP_V",
            "grad_qf_n_x_V_cm",
            "grad_qf_n_y_V_cm",
        ),
        (
            "hole_qfp_gradient",
            "hQFP_V",
            "grad_qf_p_x_V_cm",
            "grad_qf_p_y_V_cm",
        ),
    )
    for key, element in sorted(element_rows.items()):
        bias, element_id = key
        ids = [
            vertex for _, vertex in sorted(local_vertices[key])
        ]
        node_rows = [vertices[(bias, vertex)] for vertex in ids]
        points = [
            (float(row["x_um"]), float(row["y_um"]))
            for row in node_rows
        ]
        for quantity, field, x_column, y_column in vector_specs:
            raw_gradient = p1_gradient(
                points, [float(row[field]) for row in node_rows]
            )
            candidate = (
                -raw_gradient[0] * 1.0e4,
                -raw_gradient[1] * 1.0e4,
            )
            reference = (
                float(element[x_column]),
                float(element[y_column]),
            )
            reference_mag = math.hypot(*reference)
            absolute_residual = math.hypot(
                candidate[0] - reference[0],
                candidate[1] - reference[1],
            )
            vector_rows.append(
                {
                    "case": case_name,
                    "bias_V": bias,
                    "element": element_id,
                    "element_class": (
                        "contact"
                        if geometry_by_element[element_id][
                            "contact_adjacent"
                        ]
                        else "interior"
                    ),
                    "quantity": quantity,
                    "sign_convention": (
                        "-grad(psi)" if quantity == "electric_field"
                        else "-grad(qfp)"
                    ),
                    "vela_x_V_cm": candidate[0],
                    "vela_y_V_cm": candidate[1],
                    "sentaurus_x_V_cm": reference[0],
                    "sentaurus_y_V_cm": reference[1],
                    "absolute_vector_residual_V_cm": absolute_residual,
                    "relative_vector_residual": (
                        absolute_residual / reference_mag
                        if reference_mag > 0.0 else ""
                    ),
                    "angle_error_deg": (
                        ""
                        if vector_angle_deg(candidate, reference) is None
                        else vector_angle_deg(candidate, reference)
                    ),
                    "active_vector": int(reference_mag >= 100.0),
                }
            )

        field_by_carrier = {
            "electron": math.hypot(
                float(element["grad_qf_n_x_V_cm"]),
                float(element["grad_qf_n_y_V_cm"]),
            ),
            "hole": math.hypot(
                float(element["grad_qf_p_x_V_cm"]),
                float(element["grad_qf_p_y_V_cm"]),
            ),
        }
        electric_field = math.hypot(
            float(element["efield_x_V_cm"]),
            float(element["efield_y_V_cm"]),
        )
        for carrier in ("electron", "hole"):
            native_mobility = float(
                element[
                    "mu_n_cm2_Vs"
                    if carrier == "electron" else "mu_p_cm2_Vs"
                ]
            )
            low_values: list[float] = []
            high_values: list[float] = []
            for vertex in ids:
                x_um = float(vertices[(bias, vertex)]["x_um"])
                net_doping, _ = node_doping(x_um, junction_um)
                low = masetti(net_doping, carrier)
                low_values.append(low)
                high_values.append(
                    field_limit(low, field_by_carrier[carrier], carrier)
                )
            cell_low = statistics.mean(low_values)
            cell_high = statistics.mean(high_values)
            mobility_rows.append(
                {
                    "case": case_name,
                    "bias_V": bias,
                    "element": element_id,
                    "local_edge": "",
                    "carrier": carrier,
                    "support": "cell_vertex_average",
                    "driving_force": "native_element_qfp_gradient",
                    "field_V_cm": field_by_carrier[carrier],
                    "vela_low_field_mobility_cm2_Vs": cell_low,
                    "vela_final_mobility_cm2_Vs": cell_high,
                    "sentaurus_element_mobility_cm2_Vs": native_mobility,
                    "low_field_absolute_error_dex": abs_dex(
                        cell_low, native_mobility
                    ),
                    "final_absolute_error_dex": abs_dex(
                        cell_high, native_mobility
                    ),
                }
            )
            qfp_column = (
                "eQFP_V" if carrier == "electron" else "hQFP_V"
            )
            for edge in edge_rows[key]:
                node0 = int(edge["start"])
                node1 = int(edge["end"])
                qfp0 = float(vertices[(bias, node0)][qfp_column])
                qfp1 = float(vertices[(bias, node1)][qfp_column])
                edge_field = abs(qfp1 - qfp0) / float(edge["length_um"]) * 1.0e4
                endpoint_low: list[float] = []
                endpoint_high: list[float] = []
                for vertex in (node0, node1):
                    x_um = float(vertices[(bias, vertex)]["x_um"])
                    net_doping, _ = node_doping(x_um, junction_um)
                    low = masetti(net_doping, carrier)
                    endpoint_low.append(low)
                    endpoint_high.append(
                        field_limit(low, edge_field, carrier)
                    )
                edge_low = statistics.mean(endpoint_low)
                edge_high = statistics.mean(endpoint_high)
                mobility_rows.append(
                    {
                        "case": case_name,
                        "bias_V": bias,
                        "element": element_id,
                        "local_edge": int(edge["local_edge"]),
                        "carrier": carrier,
                        "support": "element_local_edge_endpoint_average",
                        "driving_force": "endpoint_qfp_difference",
                        "field_V_cm": edge_field,
                        "vela_low_field_mobility_cm2_Vs": edge_low,
                        "vela_final_mobility_cm2_Vs": edge_high,
                        "sentaurus_element_mobility_cm2_Vs": native_mobility,
                        "low_field_absolute_error_dex": abs_dex(
                            edge_low, native_mobility
                        ),
                        "final_absolute_error_dex": abs_dex(
                            edge_high, native_mobility
                        ),
                    }
                )

            forced = {
                "electric_field": electric_field,
                "quasi_fermi_gradient": field_by_carrier[carrier],
            }
            native_vertex_alpha = statistics.mean(
                float(
                    row[
                        "alpha_n_cm_inv"
                        if carrier == "electron" else "alpha_p_cm_inv"
                    ]
                )
                for row in node_rows
            )
            for driver, field in forced.items():
                candidate_alpha = alpha(field, carrier)
                alpha_rows.append(
                    {
                        "case": case_name,
                        "bias_V": bias,
                        "element": element_id,
                        "element_class": (
                            "contact"
                            if geometry_by_element[element_id][
                                "contact_adjacent"
                            ]
                            else "interior"
                        ),
                        "carrier": carrier,
                        "forced_driver": driver,
                        "field_V_cm": field,
                        "vela_alpha_cm_inv": candidate_alpha,
                        "sentaurus_vertex_mean_alpha_cm_inv": (
                            native_vertex_alpha
                        ),
                        "support_mismatched_alpha_error_dex": abs_dex(
                            candidate_alpha, native_vertex_alpha
                        ),
                        "observation_label": (
                            "unsupported_native_element_alpha"
                        ),
                    }
                )

    write_csv(output_root / "node_density_comparison.csv", density_rows)
    write_csv(output_root / "element_vector_comparison.csv", vector_rows)
    write_csv(output_root / "mobility_support_comparison.csv", mobility_rows)
    write_csv(output_root / "element_alpha_replay.csv", alpha_rows)

    density_summary = {
        carrier: error_summary(
            [
                row[f"{carrier}_absolute_error_dex"]
                for row in density_rows
            ]
        )
        for carrier in ("electron", "hole")
    }
    vector_summary: dict[str, Any] = {}
    for quantity, *_ in vector_specs:
        selected = [
            row for row in vector_rows if row["quantity"] == quantity
        ]
        active = [row for row in selected if row["active_vector"]]
        vector_summary[quantity] = {
            "all_relative": error_summary(
                [
                    None
                    if row["relative_vector_residual"] == ""
                    else float(row["relative_vector_residual"])
                    for row in selected
                ]
            ),
            "active_relative": error_summary(
                [
                    None
                    if row["relative_vector_residual"] == ""
                    else float(row["relative_vector_residual"])
                    for row in active
                ]
            ),
            "active_angle_deg": error_summary(
                [
                    None
                    if row["angle_error_deg"] == ""
                    else float(row["angle_error_deg"])
                    for row in active
                ]
            ),
        }
    mobility_summary: dict[str, Any] = {}
    for support in (
        "cell_vertex_average",
        "element_local_edge_endpoint_average",
    ):
        mobility_summary[support] = {
            carrier: error_summary(
                [
                    row["final_absolute_error_dex"]
                    for row in mobility_rows
                    if row["support"] == support
                    and row["carrier"] == carrier
                ]
            )
            for carrier in ("electron", "hole")
        }
    summary = {
        "schema": "pn2d_general_tri3_imported_state/v1",
        "status": "valid",
        "case_name": case_name,
        "source_schema": SCHEMA_ID,
        "source_manifest_sha256": digest(raw_manifest_path),
        "source_log_sha256": digest(log_path),
        "exact_biases_V": list(EXACT_BIASES_V),
        "junction_um": junction_um,
        "density": density_summary,
        "vectors": vector_summary,
        "mobility": mobility_summary,
        "alpha_native_element_observation": (
            "insufficient_native_observation"
        ),
        "outputs": {
            name: digest(output_root / name)
            for name in (
                "node_density_comparison.csv",
                "element_vector_comparison.csv",
                "mobility_support_comparison.csv",
                "element_alpha_replay.csv",
            )
        },
    }
    (output_root / "analysis_manifest.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="ascii",
        newline="\n",
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
