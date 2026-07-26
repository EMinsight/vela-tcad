#!/usr/bin/env python3
"""Replay Vela element-local SG currents and GSS/Laux cell vectors."""

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
from scripts.diagnose_pn2d_general_tri3_imported_state import (
    VT_300K_V,
    abs_dex,
    carrier_density,
    effective_ni,
    error_summary,
    field_limit,
    masetti,
    node_doping,
    percentile,
    vector_angle_deg,
)
from scripts.pn2d_general_tri3_contract import (
    EXACT_BIASES_V,
    SCHEMA_ID,
    SENTAURUS_RELEASE,
)


Q_SENT_C = 1.6021918e-19


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


def bernoulli(value: float) -> float:
    magnitude = abs(value)
    if magnitude < 1.0e-8:
        return 1.0 - 0.5 * value + value * value / 12.0
    if value > 50.0:
        return value * math.exp(-value)
    if value < -50.0:
        return -value
    return value / math.expm1(value)


def variable_ni_flux(
    ni0: float,
    ni1: float,
    psi0: float,
    psi1: float,
    qfp0: float,
    qfp1: float,
    mobility_cm2_vs: float,
    length_cm: float,
    carrier: str,
) -> float:
    if qfp0 == qfp1:
        return 0.0
    coefficient = mobility_cm2_vs * VT_300K_V / length_cm
    if carrier == "electron":
        eta = (
            (psi1 - psi0) / VT_300K_V + math.log(ni1 / ni0)
        )
        n0 = carrier_density(ni0, psi0, qfp0, carrier)
        n1 = carrier_density(ni1, psi1, qfp1, carrier)
        return coefficient * (
            bernoulli(-eta) * n0 - bernoulli(eta) * n1
        )
    eta = (psi1 - psi0) / VT_300K_V + math.log(ni0 / ni1)
    p0 = carrier_density(ni0, psi0, qfp0, carrier)
    p1 = carrier_density(ni1, psi1, qfp1, carrier)
    return coefficient * (
        bernoulli(eta) * p0 - bernoulli(-eta) * p1
    )


def projection_pair(
    tangent_a: tuple[float, float],
    value_a: float,
    tangent_b: tuple[float, float],
    value_b: float,
) -> tuple[float, float]:
    determinant = (
        tangent_a[0] * tangent_b[1]
        - tangent_a[1] * tangent_b[0]
    )
    if abs(determinant) <= 1.0e-14:
        raise ValueError("parallel Tri3 edges")
    return (
        (value_a * tangent_b[1] - tangent_a[1] * value_b)
        / determinant,
        (tangent_a[0] * value_b - value_a * tangent_b[0])
        / determinant,
    )


def geometric_partial_volumes(
    points: list[tuple[float, float]],
) -> list[float]:
    partials: list[float] = []
    for local_edge in range(3):
        p0 = points[local_edge]
        p1 = points[(local_edge + 1) % 3]
        opposite = points[(local_edge + 2) % 3]
        u = (p0[0] - opposite[0], p0[1] - opposite[1])
        v = (p1[0] - opposite[0], p1[1] - opposite[1])
        cross = abs(u[0] * v[1] - u[1] * v[0])
        cotangent = (
            (u[0] * v[0] + u[1] * v[1]) / cross
            if cross > 0.0 else 0.0
        )
        length2 = (
            (p1[0] - p0[0]) ** 2 + (p1[1] - p0[1]) ** 2
        )
        partials.append(max(0.0, 0.25 * cotangent * length2))
    return partials


def gss_laux_vector(
    points: list[tuple[float, float]],
    signed_edge_current: list[float],
    partials: list[float],
) -> tuple[float, float]:
    tangents: list[tuple[float, float]] = []
    lengths: list[float] = []
    for local_edge in range(3):
        p0 = points[local_edge]
        p1 = points[(local_edge + 1) % 3]
        delta = (p1[0] - p0[0], p1[1] - p0[1])
        length = math.hypot(*delta)
        lengths.append(length)
        tangents.append((delta[0] / length, delta[1] / length))
    pair_vectors = (
        projection_pair(
            tangents[0], signed_edge_current[0],
            tangents[1], signed_edge_current[1],
        ),
        projection_pair(
            tangents[0], signed_edge_current[0],
            tangents[2], signed_edge_current[2],
        ),
        projection_pair(
            tangents[1], signed_edge_current[1],
            tangents[2], signed_edge_current[2],
        ),
    )

    def pair_index(first: int, second: int) -> int:
        low, high = sorted((first, second))
        return high - 1 if low == 0 else 2

    edge_vectors: list[tuple[float, float]] = []
    for target in range(3):
        others = [value for value in range(3) if value != target]
        first, second = others
        first_weight = 2.0 * partials[first] / lengths[first]
        second_weight = 2.0 * partials[second] / lengths[second]
        first_pair = pair_vectors[pair_index(target, first)]
        second_pair = pair_vectors[pair_index(target, second)]
        weight_sum = first_weight + second_weight
        if weight_sum > 0.0:
            edge_vectors.append(
                (
                    (
                        first_weight * first_pair[0]
                        + second_weight * second_pair[0]
                    ) / weight_sum,
                    (
                        first_weight * first_pair[1]
                        + second_weight * second_pair[1]
                    ) / weight_sum,
                )
            )
        else:
            edge_vectors.append(
                (
                    0.5 * (first_pair[0] + second_pair[0]),
                    0.5 * (first_pair[1] + second_pair[1]),
                )
            )
    total = sum(partials)
    if total <= 0.0:
        raise ValueError("Tri3 has no positive partial volume")
    return (
        sum(
            partials[index] * edge_vectors[index][0]
            for index in range(3)
        ) / total,
        sum(
            partials[index] * edge_vectors[index][1]
            for index in range(3)
        ) / total,
    )


def main() -> int:
    args = parse_args()
    raw_root = args.raw_root.resolve()
    output_root = args.output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    manifest_path = raw_root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="ascii"))
    if manifest.get("schema") != SCHEMA_ID:
        raise ValueError("raw schema mismatch")
    if manifest.get("status") != "passed":
        raise ValueError("raw root is incomplete")
    if manifest.get("sentaurus_release") != SENTAURUS_RELEASE:
        raise ValueError("Sentaurus release mismatch")
    case_names = tuple(manifest["cases"])
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
    for row in groups["measures"]:
        local_vertices[
            (float(row["bias_V"]), int(row["element"]))
        ].append((int(row["local_vertex"]), int(row["vertex"])))
    raw_edges: dict[tuple[float, int], list[dict[str, Any]]] = defaultdict(list)
    for row in groups["edges"]:
        raw_edges[
            (float(row["bias_V"]), int(row["element"]))
        ].append(row)
    elements = {
        (float(row["bias_V"]), int(row["element"])): row
        for row in groups["elements"]
    }
    used_x = [
        float(vertices[(bias, vertex)]["x_um"])
        for (bias, _), values in local_vertices.items()
        for _, vertex in values
    ]
    junction_um = 0.5 * (min(used_x) + max(used_x))

    edge_output: list[dict[str, Any]] = []
    vector_output: list[dict[str, Any]] = []
    for key, element in sorted(elements.items()):
        bias, element_id = key
        ids = [vertex for _, vertex in sorted(local_vertices[key])]
        node_rows = [vertices[(bias, vertex)] for vertex in ids]
        points = [
            (float(row["x_um"]), float(row["y_um"]))
            for row in node_rows
        ]
        partials = geometric_partial_volumes(points)
        sentaurus_signed: dict[str, list[float]] = {
            "electron": [],
            "hole": [],
        }
        vela_signed: dict[str, list[float]] = {
            "electron": [],
            "hole": [],
        }
        for local_edge in range(3):
            node0 = ids[local_edge]
            node1 = ids[(local_edge + 1) % 3]
            raw = next(
                row for row in raw_edges[key]
                if {
                    int(row["start"]), int(row["end"])
                } == {node0, node1}
            )
            orientation = (
                1
                if (int(raw["start"]), int(raw["end"]))
                == (node0, node1)
                else -1
            )
            length_um = math.dist(
                points[local_edge], points[(local_edge + 1) % 3]
            )
            length_cm = length_um * 1.0e-4
            endpoint_rows = (
                vertices[(bias, node0)],
                vertices[(bias, node1)],
            )
            intrinsic: list[float] = []
            net_doping: list[float] = []
            for row in endpoint_rows:
                net, total = node_doping(
                    float(row["x_um"]), junction_um
                )
                net_doping.append(net)
                intrinsic.append(effective_ni(total)[0])
            for carrier in ("electron", "hole"):
                qfp_column = (
                    "eQFP_V" if carrier == "electron" else "hQFP_V"
                )
                field = abs(
                    float(endpoint_rows[1][qfp_column])
                    - float(endpoint_rows[0][qfp_column])
                ) / length_cm
                mobility = statistics.mean(
                    field_limit(
                        masetti(net, carrier), field, carrier
                    )
                    for net in net_doping
                )
                flux = variable_ni_flux(
                    intrinsic[0],
                    intrinsic[1],
                    float(endpoint_rows[0]["psi_V"]),
                    float(endpoint_rows[1]["psi_V"]),
                    float(endpoint_rows[0][qfp_column]),
                    float(endpoint_rows[1][qfp_column]),
                    mobility,
                    length_cm,
                    carrier,
                )
                # The production hole routine returns the continuity-residual
                # flux. The Sentaurus Tcl reconstruction uses conventional
                # hole current, which has the opposite sign. Its electron
                # reconstruction uses the production continuity convention.
                vela_j = (
                    Q_SENT_C * flux
                    if carrier == "electron"
                    else -Q_SENT_C * flux
                )
                raw_column = (
                    "sg_jn_A_cm2"
                    if carrier == "electron" else "sg_jp_A_cm2"
                )
                sentaurus_j = orientation * float(raw[raw_column])
                vela_signed[carrier].append(vela_j)
                sentaurus_signed[carrier].append(sentaurus_j)
                edge_output.append(
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
                        "angle_class": geometry_by_element[element_id][
                            "angle_class"
                        ],
                        "carrier": carrier,
                        "local_edge": local_edge,
                        "node0": node0,
                        "node1": node1,
                        "orientation_vs_probe": orientation,
                        "length_um": length_um,
                        # ReadCoefficient is an unoriented scalar support. Only
                        # the reconstructed edge current changes sign when the
                        # probe edge is reversed.
                        "read_coefficient": float(raw["kappa"]),
                        "vela_truncated_partial_volume_um2": partials[local_edge],
                        "vela_mobility_cm2_Vs": mobility,
                        "sentaurus_element_mobility_cm2_Vs": float(
                            element[
                                "mu_n_cm2_Vs"
                                if carrier == "electron"
                                else "mu_p_cm2_Vs"
                            ]
                        ),
                        "vela_sg_current_A_cm2": vela_j,
                        "sentaurus_box_operator_sg_current_A_cm2": (
                            sentaurus_j
                        ),
                        "absolute_error_dex": abs_dex(
                            abs(vela_j), abs(sentaurus_j)
                        ),
                        "sign_match": int(
                            vela_j == 0.0
                            or sentaurus_j == 0.0
                            or math.copysign(1.0, vela_j)
                            == math.copysign(1.0, sentaurus_j)
                        ),
                        "observation_label": "box_operator_reconstruction",
                    }
                )

        for carrier in ("electron", "hole"):
            vela_vector = gss_laux_vector(
                points, vela_signed[carrier], partials
            )
            sentaurus_box_vector = gss_laux_vector(
                points, sentaurus_signed[carrier], partials
            )
            native_vector = (
                float(
                    element[
                        "current_n_x_A_cm2"
                        if carrier == "electron"
                        else "current_p_x_A_cm2"
                    ]
                ),
                float(
                    element[
                        "current_n_y_A_cm2"
                        if carrier == "electron"
                        else "current_p_y_A_cm2"
                    ]
                ),
            )
            native_gradient = (
                float(
                    element[
                        "grad_qf_n_x_V_cm"
                        if carrier == "electron"
                        else "grad_qf_p_x_V_cm"
                    ]
                ),
                float(
                    element[
                        "grad_qf_n_y_V_cm"
                        if carrier == "electron"
                        else "grad_qf_p_y_V_cm"
                    ]
                ),
            )
            for candidate, vector, label in (
                (
                    "vela_recomputed",
                    vela_vector,
                    "vela_recomputation",
                ),
                (
                    "sentaurus_box_operator",
                    sentaurus_box_vector,
                    "box_operator_reconstruction",
                ),
            ):
                candidate_mag = math.hypot(*vector)
                native_mag = math.hypot(*native_vector)
                vector_output.append(
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
                        "angle_class": geometry_by_element[element_id][
                            "angle_class"
                        ],
                        "carrier": carrier,
                        "candidate": candidate,
                        "candidate_x_A_cm2": vector[0],
                        "candidate_y_A_cm2": vector[1],
                        "candidate_magnitude_A_cm2": candidate_mag,
                        "native_x_A_cm2": native_vector[0],
                        "native_y_A_cm2": native_vector[1],
                        "native_magnitude_A_cm2": native_mag,
                        "magnitude_absolute_error_dex": abs_dex(
                            candidate_mag, native_mag
                        ),
                        "vector_angle_error_deg": (
                            ""
                            if vector_angle_deg(
                                vector, native_vector
                            ) is None
                            else vector_angle_deg(vector, native_vector)
                        ),
                        "current_vs_native_qfp_gradient_angle_deg": (
                            ""
                            if vector_angle_deg(
                                vector, native_gradient
                            ) is None
                            else vector_angle_deg(vector, native_gradient)
                        ),
                        "observation_label": label,
                    }
                )

    write_csv(output_root / "element_local_edges.csv", edge_output)
    write_csv(output_root / "cell_current_vectors.csv", vector_output)
    edge_summary = {
        carrier: {
            "error_dex": error_summary(
                [
                    row["absolute_error_dex"]
                    for row in edge_output
                    if row["carrier"] == carrier
                ]
            ),
            "nonzero_sign_fraction": (
                sum(
                    row["sign_match"]
                    for row in edge_output
                    if row["carrier"] == carrier
                )
                / sum(
                    1 for row in edge_output
                    if row["carrier"] == carrier
                )
            ),
        }
        for carrier in ("electron", "hole")
    }
    vector_summary: dict[str, Any] = {}
    for candidate in ("vela_recomputed", "sentaurus_box_operator"):
        vector_summary[candidate] = {
            carrier: error_summary(
                [
                    row["magnitude_absolute_error_dex"]
                    for row in vector_output
                    if row["candidate"] == candidate
                    and row["carrier"] == carrier
                ]
            )
            for carrier in ("electron", "hole")
        }
    summary = {
        "schema": "pn2d_general_tri3_element_edge_current/v1",
        "status": "valid",
        "case_name": case_name,
        "source_manifest_sha256": digest(manifest_path),
        "source_log_sha256": digest(log_path),
        "exact_biases_V": list(EXACT_BIASES_V),
        "edge_current": edge_summary,
        "cell_current_vector": vector_summary,
        "terminal_current": "pending_same_support_aggregation",
        "internal_kcl": "pending_same_support_aggregation",
        "outputs": {
            name: digest(output_root / name)
            for name in (
                "element_local_edges.csv",
                "cell_current_vectors.csv",
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
