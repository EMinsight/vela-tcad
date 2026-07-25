#!/usr/bin/env python3
"""Replay Minimal6 element avalanche generation from Sentaurus edge SG data."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from collections import defaultdict
from pathlib import Path
from statistics import median
from typing import Any, Iterable


Q_LEGACY_C = 1.6021918e-19
UM2_TO_CM2 = 1.0e-8
AMPERE_PER_CM_TO_AMPERE_PER_UM = 1.0e-4
SOURCE_INTEGRAL_TO_AMPERE_PER_UM = (
    Q_LEGACY_C * UM2_TO_CM2 * AMPERE_PER_CM_TO_AMPERE_PER_UM
)
TARGET_BIASES = (-1.0, -10.0, -20.0)
TOPOLOGIES = ("mirror", "sketch")
CARRIERS = ("electron", "hole")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-root", type=Path, required=True)
    parser.add_argument("--vela-factorization", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_tokens(line: str) -> dict[str, str]:
    return dict(re.findall(r"(\w+)=([^\s]+)", line))


def typed_row(tokens: dict[str, str], integer_keys: set[str]) -> dict[str, Any]:
    row: dict[str, Any] = {}
    for key, value in tokens.items():
        if key in integer_keys:
            row[key] = int(value)
        else:
            row[key] = float(value)
    return row


def parse_log(path: Path) -> dict[str, list[dict[str, Any]]]:
    groups: dict[str, list[dict[str, Any]]] = defaultdict(list)
    prefixes = {
        "AVAL_PROBE_VERTEX ": ("vertices", {"bias_V", "vertex"}),
        "AVAL_PROBE_ELEMENT ": ("elements", {"bias_V", "element"}),
        "AVAL_PROBE_MEASURE ": (
            "measures",
            {"bias_V", "element", "local_vertex", "vertex"},
        ),
        "AVAL_PROBE_EDGE ": (
            "edges",
            {
                "bias_V",
                "element",
                "local_edge",
                "edge",
                "start",
                "end",
            },
        ),
        "AVAL_PROBE_INTEGRAL ": ("integrals", {"bias_V"}),
    }
    text = path.read_text(encoding="ascii", errors="strict")
    for line in text.splitlines():
        for prefix, (group, integer_keys) in prefixes.items():
            if line.startswith(prefix):
                groups[group].append(
                    typed_row(parse_tokens(line[len(prefix) :]), integer_keys)
                )
                break
    expected = {
        "vertices": 3 * 10,
        "elements": 3 * 4,
        "measures": 3 * 12,
        "edges": 3 * 12,
        "integrals": 3,
    }
    for group, count in expected.items():
        if len(groups[group]) != count:
            raise ValueError(
                f"{path}: expected {count} {group}, got {len(groups[group])}"
            )
    return groups


def parse_plt(path: Path) -> tuple[list[str], list[dict[str, float]]]:
    text = path.read_text(encoding="ascii", errors="strict")
    info, data = text.split("Data {", 1)
    match = re.search(r"datasets\s*=\s*\[(.*?)\]", info, re.S)
    if match is None:
        raise ValueError(f"{path}: missing datasets")
    names = re.findall(r'"([^"]+)"', match.group(1))
    values = [
        float(token)
        for token in re.findall(
            r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?", data
        )
    ]
    if len(values) % len(names):
        raise ValueError(f"{path}: non-rectangular CurrentPlot")
    rows = [
        dict(zip(names, values[index : index + len(names)], strict=True))
        for index in range(0, len(values), len(names))
    ]
    return names, rows


def currentplot_targets(path: Path) -> list[dict[str, float]]:
    names, rows = parse_plt(path)
    voltage_name = next(
        name
        for name in names
        if name.endswith("Anode OuterVoltage") or name == "Anode OuterVoltage"
    )
    result = []
    for bias in TARGET_BIASES:
        match = min(rows, key=lambda row: abs(row[voltage_name] - bias))
        if abs(match[voltage_name] - bias) > 1.0e-8:
            raise ValueError(f"{path}: missing CurrentPlot bias {bias:g}")
        result.append({"bias_V": bias, **match})
    return result


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty CSV {path}")
    fields = list(rows[0])
    with path.open("w", newline="", encoding="ascii") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def solve_pair(
    tangent_a: tuple[float, float],
    value_a: float,
    tangent_b: tuple[float, float],
    value_b: float,
) -> tuple[float, float]:
    ax, ay = tangent_a
    bx, by = tangent_b
    determinant = ax * by - ay * bx
    if abs(determinant) <= 1.0e-14:
        raise ValueError("parallel edge tangents cannot reconstruct a vector")
    return (
        (value_a * by - ay * value_b) / determinant,
        (ax * value_b - value_a * bx) / determinant,
    )


def least_squares(
    rows: Iterable[tuple[tuple[float, float], float, float]]
) -> tuple[float, float]:
    a11 = a12 = a22 = b1 = b2 = 0.0
    for (tx, ty), value, weight in rows:
        a11 += weight * tx * tx
        a12 += weight * tx * ty
        a22 += weight * ty * ty
        b1 += weight * tx * value
        b2 += weight * ty * value
    determinant = a11 * a22 - a12 * a12
    if abs(determinant) <= 1.0e-30:
        raise ValueError("singular weighted edge reconstruction")
    return (
        (b1 * a22 - b2 * a12) / determinant,
        (a11 * b2 - a12 * b1) / determinant,
    )


def gss_laux_vector(
    edges: list[dict[str, Any]], current_key: str
) -> tuple[float, float]:
    pair_vectors: dict[tuple[int, int], tuple[float, float]] = {}
    for first in range(3):
        for second in range(first + 1, 3):
            pair_vectors[(first, second)] = solve_pair(
                (edges[first]["tangent_x"], edges[first]["tangent_y"]),
                edges[first][current_key],
                (edges[second]["tangent_x"], edges[second]["tangent_y"]),
                edges[second][current_key],
            )

    edge_vectors: list[tuple[tuple[float, float], float]] = []
    for target in range(3):
        others = [index for index in range(3) if index != target]
        first, second = others
        first_weight = edges[first]["kappa"] * edges[first]["length_um"]
        second_weight = edges[second]["kappa"] * edges[second]["length_um"]
        first_pair = pair_vectors[tuple(sorted((target, first)))]
        second_pair = pair_vectors[tuple(sorted((target, second)))]
        weight_sum = first_weight + second_weight
        if weight_sum <= 0.0:
            vector = (
                0.5 * (first_pair[0] + second_pair[0]),
                0.5 * (first_pair[1] + second_pair[1]),
            )
        else:
            vector = (
                (first_weight * first_pair[0] + second_weight * second_pair[0])
                / weight_sum,
                (first_weight * first_pair[1] + second_weight * second_pair[1])
                / weight_sum,
            )
        partial_area = (
            0.5
            * edges[target]["kappa"]
            * edges[target]["length_um"]
            * edges[target]["length_um"]
        )
        edge_vectors.append((vector, partial_area))
    total_area = sum(weight for _, weight in edge_vectors)
    if total_area <= 0.0:
        raise ValueError("triangle has no positive box partial area")
    return (
        sum(vector[0] * weight for vector, weight in edge_vectors) / total_area,
        sum(vector[1] * weight for vector, weight in edge_vectors) / total_area,
    )


def charon_whitney_vector(
    edges: list[dict[str, Any]],
    local_vertices: list[int],
    vertex_by_id: dict[int, dict[str, Any]],
    current_key: str,
) -> tuple[float, float]:
    if len(local_vertices) != 3:
        raise ValueError("triangle must have three local vertices")
    points = [
        (
            vertex_by_id[vertex]["x_um"] * 1.0e-4,
            vertex_by_id[vertex]["y_um"] * 1.0e-4,
        )
        for vertex in local_vertices
    ]
    reference_edges = (
        (local_vertices[0], local_vertices[1]),
        (local_vertices[1], local_vertices[2]),
        (local_vertices[2], local_vertices[0]),
    )
    edge_by_nodes = {
        tuple(sorted((int(edge["start"]), int(edge["end"])))): edge
        for edge in edges
    }
    dofs = []
    for start, end in reference_edges:
        edge = edge_by_nodes[tuple(sorted((start, end)))]
        orientation = (
            1.0
            if (int(edge["start"]), int(edge["end"])) == (start, end)
            else -1.0
        )
        length_cm = edge["length_um"] * 1.0e-4
        dofs.append(orientation * edge[current_key] * length_cm)

    # Lowest-order Whitney edge basis at the reference-triangle centroid.
    basis = ((2.0 / 3.0, 1.0 / 3.0), (-1.0 / 3.0, 1.0 / 3.0),
             (-1.0 / 3.0, -2.0 / 3.0))
    ref_x = sum(dof * item[0] for dof, item in zip(dofs, basis, strict=True))
    ref_y = sum(dof * item[1] for dof, item in zip(dofs, basis, strict=True))
    ax = points[1][0] - points[0][0]
    ay = points[1][1] - points[0][1]
    bx = points[2][0] - points[0][0]
    by = points[2][1] - points[0][1]
    determinant = ax * by - ay * bx
    if abs(determinant) <= 1.0e-30:
        raise ValueError("degenerate physical triangle")
    # Covariant Piola transform: physical = A^{-T} reference.
    return (
        (by * ref_x - ay * ref_y) / determinant,
        (-bx * ref_x + ax * ref_y) / determinant,
    )


def vector_magnitude(vector: tuple[float, float]) -> float:
    return math.hypot(vector[0], vector[1])


def abs_dex(value: float, reference: float) -> float | None:
    if value <= 0.0 or reference <= 0.0:
        return None
    return abs(math.log10(value / reference))


def format_value(value: Any) -> Any:
    if value is None:
        return ""
    if isinstance(value, float):
        return format(value, ".17g")
    return value


def candidate_vectors(
    edges: list[dict[str, Any]],
    local_vertices: list[int],
    vertex_by_id: dict[int, dict[str, Any]],
    element: dict[str, Any],
    current_key: str,
    native_prefix: str,
) -> dict[str, tuple[float, float]]:
    all_rows = [
        (
            (edge["tangent_x"], edge["tangent_y"]),
            edge[current_key],
            1.0,
        )
        for edge in edges
    ]
    active_rows = [
        (
            (edge["tangent_x"], edge["tangent_y"]),
            edge[current_key],
            max(edge["kappa"] * edge["length_um"], 0.0),
        )
        for edge in edges
        if edge["kappa"] > 0.0
    ]
    return {
        "gss_laux_edge_volume_weighted": gss_laux_vector(edges, current_key),
        "charon_whitney_hcurl_cell_average": charon_whitney_vector(
            edges, local_vertices, vertex_by_id, current_key
        ),
        "genius_least_squares_tangent": least_squares(all_rows),
        "box_active_edge_exact": least_squares(active_rows),
        "native_element_vector_control": (
            element[f"{native_prefix}_x_A_cm2"],
            element[f"{native_prefix}_y_A_cm2"],
        ),
    }


def find_integral_name(names: Iterable[str], token: str) -> str:
    matches = [name for name in names if token in name]
    if len(matches) != 1:
        raise ValueError(f"expected one CurrentPlot name containing {token!r}")
    return matches[0]


def load_vela_state_sources(path: Path) -> dict[tuple[str, float], float]:
    result = {}
    with path.open(newline="", encoding="ascii") as handle:
        for row in csv.DictReader(handle):
            bias = float(row["bias_V"])
            if row["topology"] in TOPOLOGIES and bias in TARGET_BIASES:
                result[(row["topology"], bias)] = float(
                    row["vela_candidate_source_per_cm_s"]
                )
    if len(result) != len(TOPOLOGIES) * len(TARGET_BIASES):
        raise ValueError("Vela factorization lacks selected topology/bias states")
    return result


def run(raw_root: Path, vela_factorization: Path, output: Path) -> dict[str, Any]:
    raw_root = raw_root.resolve()
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    vela_sources = load_vela_state_sources(vela_factorization)

    all_vertices: list[dict[str, Any]] = []
    all_elements: list[dict[str, Any]] = []
    all_edges: list[dict[str, Any]] = []
    all_measures: list[dict[str, Any]] = []
    all_integrals: list[dict[str, Any]] = []
    current_rows: list[dict[str, Any]] = []
    reconstructions: list[dict[str, Any]] = []
    predicted_nodes: list[dict[str, Any]] = []
    state_summary: list[dict[str, Any]] = []
    raw_hashes: dict[str, str] = {}

    for topology in TOPOLOGIES:
        log_path = raw_root / topology / "default" / "run_default.out"
        plt_path = (
            raw_root
            / topology
            / "default"
            / "runtime_element_avalanche_probe_default.plt"
        )
        groups = parse_log(log_path)
        plt_targets = currentplot_targets(plt_path)
        names, _ = parse_plt(plt_path)
        e_integral_name = find_integral_name(names, "eAvalancheIntegral")
        h_integral_name = find_integral_name(names, "hAvalancheIntegral")
        total_integral_name = find_integral_name(
            [
                name
                for name in names
                if "eAvalancheIntegral" not in name
                and "hAvalancheIntegral" not in name
            ],
            "AvalancheIntegral",
        )
        for path in (log_path, plt_path):
            raw_hashes[path.relative_to(raw_root).as_posix()] = sha256(path)

        for group_name, target in (
            ("vertices", all_vertices),
            ("elements", all_elements),
            ("edges", all_edges),
            ("measures", all_measures),
            ("integrals", all_integrals),
        ):
            for row in groups[group_name]:
                target.append({"topology": topology, **row})

        current_by_bias = {}
        for row in plt_targets:
            normalized = {
                "topology": topology,
                "bias_V": row["bias_V"],
                "e_avalanche_integral_um2_cm3_s": row[e_integral_name],
                "h_avalanche_integral_um2_cm3_s": row[h_integral_name],
                "total_avalanche_integral_um2_cm3_s": row[
                    total_integral_name
                ],
            }
            for carrier, label in (
                ("electron", "eCurrent"),
                ("hole", "hCurrent"),
                ("total", "TotalCurrent"),
            ):
                for contact in ("Anode", "Cathode"):
                    key = next(
                        name
                        for name in names
                        if name.endswith(f"{contact} {label}")
                        or name == f"{contact} {label}"
                    )
                    normalized[f"{contact}_{carrier}_A_um"] = row[key]
            current_rows.append(normalized)
            current_by_bias[row["bias_V"]] = normalized

        for bias in TARGET_BIASES:
            vertices = [
                row for row in groups["vertices"] if row["bias_V"] == int(bias)
            ]
            physical_vertices = [row for row in vertices if row["vertex"] < 6]
            vertex_by_id = {row["vertex"]: row for row in physical_vertices}
            elements = [
                row for row in groups["elements"] if row["bias_V"] == int(bias)
            ]
            measures = [
                row for row in groups["measures"] if row["bias_V"] == int(bias)
            ]
            edges = [
                row for row in groups["edges"] if row["bias_V"] == int(bias)
            ]
            measure_by_element = defaultdict(list)
            for row in measures:
                measure_by_element[row["element"]].append(row)
            edge_by_element = defaultdict(list)
            for row in edges:
                edge_by_element[row["element"]].append(row)

            node_predictions: dict[
                tuple[str, str, int], dict[str, float]
            ] = defaultdict(lambda: {"weighted_source": 0.0, "measure": 0.0})
            candidate_integrals: dict[tuple[str, str], float] = defaultdict(float)

            for element in elements:
                element_id = element["element"]
                element_edges = sorted(
                    edge_by_element[element_id],
                    key=lambda row: row["local_edge"],
                )
                element_measures = sorted(
                    measure_by_element[element_id],
                    key=lambda row: row["local_vertex"],
                )
                local_vertices = [row["vertex"] for row in element_measures]
                if len(element_edges) != 3 or len(local_vertices) != 3:
                    raise ValueError("incomplete triangle topology in runtime log")
                best_measure = max(
                    element_measures, key=lambda row: row["measure_um2"]
                )
                best_vertex = vertex_by_id[best_measure["vertex"]]
                for carrier in CARRIERS:
                    if carrier == "electron":
                        current_key = "sg_jn_A_cm2"
                        native_prefix = "current_n"
                        alpha = best_vertex["alpha_n_cm_inv"]
                        native_generation_key = "generation_n_cm3_s"
                    else:
                        current_key = "sg_jp_A_cm2"
                        native_prefix = "current_p"
                        alpha = best_vertex["alpha_p_cm_inv"]
                        native_generation_key = "generation_p_cm3_s"
                    vectors = candidate_vectors(
                        element_edges,
                        local_vertices,
                        vertex_by_id,
                        element,
                        current_key,
                        native_prefix,
                    )
                    for candidate, vector in vectors.items():
                        magnitude = vector_magnitude(vector)
                        generation = alpha * magnitude / Q_LEGACY_C
                        row = {
                            "topology": topology,
                            "bias_V": bias,
                            "element": element_id,
                            "carrier": carrier,
                            "candidate": candidate,
                            "vector_x_A_cm2": vector[0],
                            "vector_y_A_cm2": vector[1],
                            "magnitude_A_cm2": magnitude,
                            "alpha_cm_inv": alpha,
                            "generation_cm3_s": generation,
                            "native_element_vector_magnitude_A_cm2": math.hypot(
                                element[f"{native_prefix}_x_A_cm2"],
                                element[f"{native_prefix}_y_A_cm2"],
                            ),
                        }
                        reconstructions.append(row)
                        for measure_row in element_measures:
                            vertex_id = measure_row["vertex"]
                            weight = measure_row["measure_um2"]
                            bucket = node_predictions[
                                (carrier, candidate, vertex_id)
                            ]
                            bucket["weighted_source"] += generation * weight
                            bucket["measure"] += weight
                            candidate_integrals[(carrier, candidate)] += (
                                generation * weight
                            )

                for carrier, native_key in (
                    ("electron", "generation_n_cm3_s"),
                    ("hole", "generation_p_cm3_s"),
                ):
                    for measure_row in element_measures:
                        vertex = vertex_by_id[measure_row["vertex"]]
                        if native_key not in vertex:
                            raise ValueError(native_generation_key)

            candidate_names = sorted(
                {
                    candidate
                    for _, candidate, _ in node_predictions
                }
            )
            native_integrals = {
                "electron": current_by_bias[bias][
                    "e_avalanche_integral_um2_cm3_s"
                ],
                "hole": current_by_bias[bias][
                    "h_avalanche_integral_um2_cm3_s"
                ],
            }
            for carrier in CARRIERS:
                native_key = (
                    "generation_n_cm3_s"
                    if carrier == "electron"
                    else "generation_p_cm3_s"
                )
                for candidate in candidate_names:
                    node_errors = []
                    for vertex_id, vertex in sorted(vertex_by_id.items()):
                        bucket = node_predictions[
                            (carrier, candidate, vertex_id)
                        ]
                        predicted = (
                            bucket["weighted_source"] / bucket["measure"]
                        )
                        native = vertex[native_key]
                        error = abs_dex(predicted, native)
                        if error is not None:
                            node_errors.append(error)
                        predicted_nodes.append(
                            {
                                "topology": topology,
                                "bias_V": bias,
                                "vertex": vertex_id,
                                "carrier": carrier,
                                "candidate": candidate,
                                "predicted_generation_cm3_s": predicted,
                                "native_generation_cm3_s": native,
                                "absolute_error_dex": error,
                            }
                        )
                    predicted_integral = candidate_integrals[
                        (carrier, candidate)
                    ]
                    native_integral = native_integrals[carrier]
                    state_summary.append(
                        {
                            "topology": topology,
                            "bias_V": bias,
                            "carrier": carrier,
                            "candidate": candidate,
                            "predicted_integral_um2_cm3_s": predicted_integral,
                            "native_integral_um2_cm3_s": native_integral,
                            "integral_absolute_error_dex": abs_dex(
                                predicted_integral, native_integral
                            ),
                            "node_median_absolute_error_dex": (
                                median(node_errors) if node_errors else None
                            ),
                            "node_max_absolute_error_dex": (
                                max(node_errors) if node_errors else None
                            ),
                        }
                    )

            native_total = current_by_bias[bias][
                "total_avalanche_integral_um2_cm3_s"
            ]
            vela_source = vela_sources[(topology, bias)]
            state_summary.append(
                {
                    "topology": topology,
                    "bias_V": bias,
                    "carrier": "total",
                    "candidate": "vela_triangle_proxy_existing",
                    "predicted_integral_um2_cm3_s": vela_source,
                    "native_integral_um2_cm3_s": native_total,
                    "integral_absolute_error_dex": abs_dex(
                        vela_source, native_total
                    ),
                    "node_median_absolute_error_dex": None,
                    "node_max_absolute_error_dex": None,
                }
            )

    for collection in (
        all_vertices,
        all_elements,
        all_edges,
        all_measures,
        all_integrals,
        current_rows,
        reconstructions,
        predicted_nodes,
        state_summary,
    ):
        for row in collection:
            for key, value in list(row.items()):
                row[key] = format_value(value)

    outputs = {
        "vertices.csv": all_vertices,
        "elements.csv": all_elements,
        "element_edges.csv": all_edges,
        "element_vertex_measures.csv": all_measures,
        "runtime_integrals.csv": all_integrals,
        "currentplot_targets.csv": current_rows,
        "element_reconstructions.csv": reconstructions,
        "node_generation_replay.csv": predicted_nodes,
        "state_source_summary.csv": state_summary,
    }
    output_hashes = {}
    for name, rows in outputs.items():
        path = output / name
        write_csv(path, rows)
        output_hashes[name] = sha256(path)

    summary_rows = [
        row
        for row in state_summary
        if row["candidate"] != "native_element_vector_control"
    ]
    rankings = []
    candidates = sorted(
        {
            row["candidate"]
            for row in summary_rows
            if row["candidate"] != "vela_triangle_proxy_existing"
        }
    )
    for candidate in candidates:
        values = [
            float(row["integral_absolute_error_dex"])
            for row in summary_rows
            if row["candidate"] == candidate
            and row["integral_absolute_error_dex"] != ""
        ]
        rankings.append(
            {
                "candidate": candidate,
                "sample_count": len(values),
                "median_integral_absolute_error_dex": median(values),
                "max_integral_absolute_error_dex": max(values),
            }
        )
    rankings.sort(key=lambda row: row["median_integral_absolute_error_dex"])

    manifest = {
        "schema_version": 1,
        "status": "valid_diagnostic_replay",
        "experiment": "pn2d_minimal6_element_avalanche_replay",
        "scope": {
            "topologies": list(TOPOLOGIES),
            "biases_V": list(TARGET_BIASES),
            "production_formula_changed": False,
            "native_edge_current_observed": False,
            "element_edge_sg_status": "documented_operator_reconstruction",
            "native_element_generation_status": (
                "vertex_output_redistributed_with_ReadMeasure"
            ),
        },
        "raw_sha256": raw_hashes,
        "output_sha256": output_hashes,
        "rankings": rankings,
    }
    manifest_path = output / "manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="ascii",
    )
    return manifest


def main() -> int:
    args = parse_args()
    manifest = run(args.raw_root, args.vela_factorization, args.output)
    print(json.dumps(manifest["rankings"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
