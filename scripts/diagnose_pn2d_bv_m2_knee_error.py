#!/usr/bin/env python3
"""Read-only localization of PN2D M2 BV knee-region discrepancies."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path
from typing import Any


Q_C = 1.602176634e-19
PARTICLE_FLUX_TO_PER_CM2_S = 1.0e6
LINE_SOURCE_TO_A_PER_UM = 1.0e-12
BIAS_TOLERANCE_V = 1.0e-8
SOURCE_ACTIVE_FRACTION = 1.0e-3


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vela-manifest", type=Path, required=True)
    parser.add_argument("--vela-probe", type=Path, required=True)
    parser.add_argument("--vela-run-a-root", type=Path, required=True)
    parser.add_argument("--vela-run-b-root", type=Path, required=True)
    parser.add_argument("--sentaurus-manifest", type=Path, required=True)
    parser.add_argument("--parity", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def magnitude(values: list[float]) -> float:
    return math.sqrt(sum(float(value) ** 2 for value in values))


def log_ratio(candidate: float, reference: float) -> float:
    if candidate <= 0.0 or reference <= 0.0:
        raise ValueError("log-ratio inputs must be positive")
    return math.log10(candidate / reference)


def weighted_mean(values: list[tuple[float, float]]) -> float:
    total = sum(weight for _, weight in values)
    if total <= 0.0:
        raise ValueError("weighted mean has no positive support")
    return sum(value * weight for value, weight in values) / total


def weighted_ratio_log(values: list[tuple[float, float]]) -> float:
    """Log of a source-weighted arithmetic ratio for integral counterfactuals."""
    total = sum(weight for _, weight in values)
    if total <= 0.0:
        raise ValueError("weighted ratio has no positive support")
    ratio = sum(value * weight for value, weight in values) / total
    return math.log10(ratio)


def percentile(values: list[float], fraction: float) -> float:
    if not values:
        raise ValueError("percentile has no values")
    ordered = sorted(values)
    index = min(len(ordered) - 1, int(fraction * (len(ordered) - 1)))
    return ordered[index]


def pearson(left: list[float], right: list[float]) -> float:
    if len(left) != len(right) or len(left) < 2:
        raise ValueError("Pearson correlation requires paired values")
    mean_left = statistics.mean(left)
    mean_right = statistics.mean(right)
    numerator = sum(
        (a - mean_left) * (b - mean_right) for a, b in zip(left, right)
    )
    denominator = math.sqrt(
        sum((a - mean_left) ** 2 for a in left)
        * sum((b - mean_right) ** 2 for b in right)
    )
    return numerator / denominator if denominator else 0.0


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty CSV {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def nearest_bias(value: float, biases: list[float]) -> float | None:
    closest = min(biases, key=lambda bias: abs(value - bias))
    return closest if abs(value - closest) <= BIAS_TOLERANCE_V else None


def knee_biases(sentaurus: dict[str, Any], vela: dict[str, Any]) -> list[float]:
    def observed(manifest: dict[str, Any]) -> set[float]:
        return {
            float(row["requested_bias_V"])
            for row in manifest["aggregate_records"]
            if row["branch"] == "avalanche_on"
            and row["quantity"] == "terminal_current"
            and row["carrier"] == "total"
        }

    common = observed(sentaurus) & observed(vela)
    selected = sorted(
        (bias for bias in common if -20.0 - BIAS_TOLERANCE_V <= bias <= -18.0 + BIAS_TOLERANCE_V),
        reverse=True,
    )
    if len(selected) < 3 or selected[0] != -18.0 or selected[-1] != -20.0:
        raise ValueError(f"incomplete -18 V to -20 V knee lattice: {selected}")
    return selected


def aggregate_index(
    manifest: dict[str, Any], biases: list[float]
) -> dict[tuple[float, str, str, str], float]:
    result: dict[tuple[float, str, str, str], float] = {}
    for row in manifest["aggregate_records"]:
        if row["branch"] != "avalanche_on":
            continue
        bias = nearest_bias(float(row["requested_bias_V"]), biases)
        if bias is None:
            continue
        result[(bias, row["quantity"], row["carrier"], row["provenance"])] = abs(
            float(row["value"])
        )
    return result


def field_indexes(
    manifest: dict[str, Any], biases: list[float]
) -> tuple[
    dict[tuple[float, str, str, str, str], tuple[float, tuple[int, ...]]],
    dict[int, tuple[int, ...]],
]:
    result: dict[
        tuple[float, str, str, str, str], tuple[float, tuple[int, ...]]
    ] = {}
    cells: dict[int, tuple[int, ...]] = {}
    selected_quantities = {
        "current_density",
        "quasi_fermi_gradient",
        "mobility",
        "avalanche_alpha",
        "avalanche_generation",
        "integrated_source",
        "density",
    }
    for row in manifest["field_records"]:
        if row["branch"] != "avalanche_on" or row["quantity"] not in selected_quantities:
            continue
        bias = nearest_bias(float(row["requested_bias_V"]), biases)
        if bias is None:
            continue
        support = str(row["support_key"])
        connectivity = tuple(int(value) for value in row.get("connectivity", []))
        result[
            (
                bias,
                row["quantity"],
                row["carrier"],
                row["provenance"],
                support,
            )
        ] = (magnitude(row["values"]), connectivity)
        if (
            row["quantity"] == "current_density"
            and row["support_kind"] == "cell"
            and connectivity
        ):
            cells[int(support.split(":")[1])] = tuple(sorted(connectivity))
    return result, cells


def load_probe(
    path: Path, biases: list[float]
) -> tuple[
    dict[tuple[float, int, str], tuple[float, float, float]],
    dict[tuple[float, int, str], list[dict[str, Any]]],
    dict[tuple[float, int, str], float],
    dict[tuple[float, int, int, str], float],
]:
    cell: dict[tuple[float, int, str], tuple[float, float, float]] = {}
    edges: dict[tuple[float, int, str], list[dict[str, Any]]] = defaultdict(list)
    alpha: dict[tuple[float, int, str], float] = {}
    measures: dict[tuple[float, int, int, str], float] = {}
    with path.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            bias = nearest_bias(float(row["bias_V"]), biases)
            if bias is None:
                continue
            cell_id = int(row["cell_id"])
            carrier = row["carrier"]
            if row["support_kind"] == "element_edge_gss_laux":
                edges[(bias, cell_id, carrier)].append(
                    {
                        "connectivity": tuple(
                            sorted((int(row["node0"]), int(row["node1"])))
                        ),
                        "current_A_cm2": abs(float(row["selected_flux_magnitude"]))
                        * Q_C
                        * PARTICLE_FLUX_TO_PER_CM2_S,
                        "mobility_cm2_Vs": float(row["final_mobility"]),
                        "drive_V_cm": abs(float(row["high_field_drive"])),
                    }
                )
            elif row["support_kind"] == "element_vertex_gss_laux":
                current = math.hypot(
                    float(row["current_vector_x"]), float(row["current_vector_y"])
                ) * Q_C * PARTICLE_FLUX_TO_PER_CM2_S
                drive = math.hypot(
                    float(row["qf_gradient_x"]), float(row["qf_gradient_y"])
                )
                cell.setdefault(
                    (bias, cell_id, carrier),
                    (current, drive, abs(float(row["impact_field"]))),
                )
                alpha[(bias, cell_id, carrier)] = abs(float(row["alpha"]))
                measures[(bias, cell_id, int(row["node0"]), carrier)] = float(
                    row["source_measure"]
                )
    return cell, edges, alpha, measures


def overlap(left: dict[str, float], right: dict[str, float]) -> float:
    left_total = sum(left.values())
    right_total = sum(right.values())
    if left_total <= 0.0 or right_total <= 0.0:
        raise ValueError("source-map overlap requires positive totals")
    return sum(
        min(left.get(key, 0.0) / left_total, right.get(key, 0.0) / right_total)
        for key in set(left) | set(right)
    )


def verify_determinism(root_a: Path, root_b: Path) -> dict[str, Any]:
    names = (
        "avalanche_on/iv.csv",
        "avalanche_on/process_probe.csv",
        "avalanche_on/newton_attempts.csv",
        "avalanche_on/newton_history.csv",
    )
    files: list[dict[str, Any]] = []
    for name in names:
        path_a = root_a / name
        path_b = root_b / name
        hash_a = sha256(path_a)
        hash_b = sha256(path_b)
        files.append(
            {"path": name, "run_a_sha256": hash_a, "run_b_sha256": hash_b, "equal": hash_a == hash_b}
        )
    return {"all_equal": all(row["equal"] for row in files), "files": files}


def build_analysis(args: argparse.Namespace) -> tuple[
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[dict[str, Any]],
    list[dict[str, Any]],
    dict[str, Any],
]:
    sentaurus = load_json(args.sentaurus_manifest)
    vela = load_json(args.vela_manifest)
    parity = load_json(args.parity)
    if sentaurus.get("simulator") != "sentaurus" or vela.get("simulator") != "vela":
        raise ValueError("simulator manifests are reversed or invalid")
    biases = knee_biases(sentaurus, vela)
    sent_aggregate = aggregate_index(sentaurus, biases)
    vela_aggregate = aggregate_index(vela, biases)
    sent_fields, cells = field_indexes(sentaurus, biases)
    vela_fields, _ = field_indexes(vela, biases)
    vela_cell, vela_edges, vela_alpha, vela_measures = load_probe(args.vela_probe, biases)

    sent_edge: dict[tuple[float, tuple[int, ...], str, str], float] = {}
    for (bias, quantity, carrier, provenance, support), (value, connectivity) in sent_fields.items():
        if quantity == "current_density" and "/local_edge:" in support:
            sent_edge[(bias, tuple(sorted(connectivity)), carrier, provenance)] = value

    bias_rows: list[dict[str, Any]] = []
    carrier_rows: list[dict[str, Any]] = []
    edge_rows: list[dict[str, Any]] = []
    mapping_rows: list[dict[str, Any]] = []

    for bias in biases:
        sent_source_maps = {
            carrier: {
                support: value
                for (row_bias, quantity, row_carrier, provenance, support), (value, _) in sent_fields.items()
                if row_bias == bias
                and quantity == "integrated_source"
                and row_carrier == carrier
                and provenance == "operator_replay"
            }
            for carrier in ("electron", "hole", "total")
        }
        vela_source_maps = {
            carrier: {
                support: value
                for (row_bias, quantity, row_carrier, provenance, support), (value, _) in vela_fields.items()
                if row_bias == bias
                and quantity == "integrated_source"
                and row_carrier == carrier
                and provenance == "solver_used"
            }
            for carrier in ("electron", "hole", "total")
        }
        sent_source_total = sum(sent_source_maps["total"].values())
        vela_source_total = sum(vela_source_maps["total"].values())
        sent_cells: dict[int, float] = defaultdict(float)
        vela_cells: dict[int, float] = defaultdict(float)
        for support, value in sent_source_maps["total"].items():
            sent_cells[int(support.split("/")[0].split(":")[1])] += value
        for support, value in vela_source_maps["total"].items():
            vela_cells[int(support.split("/")[0].split(":")[1])] += value
        active_threshold = max(sent_cells.values()) * SOURCE_ACTIVE_FRACTION
        active_cells = {
            cell_id: value for cell_id, value in sent_cells.items() if value >= active_threshold
        }

        stage_values: dict[str, float] = {}
        for carrier in ("electron", "hole"):
            carrier_sent_source = sum(sent_source_maps[carrier].values())
            carrier_vela_source = sum(vela_source_maps[carrier].values())
            carrier_rows.append(
                {
                    "bias_V": bias,
                    "reverse_bias_magnitude_V": abs(bias),
                    "carrier": carrier,
                    "sentaurus_source_A_um": carrier_sent_source,
                    "vela_source_A_um": carrier_vela_source,
                    "log10_vela_over_sentaurus_dex": log_ratio(
                        carrier_vela_source, carrier_sent_source
                    ),
                    "sentaurus_fraction": carrier_sent_source / sent_source_total,
                    "vela_fraction": carrier_vela_source / vela_source_total,
                }
            )
            for stage in ("current", "drive", "mobility"):
                comparisons: list[tuple[float, float]] = []
                for cell_id, weight in active_cells.items():
                    if stage == "current":
                        candidate = vela_cell[(bias, cell_id, carrier)][0]
                        reference = sent_fields[
                            (bias, "current_density", carrier, "native", f"cell:{cell_id}")
                        ][0]
                    elif stage == "drive":
                        candidate = vela_cell[(bias, cell_id, carrier)][1]
                        reference = sent_fields[
                            (
                                bias,
                                "quasi_fermi_gradient",
                                carrier,
                                "native",
                                f"cell:{cell_id}",
                            )
                        ][0]
                    else:
                        dominant = max(
                            vela_edges[(bias, cell_id, carrier)],
                            key=lambda row: row["current_A_cm2"],
                        )
                        candidate = dominant["mobility_cm2_Vs"]
                        reference = sent_fields[
                            (bias, "mobility", carrier, "native", f"cell:{cell_id}")
                        ][0]
                    if candidate > 0.0 and reference > 0.0:
                        comparisons.append((log_ratio(candidate, reference), weight))
                stage_values[f"{carrier}_{stage}_log_ratio_dex"] = weighted_mean(comparisons)

            alpha_comparisons: list[tuple[float, float]] = []
            density_comparisons: list[tuple[float, float]] = []
            for support, weight in sent_source_maps[carrier].items():
                if weight <= 0.0:
                    continue
                sent_record = sent_fields[
                    (bias, "integrated_source", carrier, "operator_replay", support)
                ]
                node = int(sent_record[1][0])
                cell_id = int(support.split("/")[0].split(":")[1])
                sent_alpha = sent_fields[
                    (bias, "avalanche_alpha", carrier, "native", f"node:{node}")
                ][0]
                candidate_alpha = vela_alpha[(bias, cell_id, carrier)]
                if sent_alpha > 0.0 and candidate_alpha > 0.0:
                    alpha_comparisons.append(
                        (candidate_alpha / sent_alpha, weight)
                    )
                sent_density = sent_fields[
                    (bias, "density", carrier, "native", f"node:{node}")
                ][0]
                candidate_density_record = vela_fields.get(
                    (bias, "density", carrier, "solver_used", f"node:{node}")
                )
                if candidate_density_record and sent_density > 0.0 and candidate_density_record[0] > 0.0:
                    density_comparisons.append(
                        (log_ratio(candidate_density_record[0], sent_density), weight)
                    )
            stage_values[f"{carrier}_alpha_log_ratio_dex"] = weighted_ratio_log(
                alpha_comparisons
            )
            stage_values[f"{carrier}_density_log_ratio_dex"] = weighted_mean(
                density_comparisons
            )

        native_errors: list[float] = []
        replay_raw_errors: list[float] = []
        replay_scaled_errors: list[float] = []
        for cell_id in active_cells:
            for carrier in ("electron", "hole"):
                for edge in vela_edges[(bias, cell_id, carrier)]:
                    candidate = edge["current_A_cm2"]
                    native = sent_edge.get(
                        (bias, edge["connectivity"], carrier, "reconstructed"), 0.0
                    )
                    replay = sent_edge.get(
                        (bias, edge["connectivity"], carrier, "operator_replay"), 0.0
                    )
                    if candidate > 0.0 and native > 0.0:
                        native_errors.append(abs(log_ratio(candidate, native)))
                    if candidate > 0.0 and replay > 0.0:
                        replay_raw_errors.append(abs(log_ratio(candidate, replay)))
                        replay_scaled_errors.append(
                            abs(log_ratio(candidate, replay * PARTICLE_FLUX_TO_PER_CM2_S))
                        )
        edge_rows.append(
            {
                "bias_V": bias,
                "reverse_bias_magnitude_V": abs(bias),
                "active_cell_count": len(active_cells),
                "native_projection_median_abs_error_dex": statistics.median(native_errors),
                "native_projection_p95_abs_error_dex": percentile(native_errors, 0.95),
                "operator_replay_raw_median_abs_error_dex": statistics.median(
                    replay_raw_errors
                ),
                "operator_replay_x1e6_median_abs_error_dex": statistics.median(
                    replay_scaled_errors
                ),
                "operator_replay_scale_inference": "manifest values behave as A/um^2 despite A/cm^2 label",
            }
        )

        geometry_errors: list[float] = []
        for carrier in ("electron", "hole"):
            for support, source in sent_source_maps[carrier].items():
                if source <= 0.0:
                    continue
                source_record = sent_fields[
                    (bias, "integrated_source", carrier, "operator_replay", support)
                ]
                node = int(source_record[1][0])
                cell_id = int(support.split("/")[0].split(":")[1])
                generation = sent_fields[
                    (
                        bias,
                        "avalanche_generation",
                        carrier,
                        "native",
                        f"node:{node}",
                    )
                ][0]
                if generation <= 0.0:
                    continue
                sent_measure = source / (Q_C * generation * LINE_SOURCE_TO_A_PER_UM)
                vela_measure = vela_measures[(bias, cell_id, node, carrier)]
                geometry_errors.append(abs(vela_measure - sent_measure) / sent_measure)
        vertex_overlap = overlap(sent_source_maps["total"], vela_source_maps["total"])
        cell_overlap = overlap(
            {f"cell:{key}": value for key, value in sent_cells.items()},
            {f"cell:{key}": value for key, value in vela_cells.items()},
        )
        mapping_rows.append(
            {
                "bias_V": bias,
                "reverse_bias_magnitude_V": abs(bias),
                "vertex_source_overlap": vertex_overlap,
                "cell_source_overlap": cell_overlap,
                "same_hotspot_vertex": int(
                    max(sent_source_maps["total"], key=sent_source_maps["total"].get)
                    == max(vela_source_maps["total"], key=vela_source_maps["total"].get)
                ),
                "same_hotspot_cell": int(
                    max(sent_cells, key=sent_cells.get)
                    == max(vela_cells, key=vela_cells.get)
                ),
                "source_measure_max_relative_error": max(geometry_errors),
                "source_measure_median_relative_error": statistics.median(geometry_errors),
            }
        )

        sent_current = sent_aggregate[(bias, "terminal_current", "total", "native")]
        vela_current = vela_aggregate[(bias, "terminal_current", "total", "solver_used")]
        bias_rows.append(
            {
                "bias_V": bias,
                "reverse_bias_magnitude_V": abs(bias),
                "sentaurus_terminal_current_A_um": sent_current,
                "vela_terminal_current_A_um": vela_current,
                "terminal_abs_log_error_dex": abs(log_ratio(vela_current, sent_current)),
                "sentaurus_integrated_source_A_um": sent_source_total,
                "vela_integrated_source_A_um": vela_source_total,
                "source_signed_log_ratio_dex": log_ratio(
                    vela_source_total, sent_source_total
                ),
                "source_abs_log_error_dex": abs(
                    log_ratio(vela_source_total, sent_source_total)
                ),
                "source_terminal_error_difference_dex": abs(
                    log_ratio(vela_source_total, sent_source_total)
                )
                - abs(log_ratio(vela_current, sent_current)),
                "active_cell_count": len(active_cells),
                **stage_values,
                "vertex_source_overlap": vertex_overlap,
                "cell_source_overlap": cell_overlap,
                "source_measure_max_relative_error": max(geometry_errors),
            }
        )

    terminal_errors = [row["terminal_abs_log_error_dex"] for row in bias_rows]
    source_errors = [row["source_abs_log_error_dex"] for row in bias_rows]
    electron_density_errors = [
        abs(row["electron_density_log_ratio_dex"]) for row in bias_rows
    ]
    hole_density_errors = [abs(row["hole_density_log_ratio_dex"]) for row in bias_rows]
    determinism = verify_determinism(args.vela_run_a_root, args.vela_run_b_root)
    first = bias_rows[0]
    last = bias_rows[-1]
    carrier_fraction_max_difference = max(
        abs(row["vela_fraction"] - row["sentaurus_fraction"])
        for row in carrier_rows
    )
    summary = {
        "schema": "vela.pn2d_bv_m2_knee_readonly_diagnostic.v1",
        "status": "passed" if determinism["all_equal"] else "failed",
        "outcome": "self_consistent_carrier_current_amplitude_localized",
        "observation_only": True,
        "acceptance_thresholds_modified": False,
        "biases_V": biases,
        "input_artifacts": {
            "vela_manifest": {"path": str(args.vela_manifest.resolve()), "sha256": sha256(args.vela_manifest)},
            "vela_probe": {"path": str(args.vela_probe.resolve()), "sha256": sha256(args.vela_probe)},
            "sentaurus_manifest": {"path": str(args.sentaurus_manifest.resolve()), "sha256": sha256(args.sentaurus_manifest)},
            "parity": {"path": str(args.parity.resolve()), "sha256": sha256(args.parity)},
        },
        "determinism": determinism,
        "correlations": {
            "terminal_vs_integrated_source_error": pearson(terminal_errors, source_errors),
            "terminal_vs_electron_density_error": pearson(
                terminal_errors, electron_density_errors
            ),
            "terminal_vs_hole_density_error": pearson(
                terminal_errors, hole_density_errors
            ),
        },
        "knee_error_growth": {
            "terminal_error_dex": last["terminal_abs_log_error_dex"]
            - first["terminal_abs_log_error_dex"],
            "integrated_source_error_dex": last["source_abs_log_error_dex"]
            - first["source_abs_log_error_dex"],
            "electron_cell_current_deficit_dex": abs(last["electron_current_log_ratio_dex"])
            - abs(first["electron_current_log_ratio_dex"]),
            "hole_cell_current_deficit_dex": abs(last["hole_current_log_ratio_dex"])
            - abs(first["hole_current_log_ratio_dex"]),
            "electron_density_deficit_dex": abs(last["electron_density_log_ratio_dex"])
            - abs(first["electron_density_log_ratio_dex"]),
            "hole_density_deficit_dex": abs(last["hole_density_log_ratio_dex"])
            - abs(first["hole_density_log_ratio_dex"]),
        },
        "controls": {
            "maximum_carrier_source_fraction_difference": carrier_fraction_max_difference,
            "maximum_qfp_drive_abs_log_ratio_dex": max(
                abs(row[f"{carrier}_drive_log_ratio_dex"])
                for row in bias_rows
                for carrier in ("electron", "hole")
            ),
            "maximum_dominant_edge_mobility_abs_log_ratio_dex": max(
                abs(row[f"{carrier}_mobility_log_ratio_dex"])
                for row in bias_rows
                for carrier in ("electron", "hole")
            ),
            "maximum_source_weighted_alpha_abs_log_ratio_dex": max(
                abs(row[f"{carrier}_alpha_log_ratio_dex"])
                for row in bias_rows
                for carrier in ("electron", "hole")
            ),
            "maximum_source_measure_relative_error": max(
                row["source_measure_max_relative_error"] for row in mapping_rows
            ),
            "minimum_vertex_source_overlap": min(
                row["vertex_source_overlap"] for row in mapping_rows
            ),
            "minimum_cell_source_overlap": min(
                row["cell_source_overlap"] for row in mapping_rows
            ),
            "all_hotspot_vertices_equal": all(
                row["same_hotspot_vertex"] for row in mapping_rows
            ),
            "all_hotspot_cells_equal": all(row["same_hotspot_cell"] for row in mapping_rows),
        },
        "causal_interpretation": {
            "localized_stage": "self-consistent carrier density and SG/Laux current amplitude",
            "excluded_as_growing_primary_driver": [
                "ionization coefficient",
                "quasi-Fermi-gradient drive",
                "dominant-edge mobility",
                "carrier partition",
                "source geometry measure",
                "hotspot relocation",
            ],
            "remaining_ambiguity": (
                "Self-consistent snapshots cannot order whether the density deficit or "
                "the native-current versus SG/Laux support difference initiates feedback."
            ),
            "next_read_only_control": (
                "Replay the M2 Sentaurus states at -18, -19.5, -19.7, and -20 V "
                "through Vela postprocess-only SG/Laux on the exact M2 mesh."
            ),
        },
        "parity_contract_status": parity.get("status"),
        "parity_contract_outcome": parity.get("outcome"),
    }
    return bias_rows, carrier_rows, edge_rows, mapping_rows, summary


def main() -> int:
    args = parse_args()
    args.output_root.mkdir(parents=True, exist_ok=True)
    bias_rows, carrier_rows, edge_rows, mapping_rows, summary = build_analysis(args)
    write_csv(args.output_root / "bias_summary.csv", bias_rows)
    write_csv(args.output_root / "carrier_source_components.csv", carrier_rows)
    write_csv(args.output_root / "edge_current_summary.csv", edge_rows)
    write_csv(args.output_root / "source_mapping_summary.csv", mapping_rows)
    (args.output_root / "diagnostic.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True, ensure_ascii=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0 if summary["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
