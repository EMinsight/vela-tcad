#!/usr/bin/env python3
"""Score the Task-7 frozen SG-vector candidate against two controls."""

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


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def write_rows(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise RuntimeError(f"refusing to write empty CSV: {path}")
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def percentile(values: list[float], fraction: float) -> float:
    if not values:
        return 0.0
    ordered = sorted(values)
    position = fraction * (len(ordered) - 1)
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def active_dex_metrics(
    pairs: list[tuple[float, float]],
    floor_fraction: float = 1.0e-12,
) -> dict[str, float | int]:
    peak = max((abs(reference) for _, reference in pairs), default=0.0)
    floor = max(peak * floor_fraction, 1.0e-300)
    errors = [
        abs(math.log10(abs(candidate) / abs(reference)))
        for candidate, reference in pairs
        if abs(candidate) > 0.0 and abs(reference) > floor
    ]
    return {
        "active_count": len(errors),
        "median_abs_error_dex": statistics.median(errors) if errors else 0.0,
        "p95_abs_error_dex": percentile(errors, 0.95),
        "maximum_abs_error_dex": max(errors, default=0.0),
    }


def field(path: Path) -> dict[int, tuple[float, ...]]:
    values: dict[int, tuple[float, ...]] = {}
    for row in read_rows(path):
        components = tuple(
            float(row[name])
            for name in ("component0", "component1", "component2")
            if name in row and row[name] != ""
        )
        values[int(row["node_id"])] = components
    return values


def vector_norm(values: tuple[float, ...]) -> float:
    return math.sqrt(sum(value * value for value in values))


def source_pairs_from_triangle(
    external_rows: list[dict[str, str]],
    candidate_column: str,
) -> dict[float, list[tuple[float, float]]]:
    pairs: dict[float, list[tuple[float, float]]] = defaultdict(list)
    for row in external_rows:
        pairs[float(row["bias_V"])].append((
            float(row[candidate_column]),
            float(row["sentaurus_generation_vela_geometry_source_integral_per_m_s"]),
        ))
    return pairs


def candidate_metrics(
    root: Path,
    result: dict[str, Any],
    active_source_floor_fraction: float,
) -> tuple[
    dict[float, dict[str, dict[str, float | int]]],
    dict[float, dict[str, float | int]],
    dict[float, dict[str, float]],
]:
    node_sources: dict[
        tuple[float, int], dict[str, float]
    ] = defaultdict(lambda: {"candidate": 0.0, "measure": 0.0})
    node_vectors: dict[
        tuple[float, str, int], dict[str, float]
    ] = defaultdict(
        lambda: {"weighted_x": 0.0, "weighted_y": 0.0, "measure": 0.0}
    )
    sent_generation_by_bias: dict[float, dict[int, tuple[float, ...]]] = {}
    sent_current_by_bias: dict[
        tuple[float, str], dict[int, tuple[float, ...]]
    ] = {}

    for case in result["cases"]:
        bias = float(case["bias_V"])
        sentaurus_case = Path(case["sentaurus_case"])
        sent_generation = field(
            sentaurus_case / "fields" / "ImpactIonization_region0.csv"
        )
        sent_generation_by_bias[bias] = sent_generation
        sent_current = {
            "electron": field(
                sentaurus_case / "fields" / "eCurrentDensity_region0.csv"
            ),
            "hole": field(
                sentaurus_case / "fields" / "hCurrentDensity_region0.csv"
            ),
        }
        for carrier, values in sent_current.items():
            sent_current_by_bias[(bias, carrier)] = values
        process_rows = read_rows(Path(case["outputs"]["sg_vector_process"]["path"]))
        seen_vectors: set[tuple[str, int, int]] = set()
        cell_node_measures: dict[tuple[int, int], float] = {}
        for row in process_rows:
            if row["support_kind"] != "element_vertex_gss_laux":
                continue
            carrier = row["carrier"]
            cell = int(row["cell_id"])
            node = int(row["node0"])
            vector_key = (carrier, cell, node)
            if vector_key in seen_vectors:
                raise RuntimeError(f"duplicate SG vector record: {vector_key}")
            seen_vectors.add(vector_key)
            measure = float(row["source_measure_m2"])
            measure_key = (cell, node)
            previous_measure = cell_node_measures.setdefault(measure_key, measure)
            if previous_measure != measure:
                raise RuntimeError(
                    f"inconsistent carrier measure for cell/node {measure_key}"
                )
            node_sources[(bias, node)]["candidate"] += float(
                row["source_integral_per_m_s"]
            )
            vela_sign = -1.0 if carrier == "electron" else 1.0
            vector = node_vectors[(bias, carrier, node)]
            vector["weighted_x"] += (
                vela_sign * float(row["current_vector_x_per_m2_s"]) * measure
            )
            vector["weighted_y"] += (
                vela_sign * float(row["current_vector_y_per_m2_s"]) * measure
            )
            vector["measure"] += measure
        for (_, node), measure in cell_node_measures.items():
            node_sources[(bias, node)]["measure"] += measure

    for (bias, node), group in node_sources.items():
        group["reference"] = (
            sent_generation_by_bias[bias][node][0] * 1.0e6 * group["measure"]
        )

    active_keys: set[tuple[float, int]] = set()
    for bias in {key[0] for key in node_sources}:
        peak = max(
            group["reference"]
            for key, group in node_sources.items()
            if key[0] == bias
        )
        floor = peak * active_source_floor_fraction
        active_keys.update(
            key for key, group in node_sources.items()
            if key[0] == bias and group["reference"] >= floor
        )

    current_by_bias_carrier: dict[
        tuple[float, str], list[tuple[float, float]]
    ] = defaultdict(list)
    angle_by_bias_carrier: dict[tuple[float, str], list[float]] = defaultdict(list)
    sign_by_bias_carrier: dict[tuple[float, str], list[bool]] = defaultdict(list)
    for (bias, carrier, node), vector in node_vectors.items():
        if (bias, node) not in active_keys or vector["measure"] <= 0.0:
            continue
        vela = (
            vector["weighted_x"] / vector["measure"],
            vector["weighted_y"] / vector["measure"],
        )
        sent = sent_current_by_bias[(bias, carrier)][node][:2]
        vela_norm = vector_norm(vela)
        sent_norm = vector_norm(sent)
        current_by_bias_carrier[(bias, carrier)].append((
            Q_C * vela_norm,
            1.0e4 * sent_norm,
        ))
        if vela_norm > 0.0 and sent_norm > 0.0:
            cosine = max(-1.0, min(
                1.0,
                sum(a * b for a, b in zip(vela, sent))
                / (vela_norm * sent_norm),
            ))
            angle_by_bias_carrier[(bias, carrier)].append(
                math.degrees(math.acos(cosine))
            )
            sign_by_bias_carrier[(bias, carrier)].append(
                cosine > 0.0
            )

    current_metrics: dict[
        float, dict[str, dict[str, float | int]]
    ] = defaultdict(dict)
    direction_metrics: dict[float, dict[str, float]] = defaultdict(dict)
    for (bias, carrier), pairs in current_by_bias_carrier.items():
        current_metrics[bias][carrier] = active_dex_metrics(pairs)
        angles = angle_by_bias_carrier[(bias, carrier)]
        signs = sign_by_bias_carrier[(bias, carrier)]
        direction_metrics[bias][f"{carrier}_median_angle_deg"] = (
            statistics.median(angles) if angles else 0.0
        )
        direction_metrics[bias][f"{carrier}_p95_angle_deg"] = percentile(
            angles, 0.95
        )
        direction_metrics[bias][f"{carrier}_sign_agreement_fraction"] = (
            sum(signs) / len(signs) if signs else 1.0
        )

    source_pairs: dict[float, list[tuple[float, float]]] = defaultdict(list)
    for (bias, _), group in node_sources.items():
        source_pairs[bias].append((group["candidate"], group["reference"]))
    source_metrics = {
        bias: active_dex_metrics(pairs, active_source_floor_fraction)
        for bias, pairs in source_pairs.items()
    }
    return current_metrics, source_metrics, direction_metrics


def branch_current_pairs(
    rows: list[dict[str, str]],
    candidate_column: str,
) -> dict[tuple[float, str], list[tuple[float, float]]]:
    result: dict[tuple[float, str], list[tuple[float, float]]] = defaultdict(list)
    for row in rows:
        result[(float(row["bias_V"]), row["carrier"])].append((
            float(row[candidate_column]),
            float(row["sentaurus_current_magnitude_A_per_m2"]),
        ))
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-a", type=Path, required=True)
    parser.add_argument("--run-b", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument(
        "--active-source-floor-fraction", type=float, default=1.0e-6
    )
    args = parser.parse_args()
    if not 0.0 < args.active_source_floor_fraction < 1.0:
        raise ValueError("--active-source-floor-fraction must be in (0,1)")

    result_a_path = args.run_a / "result.json"
    result_b_path = args.run_b / "result.json"
    result_a = json.loads(result_a_path.read_text(encoding="utf-8-sig"))
    result_b = json.loads(result_b_path.read_text(encoding="utf-8-sig"))
    for result in (result_a, result_b):
        if (
            not result.get("observation_only", False)
            or result.get("state_advanced", True)
            or result.get("continuity_feedback_enabled", True)
        ):
            raise RuntimeError("candidate input violates observation-only contract")

    gss_a_path = args.run_a / "gss_aux2" / "result.json"
    gss_b_path = args.run_b / "gss_aux2" / "result.json"
    gss_a = json.loads(gss_a_path.read_text(encoding="utf-8-sig"))
    gss_b = json.loads(gss_b_path.read_text(encoding="utf-8-sig"))
    factor_a = read_rows(args.run_a / "current_proxy_factorization.csv")
    gss_details_a = read_rows(
        args.run_a / "gss_aux2" / "gss_aux2_ownership_details.csv"
    )
    external_a = read_rows(args.run_a / "external_current_substitution.csv")

    deterministic_files = [
        "stage_summary.csv",
        "support_comparison.csv",
        "external_current_substitution.csv",
        "current_proxy_factorization.csv",
        "sg_vector_current_control.csv",
        "gss_aux2/gss_aux2_ownership_details.csv",
        "gss_aux2/gss_aux2_ownership_summary.csv",
    ]
    determinism_rows: list[dict[str, Any]] = []
    for relative in deterministic_files:
        path_a = args.run_a / relative
        path_b = args.run_b / relative
        hash_a = sha256(path_a)
        hash_b = sha256(path_b)
        determinism_rows.append({
            "artifact": relative,
            "run_a_sha256": hash_a,
            "run_b_sha256": hash_b,
            "identical": int(hash_a == hash_b),
        })
    for case_a, case_b in zip(result_a["cases"], result_b["cases"], strict=True):
        if float(case_a["bias_V"]) != float(case_b["bias_V"]):
            raise RuntimeError("run bias ordering mismatch")
        for name in (
            "node", "edge", "triangle", "element", "process",
            "sg_vector_node", "sg_vector_edge", "sg_vector_triangle",
            "sg_vector_element", "sg_vector_process",
        ):
            path_a = Path(case_a["outputs"][name]["path"])
            path_b = Path(case_b["outputs"][name]["path"])
            hash_a = sha256(path_a)
            hash_b = sha256(path_b)
            determinism_rows.append({
                "artifact": f"{case_a['bias_V']}/{name}",
                "run_a_sha256": hash_a,
                "run_b_sha256": hash_b,
                "identical": int(hash_a == hash_b),
            })
    deterministic = all(row["identical"] for row in determinism_rows)

    candidate_current, candidate_source, candidate_direction = candidate_metrics(
        args.run_a, result_a, args.active_source_floor_fraction
    )
    baseline_current_pairs = branch_current_pairs(
        factor_a, "production_proxy_current_A_per_m2"
    )
    negative_current_pairs = branch_current_pairs(
        gss_details_a, "gss_reference_midpoint_proxy_current_A_per_m2"
    )
    baseline_source_pairs = source_pairs_from_triangle(
        external_a, "vela_baseline_source_integral_per_m_s"
    )
    negative_source_by_key: dict[
        tuple[float, int, int, int, int, int], float
    ] = defaultdict(float)
    for row in gss_details_a:
        key = (
            float(row["bias_V"]),
            int(row["cell_id"]),
            int(row["local_edge"]),
            int(row["edge_id"]),
            int(row["node0"]),
            int(row["node1"]),
        )
        negative_source_by_key[key] += float(
            row["gss_reference_midpoint_source_integral_per_m_s"]
        )
    negative_source_pairs: dict[
        float, list[tuple[float, float]]
    ] = defaultdict(list)
    for row in external_a:
        bias = float(row["bias_V"])
        key = (
            bias,
            int(row["cell_id"]),
            int(row["local_edge"]),
            int(row["edge_id"]),
            int(row["node0"]),
            int(row["node1"]),
        )
        negative_source_pairs[bias].append((
            negative_source_by_key[key],
            float(row[
                "sentaurus_generation_vela_geometry_source_integral_per_m_s"
            ]),
        ))

    gss_by_bias = {
        float(case["bias_V"]): case for case in gss_a["cases"]
    }
    scorecard: list[dict[str, Any]] = []
    for case in result_a["cases"]:
        bias = float(case["bias_V"])
        sent_total = float(case["sentaurus_total_source_integral_per_m_s"])
        branch_data = {
            "production_triangle_baseline": {
                "total": float(case["vela_total_source_integral_per_m_s"]),
                "source": active_dex_metrics(
                    baseline_source_pairs[bias],
                    args.active_source_floor_fraction,
                ),
                "current_pairs": baseline_current_pairs,
            },
            "sign_correct_midpoint_negative_control": {
                "total": float(gss_by_bias[bias][
                    "gss_reference_midpoint_proxy_total_source_integral_per_m_s"
                ]),
                "source": active_dex_metrics(
                    negative_source_pairs[bias],
                    args.active_source_floor_fraction,
                ),
                "current_pairs": negative_current_pairs,
            },
            "complete_element_edge_sg_gss_laux": {
                "total": float(case["sg_vector_total_source_integral_per_m_s"]),
                "source": candidate_source[bias],
                "current_pairs": None,
            },
        }
        for branch, data in branch_data.items():
            ratio = data["total"] / sent_total
            row: dict[str, Any] = {
                "bias_V": bias,
                "branch": branch,
                "total_source_integral_per_m_s": data["total"],
                "sentaurus_total_source_integral_per_m_s": sent_total,
                "integrated_source_ratio": ratio,
                "integrated_abs_error_dex": abs(math.log10(ratio)),
                "local_source_active_count": data["source"]["active_count"],
                "local_source_median_abs_error_dex": data["source"][
                    "median_abs_error_dex"
                ],
                "local_source_p95_abs_error_dex": data["source"][
                    "p95_abs_error_dex"
                ],
                "local_source_maximum_abs_error_dex": data["source"][
                    "maximum_abs_error_dex"
                ],
            }
            for carrier in ("electron", "hole"):
                metrics = (
                    candidate_current[bias][carrier]
                    if data["current_pairs"] is None
                    else active_dex_metrics(
                        data["current_pairs"][(bias, carrier)]
                    )
                )
                row[f"{carrier}_current_active_count"] = metrics["active_count"]
                row[f"{carrier}_current_median_abs_error_dex"] = metrics[
                    "median_abs_error_dex"
                ]
                row[f"{carrier}_current_p95_abs_error_dex"] = metrics[
                    "p95_abs_error_dex"
                ]
                row[f"{carrier}_current_maximum_abs_error_dex"] = metrics[
                    "maximum_abs_error_dex"
                ]
                row[f"{carrier}_current_sign_agreement_fraction"] = (
                    candidate_direction[bias][
                        f"{carrier}_sign_agreement_fraction"
                    ]
                    if data["current_pairs"] is None else ""
                )
                row[f"{carrier}_current_median_angle_deg"] = (
                    candidate_direction[bias][f"{carrier}_median_angle_deg"]
                    if data["current_pairs"] is None else ""
                )
                row[f"{carrier}_current_p95_angle_deg"] = (
                    candidate_direction[bias][f"{carrier}_p95_angle_deg"]
                    if data["current_pairs"] is None else ""
                )
            scorecard.append(row)

    candidate_rows = [
        row for row in scorecard
        if row["branch"] == "complete_element_edge_sg_gss_laux"
    ]
    fixed_state_gates = {
        "observation_only": True,
        "duplicate_run_determinism": deterministic,
        "integrated_source_within_0p02_relative": all(
            abs(float(row["integrated_source_ratio"]) - 1.0) <= 0.02
            for row in candidate_rows
        ),
        "matching_current_median_p95": all(
            float(row[f"{carrier}_current_median_abs_error_dex"]) <= 0.05
            and float(row[f"{carrier}_current_p95_abs_error_dex"]) <= 0.15
            for row in candidate_rows
            for carrier in ("electron", "hole")
        ),
        "active_source_median_maximum": all(
            float(row["local_source_median_abs_error_dex"]) <= 0.10
            and float(row["local_source_maximum_abs_error_dex"]) <= 0.30
            for row in candidate_rows
        ),
        "nonzero_vector_direction_agreement": all(
            float(row[f"{carrier}_current_sign_agreement_fraction"]) == 1.0
            for row in candidate_rows
            for carrier in ("electron", "hole")
        ),
        "sign_only_negative_control_rejected": all(
            0.45 <= float(row["integrated_source_ratio"]) <= 0.55
            for row in scorecard
            if row["branch"] == "sign_correct_midpoint_negative_control"
        ),
    }
    fixed_state_prequalified = all(fixed_state_gates.values())
    typed_outcome = (
        "complete_sg_vector_fixed_state_prequalified"
        if fixed_state_prequalified
        else "complete_sg_vector_fixed_state_gate_failed"
    )

    args.out_dir.mkdir(parents=True, exist_ok=True)
    scorecard_path = args.out_dir / "task7_frozen_candidate_scorecard.csv"
    determinism_path = args.out_dir / "determinism.csv"
    write_rows(scorecard_path, scorecard)
    write_rows(determinism_path, determinism_rows)
    output = {
        "schema": "vela.task7_frozen_sg_candidate.v1",
        "typed_outcome": typed_outcome,
        "fixed_state_prequalified": fixed_state_prequalified,
        "task8_authorized": False,
        "production_default_change_authorized": False,
        "observation_only": True,
        "state_advanced": False,
        "continuity_feedback_enabled": False,
        "active_source_floor_fraction": args.active_source_floor_fraction,
        "fixed_state_gates": fixed_state_gates,
        "run_a": {
            "path": str(args.run_a.resolve()),
            "result_sha256": sha256(result_a_path),
            "gss_audit_sha256": sha256(gss_a_path),
        },
        "run_b": {
            "path": str(args.run_b.resolve()),
            "result_sha256": sha256(result_b_path),
            "gss_audit_sha256": sha256(gss_b_path),
        },
        "artifacts": {
            "scorecard": {
                "path": str(scorecard_path.resolve()),
                "sha256": sha256(scorecard_path),
            },
            "determinism": {
                "path": str(determinism_path.resolve()),
                "sha256": sha256(determinism_path),
            },
        },
    }
    (args.out_dir / "result.json").write_text(
        json.dumps(output, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
