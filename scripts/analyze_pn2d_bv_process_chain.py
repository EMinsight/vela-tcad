#!/usr/bin/env python3
"""Fail-closed paired PN2D BV process-chain difference localization."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable


CHAIN_SCHEMA = "vela.pn2d_bv_process_chain_input.v1"
PROCESS_SCHEMA = "vela.pn2d_bv_process_run.v1"
OUTPUT_SCHEMA = "vela.pn2d_bv_process_chain_analysis.v1"
STAGES = (
    "state",
    "density",
    "drive",
    "mobility",
    "current",
    "alpha",
    "generation",
    "geometric_source",
    "residual_jacobian",
    "newton_update",
    "terminal_current",
)
STAGE_INDEX = {stage: index for index, stage in enumerate(STAGES)}
QUANTITY_STAGE = {
    "potential": "state",
    "quasi_fermi": "state",
    "density": "density",
    "electric_field": "drive",
    "quasi_fermi_gradient": "drive",
    "mobility": "mobility",
    "velocity": "current",
    "current_density": "current",
    "avalanche_alpha": "alpha",
    "avalanche_generation": "generation",
    "integrated_source": "geometric_source",
    "residual": "residual_jacobian",
    "jacobian": "residual_jacobian",
    "newton_update": "newton_update",
    "terminal_current": "terminal_current",
}
BRANCHES = ("avalanche_off", "iic_postprocess", "avalanche_on")
RELATIVE_THRESHOLD = 0.05
LOG_THRESHOLD_DEX = 0.05
ABSOLUTE_FLOOR = 1.0e-30
TAIL_FRACTION = 1.0e-12
ACTIVE_SOURCE_FRACTION = 1.0e-3
KNEE_BIASES_V = (-19.7, -19.8, -19.85, -19.9, -19.95, -20.0)
PROVENANCE_PRIORITY = {
    "solver_used": 5,
    "native": 4,
    "operator_replay": 3,
    "reconstructed": 2,
    "postprocessed": 1,
    "unknown": 0,
}

STAGE_SUMMARY_FIELDS = (
    "comparison_id", "bias_V", "stage", "matched_records",
    "missing_left", "missing_right", "active_records",
    "max_absolute_difference", "max_relative_difference",
    "max_log_error_dex", "max_vector_angle_deg", "worst_support_key",
    "departed",
)
SUPPORT_SUMMARY_FIELDS = (
    "comparison_id", "bias_V", "left_hotspot", "right_hotspot",
    "active_support_overlap", "left_centroid_x_um", "left_centroid_y_um",
    "right_centroid_x_um", "right_centroid_y_um",
    "left_cumulative_10_support", "left_cumulative_50_support",
    "left_cumulative_90_support", "right_cumulative_10_support",
    "right_cumulative_50_support", "right_cumulative_90_support",
)
HOTSPOT_FIELDS = (
    "comparison_id", "bias_V", "stage", "carrier", "quantity",
    "support_key", "left_magnitude", "right_magnitude",
    "signed_difference", "relative_difference", "log_error_dex",
)
CLOSURE_FIELDS = (
    "simulator", "branch", "bias_V", "source_native",
    "source_reintegrated", "source_relative_error", "terminal_source",
    "terminal_current", "terminal_relative_error", "passed",
)
NEWTON_FIELDS = (
    "comparison_id", "bias_V", "support_key", "carrier", "quantity",
    "left_value", "right_value", "difference",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def finite(value: Any, label: str) -> float:
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        raise ValueError(f"{label} must be numeric")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{label} must be finite")
    return result


def normalized_record(record: dict[str, Any], simulator: str) -> dict[str, Any]:
    stage = record.get("stage") or QUANTITY_STAGE.get(str(record.get("quantity")))
    if stage not in STAGE_INDEX:
        raise ValueError(f"unknown process stage {stage!r}")
    values = record.get("values")
    if not isinstance(values, list) or not values:
        raise ValueError("record values must be a non-empty list")
    coordinates = record.get("coordinates_um", [])
    if coordinates and len(coordinates) < 2:
        raise ValueError("coordinates_um must contain x and y")
    return {
        "simulator": simulator,
        "branch": str(record["branch"]),
        "bias_V": finite(record["bias_V"], "bias_V"),
        "stage": stage,
        "quantity": str(record["quantity"]),
        "carrier": str(record.get("carrier", "none")),
        "support_kind": str(record.get("support_kind", "aggregate")),
        "support_key": str(record.get("support_key", "device")),
        "values": [finite(value, "values") for value in values],
        "unit": str(record.get("unit", "1")),
        "coordinates_um": [finite(value, "coordinates_um") for value in coordinates],
        "provenance": str(record.get("provenance", "unknown")),
    }


def from_process_run(payload: dict[str, Any]) -> dict[str, Any]:
    simulator = str(payload["simulator"])
    records: list[dict[str, Any]] = []
    for field in payload.get("field_records", []):
        record = {
            "branch": field["branch"],
            "bias_V": field["requested_bias_V"],
            "quantity": field["quantity"],
            "carrier": field["carrier"],
            "support_kind": field["support_kind"],
            "support_key": field["support_key"],
            "values": field["values"],
            "unit": field["unit"],
            "coordinates_um": field.get("coordinates_um", []),
            "provenance": field["provenance"],
        }
        if record["quantity"] in QUANTITY_STAGE:
            records.append(normalized_record(record, simulator))
    for aggregate in payload.get("aggregate_records", []):
        record = {
            "branch": aggregate["branch"],
            "bias_V": aggregate["requested_bias_V"],
            "quantity": aggregate["quantity"],
            "carrier": aggregate["carrier"],
            "support_kind": "aggregate",
            "support_key": "device",
            "values": [aggregate["value"]],
            "unit": aggregate["unit"],
            "provenance": aggregate["provenance"],
        }
        if record["quantity"] in QUANTITY_STAGE:
            records.append(normalized_record(record, simulator))
    closures: list[dict[str, Any]] = []
    grouped_aggregates: dict[tuple[str, float], list[dict[str, Any]]] = defaultdict(list)
    for aggregate in payload.get("aggregate_records", []):
        grouped_aggregates[
            (str(aggregate["branch"]), float(aggregate["requested_bias_V"]))
        ].append(aggregate)
    for (branch, bias), aggregates in sorted(grouped_aggregates.items()):
        source_totals = [
            aggregate for aggregate in aggregates
            if aggregate["quantity"] == "integrated_source"
            and aggregate["carrier"] == "total"
        ]
        terminals = {
            aggregate["carrier"]: float(aggregate["value"])
            for aggregate in aggregates
            if aggregate["quantity"] == "terminal_current"
        }
        native = next(
            (
                float(aggregate["value"])
                for aggregate in source_totals
                if aggregate["provenance"] in {"native", "solver_used", "postprocessed"}
            ),
            None,
        )
        replay = next(
            (
                float(aggregate["value"])
                for aggregate in source_totals
                if aggregate["provenance"] == "operator_replay"
            ),
            None,
        )
        if (
            native is not None
            and replay is not None
            and {"electron", "hole", "total"} <= set(terminals)
        ):
            closures.append(
                {
                    "branch": branch,
                    "bias_V": bias,
                    "source_native": native,
                    "source_reintegrated": replay,
                    "terminal_source": (
                        terminals["electron"] - terminals["hole"]
                        if simulator == "vela"
                        else terminals["electron"] + terminals["hole"]
                    ),
                    "terminal_current": terminals["total"],
                }
            )
    return {
        "schema": CHAIN_SCHEMA,
        "simulator": simulator,
        "records": records,
        "closures": closures,
        "newton_updates": [],
    }


def load_input(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("schema") == PROCESS_SCHEMA:
        return from_process_run(payload)
    if payload.get("schema") != CHAIN_SCHEMA:
        raise ValueError(f"{path}: unsupported schema {payload.get('schema')!r}")
    simulator = str(payload["simulator"])
    records = [
        normalized_record(dict(record), simulator)
        for record in payload.get("records", [])
    ]
    return {
        "schema": CHAIN_SCHEMA,
        "simulator": simulator,
        "records": records,
        "closures": list(payload.get("closures", [])),
        "newton_updates": list(payload.get("newton_updates", [])),
    }


def record_key(record: dict[str, Any]) -> tuple[Any, ...]:
    return (
        record["branch"], record["bias_V"], record["stage"],
        record["quantity"], record["carrier"], record["support_kind"],
        record["support_key"], record["unit"],
    )


def magnitude(values: Iterable[float]) -> float:
    return math.sqrt(sum(float(value) ** 2 for value in values))


def vector_angle(left: list[float], right: list[float]) -> float:
    if len(left) < 2 or len(right) < 2:
        return 0.0
    left_norm = magnitude(left)
    right_norm = magnitude(right)
    if left_norm == 0.0 or right_norm == 0.0:
        return 0.0
    cosine = sum(a * b for a, b in zip(left, right)) / (left_norm * right_norm)
    return math.degrees(math.acos(max(-1.0, min(1.0, cosine))))


def comparison_specs(
    sentaurus: dict[str, Any], vela: dict[str, Any]
) -> list[tuple[str, dict[str, Any], str, dict[str, Any], str]]:
    result = [
        (f"sentaurus_vs_vela_{branch}", sentaurus, branch, vela, branch)
        for branch in BRANCHES
    ]
    result.extend(
        [
            (
                "sentaurus_iic_vs_on", sentaurus, "iic_postprocess",
                sentaurus, "avalanche_on",
            ),
            (
                "vela_iic_vs_on", vela, "iic_postprocess",
                vela, "avalanche_on",
            ),
        ]
    )
    return result


def selected(
    dataset: dict[str, Any], branch: str
) -> dict[tuple[Any, ...], dict[str, Any]]:
    result: dict[tuple[Any, ...], dict[str, Any]] = {}
    for record in dataset["records"]:
        if record["branch"] != branch:
            continue
        key = record_key({**record, "branch": "_"})
        if key in result:
            previous = result[key]
            old_priority = PROVENANCE_PRIORITY.get(previous["provenance"], -1)
            new_priority = PROVENANCE_PRIORITY.get(record["provenance"], -1)
            if new_priority == old_priority:
                raise ValueError(
                    f"duplicate normalized record at equal provenance priority {key}"
                )
            if new_priority < old_priority:
                continue
        result[key] = record
    return result


def stage_rows_for_comparison(
    comparison_id: str,
    left_data: dict[str, Any],
    left_branch: str,
    right_data: dict[str, Any],
    right_branch: str,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    left = selected(left_data, left_branch)
    right = selected(right_data, right_branch)
    biases = sorted(
        {float(key[1]) for key in left} | {float(key[1]) for key in right},
        reverse=True,
    )
    summaries: list[dict[str, Any]] = []
    details: list[dict[str, Any]] = []
    for bias in biases:
        for stage in STAGES:
            left_stage = {
                key: value for key, value in left.items()
                if key[1] == bias and key[2] == stage
            }
            right_stage = {
                key: value for key, value in right.items()
                if key[1] == bias and key[2] == stage
            }
            common = sorted(set(left_stage) & set(right_stage), key=str)
            stage_scale = max(
                [
                    magnitude(record["values"])
                    for record in (*left_stage.values(), *right_stage.values())
                ]
                or [0.0]
            )
            active_floor = max(ABSOLUTE_FLOOR, stage_scale * TAIL_FRACTION)
            max_abs = max_rel = max_log = max_angle = 0.0
            worst_support = ""
            active = 0
            for key in common:
                lrec, rrec = left_stage[key], right_stage[key]
                lmag, rmag = magnitude(lrec["values"]), magnitude(rrec["values"])
                if max(lmag, rmag) < active_floor:
                    continue
                active += 1
                difference = magnitude(
                    [a - b for a, b in zip(lrec["values"], rrec["values"])]
                )
                relative = difference / max(lmag, rmag, active_floor)
                log_error = (
                    abs(math.log10(rmag / lmag))
                    if lmag > active_floor and rmag > active_floor else 0.0
                )
                angle = vector_angle(lrec["values"], rrec["values"])
                if relative > max_rel:
                    worst_support = lrec["support_key"]
                max_abs = max(max_abs, difference)
                max_rel = max(max_rel, relative)
                max_log = max(max_log, log_error)
                max_angle = max(max_angle, angle)
                details.append(
                    {
                        "comparison_id": comparison_id,
                        "bias_V": bias,
                        "stage": stage,
                        "carrier": lrec["carrier"],
                        "quantity": lrec["quantity"],
                        "support_key": lrec["support_key"],
                        "left_magnitude": lmag,
                        "right_magnitude": rmag,
                        "signed_difference": rrec["values"][0] - lrec["values"][0],
                        "relative_difference": relative,
                        "log_error_dex": log_error,
                    }
                )
            missing_left = len(set(right_stage) - set(left_stage))
            missing_right = len(set(left_stage) - set(right_stage))
            departed = int(
                active > 0
                and max_abs > active_floor
                and (max_rel > RELATIVE_THRESHOLD or max_log > LOG_THRESHOLD_DEX)
            )
            summaries.append(
                {
                    "comparison_id": comparison_id,
                    "bias_V": bias,
                    "stage": stage,
                    "matched_records": len(common),
                    "missing_left": missing_left,
                    "missing_right": missing_right,
                    "active_records": active,
                    "max_absolute_difference": max_abs,
                    "max_relative_difference": max_rel,
                    "max_log_error_dex": max_log,
                    "max_vector_angle_deg": max_angle,
                    "worst_support_key": worst_support,
                    "departed": departed,
                }
            )
    return summaries, details


def source_records(
    dataset: dict[str, Any], branch: str, bias: float
) -> list[dict[str, Any]]:
    return [
        record for record in dataset["records"]
        if record["branch"] == branch
        and record["bias_V"] == bias
        and record["stage"] in {"generation", "geometric_source"}
        and record["support_kind"] != "aggregate"
    ]


def active_support(records: list[dict[str, Any]]) -> tuple[set[str], str, tuple[float, float], list[str]]:
    if not records:
        return set(), "", (math.nan, math.nan), []
    by_support: dict[str, tuple[float, list[float]]] = {}
    for record in records:
        value = magnitude(record["values"])
        old = by_support.get(record["support_key"], (0.0, record["coordinates_um"]))
        by_support[record["support_key"]] = (old[0] + value, old[1])
    peak_key = max(by_support, key=lambda key: (by_support[key][0], key))
    peak = by_support[peak_key][0]
    active = {
        key for key, (value, _) in by_support.items()
        if peak > 0.0 and value >= peak * ACTIVE_SOURCE_FRACTION
    }
    total = sum(value for value, _ in by_support.values())
    if total > 0.0:
        x = sum(
            value * coords[0] for value, coords in by_support.values()
            if len(coords) >= 2
        ) / total
        y = sum(
            value * coords[1] for value, coords in by_support.values()
            if len(coords) >= 2
        ) / total
    else:
        x = y = math.nan
    ordered = sorted(by_support, key=lambda key: (-by_support[key][0], key))
    cumulative: list[str] = []
    running = 0.0
    targets = (0.1, 0.5, 0.9)
    target_index = 0
    for key in ordered:
        running += by_support[key][0]
        while target_index < len(targets) and (
            total == 0.0 or running / total >= targets[target_index]
        ):
            cumulative.append(key)
            target_index += 1
    while len(cumulative) < 3:
        cumulative.append("")
    return active, peak_key, (x, y), cumulative


def support_summary(
    comparison_id: str,
    left: dict[str, Any],
    left_branch: str,
    right: dict[str, Any],
    right_branch: str,
) -> list[dict[str, Any]]:
    biases = sorted(
        {
            record["bias_V"] for record in left["records"]
            if record["branch"] == left_branch
        }
        | {
            record["bias_V"] for record in right["records"]
            if record["branch"] == right_branch
        },
        reverse=True,
    )
    result = []
    for bias in biases:
        la, lh, lc, lq = active_support(source_records(left, left_branch, bias))
        ra, rh, rc, rq = active_support(source_records(right, right_branch, bias))
        union = la | ra
        overlap = len(la & ra) / len(union) if union else math.nan
        result.append(
            {
                "comparison_id": comparison_id,
                "bias_V": bias,
                "left_hotspot": lh,
                "right_hotspot": rh,
                "active_support_overlap": overlap,
                "left_centroid_x_um": lc[0],
                "left_centroid_y_um": lc[1],
                "right_centroid_x_um": rc[0],
                "right_centroid_y_um": rc[1],
                "left_cumulative_10_support": lq[0],
                "left_cumulative_50_support": lq[1],
                "left_cumulative_90_support": lq[2],
                "right_cumulative_10_support": rq[0],
                "right_cumulative_50_support": rq[1],
                "right_cumulative_90_support": rq[2],
            }
        )
    return result


def first_departures(stage_rows: list[dict[str, Any]]) -> dict[str, Any]:
    by_comparison: dict[str, dict[float, str]] = defaultdict(dict)
    for row in stage_rows:
        if int(row["departed"]):
            bias = float(row["bias_V"])
            previous = by_comparison[row["comparison_id"]].get(bias)
            if previous is None or STAGE_INDEX[row["stage"]] < STAGE_INDEX[previous]:
                by_comparison[row["comparison_id"]][bias] = row["stage"]
    comparisons: dict[str, Any] = {}
    causal_candidates: list[tuple[int, str, str, list[float]]] = []
    for comparison, bias_stages in sorted(by_comparison.items()):
        ordered_biases = [
            bias for bias in KNEE_BIASES_V
            if bias in bias_stages
        ]
        evidence: list[tuple[str, list[float]]] = []
        for first, second in zip(ordered_biases, ordered_biases[1:]):
            if bias_stages[first] == bias_stages[second]:
                evidence.append((bias_stages[first], [first, second]))
        if evidence:
            stage, biases = min(evidence, key=lambda item: STAGE_INDEX[item[0]])
            causal_candidates.append((STAGE_INDEX[stage], comparison, stage, biases))
        comparisons[comparison] = {
            "by_bias": {f"{bias:.17g}": bias_stages[bias] for bias in ordered_biases},
            "adjacent_bias_evidence": [
                {"stage": stage, "biases_V": biases} for stage, biases in evidence
            ],
        }
    if not causal_candidates:
        return {
            "comparisons": comparisons,
            "causal_stage": None,
            "comparison_id": None,
            "adjacent_biases_V": [],
        }
    _, comparison, stage, biases = min(causal_candidates)
    return {
        "comparisons": comparisons,
        "causal_stage": stage,
        "comparison_id": comparison,
        "adjacent_biases_V": biases,
    }


def classify(first: dict[str, Any], missing_observation: bool) -> str:
    if missing_observation or first["causal_stage"] is None:
        return "insufficient_observation"
    stage = str(first["causal_stage"])
    comparison = str(first["comparison_id"])
    if stage in {"residual_jacobian", "newton_update"}:
        return "continuation_solver_path_cause"
    if comparison.endswith("iic_vs_on"):
        if stage in {"state", "density"}:
            return "density_qfp_feedback_cause"
        if stage in {"mobility", "current"}:
            return "mobility_current_feedback_cause"
        if stage in {"generation", "geometric_source"}:
            return "source_support_feedback_cause"
    return "fixed_state_operator_cause"


def closure_rows(datasets: Iterable[dict[str, Any]]) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    for dataset in datasets:
        for row in dataset.get("closures", []):
            native = finite(row.get("source_native", 0.0), "source_native")
            replay = finite(row.get("source_reintegrated", 0.0), "source_reintegrated")
            terminal_source = finite(row.get("terminal_source", 0.0), "terminal_source")
            terminal = finite(row.get("terminal_current", 0.0), "terminal_current")
            source_error = abs(native - replay) / max(abs(native), abs(replay), ABSOLUTE_FLOOR)
            terminal_error = abs(terminal_source - terminal) / max(
                abs(terminal_source), abs(terminal), ABSOLUTE_FLOOR
            )
            result.append(
                {
                    "simulator": dataset["simulator"],
                    "branch": row["branch"],
                    "bias_V": row["bias_V"],
                    "source_native": native,
                    "source_reintegrated": replay,
                    "source_relative_error": source_error,
                    "terminal_source": terminal_source,
                    "terminal_current": terminal,
                    "terminal_relative_error": terminal_error,
                    "passed": int(source_error <= 1.0e-12 and terminal_error <= 1.0e-12),
                }
            )
    return result


def write_csv(path: Path, fields: tuple[str, ...], rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field) for field in fields})


def write_svg(path: Path, title: str, labels: list[str], values: list[float]) -> None:
    width, height = 960, 420
    bars = max(1, len(labels))
    bar_width = 820 / bars
    maximum = max(values or [1.0]) or 1.0
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        f'<text x="40" y="35" font-family="sans-serif" font-size="20">{title}</text>',
    ]
    for index, (label, value) in enumerate(zip(labels, values)):
        x = 80 + index * bar_width
        h = 300 * value / maximum
        parts.append(
            f'<rect x="{x:.2f}" y="{360-h:.2f}" width="{bar_width*0.75:.2f}" '
            f'height="{h:.2f}" fill="#2f6f9f"/>'
        )
        parts.append(
            f'<text x="{x:.2f}" y="382" font-family="sans-serif" '
            f'font-size="10" transform="rotate(35 {x:.2f} 382)">{label}</text>'
        )
    parts.append("</svg>")
    path.write_text("\n".join(parts) + "\n", encoding="utf-8")


def analyze(
    sentaurus: dict[str, Any] | None,
    vela: dict[str, Any] | None,
    output_root: Path,
    *,
    input_artifacts: dict[str, dict[str, str] | None] | None = None,
) -> dict[str, Any]:
    output_root.mkdir(parents=True, exist_ok=True)
    if sentaurus is None or vela is None:
        stage_rows: list[dict[str, Any]] = []
        detail_rows: list[dict[str, Any]] = []
        support_rows: list[dict[str, Any]] = []
        first = {
            "comparisons": {}, "causal_stage": None,
            "comparison_id": None, "adjacent_biases_V": [],
        }
        missing = True
        closures: list[dict[str, Any]] = []
        missing_stage_observations: dict[str, list[str]] = {
            "all": list(STAGES)
        }
    else:
        stage_rows, detail_rows, support_rows = [], [], []
        specs = comparison_specs(sentaurus, vela)
        for comparison_id, left, left_branch, right, right_branch in specs:
            summary, details = stage_rows_for_comparison(
                comparison_id, left, left_branch, right, right_branch
            )
            stage_rows.extend(summary)
            detail_rows.extend(details)
            support_rows.extend(
                support_summary(
                    comparison_id, left, left_branch, right, right_branch
                )
            )
        first = first_departures(stage_rows)
        required_comparisons = {
            *(f"sentaurus_vs_vela_{branch}" for branch in BRANCHES),
            "sentaurus_iic_vs_on",
            "vela_iic_vs_on",
        }
        covered = {
            (row["comparison_id"], row["stage"])
            for row in stage_rows
            if int(row["matched_records"]) > 0
        }
        missing = any(
            (comparison, stage) not in covered
            for comparison in required_comparisons
            for stage in STAGES
        )
        missing_stage_observations = {
            comparison: [
                stage for stage in STAGES
                if (comparison, stage) not in covered
            ]
            for comparison in sorted(required_comparisons)
        }
        missing_stage_observations = {
            comparison: stages
            for comparison, stages in missing_stage_observations.items()
            if stages
        }
        closures = closure_rows((sentaurus, vela))
        closure_simulators = {row["simulator"] for row in closures}
        if (
            closure_simulators != {sentaurus["simulator"], vela["simulator"]}
            or not closures
            or not all(int(row["passed"]) for row in closures)
        ):
            missing = True

    outcome = classify(first, missing)
    accepted_causal_stage = (
        first["causal_stage"]
        if outcome != "insufficient_observation"
        else None
    )
    accepted_comparison = (
        first["comparison_id"]
        if outcome != "insufficient_observation"
        else None
    )
    accepted_biases = (
        first["adjacent_biases_V"]
        if outcome != "insufficient_observation"
        else []
    )
    write_csv(output_root / "stage_summary.csv", STAGE_SUMMARY_FIELDS, stage_rows)
    write_csv(output_root / "support_summary.csv", SUPPORT_SUMMARY_FIELDS, support_rows)
    hotspot = [
        row for row in detail_rows
        if any(
            row["support_key"] == support.get("left_hotspot")
            for support in support_rows
            if support["comparison_id"] == row["comparison_id"]
            and float(support["bias_V"]) == float(row["bias_V"])
        )
    ]
    write_csv(output_root / "hotspot_chain.csv", HOTSPOT_FIELDS, hotspot)
    write_csv(output_root / "source_terminal_closure.csv", CLOSURE_FIELDS, closures)
    newton_rows = [
        {
            "comparison_id": row["comparison_id"],
            "bias_V": row["bias_V"],
            "support_key": row["support_key"],
            "carrier": row["carrier"],
            "quantity": row["quantity"],
            "left_value": row["left_magnitude"],
            "right_value": row["right_magnitude"],
            "difference": row["signed_difference"],
        }
        for row in detail_rows if row["stage"] == "newton_update"
    ]
    write_csv(output_root / "newton_first_update.csv", NEWTON_FIELDS, newton_rows)
    first_payload = {
        "schema": OUTPUT_SCHEMA,
        "outcome": outcome,
        "causal_stage": accepted_causal_stage,
        "comparison_id": accepted_comparison,
        "adjacent_biases_V": accepted_biases,
        "candidate_causal_stage": first["causal_stage"],
        "candidate_comparison_id": first["comparison_id"],
        "candidate_adjacent_biases_V": first["adjacent_biases_V"],
        "comparisons": first["comparisons"],
        "missing_stage_observations": missing_stage_observations,
    }
    (output_root / "first_departure.json").write_text(
        json.dumps(first_payload, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    stage_max = {
        stage: max(
            [
                float(row["max_relative_difference"])
                for row in stage_rows if row["stage"] == stage
            ]
            or [0.0]
        )
        for stage in STAGES
    }
    write_svg(
        output_root / "process_chain.svg",
        "PN2D BV process-chain maximum relative difference",
        list(STAGES),
        [stage_max[stage] for stage in STAGES],
    )
    write_svg(
        output_root / "hotspot.svg",
        "PN2D BV hotspot-chain relative difference",
        [row["stage"] for row in hotspot[:20]],
        [float(row["relative_difference"]) for row in hotspot[:20]],
    )
    acceptance = {
        "schema": OUTPUT_SCHEMA,
        "outcome": outcome,
        "status": "passed" if outcome != "insufficient_observation" else "failed",
        "missing_observation": missing,
        "causal_stage": accepted_causal_stage,
        "comparison_id": accepted_comparison,
        "adjacent_biases_V": accepted_biases,
        "candidate_causal_stage": first["causal_stage"],
        "candidate_comparison_id": first["comparison_id"],
        "candidate_adjacent_biases_V": first["adjacent_biases_V"],
        "missing_stage_observations": missing_stage_observations,
        "failed_closure_rows": sum(
            1 for row in closures if not int(row["passed"])
        ),
        "thresholds": {
            "relative": RELATIVE_THRESHOLD,
            "log_error_dex": LOG_THRESHOLD_DEX,
            "tail_fraction": TAIL_FRACTION,
            "source_active_fraction": ACTIVE_SOURCE_FRACTION,
            "closure_relative": 1.0e-12,
        },
        "input_artifacts": input_artifacts or {},
        "row_counts": {
            "stage_summary": len(stage_rows),
            "support_summary": len(support_rows),
            "hotspot_chain": len(hotspot),
            "source_terminal_closure": len(closures),
            "newton_first_update": len(newton_rows),
        },
    }
    (output_root / "acceptance.json").write_text(
        json.dumps(acceptance, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return acceptance


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sentaurus", type=Path)
    parser.add_argument("--vela", type=Path)
    parser.add_argument("--output-root", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    artifacts: dict[str, dict[str, str] | None] = {}
    datasets: dict[str, dict[str, Any] | None] = {}
    for name, path in (("sentaurus", args.sentaurus), ("vela", args.vela)):
        if path is None:
            datasets[name] = None
            artifacts[name] = None
            continue
        resolved = path.resolve()
        if not resolved.is_file():
            raise FileNotFoundError(resolved)
        datasets[name] = load_input(resolved)
        artifacts[name] = {"path": str(resolved), "sha256": sha256(resolved)}
    result = analyze(
        datasets["sentaurus"],
        datasets["vela"],
        args.output_root.resolve(),
        input_artifacts=artifacts,
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if result["outcome"] != "insufficient_observation" else 2


if __name__ == "__main__":
    raise SystemExit(main())
