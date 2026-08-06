#!/usr/bin/env python3
"""Extract Sentaurus BVmethods paths and compare them with Vela path traces.

The Sentaurus ion-integral fields are piecewise-constant on the node support of
each breakdown path.  This script treats each distinct positive
MeanIonIntegral plateau as one path, exports that raw support without geometric
guesswork, then assigns the three Vela paths to the three Sentaurus supports by
minimum symmetric nearest-neighbour distance.  When a WriteAll log is
available, ranking values come from the log's electron and hole path integrals;
the plotted MeanIonIntegral plateau is retained only as a geometry label.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
import re
from collections import defaultdict
from pathlib import Path
from statistics import fmean, median


REPO = Path(__file__).resolve().parents[1]
DEFAULT_SENT = (
    REPO / "build-release/reference_tcad/bvmethods_sentaurus2018/run01"
    / "sentaurus_path_export_20260804/neutral"
)
DEFAULT_VELA = (
    REPO / "build-release/reference_tcad/bvmethods_sentaurus2018/run01"
    / "vela_validation/eparallel_extend_7_10p45_20260804/postprocess_only"
)
DEFAULT_OUT = (
    REPO / "build-release/reference_tcad/bvmethods_sentaurus2018/run01"
    / "vela_validation/eparallel_path_compare_20260804"
)
DEFAULT_SENT_LOG = (
    REPO / "build-release/reference_tcad/bvmethods_sentaurus2018/run01"
    / "sentaurus_mean_ionization_controls_20260805/extracted"
    / "mean_control_baseline.log_des.log"
)


_PATH_RE = re.compile(r"^Path number\s+(\d+)\s*$")
_FIELD_RE = re.compile(r"^(Maximum Field|Electron|Hole):\s+([+\-0-9.eE]+)\s*$")


def final_write_all_path_integrals(path: Path) -> list[dict[str, float | int]]:
    """Return distinct final WriteAll paths ranked by arithmetic carrier mean.

    Sentaurus writes many bias points to one log.  A WriteAll inventory is the
    numbered block immediately preceding ``Best Path``.  Coincident local-peak
    aliases can have identical carrier integrals, so they are retained through
    a ``multiplicity`` count while the numeric rank uses one distinct value.
    """
    inventories: list[list[dict[str, float | int]]] = []
    current: list[dict[str, float | int]] = []
    active: dict[str, float | int] | None = None
    for raw_line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw_line.strip()
        path_match = _PATH_RE.match(line)
        if path_match:
            active = {"path_number": int(path_match.group(1))}
            current.append(active)
            continue
        if line == "Best Path":
            complete = [
                row for row in current
                if "electron_ionization_integral" in row
                and "hole_ionization_integral" in row
            ]
            if complete:
                inventories.append(complete)
            current = []
            active = None
            continue
        field_match = _FIELD_RE.match(line)
        if active is None or field_match is None:
            continue
        label, raw_value = field_match.groups()
        key = {
            "Maximum Field": "max_electric_field_V_per_cm",
            "Electron": "electron_ionization_integral",
            "Hole": "hole_ionization_integral",
        }[label]
        active[key] = float(raw_value)

    if not inventories:
        raise RuntimeError(f"No complete Sentaurus WriteAll inventory in {path}")

    distinct: dict[tuple[float, float], dict[str, float | int]] = {}
    for row in inventories[-1]:
        electron = float(row["electron_ionization_integral"])
        hole = float(row["hole_ionization_integral"])
        key = (electron, hole)
        if key not in distinct:
            distinct[key] = dict(row, multiplicity=1)
        else:
            distinct[key]["multiplicity"] = int(distinct[key]["multiplicity"]) + 1
    ranked = list(distinct.values())
    for row in ranked:
        row["mean_ionization_integral"] = 0.5 * (
            float(row["electron_ionization_integral"])
            + float(row["hole_ionization_integral"])
        )
    ranked.sort(
        key=lambda row: (
            -float(row["mean_ionization_integral"]),
            int(row["path_number"]),
        )
    )
    return ranked


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def scalar_field(root: Path, name: str) -> dict[int, float]:
    rows = read_rows(root / "fields" / f"{name}_region3.csv")
    return {int(row["node_id"]): float(row["component0"]) for row in rows}


def optional_scalar_field(root: Path, name: str) -> dict[int, float]:
    path = root / "fields" / f"{name}_region3.csv"
    return scalar_field(root, name) if path.exists() else {}


def vector_magnitude_field(root: Path, name: str) -> dict[int, float]:
    rows = read_rows(root / "fields" / f"{name}_region3.csv")
    return {
        int(row["node_id"]): math.hypot(
            float(row["component0"]), float(row["component1"])
        )
        for row in rows
    }


def distance(a: tuple[float, float], b: tuple[float, float]) -> float:
    return math.hypot(a[0] - b[0], a[1] - b[1])


def nearest_metrics(
    left: list[tuple[float, float]], right: list[tuple[float, float]]
) -> tuple[float, float, float]:
    left_nearest = [min(distance(a, b) for b in right) for a in left]
    right_nearest = [min(distance(b, a) for a in left) for b in right]
    return (
        0.5 * (fmean(left_nearest) + fmean(right_nearest)),
        max(max(left_nearest), max(right_nearest)),
        math.sqrt(fmean(value * value for value in left_nearest)),
    )


def safe_ratio(value: float, reference: float) -> float:
    if reference == 0.0:
        return math.inf if value != 0.0 else 1.0
    return value / reference


def finite_median(values: list[float]) -> float | None:
    finite = [value for value in values if math.isfinite(value)]
    return median(finite) if finite else None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sentaurus-neutral", type=Path, default=DEFAULT_SENT)
    parser.add_argument("--sentaurus-log", type=Path, default=DEFAULT_SENT_LOG)
    parser.add_argument("--vela-dir", type=Path, default=DEFAULT_VELA)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--top-paths", type=int, default=3)
    parser.add_argument("--sentaurus-bias", type=float, default=10.448266730833666)
    args = parser.parse_args()

    sent = args.sentaurus_neutral.resolve()
    vela = args.vela_dir.resolve()
    out = args.out_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)

    nodes = {
        int(row["id"]): (float(row["x_um"]), float(row["y_um"]))
        for row in read_rows(sent / "nodes.csv")
    }
    e_integral = scalar_field(sent, "eIonIntegral")
    h_integral = scalar_field(sent, "hIonIntegral")
    mean_integral = scalar_field(sent, "MeanIonIntegral")
    e_alpha_cm = scalar_field(sent, "eAlphaAvalanche")
    h_alpha_cm = scalar_field(sent, "hAlphaAvalanche")
    efield_v_cm = vector_magnitude_field(sent, "ElectricField")
    eparallel_v_cm = optional_scalar_field(sent, "eEparallel")
    hparallel_v_cm = optional_scalar_field(sent, "hEparallel")

    plateau_nodes: dict[float, list[int]] = defaultdict(list)
    for node_id, value in mean_integral.items():
        if value > 0.0:
            plateau_nodes[value].append(node_id)
    plateau_values = sorted(plateau_nodes, reverse=True)[: args.top_paths]
    if len(plateau_values) < args.top_paths:
        raise RuntimeError("Sentaurus TDR contains fewer positive path plateaus than requested")

    sent_paths: dict[int, list[int]] = {}
    sent_support_rows: list[dict[str, object]] = []
    sent_summary: dict[int, dict[str, object]] = {}
    for rank, plateau in enumerate(plateau_values, start=1):
        support = plateau_nodes[plateau]
        sent_paths[rank] = support
        xs = [nodes[node_id][0] for node_id in support]
        ys = [nodes[node_id][1] for node_id in support]
        for node_id in support:
            x_um, y_um = nodes[node_id]
            sent_support_rows.append({
                "sentaurus_path_rank": rank,
                "node_id": node_id,
                "x_um": x_um,
                "y_um": y_um,
                "electric_field_V_per_cm": efield_v_cm[node_id],
                "electron_eparallel_V_per_cm": eparallel_v_cm.get(node_id),
                "hole_eparallel_V_per_cm": hparallel_v_cm.get(node_id),
                "electron_alpha_cm_inv": e_alpha_cm[node_id],
                "hole_alpha_cm_inv": h_alpha_cm[node_id],
                "electron_ionization_integral": e_integral[node_id],
                "hole_ionization_integral": h_integral[node_id],
                "mean_ionization_integral": mean_integral[node_id],
            })
        sent_summary[rank] = {
            "sentaurus_path_rank": rank,
            "support_node_count": len(support),
            "x_min_um": min(xs),
            "x_max_um": max(xs),
            "y_min_um": min(ys),
            "y_max_um": max(ys),
            "electron_ionization_integral": median(e_integral[i] for i in support),
            "hole_ionization_integral": median(h_integral[i] for i in support),
            "mean_ionization_integral": plateau,
            "tdr_mean_ionization_integral_plateau": plateau,
            "mean_ionization_integral_source": "tdr_plateau_fallback",
            "max_electric_field_V_per_m": max(efield_v_cm[i] for i in support) * 100.0,
            "max_electron_eparallel_V_per_m": (
                max(eparallel_v_cm[i] for i in support) * 100.0
                if eparallel_v_cm else None
            ),
            "max_hole_eparallel_V_per_m": (
                max(hparallel_v_cm[i] for i in support) * 100.0
                if hparallel_v_cm else None
            ),
            "max_electron_alpha_m_inv": max(e_alpha_cm[i] for i in support) * 100.0,
            "max_hole_alpha_m_inv": max(h_alpha_cm[i] for i in support) * 100.0,
        }
    sent_log = args.sentaurus_log.resolve()
    if sent_log.is_file():
        log_paths = final_write_all_path_integrals(sent_log)
        if len(log_paths) < args.top_paths:
            raise RuntimeError(
                "Sentaurus WriteAll log contains fewer distinct paths than requested"
            )
        for rank, log_path in enumerate(log_paths[: args.top_paths], start=1):
            sent_summary[rank].update({
                "write_all_path_number": int(log_path["path_number"]),
                "write_all_path_multiplicity": int(log_path["multiplicity"]),
                "electron_ionization_integral": float(
                    log_path["electron_ionization_integral"]
                ),
                "hole_ionization_integral": float(
                    log_path["hole_ionization_integral"]
                ),
                "mean_ionization_integral": float(
                    log_path["mean_ionization_integral"]
                ),
                "mean_ionization_integral_source": (
                    "WriteAll arithmetic mean of electron and hole path integrals"
                ),
            })
    write_rows(out / "sentaurus_top3_path_support_nodes.csv", sent_support_rows)

    vela_summary_rows = read_rows(vela / "path_ionization_integrals.csv")
    vela_segment_rows = read_rows(vela / "path_ionization_integral_segments.csv")
    if not vela_summary_rows or not vela_segment_rows:
        raise RuntimeError("Vela path summary or segment trace is empty")
    vela_bias = max(float(row["bias_V"]) for row in vela_summary_rows)
    vela_summary_rows = [
        row for row in vela_summary_rows
        if math.isclose(float(row["bias_V"]), vela_bias, abs_tol=1.0e-10)
        and int(row["path_rank"]) <= args.top_paths
    ]
    vela_segment_rows = [
        row for row in vela_segment_rows
        if math.isclose(float(row["bias_V"]), vela_bias, abs_tol=1.0e-10)
        and int(row["path_rank"]) <= args.top_paths
    ]
    vela_segments: dict[int, list[dict[str, str]]] = defaultdict(list)
    for row in vela_segment_rows:
        vela_segments[int(row["path_rank"])].append(row)
    for rows in vela_segments.values():
        rows.sort(key=lambda row: int(row["segment_index"]))

    vela_points = {
        rank: [
            (
                0.5 * (float(row["x0_um"]) + float(row["x1_um"])),
                0.5 * (float(row["y0_um"]) + float(row["y1_um"])),
            )
            for row in rows
        ]
        for rank, rows in vela_segments.items()
    }
    sent_points = {
        rank: [nodes[node_id] for node_id in support]
        for rank, support in sent_paths.items()
    }
    ranks = list(range(1, args.top_paths + 1))
    assignments = []
    for permutation in itertools.permutations(ranks):
        metrics = {
            vela_rank: nearest_metrics(
                vela_points[vela_rank], sent_points[sent_rank]
            )
            for vela_rank, sent_rank in zip(ranks, permutation)
        }
        assignments.append((sum(value[0] for value in metrics.values()), permutation, metrics))
    _, permutation, geometry_metrics = min(assignments, key=lambda item: item[0])
    mapping = dict(zip(ranks, permutation))

    matched_segment_rows: list[dict[str, object]] = []
    compare_rows: list[dict[str, object]] = []
    summary_by_rank = {int(row["path_rank"]): row for row in vela_summary_rows}
    rankwise_rows: list[dict[str, object]] = []
    for rank in ranks:
        segment_rows = vela_segments[rank]
        vela_row = summary_by_rank[rank]
        sent_row = sent_summary[rank]
        vela_x = [
            float(row[key]) for row in segment_rows for key in ("x0_um", "x1_um")
        ]
        vela_y = [
            float(row[key]) for row in segment_rows for key in ("y0_um", "y1_um")
        ]
        vela_max_field = max(float(row["electric_field_V_per_m"]) for row in segment_rows)
        vela_max_eparallel = max(
            float(row.get("electron_driving_field_V_per_m", 0.0))
            for row in segment_rows
        )
        vela_max_hparallel = max(
            float(row.get("hole_driving_field_V_per_m", 0.0))
            for row in segment_rows
        )
        vela_max_e_alpha = max(float(row["electron_alpha_m_inv"]) for row in segment_rows)
        vela_max_h_alpha = max(float(row["hole_alpha_m_inv"]) for row in segment_rows)
        rankwise_rows.append({
            "path_rank": rank,
            "vela_bias_V": vela_bias,
            "sentaurus_bias_V": args.sentaurus_bias,
            "vela_x_min_um": min(vela_x),
            "vela_x_max_um": max(vela_x),
            "vela_y_min_um": min(vela_y),
            "vela_y_max_um": max(vela_y),
            "sentaurus_x_min_um": sent_row["x_min_um"],
            "sentaurus_x_max_um": sent_row["x_max_um"],
            "sentaurus_y_min_um": sent_row["y_min_um"],
            "sentaurus_y_max_um": sent_row["y_max_um"],
            "vela_path_length_um": float(vela_row["path_length_m"]) * 1.0e6,
            "sentaurus_support_node_count": sent_row["support_node_count"],
            "vela_max_electric_field_V_per_m": vela_max_field,
            "sentaurus_max_electric_field_V_per_m": sent_row["max_electric_field_V_per_m"],
            "max_field_ratio_vela_over_sentaurus": safe_ratio(
                vela_max_field, float(sent_row["max_electric_field_V_per_m"])
            ),
            "vela_max_electron_eparallel_V_per_m": vela_max_eparallel,
            "sentaurus_max_electron_eparallel_V_per_m": sent_row["max_electron_eparallel_V_per_m"],
            "max_electron_eparallel_ratio_vela_over_sentaurus": (
                safe_ratio(vela_max_eparallel, float(sent_row["max_electron_eparallel_V_per_m"]))
                if sent_row["max_electron_eparallel_V_per_m"] is not None else None
            ),
            "vela_max_hole_eparallel_V_per_m": vela_max_hparallel,
            "sentaurus_max_hole_eparallel_V_per_m": sent_row["max_hole_eparallel_V_per_m"],
            "max_hole_eparallel_ratio_vela_over_sentaurus": (
                safe_ratio(vela_max_hparallel, float(sent_row["max_hole_eparallel_V_per_m"]))
                if sent_row["max_hole_eparallel_V_per_m"] is not None else None
            ),
            "vela_max_electron_alpha_m_inv": vela_max_e_alpha,
            "sentaurus_max_electron_alpha_m_inv": sent_row["max_electron_alpha_m_inv"],
            "max_electron_alpha_ratio_vela_over_sentaurus": safe_ratio(
                vela_max_e_alpha, float(sent_row["max_electron_alpha_m_inv"])
            ),
            "vela_max_hole_alpha_m_inv": vela_max_h_alpha,
            "sentaurus_max_hole_alpha_m_inv": sent_row["max_hole_alpha_m_inv"],
            "max_hole_alpha_ratio_vela_over_sentaurus": safe_ratio(
                vela_max_h_alpha, float(sent_row["max_hole_alpha_m_inv"])
            ),
            "vela_electron_ionization_integral": float(vela_row["electron_ionization_integral"]),
            "sentaurus_electron_ionization_integral": sent_row["electron_ionization_integral"],
            "vela_hole_ionization_integral": float(vela_row["hole_ionization_integral"]),
            "sentaurus_hole_ionization_integral": sent_row["hole_ionization_integral"],
            "vela_mean_ionization_integral": float(vela_row["mean_ionization_integral"]),
            "sentaurus_mean_ionization_integral": sent_row["mean_ionization_integral"],
            "mean_integral_ratio_vela_over_sentaurus": safe_ratio(
                float(vela_row["mean_ionization_integral"]),
                float(sent_row["mean_ionization_integral"]),
            ),
        })
    for vela_rank in ranks:
        sent_rank = mapping[vela_rank]
        support = sent_paths[sent_rank]
        sample_e_ratios: list[float] = []
        sample_h_ratios: list[float] = []
        sample_field_ratios: list[float] = []
        nearest_distances: list[float] = []
        vela_fields: list[float] = []
        sent_fields: list[float] = []
        vela_e_alphas: list[float] = []
        sent_e_alphas: list[float] = []
        vela_h_alphas: list[float] = []
        sent_h_alphas: list[float] = []
        for row in vela_segments[vela_rank]:
            midpoint = (
                0.5 * (float(row["x0_um"]) + float(row["x1_um"])),
                0.5 * (float(row["y0_um"]) + float(row["y1_um"])),
            )
            nearest = min(support, key=lambda node_id: distance(midpoint, nodes[node_id]))
            nearest_distance = distance(midpoint, nodes[nearest])
            vela_field = float(row["electric_field_V_per_m"])
            sent_field = efield_v_cm[nearest] * 100.0
            vela_e_alpha = float(row["electron_alpha_m_inv"])
            vela_h_alpha = float(row["hole_alpha_m_inv"])
            sent_e_alpha = e_alpha_cm[nearest] * 100.0
            sent_h_alpha = h_alpha_cm[nearest] * 100.0
            field_ratio = safe_ratio(vela_field, sent_field)
            e_ratio = safe_ratio(vela_e_alpha, sent_e_alpha)
            h_ratio = safe_ratio(vela_h_alpha, sent_h_alpha)
            sample_field_ratios.append(field_ratio)
            sample_e_ratios.append(e_ratio)
            sample_h_ratios.append(h_ratio)
            nearest_distances.append(nearest_distance)
            vela_fields.append(vela_field)
            sent_fields.append(sent_field)
            vela_e_alphas.append(vela_e_alpha)
            sent_e_alphas.append(sent_e_alpha)
            vela_h_alphas.append(vela_h_alpha)
            sent_h_alphas.append(sent_h_alpha)
            matched_segment_rows.append({
                "vela_path_rank": vela_rank,
                "sentaurus_path_rank": sent_rank,
                "segment_index": int(row["segment_index"]),
                "vela_mid_x_um": midpoint[0],
                "vela_mid_y_um": midpoint[1],
                "sentaurus_nearest_node": nearest,
                "sentaurus_x_um": nodes[nearest][0],
                "sentaurus_y_um": nodes[nearest][1],
                "nearest_distance_um": nearest_distance,
                "vela_electric_field_V_per_m": vela_field,
                "sentaurus_electric_field_V_per_m": sent_field,
                "vela_over_sentaurus_field": field_ratio,
                "vela_electron_alpha_m_inv": vela_e_alpha,
                "sentaurus_electron_alpha_m_inv": sent_e_alpha,
                "vela_over_sentaurus_electron_alpha": e_ratio,
                "vela_hole_alpha_m_inv": vela_h_alpha,
                "sentaurus_hole_alpha_m_inv": sent_h_alpha,
                "vela_over_sentaurus_hole_alpha": h_ratio,
            })
        vela_row = summary_by_rank[vela_rank]
        symmetric_mean, hausdorff, vela_rms = geometry_metrics[vela_rank]
        compare_rows.append({
            "vela_bias_V": vela_bias,
            "sentaurus_bias_V": args.sentaurus_bias,
            "vela_path_rank": vela_rank,
            "sentaurus_path_rank": sent_rank,
            "geometry_symmetric_mean_distance_um": symmetric_mean,
            "geometry_hausdorff_distance_um": hausdorff,
            "vela_midpoint_rms_nearest_distance_um": vela_rms,
            "vela_path_length_um": float(vela_row["path_length_m"]) * 1.0e6,
            "sentaurus_support_node_count": sent_summary[sent_rank]["support_node_count"],
            "median_finite_nearest_field_ratio_vela_over_sentaurus": finite_median(sample_field_ratios),
            "median_finite_nearest_electron_alpha_ratio_vela_over_sentaurus": finite_median(sample_e_ratios),
            "median_finite_nearest_hole_alpha_ratio_vela_over_sentaurus": finite_median(sample_h_ratios),
            "vela_max_electric_field_V_per_m": max(vela_fields),
            "sentaurus_max_matched_electric_field_V_per_m": max(sent_fields),
            "max_field_ratio_vela_over_sentaurus": safe_ratio(max(vela_fields), max(sent_fields)),
            "vela_max_electron_alpha_m_inv": max(vela_e_alphas),
            "sentaurus_max_matched_electron_alpha_m_inv": max(sent_e_alphas),
            "max_electron_alpha_ratio_vela_over_sentaurus": safe_ratio(max(vela_e_alphas), max(sent_e_alphas)),
            "vela_max_hole_alpha_m_inv": max(vela_h_alphas),
            "sentaurus_max_matched_hole_alpha_m_inv": max(sent_h_alphas),
            "max_hole_alpha_ratio_vela_over_sentaurus": safe_ratio(max(vela_h_alphas), max(sent_h_alphas)),
            "vela_electron_ionization_integral": float(vela_row["electron_ionization_integral"]),
            "sentaurus_electron_ionization_integral": sent_summary[sent_rank]["electron_ionization_integral"],
            "vela_hole_ionization_integral": float(vela_row["hole_ionization_integral"]),
            "sentaurus_hole_ionization_integral": sent_summary[sent_rank]["hole_ionization_integral"],
            "vela_mean_ionization_integral": float(vela_row["mean_ionization_integral"]),
            "sentaurus_mean_ionization_integral": sent_summary[sent_rank]["mean_ionization_integral"],
            "mean_integral_ratio_vela_over_sentaurus": safe_ratio(
                float(vela_row["mean_ionization_integral"]),
                float(sent_summary[sent_rank]["mean_ionization_integral"]),
            ),
        })

    write_rows(out / "sentaurus_top3_path_summary.csv", list(sent_summary.values()))
    write_rows(out / "rankwise_physics_compare.csv", rankwise_rows)
    write_rows(out / "vela_sentaurus_matched_segments.csv", matched_segment_rows)
    write_rows(out / "vela_sentaurus_path_summary_compare.csv", compare_rows)
    payload = {
        "method": {
            "sentaurus_path_extraction": (
                "distinct positive MeanIonIntegral node plateaus for geometry; "
                "final WriteAll carrier integrals for numeric ranking"
                if sent_log.is_file()
                else "distinct positive MeanIonIntegral node plateaus"
            ),
            "sentaurus_mean_definition": (
                "arithmetic mean of final WriteAll electron and hole path integrals"
                if sent_log.is_file()
                else "TDR MeanIonIntegral plateau fallback"
            ),
            "path_assignment": "minimum total symmetric nearest-neighbour distance",
            "alpha_and_field_sampling": "nearest Sentaurus path-support node to each Vela segment midpoint",
            "geometry_warning": (
                "Sentaurus coordinates are raw plateau support nodes, not an inferred centerline"
            ),
        },
        "vela_bias_V": vela_bias,
        "sentaurus_bias_V": args.sentaurus_bias,
        "assignment_vela_to_sentaurus": mapping,
        "sentaurus_positive_path_plateau_count": len(plateau_nodes),
        "comparison": compare_rows,
        "rankwise_physics_comparison": rankwise_rows,
    }
    (out / "summary.json").write_text(
        json.dumps(payload, indent=2, allow_nan=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(payload, indent=2, allow_nan=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
