#!/usr/bin/env python3
"""Parse and classify general-Tri3 Sentaurus avalanche control evidence."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.pn2d_general_tri3_contract import (
    EXACT_BIASES_V,
    SCHEMA_ID,
    SENTAURUS_RELEASE,
    validate_source_manifests,
)


PREFIXES = {
    "AVAL_PROBE_BEGIN": "begins",
    "AVAL_PROBE_VERTEX": "vertices",
    "AVAL_PROBE_ELEMENT": "elements",
    "AVAL_PROBE_MEASURE": "measures",
    "AVAL_PROBE_EDGE": "edges",
    "AVAL_PROBE_INTEGRAL": "integrals",
    "AVAL_PROBE_END": "ends",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_scalar(value: str) -> Any:
    if re.fullmatch(r"[-+]?\d+", value):
        return int(value)
    try:
        parsed = float(value)
    except ValueError:
        return value
    return parsed if math.isfinite(parsed) else value


def parse_probe_line(line: str) -> tuple[str, dict[str, Any]] | None:
    for prefix, group in PREFIXES.items():
        if line.startswith(prefix + " "):
            row: dict[str, Any] = {}
            for token in line[len(prefix) + 1 :].split():
                if "=" not in token:
                    continue
                key, value = token.split("=", 1)
                row[key] = parse_scalar(value)
            return group, row
    return None


def parse_log(path: Path) -> dict[str, list[dict[str, Any]]]:
    groups: dict[str, list[dict[str, Any]]] = {
        value: [] for value in PREFIXES.values()
    }
    with path.open(encoding="ascii", errors="strict") as handle:
        for raw_line in handle:
            parsed = parse_probe_line(raw_line.strip())
            if parsed is not None:
                group, row = parsed
                groups[group].append(row)
    if not groups["begins"]:
        raise ValueError(f"no avalanche probe records in {path}")
    return groups


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    columns: list[str] = []
    for row in rows:
        for key in row:
            if key not in columns:
                columns.append(key)
    with path.open("w", newline="", encoding="ascii") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, value: dict[str, Any]) -> None:
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n",
        encoding="ascii",
        newline="\n",
    )


def triangle_angles(points: list[tuple[float, float]]) -> list[float]:
    if len(points) != 3:
        raise ValueError(f"Tri3 requires three unique points, got {points}")
    angles: list[float] = []
    for index in range(3):
        origin = points[index]
        left = points[(index + 1) % 3]
        right = points[(index + 2) % 3]
        ux, uy = left[0] - origin[0], left[1] - origin[1]
        vx, vy = right[0] - origin[0], right[1] - origin[1]
        denominator = math.hypot(ux, uy) * math.hypot(vx, vy)
        if denominator <= 0.0:
            raise ValueError("degenerate triangle edge")
        cosine = max(-1.0, min(1.0, (ux * vx + uy * vy) / denominator))
        angles.append(math.degrees(math.acos(cosine)))
    return sorted(angles)


def signed_area(points: list[tuple[float, float]]) -> float:
    if len(points) != 3:
        raise ValueError(f"Tri3 requires three ordered points, got {points}")
    return 0.5 * sum(
        points[index][0] * points[(index + 1) % 3][1]
        - points[(index + 1) % 3][0] * points[index][1]
        for index in range(3)
    )


def geometry_rows(
    groups: dict[str, list[dict[str, Any]]],
    bias: float,
) -> list[dict[str, Any]]:
    selected = [
        row for row in groups["edges"] if float(row["bias_V"]) == bias
    ]
    if not selected:
        raise ValueError(f"no edge rows at bias {bias}")
    vertices = {
        int(row["vertex"]): (float(row["x_um"]), float(row["y_um"]))
        for row in groups["vertices"]
        if float(row["bias_V"]) == bias
    }
    if not vertices:
        raise ValueError(f"no vertex rows at bias {bias}")
    x_values = [point[0] for point in vertices.values()]
    xmin, xmax = min(x_values), max(x_values)
    by_element: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for row in selected:
        by_element[int(row["element"])].append(row)
    local_vertices: dict[int, list[tuple[int, int]]] = defaultdict(list)
    for row in groups["measures"]:
        if float(row["bias_V"]) == bias:
            local_vertices[int(row["element"])].append(
                (int(row["local_vertex"]), int(row["vertex"]))
            )

    output: list[dict[str, Any]] = []
    for element, edges in sorted(by_element.items()):
        permutation = [
            vertex
            for _, vertex in sorted(local_vertices[element])
        ]
        if len(permutation) != 3 or len(set(permutation)) != 3:
            raise ValueError(
                f"element {element}: invalid Tri3 permutation {permutation}"
            )
        points = [vertices[vertex] for vertex in permutation]
        angles = triangle_angles(points)
        area = signed_area(points)
        contact = any(
            abs(point[0] - xmin) <= 1.0e-12
            or abs(point[0] - xmax) <= 1.0e-12
            for point in points
        )
        max_angle = max(angles)
        if abs(max_angle - 90.0) <= 1.0e-8:
            angle_class = "right"
        elif max_angle < 90.0:
            angle_class = "acute"
        else:
            angle_class = "obtuse"
        lengths = sorted(
            math.hypot(
                points[(index + 1) % 3][0] - points[index][0],
                points[(index + 1) % 3][1] - points[index][1],
            )
            for index in range(3)
        )
        scalene = (
            abs(lengths[0] - lengths[1]) > 1.0e-10
            and abs(lengths[1] - lengths[2]) > 1.0e-10
        )
        kappas = [float(row["kappa"]) for row in edges]
        output.append(
            {
                "element": element,
                "cell_vertex_permutation": ";".join(
                    str(vertex) for vertex in permutation
                ),
                "node_coordinates": ";".join(
                    f"{vertex}:{point[0]:.17g}:{point[1]:.17g}"
                    for vertex, point in zip(permutation, points)
                ),
                "angle0_deg": angles[0],
                "angle1_deg": angles[1],
                "angle2_deg": angles[2],
                "angle_sum_error_deg": abs(sum(angles) - 180.0),
                "angle_class": angle_class,
                "scalene": int(scalene),
                "signed_area_um2": area,
                "orientation": "ccw" if area > 0.0 else "cw",
                "contact_adjacent": int(contact),
                "interior_element": int(not contact),
                "positive_support_count": sum(value > 1.0e-14 for value in kappas),
                "zero_support_count": sum(abs(value) <= 1.0e-14 for value in kappas),
                "negative_support_count": sum(value < -1.0e-14 for value in kappas),
                "kappa_values": ";".join(f"{value:.17g}" for value in kappas),
            }
        )
    return output


def state_count_rows(
    variant: str,
    groups: dict[str, list[dict[str, Any]]],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for begin in groups["begins"]:
        bias = float(begin["bias_V"])
        row: dict[str, Any] = {
            "variant": variant,
            "bias_V": bias,
            "declared_vertices": int(begin["vertices"]),
            "declared_edges": int(begin["edges"]),
            "declared_elements": int(begin["elements"]),
            "declared_element_vertices": int(begin["element_vertices"]),
        }
        for group in ("vertices", "elements", "measures", "edges", "integrals"):
            row[f"observed_{group}"] = sum(
                float(item["bias_V"]) == bias for item in groups[group]
            )
        rows.append(row)
    return rows


def integral_rows(
    variant: str,
    groups: dict[str, list[dict[str, Any]]],
) -> list[dict[str, Any]]:
    return [
        {"variant": variant, **row}
        for row in groups["integrals"]
    ]


def element_source_rows(
    variant: str,
    groups: dict[str, list[dict[str, Any]]],
) -> list[dict[str, Any]]:
    totals: dict[tuple[float, int], dict[str, float]] = defaultdict(
        lambda: {
            "qg_n_A_um": 0.0,
            "qg_p_A_um": 0.0,
            "qg_total_A_um": 0.0,
        }
    )
    for row in groups["measures"]:
        key = (float(row["bias_V"]), int(row["element"]))
        for column in ("qg_n_A_um", "qg_p_A_um", "qg_total_A_um"):
            totals[key][column] += float(row[column])
    return [
        {"variant": variant, "bias_V": bias, "element": element, **values}
        for (bias, element), values in sorted(totals.items())
    ]


def driver_element_comparison_rows(
    element_sources: list[dict[str, Any]],
    geometry: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    by_key = {
        (str(row["variant"]), float(row["bias_V"]), int(row["element"])): row
        for row in element_sources
    }
    contact = {
        int(row["element"]): bool(row["contact_adjacent"])
        for row in geometry
    }
    references = (
        ("implicit_default", "explicit_grad_qf"),
        ("explicit_electric_field", "explicit_grad_qf"),
        ("grad_qf_use_qf_contacts", "explicit_grad_qf"),
        ("electric_field_use_qf_contacts", "grad_qf_use_qf_contacts"),
        (
            "lowfield_mobility_avalanche_grad_qf",
            "lowfield_mobility_avalanche_electric_field",
        ),
        ("grad_qf_aval_dens_grad_qf", "explicit_grad_qf"),
    )
    rows: list[dict[str, Any]] = []
    for candidate, reference in references:
        for bias in EXACT_BIASES_V:
            for element in sorted(contact):
                cand = by_key[(candidate, bias, element)]
                ref = by_key[(reference, bias, element)]
                row: dict[str, Any] = {
                    "candidate": candidate,
                    "reference": reference,
                    "bias_V": bias,
                    "element": element,
                    "element_class": "contact" if contact[element] else "interior",
                }
                for carrier, column in (
                    ("electron", "qg_n_A_um"),
                    ("hole", "qg_p_A_um"),
                    ("total", "qg_total_A_um"),
                ):
                    candidate_value = float(cand[column])
                    reference_value = float(ref[column])
                    row[f"{carrier}_candidate_A_um"] = candidate_value
                    row[f"{carrier}_reference_A_um"] = reference_value
                    row[f"{carrier}_absolute_difference_A_um"] = abs(
                        candidate_value - reference_value
                    )
                    value = abs_dex(candidate_value, reference_value)
                    row[f"{carrier}_absolute_error_dex"] = (
                        "" if value is None else value
                    )
                rows.append(row)
    return rows


def driver_contact_class_summary_rows(
    rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, str, float, str], list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        grouped[(str(row["candidate"]), str(row["reference"]), float(row["bias_V"]), str(row["element_class"]))].append(row)
    output: list[dict[str, Any]] = []
    for key, group in sorted(grouped.items()):
        candidate, reference, bias, element_class = key
        row: dict[str, Any] = {
            "candidate": candidate,
            "reference": reference,
            "bias_V": bias,
            "element_class": element_class,
            "element_count": len(group),
        }
        for carrier in ("electron", "hole", "total"):
            dex_values = [float(item[f"{carrier}_absolute_error_dex"]) for item in group if item[f"{carrier}_absolute_error_dex"] != ""]
            absolute_values = [float(item[f"{carrier}_absolute_difference_A_um"]) for item in group]
            row[f"{carrier}_finite_dex_count"] = len(dex_values)
            row[f"{carrier}_median_error_dex"] = "" if not dex_values else statistics.median(dex_values)
            row[f"{carrier}_maximum_error_dex"] = "" if not dex_values else max(dex_values)
            row[f"{carrier}_maximum_absolute_difference_A_um"] = max(absolute_values)
        output.append(row)
    return output

def validate_counts(row: dict[str, Any]) -> None:
    expected = {
        "observed_vertices": row["declared_vertices"],
        "observed_elements": row["declared_elements"],
        "observed_measures": row["declared_element_vertices"],
        "observed_edges": 3 * row["declared_elements"],
        "observed_integrals": 1,
    }
    for key, value in expected.items():
        if int(row[key]) != int(value):
            raise ValueError(
                f"{row['variant']} {row['bias_V']}: "
                f"{key}={row[key]}, expected {value}"
            )


def abs_dex(candidate: float, reference: float) -> float | None:
    if candidate == 0.0 or reference == 0.0:
        return None
    return abs(math.log10(abs(candidate) / abs(reference)))


def driver_comparison_rows(
    integrals: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    by_key = {
        (str(row["variant"]), float(row["bias_V"])): row
        for row in integrals
    }
    references = (
        ("implicit_default", "explicit_grad_qf"),
        ("explicit_electric_field", "explicit_grad_qf"),
        ("grad_qf_use_qf_contacts", "explicit_grad_qf"),
        ("electric_field_use_qf_contacts", "grad_qf_use_qf_contacts"),
        (
            "lowfield_mobility_avalanche_grad_qf",
            "lowfield_mobility_avalanche_electric_field",
        ),
        ("grad_qf_aval_dens_grad_qf", "explicit_grad_qf"),
    )
    rows: list[dict[str, Any]] = []
    for candidate, reference in references:
        for bias in EXACT_BIASES_V:
            cand = by_key[(candidate, bias)]
            ref = by_key[(reference, bias)]
            row: dict[str, Any] = {
                "candidate": candidate,
                "reference": reference,
                "bias_V": bias,
            }
            for carrier, column in (
                ("electron", "qg_n_A_um"),
                ("hole", "qg_p_A_um"),
                ("total", "qg_total_A_um"),
            ):
                candidate_value = float(cand[column])
                reference_value = float(ref[column])
                row[f"{carrier}_candidate_A_um"] = candidate_value
                row[f"{carrier}_reference_A_um"] = reference_value
                value = abs_dex(candidate_value, reference_value)
                row[f"{carrier}_absolute_error_dex"] = (
                    "" if value is None else value
                )
            rows.append(row)
    return rows


def normalized_manifest(
    raw_manifest: dict[str, Any],
    case_name: str,
) -> dict[str, dict[str, Any]]:
    return {
        variant: {
            "schema": raw_manifest["schema"],
            "status": result["status"],
            "sentaurus_release": raw_manifest["sentaurus_release"],
            "exact_biases_V": raw_manifest["exact_biases_V"],
            "case_hashes": raw_manifest["case_hashes"],
        }
        for variant, result in raw_manifest["cases"][case_name].items()
    }


def main() -> int:
    args = parse_args()
    raw_root = args.raw_root.resolve()
    output_root = args.output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    raw_manifest_path = raw_root / "manifest.json"
    raw_manifest = json.loads(raw_manifest_path.read_text(encoding="ascii"))
    if raw_manifest.get("schema") != SCHEMA_ID:
        raise ValueError("raw schema mismatch")
    if raw_manifest.get("sentaurus_release") != SENTAURUS_RELEASE:
        raise ValueError("raw Sentaurus release mismatch")
    if tuple(raw_manifest.get("exact_biases_V", ())) != EXACT_BIASES_V:
        raise ValueError("raw bias matrix mismatch")
    case_names = tuple(raw_manifest.get("cases", {}))
    if len(case_names) != 1:
        raise ValueError(f"expected one case, got {case_names}")
    case_name = case_names[0]
    validate_source_manifests(
        normalized_manifest(raw_manifest, case_name)
    )

    parsed: dict[str, dict[str, list[dict[str, Any]]]] = {}
    counts: list[dict[str, Any]] = []
    integrals: list[dict[str, Any]] = []
    element_sources: list[dict[str, Any]] = []
    for variant, result in raw_manifest["cases"][case_name].items():
        if result["status"] != "passed":
            raise ValueError(f"{variant}: raw status is not passed")
        path = (
            raw_root
            / case_name
            / variant
            / "fetched"
            / f"run_{variant}.out"
        )
        groups = parse_log(path)
        parsed[variant] = groups
        variant_counts = state_count_rows(variant, groups)
        for row in variant_counts:
            validate_counts(row)
        counts.extend(variant_counts)
        integrals.extend(integral_rows(variant, groups))
        element_sources.extend(element_source_rows(variant, groups))

    geometry = geometry_rows(
        parsed["implicit_default"],
        EXACT_BIASES_V[0],
    )
    classes = Counter(row["angle_class"] for row in geometry)
    contact_counts = Counter(
        "contact" if row["contact_adjacent"] else "interior"
        for row in geometry
    )
    comparisons = driver_comparison_rows(integrals)
    element_comparisons = driver_element_comparison_rows(
        element_sources, geometry
    )
    contact_summaries = driver_contact_class_summary_rows(
        element_comparisons
    )

    write_csv(output_root / "mesh_classification.csv", geometry)
    write_csv(output_root / "state_counts.csv", counts)
    write_csv(output_root / "source_integrals.csv", integrals)
    write_csv(output_root / "driver_control_summary.csv", comparisons)
    write_csv(
        output_root / "driver_element_summary.csv", element_comparisons
    )
    write_csv(
        output_root / "driver_contact_class_summary.csv", contact_summaries
    )
    analysis_manifest = {
        "schema": SCHEMA_ID,
        "status": "valid",
        "case_name": case_name,
        "source_manifest": str(raw_manifest_path),
        "source_manifest_sha256": sha256(raw_manifest_path),
        "sentaurus_release": raw_manifest["sentaurus_release"],
        "exact_biases_V": raw_manifest["exact_biases_V"],
        "variant_count": len(parsed),
        "state_count": len(counts),
        "element_count": len(geometry),
        "angle_class_counts": dict(sorted(classes.items())),
        "contact_class_counts": dict(sorted(contact_counts.items())),
        "has_interior_element": contact_counts["interior"] > 0,
        "has_acute_scalene_element": any(
            row["angle_class"] == "acute"
            and row["scalene"]
            and row["positive_support_count"] == 3
            for row in geometry
        ),
        "has_obtuse_element": classes["obtuse"] > 0,
        "outputs": {
            name: sha256(output_root / name)
            for name in (
                "mesh_classification.csv",
                "state_counts.csv",
                "source_integrals.csv",
                "driver_control_summary.csv",
                "driver_element_summary.csv",
                "driver_contact_class_summary.csv",
            )
        },
    }
    write_json(output_root / "analysis_manifest.json", analysis_manifest)
    print(json.dumps(analysis_manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
