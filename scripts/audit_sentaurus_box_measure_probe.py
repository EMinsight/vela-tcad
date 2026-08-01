#!/usr/bin/env python3
"""Audit Sentaurus box measures against geometric and Vela local-volume rules."""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import defaultdict
from pathlib import Path


MEASURE_LINE = re.compile(
    r"^\s*(?P<grd>\d+)\s+(?P<des>-?\d+)\s+(?P<type>\d+)\s+"
    r"(?P<values>[^#]+?)\s*$"
)
PROBE_LINE = re.compile(
    r"AVAL_PROBE_MEASURE .*?element=(?P<element>\d+) "
    r"local_vertex=(?P<local>\d+) vertex=(?P<vertex>\d+) "
    r"measure_um2=(?P<measure>[-+0-9.eE]+)"
)
REGION_SUMMARY = re.compile(
    r"^\s*R\.Si\s+(?P<volume>[-+0-9.eE]+)\s+"
    r"(?P<box>[-+0-9.eE]+)\s+(?P<delta>[-+0-9.eE]+)",
    re.MULTILINE,
)
SIGNED_INTERSECTION = re.compile(
    r"maximum CoeffIntersection\s*=\s*(?P<value>[-+0-9.eE]+)"
)


def parse_vertices(path: Path) -> list[tuple[float, float]]:
    text = path.read_text(encoding="utf-8")
    match = re.search(r"Vertices\s*\((\d+)\)\s*\{(.*?)\n\s*\}", text, re.DOTALL)
    if match is None:
        raise ValueError(f"Vertices block not found in {path}")
    count = int(match.group(1))
    vertices = [
        tuple(float(value) for value in line.split()[:2])
        for line in match.group(2).splitlines()
        if line.strip()
    ]
    if len(vertices) != count:
        raise ValueError(f"declared {count} vertices but parsed {len(vertices)}")
    return vertices  # type: ignore[return-value]


def parse_measures(path: Path) -> dict[int, list[float]]:
    text = path.read_text(encoding="utf-8")
    match = re.search(r"\n\s*Measure\s*\{.*?\n(.*?)\n\s*\}", text, re.DOTALL)
    if match is None:
        raise ValueError(f"Measure section not found in {path}")
    result = {}
    for line in match.group(1).splitlines():
        parsed = MEASURE_LINE.match(line)
        if parsed is None or int(parsed.group("type")) != 2:
            continue
        result[int(parsed.group("des"))] = [
            float(value) for value in parsed.group("values").split()
        ]
    return result


def parse_element_vertices(path: Path) -> dict[int, list[int]]:
    mapped: dict[int, dict[int, int]] = defaultdict(dict)
    for match in PROBE_LINE.finditer(path.read_text(encoding="utf-8")):
        mapped[int(match.group("element"))][int(match.group("local"))] = int(
            match.group("vertex")
        )
    return {
        element: [local[index] for index in sorted(local)]
        for element, local in mapped.items()
    }


def parse_dataset(path: Path, name: str) -> list[float]:
    text = path.read_text(encoding="utf-8")
    match = re.search(
        rf'Dataset \("{re.escape(name)}"\)\s*\{{.*?'
        r"Values\s*\(\d+\)\s*\{(?P<values>.*?)\n\s*\}\s*\n\s*\}",
        text,
        re.DOTALL,
    )
    if match is None:
        raise ValueError(f"dataset {name} or its Values block not found in {path}")
    return [float(value) for value in match.group("values").split()]


def geometry(points: list[tuple[float, float]]) -> tuple[float, list[float]]:
    area = 0.5 * abs(
        (points[1][0] - points[0][0]) * (points[2][1] - points[0][1])
        - (points[2][0] - points[0][0]) * (points[1][1] - points[0][1])
    )
    angles = []
    for index in range(3):
        center = points[index]
        left = points[(index + 1) % 3]
        right = points[(index + 2) % 3]
        u = (left[0] - center[0], left[1] - center[1])
        v = (right[0] - center[0], right[1] - center[1])
        denom = math.hypot(*u) * math.hypot(*v)
        cosine = max(-1.0, min(1.0, (u[0] * v[0] + u[1] * v[1]) / denom))
        angles.append(math.degrees(math.acos(cosine)))
    return area, angles


def raw_circumcentric_shares(points: list[tuple[float, float]]) -> list[float]:
    shares = []
    for index in range(3):
        j = (index + 1) % 3
        k = (index + 2) % 3

        def cotangent(a: int, b: int, opposite: int) -> float:
            u = (
                points[a][0] - points[opposite][0],
                points[a][1] - points[opposite][1],
            )
            v = (
                points[b][0] - points[opposite][0],
                points[b][1] - points[opposite][1],
            )
            return (u[0] * v[0] + u[1] * v[1]) / abs(u[0] * v[1] - u[1] * v[0])

        distance_ik2 = sum((points[index][axis] - points[k][axis]) ** 2 for axis in (0, 1))
        distance_ij2 = sum((points[index][axis] - points[j][axis]) ** 2 for axis in (0, 1))
        shares.append(
            0.125
            * (
                distance_ik2 * cotangent(index, k, j)
                + distance_ij2 * cotangent(index, j, k)
            )
        )
    return shares


def vela_mixed_shares(
    points: list[tuple[float, float]], area: float, angles: list[float]
) -> list[float]:
    max_index = max(range(3), key=angles.__getitem__)
    if angles[max_index] > 90.0:
        result = [0.25 * area] * 3
        result[max_index] = 0.5 * area
        return result
    return [max(0.0, value) for value in raw_circumcentric_shares(points)]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grid", type=Path, required=True)
    parser.add_argument("--average-debug", type=Path, required=True)
    parser.add_argument("--average-log", type=Path, required=True)
    parser.add_argument("--mix-debug", type=Path, required=True)
    parser.add_argument("--bm-data", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    vertices = parse_vertices(args.grid)
    average = parse_measures(args.average_debug)
    mixed = parse_measures(args.mix_debug)
    element_vertices = parse_element_vertices(args.average_log)
    bm_coeff = parse_dataset(args.bm_data, "BM_CoeffIntersectionNonDelaunayElements")
    rows = []
    node_totals: dict[str, dict[int, float]] = {
        key: defaultdict(float)
        for key in ("sentaurus_average", "sentaurus_mix", "vela_mixed")
    }
    for element in sorted(average):
        node_ids = element_vertices[element]
        points = [vertices[node] for node in node_ids]
        area, angles = geometry(points)
        raw = raw_circumcentric_shares(points)
        vela = vela_mixed_shares(points, area, angles)
        for local, node in enumerate(node_ids):
            for key, values in (
                ("sentaurus_average", average[element]),
                ("sentaurus_mix", mixed[element]),
                ("vela_mixed", vela),
            ):
                node_totals[key][node] += values[local]
            rows.append(
                {
                    "element": element,
                    "local_vertex": local,
                    "vertex": node,
                    "area_um2": area,
                    "angle_deg": angles[local],
                    "is_obtuse_element": max(angles) > 90.0 + 1e-10,
                    "sentaurus_average_measure_um2": average[element][local],
                    "sentaurus_mix_measure_um2": mixed[element][local],
                    "raw_circumcentric_measure_um2": raw[local],
                    "vela_mixed_measure_um2": vela[local],
                    "sent_mix_minus_vela_mixed_um2": mixed[element][local] - vela[local],
                    "bm_coeff_intersection": bm_coeff[element],
                }
            )

    log_text = args.average_log.read_text(encoding="utf-8")
    region = REGION_SUMMARY.search(log_text)
    signed = SIGNED_INTERSECTION.search(log_text)
    if region is None or signed is None:
        raise ValueError("box summary missing from AverageBox run log")
    obtuse_elements = sorted(
        {int(row["element"]) for row in rows if row["is_obtuse_element"]}
    )
    summary = {
        "schema": "vela.sentaurus_box_measure_audit.v1",
        "geometry_area_um2": sum(
            float(row["area_um2"])
            for row in rows
            if int(row["local_vertex"]) == 0
        ),
        "sentaurus_average_measure_sum_um2": sum(sum(values) for values in average.values()),
        "sentaurus_mix_measure_sum_um2": sum(sum(values) for values in mixed.values()),
        "sentaurus_average_negative_measure_count": sum(
            value < 0.0 for values in average.values() for value in values
        ),
        "sentaurus_mix_negative_measure_count": sum(
            value < 0.0 for values in mixed.values() for value in values
        ),
        "raw_circumcentric_negative_measure_count": sum(
            float(row["raw_circumcentric_measure_um2"]) < 0.0 for row in rows
        ),
        "obtuse_elements": obtuse_elements,
        "sentaurus_log_geometry_volume_um2": float(region.group("volume")),
        "sentaurus_log_average_box_volume_um2": float(region.group("box")),
        "sentaurus_log_average_delta_percent": float(region.group("delta")),
        "sentaurus_log_signed_coeff_intersection": float(signed.group("value")),
        "bm_plot_max_coeff_intersection": max(bm_coeff),
        "sent_mix_vs_vela_mixed_max_abs_local_um2": max(
            abs(float(row["sent_mix_minus_vela_mixed_um2"])) for row in rows
        ),
        "sent_mix_vs_vela_mixed_l1_node_um2": sum(
            abs(node_totals["sentaurus_mix"][node] - node_totals["vela_mixed"][node])
            for node in range(len(vertices))
        ),
    }

    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    with (output / "element_local_box_measure_compare.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    (output / "box_measure_audit.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
