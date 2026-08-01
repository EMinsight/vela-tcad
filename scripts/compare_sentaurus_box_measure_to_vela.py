#!/usr/bin/env python3
"""Compare direct Sentaurus element-vertex Measure with Vela mixed Voronoi."""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

try:
    from scripts.audit_sentaurus_box_measure_probe import (
        geometry,
        parse_element_vertices,
        parse_measures,
        parse_vertices,
        vela_mixed_shares,
    )
except ModuleNotFoundError:
    from audit_sentaurus_box_measure_probe import (  # type: ignore[no-redef]
        geometry,
        parse_element_vertices,
        parse_measures,
        parse_vertices,
        vela_mixed_shares,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grid", type=Path, required=True)
    parser.add_argument("--debug", type=Path, required=True)
    parser.add_argument("--log", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    vertices = parse_vertices(args.grid)
    measures = parse_measures(args.debug)
    element_vertices = parse_element_vertices(args.log)
    sentaurus_nodes: dict[int, float] = defaultdict(float)
    vela_nodes: dict[int, float] = defaultdict(float)
    barycentric_nodes: dict[int, float] = defaultdict(float)
    local_differences = []
    max_angle = 0.0
    obtuse_elements = []
    geometry_area = 0.0
    for element in sorted(measures):
        node_ids = element_vertices[element]
        points = [vertices[node] for node in node_ids]
        area, angles = geometry(points)
        geometry_area += area
        max_angle = max(max_angle, max(angles))
        if max(angles) > 90.0 + 1e-10:
            obtuse_elements.append(element)
        vela = vela_mixed_shares(points, area, angles)
        for local, node in enumerate(node_ids):
            sentaurus = measures[element][local]
            sentaurus_nodes[node] += sentaurus
            vela_nodes[node] += vela[local]
            barycentric_nodes[node] += area / 3.0
            local_differences.append(sentaurus - vela[local])

    node_differences = [
        sentaurus_nodes[node] - vela_nodes[node] for node in range(len(vertices))
    ]
    barycentric_differences = [
        sentaurus_nodes[node] - barycentric_nodes[node] for node in range(len(vertices))
    ]
    summary = {
        "schema": "vela.sentaurus_box_measure_direct_compare.v1",
        "vertex_count": len(vertices),
        "triangle_count": len(measures),
        "geometry_area_um2": geometry_area,
        "sentaurus_measure_sum_um2": sum(sum(values) for values in measures.values()),
        "vela_mixed_measure_sum_um2": sum(vela_nodes.values()),
        "max_angle_deg": max_angle,
        "obtuse_element_count": len(obtuse_elements),
        "obtuse_elements": obtuse_elements,
        "sentaurus_negative_measure_count": sum(
            value < 0.0 for values in measures.values() for value in values
        ),
        "max_abs_element_local_difference_um2": max(map(abs, local_differences)),
        "max_abs_node_volume_difference_um2": max(map(abs, node_differences)),
        "l1_node_volume_difference_um2": sum(map(abs, node_differences)),
        "max_abs_barycentric_node_volume_difference_um2": max(
            map(abs, barycentric_differences)
        ),
        "l1_barycentric_node_volume_difference_um2": sum(
            map(abs, barycentric_differences)
        ),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
