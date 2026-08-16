#!/usr/bin/env python3
"""Identify the Sentaurus nodal-gradient reconstruction from exported fields."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

from audit_singledevice_eq231_box_formula0 import mixed_voronoi_shares, triangle_area


def read_scalar(path: Path) -> dict[int, float]:
    with path.open(newline="", encoding="utf-8") as stream:
        return {
            int(row["node_id"]): float(row["component0"])
            for row in csv.DictReader(stream)
        }


def read_vector(path: Path) -> dict[int, tuple[float, float]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return {
            int(row["node_id"]):
                (float(row["component0"]), float(row["component1"]))
            for row in csv.DictReader(stream)
        }


def field_entries(manifest: dict, name: str, region: int) -> list[dict]:
    entries = [
        entry for entry in manifest["fields"]
        if entry["name"] == name
        and (region < 0 or int(entry["region"]) == region)
    ]
    if not entries:
        raise ValueError(f"field {name!r} has no entries for region {region}")
    return entries


def merged_scalar(export_dir: Path, entries: list[dict]) -> dict[int, float]:
    result: dict[int, float] = {}
    for entry in entries:
        result.update(read_scalar(export_dir / "fields" / entry["csv_file"]))
    return result


def merged_vector(
    export_dir: Path, entries: list[dict]
) -> dict[int, tuple[float, float]]:
    result: dict[int, tuple[float, float]] = {}
    for entry in entries:
        result.update(read_vector(export_dir / "fields" / entry["csv_file"]))
    return result


def solve_2x2(
    a00: float, a01: float, a11: float, b0: float, b1: float
) -> tuple[float, float] | None:
    determinant = a00 * a11 - a01 * a01
    if abs(determinant) < 1.0e-30:
        return None
    return (
        (a11 * b0 - a01 * b1) / determinant,
        (-a01 * b0 + a00 * b1) / determinant,
    )


def summarize(
    predicted: dict[int, tuple[float, float]],
    reference: dict[int, tuple[float, float]],
) -> dict[str, float | int]:
    errors = []
    references = []
    for node in sorted(reference.keys() & predicted.keys()):
        px, py = predicted[node]
        rx, ry = reference[node]
        errors.append(math.hypot(px - rx, py - ry))
        references.append(math.hypot(rx, ry))
    squared = sum(error * error for error in errors)
    reference_squared = sum(value * value for value in references)
    return {
        "nodes": len(errors),
        "mean_absolute_error_V_per_um": sum(errors) / len(errors),
        "rms_error_V_per_um": math.sqrt(squared / len(errors)),
        "relative_l2": math.sqrt(squared / reference_squared),
        "maximum_error_V_per_um": max(errors),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--export-dir", type=Path, required=True)
    parser.add_argument(
        "--region", type=int, required=True,
        help="TDR region id; use -1 for a global all-region reconstruction",
    )
    parser.add_argument("--scalar", default="ElectrostaticPotential")
    parser.add_argument("--vector", default="ElectricField")
    parser.add_argument(
        "--reference-scale",
        type=float,
        default=-1.0e-4,
        help="convert vector CSV V/cm to the target scalar gradient V/um",
    )
    args = parser.parse_args()

    manifest = json.loads(
        (args.export_dir / "field_manifest.json").read_text(encoding="utf-8")
    )
    scalar_entries = field_entries(manifest, args.scalar, args.region)
    vector_entries = field_entries(manifest, args.vector, args.region)
    region_names = {entry["region_name"] for entry in scalar_entries}
    region_name = next(iter(region_names)) if args.region >= 0 else "<all>"
    scalar = merged_scalar(args.export_dir, scalar_entries)
    vector_native = merged_vector(args.export_dir, vector_entries)
    reference = {
        node: (
            args.reference_scale * value[0],
            args.reference_scale * value[1],
        )
        for node, value in vector_native.items()
    }

    with (args.export_dir / "nodes.csv").open(newline="", encoding="utf-8") as stream:
        coordinates = {
            int(row["id"]): {"x": float(row["x_um"]), "y": float(row["y_um"])}
            for row in csv.DictReader(stream)
        }
    with (args.export_dir / "elements.csv").open(newline="", encoding="utf-8") as stream:
        triangles = [
            (int(row["node0"]), int(row["node1"]), int(row["node2"]))
            for row in csv.DictReader(stream)
            if args.region < 0 or row["region"] == region_name
        ]

    gradient_sums = {
        name: defaultdict(lambda: [0.0, 0.0, 0.0])
        for name in ("unweighted_cell", "area_cell", "box_cell")
    }
    neighbours: dict[int, set[int]] = defaultdict(set)
    for triangle in triangles:
        points = [coordinates[node] for node in triangle]
        area = triangle_area(*points)
        x = [point["x"] for point in points]
        y = [point["y"] for point in points]
        b = [y[1] - y[2], y[2] - y[0], y[0] - y[1]]
        c = [x[2] - x[1], x[0] - x[2], x[1] - x[0]]
        gx = sum(scalar[triangle[i]] * b[i] for i in range(3)) / (2.0 * area)
        gy = sum(scalar[triangle[i]] * c[i] for i in range(3)) / (2.0 * area)
        shares = mixed_voronoi_shares(points)
        for local, node in enumerate(triangle):
            for scheme, weight in (
                ("unweighted_cell", 1.0),
                ("area_cell", area),
                ("box_cell", shares[local]),
            ):
                target = gradient_sums[scheme][node]
                target[0] += weight * gx
                target[1] += weight * gy
                target[2] += weight
            for other in triangle:
                if other != node:
                    neighbours[node].add(other)

    predictions = {
        name: {
            node: (value[0] / value[2], value[1] / value[2])
            for node, value in values.items()
            if value[2] > 0.0
        }
        for name, values in gradient_sums.items()
    }
    for exponent in (0, -1, -2):
        prediction = {}
        for node, adjacent in neighbours.items():
            x0 = coordinates[node]["x"]
            y0 = coordinates[node]["y"]
            f0 = scalar[node]
            a00 = a01 = a11 = b0 = b1 = 0.0
            for other in adjacent:
                dx = coordinates[other]["x"] - x0
                dy = coordinates[other]["y"] - y0
                distance = math.hypot(dx, dy)
                weight = distance ** exponent
                delta = scalar[other] - f0
                a00 += weight * dx * dx
                a01 += weight * dx * dy
                a11 += weight * dy * dy
                b0 += weight * dx * delta
                b1 += weight * dy * delta
            result = solve_2x2(a00, a01, a11, b0, b1)
            if result is not None:
                prediction[node] = result
        predictions[f"least_squares_length_pow_{exponent}"] = prediction

    print(json.dumps({
        "region": args.region,
        "region_name": region_name,
        "scalar": args.scalar,
        "reference_vector": args.vector,
        "schemes": {
            name: summarize(prediction, reference)
            for name, prediction in predictions.items()
        },
    }, indent=2))


if __name__ == "__main__":
    main()
