#!/usr/bin/env python3
"""Independently verify the Minimal6 element-edge fixed-state replay."""

from __future__ import annotations

import argparse
import csv
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--root",
        type=Path,
        default=Path(
            "build-release/"
            "pn2d-minimal6-element-edge-gss-laux-fixed-state-20260725"
        ),
    )
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="ascii") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def dex(candidate: float, reference: float) -> float | None:
    if candidate <= 0.0 or reference <= 0.0:
        return None
    return abs(math.log10(candidate / reference))


def finite_summary(values: list[float | None]) -> dict[str, float | int]:
    finite = sorted(
        value for value in values
        if value is not None and math.isfinite(value)
    )
    if not finite:
        raise RuntimeError("summary has no finite values")
    return {
        "count": len(finite),
        "median_dex": statistics.median(finite),
        "max_dex": finite[-1],
    }


def main() -> int:
    root = parse_args().root.resolve()
    with (root / "manifest.json").open(encoding="ascii") as stream:
        manifest = json.load(stream)
    topologies = tuple(manifest["scope"]["topologies"])
    biases = tuple(int(value) for value in manifest["scope"]["biases_V"])
    rows = read_csv(root / "fixed_state_comparison.csv")
    failures: list[str] = []
    expected_row_count = len(topologies) * len(biases) * 12
    if len(rows) != expected_row_count:
        failures.append(
            f"expected {expected_row_count} element-vertex rows, got {len(rows)}"
        )

    states = defaultdict(list)
    for row in rows:
        states[(row["topology"], int(row["bias_V"]))].append(row)
    expected_states = {
        (topology, bias)
        for topology in topologies
        for bias in biases
    }
    if set(states) != expected_states:
        failures.append("state lattice does not match manifest scope")
    if any(len(state_rows) != 12 for state_rows in states.values()):
        failures.append("a state does not have exactly 12 element-vertex rows")

    zero_partial_volume_count = sum(
        float(row["edge_partial_volume_m2"]) == 0.0 for row in rows
    )
    expected_zero_partial_volumes = 4 * len(expected_states)
    if zero_partial_volume_count != expected_zero_partial_volumes:
        failures.append(
            f"expected {expected_zero_partial_volumes} zero diagonal "
            f"partial volumes, got "
            f"{zero_partial_volume_count}"
        )

    metrics: dict[str, dict[str, float | int]] = {}
    for carrier in ("electron", "hole"):
        for quantity in ("edge", "vector", "alpha"):
            metrics[f"{carrier}_{quantity}"] = finite_summary(
                [
                    float(row[f"{carrier}_{quantity}_error_dex"])
                    if row[f"{carrier}_{quantity}_error_dex"]
                    else None
                    for row in rows
                ]
            )

    source_identity_max_relative = 0.0
    node_buckets: dict[
        tuple[str, int, int, str], list[float]
    ] = defaultdict(lambda: [0.0, 0.0])
    for row in rows:
        topology = row["topology"]
        bias = int(row["bias_V"])
        node = int(row["node_id"])
        measure_m2 = float(row["vertex_measure_m2"])
        for carrier in ("electron", "hole"):
            candidate = float(row[f"{carrier}_cpp_qg_A_um"])
            expected = (
                float(row[f"{carrier}_cpp_alpha_cm_inv"])
                * float(row[f"{carrier}_cpp_vector_A_cm2"])
                * measure_m2
            )
            scale = max(abs(candidate), abs(expected), 1.0e-300)
            source_identity_max_relative = max(
                source_identity_max_relative,
                abs(candidate - expected) / scale,
            )
            bucket = node_buckets[(topology, bias, node, carrier)]
            bucket[0] += candidate
            bucket[1] += float(
                row[f"{carrier}_sentaurus_qg_A_um"]
            )
    if source_identity_max_relative > 5.0e-14:
        failures.append(
            "element source identity relative error exceeds 5e-14"
        )

    node_rows: list[dict[str, object]] = []
    for (topology, bias, node, carrier), (candidate, reference) in sorted(
        node_buckets.items()
    ):
        node_rows.append(
            {
                "topology": topology,
                "bias_V": bias,
                "node_id": node,
                "carrier": carrier,
                "cpp_accumulated_qg_A_um": format(candidate, ".17g"),
                "sentaurus_accumulated_qg_A_um": format(reference, ".17g"),
                "absolute_error_dex": (
                    "" if dex(candidate, reference) is None
                    else format(dex(candidate, reference), ".17g")
                ),
            }
        )
    write_csv(root / "node_source_comparison.csv", node_rows)

    state_rows: list[dict[str, object]] = []
    for topology, bias in sorted(states):
        for carrier in ("electron", "hole"):
            selected = [
                row for row in node_rows
                if row["topology"] == topology
                and row["bias_V"] == bias
                and row["carrier"] == carrier
            ]
            candidate = sum(
                float(row["cpp_accumulated_qg_A_um"])
                for row in selected
            )
            reference = sum(
                float(row["sentaurus_accumulated_qg_A_um"])
                for row in selected
            )
            state_rows.append(
                {
                    "topology": topology,
                    "bias_V": bias,
                    "carrier": carrier,
                    "cpp_integral_A_um": format(candidate, ".17g"),
                    "sentaurus_integral_A_um": format(reference, ".17g"),
                    "absolute_error_dex": format(
                        dex(candidate, reference), ".17g"
                    ),
                }
            )
    write_csv(root / "state_source_summary.csv", state_rows)

    for carrier in ("electron", "hole"):
        metrics[f"{carrier}_node_source"] = finite_summary(
            [
                float(row["absolute_error_dex"])
                for row in node_rows
                if row["carrier"] == carrier and row["absolute_error_dex"]
            ]
        )
        metrics[f"{carrier}_integral_source"] = finite_summary(
            [
                float(row["absolute_error_dex"])
                for row in state_rows
                if row["carrier"] == carrier
            ]
        )

    gates = {
        "electron_alpha_max_dex": (metrics["electron_alpha"]["max_dex"], 1e-7),
        "hole_alpha_max_dex": (metrics["hole_alpha"]["max_dex"], 1e-7),
        "electron_edge_median_dex": (
            metrics["electron_edge"]["median_dex"], 0.1
        ),
        "hole_edge_median_dex": (
            metrics["hole_edge"]["median_dex"], 0.1
        ),
        "electron_vector_median_dex": (
            metrics["electron_vector"]["median_dex"], 0.1
        ),
        "hole_vector_median_dex": (
            metrics["hole_vector"]["median_dex"], 0.1
        ),
        "electron_node_source_max_dex": (
            metrics["electron_node_source"]["max_dex"], 0.25
        ),
        "hole_node_source_max_dex": (
            metrics["hole_node_source"]["max_dex"], 0.1
        ),
        "electron_integral_source_max_dex": (
            metrics["electron_integral_source"]["max_dex"], 0.21
        ),
        "hole_integral_source_max_dex": (
            metrics["hole_integral_source"]["max_dex"], 0.08
        ),
    }
    for name, (value, limit) in gates.items():
        if float(value) > limit:
            failures.append(f"{name}={value} exceeds {limit}")

    result = {
        "schema_version": 1,
        "status": "passed" if not failures else "failed",
        "failures": failures,
        "row_count": len(rows),
        "state_count": len(states),
        "zero_partial_volume_count": zero_partial_volume_count,
        "source_identity_max_relative_error": source_identity_max_relative,
        "metrics": metrics,
        "gates": {
            name: {"value": value, "limit": limit}
            for name, (value, limit) in gates.items()
        },
    }
    with (root / "independent_verification.json").open(
        "w", encoding="ascii", newline="\n"
    ) as stream:
        json.dump(result, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0 if not failures else 1


if __name__ == "__main__":
    raise SystemExit(main())
