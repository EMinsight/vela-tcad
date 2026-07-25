#!/usr/bin/env python3
"""Close Minimal6 box current and avalanche integrals for selected states."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from collections import defaultdict
from pathlib import Path
from typing import Any


Q_LEGACY_C = 1.6021918e-19
SOURCE_INTEGRAL_TO_A_UM = Q_LEGACY_C * 1.0e-12
TOPOLOGIES = ("mirror", "sketch")
TARGET_BIASES = (-1.0, -10.0, -20.0)
CARRIERS = ("electron", "hole")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis", type=Path, required=True)
    parser.add_argument("--vela-factorization", type=Path, required=True)
    return parser.parse_args()


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="ascii") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, values: list[dict[str, Any]]) -> None:
    if not values:
        raise ValueError(f"refusing to write empty CSV {path}")
    with path.open("w", newline="", encoding="ascii") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(values[0]))
        writer.writeheader()
        writer.writerows(values)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def abs_dex(value: float, reference: float) -> float | None:
    if value <= 0.0 or reference <= 0.0:
        return None
    return abs(math.log10(value / reference))


def f(value: Any) -> float:
    return float(str(value))


def parse_tokens(line: str) -> dict[str, str]:
    return dict(re.findall(r"(\w+)=([^\s]+)", line))


def runtime_integrals(log: Path) -> dict[float, dict[str, float]]:
    result = {}
    for line in log.read_text(encoding="ascii").splitlines():
        if not line.startswith("AVAL_PROBE_INTEGRAL "):
            continue
        values = parse_tokens(line)
        result[float(values["bias_V"])] = {
            key: float(value)
            for key, value in values.items()
            if key != "bias_V"
        }
    if set(result) != set(TARGET_BIASES):
        raise ValueError(f"{log}: incomplete runtime integrals")
    return result


def main() -> int:
    args = parse_args()
    analysis = args.analysis.resolve()
    edge_rows = rows(analysis / "element_edges.csv")
    vertex_rows = rows(analysis / "vertices.csv")
    current_rows = rows(analysis / "currentplot_targets.csv")
    state_rows = rows(analysis / "state_source_summary.csv")
    current_by_state = {
        (row["topology"], f(row["bias_V"])): row for row in current_rows
    }
    vertices_by_state = defaultdict(dict)
    for row in vertex_rows:
        vertex = int(row["vertex"])
        if vertex < 6:
            vertices_by_state[(row["topology"], f(row["bias_V"]))][
                vertex
            ] = row

    closure_rows: list[dict[str, Any]] = []
    integral_rows: list[dict[str, Any]] = []
    corrected_state_rows: list[dict[str, Any]] = []

    for topology in TOPOLOGIES:
        log = (
            analysis.parent
            / "raw"
            / topology
            / "default"
            / "run_default.out"
        )
        runtime = runtime_integrals(log)
        for bias in TARGET_BIASES:
            state = (topology, bias)
            selected_edges = [
                row
                for row in edge_rows
                if row["topology"] == topology and f(row["bias_V"]) == bias
            ]
            vertices = vertices_by_state[state]
            current = current_by_state[state]
            min_x = min(f(row["x_um"]) for row in vertices.values())
            max_x = max(f(row["x_um"]) for row in vertices.values())
            contact_nodes = {
                "Anode": {
                    node
                    for node, row in vertices.items()
                    if abs(f(row["x_um"]) - min_x) <= 1.0e-12
                },
                "Cathode": {
                    node
                    for node, row in vertices.items()
                    if abs(f(row["x_um"]) - max_x) <= 1.0e-12
                },
            }
            endpoints = {}
            fluxes = {
                carrier: defaultdict(float) for carrier in CARRIERS
            }
            for row in selected_edges:
                edge = int(row["edge"])
                endpoints[edge] = (int(row["start"]), int(row["end"]))
                fluxes["electron"][edge] += f(row["box_flux_n_A_um"])
                fluxes["hole"][edge] += f(row["box_flux_p_A_um"])

            balances = {}
            for carrier in CARRIERS:
                balance = {node: 0.0 for node in vertices}
                for edge, flux in fluxes[carrier].items():
                    start, end = endpoints[edge]
                    balance[start] -= flux
                    balance[end] += flux
                balances[carrier] = balance
                label = "electron" if carrier == "electron" else "hole"
                for contact, nodes in contact_nodes.items():
                    predicted = sum(balance[node] for node in nodes)
                    reference = f(current[f"{contact}_{label}_A_um"])
                    closure_rows.append(
                        {
                            "topology": topology,
                            "bias_V": format(bias, ".17g"),
                            "carrier": carrier,
                            "location": contact,
                            "predicted_A_um": format(predicted, ".17g"),
                            "reference_A_um": format(reference, ".17g"),
                            "absolute_error_A_um": format(
                                abs(predicted - reference), ".17g"
                            ),
                            "relative_error": format(
                                abs(predicted - reference)
                                / max(abs(reference), 1.0e-300),
                                ".17g",
                            ),
                        }
                    )
            terminal_scale = max(
                abs(f(current["Anode_total_A_um"])),
                abs(f(current["Cathode_total_A_um"])),
                1.0e-300,
            )
            internal_nodes = (
                set(vertices)
                - contact_nodes["Anode"]
                - contact_nodes["Cathode"]
            )
            for node in sorted(internal_nodes):
                residual = sum(
                    balances[carrier][node] for carrier in CARRIERS
                )
                closure_rows.append(
                    {
                        "topology": topology,
                        "bias_V": format(bias, ".17g"),
                        "carrier": "total",
                        "location": f"internal_vertex_{node}",
                        "predicted_A_um": format(residual, ".17g"),
                        "reference_A_um": "0",
                        "absolute_error_A_um": format(
                            abs(residual), ".17g"
                        ),
                        "relative_error": format(
                            abs(residual) / terminal_scale, ".17g"
                        ),
                    }
                )

            for carrier, runtime_key, current_key in (
                (
                    "electron",
                    "qg_n_A_um",
                    "e_avalanche_integral_um2_cm3_s",
                ),
                (
                    "hole",
                    "qg_p_A_um",
                    "h_avalanche_integral_um2_cm3_s",
                ),
                (
                    "total",
                    "qg_total_A_um",
                    "total_avalanche_integral_um2_cm3_s",
                ),
            ):
                predicted = runtime[bias][runtime_key]
                reference = f(current[current_key]) * SOURCE_INTEGRAL_TO_A_UM
                integral_rows.append(
                    {
                        "topology": topology,
                        "bias_V": format(bias, ".17g"),
                        "carrier": carrier,
                        "readmeasure_qg_A_um": format(predicted, ".17g"),
                        "currentplot_qg_A_um": format(reference, ".17g"),
                        "absolute_error_A_um": format(
                            abs(predicted - reference), ".17g"
                        ),
                        "relative_error": format(
                            abs(predicted - reference)
                            / max(abs(reference), 1.0e-300),
                            ".17g",
                        ),
                    }
                )

    for row in state_rows:
        corrected = dict(row)
        if row["candidate"] == "vela_triangle_proxy_existing":
            predicted = f(row["predicted_integral_um2_cm3_s"]) / 1.0e-8
            reference = f(row["native_integral_um2_cm3_s"])
            corrected["predicted_integral_um2_cm3_s"] = format(
                predicted, ".17g"
            )
            error = abs_dex(predicted, reference)
            corrected["integral_absolute_error_dex"] = (
                "" if error is None else format(error, ".17g")
            )
        corrected_state_rows.append(corrected)

    closure_path = analysis / "box_current_closure.csv"
    integral_path = analysis / "source_integral_closure.csv"
    state_path = analysis / "state_source_summary_corrected.csv"
    write_csv(closure_path, closure_rows)
    write_csv(integral_path, integral_rows)
    write_csv(state_path, corrected_state_rows)

    manifest = json.loads(
        (analysis / "manifest.json").read_text(encoding="ascii")
    )
    manifest["status"] = "valid_diagnostic_replay_with_closure"
    manifest["output_sha256"].update(
        {
            closure_path.name: sha256(closure_path),
            integral_path.name: sha256(integral_path),
            state_path.name: sha256(state_path),
        }
    )
    manifest["closure"] = {
        "max_contact_relative_error": max(
            f(row["relative_error"])
            for row in closure_rows
            if row["location"] in {"Anode", "Cathode"}
        ),
        "max_internal_kcl_relative_error": max(
            f(row["relative_error"])
            for row in closure_rows
            if row["location"].startswith("internal_vertex_")
        ),
        "max_readmeasure_currentplot_relative_error": max(
            f(row["relative_error"]) for row in integral_rows
        ),
    }
    (analysis / "manifest_closed.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="ascii",
    )
    print(json.dumps(manifest["closure"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
