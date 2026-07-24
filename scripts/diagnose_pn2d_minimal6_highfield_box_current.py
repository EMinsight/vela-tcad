#!/usr/bin/env python3
"""Replay Minimal6 box-edge current with high-field mobility branches."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from pathlib import Path

from pn2d_minimal6_diagnostics.highfield_box_replay import (
    BRANCH_ID,
    CELL_MAPPING,
    FIELD,
    RECONSTRUCTION_LABEL,
    coefficient_weighted_mobility,
    validate_sample_record,
)


BRANCH_FIELDS = {
    "sentaurus_native_final": "sentaurus_final_m2_per_Vs",
    BRANCH_ID: "electric_replay_m2_per_Vs",
    "sentaurus_lowfield_element_triangle_qfp": (
        "triangle_qf_replay_m2_per_Vs"
    ),
    "sentaurus_lowfield_element_native_qfp": (
        "native_qf_replay_m2_per_Vs"
    ),
    "sentaurus_lowfield_only": "sentaurus_low_field_m2_per_Vs",
}
PRODUCTION_BRANCH = "vela_imported_state_production_mobility"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, values: list[dict[str, object]]) -> None:
    if not values:
        raise ValueError(f"refusing to write empty CSV: {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(values[0]))
        writer.writeheader()
        writer.writerows(values)


def quantile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * fraction
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    return ordered[lower] + (
        ordered[upper] - ordered[lower]
    ) * (position - lower)


def contact_current(
    edge_values: dict[tuple[int, int], float],
    contact_nodes: set[int],
) -> float:
    outward = 0.0
    for (start, end), value in edge_values.items():
        start_in = start in contact_nodes
        end_in = end in contact_nodes
        if start_in == end_in:
            continue
        outward += value if start_in else -value
    return -outward


def node_divergence(
    edge_values: dict[tuple[int, int], float], node: int
) -> float:
    result = 0.0
    for (start, end), value in edge_values.items():
        if start == node:
            result += value
        elif end == node:
            result -= value
    return result


def run(args: argparse.Namespace) -> dict[str, object]:
    paths = {
        "highfield": args.highfield.resolve(),
        "adjacent": args.adjacent.resolve(),
        "stage_edges": args.stage_edges.resolve(),
        "terminal": args.terminal.resolve(),
    }
    highfield = {
        (
            row["topology"],
            float(row["bias_V"]),
            int(row["cell_id"]),
            row["carrier"],
        ): row
        for row in rows(paths["highfield"])
    }
    adjacent_rows = rows(paths["adjacent"])
    final_stage = "sentaurus_qfp_density_element_mobility"
    stage_rows = rows(paths["stage_edges"])
    reference = {
        (
            row["topology"],
            float(row["bias_V"]),
            row["carrier"],
            int(row["edge_id"]),
        ): row
        for row in stage_rows
        if row["stage"] == final_stage
    }
    production = {
        (
            row["topology"],
            float(row["bias_V"]),
            row["carrier"],
            int(row["edge_id"]),
        ): row
        for row in stage_rows
        if row["stage"] == "sentaurus_qfp_density"
    }
    terminal_rows = rows(paths["terminal"])
    terminal = {
        (
            row["topology"],
            float(row["bias_V"]),
            row["contact"],
            row["carrier"],
        ): float(row["sentaurus_A_per_um"])
        for row in terminal_rows
    }
    if len(highfield) != 320 or len(reference) != 720:
        raise ValueError("input lattice does not match 40 states")

    highfield_manifest = json.loads(
        (paths["highfield"].parent / "manifest.json").read_text(
            encoding="utf-8"
        )
    )
    low_hashes = highfield_manifest["low_field_members"]
    electric_hash = highfield_manifest["inputs"]["transport"]["sha256"]

    adjacent: dict[
        tuple[str, float, str, int], list[dict[str, str]]
    ] = {}
    element_rows: list[dict[str, object]] = []
    for row in adjacent_rows:
        key = (
            row["topology"],
            float(row["bias_V"]),
            row["carrier"],
            int(row["edge_id"]),
        )
        adjacent.setdefault(key, []).append(row)
        cell = int(row["cell_id"])
        field_row = highfield[
            (key[0], key[1], cell, key[2])
        ]
        low_name = (
            "eMobility_region0_cells.csv"
            if key[2] == "electron"
            else "hMobility_region0_cells.csv"
        )
        for branch, field in BRANCH_FIELDS.items():
            sample = {
                "topology": key[0],
                "bias_V": key[1],
                "carrier": key[2],
                "edge_id": key[3],
                "node0": int(row["node0"]),
                "node1": int(row["node1"]),
                "vela_triangle_id": cell,
                "sentaurus_region_cell_id": CELL_MAPPING[key[0]][cell],
                "branch": branch,
                "mobility_m2_per_Vs": float(field_row[field]),
                "low_field_source_sha256": low_hashes[
                    f"{key[0]}/{low_name}"
                ]["sha256"],
                "electric_field_source_sha256": electric_hash,
                "saturation_velocity_m_per_s": FIELD[key[2]][
                    "saturation_velocity"
                ],
                "field_beta": FIELD[key[2]]["beta"],
                "kappa": float(row["kappa"]),
                "status": (
                    "geometric_zero"
                    if float(row["kappa"]) == 0.0
                    else "valid"
                ),
                "reconstruction_label": RECONSTRUCTION_LABEL,
            }
            validate_sample_record(sample)
            element_rows.append(sample)

    element_index = {
        (
            row["topology"],
            float(row["bias_V"]),
            row["carrier"],
            int(row["edge_id"]),
            int(row["vela_triangle_id"]),
            row["branch"],
        ): row
        for row in element_rows
    }
    edge_rows: list[dict[str, object]] = []
    edge_values_by_state: dict[
        tuple[str, float, str, str],
        dict[tuple[int, int], float],
    ] = {}
    for key, local_rows in adjacent.items():
        ref_row = reference[key]
        ref_current = float(ref_row["reference_A_per_um"])
        node0 = int(ref_row["node0"])
        node1 = int(ref_row["node1"])
        ref_mobility = None
        branch_mobility: dict[str, float | None] = {}
        branch_status: dict[str, str] = {}
        for branch in BRANCH_FIELDS:
            weighted = coefficient_weighted_mobility(
                [
                    (
                        float(local["kappa"]),
                        float(
                            element_index[
                                (
                                    key[0],
                                    key[1],
                                    key[2],
                                    key[3],
                                    int(local["cell_id"]),
                                    branch,
                                )
                            ]["mobility_m2_per_Vs"]
                        ),
                    )
                    for local in local_rows
                ]
            )
            branch_status[branch] = str(weighted["status"])
            value = weighted["mobility_m2_per_Vs"]
            branch_mobility[branch] = (
                None if value is None else float(value)
            )
            if branch == "sentaurus_native_final":
                ref_mobility = branch_mobility[branch]
        for branch in BRANCH_FIELDS:
            mobility = branch_mobility[branch]
            if ref_mobility is None and mobility is None:
                candidate = 0.0
                status = "exact_zero"
            elif ref_mobility is None or mobility is None:
                candidate = 0.0
                status = "support_mismatch"
            else:
                candidate = ref_current * mobility / ref_mobility
                if ref_current == 0.0 and candidate == 0.0:
                    status = "exact_zero"
                elif ref_current == 0.0:
                    status = "reference_zero_candidate_nonzero"
                elif candidate == 0.0:
                    status = "candidate_zero"
                else:
                    status = "valid"
            error = (
                abs(math.log10(abs(candidate) / abs(ref_current)))
                if status == "valid"
                else None
            )
            sign = (
                float(
                    math.copysign(1.0, candidate)
                    == math.copysign(1.0, ref_current)
                )
                if status == "valid"
                else None
            )
            edge_rows.append(
                {
                    "topology": key[0],
                    "bias_V": key[1],
                    "carrier": key[2],
                    "edge_id": key[3],
                    "node0": node0,
                    "node1": node1,
                    "branch": branch,
                    "reference_mobility_m2_per_Vs": ref_mobility,
                    "candidate_mobility_m2_per_Vs": mobility,
                    "reference_A_per_um": ref_current,
                    "candidate_A_per_um": candidate,
                    "absolute_log10_error_dex": error,
                    "sign_agreement": sign,
                    "status": status,
                    "reconstruction_label": RECONSTRUCTION_LABEL,
                }
            )
            edge_values_by_state.setdefault(
                (key[0], key[1], branch, key[2]), {}
            )[(node0, node1)] = candidate

        production_row = production[key]
        edge_rows.append(
            {
                "topology": key[0],
                "bias_V": key[1],
                "carrier": key[2],
                "edge_id": key[3],
                "node0": node0,
                "node1": node1,
                "branch": PRODUCTION_BRANCH,
                "reference_mobility_m2_per_Vs": ref_mobility,
                "candidate_mobility_m2_per_Vs": None,
                "reference_A_per_um": ref_current,
                "candidate_A_per_um": float(
                    production_row["candidate_A_per_um"]
                ),
                "absolute_log10_error_dex": (
                    None
                    if production_row["absolute_log10_error_dex"] == ""
                    else float(
                        production_row["absolute_log10_error_dex"]
                    )
                ),
                "sign_agreement": (
                    None
                    if production_row["sign_agreement"] == ""
                    else float(production_row["sign_agreement"])
                ),
                "status": production_row["status"],
                "reconstruction_label": RECONSTRUCTION_LABEL,
            }
        )

    summary_rows: list[dict[str, object]] = []
    for branch in (*BRANCH_FIELDS, PRODUCTION_BRANCH):
        for carrier in ("electron", "hole"):
            selected = [
                row
                for row in edge_rows
                if row["branch"] == branch
                and row["carrier"] == carrier
                and row["status"] == "valid"
            ]
            errors = [
                float(row["absolute_log10_error_dex"])
                for row in selected
            ]
            signs = [float(row["sign_agreement"]) for row in selected]
            summary_rows.append(
                {
                    "branch": branch,
                    "carrier": carrier,
                    "valid_count": len(selected),
                    "median_abs_dex": statistics.median(errors),
                    "p95_abs_dex": quantile(errors, 0.95),
                    "maximum_abs_dex": max(errors),
                    "sign_agreement_fraction": statistics.mean(signs),
                }
            )

    terminal_output: list[dict[str, object]] = []
    kcl_output: list[dict[str, object]] = []
    contacts = {"Anode": {0, 4}, "Cathode": {2, 3}}
    for topology in ("mirror", "sketch"):
        for magnitude in range(1, 21):
            bias = -float(magnitude)
            observed_total = {
                contact: sum(
                    terminal[(topology, bias, contact, carrier)]
                    for carrier in ("electron", "hole")
                )
                for contact in contacts
            }
            terminal_scale = max(
                abs(observed_total["Anode"]),
                abs(observed_total["Cathode"]),
                1.0e-300,
            )
            for branch in BRANCH_FIELDS:
                candidate_by_carrier = {
                    carrier: contact_current(
                        edge_values_by_state[
                            (topology, bias, branch, carrier)
                        ],
                        contacts["Anode"],
                    )
                    for carrier in ("electron", "hole")
                }
                candidate_total = sum(candidate_by_carrier.values())
                observed = observed_total["Anode"]
                terminal_output.append(
                    {
                        "topology": topology,
                        "bias_V": bias,
                        "branch": branch,
                        "contact": "Anode",
                        "candidate_electron_A_per_um": (
                            candidate_by_carrier["electron"]
                        ),
                        "candidate_hole_A_per_um": (
                            candidate_by_carrier["hole"]
                        ),
                        "candidate_total_A_per_um": candidate_total,
                        "sentaurus_total_A_per_um": observed,
                        "relative_error": abs(candidate_total - observed)
                        / max(abs(observed), 1.0e-300),
                    }
                )
                total_edges = {
                    pair: (
                        edge_values_by_state[
                            (topology, bias, branch, "electron")
                        ][pair]
                        + edge_values_by_state[
                            (topology, bias, branch, "hole")
                        ][pair]
                    )
                    for pair in edge_values_by_state[
                        (topology, bias, branch, "electron")
                    ]
                }
                for node in (1, 5):
                    divergence = node_divergence(total_edges, node)
                    kcl_output.append(
                        {
                            "topology": topology,
                            "bias_V": bias,
                            "branch": branch,
                            "node": node,
                            "total_current_divergence_A_per_um": divergence,
                            "relative_to_terminal_total": abs(divergence)
                            / terminal_scale,
                        }
                    )

    electric_summary = {
        row["carrier"]: row
        for row in summary_rows
        if row["branch"] == BRANCH_ID
    }
    electric_terminal = [
        row for row in terminal_output if row["branch"] == BRANCH_ID
    ]
    electric_kcl = [
        row for row in kcl_output if row["branch"] == BRANCH_ID
    ]
    gates = {
        "edge_count": all(
            int(electric_summary[carrier]["valid_count"]) == 200
            for carrier in ("electron", "hole")
        ),
        "edge_median": all(
            float(electric_summary[carrier]["median_abs_dex"]) <= 0.01
            for carrier in ("electron", "hole")
        ),
        "edge_p95": all(
            float(electric_summary[carrier]["p95_abs_dex"]) <= 0.05
            for carrier in ("electron", "hole")
        ),
        "edge_sign": all(
            float(
                electric_summary[carrier]["sign_agreement_fraction"]
            )
            == 1.0
            for carrier in ("electron", "hole")
        ),
        "terminal": max(
            float(row["relative_error"]) for row in electric_terminal
        )
        <= 0.02,
        "kcl": max(
            float(row["relative_to_terminal_total"])
            for row in electric_kcl
        )
        <= 1.0e-8,
    }
    status = "valid" if all(gates.values()) else "bounded_gate_failure"

    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    outputs = {
        "element_edge_samples.csv": element_rows,
        "box_edge_samples.csv": edge_rows,
        "box_edge_summary.csv": summary_rows,
        "terminal_comparison.csv": terminal_output,
        "total_current_kcl.csv": kcl_output,
    }
    for name, values in outputs.items():
        write_csv(output / name, values)
    report = [
        "# PN2D Minimal6 high-field box-current replay",
        "",
        f"Status: `{status}`",
        "",
        "| Carrier | Median dex | P95 dex | Maximum dex | Sign |",
        "|---|---:|---:|---:|---:|",
    ]
    for carrier in ("electron", "hole"):
        row = electric_summary[carrier]
        report.append(
            f"| {carrier} | {float(row['median_abs_dex']):.6g} | "
            f"{float(row['p95_abs_dex']):.6g} | "
            f"{float(row['maximum_abs_dex']):.6g} | "
            f"{float(row['sign_agreement_fraction']):.6g} |"
        )
    report.extend(
        [
            "",
            f"Maximum terminal relative error: "
            f"`{max(float(row['relative_error']) for row in electric_terminal):.6e}`.",
            "",
            f"Maximum internal KCL relative value: "
            f"`{max(float(row['relative_to_terminal_total']) for row in electric_kcl):.6e}`.",
            "",
        ]
    )
    (output / "report.md").write_text(
        "\n".join(report), encoding="utf-8", newline="\n"
    )
    output_hashes = {
        name: sha256(output / name)
        for name in (*outputs, "report.md")
    }
    manifest = {
        "schema_version": 1,
        "status": status,
        "experiment": "pn2d_minimal6_highfield_box_current",
        "typed_branch": BRANCH_ID,
        "reconstruction_label": RECONSTRUCTION_LABEL,
        "state_count": 40,
        "active_carrier_edge_count": 400,
        "element_edge_branch_sample_count": len(element_rows),
        "gates": gates,
        "maximum_terminal_relative_error": max(
            float(row["relative_error"]) for row in electric_terminal
        ),
        "maximum_kcl_relative_to_terminal": max(
            float(row["relative_to_terminal_total"])
            for row in electric_kcl
        ),
        "production_formula_modified": False,
        "inputs": {
            name: {"path": str(path), "sha256": sha256(path)}
            for name, path in paths.items()
        },
        "outputs": output_hashes,
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    print(json.dumps(manifest, indent=2))
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--highfield", type=Path, required=True)
    parser.add_argument("--adjacent", type=Path, required=True)
    parser.add_argument("--stage-edges", type=Path, required=True)
    parser.add_argument("--terminal", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
