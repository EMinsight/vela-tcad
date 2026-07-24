#!/usr/bin/env python3
"""Attribute Minimal6 current error to low field, drive, and support."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from pathlib import Path

from pn2d_minimal6_diagnostics.factorial_attribution import (
    FACTORS,
    interaction_remainder,
    replacement_order_contributions,
    shapley_contributions,
)
from pn2d_minimal6_diagnostics.highfield_box_replay import (
    coefficient_weighted_mobility,
    field_limited_mobility,
)


LOW_SOURCES = ("vela", "sentaurus")
DRIVES = (
    "global_edge_qfp",
    "triangle_qfp",
    "native_element_electric_field",
)
SUPPORTS = ("global_edge", "native_elements")
REFERENCE_BRANCH = "sentaurus_native_final"
TARGET_BRANCH = "sentaurus_lowfield_element_electric_field"


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


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


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


def invert_low_field(
    carrier: str, final_mobility: float, field_V_per_m: float
) -> float:
    from pn2d_minimal6_diagnostics.highfield_box_replay import FIELD

    beta = FIELD[carrier]["beta"]
    velocity = FIELD[carrier]["saturation_velocity"]
    inverse_power = final_mobility ** (-beta) - (
        abs(field_V_per_m) / velocity
    ) ** beta
    if inverse_power <= 0.0:
        raise ValueError("global low-field mobility inversion is nonphysical")
    return inverse_power ** (-1.0 / beta)


def weighted(values: list[tuple[float, float]]) -> float:
    result = coefficient_weighted_mobility(values)
    if result["status"] != "valid":
        raise ValueError("factorial sample has zero geometric support")
    return float(result["mobility_m2_per_Vs"])


def summarize(values: list[float], weights: list[float]) -> dict[str, float]:
    absolute = [abs(value) for value in values]
    weight_sum = sum(weights)
    return {
        "median_abs_dex": statistics.median(absolute),
        "p95_abs_dex": quantile(absolute, 0.95),
        "maximum_abs_dex": max(absolute),
        "mean_signed_dex": statistics.mean(values),
        "current_weighted_mean_abs_dex": sum(
            weight * value for weight, value in zip(weights, absolute)
        )
        / weight_sum,
    }


def run(args: argparse.Namespace) -> dict[str, object]:
    paths = {
        "highfield": args.highfield.resolve(),
        "adjacent": args.adjacent.resolve(),
        "box_current": args.box_current.resolve(),
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
    adjacent: dict[
        tuple[str, float, str, int], list[dict[str, str]]
    ] = {}
    for row in rows(paths["adjacent"]):
        key = (
            row["topology"],
            float(row["bias_V"]),
            row["carrier"],
            int(row["edge_id"]),
        )
        adjacent.setdefault(key, []).append(row)
    box_rows = rows(paths["box_current"])
    reference = {
        (
            row["topology"],
            float(row["bias_V"]),
            row["carrier"],
            int(row["edge_id"]),
        ): row
        for row in box_rows
        if row["branch"] == REFERENCE_BRANCH
    }
    central = [
        row
        for row in box_rows
        if {int(row["node0"]), int(row["node1"])} == {1, 5}
    ]

    factorial_rows: list[dict[str, object]] = []
    order_rows: list[dict[str, object]] = []
    shapley_rows: list[dict[str, object]] = []
    interaction_rows: list[dict[str, object]] = []
    maximum_closure = 0.0

    for key, local_rows in sorted(adjacent.items()):
        ref = reference[key]
        ref_current = float(ref["reference_A_per_um"])
        ref_mobility_raw = ref["reference_mobility_m2_per_Vs"]
        edge_class = (
            "central_1_5"
            if {int(ref["node0"]), int(ref["node1"])} == {1, 5}
            else "ordinary_active"
        )
        kappa_sum = sum(float(row["kappa"]) for row in local_rows)
        if (
            kappa_sum == 0.0
            or ref_current == 0.0
            or ref_mobility_raw == ""
        ):
            continue
        ref_mobility = float(ref_mobility_raw)
        global_fields = {
            float(row["vela_triangle_local_edge_qf_field_V_per_m"])
            for row in local_rows
        }
        if len(global_fields) != 1:
            raise ValueError(f"global edge field is inconsistent for {key}")
        global_qfp_field = global_fields.pop()
        production_values = {
            row["vela_production_global_edge_mobility_m2_per_Vs"]
            for row in local_rows
            if row["vela_production_global_edge_mobility_m2_per_Vs"]
        }
        if len(production_values) != 1:
            raise ValueError(f"global production mobility missing for {key}")
        production_final = float(production_values.pop())
        vela_global_low = invert_low_field(
            key[2], production_final, global_qfp_field
        )
        element: list[dict[str, float]] = []
        for row in local_rows:
            cell = int(row["cell_id"])
            source = highfield[(key[0], key[1], cell, key[2])]
            element.append(
                {
                    "kappa": float(row["kappa"]),
                    "vela_low": float(
                        source["vela_cell_average_low_field_m2_per_Vs"]
                    ),
                    "sentaurus_low": float(
                        source["sentaurus_low_field_m2_per_Vs"]
                    ),
                    "global_edge_qfp": global_qfp_field,
                    "triangle_qfp": float(
                        source["triangle_qf_field_V_per_m"]
                    ),
                    "native_element_electric_field": float(
                        source["electric_field_V_per_m"]
                    ),
                }
            )
        global_low = {
            "vela": vela_global_low,
            "sentaurus": weighted(
                [
                    (item["kappa"], item["sentaurus_low"])
                    for item in element
                ]
            ),
        }
        global_drive = {
            drive: weighted(
                [(item["kappa"], item[drive]) for item in element]
            )
            for drive in DRIVES
        }
        branch_value: dict[tuple[str, str, str], float] = {}
        for low_source in LOW_SOURCES:
            for drive in DRIVES:
                for support in SUPPORTS:
                    if support == "global_edge":
                        mobility = field_limited_mobility(
                            key[2],
                            global_low[low_source],
                            global_drive[drive],
                        )
                    else:
                        mobility = weighted(
                            [
                                (
                                    item["kappa"],
                                    field_limited_mobility(
                                        key[2],
                                        item[f"{low_source}_low"],
                                        item[drive],
                                    ),
                                )
                                for item in element
                            ]
                        )
                    signed = math.log10(mobility / ref_mobility)
                    candidate_current = ref_current * mobility / ref_mobility
                    branch_value[(low_source, drive, support)] = signed
                    factorial_rows.append(
                        {
                            "topology": key[0],
                            "bias_V": key[1],
                            "carrier": key[2],
                            "edge_id": key[3],
                            "node0": int(ref["node0"]),
                            "node1": int(ref["node1"]),
                            "edge_class": edge_class,
                            "low_field_source": low_source,
                            "drive": drive,
                            "support": support,
                            "candidate_mobility_m2_per_Vs": mobility,
                            "reference_mobility_m2_per_Vs": ref_mobility,
                            "candidate_A_per_um": candidate_current,
                            "reference_A_per_um": ref_current,
                            "signed_log10_ratio_dex": signed,
                            "absolute_log10_error_dex": abs(signed),
                            "same_edge_paired_baseline": True,
                            "status": "valid",
                        }
                    )
        values: dict[frozenset[str], float] = {}
        for mask in range(8):
            active = frozenset(
                factor
                for bit, factor in enumerate(FACTORS)
                if mask & (1 << bit)
            )
            values[active] = branch_value[
                (
                    "sentaurus" if "low_field" in active else "vela",
                    (
                        "native_element_electric_field"
                        if "drive" in active
                        else "global_edge_qfp"
                    ),
                    (
                        "native_elements"
                        if "support" in active
                        else "global_edge"
                    ),
                )
            ]
        increments = replacement_order_contributions(values)
        for row in increments:
            order_rows.append(
                {
                    "topology": key[0],
                    "bias_V": key[1],
                    "carrier": key[2],
                    "edge_id": key[3],
                    "edge_class": edge_class,
                    "order": row["order"],
                    "step": row["step"],
                    "factor": row["factor"],
                    "before_signed_dex": row["before"],
                    "after_signed_dex": row["after"],
                    "increment_dex": row["increment"],
                    "reference_abs_current_A_per_um": abs(ref_current),
                    "same_edge_paired_baseline": True,
                }
            )
        shapley = shapley_contributions(values)
        closure = (
            values[frozenset(FACTORS)]
            - values[frozenset()]
            - sum(shapley.values())
        )
        maximum_closure = max(maximum_closure, abs(closure))
        interaction = interaction_remainder(values)
        shapley_rows.append(
            {
                "topology": key[0],
                "bias_V": key[1],
                "carrier": key[2],
                "edge_id": key[3],
                "edge_class": edge_class,
                "baseline_signed_dex": values[frozenset()],
                "target_signed_dex": values[frozenset(FACTORS)],
                "low_field_shapley_dex": shapley["low_field"],
                "drive_shapley_dex": shapley["drive"],
                "support_shapley_dex": shapley["support"],
                "interaction_remainder_dex": interaction,
                "shapley_closure_dex": closure,
                "reference_abs_current_A_per_um": abs(ref_current),
                "same_edge_paired_baseline": True,
            }
        )
        interaction_rows.append(
            {
                "topology": key[0],
                "bias_V": key[1],
                "carrier": key[2],
                "edge_id": key[3],
                "edge_class": edge_class,
                "interaction_remainder_dex": interaction,
                "reference_abs_current_A_per_um": abs(ref_current),
            }
        )

    factor_summary: list[dict[str, object]] = []
    for edge_class_name in ("ordinary_active", "central_1_5"):
        for carrier in ("electron", "hole"):
            selected = [
                row
                for row in shapley_rows
                if row["carrier"] == carrier
                and row["edge_class"] == edge_class_name
            ]
            for factor in FACTORS:
                values = [
                    float(row[f"{factor}_shapley_dex"])
                    for row in selected
                ]
                weights = [
                    float(row["reference_abs_current_A_per_um"])
                    for row in selected
                ]
                factor_summary.append(
                    {
                        "edge_class": edge_class_name,
                        "carrier": carrier,
                        "factor": factor,
                        "sample_count": len(values),
                        **summarize(values, weights),
                    }
                )
    interaction_summary: list[dict[str, object]] = []
    for edge_class_name in ("ordinary_active", "central_1_5"):
        for carrier in ("electron", "hole"):
            selected = [
                row
                for row in interaction_rows
                if row["carrier"] == carrier
                and row["edge_class"] == edge_class_name
            ]
            interaction_summary.append(
                {
                    "edge_class": edge_class_name,
                    "carrier": carrier,
                    "sample_count": len(selected),
                    **summarize(
                        [
                            float(row["interaction_remainder_dex"])
                            for row in selected
                        ],
                        [
                            float(row["reference_abs_current_A_per_um"])
                            for row in selected
                        ],
                    ),
                }
            )

    aggregate: dict[str, dict[str, float]] = {}
    for factor in FACTORS:
        factor_rows = [
            row
            for row in factor_summary
            if row["factor"] == factor
            and row["edge_class"] == "ordinary_active"
        ]
        aggregate[factor] = {
            "unweighted": statistics.mean(
                float(row["median_abs_dex"]) for row in factor_rows
            ),
            "weighted": statistics.mean(
                float(row["current_weighted_mean_abs_dex"])
                for row in factor_rows
            ),
        }
    unweighted_rank = sorted(
        FACTORS, key=lambda factor: aggregate[factor]["unweighted"], reverse=True
    )
    weighted_rank = sorted(
        FACTORS, key=lambda factor: aggregate[factor]["weighted"], reverse=True
    )
    ranking_stable = unweighted_rank == weighted_rank
    interaction_score = statistics.mean(
        float(row["current_weighted_mean_abs_dex"])
        for row in interaction_summary
        if row["edge_class"] == "ordinary_active"
    )
    top_factor = weighted_rank[0]
    if interaction_score > aggregate[top_factor]["weighted"]:
        outcome = "interaction_dominant"
    elif ranking_stable:
        outcome = {
            "low_field": "low_field_coefficient_dominant",
            "drive": "driving_force_dominant",
            "support": "support_dominant",
        }[top_factor]
    else:
        outcome = "mixed_bounded_contributions"

    central_rows = [
        {
            **row,
            "reference_abs_A_per_um": abs(float(row["reference_A_per_um"])),
            "candidate_abs_A_per_um": abs(float(row["candidate_A_per_um"])),
        }
        for row in central
    ]
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    outputs = {
        "factorial_samples.csv": factorial_rows,
        "replacement_orders.csv": order_rows,
        "shapley_samples.csv": shapley_rows,
        "factor_summary.csv": factor_summary,
        "interaction_summary.csv": interaction_summary,
        "central_tail.csv": central_rows,
    }
    for name, values in outputs.items():
        write_csv(output / name, values)
    report = [
        "# PN2D Minimal6 current-factor attribution",
        "",
        f"Typed outcome: `{outcome}`",
        "",
        f"Maximum Shapley closure: `{maximum_closure:.6e} dex`.",
        "",
        f"Unweighted rank: `{' > '.join(unweighted_rank)}`.",
        "",
        f"Current-weighted rank: `{' > '.join(weighted_rank)}`.",
        "",
        f"Ranking stable: `{str(ranking_stable).lower()}`.",
        "",
        "| Edge class | Carrier | Factor | Median abs dex | P95 abs dex | "
        "Max abs dex | Current-weighted mean abs dex |",
        "|---|---|---|---:|---:|---:|---:|",
    ]
    for row in factor_summary:
        report.append(
            f"| {row['edge_class']} | {row['carrier']} | "
            f"{row['factor']} | "
            f"{float(row['median_abs_dex']):.6g} | "
            f"{float(row['p95_abs_dex']):.6g} | "
            f"{float(row['maximum_abs_dex']):.6g} | "
            f"{float(row['current_weighted_mean_abs_dex']):.6g} |"
        )
    report.extend(
        [
            "",
            "The central 1-5 edge is reported separately in `central_tail.csv` "
            "and in the `central_1_5` summary rows. Rankings use only "
            "`ordinary_active` rows, so its near-zero current cannot dominate "
            "the unweighted or current-weighted conclusion.",
            "",
            "No fitted mobility, field scale, or edge coefficient was used.",
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
        "status": "valid" if maximum_closure <= 1.0e-12 else "failed",
        "experiment": "pn2d_minimal6_current_factorial",
        "typed_outcome": outcome,
        "active_sample_count": len(shapley_rows),
        "factorial_sample_count": len(factorial_rows),
        "replacement_increment_count": len(order_rows),
        "maximum_shapley_closure_dex": maximum_closure,
        "unweighted_rank": unweighted_rank,
        "current_weighted_rank": weighted_rank,
        "ranking_stable": ranking_stable,
        "same_edge_paired_baseline": True,
        "fitted_parameter_count": 0,
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
    parser.add_argument("--box-current", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
