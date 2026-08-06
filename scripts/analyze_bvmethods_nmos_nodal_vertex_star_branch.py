#!/usr/bin/env python3
"""Summarize the BVmethods NMOS nodal-vertex-star high-voltage branch."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any

import audit_bvmethods_nmos_hotspot_slope_and_criteria as base


Q_C = 1.602176634e-19
BRANCH = (
    base.RUN
    / "vela_validation/btbt_e2_iic_qf_vector_nodal_vertex_star_branch_6p4_7p1_20260805"
)
DEFAULT_EDGE_PATHS = [
    BRANCH / "sg_avalanche_edges.csv",
    BRANCH / "segment_6p8_6p9/sg_avalanche_edges.csv",
    BRANCH / "segment_7p0/sg_avalanche_edges.csv",
    BRANCH / "segment_7p1/sg_avalanche_edges.csv",
]
DEFAULT_TERMINAL_PATHS = [
    base.RUN
    / "vela_validation/btbt_e2_iic_qf_vector_nodal_vertex_star_fixed6p4_20260805/terminal_current_method_compare.csv",
    BRANCH / "fixed_replay_6p5/terminal_current_method_compare.csv",
    BRANCH / "fixed_replay_6p6/terminal_current_method_compare.csv",
    BRANCH / "fixed_replay_6p7/terminal_current_method_compare.csv",
    BRANCH / "segment_6p8_6p9/terminal_current_method_compare.csv",
    BRANCH / "segment_7p0/terminal_current_method_compare.csv",
    BRANCH / "segment_7p1/terminal_current_method_compare.csv",
]
DEFAULT_OUT = BRANCH / "analysis/branch_closure"


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def value(row: dict[str, str], key: str) -> float:
    return float(row.get(key, "0") or 0.0)


def write_csv(path: Path, records: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def load_edges(paths: list[Path]) -> dict[float, list[dict[str, str]]]:
    grouped: dict[float, list[dict[str, str]]] = {}
    for path in paths:
        local: dict[float, list[dict[str, str]]] = {}
        for row in rows(path):
            local.setdefault(round(value(row, "bias_V"), 9), []).append(row)
        grouped.update(local)
    return grouped


def load_drain_currents(paths: list[Path]) -> dict[float, float]:
    output: dict[float, float] = {}
    for path in paths:
        for row in rows(path):
            if row["contact"] != "drain":
                continue
            output[round(value(row, "bias_V"), 9)] = abs(
                value(row, "I_sgflux_A_per_um")
            )
    return output


def source_current(edge_rows: list[dict[str, str]]) -> float:
    native_source = sum(value(row, "edge_source_integral") for row in edge_rows)
    return Q_C * native_source * 1.0e-12


def linear_crossing(records: list[dict[str, Any]], key: str) -> float | None:
    ordered = sorted(records, key=lambda row: row["bias_V"])
    for left, right in zip(ordered, ordered[1:]):
        y0, y1 = left[key], right[key]
        if y0 <= 0.0 < y1 and y1 != y0:
            return left["bias_V"] + (right["bias_V"] - left["bias_V"]) * (
                -y0 / (y1 - y0)
            )
    return None


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--edges", type=Path, action="append")
    parser.add_argument("--terminal", type=Path, action="append")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()

    edge_paths = args.edges or DEFAULT_EDGE_PATHS
    terminal_paths = args.terminal or DEFAULT_TERMINAL_PATHS
    new_edges = load_edges(edge_paths)
    currents = load_drain_currents(terminal_paths)
    missing = sorted(new_edges.keys() - currents.keys())
    if missing:
        raise ValueError(f"missing drain currents at biases {missing}")

    old_edges = base.load_vela_edges(base.DEFAULT_VELA_EDGES)
    sent_plot = base.plot_index(base.DEFAULT_SENT_PLOT)
    records: list[dict[str, Any]] = []
    for bias in sorted(new_edges):
        drain = currents[bias]
        new_iava = source_current(new_edges[bias])
        old_iava = source_current(old_edges[bias]) if bias in old_edges else math.nan
        sent = sent_plot.get(bias, {})
        sent_id = abs(sent.get("drain TotalCurrent", math.nan))
        sent_native_source = abs(
            sent.get("IntegrSemiconductor AvalancheGeneration", math.nan)
        )
        sent_iava = Q_C * sent_native_source * 1.0e-12
        records.append({
            "bias_V": bias,
            "drain_current_A_per_um": drain,
            "nodal_vertex_star_iava_A_per_um": new_iava,
            "nodal_vertex_star_iava_over_id": new_iava / drain,
            "nodal_vertex_star_iava_minus_id_A_per_um": new_iava - drain,
            "legacy_edge_adjacent_iava_A_per_um": old_iava,
            "legacy_edge_adjacent_iava_over_id": old_iava / drain,
            "sentaurus_drain_current_A_per_um": sent_id,
            "sentaurus_iava_A_per_um": sent_iava,
            "sentaurus_iava_over_id": sent_iava / sent_id,
            "nodal_vertex_star_over_sentaurus_iava": new_iava / sent_iava,
        })

    new_crossing = linear_crossing(
        records, "nodal_vertex_star_iava_minus_id_A_per_um"
    )
    legacy_records = [
        {
            "bias_V": row["bias_V"],
            "difference": row["legacy_edge_adjacent_iava_A_per_um"]
            - row["drain_current_A_per_um"],
        }
        for row in records
        if math.isfinite(row["legacy_edge_adjacent_iava_A_per_um"])
    ]
    legacy_crossing = linear_crossing(legacy_records, "difference")
    sentaurus_dense_crossing = 6.734425890478791
    sentaurus_sparse_crossing = 6.377494277837012
    summary = {
        "biases_V": [row["bias_V"] for row in records],
        "nodal_vertex_star_current_source_crossing_V": new_crossing,
        "legacy_edge_adjacent_current_source_crossing_V": legacy_crossing,
        "sentaurus_dense_current_source_crossing_V": sentaurus_dense_crossing,
        "sentaurus_sparse_official_linear_crossing_V": sentaurus_sparse_crossing,
        "nodal_vertex_star_minus_sentaurus_dense_crossing_V": (
            new_crossing - sentaurus_dense_crossing if new_crossing else math.nan
        ),
        "nodal_vertex_star_minus_legacy_crossing_V": (
            new_crossing - legacy_crossing
            if new_crossing and legacy_crossing else math.nan
        ),
        "iava_ratio_to_sentaurus_at_6p4": next(
            row["nodal_vertex_star_over_sentaurus_iava"]
            for row in records if row["bias_V"] == 6.4
        ),
        "iava_ratio_to_sentaurus_at_7p0": next(
            row["nodal_vertex_star_over_sentaurus_iava"]
            for row in records if row["bias_V"] == 7.0
        ),
    }
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "branch_current_source_compare.csv", records)
    (args.out_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, allow_nan=True) + "\n", encoding="utf-8"
    )
    crossing_text = (
        f"`{new_crossing:.9f} V`" if new_crossing is not None
        else "not bracketed by the supplied bias range"
    )
    report = [
        "# BVmethods NMOS nodal-vertex-star high-voltage branch",
        "",
        f"- Bias coverage: `{records[0]['bias_V']:.1f}--{records[-1]['bias_V']:.1f} V`.",
        f"- New current-source crossing: {crossing_text}.",
        f"- Legacy edge-adjacent crossing on the same states: `{legacy_crossing:.9f} V`.",
        f"- Sentaurus dense crossing: `{sentaurus_dense_crossing:.9f} V`.",
        f"- New minus Sentaurus dense: `{new_crossing - sentaurus_dense_crossing:+.9f} V`.",
        f"- New minus legacy: `{new_crossing - legacy_crossing:+.9f} V`.",
        "",
        "The official sparse ABA-coupled linear interpolation is "
        f"`{sentaurus_sparse_crossing:.9f} V`; it is retained as a separate oracle "
        "and is not mixed with the dense same-bias crossing.",
    ]
    (args.out_dir / "summary.md").write_text(
        "\n".join(report) + "\n", encoding="utf-8"
    )
    print(args.out_dir / "summary.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
