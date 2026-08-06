#!/usr/bin/env python3
"""Compare Vela peak/saddle diagnostics with Sentaurus WriteAll paths.

The Sentaurus log prints one numbered path inventory followed by the terminal
contact table.  This script associates each complete inventory with that drain
voltage, collapses exact local-peak aliases, and ranks the remaining paths by
the validated arithmetic carrier mean.  Vela rows retain every numbered peak
and its merge-tree saddle diagnostics.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
DEFAULT_SENT_LOG = (
    REPO
    / "build-release/reference_tcad/bvmethods_sentaurus2018/run01"
    / "sentaurus_eparallel_vector_exact_20260804/extracted/iic_multibias_des.log"
)
DEFAULT_VELA_BRANCH = (
    REPO
    / "build-release/reference_tcad/bvmethods_sentaurus2018/run01"
    / "vela_validation/adaptive_minority_qf_branch_20260806"
)
DEFAULT_OUT = (
    REPO
    / "build-release/reference_tcad/bvmethods_sentaurus2018/run01"
    / "vela_validation/peak_saddle_retention_audit_20260806"
)

PATH_RE = re.compile(r"^Path number\s+(\d+)\s*$")
FIELD_RE = re.compile(r"^(Maximum Field|Electron|Hole):\s+([+\-0-9.eE]+)\s*$")
DRAIN_RE = re.compile(r"^drain\s+([+\-0-9.eE]+)\s+")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def sentaurus_inventories(path: Path) -> list[dict[str, object]]:
    completed: list[list[dict[str, float | int]]] = []
    pending: list[dict[str, float | int]] | None = None
    current: list[dict[str, float | int]] = []
    active: dict[str, float | int] | None = None
    output: list[dict[str, object]] = []
    sequence = 0
    for raw in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        match = PATH_RE.match(line)
        if match:
            active = {"path_number": int(match.group(1))}
            current.append(active)
            continue
        if line == "Best Path":
            complete = [
                row
                for row in current
                if "electron" in row and "hole" in row and "maximum_field_V_per_cm" in row
            ]
            pending = complete or None
            current = []
            active = None
            continue
        field = FIELD_RE.match(line)
        if field is not None and active is not None:
            label, value = field.groups()
            active[
                {
                    "Maximum Field": "maximum_field_V_per_cm",
                    "Electron": "electron",
                    "Hole": "hole",
                }[label]
            ] = float(value)
            continue
        drain = DRAIN_RE.match(line)
        if drain is None or pending is None:
            continue
        bias = float(drain.group(1))
        completed.append(pending)
        for row in pending:
            electron = float(row["electron"])
            hole = float(row["hole"])
            output.append(
                {
                    "inventory_sequence": sequence,
                    "bias_V": bias,
                    "path_number": int(row["path_number"]),
                    "maximum_field_V_per_m": 100.0
                    * float(row["maximum_field_V_per_cm"]),
                    "electron_ionization_integral": electron,
                    "hole_ionization_integral": hole,
                    "mean_ionization_integral": 0.5 * (electron + hole),
                }
            )
        sequence += 1
        pending = None
    return output


def final_inventory_per_bias(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    latest: dict[float, int] = {}
    for row in rows:
        latest[float(row["bias_V"])] = int(row["inventory_sequence"])
    return [
        row
        for row in rows
        if int(row["inventory_sequence"]) == latest[float(row["bias_V"])]
    ]


def distinct_ranked(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    by_bias: dict[float, list[dict[str, object]]] = {}
    for row in rows:
        by_bias.setdefault(float(row["bias_V"]), []).append(row)
    result: list[dict[str, object]] = []
    for bias, inventory in sorted(by_bias.items()):
        distinct: dict[tuple[float, float], dict[str, object]] = {}
        for row in inventory:
            key = (
                float(row["electron_ionization_integral"]),
                float(row["hole_ionization_integral"]),
            )
            if key not in distinct:
                distinct[key] = dict(row, multiplicity=1)
            else:
                distinct[key]["multiplicity"] = int(distinct[key]["multiplicity"]) + 1
        ranked = sorted(
            distinct.values(),
            key=lambda row: (
                -float(row["mean_ionization_integral"]),
                int(row["path_number"]),
            ),
        )
        for rank, row in enumerate(ranked, 1):
            result.append(
                {
                    "bias_V": bias,
                    "path_rank": rank,
                    "path_number": row["path_number"],
                    "multiplicity": row["multiplicity"],
                    "maximum_field_V_per_m": row["maximum_field_V_per_m"],
                    "electron_ionization_integral": row["electron_ionization_integral"],
                    "hole_ionization_integral": row["hole_ionization_integral"],
                    "mean_ionization_integral": row["mean_ionization_integral"],
                }
            )
    return result


def vela_rows(branch: Path) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    for path in sorted(branch.glob("bias_*/postprocess_only/path_ionization_integrals.csv")):
        for row in read_csv(path):
            output.append(
                {
                    "bias_V": float(row["bias_V"]),
                    "path_rank": int(row["path_rank"]),
                    "physical_path_rank": int(row.get("physical_path_rank", 0)),
                    "seed_node_id": int(row["seed_node_id"]),
                    "parent_peak_node_id": row["parent_peak_node_id"],
                    "physical_path_group_id": row.get("physical_path_group_id", ""),
                    "seed_field_V_per_m": float(row["seed_field_V_per_m"]),
                    "saddle_field_V_per_m": float(row["saddle_field_V_per_m"]),
                    "peak_prominence_V_per_m": float(row["peak_prominence_V_per_m"]),
                    "peak_prominence_ratio": float(row["peak_prominence_ratio"]),
                    "seed_electron_qf_relative_magnitude": float(
                        row.get("seed_electron_qf_relative_magnitude", 0.0)
                    ),
                    "seed_hole_qf_relative_magnitude": float(
                        row.get("seed_hole_qf_relative_magnitude", 0.0)
                    ),
                    "electron_ionization_integral": float(
                        row["electron_ionization_integral"]
                    ),
                    "hole_ionization_integral": float(row["hole_ionization_integral"]),
                    "mean_ionization_integral": float(row["mean_ionization_integral"]),
                }
            )
    return sorted(output, key=lambda row: (row["bias_V"], row["path_rank"]))


def vela_distinct_ranked(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    """Collapse numbered aliases while preserving the strongest row per corridor."""
    by_bias: dict[float, list[dict[str, object]]] = {}
    for row in rows:
        by_bias.setdefault(float(row["bias_V"]), []).append(row)
    output: list[dict[str, object]] = []
    for bias, inventory in sorted(by_bias.items()):
        groups: dict[str, dict[str, object]] = {}
        for row in inventory:
            raw_group = str(row.get("physical_path_group_id", ""))
            group = raw_group or f"raw:{int(row['path_rank'])}"
            existing = groups.get(group)
            if existing is None or float(row["mean_ionization_integral"]) > float(
                existing["mean_ionization_integral"]
            ):
                groups[group] = row
        ranked = sorted(
            groups.values(),
            key=lambda row: (
                -float(row["mean_ionization_integral"]),
                int(row["path_rank"]),
            ),
        )
        for rank, row in enumerate(ranked, 1):
            output.append(
                {
                    "bias_V": bias,
                    "physical_path_rank": rank,
                    "raw_path_rank": row["path_rank"],
                    "physical_path_group_id": row["physical_path_group_id"],
                    "seed_node_id": row["seed_node_id"],
                    "seed_field_V_per_m": row["seed_field_V_per_m"],
                    "saddle_field_V_per_m": row["saddle_field_V_per_m"],
                    "peak_prominence_ratio": row["peak_prominence_ratio"],
                    "seed_electron_qf_relative_magnitude": row[
                        "seed_electron_qf_relative_magnitude"
                    ],
                    "seed_hole_qf_relative_magnitude": row[
                        "seed_hole_qf_relative_magnitude"
                    ],
                    "electron_ionization_integral": row[
                        "electron_ionization_integral"
                    ],
                    "hole_ionization_integral": row["hole_ionization_integral"],
                    "mean_ionization_integral": row["mean_ionization_integral"],
                }
            )
    return output


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sentaurus-log", type=Path, default=DEFAULT_SENT_LOG)
    parser.add_argument("--vela-branch", type=Path, default=DEFAULT_VELA_BRANCH)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()

    all_sent = sentaurus_inventories(args.sentaurus_log)
    final_sent = final_inventory_per_bias(all_sent)
    ranked_sent = distinct_ranked(final_sent)
    vela = vela_rows(args.vela_branch)
    ranked_vela = vela_distinct_ranked(vela)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.out_dir / "sentaurus_writeall_all_inventories.csv", all_sent)
    write_csv(args.out_dir / "sentaurus_writeall_final_by_bias.csv", final_sent)
    write_csv(args.out_dir / "sentaurus_distinct_ranked_by_bias.csv", ranked_sent)
    write_csv(args.out_dir / "vela_peak_saddle_by_bias.csv", vela)
    write_csv(args.out_dir / "vela_distinct_physical_ranked_by_bias.csv", ranked_vela)

    selected_sent = [
        row
        for row in ranked_sent
        if int(row["path_rank"]) <= 4 and float(row["bias_V"]) >= 6.0
    ]
    selected_vela = [
        row
        for row in vela
        if int(row["path_rank"]) <= 5
        and int(row["seed_node_id"]) in {327, 990, 1381, 1460, 1461, 1462}
    ]
    summary = {
        "sentaurus_source": str(args.sentaurus_log),
        "vela_source": str(args.vela_branch),
        "sentaurus_inventory_count": 1
        + max(int(row["inventory_sequence"]) for row in all_sent),
        "sentaurus_distinct_bias_count": len({row["bias_V"] for row in final_sent}),
        "vela_bias_count": len({row["bias_V"] for row in vela}),
        "selected_sentaurus_top_paths": selected_sent,
        "selected_vela_peaks": selected_vela,
        "vela_rank3_by_bias": [
            row for row in ranked_vela if int(row["physical_path_rank"]) == 3
        ],
    }
    (args.out_dir / "summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
