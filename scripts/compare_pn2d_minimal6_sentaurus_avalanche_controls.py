#!/usr/bin/env python3
"""Compare default, AvalDensGradQF, and ElementVolumeAvalanche controls."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Any

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.diagnose_pn2d_minimal6_element_avalanche_replay import (
    Q_LEGACY_C,
    TARGET_BIASES,
    TOPOLOGIES,
    currentplot_targets,
    find_integral_name,
    parse_log,
    parse_plt,
)


SOURCE_INTEGRAL_TO_A_UM = Q_LEGACY_C * 1.0e-12
VARIANTS = {
    "default": (
        "run_default.out",
        "runtime_element_avalanche_probe_default.plt",
    ),
    "aval_dens_grad_qf": (
        "run_gradqf.out",
        "runtime_element_avalanche_probe_gradqf.plt",
    ),
    "element_volume_avalanche": (
        "run_elementvolume.out",
        "runtime_element_avalanche_probe_elementvolume.plt",
    ),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="", encoding="ascii") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def dex_ratio(value: float, reference: float) -> float | None:
    if value <= 0.0 or reference <= 0.0:
        return None
    return math.log10(value / reference)


def format_value(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, float):
        return format(value, ".17g")
    return str(value)


def run(raw_root: Path, output: Path) -> dict[str, Any]:
    raw_root = raw_root.resolve()
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    state_rows: list[dict[str, Any]] = []
    node_rows: list[dict[str, Any]] = []
    input_hashes = {}

    for topology in TOPOLOGIES:
        for variant, (log_name, plt_name) in VARIANTS.items():
            log_path = raw_root / topology / variant / log_name
            plt_path = raw_root / topology / variant / plt_name
            groups = parse_log(log_path)
            names, _ = parse_plt(plt_path)
            targets = {
                row["bias_V"]: row for row in currentplot_targets(plt_path)
            }
            e_integral = find_integral_name(names, "eAvalancheIntegral")
            h_integral = find_integral_name(names, "hAvalancheIntegral")
            total_integral = find_integral_name(
                [
                    name
                    for name in names
                    if "eAvalancheIntegral" not in name
                    and "hAvalancheIntegral" not in name
                ],
                "AvalancheIntegral",
            )
            for path in (log_path, plt_path):
                input_hashes[path.relative_to(raw_root).as_posix()] = sha256(path)
            for bias in TARGET_BIASES:
                current = targets[bias]
                runtime = next(
                    row
                    for row in groups["integrals"]
                    if row["bias_V"] == int(bias)
                )
                anode_total_name = next(
                    name
                    for name in names
                    if name.endswith("Anode TotalCurrent")
                    or name == "Anode TotalCurrent"
                )
                cathode_total_name = next(
                    name
                    for name in names
                    if name.endswith("Cathode TotalCurrent")
                    or name == "Cathode TotalCurrent"
                )
                state_rows.append(
                    {
                        "topology": topology,
                        "bias_V": bias,
                        "variant": variant,
                        "electron_qg_A_um": runtime["qg_n_A_um"],
                        "hole_qg_A_um": runtime["qg_p_A_um"],
                        "total_qg_A_um": runtime["qg_total_A_um"],
                        "currentplot_electron_qg_A_um": (
                            current[e_integral] * SOURCE_INTEGRAL_TO_A_UM
                        ),
                        "currentplot_hole_qg_A_um": (
                            current[h_integral] * SOURCE_INTEGRAL_TO_A_UM
                        ),
                        "currentplot_total_qg_A_um": (
                            current[total_integral] * SOURCE_INTEGRAL_TO_A_UM
                        ),
                        "anode_total_current_A_um": current[anode_total_name],
                        "cathode_total_current_A_um": current[
                            cathode_total_name
                        ],
                    }
                )
                for vertex in groups["vertices"]:
                    if (
                        vertex["bias_V"] == int(bias)
                        and int(vertex["vertex"]) < 6
                    ):
                        node_rows.append(
                            {
                                "topology": topology,
                                "bias_V": bias,
                                "variant": variant,
                                "vertex": int(vertex["vertex"]),
                                "electron_generation_cm3_s": vertex[
                                    "generation_n_cm3_s"
                                ],
                                "hole_generation_cm3_s": vertex[
                                    "generation_p_cm3_s"
                                ],
                                "total_generation_cm3_s": vertex[
                                    "generation_total_cm3_s"
                                ],
                            }
                        )

    by_state = {
        (row["topology"], row["bias_V"], row["variant"]): row
        for row in state_rows
    }
    comparison_rows: list[dict[str, Any]] = []
    for topology in TOPOLOGIES:
        for bias in TARGET_BIASES:
            baseline = by_state[(topology, bias, "default")]
            for variant in ("aval_dens_grad_qf", "element_volume_avalanche"):
                candidate = by_state[(topology, bias, variant)]
                comparison_rows.append(
                    {
                        "topology": topology,
                        "bias_V": bias,
                        "variant": variant,
                        "electron_qg_log10_ratio_to_default_dex": dex_ratio(
                            candidate["electron_qg_A_um"],
                            baseline["electron_qg_A_um"],
                        ),
                        "hole_qg_log10_ratio_to_default_dex": dex_ratio(
                            candidate["hole_qg_A_um"],
                            baseline["hole_qg_A_um"],
                        ),
                        "total_qg_log10_ratio_to_default_dex": dex_ratio(
                            candidate["total_qg_A_um"],
                            baseline["total_qg_A_um"],
                        ),
                        "anode_total_current_relative_change": (
                            candidate["anode_total_current_A_um"]
                            - baseline["anode_total_current_A_um"]
                        )
                        / max(
                            abs(baseline["anode_total_current_A_um"]),
                            1.0e-300,
                        ),
                        "cathode_total_current_relative_change": (
                            candidate["cathode_total_current_A_um"]
                            - baseline["cathode_total_current_A_um"]
                        )
                        / max(
                            abs(baseline["cathode_total_current_A_um"]),
                            1.0e-300,
                        ),
                    }
                )

    for collection in (state_rows, node_rows, comparison_rows):
        for row in collection:
            for key, value in list(row.items()):
                row[key] = format_value(value)
    outputs = {
        "control_states.csv": state_rows,
        "control_nodes.csv": node_rows,
        "control_vs_default.csv": comparison_rows,
    }
    output_hashes = {}
    for name, rows in outputs.items():
        path = output / name
        write_csv(path, rows)
        output_hashes[name] = sha256(path)

    numeric_comparison = [
        {
            key: (
                float(value)
                if key not in {"topology", "variant"} and value != ""
                else value
            )
            for key, value in row.items()
        }
        for row in comparison_rows
    ]
    manifest = {
        "schema_version": 1,
        "status": "valid_sentaurus_control_comparison",
        "experiment": "pn2d_minimal6_sentaurus_avalanche_controls",
        "state_count": len(state_rows),
        "node_count": len(node_rows),
        "input_sha256": input_hashes,
        "output_sha256": output_hashes,
        "max_absolute_total_qg_change_dex": {
            variant: max(
                abs(row["total_qg_log10_ratio_to_default_dex"])
                for row in numeric_comparison
                if row["variant"] == variant
            )
            for variant in ("aval_dens_grad_qf", "element_volume_avalanche")
        },
        "max_absolute_terminal_relative_change": {
            variant: max(
                abs(row[key])
                for row in numeric_comparison
                if row["variant"] == variant
                for key in (
                    "anode_total_current_relative_change",
                    "cathode_total_current_relative_change",
                )
            )
            for variant in ("aval_dens_grad_qf", "element_volume_avalanche")
        },
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="ascii",
    )
    return manifest


def main() -> int:
    args = parse_args()
    manifest = run(args.raw_root, args.output)
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
