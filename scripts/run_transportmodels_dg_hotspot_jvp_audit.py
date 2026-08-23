#!/usr/bin/env python3
"""Audit the TransportModels DG hotspot Jacobian with centered Jv probes."""

from __future__ import annotations

import argparse
import copy
import csv
import json
import subprocess
from pathlib import Path
from typing import Any


DEFAULT_HOTSPOT_NODES = (706, 705, 794, 342, 743, 712, 714, 711, 793, 716)
DEFAULT_AMPLITUDES_V = (1.0e-4, 1.0e-6, 1.0e-8)


def build_document(
    base: dict[str, Any],
    *,
    state_file: Path,
    output_csv: Path,
    gate_bias_V: float,
    hotspot_nodes: tuple[int, ...] = DEFAULT_HOTSPOT_NODES,
    amplitudes_V: tuple[float, ...] = DEFAULT_AMPLITUDES_V,
) -> dict[str, Any]:
    """Build a global Jv audit with perturbations localized to DG hotspots."""
    document = copy.deepcopy(base)
    document["simulation_type"] = "newton_jvp_probe"
    document["state_file"] = str(state_file.resolve())
    document["output_csv"] = str(output_csv.resolve())
    document.pop("sweep", None)
    document.pop("log_file", None)

    for contact in document["contacts"]:
        if contact["name"] == "gate":
            contact["bias"] = gate_bias_V

    solver = document.setdefault("solver", {})
    solver["quasi_fermi_reference"] = "contact_majority"

    directions: list[dict[str, Any]] = []
    for amplitude in amplitudes_V:
        for mode in ("psi", "phin", "psi_minus_phin"):
            directions.append(
                {
                    "name": f"hotspot_cluster_{mode}_h{amplitude:.0e}",
                    "mode": mode,
                    "amplitude_V": amplitude,
                    "exclude_contacts": False,
                    "node_ids": list(hotspot_nodes),
                    "node_index_base": 0,
                    "adjacent_cell_rings": 0,
                }
            )

    # Resolve whether a single hotspot dominates a cluster cancellation.
    for node in hotspot_nodes[:4]:
        for amplitude in (1.0e-6, 1.0e-8):
            directions.append(
                {
                    "name": f"node{node}_phin_h{amplitude:.0e}",
                    "mode": "phin",
                    "amplitude_V": amplitude,
                    "exclude_contacts": False,
                    "node_ids": [node],
                    "node_index_base": 0,
                    "adjacent_cell_rings": 0,
                }
            )

    document["directions"] = directions
    document["_dg_hotspot_jvp_audit"] = {
        "gate_bias_V": gate_bias_V,
        "hotspot_nodes_zero_based": list(hotspot_nodes),
        "amplitudes_V": list(amplitudes_V),
        "reference_coordinates": "contact_majority",
        "finite_difference": "centered residual difference",
        "response_scope": "global Jv with hotspot-localized perturbations",
    }
    return document


def scaled_error(error: float, analytic: float, finite_difference: float) -> float:
    return error / max(abs(analytic), abs(finite_difference), 1.0e-300)


def summarize(csv_path: Path) -> dict[str, Any]:
    rows: list[dict[str, Any]] = []
    with csv_path.open(newline="", encoding="utf-8") as handle:
        for raw in csv.DictReader(handle):
            row: dict[str, Any] = dict(raw)
            row["scaled_relative_error"] = scaled_error(
                float(raw["absolute_error"]),
                float(raw["analytic_norm"]),
                float(raw["finite_difference_norm"]),
            )
            for block in ("psi", "phin", "phip"):
                analytic = float(raw[f"analytic_{block}_norm"])
                finite_difference = float(raw[f"finite_difference_{block}_norm"])
                # Recover the absolute block error from the runner's legacy
                # max(1, fd_norm) normalization.
                legacy_relative = float(raw[f"{block}_relative_error"])
                absolute_error = legacy_relative * max(1.0, finite_difference)
                row[f"{block}_scaled_relative_error"] = scaled_error(
                    absolute_error, analytic, finite_difference
                )
            rows.append(row)

    worst = sorted(rows, key=lambda row: row["scaled_relative_error"], reverse=True)
    return {
        "row_count": len(rows),
        "max_scaled_relative_error": worst[0]["scaled_relative_error"] if worst else 0.0,
        "worst_directions": [
            {
                "direction": row["direction"],
                "mode": row["mode"],
                "amplitude_V": float(row["amplitude_V"]),
                "scaled_relative_error": row["scaled_relative_error"],
                "psi_scaled_relative_error": row["psi_scaled_relative_error"],
                "phin_scaled_relative_error": row["phin_scaled_relative_error"],
                "phip_scaled_relative_error": row["phip_scaled_relative_error"],
            }
            for row in worst[:8]
        ],
        "rows": rows,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-config", type=Path, required=True)
    parser.add_argument("--state-file", type=Path, required=True)
    parser.add_argument("--runner", type=Path, default=Path("build-release/vela_example_runner.exe"))
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--gate-bias-V", type=float, default=-0.4)
    args = parser.parse_args()

    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    output_csv = output_dir / "dg_hotspot_jvp.csv"
    output_config = output_dir / "config.json"
    with args.base_config.resolve().open(encoding="utf-8") as handle:
        base = json.load(handle)
    document = build_document(
        base,
        state_file=args.state_file,
        output_csv=output_csv,
        gate_bias_V=args.gate_bias_V,
    )
    output_config.write_text(json.dumps(document, indent=2) + "\n", encoding="utf-8")

    completed = subprocess.run(
        [str(args.runner.resolve()), "--config", str(output_config)],
        check=False,
        capture_output=True,
        text=True,
    )
    (output_dir / "stdout.txt").write_text(completed.stdout, encoding="utf-8")
    (output_dir / "stderr.txt").write_text(completed.stderr, encoding="utf-8")
    execution = {"returncode": completed.returncode, "config": str(output_config)}
    (output_dir / "execution.json").write_text(
        json.dumps(execution, indent=2) + "\n", encoding="utf-8"
    )
    if completed.returncode != 0:
        print(json.dumps(execution, indent=2))
        return completed.returncode

    audit_summary = summarize(output_csv)
    (output_dir / "summary.json").write_text(
        json.dumps(audit_summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({key: value for key, value in audit_summary.items() if key != "rows"}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
