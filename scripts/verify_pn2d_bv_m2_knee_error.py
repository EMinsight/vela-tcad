#!/usr/bin/env python3
"""Verify the M2 BV knee read-only diagnostic output."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    result = json.loads((root / "diagnostic.json").read_text(encoding="utf-8"))
    bias = rows(root / "bias_summary.csv")
    carrier = rows(root / "carrier_source_components.csv")
    edge = rows(root / "edge_current_summary.csv")
    mapping = rows(root / "source_mapping_summary.csv")

    require(result["status"] == "passed", "diagnostic did not pass")
    require(result["observation_only"], "diagnostic is not observation-only")
    require(
        not result["acceptance_thresholds_modified"],
        "acceptance thresholds were modified",
    )
    require(result["determinism"]["all_equal"], "Vela A/B outputs differ")
    require(len(bias) == 11 and len(edge) == 11 and len(mapping) == 11, "wrong knee lattice")
    require(len(carrier) == 22, "carrier source decomposition is incomplete")

    growth = result["knee_error_growth"]
    require(
        abs(growth["terminal_error_dex"] - growth["integrated_source_error_dex"])
        < 0.002,
        "source-error growth does not track terminal-error growth",
    )
    for key in (
        "electron_cell_current_deficit_dex",
        "hole_cell_current_deficit_dex",
        "electron_density_deficit_dex",
        "hole_density_deficit_dex",
    ):
        require(
            abs(growth[key] - growth["terminal_error_dex"]) < 0.01,
            f"{key} does not track terminal-error growth",
        )
    require(
        result["correlations"]["terminal_vs_integrated_source_error"] > 0.999,
        "terminal/source correlation is too weak",
    )
    controls = result["controls"]
    require(controls["maximum_qfp_drive_abs_log_ratio_dex"] < 0.005, "QFP drive differs")
    require(controls["maximum_source_weighted_alpha_abs_log_ratio_dex"] < 0.01, "alpha differs")
    require(controls["maximum_source_measure_relative_error"] < 1.0e-4, "source measure differs")
    require(controls["all_hotspot_cells_equal"], "cell hotspot moved")
    require(controls["all_hotspot_vertices_equal"], "vertex hotspot moved")
    require(
        max(float(row["operator_replay_x1e6_median_abs_error_dex"]) for row in edge)
        < 0.15,
        "unit-corrected operator replay does not close",
    )
    require(
        min(float(row["operator_replay_raw_median_abs_error_dex"]) for row in edge)
        > 5.9,
        "raw replay does not expose the expected 1e6 unit signature",
    )
    print(
        json.dumps(
            {
                "schema": "vela.pn2d_bv_m2_knee_readonly_verification.v1",
                "status": "passed",
                "bias_count": len(bias),
                "outcome": result["outcome"],
            },
            indent=2,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
