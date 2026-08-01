#!/usr/bin/env python3
"""Verify the frozen M2 carrier-QFP residual/Jacobian diagnostic artifact."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    args = parser.parse_args()

    result = json.loads((args.root / "result.json").read_text(encoding="utf-8"))
    require(
        result["schema"] == "vela.pn2d_bv_m2_qfp_carrier_jacobian_verification.v1",
        "unexpected result schema",
    )
    require(result["status"] == "passed", "experiment did not finish")
    require(result["biases_V"] == [-18.0, -19.5, -19.7, -20.0], "bias lattice changed")
    require(not result["physics_modified"], "physics was modified")
    require(not result["production_defaults_modified"], "production defaults changed")
    require(not result["acceptance_thresholds_modified"], "acceptance thresholds changed")
    require(result["determinism"]["all_byte_identical"], "independent runs differ")
    require(result["determinism"]["artifact_count"] == 148, "artifact count changed")
    require(result["maximum_term_closure"] <= 1.0e-20, "carrier-term closure failed")

    verdict = result["verdict"]
    require(verdict["source_outcome"] == "phip_dominant", "QFP carrier attribution changed")
    require(
        verdict["carrier_term_outcomes"]
        == {"electron": "flux_dominant", "hole": "flux_dominant"},
        "carrier residual-term attribution changed",
    )
    require(
        verdict["jacobian_fd_outcome"] == "analytic_fd_inconsistent",
        "formal predeclared Jacobian outcome changed",
    )
    require(
        verdict["jacobian_fd_interpretation"]
        == "formal_relative_gate_fails_only_at_srh_absolute_fd_floor",
        "finite-difference floor interpretation changed",
    )
    require(
        float(verdict["worst_non_srh_jacobian_fd_check"]["relative_difference"])
        <= float(result["finite_difference_threshold"]),
        "a non-SRH Jacobian block failed the predeclared threshold",
    )
    require(
        verdict["srh_fd_step_sensitivity"]["classified_as_absolute_fd_floor"],
        "SRH step sensitivity no longer supports an absolute finite-difference floor",
    )
    require(
        verdict["qfp_update_outcome"] == "both_qfp_updates_roll_back_from_sentaurus",
        "first-update direction changed",
    )

    expected_rows = {
        "source_carrier_substitution.csv": 16,
        "carrier_term_decomposition.csv": 480,
        "carrier_term_closure.csv": 32,
        "first_qfp_updates.csv": 32,
        "jacobian_fd_blocks.csv": 40,
        "jacobian_fd_step_sensitivity.csv": 28,
        "determinism.csv": 148,
    }
    for name, expected in expected_rows.items():
        actual = len(read_rows(args.root / name))
        require(actual == expected, f"{name}: expected {expected} rows, found {actual}")

    print(json.dumps({
        "status": "passed",
        "typed_outcome": verdict["typed_outcome"],
        "interpretation": verdict["jacobian_fd_interpretation"],
        "artifact_count": result["determinism"]["artifact_count"],
    }, indent=2))


if __name__ == "__main__":
    main()
