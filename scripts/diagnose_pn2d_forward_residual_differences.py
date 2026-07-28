#!/usr/bin/env python3
"""Quantify residual Vela/Sentaurus differences after the 2-D charge-area fix."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


ANCHORS = (1, 2, 5, 10, 15, 20)


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def rms(values: list[float]) -> float:
    return math.sqrt(sum(value * value for value in values) / len(values))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--audit-dir", type=Path, required=True)
    parser.add_argument("--sentaurus-fields", type=Path, required=True)
    parser.add_argument("--vela-iv", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    field_rows = rows(args.audit_dir / "node_field_comparison.csv")
    vela_rows = {
        round(float(row["bias_V"]), 9): row
        for row in rows(args.vela_iv)
        if row["converged"] == "1"
    }

    potential_decomposition = []
    for bias in (0, 1, 2, 5, 10, 15, 20):
        selected = [
            row for row in field_rows if abs(float(row["bias_V"]) - bias) < 1e-9
        ]
        epsi = [float(row["vela_minus_sent_psi"]) for row in selected]
        en = [float(row["vela_minus_sent_phin"]) for row in selected]
        ep = [float(row["vela_minus_sent_phip"]) for row in selected]
        common = [(a + b + c) / 3.0 for a, b, c in zip(epsi, en, ep)]
        electron_driver = [a - b for a, b in zip(epsi, en)]
        hole_driver = [c - a for a, c in zip(epsi, ep)]
        potential_decomposition.append(
            {
                "bias_V": bias,
                "psi_error_rms_V": rms(epsi),
                "common_mode_error_rms_V": rms(common),
                "electron_psi_minus_phin_error_rms_V": rms(electron_driver),
                "hole_phip_minus_psi_error_rms_V": rms(hole_driver),
                "common_mode_fraction_of_psi_rms": (
                    rms(common) / rms(epsi) if rms(epsi) else 0.0
                ),
            }
        )

    current_comparison = []
    for bias in ANCHORS:
        token = f"{bias}v"
        flux_path = (
            args.sentaurus_fields
            / token
            / "fields"
            / "ContactCurrentFlux_region2.csv"
        )
        sentaurus_current = abs(float(rows(flux_path)[0]["component0"]))
        vela_current = abs(
            float(vela_rows[float(bias)]["current_total_A_per_um"])
        )
        current_comparison.append(
            {
                "bias_V": bias,
                "sentaurus_A_per_um": sentaurus_current,
                "vela_A_per_um": vela_current,
                "vela_over_sentaurus": vela_current / sentaurus_current,
                "relative_difference_percent": (
                    100.0 * (vela_current - sentaurus_current) / sentaurus_current
                ),
            }
        )

    output = {
        "schema": "vela.pn2d.forward_residual_diagnostic.v1",
        "current_comparison_basis": (
            "Exact same-bias Sentaurus ContactCurrentFlux_region2 anchors; "
            "no interpolation of sparse PLT points."
        ),
        "current_comparison": current_comparison,
        "potential_error_decomposition": potential_decomposition,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(output, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
