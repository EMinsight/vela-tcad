#!/usr/bin/env python3
"""Infer the production SG edge mobility from each C++ baseline edge flux."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.build_pn2d_minimal6_cpp_mobility_current_input import triangle_mobility
from scripts.pn2d_minimal6_diagnostics.qfp_sg_experiment import (
    SILICON_NI_300K_M3,
    THERMAL_VOLTAGE_300K_V,
    _observation_index,
    _state,
)
from scripts.pn2d_minimal6_diagnostics.qfp_sg_replacement import qf_sg_flux


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--observations", type=Path, required=True)
    parser.add_argument("--current-edges", type=Path, required=True)
    parser.add_argument("--baseline-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    observations, _ = _observation_index(args.observations)
    state_cache: dict[tuple[str, float], dict[int, dict[str, float]]] = {}
    edge_cache: dict[tuple[str, float], dict[tuple[int, int], dict[str, str]]] = {}
    triangle_cache: dict[
        tuple[str, float], dict[tuple[int, int, str], float]
    ] = {}
    calibrated_count = 0
    fallback_count = 0
    affected_fallback_count = 0
    triangle_ratios: list[float] = []
    rows: list[dict[str, str]] = []

    with args.current_edges.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ValueError("current-edge input has no header")
        fieldnames = list(reader.fieldnames)
        for row in reader:
            topology, bias = row["topology"], float(row["bias_V"])
            state_key = (topology, bias)
            label = f"m{abs(int(bias))}V"
            if state_key not in state_cache:
                state_cache[state_key] = _state(
                    observations, "vela", topology, bias
                )
                edge_path = args.baseline_root / topology / label / "edges.csv"
                with edge_path.open(newline="", encoding="utf-8") as edge_handle:
                    edge_cache[state_key] = {
                        (int(value["node0"]), int(value["node1"])): value
                        for value in csv.DictReader(edge_handle)
                    }
                triangle_cache[state_key] = triangle_mobility(
                    args.baseline_root / topology / label / "triangles.csv"
                )
            state = state_cache[state_key]
            node0, node1 = int(row["node0"]), int(row["node1"])
            pair = tuple(sorted((node0, node1)))
            carrier = row["carrier"]
            qf_key = "phin_V" if carrier == "electron" else "phip_V"
            cpp_key = (
                "electron_raw_signed_flux_per_m2_s"
                if carrier == "electron"
                else "hole_raw_signed_flux_per_m2_s"
            )
            length = float(row["length_m"])
            unit_flux = qf_sg_flux(
                carrier,
                SILICON_NI_300K_M3,
                SILICON_NI_300K_M3,
                state[node0]["psi_V"],
                state[node1]["psi_V"],
                state[node0][qf_key],
                state[node1][qf_key],
                THERMAL_VOLTAGE_300K_V,
                THERMAL_VOLTAGE_300K_V / length,
            )
            cpp_flux = float(edge_cache[state_key][pair][cpp_key])
            triangle_value = triangle_cache[state_key][
                (pair[0], pair[1], carrier)
            ]
            if unit_flux != 0.0:
                mobility = cpp_flux / unit_flux
                if not math.isfinite(mobility) or mobility <= 0.0:
                    raise ValueError(
                        f"invalid inferred mobility for {state_key} {pair} {carrier}"
                    )
                calibrated_count += 1
                triangle_ratios.append(triangle_value / mobility)
            else:
                if cpp_flux != 0.0:
                    raise ValueError("zero unit flux has nonzero C++ baseline flux")
                mobility = triangle_value
                fallback_count += 1
                if node0 in (1, 5) or node1 in (1, 5):
                    affected_fallback_count += 1
            row["vela_masetti_native_state_mobility_m2_per_Vs"] = format(
                mobility, ".17g"
            )
            rows.append(row)

    if len(state_cache) != 40 or len(rows) != 720:
        raise ValueError(
            f"expected 40 states/720 rows, got {len(state_cache)}/{len(rows)}"
        )
    if affected_fallback_count:
        raise ValueError(
            f"{affected_fallback_count} replacement-affected edges need fallback mobility"
        )
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    manifest = {
        "schema_version": 1,
        "status": "valid",
        "semantics": (
            "production SG edge mobility inferred from baseline C++ edge flux "
            "using exact linearity in mobility; triangle mobility only for zero-QFP-flux controls"
        ),
        "state_count": len(state_cache),
        "row_count": len(rows),
        "calibrated_count": calibrated_count,
        "zero_flux_fallback_count": fallback_count,
        "replacement_affected_fallback_count": affected_fallback_count,
        "triangle_to_inferred_mobility_ratio_min": min(triangle_ratios),
        "triangle_to_inferred_mobility_ratio_max": max(triangle_ratios),
        "observations": str(args.observations.resolve()),
        "observations_sha256": sha256(args.observations),
        "input_current_edges": str(args.current_edges.resolve()),
        "input_current_edges_sha256": sha256(args.current_edges),
        "baseline_root": str(args.baseline_root.resolve()),
        "output": str(args.output.resolve()),
        "output_sha256": sha256(args.output),
    }
    args.output.with_suffix(".manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
