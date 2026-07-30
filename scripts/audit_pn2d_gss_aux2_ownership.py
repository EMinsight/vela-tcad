#!/usr/bin/env python3
"""Audit the archived GSS aux2 midpoint ownership without changing a solver."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any


Q_C = 1.602176634e-19
KB_J_PER_K = 1.380649e-23


def aux2(value: float) -> float:
    if value >= 0.0:
        return math.exp(-value) / (1.0 + math.exp(-value))
    return 1.0 / (1.0 + math.exp(value))


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def write_rows(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise RuntimeError(f"refusing to write empty audit CSV: {path}")
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def relative_error(actual: float, expected: float) -> float:
    return abs(actual - expected) / max(abs(expected), 1.0e-300)


def relative_l2(pairs: list[tuple[float, float]]) -> float:
    numerator = sum((actual - expected) ** 2 for actual, expected in pairs)
    denominator = sum(expected**2 for _, expected in pairs)
    return math.sqrt(numerator / denominator) if denominator > 0.0 else 0.0


def reference_midpoint(
    carrier: str,
    density0: float,
    density1: float,
    psi0: float,
    psi1: float,
    thermal_voltage: float,
) -> float:
    alpha = (psi0 - psi1) / (2.0 * thermal_voltage)
    if carrier == "electron":
        return density0 * aux2(alpha) + density1 * aux2(-alpha)
    if carrier == "hole":
        return density0 * aux2(-alpha) + density1 * aux2(alpha)
    raise ValueError(f"unsupported carrier: {carrier}")


def production_midpoint(
    carrier: str,
    density0: float,
    density1: float,
    psi0: float,
    psi1: float,
    thermal_voltage: float,
) -> float:
    alpha = (psi0 - psi1) / (2.0 * thermal_voltage)
    if carrier == "electron":
        return density0 * aux2(-alpha) + density1 * aux2(alpha)
    if carrier == "hole":
        return density0 * aux2(alpha) + density1 * aux2(-alpha)
    raise ValueError(f"unsupported carrier: {carrier}")


def key(row: dict[str, str]) -> tuple[str, int, int, int, int, int]:
    return (
        row["carrier"],
        int(row["cell_id"]),
        int(row["local_edge"]),
        int(row["edge_id"]),
        int(row["node0"]),
        int(row["node1"]),
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--temperature-K", type=float, default=300.0)
    parser.add_argument("--gss-guide-text", type=Path, required=True)
    parser.add_argument("--gss-jflux-header", type=Path, required=True)
    args = parser.parse_args()

    if args.temperature_K <= 0.0:
        raise ValueError("--temperature-K must be positive")
    thermal_voltage = KB_J_PER_K * args.temperature_K / Q_C
    result_path = args.root / "result.json"
    result = json.loads(result_path.read_text(encoding="utf-8-sig"))
    if not result.get("observation_only", False):
        raise RuntimeError("input comparison is not observation-only")
    factorization_path = Path(
        result["artifacts"]["current_proxy_factorization"]["path"]
    )
    factor_rows = read_rows(factorization_path)
    factor_by_bias: dict[float, list[dict[str, str]]] = {}
    for row in factor_rows:
        factor_by_bias.setdefault(float(row["bias_V"]), []).append(row)

    details: list[dict[str, Any]] = []
    summaries: list[dict[str, Any]] = []
    for case in result["cases"]:
        bias = float(case["bias_V"])
        process_path = Path(case["outputs"]["process"]["path"])
        process_by_key = {key(row): row for row in read_rows(process_path)}
        reference_edge_pairs: list[tuple[float, float]] = []
        recorded_production_pairs: list[tuple[float, float]] = []
        swap_errors: list[float] = []
        reference_source_total = 0.0
        carrier_rows: dict[str, list[dict[str, Any]]] = {
            "electron": [],
            "hole": [],
        }

        for factor in factor_by_bias[bias]:
            process = process_by_key[key(factor)]
            carrier = factor["carrier"]
            density0 = float(factor["density0_m3"])
            density1 = float(factor["density1_m3"])
            psi0 = float(process["psi0_V"])
            psi1 = float(process["psi1_V"])
            archived_midpoint = reference_midpoint(
                carrier, density0, density1, psi0, psi1, thermal_voltage
            )
            archived_swapped = reference_midpoint(
                carrier, density1, density0, psi1, psi0, thermal_voltage
            )
            reconstructed_production = production_midpoint(
                carrier, density0, density1, psi0, psi1, thermal_voltage
            )
            recorded_production = float(factor["gss_midpoint_density_m3"])
            edge_midpoint = float(
                factor["vela_edge_audit_midpoint_density_m3"]
            )
            mobility = float(factor["mobility_m2_per_V_s"])
            qf_drive = float(factor["edge_qf_drive_V_per_m"])
            archived_proxy_current = (
                Q_C * mobility * archived_midpoint * qf_drive
            )
            alpha = float(process["alpha_per_m"])
            source_measure = float(process["source_measure_m2"])
            archived_source = (
                alpha * archived_proxy_current / Q_C * source_measure
            )
            reference_source_total += archived_source
            row = {
                "bias_V": bias,
                "carrier": carrier,
                "cell_id": int(factor["cell_id"]),
                "local_edge": int(factor["local_edge"]),
                "edge_id": int(factor["edge_id"]),
                "node0": int(factor["node0"]),
                "node1": int(factor["node1"]),
                "psi0_V": psi0,
                "psi1_V": psi1,
                "density0_m3": density0,
                "density1_m3": density1,
                "production_midpoint_density_m3": recorded_production,
                "production_formula_reconstructed_midpoint_density_m3": (
                    reconstructed_production
                ),
                "production_formula_closure_relative_error": relative_error(
                    reconstructed_production, recorded_production
                ),
                "gss_reference_midpoint_density_m3": archived_midpoint,
                "gss_reference_endpoint_swap_relative_error": relative_error(
                    archived_swapped, archived_midpoint
                ),
                "edge_audit_midpoint_density_m3": edge_midpoint,
                "gss_reference_to_edge_midpoint_relative_error": relative_error(
                    archived_midpoint, edge_midpoint
                ),
                "production_over_gss_reference_midpoint": (
                    recorded_production / archived_midpoint
                    if archived_midpoint > 0.0
                    else math.inf
                ),
                "production_proxy_current_A_per_m2": float(
                    factor["production_proxy_current_A_per_m2"]
                ),
                "gss_reference_midpoint_proxy_current_A_per_m2": (
                    archived_proxy_current
                ),
                "raw_sg_transport_current_A_per_m2": float(
                    factor["vela_raw_sg_transport_current_A_per_m2"]
                ),
                "sentaurus_current_magnitude_A_per_m2": float(
                    factor["sentaurus_current_magnitude_A_per_m2"]
                ),
                "gss_reference_midpoint_source_integral_per_m_s": (
                    archived_source
                ),
            }
            details.append(row)
            carrier_rows[carrier].append(row)
            reference_edge_pairs.append((archived_midpoint, edge_midpoint))
            recorded_production_pairs.append(
                (reconstructed_production, recorded_production)
            )
            swap_errors.append(
                row["gss_reference_endpoint_swap_relative_error"]
            )

        sentaurus_total = float(
            case["sentaurus_total_source_integral_per_m_s"]
        )
        for carrier, rows in carrier_rows.items():
            summaries.append(
                {
                    "bias_V": bias,
                    "carrier": carrier,
                    "record_count": len(rows),
                    "production_formula_closure_relative_l2": relative_l2(
                        [
                            (
                                float(row[
                                    "production_formula_reconstructed_midpoint_density_m3"
                                ]),
                                float(row[
                                    "production_midpoint_density_m3"
                                ]),
                            )
                            for row in rows
                        ]
                    ),
                    "gss_reference_to_edge_midpoint_relative_l2": relative_l2(
                        [
                            (
                                float(row[
                                    "gss_reference_midpoint_density_m3"
                                ]),
                                float(row["edge_audit_midpoint_density_m3"]),
                            )
                            for row in rows
                        ]
                    ),
                    "maximum_endpoint_swap_relative_error": max(
                        float(row[
                            "gss_reference_endpoint_swap_relative_error"
                        ])
                        for row in rows
                    ),
                }
            )
        hotspot = max(
            (row for row in details if float(row["bias_V"]) == bias),
            key=lambda row: float(row[
                "production_over_gss_reference_midpoint"
            ]),
        )
        case["gss_aux2_ownership_audit"] = {
            "gss_reference_midpoint_proxy_total_source_integral_per_m_s": (
                reference_source_total
            ),
            "gss_reference_midpoint_proxy_source_ratio_to_sentaurus": (
                reference_source_total / sentaurus_total
            ),
            "actual_sg_vector_source_ratio_to_sentaurus": (
                float(case["sg_vector_total_source_integral_per_m_s"])
                / sentaurus_total
            ),
            "maximum_reference_endpoint_swap_relative_error": max(
                swap_errors
            ),
            "gss_reference_to_edge_midpoint_relative_l2": relative_l2(
                reference_edge_pairs
            ),
            "production_formula_closure_relative_l2": relative_l2(
                recorded_production_pairs
            ),
            "hotspot": hotspot,
        }

    args.out_dir.mkdir(parents=True, exist_ok=True)
    detail_path = args.out_dir / "gss_aux2_ownership_details.csv"
    summary_path = args.out_dir / "gss_aux2_ownership_summary.csv"
    write_rows(detail_path, details)
    write_rows(summary_path, summaries)
    audit = {
        "schema": "vela.gss_aux2_ownership_audit.v1",
        "observation_only": True,
        "solver_state_advanced": False,
        "continuity_feedback_enabled": False,
        "temperature_K": args.temperature_K,
        "thermal_voltage_V": thermal_voltage,
        "reference_interpretation": {
            "aux2": "1/(1+exp(x))",
            "electron_alpha_isothermal": "(psi_i-psi_j)/(2*Vt)",
            "hole_alpha_isothermal": "(psi_i-psi_j)/(2*Vt)",
            "electron_midpoint": "n_i*aux2(alpha)+n_j*aux2(-alpha)",
            "hole_midpoint": "p_i*aux2(-alpha)+p_j*aux2(alpha)",
            "ownership": (
                "midpoint density inside the complete SG drift-diffusion "
                "decomposition; not a standalone avalanche current"
            ),
        },
        "inputs": {
            "frozen_comparison_result": {
                "path": str(result_path.resolve()),
                "sha256": sha256(result_path),
            },
            "gss_guide_text": {
                "path": str(args.gss_guide_text.resolve()),
                "sha256": sha256(args.gss_guide_text),
                "equations": ["9.100", "9.103", "9.107", "9.108"],
            },
            "gss_jflux_header": {
                "path": str(args.gss_jflux_header.resolve()),
                "sha256": sha256(args.gss_jflux_header),
                "functions": ["nmid", "pmid", "In", "Ip"],
            },
        },
        "cases": [
            {
                "bias_V": float(case["bias_V"]),
                **case["gss_aux2_ownership_audit"],
            }
            for case in result["cases"]
        ],
        "artifacts": {
            "details": {
                "path": str(detail_path.resolve()),
                "sha256": sha256(detail_path),
            },
            "summary": {
                "path": str(summary_path.resolve()),
                "sha256": sha256(summary_path),
            },
        },
    }
    audit_path = args.out_dir / "result.json"
    audit_path.write_text(json.dumps(audit, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
