#!/usr/bin/env python3
"""Factor Minimal6 avalanche-source differences on one triangle-edge support."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from pathlib import Path
from typing import Any


Q_C = 1.602176634e-19
TOPOLOGIES = ("mirror", "sketch")
BIASES = tuple(range(1, 21))
CARRIERS = ("electron", "hole")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def read_scalar(path: Path) -> dict[int, float]:
    return {
        int(row["node_id"]): float(row["component0"])
        for row in read_csv(path)
    }


def read_vector(path: Path) -> dict[int, tuple[float, float]]:
    return {
        int(row["node_id"]): (
            float(row["component0"]),
            float(row["component1"]),
        )
        for row in read_csv(path)
    }


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def log_error(value: float, reference: float) -> float | None:
    if value <= 0.0 or reference <= 0.0:
        return None
    return abs(math.log10(value / reference))


def percentile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * fraction
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    return ordered[lower] + (
        ordered[upper] - ordered[lower]
    ) * (position - lower)


def summary(values: list[float]) -> dict[str, float | int]:
    return {
        "count": len(values),
        "median": statistics.median(values),
        "p95": percentile(values, 0.95),
        "maximum": max(values),
    }


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table {path}")
    with path.open("w", newline="", encoding="ascii") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def local_field(row: dict[str, str], local: int, name: str) -> float:
    return float(row[f"local_edge{local}_{name}"])


def local_node(row: dict[str, str], local: int, endpoint: int) -> int:
    return int(row[f"local_edge{local}_node{endpoint}"])


def run(
    *,
    phase_f_root: Path,
    candidate_sweep_root: Path,
    inverse_inputs_root: Path,
    output_root: Path,
) -> dict[str, Any]:
    phase_f_root = phase_f_root.resolve()
    candidate_sweep_root = candidate_sweep_root.resolve()
    inverse_inputs_root = inverse_inputs_root.resolve()
    output_root = output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    local_rows: list[dict[str, Any]] = []
    state_rows: list[dict[str, Any]] = []
    maximum_geometry_relative_error = 0.0
    maximum_node_mapping_relative_error = 0.0
    zero_volume_count = 0
    zero_volume_nonzero_source_count = 0

    for topology in TOPOLOGIES:
        for bias in BIASES:
            raw = (
                phase_f_root / "raw" / topology / f"m{bias}V"
                / "triangle.csv"
            )
            triangle_rows = read_csv(raw)
            if len(triangle_rows) != 4:
                raise ValueError(
                    f"{topology} -{bias} V requires four triangle rows"
                )
            fields = (
                inverse_inputs_root / "sentaurus" / "source" / "reimport"
                / topology / f"m{bias}V" / "fields"
            )
            sent_alpha = {
                "electron": read_scalar(
                    fields / "eAlphaAvalanche_region0.csv"
                ),
                "hole": read_scalar(
                    fields / "hAlphaAvalanche_region0.csv"
                ),
            }
            sent_current = {
                "electron": read_vector(
                    fields / "eCurrentDensity_region0.csv"
                ),
                "hole": read_vector(
                    fields / "hCurrentDensity_region0.csv"
                ),
            }
            sent_impact = read_scalar(
                fields / "ImpactIonization_region0.csv"
            )
            sent_electric = read_vector(
                fields / "ElectricField_region0.csv"
            )

            sent_particle_flux: dict[str, dict[int, float]] = {}
            sent_generation: dict[str, dict[int, float]] = {}
            for carrier in CARRIERS:
                sent_particle_flux[carrier] = {
                    node: math.hypot(*vector) / Q_C
                    for node, vector in sent_current[carrier].items()
                }
                sent_generation[carrier] = {
                    node: sent_alpha[carrier][node]
                    * sent_particle_flux[carrier][node]
                    for node in sent_alpha[carrier]
                }

            native_source = 0.0
            nodal_reconstructed_source = 0.0
            projected_sentaurus_source = 0.0
            vela_source = 0.0
            sentaurus_alpha_hybrid_source = 0.0
            sentaurus_current_hybrid_source = 0.0
            node_mapped_source = {node: 0.0 for node in sent_impact}

            for triangle in triangle_rows:
                nodes = [
                    int(triangle["node0"]),
                    int(triangle["node1"]),
                    int(triangle["node2"]),
                ]
                area_m2 = 0.5 * abs(
                    float(triangle["signed_double_area_m2"])
                )
                area_cm2 = area_m2 * 1.0e4
                native_source += (
                    area_cm2
                    * sum(sent_impact[node] for node in nodes)
                    / 3.0
                )
                nodal_reconstructed_source += area_cm2 * sum(
                    sent_generation[carrier][node]
                    for carrier in CARRIERS
                    for node in nodes
                ) / 3.0

                triangle_volume = 0.0
                candidate_electric_field = math.hypot(
                    float(triangle["grad_psi_x_V_per_m"]),
                    float(triangle["grad_psi_y_V_per_m"]),
                )
                for local in range(3):
                    node0 = local_node(triangle, local, 0)
                    node1 = local_node(triangle, local, 1)
                    volume_m2 = local_field(
                        triangle, local, "truncated_partial_volume_m2"
                    )
                    volume_cm2 = volume_m2 * 1.0e4
                    triangle_volume += volume_m2
                    if volume_m2 == 0.0:
                        zero_volume_count += 1

                    sent_electric_field = 0.5 * (
                        math.hypot(*sent_electric[node0])
                        + math.hypot(*sent_electric[node1])
                    ) * 100.0
                    candidate_total = 0.0
                    sent_alpha_total = 0.0
                    sent_current_total = 0.0
                    projected_sent_total = 0.0

                    for carrier in CARRIERS:
                        prefix = carrier
                        candidate_alpha_m_inv = local_field(
                            triangle, local, f"{prefix}_alpha_per_m"
                        )
                        candidate_flux_m2 = abs(
                            local_field(
                                triangle,
                                local,
                                f"{prefix}_flux_proxy_per_m2_s",
                            )
                        )
                        candidate_source_per_m = local_field(
                            triangle,
                            local,
                            f"{prefix}_source_integral_per_m_s",
                        )
                        sent_alpha_cm_inv = 0.5 * (
                            sent_alpha[carrier][node0]
                            + sent_alpha[carrier][node1]
                        )
                        sent_flux_cm2 = 0.5 * (
                            sent_particle_flux[carrier][node0]
                            + sent_particle_flux[carrier][node1]
                        )
                        sent_generation_cm3 = 0.5 * (
                            sent_generation[carrier][node0]
                            + sent_generation[carrier][node1]
                        )
                        candidate_alpha_cm_inv = (
                            candidate_alpha_m_inv / 100.0
                        )
                        candidate_flux_cm2 = (
                            candidate_flux_m2 / 1.0e4
                        )

                        candidate_source = candidate_source_per_m * 1.0e-2
                        sent_alpha_hybrid = (
                            sent_alpha_cm_inv
                            * candidate_flux_cm2
                            * volume_cm2
                        )
                        sent_current_hybrid = (
                            candidate_alpha_cm_inv
                            * sent_flux_cm2
                            * volume_cm2
                        )
                        projected_sent = (
                            sent_generation_cm3 * volume_cm2
                        )
                        if volume_m2 == 0.0 and any(
                            value != 0.0
                            for value in (
                                candidate_source,
                                sent_alpha_hybrid,
                                sent_current_hybrid,
                                projected_sent,
                            )
                        ):
                            zero_volume_nonzero_source_count += 1

                        candidate_total += candidate_source
                        sent_alpha_total += sent_alpha_hybrid
                        sent_current_total += sent_current_hybrid
                        projected_sent_total += projected_sent

                        alpha_error = log_error(
                            candidate_alpha_m_inv,
                            sent_alpha_cm_inv * 100.0,
                        )
                        current_error = log_error(
                            candidate_flux_cm2,
                            sent_flux_cm2,
                        )
                        qfp_field = local_field(
                            triangle,
                            local,
                            f"{prefix}_cell_qf_field_V_per_m",
                        )
                        local_rows.append(
                            {
                                "topology": topology,
                                "bias_V": -bias,
                                "cell_id": int(triangle["cell_id"]),
                                "local_edge": local,
                                "node0": node0,
                                "node1": node1,
                                "carrier": carrier,
                                "volume_m2": volume_m2,
                                "geometric_zero": int(volume_m2 == 0.0),
                                "candidate_electric_field_V_per_m":
                                    candidate_electric_field,
                                "candidate_qfp_field_V_per_m": qfp_field,
                                "sentaurus_electric_field_V_per_m":
                                    sent_electric_field,
                                "candidate_alpha_per_m":
                                    candidate_alpha_m_inv,
                                "sentaurus_endpoint_alpha_per_m":
                                    sent_alpha_cm_inv * 100.0,
                                "alpha_abs_error_dex":
                                    "" if alpha_error is None else alpha_error,
                                "candidate_flux_proxy_per_cm2_s":
                                    candidate_flux_cm2,
                                "sentaurus_endpoint_flux_per_cm2_s":
                                    sent_flux_cm2,
                                "current_abs_error_dex":
                                    "" if current_error is None else current_error,
                                "candidate_source_per_cm_s":
                                    candidate_source,
                                "sentaurus_alpha_hybrid_source_per_cm_s":
                                    sent_alpha_hybrid,
                                "sentaurus_current_hybrid_source_per_cm_s":
                                    sent_current_hybrid,
                                "projected_sentaurus_source_per_cm_s":
                                    projected_sent,
                            }
                        )

                    vela_source += candidate_total
                    sentaurus_alpha_hybrid_source += sent_alpha_total
                    sentaurus_current_hybrid_source += sent_current_total
                    projected_sentaurus_source += projected_sent_total
                    node_mapped_source[node0] += 0.5 * candidate_total
                    node_mapped_source[node1] += 0.5 * candidate_total

                geometry_denominator = max(area_m2, 1.0e-300)
                maximum_geometry_relative_error = max(
                    maximum_geometry_relative_error,
                    abs(triangle_volume - area_m2) / geometry_denominator,
                )

            mapped_total = sum(node_mapped_source.values())
            mapping_denominator = max(abs(vela_source), 1.0e-300)
            mapping_error = abs(mapped_total - vela_source) / mapping_denominator
            maximum_node_mapping_relative_error = max(
                maximum_node_mapping_relative_error, mapping_error
            )
            state_rows.append(
                {
                    "topology": topology,
                    "bias_V": -bias,
                    "sentaurus_native_source_per_cm_s": native_source,
                    "sentaurus_nodal_alpha_current_source_per_cm_s":
                        nodal_reconstructed_source,
                    "sentaurus_projected_triangle_source_per_cm_s":
                        projected_sentaurus_source,
                    "vela_candidate_source_per_cm_s": vela_source,
                    "sentaurus_alpha_hybrid_source_per_cm_s":
                        sentaurus_alpha_hybrid_source,
                    "sentaurus_current_hybrid_source_per_cm_s":
                        sentaurus_current_hybrid_source,
                    "node_mapped_candidate_source_per_cm_s": mapped_total,
                    "nodal_reconstruction_abs_error_dex":
                        log_error(nodal_reconstructed_source, native_source),
                    "triangle_projection_abs_error_dex":
                        log_error(
                            projected_sentaurus_source,
                            nodal_reconstructed_source,
                        ),
                    "candidate_source_abs_error_dex":
                        log_error(vela_source, native_source),
                    "sentaurus_alpha_hybrid_abs_error_dex":
                        log_error(
                            sentaurus_alpha_hybrid_source, native_source
                        ),
                    "sentaurus_current_hybrid_abs_error_dex":
                        log_error(
                            sentaurus_current_hybrid_source, native_source
                        ),
                    "projected_sentaurus_abs_error_dex":
                        log_error(projected_sentaurus_source, native_source),
                    "node_mapping_relative_error": mapping_error,
                }
            )

    if len(local_rows) != 960 or len(state_rows) != 40:
        raise ValueError("impact factorization row-count contract failed")

    metric_columns = {
        "alpha": "alpha_abs_error_dex",
        "current": "current_abs_error_dex",
        "nodal_reconstruction": "nodal_reconstruction_abs_error_dex",
        "triangle_projection": "triangle_projection_abs_error_dex",
        "candidate_source": "candidate_source_abs_error_dex",
        "sentaurus_alpha_hybrid": "sentaurus_alpha_hybrid_abs_error_dex",
        "sentaurus_current_hybrid": "sentaurus_current_hybrid_abs_error_dex",
        "projected_sentaurus": "projected_sentaurus_abs_error_dex",
    }
    summaries: dict[str, Any] = {}
    for name, column in metric_columns.items():
        source = local_rows if column in local_rows[0] else state_rows
        values = [
            float(row[column])
            for row in source
            if row[column] not in ("", None)
            and (name not in {"alpha", "current"} or not row["geometric_zero"])
        ]
        summaries[name] = summary(values)

    outcome = "current_support_dominant"
    if summaries["alpha"]["median"] > summaries["current"]["median"]:
        outcome = "alpha_difference_dominant"
    if summaries["triangle_projection"]["p95"] > 0.1:
        outcome = "support_projection_material"

    write_csv(output_root / "local_edge_factorization.csv", local_rows)
    write_csv(output_root / "state_source_factorization.csv", state_rows)
    report = [
        "# PN2D Minimal6 impact factorization",
        "",
        f"Typed outcome: `{outcome}`.",
        "",
        "The comparison uses 40 exact states, four triangles per state, "
        "and three local edges per triangle. Sentaurus directed edges remain "
        "a projection, not a native edge observation.",
        "",
        f"- zero-volume local edges: `{zero_volume_count}`",
        f"- zero-volume edges with nonzero source: "
        f"`{zero_volume_nonzero_source_count}`",
        f"- maximum geometry closure error: "
        f"`{maximum_geometry_relative_error:.17g}`",
        f"- maximum source-to-node mapping error: "
        f"`{maximum_node_mapping_relative_error:.17g}`",
        "",
        "## Metric medians",
        "",
        "| metric | count | median dex | P95 dex | maximum dex |",
        "|---|---:|---:|---:|---:|",
    ]
    for name, values in summaries.items():
        report.append(
            f"| {name} | {values['count']} | {values['median']:.9g} | "
            f"{values['p95']:.9g} | {values['maximum']:.9g} |"
        )
    report.extend(
        [
            "",
            "No coefficient, field scale, saturation velocity, beta, "
            "or geometric weight is fitted.",
        ]
    )
    (output_root / "report.md").write_text(
        "\n".join(report) + "\n", encoding="ascii", newline="\n"
    )
    manifest = {
        "schema_version": 1,
        "experiment": "pn2d_minimal6_impact_factorization",
        "status": "valid",
        "outcome": outcome,
        "contracts": {
            "state_count": len(state_rows),
            "local_carrier_edge_count": len(local_rows),
            "zero_volume_count": zero_volume_count,
            "zero_volume_nonzero_source_count":
                zero_volume_nonzero_source_count,
            "maximum_geometry_relative_error":
                maximum_geometry_relative_error,
            "maximum_node_mapping_relative_error":
                maximum_node_mapping_relative_error,
            "sentaurus_current_support": "endpoint_to_triangle_projection",
        },
        "summaries": summaries,
        "inputs": {
            "phase_f_manifest_sha256":
                sha256(phase_f_root / "manifest.json"),
            "candidate_sweep_manifest_sha256":
                sha256(candidate_sweep_root / "sweep_manifest.json"),
            "inverse_sentaurus_manifest_sha256":
                sha256(
                    inverse_inputs_root / "sentaurus" / "manifest.json"
                ),
        },
    }
    for name in (
        "local_edge_factorization.csv",
        "state_source_factorization.csv",
        "report.md",
    ):
        manifest.setdefault("outputs", {})[name] = sha256(output_root / name)
    (output_root / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="ascii",
        newline="\n",
    )
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase-f-root", type=Path, required=True)
    parser.add_argument("--candidate-sweep-root", type=Path, required=True)
    parser.add_argument("--inverse-inputs-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()
    result = run(
        phase_f_root=args.phase_f_root,
        candidate_sweep_root=args.candidate_sweep_root,
        inverse_inputs_root=args.inverse_inputs_root,
        output_root=args.output_root,
    )
    print(json.dumps(
        {"status": result["status"], "outcome": result["outcome"]},
        sort_keys=True,
    ))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
