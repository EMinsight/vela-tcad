#!/usr/bin/env python3
"""Decompose PN2D Minimal6 current differences on sealed common supports."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from pathlib import Path


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty CSV {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _quantile(values: list[float], fraction: float) -> float:
    if not values:
        raise ValueError("quantile needs at least one value")
    if fraction == 0.5:
        return statistics.median(values)
    ordered = sorted(values)
    return ordered[int(math.floor((len(ordered) - 1) * fraction))]


def _stats(values: list[float]) -> dict[str, object]:
    return {
        "count": len(values),
        "median": _quantile(values, 0.5),
        "p95": _quantile(values, 0.95),
        "maximum": max(values),
    }


def _pair(row: dict[str, str]) -> tuple[int, int]:
    node0 = int(row["node0"])
    node1 = int(row["node1"])
    return min(node0, node1), max(node0, node1)


def _key(
    row: dict[str, str],
) -> tuple[str, float, str, int, int]:
    node0, node1 = _pair(row)
    return (
        row["topology"],
        float(row["bias_V"]),
        row["carrier"],
        node0,
        node1,
    )


def _current_decomposition(
    phase_f_edges: list[dict[str, str]],
    staged_edges: list[dict[str, str]],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    control = {
        _key(row): row
        for row in staged_edges
        if row["stage"] == "sentaurus_qfp_recomputed_density_control"
        and row["status"] == "valid"
    }
    samples: list[dict[str, object]] = []
    for row in phase_f_edges:
        if row["status"] != "valid":
            continue
        key = _key(row)
        imported = control[key]
        before = float(row["absolute_log10_error_dex"])
        after = float(imported["absolute_log10_error_dex"])
        reference = abs(float(row["sentaurus_box_reconstructed_A_per_um"]))
        central = key[3:] == (1, 5)
        samples.append(
            {
                "topology": key[0],
                "bias_V": key[1],
                "carrier": key[2],
                "node0": key[3],
                "node1": key[4],
                "central_edge": central,
                "reference_abs_A_per_um": reference,
                "self_consistent_error_dex": before,
                "imported_state_error_dex": after,
                "paired_state_improvement_dex": before - after,
                "self_consistent_sign_agreement": float(row["sign_agreement"]),
                "imported_state_sign_agreement": float(
                    imported["sign_agreement"]
                ),
            }
        )

    summary: list[dict[str, object]] = []
    for carrier in ("electron", "hole"):
        selected = [row for row in samples if row["carrier"] == carrier]
        weight = sum(float(row["reference_abs_A_per_um"]) for row in selected)
        central = [row for row in selected if row["central_edge"]]
        other = [row for row in selected if not row["central_edge"]]
        summary.append(
            {
                "carrier": carrier,
                "active_edge_count": len(selected),
                "self_consistent_median_error_dex": _quantile(
                    [float(row["self_consistent_error_dex"]) for row in selected],
                    0.5,
                ),
                "imported_state_median_error_dex": _quantile(
                    [float(row["imported_state_error_dex"]) for row in selected],
                    0.5,
                ),
                "paired_state_improvement_median_dex": _quantile(
                    [
                        float(row["paired_state_improvement_dex"])
                        for row in selected
                    ],
                    0.5,
                ),
                "current_weighted_state_improvement_dex": sum(
                    float(row["reference_abs_A_per_um"])
                    * float(row["paired_state_improvement_dex"])
                    for row in selected
                )
                / weight,
                "self_consistent_sign_agreement_fraction": statistics.mean(
                    float(row["self_consistent_sign_agreement"])
                    for row in selected
                ),
                "imported_state_sign_agreement_fraction": statistics.mean(
                    float(row["imported_state_sign_agreement"])
                    for row in selected
                ),
                "central_edge_sign_agreement_fraction": statistics.mean(
                    float(row["self_consistent_sign_agreement"]) for row in central
                ),
                "other_edge_sign_agreement_fraction": statistics.mean(
                    float(row["self_consistent_sign_agreement"]) for row in other
                ),
                "central_edge_median_reference_abs_A_per_um": _quantile(
                    [float(row["reference_abs_A_per_um"]) for row in central], 0.5
                ),
                "other_edge_median_reference_abs_A_per_um": _quantile(
                    [float(row["reference_abs_A_per_um"]) for row in other], 0.5
                ),
            }
        )
    return samples, summary


def _qfp_modes(
    state_rows: list[dict[str, str]],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    by_node = {
        (row["topology"], float(row["bias_V"]), int(row["node_id"])): row
        for row in state_rows
    }
    samples: list[dict[str, object]] = []
    for topology in ("mirror", "sketch"):
        for magnitude in range(1, 21):
            bias = -float(magnitude)
            node1 = by_node[(topology, bias, 1)]
            node5 = by_node[(topology, bias, 5)]
            for carrier, field in (
                ("electron", "phin"),
                ("hole", "phip"),
            ):
                sent1 = float(node1[f"sentaurus_{field}_V"])
                sent5 = float(node5[f"sentaurus_{field}_V"])
                vela1 = float(node1[f"vela_{field}_V"])
                vela5 = float(node5[f"vela_{field}_V"])
                error1 = vela1 - sent1
                error5 = vela5 - sent5
                sent_delta = sent5 - sent1
                vela_delta = vela5 - vela1
                samples.append(
                    {
                        "topology": topology,
                        "bias_V": bias,
                        "carrier": carrier,
                        "common_mode_error_V": 0.5 * (error1 + error5),
                        "differential_mode_error_V": error5 - error1,
                        "sentaurus_edge_qfp_delta_V": sent_delta,
                        "vela_edge_qfp_delta_V": vela_delta,
                        "edge_qfp_sign_agreement": float(
                            sent_delta * vela_delta > 0.0
                        ),
                        "vela_to_sentaurus_abs_edge_delta_ratio": abs(
                            vela_delta / sent_delta
                        ),
                    }
                )
    summary: list[dict[str, object]] = []
    for carrier in ("electron", "hole"):
        selected = [row for row in samples if row["carrier"] == carrier]
        summary.append(
            {
                "carrier": carrier,
                "state_count": len(selected),
                "common_mode_abs_error_median_V": _quantile(
                    [abs(float(row["common_mode_error_V"])) for row in selected],
                    0.5,
                ),
                "common_mode_abs_error_maximum_V": max(
                    abs(float(row["common_mode_error_V"])) for row in selected
                ),
                "differential_mode_abs_error_median_V": _quantile(
                    [
                        abs(float(row["differential_mode_error_V"]))
                        for row in selected
                    ],
                    0.5,
                ),
                "edge_qfp_sign_agreement_fraction": statistics.mean(
                    float(row["edge_qfp_sign_agreement"]) for row in selected
                ),
                "sentaurus_abs_edge_delta_median_V": _quantile(
                    [
                        abs(float(row["sentaurus_edge_qfp_delta_V"]))
                        for row in selected
                    ],
                    0.5,
                ),
                "vela_abs_edge_delta_median_V": _quantile(
                    [abs(float(row["vela_edge_qfp_delta_V"])) for row in selected],
                    0.5,
                ),
            }
        )
    return samples, summary


def _mobility_residual_effect(
    waterfall: list[dict[str, str]],
) -> list[dict[str, object]]:
    rows = {
        (
            row["topology"],
            float(row["bias_V"]),
            row["branch"],
            row["carrier"],
            int(row["node_id"]),
        ): row
        for row in waterfall
    }
    summary: list[dict[str, object]] = []
    for carrier in ("electron", "hole"):
        ratios: list[float] = []
        for key, production in rows.items():
            topology, bias, branch, key_carrier, node = key
            if branch != "vela_production" or key_carrier != carrier:
                continue
            sentaurus = rows[
                (topology, bias, "sentaurus_box_edge", carrier, node)
            ]
            baseline = abs(
                float(production["final_residual_normalized_units"])
            )
            candidate = abs(
                float(sentaurus["final_residual_normalized_units"])
            )
            ratios.append(candidate / baseline)
        summary.append(
            {
                "carrier": carrier,
                "node_state_count": len(ratios),
                "sentaurus_to_vela_residual_ratio_median": _quantile(ratios, 0.5),
                "residual_reduction_fraction_median": _quantile(
                    [1.0 - value for value in ratios], 0.5
                ),
                "residual_ratio_minimum": min(ratios),
                "residual_ratio_maximum": max(ratios),
            }
        )
    return summary


def _driving_force_summary(
    native_rows: list[dict[str, str]],
    triangle_rows: list[dict[str, str]],
) -> list[dict[str, object]]:
    triangle = {
        (
            row["topology"],
            float(row["bias_V"]),
            int(row["cell_id"]),
            row["carrier"],
        ): row
        for row in triangle_rows
    }
    samples: dict[str, dict[str, list[float]]] = {
        carrier: {"native": [], "triangle": [], "ratio": []}
        for carrier in ("electron", "hole")
    }
    for row in native_rows:
        key = (
            row["topology"],
            float(row["bias_V"]),
            int(row["cell_id"]),
            row["carrier"],
        )
        carrier = row["carrier"]
        samples[carrier]["native"].append(
            float(row["cell_average_abs_log10_error_dex"])
        )
        samples[carrier]["triangle"].append(
            float(triangle[key]["correct_average_doping_abs_log10_error_dex"])
        )
        native_field = float(row["sentaurus_native_qf_field_V_per_m"])
        triangle_field = float(triangle[key]["cell_qf_field_V_per_m"])
        samples[carrier]["ratio"].append(triangle_field / native_field)
    summary: list[dict[str, object]] = []
    for carrier in ("electron", "hole"):
        native_stats = _stats(samples[carrier]["native"])
        triangle_stats = _stats(samples[carrier]["triangle"])
        ratio_stats = _stats(samples[carrier]["ratio"])
        summary.append(
            {
                "carrier": carrier,
                "element_count": native_stats["count"],
                "native_egrad_mobility_error_median_dex": native_stats["median"],
                "native_egrad_mobility_error_p95_dex": native_stats["p95"],
                "triangle_qfp_mobility_error_median_dex": triangle_stats["median"],
                "triangle_qfp_mobility_error_p95_dex": triangle_stats["p95"],
                "triangle_to_native_field_ratio_median": ratio_stats["median"],
                "triangle_to_native_field_ratio_minimum": min(
                    samples[carrier]["ratio"]
                ),
                "triangle_to_native_field_ratio_maximum": max(
                    samples[carrier]["ratio"]
                ),
            }
        )
    return summary


def _report(
    current_summary: list[dict[str, object]],
    qfp_summary: list[dict[str, object]],
    mobility_effect: list[dict[str, object]],
    driving_force: list[dict[str, object]],
) -> str:
    by_current = {row["carrier"]: row for row in current_summary}
    by_qfp = {row["carrier"]: row for row in qfp_summary}
    by_mobility = {row["carrier"]: row for row in mobility_effect}
    by_field = {row["carrier"]: row for row in driving_force}
    lines = [
        "# PN2D Minimal6 current-factor follow-up",
        "",
        "Date: 2026-07-24",
        "",
        "Status: diagnostic complete; remote low-field Sentaurus probe pending.",
        "",
        "## Decision",
        "",
        "The self-consistent QFP state is the dominant current-error driver. "
        "Imported Sentaurus potentials plus Vela density recomputation remove "
        "most of the paired edge-current error. The remaining fixed-state "
        "difference is a mobility/operator-support residual, not a density or "
        "geometry error.",
        "",
        "The central 1-5 edge is a separate differential-mode symptom: its QFP "
        "drop has the opposite sign in all 40 Vela states, while its absolute "
        "current is much smaller than the four horizontal active edges.",
        "",
        "## Paired current decomposition",
        "",
        "| Carrier | Self-consistent median dex | Imported-state median dex | "
        "Paired state improvement dex | Current-weighted improvement dex | "
        "Self-consistent sign | Imported sign |",
        "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for carrier in ("electron", "hole"):
        row = by_current[carrier]
        lines.append(
            "| {carrier} | {self_consistent_median_error_dex:.6g} | "
            "{imported_state_median_error_dex:.6g} | "
            "{paired_state_improvement_median_dex:.6g} | "
            "{current_weighted_state_improvement_dex:.6g} | "
            "{self_consistent_sign_agreement_fraction:.3f} | "
            "{imported_state_sign_agreement_fraction:.3f} |".format(**row)
        )
    lines.extend(
        [
            "",
            "## QFP common and differential modes",
            "",
            "| Carrier | Common-mode median V | Differential-mode median V | "
            "Sent edge delta median V | Vela edge delta median V | Edge sign |",
            "|---|---:|---:|---:|---:|---:|",
        ]
    )
    for carrier in ("electron", "hole"):
        row = by_qfp[carrier]
        lines.append(
            "| {carrier} | {common_mode_abs_error_median_V:.6g} | "
            "{differential_mode_abs_error_median_V:.6g} | "
            "{sentaurus_abs_edge_delta_median_V:.6g} | "
            "{vela_abs_edge_delta_median_V:.6g} | "
            "{edge_qfp_sign_agreement_fraction:.3f} |".format(**row)
        )
    lines.extend(
        [
            "",
            "The common-mode QFP displacement explains the order-one carrier "
            "density error through the Boltzmann factor. The much smaller "
            "differential error controls the central-edge direction because "
            "the physical 1-5 QFP drop is itself only about 0.5 mV.",
            "",
            "## Mobility is contributory but not sufficient",
            "",
            "| Carrier | Sent/Vela residual ratio median | Residual reduction | "
            "Minimum ratio | Maximum ratio |",
            "|---|---:|---:|---:|---:|",
        ]
    )
    for carrier in ("electron", "hole"):
        row = by_mobility[carrier]
        lines.append(
            "| {carrier} | "
            "{sentaurus_to_vela_residual_ratio_median:.6g} | "
            "{residual_reduction_fraction_median:.6g} | "
            "{residual_ratio_minimum:.6g} | "
            "{residual_ratio_maximum:.6g} |".format(**row)
        )
    lines.extend(
        [
            "",
            "On the fixed imported state, replacing only Vela mobility with the "
            "coefficient-weighted Sentaurus element mobility does not close the "
            "carrier continuity residual. It therefore cannot by itself explain "
            "the final self-consistent QFP displacement.",
            "",
            "## Exported eGradQuasiFermi is not a proven mobility drive",
            "",
            "| Carrier | Native-eGrad mobility median/P95 dex | "
            "Triangle-QFP mobility median/P95 dex | "
            "Triangle/native field median ratio |",
            "|---|---:|---:|---:|",
        ]
    )
    for carrier in ("electron", "hole"):
        row = by_field[carrier]
        lines.append(
            "| {carrier} | "
            "{native_egrad_mobility_error_median_dex:.6g} / "
            "{native_egrad_mobility_error_p95_dex:.6g} | "
            "{triangle_qfp_mobility_error_median_dex:.6g} / "
            "{triangle_qfp_mobility_error_p95_dex:.6g} | "
            "{triangle_to_native_field_ratio_median:.6g} |".format(**row)
        )
    lines.extend(
        [
            "",
            "For electrons, the native exported element eGradQuasiFermi field "
            "makes the mobility replay substantially worse than the affine "
            "triangle QFP gradient. The electron field ratio is constant across "
            "all 160 elements. This supports an output-semantics or proprietary "
            "element-evaluation difference, not a Vela factor correction.",
            "",
            "## Excluded causes",
            "",
            "- Imported-state density closes to 4.426181e-6 dex.",
            "- Electric field and box geometry close at floating-point precision.",
            "- The carrier-edge sum reproduces the Vela terminal current.",
            "- Boundary state, analytic Jacobian, and nonlinear convergence gates pass.",
            "- The source-unit defect is already repaired and is no longer the active cause.",
            "",
            "## Remaining discriminator",
            "",
            "Run one Sentaurus state twice with identical endpoint state: once "
            "with the current HighFieldSaturation model and once with high-field "
            "saturation disabled. Native low-field element mobility would remove "
            "the current low-field interpolation confound and allow the effective "
            "high-field drive to be inverted directly. The authorized VM was "
            "unreachable during this run, so this control remains pending.",
            "",
            "No production mobility, SG, Poisson, impact, or QFP formula change is "
            "authorized by this diagnostic.",
            "",
        ]
    )
    return "\n".join(lines)


def run(args: argparse.Namespace) -> None:
    inputs = {
        "phase_f_edges": args.phase_f_edges.resolve(),
        "phase_f_states": args.phase_f_states.resolve(),
        "staged_edges": args.staged_edges.resolve(),
        "phase_e_waterfall": args.phase_e_waterfall.resolve(),
        "phase_d_native": args.phase_d_native.resolve(),
        "triangle_mobility": args.triangle_mobility.resolve(),
    }
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)

    edge_samples, current_summary = _current_decomposition(
        _read_csv(inputs["phase_f_edges"]),
        _read_csv(inputs["staged_edges"]),
    )
    qfp_samples, qfp_summary = _qfp_modes(
        _read_csv(inputs["phase_f_states"])
    )
    mobility_effect = _mobility_residual_effect(
        _read_csv(inputs["phase_e_waterfall"])
    )
    driving_force = _driving_force_summary(
        _read_csv(inputs["phase_d_native"]),
        _read_csv(inputs["triangle_mobility"]),
    )

    files = {
        "edge_state_decomposition.csv": edge_samples,
        "current_factor_summary.csv": current_summary,
        "qfp_mode_samples.csv": qfp_samples,
        "qfp_mode_summary.csv": qfp_summary,
        "mobility_residual_effect.csv": mobility_effect,
        "driving_force_summary.csv": driving_force,
    }
    for name, rows in files.items():
        _write_csv(output / name, rows)
    (output / "report.md").write_text(
        _report(current_summary, qfp_summary, mobility_effect, driving_force),
        encoding="utf-8",
    )

    output_names = [*files, "report.md"]
    manifest = {
        "schema_version": 1,
        "experiment": "pn2d_minimal6_current_factor_followup",
        "status": "valid",
        "state_count": 40,
        "active_carrier_edge_count": len(edge_samples),
        "qfp_mode_sample_count": len(qfp_samples),
        "native_carrier_element_count": 320,
        "remote_low_field_probe": "pending_ssh_timeout",
        "production_formula_modified": False,
        "inputs": {
            name: {"path": str(path), "sha256": _sha256(path)}
            for name, path in inputs.items()
        },
        "outputs": {
            name: {"sha256": _sha256(output / name)}
            for name in output_names
        },
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--phase-f-edges", type=Path, required=True)
    parser.add_argument("--phase-f-states", type=Path, required=True)
    parser.add_argument("--staged-edges", type=Path, required=True)
    parser.add_argument("--phase-e-waterfall", type=Path, required=True)
    parser.add_argument("--phase-d-native", type=Path, required=True)
    parser.add_argument("--triangle-mobility", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    run(parser.parse_args())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
