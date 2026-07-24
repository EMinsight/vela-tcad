#!/usr/bin/env python3
"""Carry high-field mobility branches into residual and first-step audits."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path

from pn2d_minimal6_diagnostics.phase_e_continuity_residual import (
    INTERNAL_NODES,
    _deck,
    _node_divergence,
    _rows,
    _run_probe,
    _source_ratio,
)


REPLAY_BRANCHES = {
    "vela_production": None,
    "sentaurus_native_final": "sentaurus_native_final",
    "sentaurus_lowfield_element_electric_field": (
        "sentaurus_lowfield_element_electric_field"
    ),
    "sentaurus_lowfield_element_triangle_qfp": (
        "sentaurus_lowfield_element_triangle_qfp"
    ),
    "constant": "constant",
}
EXECUTABLE_BRANCHES = {
    "vela_global_qfp_config": {
        "model": "masetti_field",
        "high_field_driving_force": "quasi_fermi_gradient",
    },
    "vela_global_electric_field_config": {
        "model": "masetti_field",
        "high_field_driving_force": "electric_field",
    },
    "constant": "constant",
}


def write_csv(path: Path, values: list[dict[str, object]]) -> None:
    if not values:
        raise ValueError(f"refusing to write empty CSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(values[0]))
        writer.writeheader()
        writer.writerows(values)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def quantile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * fraction
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    return ordered[lower] + (
        ordered[upper] - ordered[lower]
    ) * (position - lower)


def branch_config(
    base: dict[str, object],
    bias: int,
    simulation_type: str,
    fields: Path,
    output: Path,
    branch: str,
    restart: Path | None = None,
) -> dict[str, object]:
    mobility = EXECUTABLE_BRANCHES[branch]
    cfg = _deck(
        base,
        bias,
        simulation_type,
        fields,
        output,
        "constant" if mobility == "constant" else "production",
        state_file=restart,
    )
    cfg["solver"]["mobility"] = mobility
    return cfg


def key(row: dict[str, str]) -> tuple[str, float, str, int]:
    return (
        row["topology"],
        float(row["bias_V"]),
        row["carrier"],
        int(row["edge_id"]),
    )


def run(args: argparse.Namespace) -> dict[str, object]:
    paths = {
        "runner": args.runner.resolve(),
        "inverse": args.inverse.resolve(),
        "phase_e": args.phase_e.resolve(),
        "box_current": args.box_current.resolve(),
        "physical_replay": args.physical_replay.resolve(),
    }
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    raw = output / "raw"
    box = _rows(paths["box_current"])
    box_by_key_branch = {
        (key(row), row["branch"]): row for row in box
    }

    residual_rows: list[dict[str, object]] = []
    state_norms: list[dict[str, object]] = []
    maximum_production_closure = 0.0
    topologies = ("mirror", "sketch")
    for topology in topologies:
        for bias in range(1, 21):
            tag = f"m{bias}V"
            phase_raw = paths["phase_e"] / "raw" / topology / tag
            production_edges = _rows(
                phase_raw / "vela_production_edges.csv"
            )
            constant_edges = _rows(phase_raw / "constant_edges.csv")
            terms = _rows(phase_raw / "vela_production_terms.csv")
            physical_edges = _rows(
                paths["physical_replay"]
                / "self_consistent_replay"
                / topology
                / tag
                / "edges.csv"
            )
            ratios: dict[str, dict[str, dict[int, float]]] = {}
            for branch, source_branch in REPLAY_BRANCHES.items():
                ratios[branch] = {}
                for carrier in ("electron", "hole"):
                    production_mu = {
                        int(row["edge_id"]): float(
                            row[f"{carrier}_mobility_m2_V_s"]
                        )
                        for row in production_edges
                    }
                    if branch == "vela_production":
                        candidate_mu = production_mu
                    elif branch == "constant":
                        candidate_mu = {
                            int(row["edge_id"]): float(
                                row[f"{carrier}_mobility_m2_V_s"]
                            )
                            for row in constant_edges
                        }
                    else:
                        candidate_mu = {}
                        for edge_id, base_mu in production_mu.items():
                            row = box_by_key_branch.get(
                                (
                                    (
                                        topology,
                                        -float(bias),
                                        carrier,
                                        edge_id,
                                    ),
                                    str(source_branch),
                                )
                            )
                            raw_mu = (
                                ""
                                if row is None
                                else row[
                                    "candidate_mobility_m2_per_Vs"
                                ]
                            )
                            candidate_mu[edge_id] = (
                                base_mu if raw_mu == "" else float(raw_mu)
                            )
                    ratios[branch][carrier] = {
                        edge_id: (
                            candidate_mu[edge_id] / value
                            if value != 0.0
                            else 1.0
                        )
                        for edge_id, value in production_mu.items()
                    }
            source_ratios = {
                branch: {
                    carrier: _source_ratio(
                        physical_edges, ratios[branch][carrier], carrier
                    )
                    for carrier in ("electron", "hole")
                }
                for branch in REPLAY_BRANCHES
            }
            production_norm: dict[str, float] = {}
            local_norms: dict[tuple[str, str], float] = {}
            for branch in REPLAY_BRANCHES:
                carrier_residuals: dict[str, list[float]] = defaultdict(list)
                for carrier in ("electron", "hole"):
                    divergence, absolute = _node_divergence(
                        production_edges,
                        carrier,
                        ratios[branch][carrier],
                    )
                    for node in INTERNAL_NODES:
                        term = terms[node]
                        impact_source = (
                            float(term["impact_electron_source"])
                            * source_ratios[branch]["electron"][node]
                            + float(term["impact_hole_source"])
                            * source_ratios[branch]["hole"][node]
                        )
                        srh = float(term[f"{carrier}_recombination"])
                        impact = -impact_source
                        gauge = float(term[f"{carrier}_gauge"])
                        boundary = float(term[f"{carrier}_boundary"])
                        final = (
                            divergence[node]
                            + srh
                            + impact
                            + gauge
                            + boundary
                        )
                        recorded = float(
                            term[f"{carrier}_residual"]
                        )
                        if branch == "vela_production":
                            closure = abs(final - recorded) / max(
                                1.0, abs(final), abs(recorded)
                            )
                            maximum_production_closure = max(
                                maximum_production_closure, closure
                            )
                        carrier_residuals[carrier].append(final)
                        residual_rows.append(
                            {
                                "topology": topology,
                                "bias_V": -bias,
                                "branch": branch,
                                "support": (
                                    "native_cpp_global_edge"
                                    if branch
                                    in ("vela_production", "constant")
                                    else "box_operator_reconstruction"
                                ),
                                "carrier": carrier,
                                "node_id": node,
                                "sg_divergence_normalized": divergence[node],
                                "sg_abs_incident_normalized": absolute[node],
                                "srh_normalized": srh,
                                "impact_normalized": impact,
                                "gauge_normalized": gauge,
                                "boundary_normalized": boundary,
                                "final_residual_normalized_units": final,
                                "recorded_production_residual": recorded,
                                "source_unit_policy": (
                                    "unchanged_production_snapshot"
                                ),
                            }
                        )
                for carrier, values in carrier_residuals.items():
                    norm = math.hypot(*values)
                    local_norms[(branch, carrier)] = norm
                    if branch == "vela_production":
                        production_norm[carrier] = norm
            for branch in REPLAY_BRANCHES:
                for carrier in ("electron", "hole"):
                    norm = local_norms[(branch, carrier)]
                    state_norms.append(
                        {
                            "topology": topology,
                            "bias_V": -bias,
                            "branch": branch,
                            "carrier": carrier,
                            "residual_l2": norm,
                            "production_residual_l2": production_norm[
                                carrier
                            ],
                            "ratio_to_production": norm
                            / max(production_norm[carrier], 1.0e-300),
                        }
                    )

    residual_summary: list[dict[str, object]] = []
    for branch in REPLAY_BRANCHES:
        for carrier in ("electron", "hole"):
            selected = [
                float(row["ratio_to_production"])
                for row in state_norms
                if row["branch"] == branch
                and row["carrier"] == carrier
            ]
            residual_summary.append(
                {
                    "branch": branch,
                    "carrier": carrier,
                    "state_count": len(selected),
                    "median_ratio_to_production": statistics.median(
                        selected
                    ),
                    "p95_ratio_to_production": quantile(selected, 0.95),
                    "maximum_ratio_to_production": max(selected),
                    "improved_state_fraction": statistics.mean(
                        float(value < 1.0) for value in selected
                    ),
                }
            )

    update_rows: list[dict[str, object]] = []
    jacobian_rows: list[dict[str, object]] = []
    for topology in topologies:
        for bias in (1, 10, 20):
            tag = f"m{bias}V"
            phase_raw = paths["phase_e"] / "raw" / topology / tag
            fields = phase_raw / "imported_fields"
            restart = phase_raw / "imported_restart.csv"
            base_path = (
                paths["inverse"]
                / "vela"
                / "source"
                / "decks"
                / topology
                / f"{tag}.json"
            )
            base = json.loads(base_path.read_text(encoding="utf-8"))
            work = raw / topology / tag
            for branch in EXECUTABLE_BRANCHES:
                carrier_output = work / f"{branch}_carrier.csv"
                carrier_cfg = branch_config(
                    base,
                    bias,
                    "newton_block_step_probe",
                    fields,
                    carrier_output,
                    branch,
                )
                carrier_cfg["block_modes"] = ["carrier_only"]
                _run_probe(
                    paths["runner"],
                    carrier_cfg,
                    work / f"{branch}_carrier.json",
                )
                coupled_output = work / f"{branch}_coupled.csv"
                coupled_cfg = branch_config(
                    base,
                    bias,
                    "newton_step_probe",
                    fields,
                    coupled_output,
                    branch,
                )
                _run_probe(
                    paths["runner"],
                    coupled_cfg,
                    work / f"{branch}_coupled.json",
                )
                for mode, source in (
                    ("carrier_only", carrier_output),
                    ("coupled", coupled_output),
                ):
                    for row in _rows(source):
                        if int(row["node_id"]) not in INTERNAL_NODES:
                            continue
                        for carrier, delta, before, after in (
                            (
                                "electron",
                                "delta_phin_V",
                                "phin_residual",
                                "trial_phin_residual",
                            ),
                            (
                                "hole",
                                "delta_phip_V",
                                "phip_residual",
                                "trial_phip_residual",
                            ),
                        ):
                            update_rows.append(
                                {
                                    "topology": topology,
                                    "bias_V": -bias,
                                    "branch": branch,
                                    "support": "native_cpp_global_edge",
                                    "mode": mode,
                                    "carrier": carrier,
                                    "node_id": int(row["node_id"]),
                                    "delta_qfp_V": float(row[delta]),
                                    "absolute_delta_qfp_V": abs(
                                        float(row[delta])
                                    ),
                                    "residual_before": float(row[before]),
                                    "residual_after": float(row[after]),
                                }
                            )
                jac_output = work / f"{branch}_jacobian.csv"
                jac_cfg = branch_config(
                    base,
                    bias,
                    "newton_jacobian_block_probe",
                    fields,
                    jac_output,
                    branch,
                    restart,
                )
                jac_cfg["finite_difference_step"] = 1.0e-7
                _run_probe(
                    paths["runner"],
                    jac_cfg,
                    work / f"{branch}_jacobian.json",
                )
                for row in _rows(jac_output):
                    jacobian_rows.append(
                        {
                            "topology": topology,
                            "bias_V": -bias,
                            "branch": branch,
                            **row,
                        }
                    )

    update_summary: list[dict[str, object]] = []
    for branch in EXECUTABLE_BRANCHES:
        for mode in ("carrier_only", "coupled"):
            for carrier in ("electron", "hole"):
                selected = [
                    float(row["absolute_delta_qfp_V"])
                    for row in update_rows
                    if row["branch"] == branch
                    and row["mode"] == mode
                    and row["carrier"] == carrier
                ]
                update_summary.append(
                    {
                        "branch": branch,
                        "mode": mode,
                        "carrier": carrier,
                        "sample_count": len(selected),
                        "median_abs_delta_qfp_V": statistics.median(
                            selected
                        ),
                        "p95_abs_delta_qfp_V": quantile(selected, 0.95),
                        "maximum_abs_delta_qfp_V": max(selected),
                    }
                )
    production_update = {
        (row["mode"], row["carrier"]): float(
            row["median_abs_delta_qfp_V"]
        )
        for row in update_summary
        if row["branch"] == "vela_global_qfp_config"
    }
    electric_update = {
        (row["mode"], row["carrier"]): float(
            row["median_abs_delta_qfp_V"]
        )
        for row in update_summary
        if row["branch"] == "vela_global_electric_field_config"
    }
    electric_first_update_improves = all(
        electric_update[key_] < production_update[key_]
        for key_ in production_update
    )
    maximum_jacobian = max(
        float(row["rel_diff"]) for row in jacobian_rows
    )
    outcome = (
        "mobility_candidate_causal"
        if electric_first_update_improves
        else "current_coefficient_improvement_without_qfp_causality"
    )

    outputs = {
        "residual_waterfall.csv": residual_rows,
        "state_residual_norms.csv": state_norms,
        "residual_summary.csv": residual_summary,
        "first_update.csv": update_rows,
        "first_update_summary.csv": update_summary,
        "jacobian_audit.csv": jacobian_rows,
    }
    for name, values in outputs.items():
        write_csv(output / name, values)
    report = [
        "# PN2D Minimal6 high-field residual and first-step follow-up",
        "",
        f"Typed outcome: `{outcome}`.",
        "",
        f"Production edge-to-node replay closure: "
        f"`{maximum_production_closure:.6e}`.",
        "",
        f"Maximum analytic/FD Jacobian block difference: "
        f"`{maximum_jacobian:.6e}`.",
        "",
        f"Global electric-field first update improves every mode/carrier "
        f"median: `{str(electric_first_update_improves).lower()}`.",
        "",
        "The element-electric, native-element, and triangle branches are "
        "fixed-state box-operator reconstructions. First-step and Jacobian "
        "rows use executable native global-edge configurations only.",
        "",
        "All branches retain the same production source-unit policy; no SRH "
        "or impact scaling change is combined with mobility.",
        "",
    ]
    (output / "report.md").write_text(
        "\n".join(report), encoding="utf-8", newline="\n"
    )
    output_hashes = {
        name: sha256(output / name)
        for name in (*outputs, "report.md")
    }
    manifest = {
        "schema_version": 1,
        "status": (
            "valid"
            if maximum_production_closure <= 1.0e-12
            and maximum_jacobian <= 1.0e-8
            else "failed"
        ),
        "experiment": "pn2d_minimal6_phase_e_highfield_followup",
        "typed_outcome": outcome,
        "state_count": 40,
        "controlled_state_count": 6,
        "replay_branches": list(REPLAY_BRANCHES),
        "executable_branches": list(EXECUTABLE_BRANCHES),
        "maximum_production_edge_to_node_closure": (
            maximum_production_closure
        ),
        "maximum_jacobian_relative_difference": maximum_jacobian,
        "electric_first_update_improves_all_mode_carrier_medians": (
            electric_first_update_improves
        ),
        "source_unit_policy": "unchanged_production_snapshot",
        "production_formula_modified": False,
        "inputs": {
            name: {"path": str(path), "sha256": sha256(path)}
            for name, path in paths.items()
            if path.is_file()
        },
        "input_roots": {
            name: str(path)
            for name, path in paths.items()
            if path.is_dir()
        },
        "outputs": output_hashes,
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    print(json.dumps(manifest, indent=2))
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--inverse", type=Path, required=True)
    parser.add_argument("--phase-e", type=Path, required=True)
    parser.add_argument("--box-current", type=Path, required=True)
    parser.add_argument("--physical-replay", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
