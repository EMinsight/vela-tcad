#!/usr/bin/env python3
"""Verify carrier-resolved QFP residual and Jacobian behavior on M2 states."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
import sys
from pathlib import Path
from typing import Any

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.diagnose_pn2d_bv_predictor_first_step_audit import make_probe_config
from scripts.run_pn2d_bv_m2_sentaurus_frozen_sg_laux import (
    bias_key,
    bias_tag,
    run_audit,
    summarize_audit,
)
from scripts.run_pn2d_bv_m2_single_family_state_substitution import (
    contact_nodes,
    exact_bias_record,
    normalize_sentaurus_state,
    normalize_vela_state,
    read_json,
    read_rows,
    write_external_fields,
    write_json,
    write_rows,
)


DEFAULT_BIASES = (-18.0, -19.5, -19.7, -20.0)
VARIANTS = (
    "vela_baseline",
    "sent_phin_only",
    "sent_phip_only",
    "sent_qfp_only",
)
JACOBIAN_STATES = ("vela_baseline", "sent_qfp_only")
TERMS = ("flux", "recombination", "impact", "gauge", "boundary")
CARRIERS = ("electron", "hole")
FD_THRESHOLD = 5.0e-5
FD_SENSITIVITY_BIASES = (-19.5, -20.0)
FD_SENSITIVITY_STEPS = (1.0e-5, 3.0e-6, 1.0e-6, 3.0e-7, 1.0e-7, 3.0e-8, 1.0e-8)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def mixed_qfp_state(
    vela: list[dict[str, object]],
    sentaurus: list[dict[str, object]],
    variant: str,
) -> list[dict[str, object]]:
    replacements = {
        "vela_baseline": set(),
        "sent_phin_only": {"phin_V"},
        "sent_phip_only": {"phip_V"},
        "sent_qfp_only": {"phin_V", "phip_V"},
    }[variant]
    if len(vela) != len(sentaurus):
        raise RuntimeError("Vela/Sentaurus state node counts do not match")
    output = []
    for left, right in zip(vela, sentaurus):
        if left["node_id"] != right["node_id"]:
            raise RuntimeError("Vela/Sentaurus state node IDs do not match")
        output.append(
            {
                "node_id": left["node_id"],
                **{
                    column: right[column] if column in replacements else left[column]
                    for column in ("psi_V", "phin_V", "phip_V", "n_m3", "p_m3")
                },
            }
        )
    return output


def write_restart(path: Path, state: list[dict[str, object]]) -> None:
    write_rows(
        path,
        [
            {
                "node_id": row["node_id"],
                "psi": row["psi_V"],
                "phin": row["phin_V"],
                "phip": row["phip_V"],
                "electrons_m3": row["n_m3"],
                "holes_m3": row["p_m3"],
            }
            for row in state
        ],
    )


def run_runner(runner: Path, config: Path) -> dict[str, Any]:
    completed = subprocess.run(
        [str(runner.resolve()), "--config", str(config.resolve())],
        cwd=config.parent,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    (config.parent / f"{config.stem}.stdout.log").write_text(
        completed.stdout, encoding="utf-8"
    )
    (config.parent / f"{config.stem}.stderr.log").write_text(
        completed.stderr, encoding="utf-8"
    )
    if completed.returncode:
        raise RuntimeError(
            f"runner failed for {config}: "
            f"{completed.stderr.strip() or completed.stdout.strip()}"
        )
    return json.loads(completed.stdout)


def source_error(source: float, golden: float) -> float:
    return abs(math.log10(abs(source) / abs(golden)))


def source_rows_for_bias(
    bias: float,
    sources: dict[str, float],
    golden: float,
) -> list[dict[str, object]]:
    baseline_error = source_error(sources["vela_baseline"], golden)
    qfp_error = source_error(sources["sent_qfp_only"], golden)
    qfp_removal = baseline_error - qfp_error
    output = []
    for variant in VARIANTS:
        error = source_error(sources[variant], golden)
        removal = baseline_error - error
        output.append(
            {
                "bias_V": bias,
                "variant": variant,
                "source_A_per_um": sources[variant],
                "sentaurus_source_A_per_um": golden,
                "source_to_sentaurus_ratio": sources[variant] / golden,
                "abs_log10_error_dex": error,
                "error_removal_from_vela_dex": removal,
                "fraction_of_qfp_error_removal": (
                    removal / qfp_removal if qfp_removal else None
                ),
            }
        )
    return output


def l2(values: list[float]) -> float:
    return math.sqrt(sum(value * value for value in values))


def term_rows_for_bias(
    *,
    bias: float,
    by_variant: dict[str, list[dict[str, str]]],
    interior: set[int],
    contacts: set[int],
) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    baseline = {
        int(row["node_id"]): row for row in by_variant["vela_baseline"]
    }
    output: list[dict[str, object]] = []
    closure: list[dict[str, object]] = []
    scopes = {"interior": interior, "contacts": contacts, "all_nodes": set(baseline)}
    for variant, rows in by_variant.items():
        selected = {int(row["node_id"]): row for row in rows}
        if set(selected) != set(baseline):
            raise RuntimeError(f"{bias:g} V {variant}: carrier-term node set changed")
        for carrier in CARRIERS:
            maximum_closure = 0.0
            for node, row in selected.items():
                term_sum = sum(float(row[f"{carrier}_{term}"]) for term in TERMS)
                residual = float(row[f"{carrier}_residual"])
                maximum_closure = max(maximum_closure, abs(term_sum - residual))
            closure.append(
                {
                    "bias_V": bias,
                    "variant": variant,
                    "carrier": carrier,
                    "maximum_absolute_term_sum_closure": maximum_closure,
                }
            )
            for scope, nodes in scopes.items():
                residual_delta = [
                    float(selected[node][f"{carrier}_residual"])
                    - float(baseline[node][f"{carrier}_residual"])
                    for node in sorted(nodes)
                ]
                term_deltas: dict[str, list[float]] = {
                    term: [
                        float(selected[node][f"{carrier}_{term}"])
                        - float(baseline[node][f"{carrier}_{term}"])
                        for node in sorted(nodes)
                    ]
                    for term in TERMS
                }
                denominator = sum(l2(values) for values in term_deltas.values())
                for term, values in term_deltas.items():
                    magnitude = l2(values)
                    hotspot_node = (
                        sorted(nodes)[max(range(len(values)), key=lambda i: abs(values[i]))]
                        if values else -1
                    )
                    output.append(
                        {
                            "bias_V": bias,
                            "variant": variant,
                            "carrier": carrier,
                            "scope": scope,
                            "term": term,
                            "term_delta_l2": magnitude,
                            "term_delta_share": magnitude / denominator if denominator else 0.0,
                            "residual_delta_l2": l2(residual_delta),
                            "maximum_absolute_term_delta": (
                                max((abs(value) for value in values), default=0.0)
                            ),
                            "hotspot_node_id": hotspot_node,
                        }
                    )
    return output, closure


def update_projection(
    update: list[float], target: list[float]
) -> tuple[float | None, float, float]:
    norm_sq = sum(value * value for value in target)
    update_norm = l2(update)
    target_norm = math.sqrt(norm_sq)
    projection = (
        sum(left * right for left, right in zip(update, target)) / norm_sq
        if norm_sq else None
    )
    return projection, update_norm, target_norm


def first_update_rows(
    *,
    bias: float,
    variant: str,
    state: list[dict[str, object]],
    vela: list[dict[str, object]],
    sentaurus: list[dict[str, object]],
    probe_rows: list[dict[str, str]],
    status: dict[str, Any],
    interior: list[int],
) -> list[dict[str, object]]:
    by_node = {int(row["node_id"]): row for row in probe_rows}
    output = []
    for carrier, field, delta_column, trial_column, residual_name in (
        ("electron", "phin_V", "delta_phin_V", "trial_phin", "phin"),
        ("hole", "phip_V", "delta_phip_V", "trial_phip", "phip"),
    ):
        update = [float(by_node[node][delta_column]) for node in interior]
        target = [
            float(sentaurus[node][field]) - float(vela[node][field])
            for node in interior
        ]
        projection, update_norm, target_norm = update_projection(update, target)
        trial_distance = l2(
            [
                float(by_node[node][trial_column]) - float(sentaurus[node][field])
                for node in interior
            ]
        )
        initial_distance = l2(
            [
                float(state[node][field]) - float(sentaurus[node][field])
                for node in interior
            ]
        )
        output.append(
            {
                "bias_V": bias,
                "variant": variant,
                "carrier": carrier,
                "selected_field_is_sentaurus": int(
                    (carrier == "electron" and variant in {"sent_phin_only", "sent_qfp_only"})
                    or (carrier == "hole" and variant in {"sent_phip_only", "sent_qfp_only"})
                ),
                "initial_residual": float(status["block_residuals"][residual_name]),
                "trial_residual": float(status["trial_block_residuals"][residual_name]),
                "update_l2_V": update_norm,
                "vela_to_sentaurus_target_l2_V": target_norm,
                "update_projection_on_vela_to_sentaurus_target": projection,
                "initial_distance_to_sentaurus_l2_V": initial_distance,
                "trial_distance_to_sentaurus_l2_V": trial_distance,
                "trial_to_vela_target_distance_ratio": (
                    trial_distance / target_norm if target_norm else None
                ),
            }
        )
    return output


def jacobian_rows(
    bias: float, state_variant: str, path: Path
) -> list[dict[str, object]]:
    return [
        {
            "bias_V": bias,
            "state_variant": state_variant,
            **row,
        }
        for row in read_rows(path)
    ]


def classify(
    source_rows: list[dict[str, object]],
    term_rows: list[dict[str, object]],
    update_rows: list[dict[str, object]],
    jacobian: list[dict[str, object]],
    sensitivity: list[dict[str, object]],
) -> dict[str, object]:
    source_variants = ("sent_phin_only", "sent_phip_only")
    source_wins = {variant: 0 for variant in source_variants}
    source_by_bias: dict[str, dict[str, float]] = {}
    for bias in DEFAULT_BIASES:
        selected = {
            str(row["variant"]): float(row["fraction_of_qfp_error_removal"])
            for row in source_rows
            if float(row["bias_V"]) == bias and row["variant"] in source_variants
        }
        winner = max(source_variants, key=lambda name: abs(selected[name]))
        source_wins[winner] += 1
        source_by_bias[bias_key(bias)] = selected
    at_m20 = source_by_bias["-20"]
    source_driver = max(source_variants, key=lambda name: abs(at_m20[name]))
    source_fraction = at_m20[source_driver]
    source_outcome = (
        source_driver.removeprefix("sent_").removesuffix("_only") + "_dominant"
        if abs(source_fraction) >= 0.60 and source_wins[source_driver] >= 3
        else "phin_phip_interaction"
    )

    term_wins: dict[str, dict[str, int]] = {
        carrier: {term: 0 for term in TERMS} for carrier in CARRIERS
    }
    term_at_m20: dict[str, dict[str, float]] = {}
    for carrier in CARRIERS:
        for bias in DEFAULT_BIASES:
            selected = {
                str(row["term"]): float(row["term_delta_share"])
                for row in term_rows
                if row["variant"] == "sent_qfp_only"
                and row["scope"] == "interior"
                and row["carrier"] == carrier
                and float(row["bias_V"]) == bias
            }
            winner = max(TERMS, key=lambda term: selected[term])
            term_wins[carrier][winner] += 1
            if bias == -20.0:
                term_at_m20[carrier] = selected
    carrier_term_outcomes: dict[str, str] = {}
    for carrier in CARRIERS:
        winner = max(TERMS, key=lambda term: term_at_m20[carrier][term])
        share = term_at_m20[carrier][winner]
        carrier_term_outcomes[carrier] = (
            f"{winner}_dominant"
            if share >= 0.60 and term_wins[carrier][winner] >= 3
            else "multi_term"
        )

    relevant_fields = {
        "poisson": ("rel_phin_column_diff", "rel_phip_column_diff"),
        "transport": (
            "rel_electron_phin_diff", "rel_electron_phip_diff",
            "rel_hole_phin_diff", "rel_hole_phip_diff",
        ),
        "srh_auger": (
            "rel_electron_phin_diff", "rel_electron_phip_diff",
            "rel_hole_phin_diff", "rel_hole_phip_diff",
        ),
        "sg_avalanche": (
            "rel_electron_phin_diff", "rel_electron_phip_diff",
            "rel_hole_phin_diff", "rel_hole_phip_diff",
        ),
        "dirichlet_or_gauge": (
            "rel_electron_phin_diff", "rel_electron_phip_diff",
            "rel_hole_phin_diff", "rel_hole_phip_diff",
        ),
    }
    fd_checks = []
    for row in jacobian:
        for field in relevant_fields[str(row["block"])]:
            fd_checks.append(
                {
                    "bias_V": float(row["bias_V"]),
                    "state_variant": row["state_variant"],
                    "block": row["block"],
                    "subblock": field.removeprefix("rel_").removesuffix("_diff"),
                    "relative_difference": float(row[field]),
                }
            )
    worst_fd = max(fd_checks, key=lambda row: row["relative_difference"])
    fd_outcome = (
        "analytic_fd_consistent"
        if worst_fd["relative_difference"] <= FD_THRESHOLD
        else "analytic_fd_inconsistent"
    )
    non_srh_checks = [row for row in fd_checks if row["block"] != "srh_auger"]
    worst_non_srh = max(non_srh_checks, key=lambda row: row["relative_difference"])
    sensitivity_fields = (
        "electron_phin", "electron_phip", "hole_phin", "hole_phip"
    )
    srh_absolute_checks = []
    for row in sensitivity:
        for field in sensitivity_fields:
            srh_absolute_checks.append(
                {
                    "bias_V": float(row["bias_V"]),
                    "state_variant": row["state_variant"],
                    "finite_difference_step": float(row["finite_difference_step"]),
                    "subblock": field,
                    "analytic_norm": float(row[f"analytic_{field}_norm"]),
                    "finite_difference_norm": float(row[f"fd_{field}_norm"]),
                    "absolute_difference": float(row[f"diff_{field}_norm"]),
                    "relative_difference": float(row[f"rel_{field}_diff"]),
                }
            )
    worst_srh_absolute = max(
        srh_absolute_checks, key=lambda row: row["absolute_difference"]
    )
    best_srh_absolute = min(
        srh_absolute_checks, key=lambda row: row["absolute_difference"]
    )
    srh_at_fd_floor = (
        max(row["analytic_norm"] for row in srh_absolute_checks) <= 1.0e-13
        and worst_srh_absolute["absolute_difference"] <= 1.0e-13
        and worst_non_srh["relative_difference"] <= FD_THRESHOLD
    )
    fd_interpretation = (
        "formal_relative_gate_fails_only_at_srh_absolute_fd_floor"
        if fd_outcome == "analytic_fd_inconsistent" and srh_at_fd_floor
        else fd_outcome
    )

    selected_updates = [
        row for row in update_rows
        if row["variant"] == "sent_qfp_only" and float(row["bias_V"]) == -20.0
    ]
    rollback = {
        str(row["carrier"]): row["update_projection_on_vela_to_sentaurus_target"]
        for row in selected_updates
    }
    rollback_outcome = (
        "both_qfp_updates_roll_back_from_sentaurus"
        if all(value is not None and float(value) < 0.0 for value in rollback.values())
        else "qfp_update_direction_mixed"
    )
    typed = "__".join(
        [
            source_outcome,
            f"electron_{carrier_term_outcomes['electron']}",
            f"hole_{carrier_term_outcomes['hole']}",
            fd_outcome,
            rollback_outcome,
        ]
    )
    return {
        "typed_outcome": typed,
        "source_outcome": source_outcome,
        "source_fraction_at_minus20": source_fraction,
        "source_win_count": source_wins,
        "source_fraction_by_bias": source_by_bias,
        "carrier_term_outcomes": carrier_term_outcomes,
        "term_win_count": term_wins,
        "term_share_at_minus20": term_at_m20,
        "jacobian_fd_outcome": fd_outcome,
        "jacobian_fd_interpretation": fd_interpretation,
        "jacobian_fd_threshold": FD_THRESHOLD,
        "worst_jacobian_fd_check": worst_fd,
        "worst_non_srh_jacobian_fd_check": worst_non_srh,
        "srh_fd_step_sensitivity": {
            "biases_V": list(FD_SENSITIVITY_BIASES),
            "steps": list(FD_SENSITIVITY_STEPS),
            "best_absolute_check": best_srh_absolute,
            "worst_absolute_check": worst_srh_absolute,
            "absolute_floor_threshold": 1.0e-13,
            "classified_as_absolute_fd_floor": srh_at_fd_floor,
        },
        "minus20_qfp_update_projection": rollback,
        "qfp_update_outcome": rollback_outcome,
        "decision_rules": {
            "source_carrier_dominant": (
                "absolute recovery fraction >= 0.60 at -20 V and larger absolute "
                "contribution at least three of four biases"
            ),
            "residual_term_dominant": (
                "interior term-delta L2 share >= 0.60 at -20 V and wins at "
                "least three of four biases"
            ),
            "jacobian_fd_consistent": "all selected relative differences <= 5e-5",
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--doping", type=Path, required=True)
    parser.add_argument("--base-config", type=Path, required=True)
    parser.add_argument("--sentaurus-frozen-root", type=Path, required=True)
    parser.add_argument("--vela-manifest", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--repeats", type=int, default=2)
    args = parser.parse_args()
    if args.repeats < 2:
        raise RuntimeError("at least two repeats are required")
    args.output_root.mkdir(parents=True, exist_ok=True)

    mesh = read_json(args.mesh)
    contacts = contact_nodes(mesh)
    interior = [
        int(node["id"]) for node in mesh["nodes"]
        if int(node["id"]) not in contacts
    ]
    manifest = read_json(args.vela_manifest)
    manifest_root = args.vela_manifest.parent
    prior_source = read_rows(args.sentaurus_frozen_root / "source_comparison.csv")
    golden = {
        float(row["bias_V"]): float(row["sentaurus_source_A_per_um"])
        for row in prior_source if row["carrier"] == "total"
    }

    base_config = read_json(args.base_config)
    impact = base_config["solver"]["impact_ionization"]
    required = {
        "current_approximation": "element_edge_sg_gss_laux",
        "source_mapping_mode": "element_vertex_box_measure",
        "driving_force": "quasi_fermi_gradient",
        "generation": "current_density",
    }
    for key, value in required.items():
        if impact.get(key) != value:
            raise RuntimeError(f"impact_ionization.{key} must be {value}")
    fixed_config = json.loads(json.dumps(base_config))
    fixed_config["solver"]["impact_ionization"]["coupling_mode"] = "postprocess_only"
    fixed_config_path = args.output_root / "postprocess_only_sg_laux.json"
    write_json(fixed_config_path, fixed_config)
    probe_config = json.loads(json.dumps(fixed_config))
    probe_config["solver"]["impact_ionization"]["coupling_mode"] = "self_consistent"
    probe_config_path = args.output_root / "self_consistent_sg_laux_probe.json"
    write_json(probe_config_path, probe_config)

    states: dict[float, dict[str, list[dict[str, object]]]] = {}
    for bias in DEFAULT_BIASES:
        record = exact_bias_record(manifest, bias)
        vela_path = manifest_root / record["snapshot_tdr"]["path"]
        sent_path = args.sentaurus_frozen_root / f"sentaurus_state_{bias_tag(bias)}.csv"
        vela = normalize_vela_state(vela_path)
        sentaurus = normalize_sentaurus_state(sent_path)
        states[bias] = {
            variant: mixed_qfp_state(vela, sentaurus, variant)
            for variant in VARIANTS
        }
        input_root = args.output_root / "inputs" / bias_tag(bias)
        for variant, state in states[bias].items():
            write_rows(input_root / f"{variant}.csv", state)
            write_external_fields(input_root / f"{variant}_fields", state)
            write_restart(input_root / f"{variant}_restart.csv", state)

    repeat_source: dict[str, list[dict[str, object]]] = {}
    repeat_terms: dict[str, list[dict[str, object]]] = {}
    repeat_closure: dict[str, list[dict[str, object]]] = {}
    repeat_updates: dict[str, list[dict[str, object]]] = {}
    repeat_jacobian: dict[str, list[dict[str, object]]] = {}
    repeat_sensitivity: dict[str, list[dict[str, object]]] = {}
    hashes: dict[str, dict[str, str]] = {}
    for repeat in range(args.repeats):
        label = chr(ord("a") + repeat)
        source_output: list[dict[str, object]] = []
        term_output: list[dict[str, object]] = []
        closure_output: list[dict[str, object]] = []
        update_output: list[dict[str, object]] = []
        jacobian_output: list[dict[str, object]] = []
        sensitivity_output: list[dict[str, object]] = []
        run_hashes: dict[str, str] = {}
        for bias in DEFAULT_BIASES:
            tag = bias_tag(bias)
            source_values: dict[str, float] = {}
            term_by_variant: dict[str, list[dict[str, str]]] = {}
            for variant in VARIANTS:
                case = args.output_root / f"run-{label}" / tag / variant
                audit_paths = run_audit(
                    args.audit,
                    args.mesh,
                    args.doping,
                    args.output_root / "inputs" / tag / f"{variant}.csv",
                    fixed_config_path,
                    case / "fixed_source",
                )
                audit_summary = summarize_audit(audit_paths, states[bias][variant])
                source_values[variant] = float(audit_summary["source_A_per_um"]["total"])
                for name, path in audit_paths.items():
                    run_hashes[f"{tag}/{variant}/fixed/{name}"] = sha256(path)

                fields = args.output_root / "inputs" / tag / f"{variant}_fields"
                term_csv = case / "carrier_terms.csv"
                term_cfg_path = case / "carrier_terms.json"
                term_cfg = make_probe_config(
                    probe_config_path,
                    term_csv,
                    fields,
                    "newton_carrier_term_probe",
                    bias,
                    "Anode",
                    "Cathode",
                )
                write_json(term_cfg_path, term_cfg)
                write_json(case / "carrier_terms_status.json", run_runner(args.runner, term_cfg_path))
                term_by_variant[variant] = read_rows(term_csv)
                run_hashes[f"{tag}/{variant}/terms"] = sha256(term_csv)

                step_csv = case / "first_update.csv"
                step_cfg_path = case / "first_update.json"
                step_cfg = make_probe_config(
                    probe_config_path,
                    step_csv,
                    fields,
                    "newton_step_probe",
                    bias,
                    "Anode",
                    "Cathode",
                )
                write_json(step_cfg_path, step_cfg)
                step_status = run_runner(args.runner, step_cfg_path)
                write_json(case / "first_update_status.json", step_status)
                update_output.extend(
                    first_update_rows(
                        bias=bias,
                        variant=variant,
                        state=states[bias][variant],
                        vela=states[bias]["vela_baseline"],
                        sentaurus=mixed_qfp_state(
                            states[bias]["vela_baseline"],
                            normalize_sentaurus_state(
                                args.sentaurus_frozen_root
                                / f"sentaurus_state_{bias_tag(bias)}.csv"
                            ),
                            "sent_qfp_only",
                        ),
                        probe_rows=read_rows(step_csv),
                        status=step_status,
                        interior=interior,
                    )
                )
                run_hashes[f"{tag}/{variant}/step"] = sha256(step_csv)

                if variant in JACOBIAN_STATES:
                    jac_csv = case / "jacobian_blocks.csv"
                    jac_cfg_path = case / "jacobian_blocks.json"
                    jac_cfg = make_probe_config(
                        probe_config_path,
                        jac_csv,
                        fields,
                        "newton_jacobian_block_probe",
                        bias,
                        "Anode",
                        "Cathode",
                    )
                    jac_cfg["state_file"] = str(
                        (args.output_root / "inputs" / tag / f"{variant}_restart.csv").resolve()
                    )
                    jac_cfg["finite_difference_step"] = 1.0e-7
                    jac_cfg["finite_difference_mode"] = "double_symmetric"
                    jac_cfg["blocks"] = [
                        "poisson",
                        "transport",
                        "srh_auger",
                        "sg_avalanche",
                        "dirichlet_or_gauge",
                    ]
                    write_json(jac_cfg_path, jac_cfg)
                    write_json(case / "jacobian_blocks_status.json", run_runner(args.runner, jac_cfg_path))
                    jacobian_output.extend(jacobian_rows(bias, variant, jac_csv))
                    run_hashes[f"{tag}/{variant}/jacobian"] = sha256(jac_csv)

            source_output.extend(source_rows_for_bias(bias, source_values, golden[bias]))
            terms, closure = term_rows_for_bias(
                bias=bias,
                by_variant=term_by_variant,
                interior=set(interior),
                contacts=contacts,
            )
            term_output.extend(terms)
            closure_output.extend(closure)
        for bias in FD_SENSITIVITY_BIASES:
            tag = bias_tag(bias)
            for variant in JACOBIAN_STATES:
                fields = args.output_root / "inputs" / tag / f"{variant}_fields"
                for step in FD_SENSITIVITY_STEPS:
                    step_tag = f"{step:.0e}".replace("-", "m")
                    case = (
                        args.output_root / f"run-{label}" / tag / variant
                        / "srh_fd_sensitivity" / step_tag
                    )
                    jac_csv = case / "jacobian_blocks.csv"
                    jac_cfg_path = case / "jacobian_blocks.json"
                    jac_cfg = make_probe_config(
                        probe_config_path,
                        jac_csv,
                        fields,
                        "newton_jacobian_block_probe",
                        bias,
                        "Anode",
                        "Cathode",
                    )
                    jac_cfg["state_file"] = str(
                        (args.output_root / "inputs" / tag / f"{variant}_restart.csv").resolve()
                    )
                    jac_cfg["finite_difference_step"] = step
                    jac_cfg["finite_difference_mode"] = "double_symmetric"
                    jac_cfg["blocks"] = ["srh_auger"]
                    write_json(jac_cfg_path, jac_cfg)
                    write_json(
                        case / "jacobian_blocks_status.json",
                        run_runner(args.runner, jac_cfg_path),
                    )
                    for row in jacobian_rows(bias, variant, jac_csv):
                        row["finite_difference_step"] = step
                        sensitivity_output.append(row)
                    run_hashes[
                        f"{tag}/{variant}/srh_fd_sensitivity/{step_tag}"
                    ] = sha256(jac_csv)
        repeat_source[label] = source_output
        repeat_terms[label] = term_output
        repeat_closure[label] = closure_output
        repeat_updates[label] = update_output
        repeat_jacobian[label] = jacobian_output
        repeat_sensitivity[label] = sensitivity_output
        hashes[label] = run_hashes

    write_rows(args.output_root / "source_carrier_substitution.csv", repeat_source["a"])
    write_rows(args.output_root / "carrier_term_decomposition.csv", repeat_terms["a"])
    write_rows(args.output_root / "carrier_term_closure.csv", repeat_closure["a"])
    write_rows(args.output_root / "first_qfp_updates.csv", repeat_updates["a"])
    write_rows(args.output_root / "jacobian_fd_blocks.csv", repeat_jacobian["a"])
    write_rows(
        args.output_root / "jacobian_fd_step_sensitivity.csv",
        repeat_sensitivity["a"],
    )
    determinism = []
    for key in sorted(hashes["a"]):
        values = {hashes[label][key] for label in hashes}
        determinism.append(
            {
                "artifact": key,
                "repeat_count": args.repeats,
                "unique_hash_count": len(values),
                "byte_identical": int(len(values) == 1),
                "sha256": hashes["a"][key],
            }
        )
    write_rows(args.output_root / "determinism.csv", determinism)

    verdict = classify(
        repeat_source["a"],
        repeat_terms["a"],
        repeat_updates["a"],
        repeat_jacobian["a"],
        repeat_sensitivity["a"],
    )
    result = {
        "schema": "vela.pn2d_bv_m2_qfp_carrier_jacobian_verification.v1",
        "status": "passed",
        "biases_V": list(DEFAULT_BIASES),
        "variants": list(VARIANTS),
        "jacobian_states": list(JACOBIAN_STATES),
        "finite_difference_step": 1.0e-7,
        "finite_difference_threshold": FD_THRESHOLD,
        "finite_difference_sensitivity_steps": list(FD_SENSITIVITY_STEPS),
        "physics_modified": False,
        "production_defaults_modified": False,
        "acceptance_thresholds_modified": False,
        "verdict": verdict,
        "determinism": {
            "repeat_count": args.repeats,
            "artifact_count": len(determinism),
            "all_byte_identical": all(row["byte_identical"] == 1 for row in determinism),
        },
        "maximum_term_closure": max(
            float(row["maximum_absolute_term_sum_closure"])
            for row in repeat_closure["a"]
        ),
        "input_hashes": {
            "audit": sha256(args.audit),
            "runner": sha256(args.runner),
            "mesh": sha256(args.mesh),
            "doping": sha256(args.doping),
            "base_config": sha256(args.base_config),
            "vela_manifest": sha256(args.vela_manifest),
            "prior_frozen_result": sha256(args.sentaurus_frozen_root / "result.json"),
        },
    }
    write_json(args.output_root / "result.json", result)
    print(json.dumps(verdict, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
