#!/usr/bin/env python3
"""Run M2 single-family state substitutions and first-Newton-step probes.

The fixed-state half holds four of ``psi``, QFP, and ``n/p`` at the converged
Vela avalanche-on state while replacing one family from Sentaurus.  The Newton
half distinguishes independent-variable replacements (psi/QFP) from frozen
residual-operator feedback replacements (density/QFP).  No continuation
state is advanced and no production default is changed.
"""

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


DEFAULT_BIASES = (-18.0, -19.5, -19.7, -20.0)
STATE_COLUMNS = ("psi_V", "phin_V", "phip_V", "n_m3", "p_m3")
FIXED_VARIANTS = (
    "vela_baseline",
    "sent_psi_only",
    "sent_qfp_only",
    "sent_density_only",
    "sent_all",
)
INDEPENDENT_NEWTON_VARIANTS = (
    "vela_baseline",
    "sent_psi_only",
    "sent_qfp_only",
    "sent_all",
)
FEEDBACK_VARIANTS = ("baseline", "density_only", "qfp_only", "density_qfp")


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise RuntimeError(f"refusing to write empty CSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def exact_bias_record(manifest: dict[str, Any], bias: float) -> dict[str, Any]:
    branch = next(
        item for item in manifest["branch_records"]
        if item["branch"] == "avalanche_on"
    )
    matches = [
        item for item in branch["bias_records"]
        if math.isclose(
            float(item["requested_bias_V"]), bias, rel_tol=0.0, abs_tol=1.0e-10
        )
    ]
    if len(matches) != 1:
        raise RuntimeError(f"expected one Vela avalanche-on state at {bias:g} V")
    return matches[0]


def normalize_vela_state(path: Path) -> list[dict[str, object]]:
    rows = read_rows(path)
    expected = list(range(len(rows)))
    if [int(row["node_id"]) for row in rows] != expected:
        raise RuntimeError(f"{path}: Vela state node IDs are not canonical")
    return [
        {
            "node_id": int(row["node_id"]),
            "psi_V": float(row["psi"]),
            "phin_V": float(row["phin"]),
            "phip_V": float(row["phip"]),
            "n_m3": float(row["electrons_m3"]),
            "p_m3": float(row["holes_m3"]),
        }
        for row in rows
    ]


def normalize_sentaurus_state(path: Path) -> list[dict[str, object]]:
    rows = read_rows(path)
    expected = list(range(len(rows)))
    if [int(row["node_id"]) for row in rows] != expected:
        raise RuntimeError(f"{path}: Sentaurus state node IDs are not canonical")
    return [
        {
            "node_id": int(row["node_id"]),
            **{column: float(row[column]) for column in STATE_COLUMNS},
        }
        for row in rows
    ]


def mixed_state(
    vela: list[dict[str, object]],
    sentaurus: list[dict[str, object]],
    variant: str,
) -> list[dict[str, object]]:
    if len(vela) != len(sentaurus):
        raise RuntimeError("Vela/Sentaurus state node counts do not match")
    replacements = {
        "vela_baseline": set(),
        "sent_psi_only": {"psi_V"},
        "sent_qfp_only": {"phin_V", "phip_V"},
        "sent_density_only": {"n_m3", "p_m3"},
        "sent_all": set(STATE_COLUMNS),
    }[variant]
    output = []
    for left, right in zip(vela, sentaurus):
        if left["node_id"] != right["node_id"]:
            raise RuntimeError("Vela/Sentaurus state node IDs do not match")
        output.append(
            {
                "node_id": left["node_id"],
                **{
                    column: right[column] if column in replacements else left[column]
                    for column in STATE_COLUMNS
                },
            }
        )
    return output


def write_external_fields(path: Path, state: list[dict[str, object]]) -> None:
    mapping = {
        "ElectrostaticPotential": "psi_V",
        "eQuasiFermiPotential": "phin_V",
        "hQuasiFermiPotential": "phip_V",
    }
    for name, column in mapping.items():
        write_rows(
            path / f"{name}_region0.csv",
            [
                {"node_id": row["node_id"], "component0": row[column]}
                for row in state
            ],
        )


def write_feedback_fields(path: Path, state: list[dict[str, object]]) -> None:
    mapping = {
        "eQuasiFermiPotential": "phin_V",
        "hQuasiFermiPotential": "phip_V",
        "eDensity_m3": "n_m3",
        "hDensity_m3": "p_m3",
    }
    for name, column in mapping.items():
        write_rows(
            path / f"{name}_region0.csv",
            [
                {"node_id": row["node_id"], "component0": row[column]}
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


def contact_nodes(mesh: dict[str, Any]) -> set[int]:
    return {
        int(node)
        for contact in mesh.get("contacts", [])
        for node in contact.get("node_ids", [])
    }


def vector_metrics(
    initial: list[float],
    trial: list[float],
    target: list[float],
    selected: list[int],
) -> dict[str, float | None]:
    difference = [target[i] - initial[i] for i in selected]
    update = [trial[i] - initial[i] for i in selected]
    remaining = [target[i] - trial[i] for i in selected]
    target_norm = math.sqrt(sum(value * value for value in difference))
    update_norm = math.sqrt(sum(value * value for value in update))
    remaining_norm = math.sqrt(sum(value * value for value in remaining))
    if target_norm == 0.0:
        return {
            "initial_distance_l2_V": 0.0,
            "trial_distance_l2_V": remaining_norm,
            "trial_to_initial_distance_ratio": None,
            "update_target_cosine": None,
            "target_projection_fraction": None,
        }
    dot = sum(a * b for a, b in zip(update, difference))
    return {
        "initial_distance_l2_V": target_norm,
        "trial_distance_l2_V": remaining_norm,
        "trial_to_initial_distance_ratio": remaining_norm / target_norm,
        "update_target_cosine": (
            dot / (update_norm * target_norm) if update_norm else 0.0
        ),
        "target_projection_fraction": dot / (target_norm * target_norm),
    }


def update_direction_metrics(
    initial: list[dict[str, object]],
    target: list[dict[str, object]],
    trial_by_node: dict[int, tuple[float, float, float]],
    interior: list[int],
) -> dict[str, float | None]:
    result: dict[str, float | None] = {}
    families = {
        "psi": ("psi_V",),
        "qfp": ("phin_V", "phip_V"),
        "combined": ("psi_V", "phin_V", "phip_V"),
    }
    trial_column = {"psi_V": 0, "phin_V": 1, "phip_V": 2}
    for family, columns in families.items():
        initial_values = []
        target_values = []
        trial_values = []
        for node in interior:
            for column in columns:
                initial_values.append(float(initial[node][column]))
                target_values.append(float(target[node][column]))
                trial_values.append(trial_by_node[node][trial_column[column]])
        metrics = vector_metrics(
            initial_values,
            trial_values,
            target_values,
            list(range(len(initial_values))),
        )
        result.update({f"{family}_{key}": value for key, value in metrics.items()})
    return result


def block_values(status: dict[str, Any], key: str) -> dict[str, float]:
    return {
        name: float(status[key][name])
        for name in ("psi", "phin", "phip", "combined")
    }


def summarize_independent_step(
    *,
    bias: float,
    variant: str,
    state: list[dict[str, object]],
    sentaurus: list[dict[str, object]],
    output_csv: Path,
    status: dict[str, Any],
    interior: list[int],
) -> dict[str, object]:
    output = read_rows(output_csv)
    trial = {
        int(row["node_id"]): (
            float(row["trial_psi"]),
            float(row["trial_phin"]),
            float(row["trial_phip"]),
        )
        for row in output
    }
    initial_residual = block_values(status, "block_residuals")
    trial_residual = block_values(status, "trial_block_residuals")
    return {
        "bias_V": bias,
        "variant": variant,
        "intervention_kind": "independent_state",
        "raw_step_norm": float(status["raw_step_norm"]),
        "step_norm": float(status["step_norm"]),
        **{f"initial_{key}_residual": value for key, value in initial_residual.items()},
        **{f"trial_{key}_residual": value for key, value in trial_residual.items()},
        "trial_combined_to_initial_ratio": (
            trial_residual["combined"] / initial_residual["combined"]
            if initial_residual["combined"] else None
        ),
        **update_direction_metrics(state, sentaurus, trial, interior),
    }


def summarize_feedback_steps(
    *,
    bias: float,
    baseline: list[dict[str, object]],
    sentaurus: list[dict[str, object]],
    output_csv: Path,
    status: dict[str, Any],
    interior: list[int],
) -> list[dict[str, object]]:
    rows = read_rows(output_csv)
    status_by_variant = {
        item["variant"]: item for item in status["variants"]
    }
    output = []
    for variant in FEEDBACK_VARIANTS:
        selected_rows = [row for row in rows if row["variant"] == variant]
        trial = {
            int(row["node_id"]): (
                float(baseline[int(row["node_id"])]["psi_V"])
                + float(row["delta_psi_V"]),
                float(row["trial_phin_V"]),
                float(row["trial_phip_V"]),
            )
            for row in selected_rows
        }
        item = status_by_variant[variant]
        initial_residual = {
            key: float(item["block_residuals"][key])
            for key in ("psi", "phin", "phip", "combined")
        }
        trial_residual = {
            key: float(item["production_trial_block_residuals"][key])
            for key in ("psi", "phin", "phip", "combined")
        }
        output.append(
            {
                "bias_V": bias,
                "variant": f"feedback_{variant}",
                "intervention_kind": "frozen_operator_feedback",
                "raw_step_norm": float(item["raw_step_norm"]),
                "step_norm": float(item["step_norm"]),
                **{
                    f"initial_{key}_residual": value
                    for key, value in initial_residual.items()
                },
                **{
                    f"trial_{key}_residual": value
                    for key, value in trial_residual.items()
                },
                "trial_combined_to_initial_ratio": (
                    trial_residual["combined"] / initial_residual["combined"]
                    if initial_residual["combined"] else None
                ),
                **update_direction_metrics(baseline, sentaurus, trial, interior),
            }
        )
    return output


def source_error(source: float, golden: float) -> float:
    return abs(math.log10(abs(source) / abs(golden)))


def source_rows_for_bias(
    bias: float,
    sources: dict[str, float],
    golden: float,
) -> list[dict[str, object]]:
    baseline_error = source_error(sources["vela_baseline"], golden)
    all_error = source_error(sources["sent_all"], golden)
    full_removal = baseline_error - all_error
    output = []
    for variant in FIXED_VARIANTS:
        error = source_error(sources[variant], golden)
        reduction = baseline_error - error
        output.append(
            {
                "bias_V": bias,
                "variant": variant,
                "source_A_per_um": sources[variant],
                "sentaurus_source_A_per_um": golden,
                "source_to_sentaurus_ratio": sources[variant] / golden,
                "abs_log10_error_dex": error,
                "error_reduction_from_vela_dex": reduction,
                "fraction_of_all_sent_error_removal": (
                    reduction / full_removal if full_removal else None
                ),
            }
        )
    return output


def classify(
    source_rows: list[dict[str, object]],
    newton_rows: list[dict[str, object]],
) -> dict[str, object]:
    family_variants = ("sent_psi_only", "sent_qfp_only", "sent_density_only")
    per_bias: dict[str, dict[str, float]] = {}
    wins = {variant: 0 for variant in family_variants}
    for bias in DEFAULT_BIASES:
        selected = {
            str(row["variant"]): float(row["fraction_of_all_sent_error_removal"])
            for row in source_rows
            if float(row["bias_V"]) == bias
        }
        winner = max(family_variants, key=lambda variant: selected[variant])
        wins[winner] += 1
        per_bias[bias_key(bias)] = {variant: selected[variant] for variant in family_variants}
    at_m20 = per_bias["-20"]
    dominant = max(family_variants, key=lambda variant: at_m20[variant])
    dominant_fraction = at_m20[dominant]
    if dominant_fraction >= 0.60 and wins[dominant] >= 3:
        source_outcome = dominant.replace("sent_", "").replace("_only", "") + "_dominant"
    else:
        source_outcome = "multi_family_interaction"

    density_steps = [
        row for row in newton_rows
        if row["variant"] == "feedback_density_only"
    ]
    density_m20 = next(row for row in density_steps if float(row["bias_V"]) == -20.0)
    qfp_projection = density_m20["qfp_target_projection_fraction"]
    if qfp_projection is None:
        update_outcome = "density_update_direction_undefined"
    elif float(qfp_projection) > 0.0:
        update_outcome = "density_feedback_moves_qfp_toward_sentaurus"
    elif float(qfp_projection) < 0.0:
        update_outcome = "density_feedback_moves_qfp_away_from_sentaurus"
    else:
        update_outcome = "density_feedback_qfp_neutral"
    return {
        "typed_outcome": f"{source_outcome}__{update_outcome}",
        "source_outcome": source_outcome,
        "first_update_outcome": update_outcome,
        "source_recovery_fraction_by_bias": per_bias,
        "source_family_win_count": wins,
        "minus20_dominant_family": dominant,
        "minus20_dominant_recovery_fraction": dominant_fraction,
        "minus20_density_feedback_qfp_projection_fraction": qfp_projection,
        "decision_rule": {
            "source_dominant": (
                "one family recovers at least 60% of the all-Sentaurus error "
                "removal at -20 V and wins at least three of four biases"
            ),
            "first_update_direction": (
                "sign of density-only first-step projection onto the Vela-to-"
                "Sentaurus QFP target at -20 V"
            ),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--doping", type=Path, required=True)
    parser.add_argument("--base-config", type=Path, required=True)
    parser.add_argument("--sentaurus-frozen-root", type=Path, required=True)
    parser.add_argument("--vela-manifest", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--bias", type=float, action="append", default=[])
    parser.add_argument("--repeats", type=int, default=2)
    args = parser.parse_args()
    biases = tuple(args.bias) if args.bias else DEFAULT_BIASES
    if biases != DEFAULT_BIASES:
        raise RuntimeError(f"bias lattice must be {DEFAULT_BIASES}")
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
    previous = read_json(args.sentaurus_frozen_root / "result.json")
    previous_comparison = read_rows(
        args.sentaurus_frozen_root / "source_comparison.csv"
    )
    golden = {
        (float(row["bias_V"]), row["carrier"]): float(
            row["sentaurus_source_A_per_um"]
        )
        for row in previous_comparison
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

    states: dict[float, dict[str, list[dict[str, object]]]] = {}
    for bias in biases:
        record = exact_bias_record(manifest, bias)
        vela_path = manifest_root / record["snapshot_tdr"]["path"]
        sent_path = args.sentaurus_frozen_root / f"sentaurus_state_{bias_tag(bias)}.csv"
        vela = normalize_vela_state(vela_path)
        sentaurus = normalize_sentaurus_state(sent_path)
        states[bias] = {
            variant: mixed_state(vela, sentaurus, variant)
            for variant in FIXED_VARIANTS
        }
        input_root = args.output_root / "inputs" / bias_tag(bias)
        for variant, state in states[bias].items():
            write_rows(input_root / f"{variant}.csv", state)
            if variant in INDEPENDENT_NEWTON_VARIANTS:
                write_external_fields(input_root / f"{variant}_fields", state)
        write_feedback_fields(input_root / "sentaurus_feedback_fields", sentaurus)

    run_source_rows: dict[str, list[dict[str, object]]] = {}
    run_newton_rows: dict[str, list[dict[str, object]]] = {}
    artifact_hashes: dict[str, dict[str, str]] = {}
    for repeat in range(args.repeats):
        run_label = chr(ord("a") + repeat)
        source_output: list[dict[str, object]] = []
        newton_output: list[dict[str, object]] = []
        run_hashes: dict[str, str] = {}
        for bias in biases:
            tag = bias_tag(bias)
            source_values: dict[str, float] = {}
            for variant in FIXED_VARIANTS:
                audit_root = args.output_root / f"run-{run_label}" / tag / "fixed" / variant
                state_path = args.output_root / "inputs" / tag / f"{variant}.csv"
                paths = run_audit(
                    args.audit,
                    args.mesh,
                    args.doping,
                    state_path,
                    fixed_config_path,
                    audit_root,
                )
                summary = summarize_audit(paths, states[bias][variant])
                source_values[variant] = float(summary["source_A_per_um"]["total"])
                for name, path in paths.items():
                    run_hashes[f"{tag}/fixed/{variant}/{name}"] = sha256(path)
            source_output.extend(
                source_rows_for_bias(bias, source_values, golden[(bias, "total")])
            )

            for variant in INDEPENDENT_NEWTON_VARIANTS:
                case = args.output_root / f"run-{run_label}" / tag / "newton" / variant
                case.mkdir(parents=True, exist_ok=True)
                output_csv = case / "newton_step_probe.csv"
                config_path = case / "newton_step_probe.json"
                fields = args.output_root / "inputs" / tag / f"{variant}_fields"
                config = make_probe_config(
                    args.base_config,
                    output_csv,
                    fields,
                    "newton_step_probe",
                    bias,
                    "Anode",
                    "Cathode",
                )
                write_json(config_path, config)
                status = run_runner(args.runner, config_path)
                write_json(case / "status.json", status)
                newton_output.append(
                    summarize_independent_step(
                        bias=bias,
                        variant=variant,
                        state=states[bias][variant],
                        sentaurus=states[bias]["sent_all"],
                        output_csv=output_csv,
                        status=status,
                        interior=interior,
                    )
                )
                run_hashes[f"{tag}/newton/{variant}"] = sha256(output_csv)

            feedback_case = args.output_root / f"run-{run_label}" / tag / "feedback"
            feedback_case.mkdir(parents=True, exist_ok=True)
            feedback_csv = feedback_case / "feedback_state_substitution.csv"
            feedback_config_path = feedback_case / "feedback_state_substitution.json"
            baseline_fields = args.output_root / "inputs" / tag / "vela_baseline_fields"
            feedback_fields = args.output_root / "inputs" / tag / "sentaurus_feedback_fields"
            feedback_config = make_probe_config(
                args.base_config,
                feedback_csv,
                baseline_fields,
                "newton_feedback_substitution_probe",
                bias,
                "Anode",
                "Cathode",
            )
            feedback_config["feedback_state_fields_dir"] = str(
                feedback_fields.resolve()
            )
            write_json(feedback_config_path, feedback_config)
            feedback_status = run_runner(args.runner, feedback_config_path)
            write_json(feedback_case / "status.json", feedback_status)
            newton_output.extend(
                summarize_feedback_steps(
                    bias=bias,
                    baseline=states[bias]["vela_baseline"],
                    sentaurus=states[bias]["sent_all"],
                    output_csv=feedback_csv,
                    status=feedback_status,
                    interior=interior,
                )
            )
            run_hashes[f"{tag}/feedback"] = sha256(feedback_csv)

        run_source_rows[run_label] = source_output
        run_newton_rows[run_label] = newton_output
        artifact_hashes[run_label] = run_hashes

    write_rows(args.output_root / "source_substitution.csv", run_source_rows["a"])
    write_rows(args.output_root / "newton_first_update.csv", run_newton_rows["a"])
    keys = sorted(artifact_hashes["a"])
    determinism_rows = [
        {
            "artifact": key,
            "repeat_count": args.repeats,
            "unique_hash_count": len(
                {artifact_hashes[label][key] for label in artifact_hashes}
            ),
            "byte_identical": int(
                len({artifact_hashes[label][key] for label in artifact_hashes}) == 1
            ),
            "sha256": artifact_hashes["a"][key],
        }
        for key in keys
    ]
    write_rows(args.output_root / "determinism.csv", determinism_rows)

    verdict = classify(run_source_rows["a"], run_newton_rows["a"])
    result = {
        "schema": "vela.pn2d_bv_m2_single_family_state_substitution.v1",
        "status": "passed",
        "biases_V": list(biases),
        "physics_modified": False,
        "production_defaults_modified": False,
        "acceptance_thresholds_modified": False,
        "fixed_state_coupling_mode": "postprocess_only",
        "newton_contract": {
            "independent_state_variants": list(INDEPENDENT_NEWTON_VARIANTS),
            "feedback_variants": list(FEEDBACK_VARIANTS),
            "density_semantics": (
                "density is a frozen residual-operator input under the single "
                "production baseline Jacobian, including Poisson charge, mobility, "
                "recombination, and avalanche paths; it is not an independent Newton unknown"
            ),
            "contact_rows": "production boundary rows preserved",
        },
        "source_rows": run_source_rows["a"],
        "newton_rows": run_newton_rows["a"],
        "verdict": verdict,
        "determinism": {
            "repeat_count": args.repeats,
            "artifact_count": len(determinism_rows),
            "all_byte_identical": all(
                row["byte_identical"] == 1 for row in determinism_rows
            ),
        },
        "input_hashes": {
            "audit": sha256(args.audit),
            "runner": sha256(args.runner),
            "mesh": sha256(args.mesh),
            "doping": sha256(args.doping),
            "base_config": sha256(args.base_config),
            "vela_manifest": sha256(args.vela_manifest),
            "prior_frozen_result": sha256(args.sentaurus_frozen_root / "result.json"),
            "prior_source_comparison": sha256(
                args.sentaurus_frozen_root / "source_comparison.csv"
            ),
        },
        "prior_frozen_contract": {
            "schema": previous["schema"],
            "typed_outcome": previous["verdict"]["typed_outcome"],
        },
    }
    write_json(args.output_root / "result.json", result)
    print(json.dumps(verdict, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
