#!/usr/bin/env python3
"""Replay selected M2 Sentaurus BV states through Vela SG/Laux.

This is an observation-only discriminating experiment.  The supplied
Sentaurus potential, quasi-Fermi potentials, and carrier densities are written
directly to Vela's fixed-state audit format.  Impact ionization is forced to
``postprocess_only`` so neither carrier continuity equation is solved or fed
back into a later bias point.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
from pathlib import Path


Q_C = 1.602176634e-19
DEFAULT_BIASES = (-18.0, -19.5, -19.7, -20.0)
STATE_FIELDS = {
    "psi_V": ("potential", "none", "V", 1.0),
    "phin_V": ("quasi_fermi", "electron", "V", 1.0),
    "phip_V": ("quasi_fermi", "hole", "V", 1.0),
    "n_m3": ("density", "electron", "cm^-3", 1.0e6),
    "p_m3": ("density", "hole", "cm^-3", 1.0e6),
}


def bias_key(value: float) -> str:
    return f"{value:.12g}"


def bias_tag(value: float) -> str:
    sign = "m" if value < 0.0 else "p"
    return sign + f"{abs(value):g}".replace(".", "p")


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8-sig"))


def write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2, sort_keys=True), encoding="utf-8")


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise RuntimeError(f"refusing to write empty CSV: {path}")
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


def select_aggregate(
    records: list[dict],
    *,
    branch: str,
    bias: float,
    quantity: str,
    carrier: str,
    provenance: str,
) -> float:
    matches = [
        record
        for record in records
        if record.get("branch") == branch
        and math.isclose(
            float(record.get("requested_bias_V", math.nan)),
            bias,
            rel_tol=0.0,
            abs_tol=1.0e-10,
        )
        and record.get("quantity") == quantity
        and record.get("carrier") == carrier
        and record.get("provenance") == provenance
        and record.get("unit") == "A/um"
    ]
    if len(matches) != 1:
        raise RuntimeError(
            "expected one aggregate record for "
            f"{branch} {bias:g} V {quantity}/{carrier}/{provenance}, "
            f"found {len(matches)}"
        )
    return float(matches[0]["value"])


def extract_sentaurus_inputs(
    manifest: dict,
    mesh: dict,
    biases: tuple[float, ...],
) -> tuple[dict[str, list[dict[str, object]]], dict[str, dict[str, float]], dict]:
    mesh_coordinates = {
        int(node["id"]): (float(node["x"]), float(node["y"]))
        for node in mesh["nodes"]
    }
    requested = {bias_key(bias): bias for bias in biases}
    state_records: dict[
        tuple[str, str, str], dict[int, dict]
    ] = {}
    physical_node_coordinates: dict[str, dict[int, tuple[float, float]]] = {
        key: {} for key in requested
    }

    for record in manifest["field_records"]:
        if record.get("branch") != "avalanche_on":
            continue
        key = bias_key(float(record.get("requested_bias_V", math.nan)))
        if key not in requested or record.get("support_kind") != "physical_node":
            continue
        quantity = str(record.get("quantity"))
        carrier = str(record.get("carrier"))
        if (quantity, carrier) not in {
            (source[0], source[1]) for source in STATE_FIELDS.values()
        }:
            continue
        support = str(record["support_key"])
        if not support.startswith("node:"):
            raise RuntimeError(f"unexpected physical-node support key: {support}")
        node = int(support.split(":", 1)[1])
        state_records.setdefault((key, quantity, carrier), {})[node] = record
        physical_node_coordinates[key][node] = tuple(
            float(value) for value in record["coordinates_um"][:2]
        )

    states: dict[str, list[dict[str, object]]] = {}
    mappings: dict[str, dict[str, float]] = {}
    for key, bias in requested.items():
        columns: dict[str, dict[int, float]] = {}
        max_coordinate_error = 0.0
        for output_name, (quantity, carrier, unit, scale) in STATE_FIELDS.items():
            records = state_records.get((key, quantity, carrier), {})
            if not records:
                raise RuntimeError(
                    f"missing Sentaurus {quantity}/{carrier} at {bias:g} V"
                )
            units = {record.get("unit") for record in records.values()}
            if units != {unit}:
                raise RuntimeError(
                    f"unexpected units for {quantity}/{carrier} at {bias:g} V: {units}"
                )
            missing = sorted(set(mesh_coordinates) - set(records))
            if missing:
                raise RuntimeError(
                    f"Sentaurus state at {bias:g} V is missing Vela nodes {missing}"
                )
            columns[output_name] = {
                node: float(records[node]["values"][0]) * scale
                for node in mesh_coordinates
            }
            for node, coordinate in mesh_coordinates.items():
                imported_coordinate = tuple(
                    float(value) for value in records[node]["coordinates_um"][:2]
                )
                max_coordinate_error = max(
                    max_coordinate_error,
                    math.dist(coordinate, imported_coordinate),
                )
        if max_coordinate_error > 1.0e-12:
            raise RuntimeError(
                f"maximum M2 coordinate mismatch at {bias:g} V is "
                f"{max_coordinate_error:.17g} um"
            )
        extra_nodes = sorted(
            set(physical_node_coordinates[key]) - set(mesh_coordinates)
        )
        duplicate_extras = 0
        mesh_coordinate_set = set(mesh_coordinates.values())
        for node in extra_nodes:
            if physical_node_coordinates[key][node] in mesh_coordinate_set:
                duplicate_extras += 1
        states[key] = [
            {
                "node_id": node,
                **{name: columns[name][node] for name in STATE_FIELDS},
            }
            for node in sorted(mesh_coordinates)
        ]
        mappings[key] = {
            "bias_V": bias,
            "vela_node_count": len(mesh_coordinates),
            "sentaurus_physical_record_count": len(
                physical_node_coordinates[key]
            ),
            "excluded_extra_node_count": len(extra_nodes),
            "excluded_coordinate_duplicate_count": duplicate_extras,
            "maximum_coordinate_mismatch_um": max_coordinate_error,
        }

    sentaurus_sources: dict[str, dict[str, float]] = {}
    for key, bias in requested.items():
        sentaurus_sources[key] = {
            carrier: select_aggregate(
                manifest["aggregate_records"],
                branch="avalanche_on",
                bias=bias,
                quantity="integrated_source",
                carrier=carrier,
                provenance="native",
            )
            for carrier in ("electron", "hole", "total")
        }

    mapping_summary = {
        "contact_support_policy": (
            "use physical nodes whose IDs and coordinates match the Vela M2 mesh; "
            "exclude duplicated Sentaurus contact-support vertices"
        ),
        "per_bias": mappings,
    }
    return states, sentaurus_sources, mapping_summary


def extract_vela_sources(
    manifest: dict, biases: tuple[float, ...]
) -> dict[str, dict[str, float]]:
    return {
        bias_key(bias): {
            carrier: select_aggregate(
                manifest["aggregate_records"],
                branch="avalanche_on",
                bias=bias,
                quantity="integrated_source",
                carrier=carrier,
                provenance="solver_used",
            )
            for carrier in ("electron", "hole", "total")
        }
        for bias in biases
    }


def run_audit(
    audit: Path,
    mesh: Path,
    doping: Path,
    state: Path,
    config: Path,
    output: Path,
) -> dict[str, Path]:
    output.mkdir(parents=True, exist_ok=True)
    paths = {
        "node": output / "node_state.csv",
        "edge": output / "edge_audit.csv",
        "triangle": output / "triangle_audit.csv",
        "element": output / "element_edge_sg_gss_laux.csv",
        "process": output / "bv_process_probe.csv",
    }
    command = [
        str(audit.resolve()),
        "--mesh", str(mesh.resolve()),
        "--doping", str(doping.resolve()),
        "--state", str(state.resolve()),
        "--config", str(config.resolve()),
        "--node-out", str(paths["node"].resolve()),
        "--edge-out", str(paths["edge"].resolve()),
        "--triangle-out", str(paths["triangle"].resolve()),
        "--element-out", str(paths["element"].resolve()),
        "--process-out", str(paths["process"].resolve()),
        "--scope", "general_tri3",
    ]
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    write_json(output / "command.json", command)
    (output / "stdout.log").write_text(completed.stdout, encoding="utf-8")
    (output / "stderr.log").write_text(completed.stderr, encoding="utf-8")
    if completed.returncode:
        raise RuntimeError(
            f"fixed-state audit failed with exit code {completed.returncode}: "
            f"{completed.stderr[-2000:]}"
        )
    return paths


def summarize_audit(
    paths: dict[str, Path], state_rows: list[dict[str, object]]
) -> dict[str, object]:
    expected = {int(row["node_id"]): row for row in state_rows}
    roundtrip = read_rows(paths["node"])
    if {int(row["node_id"]) for row in roundtrip} != set(expected):
        raise RuntimeError("fixed-state audit changed the node ID set")
    maximum_state_relative_error = 0.0
    maximum_state_absolute_error = 0.0
    for row in roundtrip:
        node = int(row["node_id"])
        for column in STATE_FIELDS:
            observed = float(row[column])
            reference = float(expected[node][column])
            maximum_state_absolute_error = max(
                maximum_state_absolute_error, abs(observed - reference)
            )
            if reference != 0.0:
                maximum_state_relative_error = max(
                    maximum_state_relative_error,
                    abs(observed - reference) / abs(reference),
                )

    process = read_rows(paths["process"])
    source_per_m_s = {"electron": 0.0, "hole": 0.0}
    qg_A_per_m = {"electron": 0.0, "hole": 0.0}
    alpha_peak_per_m = {"electron": 0.0, "hole": 0.0}
    hotspots: dict[str, dict[str, object] | None] = {
        "electron": None,
        "hole": None,
    }
    solver_coupled_count = 0
    residual_feedback_count = 0
    for row in process:
        carrier = row["carrier"]
        solver_coupled_count += int(row["solver_coupled"] or 0)
        for column in (
            "electron_residual_contributions_per_m_s",
            "hole_residual_contributions_per_m_s",
        ):
            if any(
                float(token) != 0.0
                for token in row[column].split(";")
                if token
            ):
                residual_feedback_count += 1
        if carrier not in source_per_m_s:
            continue
        alpha_peak_per_m[carrier] = max(
            alpha_peak_per_m[carrier], float(row["alpha_per_m"] or 0.0)
        )
        if not row["source_integral_per_m_s"]:
            continue
        source = float(row["source_integral_per_m_s"])
        qg = float(row["qG_contribution_A_per_m"])
        source_per_m_s[carrier] += source
        qg_A_per_m[carrier] += qg
        if hotspots[carrier] is None or qg > float(hotspots[carrier]["qG_A_per_m"]):
            hotspots[carrier] = {
                "cell_id": int(row["cell_id"]),
                "node_id": int(row["node0"]),
                "alpha_per_m": float(row["alpha_per_m"]),
                "selected_flux_per_m2_s": float(
                    row["selected_flux_magnitude_per_m2_s"]
                ),
                "generation_rate_per_m3_s": float(
                    row["generation_rate_per_m3_s"]
                ),
                "qG_A_per_m": qg,
            }

    source_A_per_um = {
        carrier: value / 1.0e6 for carrier, value in qg_A_per_m.items()
    }
    source_A_per_um["total"] = sum(source_A_per_um.values())
    qg_closure = max(
        abs(qg_A_per_m[carrier] - Q_C * source_per_m_s[carrier])
        / max(abs(qg_A_per_m[carrier]), 1.0e-300)
        for carrier in source_per_m_s
    )
    return {
        "state_node_count": len(roundtrip),
        "maximum_state_relative_error": maximum_state_relative_error,
        "maximum_state_absolute_error": maximum_state_absolute_error,
        "process_record_count": len(process),
        "solver_coupled_record_count": solver_coupled_count,
        "nonzero_residual_feedback_record_count": residual_feedback_count,
        "source_A_per_um": source_A_per_um,
        "source_particles_per_m_s": source_per_m_s,
        "maximum_qG_closure_relative_error": qg_closure,
        "alpha_peak_per_m": alpha_peak_per_m,
        "hotspots": hotspots,
        "configuration_fingerprints": sorted(
            {row["configuration_fingerprint"] for row in process}
        ),
    }


def log_error(candidate: float, reference: float) -> float:
    return abs(math.log10(abs(candidate) / abs(reference)))


def classify(rows: list[dict[str, object]]) -> dict[str, object]:
    totals = [row for row in rows if row["carrier"] == "total"]
    frozen_errors = [float(row["frozen_abs_log10_error_dex"]) for row in totals]
    self_errors = [float(row["self_consistent_abs_log10_error_dex"]) for row in totals]
    improvements = [before - after for before, after in zip(self_errors, frozen_errors)]
    mean_frozen = sum(frozen_errors) / len(frozen_errors)
    mean_self = sum(self_errors) / len(self_errors)
    mean_improvement = sum(improvements) / len(improvements)
    if mean_frozen <= 0.02 and mean_improvement >= 0.02:
        outcome = "state_feedback_dominant"
    elif mean_improvement >= 0.02:
        outcome = "mixed_state_and_operator"
    elif mean_frozen >= mean_self - 0.01:
        outcome = "operator_or_support_dominant"
    else:
        outcome = "inconclusive_partial_state_effect"
    return {
        "typed_outcome": outcome,
        "mean_frozen_abs_log10_error_dex": mean_frozen,
        "mean_self_consistent_abs_log10_error_dex": mean_self,
        "mean_error_reduction_dex": mean_improvement,
        "maximum_frozen_abs_log10_error_dex": max(frozen_errors),
        "maximum_self_consistent_abs_log10_error_dex": max(self_errors),
        "decision_rule": {
            "state_feedback_dominant": (
                "mean frozen-state error <= 0.02 dex and mean reduction >= 0.02 dex"
            ),
            "mixed_state_and_operator": (
                "mean reduction >= 0.02 dex but frozen-state error remains > 0.02 dex"
            ),
            "operator_or_support_dominant": (
                "frozen-state mean error is no more than 0.01 dex below the "
                "self-consistent mean error"
            ),
        },
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--doping", type=Path, required=True)
    parser.add_argument("--baseline-config", type=Path, required=True)
    parser.add_argument("--sentaurus-manifest", type=Path, required=True)
    parser.add_argument("--vela-manifest", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--bias", type=float, action="append", default=[])
    parser.add_argument("--repeats", type=int, default=2)
    args = parser.parse_args()

    biases = tuple(args.bias) if args.bias else DEFAULT_BIASES
    if args.repeats < 2:
        raise RuntimeError("at least two repeats are required for determinism evidence")
    args.out_dir.mkdir(parents=True, exist_ok=True)

    mesh = read_json(args.mesh)
    sentaurus_manifest = read_json(args.sentaurus_manifest)
    states, sentaurus_sources, mapping = extract_sentaurus_inputs(
        sentaurus_manifest, mesh, biases
    )
    del sentaurus_manifest
    vela_sources = extract_vela_sources(read_json(args.vela_manifest), biases)

    baseline = read_json(args.baseline_config)
    config = json.loads(json.dumps(baseline))
    impact = config["solver"]["impact_ionization"]
    impact["coupling_mode"] = "postprocess_only"
    required_operator = {
        "model": "van_overstraeten",
        "generation": "current_density",
        "driving_force": "quasi_fermi_gradient",
        "current_approximation": "element_edge_sg_gss_laux",
        "source_mapping_mode": "element_vertex_box_measure",
    }
    for key, expected in required_operator.items():
        if impact.get(key) != expected:
            raise RuntimeError(
                f"baseline impact_ionization.{key}={impact.get(key)!r}; "
                f"expected {expected!r}"
            )
    config_path = args.out_dir / "postprocess_only_sg_laux.json"
    write_json(config_path, config)

    mapping_rows = [mapping["per_bias"][bias_key(bias)] for bias in biases]
    write_rows(args.out_dir / "state_mapping.csv", mapping_rows)
    run_summaries: dict[str, list[dict[str, object]]] = {}
    determinism_rows: list[dict[str, object]] = []

    for bias in biases:
        key = bias_key(bias)
        tag = bias_tag(bias)
        state_path = args.out_dir / f"sentaurus_state_{tag}.csv"
        write_rows(state_path, states[key])
        per_bias_runs = []
        run_hashes = []
        for repeat in range(args.repeats):
            label = chr(ord("a") + repeat)
            paths = run_audit(
                args.audit,
                args.mesh,
                args.doping,
                state_path,
                config_path,
                args.out_dir / f"run-{label}" / tag,
            )
            summary = summarize_audit(paths, states[key])
            summary["run"] = label
            per_bias_runs.append(summary)
            run_hashes.append({name: sha256(path) for name, path in paths.items()})
        run_summaries[key] = per_bias_runs
        reference_hashes = run_hashes[0]
        for name in reference_hashes:
            hashes = [values[name] for values in run_hashes]
            determinism_rows.append(
                {
                    "bias_V": bias,
                    "artifact": name,
                    "repeat_count": args.repeats,
                    "unique_hash_count": len(set(hashes)),
                    "byte_identical": int(len(set(hashes)) == 1),
                    "sha256": hashes[0] if len(set(hashes)) == 1 else ";".join(hashes),
                }
            )

    comparison_rows: list[dict[str, object]] = []
    for bias in biases:
        key = bias_key(bias)
        frozen = run_summaries[key][0]["source_A_per_um"]
        for carrier in ("electron", "hole", "total"):
            sentaurus = sentaurus_sources[key][carrier]
            self_consistent = vela_sources[key][carrier]
            replay = float(frozen[carrier])
            comparison_rows.append(
                {
                    "bias_V": bias,
                    "carrier": carrier,
                    "sentaurus_source_A_per_um": sentaurus,
                    "vela_self_consistent_source_A_per_um": self_consistent,
                    "vela_frozen_sentaurus_state_source_A_per_um": replay,
                    "self_consistent_to_sentaurus_ratio": self_consistent / sentaurus,
                    "frozen_to_sentaurus_ratio": replay / sentaurus,
                    "self_consistent_abs_log10_error_dex": log_error(
                        self_consistent, sentaurus
                    ),
                    "frozen_abs_log10_error_dex": log_error(replay, sentaurus),
                    "frozen_error_reduction_dex": log_error(
                        self_consistent, sentaurus
                    ) - log_error(replay, sentaurus),
                }
            )
    write_rows(args.out_dir / "source_comparison.csv", comparison_rows)
    write_rows(args.out_dir / "determinism.csv", determinism_rows)

    verdict = classify(comparison_rows)
    result = {
        "schema": "vela.pn2d_bv_m2_sentaurus_frozen_sg_laux.v1",
        "observation_only": True,
        "state_advanced": False,
        "continuity_feedback_enabled": False,
        "biases_V": list(biases),
        "operator_contract": required_operator,
        "coupling_mode": impact["coupling_mode"],
        "state_mapping": mapping,
        "runs": run_summaries,
        "source_comparison": comparison_rows,
        "determinism": {
            "repeat_count": args.repeats,
            "all_artifacts_byte_identical": all(
                row["byte_identical"] == 1 for row in determinism_rows
            ),
            "records": determinism_rows,
        },
        "verdict": verdict,
        "input_hashes": {
            "audit": sha256(args.audit),
            "mesh": sha256(args.mesh),
            "doping": sha256(args.doping),
            "baseline_config": sha256(args.baseline_config),
            "postprocess_only_config": sha256(config_path),
            "sentaurus_manifest": sha256(args.sentaurus_manifest),
            "vela_manifest": sha256(args.vela_manifest),
        },
    }
    write_json(args.out_dir / "result.json", result)
    print(json.dumps(verdict, indent=2))


if __name__ == "__main__":
    main()
