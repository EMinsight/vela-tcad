"""Run the deterministic Minimal6 internal-node QFP replacement experiment."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import statistics
import subprocess
from collections.abc import Mapping
from pathlib import Path

from .qfp_sg_replacement import (
    INTERNAL_NODE_IDS,
    absolute_log10_error,
    continuity_flux_from_current_proxy,
    density_sg_flux,
    qf_sg_flux,
    replace_internal_qfp,
    symmetric_relative_residual,
)


THERMAL_VOLTAGE_300K_V = 1.380649e-23 * 300.0 / 1.602176634e-19
SILICON_NI_300K_M3 = 1.0e16


def _quantile(values: list[float], fraction: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    position = (len(ordered) - 1) * fraction
    lower, upper = math.floor(position), math.ceil(position)
    if lower == upper:
        return ordered[lower]
    return ordered[lower] + (ordered[upper] - ordered[lower]) * (position - lower)


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
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


def _observation_index(
    path: Path,
) -> tuple[dict[tuple[str, str, float, int, str], float], tuple[tuple[str, float], ...]]:
    index: dict[tuple[str, str, float, int, str], float] = {}
    states: set[tuple[str, float]] = set()
    quantities = {
        "ElectrostaticPotential",
        "eQuasiFermiPotential",
        "hQuasiFermiPotential",
        "eDensity",
        "hDensity",
    }
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if (
                row["support_kind"] != "node"
                or row["component"] != "component0"
                or row["status"] != "valid"
                or row["quantity"] not in quantities
            ):
                continue
            key = (
                row["solver"],
                row["topology"],
                float(row["bias_V"]),
                int(row["support_id"]),
                row["quantity"],
            )
            if key in index:
                raise ValueError(f"duplicate observation {key}")
            index[key] = float(row["value_si"])
            states.add((row["topology"], float(row["bias_V"])))
    expected = {
        (topology, float(-bias))
        for topology in ("sketch", "mirror")
        for bias in range(1, 21)
    }
    if states != expected:
        raise ValueError("observation states differ from the exact 40-state contract")
    return index, tuple(sorted(states))


def _state(
    index: Mapping[tuple[str, str, float, int, str], float],
    solver: str,
    topology: str,
    bias: float,
) -> dict[int, dict[str, float]]:
    fields = {
        "psi_V": "ElectrostaticPotential",
        "phin_V": "eQuasiFermiPotential",
        "phip_V": "hQuasiFermiPotential",
        "n_m3": "eDensity",
        "p_m3": "hDensity",
    }
    result: dict[int, dict[str, float]] = {}
    for node in range(6):
        result[node] = {}
        for output, quantity in fields.items():
            key = (solver, topology, bias, node, quantity)
            if key not in index:
                raise ValueError(f"missing observation {key}")
            result[node][output] = index[key]
    return result


def _write_state(path: Path, state: Mapping[int, Mapping[str, float]]) -> None:
    rows = []
    for node in range(6):
        values = state[node]
        rows.append(
            {
                "node_id": node,
                "psi_V": format(values["psi_V"], ".17g"),
                "phin_V": format(values["phin_V"], ".17g"),
                "phip_V": format(values["phip_V"], ".17g"),
                "n_m3": format(values["n_m3"], ".17g"),
                "p_m3": format(values["p_m3"], ".17g"),
            }
        )
    _write_csv(path, rows)


def _run_baseline_audit(
    executable: Path,
    inverse_root: Path,
    output: Path,
    topology: str,
    bias: float,
    state: Mapping[int, Mapping[str, float]],
) -> dict[int, dict[str, str]]:
    label = f"m{abs(int(bias))}V"
    state_root = output / "baseline_replay" / topology / label
    state_root.mkdir(parents=True, exist_ok=True)
    state_path = state_root / "state.csv"
    node_path = state_root / "nodes.csv"
    edge_path = state_root / "edges.csv"
    triangle_path = state_root / "triangles.csv"
    _write_state(state_path, state)
    source = inverse_root / "vela" / "source"
    command = [
        str(executable),
        "--mesh",
        str(source / "topologies" / topology / "mesh.json"),
        "--doping",
        str(source / "topologies" / topology / "doping.csv"),
        "--state",
        str(state_path),
        "--config",
        str(source / "decks" / topology / f"{label}.json"),
        "--node-out",
        str(node_path),
        "--edge-out",
        str(edge_path),
        "--triangle-out",
        str(triangle_path),
    ]
    completed = subprocess.run(command, capture_output=True, text=True, check=False)
    if completed.returncode != 0:
        raise RuntimeError(
            f"operator audit failed for {topology} {bias:g} V: "
            f"{completed.stderr.strip()}"
        )
    with edge_path.open(newline="", encoding="utf-8") as handle:
        rows = {
            (int(row["node0"]), int(row["node1"])): row
            for row in csv.DictReader(handle)
        }
    if len(rows) != 9:
        raise ValueError(f"operator audit {topology} {bias:g} V did not emit 9 edges")
    return rows


def _load_current_edges(path: Path) -> dict[tuple[str, float, str, int], dict[str, str]]:
    result: dict[tuple[str, float, str, int], dict[str, str]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            key = (
                row["topology"],
                float(row["bias_V"]),
                row["carrier"],
                int(row["edge_id"]),
            )
            if key in result:
                raise ValueError(f"duplicate current edge {key}")
            result[key] = row
    if len(result) != 720:
        raise ValueError(f"current proxy input expected 720 rows, got {len(result)}")
    return result


def _summaries(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    output: list[dict[str, object]] = []
    variants = (
        "baseline",
        "electron_qfp",
        "hole_qfp",
        "both_qfp",
        "strict_frozen_density",
    )
    for carrier in ("electron", "hole"):
        carrier_rows = [
            row for row in rows if row["carrier"] == carrier and row["affected_edge"]
        ]
        baseline_by_key = {
            (row["topology"], row["bias_V"], row["edge_id"]): row
            for row in carrier_rows
            if row["variant"] == "baseline"
        }
        for variant in variants:
            selected = [row for row in carrier_rows if row["variant"] == variant]
            log_errors = [
                float(row["abs_log10_error"])
                for row in selected
                if row["abs_log10_error"] != ""
            ]
            residuals = [float(row["symmetric_relative_residual"]) for row in selected]
            signs = [
                float(row["sign_agreement"])
                for row in selected
                if row["sign_agreement"] != ""
            ]
            improvements = []
            for row in selected:
                baseline = baseline_by_key[
                    (row["topology"], row["bias_V"], row["edge_id"])
                ]
                if row["abs_log10_error"] != "" and baseline["abs_log10_error"] != "":
                    improvements.append(
                        float(baseline["abs_log10_error"])
                        - float(row["abs_log10_error"])
                    )
            output.append(
                {
                    "carrier": carrier,
                    "variant": variant,
                    "affected_edge_count": len(selected),
                    "median_abs_log10_error_dex": _quantile(log_errors, 0.5),
                    "p95_abs_log10_error_dex": _quantile(log_errors, 0.95),
                    "median_symmetric_relative_residual": _quantile(residuals, 0.5),
                    "p95_symmetric_relative_residual": _quantile(residuals, 0.95),
                    "sign_agreement_fraction": statistics.fmean(signs) if signs else None,
                    "median_paired_log_error_improvement_dex": _quantile(
                        improvements, 0.5
                    ),
                    "p95_paired_log_error_improvement_dex": _quantile(
                        improvements, 0.95
                    ),
                }
            )
    return output


def _markdown(manifest: Mapping[str, object], summaries: list[dict[str, object]]) -> str:
    replay = manifest["baseline_cpp_replay"]
    assert isinstance(replay, Mapping)
    lines = [
        "# PN2D Minimal6 internal-QFP SG replacement",
        "",
        "## Contract",
        "",
        "- Exact states: 40 (2 topologies x reverse biases -1..-20 V).",
        "- Replaced: Sentaurus electron and/or hole QFP at internal nodes 1 and 5.",
        "- Frozen: Vela electrostatic potential, stored n/p, baseline edge mobility, mesh, 300 K, and intrinsic density.",
        "- Reference: endpoint-mean Sentaurus node current vector projected on the canonical edge tangent. This is a support-mapped proxy, not a native Sentaurus edge flux.",
        "- Negative control: density-form SG with frozen Vela n/p; QFP replacement cannot affect it by construction.",
        "",
        "## Baseline replay gate",
        "",
        f"- Samples: {replay['sample_count']}",
        f"- Maximum relative difference from C++: {float(replay['max_relative_error']):.6g}",
        f"- Gate: {replay['status']} at <= {float(replay['tolerance_max_relative_error']):.6g}",
        "",
        "## Affected-edge results",
        "",
        "| Carrier | Variant | n | Median error (dex) | p95 error (dex) | Median bounded residual | Sign agreement | Paired median improvement (dex) |",
        "|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in summaries:
        lines.append(
            "| {carrier} | {variant} | {affected_edge_count} | "
            "{median_abs_log10_error_dex:.6g} | {p95_abs_log10_error_dex:.6g} | "
            "{median_symmetric_relative_residual:.6g} | {sign_agreement_fraction:.6g} | "
            "{median_paired_log_error_improvement_dex:.6g} |".format(**row)
        )
    return "\n".join(lines) + "\n"


def run_qfp_replacement_experiment(
    *,
    observations_csv: str | Path,
    current_edges_csv: str | Path,
    inverse_inputs_root: str | Path,
    operator_audit_executable: str | Path,
    output_root: str | Path,
) -> dict[str, object]:
    observations_path = Path(observations_csv).resolve()
    current_path = Path(current_edges_csv).resolve()
    inverse_root = Path(inverse_inputs_root).resolve()
    executable = Path(operator_audit_executable).resolve()
    output = Path(output_root).resolve()
    output.mkdir(parents=True, exist_ok=True)
    observations, states = _observation_index(observations_path)
    current_edges = _load_current_edges(current_path)
    flags = {
        "baseline": (False, False),
        "electron_qfp": (True, False),
        "hole_qfp": (False, True),
        "both_qfp": (True, True),
    }
    rows: list[dict[str, object]] = []
    replay_relative_errors: list[float] = []
    replay_absolute_errors: list[float] = []

    for topology, bias in states:
        vela = _state(observations, "vela", topology, bias)
        sentaurus = _state(observations, "sentaurus", topology, bias)
        baseline_edges = _run_baseline_audit(
            executable, inverse_root, output, topology, bias, vela
        )
        branches = {
            name: replace_internal_qfp(
                vela,
                sentaurus,
                replace_electron=value[0],
                replace_hole=value[1],
            )
            for name, value in flags.items()
        }
        for carrier in ("electron", "hole"):
            qf_key = "phin_V" if carrier == "electron" else "phip_V"
            density_key = "n_m3" if carrier == "electron" else "p_m3"
            cpp_key = (
                "electron_raw_signed_flux_per_m2_s"
                if carrier == "electron"
                else "hole_raw_signed_flux_per_m2_s"
            )
            for edge_id in range(9):
                current = current_edges[(topology, bias, carrier, edge_id)]
                node0, node1 = int(current["node0"]), int(current["node1"])
                length = float(current["length_m"])
                mobility = float(
                    current["vela_masetti_native_state_mobility_m2_per_Vs"]
                )
                coefficient = mobility * THERMAL_VOLTAGE_300K_V / length
                reference = continuity_flux_from_current_proxy(
                    carrier, float(current["current_tangent_A_per_m2"])
                )
                replay = qf_sg_flux(
                    carrier,
                    SILICON_NI_300K_M3,
                    SILICON_NI_300K_M3,
                    vela[node0]["psi_V"],
                    vela[node1]["psi_V"],
                    vela[node0][qf_key],
                    vela[node1][qf_key],
                    THERMAL_VOLTAGE_300K_V,
                    coefficient,
                )
                cpp_flux = float(baseline_edges[(node0, node1)][cpp_key])
                absolute = abs(replay - cpp_flux)
                replay_absolute_errors.append(absolute)
                if cpp_flux != 0.0:
                    replay_relative_errors.append(absolute / abs(cpp_flux))
                elif absolute != 0.0:
                    raise ValueError(
                        f"nonzero replay for zero C++ flux: {topology} {bias} "
                        f"{carrier} edge {edge_id}"
                    )
                branch_fluxes = {
                    name: qf_sg_flux(
                        carrier,
                        SILICON_NI_300K_M3,
                        SILICON_NI_300K_M3,
                        branch[node0]["psi_V"],
                        branch[node1]["psi_V"],
                        branch[node0][qf_key],
                        branch[node1][qf_key],
                        THERMAL_VOLTAGE_300K_V,
                        coefficient,
                    )
                    for name, branch in branches.items()
                }
                branch_fluxes["strict_frozen_density"] = density_sg_flux(
                    carrier,
                    vela[node0][density_key],
                    vela[node1][density_key],
                    vela[node0]["psi_V"],
                    vela[node1]["psi_V"],
                    THERMAL_VOLTAGE_300K_V,
                    coefficient,
                )
                for variant, candidate in branch_fluxes.items():
                    log_error = absolute_log10_error(candidate, reference)
                    sign = (
                        ""
                        if candidate == 0.0 or reference == 0.0
                        else float(
                            math.copysign(1.0, candidate)
                            == math.copysign(1.0, reference)
                        )
                    )
                    rows.append(
                        {
                            "topology": topology,
                            "bias_V": bias,
                            "carrier": carrier,
                            "edge_id": edge_id,
                            "node0": node0,
                            "node1": node1,
                            "affected_edge": int(
                                node0 in INTERNAL_NODE_IDS or node1 in INTERNAL_NODE_IDS
                            ),
                            "variant": variant,
                            "mobility_frozen_m2_per_Vs": mobility,
                            "candidate_continuity_flux_per_m2_s": candidate,
                            "sentaurus_endpoint_mean_current_proxy_continuity_flux_per_m2_s": reference,
                            "symmetric_relative_residual": symmetric_relative_residual(
                                candidate, reference
                            ),
                            "abs_log10_error": "" if log_error is None else log_error,
                            "sign_agreement": sign,
                        }
                    )

    maximum_relative = max(replay_relative_errors, default=0.0)
    maximum_absolute = max(replay_absolute_errors, default=0.0)
    tolerance = 5.0e-11
    if maximum_relative > tolerance:
        raise ValueError(
            f"offline baseline SG replay failed C++ agreement: {maximum_relative:.6g}"
        )
    summaries = _summaries(rows)
    edge_path = output / "qfp_replacement_edge_samples.csv"
    summary_path = output / "qfp_replacement_summary.csv"
    _write_csv(edge_path, rows)
    _write_csv(summary_path, summaries)
    manifest: dict[str, object] = {
        "schema_version": 1,
        "status": "valid",
        "experiment": "minimal6_internal_node_qfp_sg_replacement",
        "state_contract": {
            "state_count": len(states),
            "topologies": ["mirror", "sketch"],
            "biases_V": list(range(-20, 0)),
            "internal_node_ids": list(INTERNAL_NODE_IDS),
        },
        "frozen_fields": [
            "psi",
            "n",
            "p",
            "edge_mobility",
            "mesh",
            "temperature",
            "intrinsic_density",
        ],
        "replacement_fields": ["phin at nodes 1,5", "phip at nodes 1,5"],
        "reference_semantics": (
            "Sentaurus endpoint-mean node current-density vector projected onto "
            "the canonical edge tangent; support-mapped proxy, not native edge flux"
        ),
        "strict_frozen_density_control": (
            "density-form SG uses frozen Vela n,p and is QFP-independent by construction"
        ),
        "baseline_cpp_replay": {
            "sample_count": len(replay_absolute_errors),
            "max_relative_error": maximum_relative,
            "max_absolute_error_per_m2_s": maximum_absolute,
            "tolerance_max_relative_error": tolerance,
            "status": "passed",
        },
        "inputs": {
            "observations_csv": str(observations_path),
            "observations_sha256": _sha256(observations_path),
            "current_edges_csv": str(current_path),
            "current_edges_sha256": _sha256(current_path),
            "inverse_inputs_root": str(inverse_root),
            "operator_audit_executable": str(executable),
            "operator_audit_sha256": _sha256(executable),
        },
        "outputs": {
            "edge_samples_csv": edge_path.name,
            "summary_csv": summary_path.name,
        },
    }
    report_path = output / "report.md"
    report_path.write_text(_markdown(manifest, summaries), encoding="utf-8")
    manifest["outputs"] = {
        "edge_samples_csv": edge_path.name,
        "edge_samples_sha256": _sha256(edge_path),
        "summary_csv": summary_path.name,
        "summary_sha256": _sha256(summary_path),
        "report_md": report_path.name,
        "report_sha256": _sha256(report_path),
    }
    manifest_path = output / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest
