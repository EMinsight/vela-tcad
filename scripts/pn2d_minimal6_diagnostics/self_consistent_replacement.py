"""Self-consistent Minimal6 potential replacement and production-operator replay."""

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
    absolute_log10_error,
    continuity_flux_from_current_proxy,
    qf_sg_flux,
    symmetric_relative_residual,
)


THERMAL_VOLTAGE_300K_V = 1.380649e-23 * 300.0 / 1.602176634e-19
STAGE_ORDER = (
    "replace_psi_phin_phip",
    "recompute_n_p_from_vela_bgn",
    "recompute_mobility",
    "recompute_sg_current",
    "recompute_alpha_and_avalanche_source",
)


def _limited_exp(value: float) -> float:
    return math.exp(max(-500.0, min(500.0, float(value))))


def _effective_intrinsic_density(
    values: Mapping[str, float], carrier: str
) -> float:
    psi = float(values["psi_V"])
    if carrier == "electron":
        density = float(values["n_m3"])
        exponent = (psi - float(values["phin_V"])) / THERMAL_VOLTAGE_300K_V
    elif carrier == "hole":
        density = float(values["p_m3"])
        exponent = (float(values["phip_V"]) - psi) / THERMAL_VOLTAGE_300K_V
    else:
        raise ValueError(f"unsupported carrier {carrier!r}")
    if not math.isfinite(density) or density <= 0.0:
        raise ValueError("carrier density must be finite and positive")
    result = math.exp(math.log(density) - exponent)
    if not math.isfinite(result) or result <= 0.0:
        raise ValueError("effective intrinsic density must be finite and positive")
    return result


def replace_potentials_and_recompute_carriers(
    vela_state: Mapping[int, Mapping[str, float]],
    sentaurus_state: Mapping[int, Mapping[str, float]],
) -> tuple[dict[int, dict[str, float]], dict[int, float], float]:
    """Replace all three potentials and derive n/p from the frozen Vela BGN ni."""
    if set(vela_state) != set(range(6)) or set(sentaurus_state) != set(range(6)):
        raise ValueError("self-consistent replacement requires exactly nodes 0..5")
    output: dict[int, dict[str, float]] = {}
    effective_ni: dict[int, float] = {}
    ni_gaps_dex: list[float] = []
    for node in range(6):
        ni_e = _effective_intrinsic_density(vela_state[node], "electron")
        ni_h = _effective_intrinsic_density(vela_state[node], "hole")
        ni_gaps_dex.append(abs(math.log10(ni_e / ni_h)))
        ni = math.sqrt(ni_e * ni_h)
        effective_ni[node] = ni
        sent = sentaurus_state[node]
        psi = float(sent["psi_V"])
        phin = float(sent["phin_V"])
        phip = float(sent["phip_V"])
        output[node] = {
            "psi_V": psi,
            "phin_V": phin,
            "phip_V": phip,
            "n_m3": ni * _limited_exp((psi - phin) / THERMAL_VOLTAGE_300K_V),
            "p_m3": ni * _limited_exp((phip - psi) / THERMAL_VOLTAGE_300K_V),
        }
    return output, effective_ni, max(ni_gaps_dex, default=0.0)


def infer_edge_mobility_m2_per_Vs(
    *,
    carrier: str,
    node0: int,
    node1: int,
    length_m: float,
    signed_flux_per_m2_s: float,
    state: Mapping[int, Mapping[str, float]],
    effective_ni_m3: Mapping[int, float],
) -> float | None:
    """Invert the linear mobility factor in the production variable-ni SG flux."""
    qf_key = "phin_V" if carrier == "electron" else "phip_V"
    unit_flux = qf_sg_flux(
        carrier,
        effective_ni_m3[node0],
        effective_ni_m3[node1],
        state[node0]["psi_V"],
        state[node1]["psi_V"],
        state[node0][qf_key],
        state[node1][qf_key],
        THERMAL_VOLTAGE_300K_V,
        THERMAL_VOLTAGE_300K_V / float(length_m),
    )
    flux = float(signed_flux_per_m2_s)
    scale = max(abs(unit_flux), abs(flux), 1.0)
    if abs(unit_flux) <= 1.0e-14 * scale:
        if abs(flux) <= 1.0e-14 * scale:
            return None
        raise ValueError("nonzero SG flux has zero unit-mobility response")
    mobility = flux / unit_flux
    if not math.isfinite(mobility) or mobility <= 0.0:
        raise ValueError("inferred edge mobility must be finite and positive")
    return mobility


def _quantile(values: list[float], fraction: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    position = (len(ordered) - 1) * fraction
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    return ordered[lower] + (ordered[upper] - ordered[lower]) * (
        position - lower
    )


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty CSV {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _load_observations(
    path: Path,
) -> tuple[
    dict[tuple[str, str, float, int, str, str], float],
    tuple[tuple[str, float], ...],
]:
    index: dict[tuple[str, str, float, int, str, str], float] = {}
    states: set[tuple[str, float]] = set()
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row["support_kind"] != "node" or row["status"] != "valid":
                continue
            try:
                value = float(row["value_si"])
            except ValueError:
                continue
            key = (
                row["solver"],
                row["topology"],
                float(row["bias_V"]),
                int(row["support_id"]),
                row["quantity"],
                row["component"],
            )
            if key in index:
                raise ValueError(f"duplicate observation {key}")
            index[key] = value
            states.add((row["topology"], float(row["bias_V"])))
    expected = {
        (topology, float(-bias))
        for topology in ("mirror", "sketch")
        for bias in range(1, 21)
    }
    if states != expected:
        raise ValueError("observations differ from the exact 40-state contract")
    return index, tuple(sorted(states))


def _state(
    index: Mapping[tuple[str, str, float, int, str, str], float],
    solver: str,
    topology: str,
    bias: float,
) -> dict[int, dict[str, float]]:
    quantities = {
        "psi_V": "ElectrostaticPotential",
        "phin_V": "eQuasiFermiPotential",
        "phip_V": "hQuasiFermiPotential",
        "n_m3": "eDensity",
        "p_m3": "hDensity",
    }
    output: dict[int, dict[str, float]] = {}
    for node in range(6):
        output[node] = {}
        for field, quantity in quantities.items():
            key = (solver, topology, bias, node, quantity, "component0")
            if key not in index:
                raise ValueError(f"missing observation {key}")
            output[node][field] = index[key]
    return output


def _write_state(path: Path, state: Mapping[int, Mapping[str, float]]) -> None:
    rows: list[dict[str, object]] = []
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


def _run_operator_audit(
    *,
    executable: Path,
    inverse_root: Path,
    output: Path,
    branch: str,
    topology: str,
    bias: float,
    state: Mapping[int, Mapping[str, float]],
) -> tuple[
    dict[int, dict[str, str]],
    list[dict[str, str]],
    list[dict[str, str]],
]:
    label = f"m{abs(int(bias))}V"
    state_root = output / branch / topology / label
    state_root.mkdir(parents=True, exist_ok=True)
    state_path = state_root / "state.csv"
    node_path = state_root / "nodes.csv"
    edge_path = state_root / "edges.csv"
    triangle_path = state_root / "triangles.csv"
    _write_state(state_path, state)
    source = inverse_root / "vela" / "source"
    deck_path = source / "decks" / topology / f"{label}.json"
    deck = json.loads(deck_path.read_text(encoding="utf-8"))
    if "materials_file" in deck:
        materials = source / "topologies" / topology / "materials.json"
        if not materials.is_file():
            raise ValueError(f"missing sealed materials for {topology} {bias:g} V")
        deck["materials_file"] = str(materials.resolve())
    config_path = state_root / "audit_config.json"
    config_path.write_text(
        json.dumps(deck, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    command = [
        str(executable),
        "--mesh",
        str(source / "topologies" / topology / "mesh.json"),
        "--doping",
        str(source / "topologies" / topology / "doping.csv"),
        "--state",
        str(state_path),
        "--config",
        str(config_path),
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
            f"operator audit failed for {branch} {topology} {bias:g} V: "
            f"{completed.stderr.strip()}"
        )
    with edge_path.open(newline="", encoding="utf-8") as handle:
        edges = {int(row["edge_id"]): row for row in csv.DictReader(handle)}
    with node_path.open(newline="", encoding="utf-8") as handle:
        nodes = list(csv.DictReader(handle))
    with triangle_path.open(newline="", encoding="utf-8") as handle:
        triangles = list(csv.DictReader(handle))
    if len(edges) != 9 or len(nodes) != 6 or len(triangles) != 4:
        raise ValueError("operator audit did not emit the exact 6/9/4 contract")
    return edges, nodes, triangles


def _current_edge_key(
    topology: str,
    bias: float,
    carrier: str,
    node0: int,
    node1: int,
) -> tuple[str, float, str, int, int]:
    return topology, bias, carrier, min(node0, node1), max(node0, node1)


def _load_current_edges(
    path: Path,
) -> dict[tuple[str, float, str, int, int], dict[str, str]]:
    output: dict[tuple[str, float, str, int, int], dict[str, str]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            key = _current_edge_key(
                row["topology"],
                float(row["bias_V"]),
                row["carrier"],
                int(row["node0"]),
                int(row["node1"]),
            )
            if key in output:
                raise ValueError(f"duplicate current-edge reference {key}")
            output[key] = row
    if len(output) != 720:
        raise ValueError("current-edge reference must contain 720 carrier-edge rows")
    return output

def _metric_row(
    *,
    stage: str,
    metric: str,
    carrier: str,
    support: str,
    unit: str,
    values: list[float],
) -> dict[str, object]:
    finite = [value for value in values if math.isfinite(value)]
    return {
        "stage": stage,
        "metric": metric,
        "carrier": carrier,
        "support": support,
        "unit": unit,
        "sample_count": len(values),
        "finite_count": len(finite),
        "median": "" if not finite else statistics.median(finite),
        "p95": "" if not finite else _quantile(finite, 0.95),
        "maximum": "" if not finite else max(finite),
    }


def _markdown(
    manifest: Mapping[str, object], summaries: list[dict[str, object]]
) -> str:
    lines = [
        "# PN2D Minimal6 self-consistent potential replacement",
        "",
        "## Contract",
        "",
        "- Exact states: 40 (mirror/sketch, reverse biases -1 through -20 V).",
        "- Replace Sentaurus psi, electron QFP, and hole QFP at all six nodes.",
        "- Recover the node effective intrinsic density from the baseline Vela BGN state.",
        "- Recompute n/p from the replaced potentials; Sentaurus n/p are comparison targets only.",
        "- Re-run the C++ fixed-state production operators for mobility, SG current, alpha, and source.",
        "- Sentaurus edge current and alpha references remain endpoint support mappings, not native directed-edge quantities.",
        "",
        "## Ordered stages",
        "",
    ]
    for number, stage in enumerate(STAGE_ORDER, 1):
        lines.append(f"{number}. `{stage}`")
    lines.extend(
        [
            "",
            "## Summary",
            "",
            "| Stage | Metric | Carrier | Support | n | Median | p95 | Maximum | Unit |",
            "|---|---|---|---|---:|---:|---:|---:|---|",
        ]
    )
    for row in summaries:
        lines.append(
            "| {stage} | {metric} | {carrier} | {support} | {finite_count} | "
            "{median} | {p95} | {maximum} | {unit} |".format(**row)
        )
    lines.extend(
        [
            "",
            "## Interpretation guard",
            "",
            "This is a fixed-state counterfactual. It does not solve the nonlinear Vela equations from the replaced state. "
            "Current and alpha edge references are support-mapped diagnostics because native Sentaurus directed-edge flux is unavailable.",
            "",
            f"- Maximum baseline electron/hole recovered-ni discrepancy: {float(manifest['maximum_recovered_ni_gap_dex']):.6g} dex.",
        ]
    )
    return "\n".join(lines) + "\n"


def run_self_consistent_replacement_experiment(
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
    observations, states = _load_observations(observations_path)
    current_edges = _load_current_edges(current_path)

    node_rows: list[dict[str, object]] = []
    edge_rows: list[dict[str, object]] = []
    generation_rows: list[dict[str, object]] = []
    maximum_ni_gap = 0.0
    density_errors = {"electron_all": [], "electron_internal": [],
                      "hole_all": [], "hole_internal": []}
    mobility_errors = {"electron": [], "hole": []}
    current_errors = {"electron": [], "hole": []}
    current_sign = {"electron": [], "hole": []}
    alpha_errors = {"electron": [], "hole": []}
    generation_errors = {"baseline": [], "self_consistent": []}

    for topology, bias in states:
        vela = _state(observations, "vela", topology, bias)
        sentaurus = _state(observations, "sentaurus", topology, bias)
        replaced, effective_ni, ni_gap = replace_potentials_and_recompute_carriers(
            vela, sentaurus
        )
        maximum_ni_gap = max(maximum_ni_gap, ni_gap)
        baseline_edges, _, _ = _run_operator_audit(
            executable=executable,
            inverse_root=inverse_root,
            output=output,
            branch="baseline_replay",
            topology=topology,
            bias=bias,
            state=vela,
        )
        replaced_edges, _, replaced_triangles = _run_operator_audit(
            executable=executable,
            inverse_root=inverse_root,
            output=output,
            branch="self_consistent_replay",
            topology=topology,
            bias=bias,
            state=replaced,
        )

        for node in range(6):
            ni_e = _effective_intrinsic_density(vela[node], "electron")
            ni_h = _effective_intrinsic_density(vela[node], "hole")
            e_error = absolute_log10_error(
                replaced[node]["n_m3"], sentaurus[node]["n_m3"]
            )
            h_error = absolute_log10_error(
                replaced[node]["p_m3"], sentaurus[node]["p_m3"]
            )
            assert e_error is not None and h_error is not None
            density_errors["electron_all"].append(e_error)
            density_errors["hole_all"].append(h_error)
            if node in (1, 5):
                density_errors["electron_internal"].append(e_error)
                density_errors["hole_internal"].append(h_error)
            node_rows.append(
                {
                    "topology": topology,
                    "bias_V": bias,
                    "node_id": node,
                    "vela_effective_ni_electron_m3": ni_e,
                    "vela_effective_ni_hole_m3": ni_h,
                    "vela_effective_ni_used_m3": effective_ni[node],
                    "replaced_psi_V": replaced[node]["psi_V"],
                    "replaced_phin_V": replaced[node]["phin_V"],
                    "replaced_phip_V": replaced[node]["phip_V"],
                    "recomputed_n_m3": replaced[node]["n_m3"],
                    "sentaurus_n_m3": sentaurus[node]["n_m3"],
                    "electron_density_abs_log10_error_dex": e_error,
                    "recomputed_p_m3": replaced[node]["p_m3"],
                    "sentaurus_p_m3": sentaurus[node]["p_m3"],
                    "hole_density_abs_log10_error_dex": h_error,
                }
            )

        baseline_total = 0.0
        replaced_total = 0.0
        for carrier in ("electron", "hole"):
            flux_key = f"{carrier}_raw_signed_flux_per_m2_s"
            alpha_key = f"{carrier}_alpha_per_m"
            qf_quantity = (
                "eAlphaAvalanche" if carrier == "electron" else "hAlphaAvalanche"
            )
            for edge_id in range(9):
                base_edge = baseline_edges[edge_id]
                new_edge = replaced_edges[edge_id]
                node0 = int(new_edge["node0"])
                node1 = int(new_edge["node1"])
                reference = current_edges[
                    _current_edge_key(topology, bias, carrier, node0, node1)
                ]
                length = float(new_edge["length_m"])
                base_flux = float(base_edge[flux_key])
                new_flux = float(new_edge[flux_key])
                base_mobility = infer_edge_mobility_m2_per_Vs(
                    carrier=carrier,
                    node0=node0,
                    node1=node1,
                    length_m=length,
                    signed_flux_per_m2_s=base_flux,
                    state=vela,
                    effective_ni_m3=effective_ni,
                )
                new_mobility = infer_edge_mobility_m2_per_Vs(
                    carrier=carrier,
                    node0=node0,
                    node1=node1,
                    length_m=length,
                    signed_flux_per_m2_s=new_flux,
                    state=replaced,
                    effective_ni_m3=effective_ni,
                )
                sent_mobility = float(
                    reference["sentaurus_exported_mobility_m2_per_Vs"]
                )
                mobility_error = (
                    None
                    if new_mobility is None
                    else absolute_log10_error(new_mobility, sent_mobility)
                )
                if mobility_error is not None:
                    mobility_errors[carrier].append(mobility_error)
                sent_flux = continuity_flux_from_current_proxy(
                    carrier, float(reference["current_tangent_A_per_m2"])
                )
                current_error = absolute_log10_error(new_flux, sent_flux)
                if current_error is not None:
                    current_errors[carrier].append(current_error)
                sign_agreement = (
                    None
                    if new_flux == 0.0 or sent_flux == 0.0
                    else float(math.copysign(1.0, new_flux) == math.copysign(1.0, sent_flux))
                )
                if sign_agreement is not None:
                    current_sign[carrier].append(sign_agreement)
                sent_alpha0 = observations[
                    (
                        "sentaurus",
                        topology,
                        bias,
                        node0,
                        qf_quantity,
                        "component0",
                    )
                ]
                sent_alpha1 = observations[
                    (
                        "sentaurus",
                        topology,
                        bias,
                        node1,
                        qf_quantity,
                        "component0",
                    )
                ]
                sent_alpha = 0.5 * (sent_alpha0 + sent_alpha1)
                base_alpha = float(base_edge[alpha_key])
                new_alpha = float(new_edge[alpha_key])
                alpha_error = absolute_log10_error(new_alpha, sent_alpha)
                if alpha_error is not None:
                    alpha_errors[carrier].append(alpha_error)
                area = float(new_edge["edge_area_m2"])
                base_source = base_alpha * abs(base_flux) * area
                new_source = new_alpha * abs(new_flux) * area
                baseline_total += base_source
                replaced_total += new_source
                edge_rows.append(
                    {
                        "topology": topology,
                        "bias_V": bias,
                        "carrier": carrier,
                        "edge_id": edge_id,
                        "node0": node0,
                        "node1": node1,
                        "edge_length_m": length,
                        "edge_area_m2": area,
                        "baseline_mobility_m2_per_Vs": ""
                        if base_mobility is None
                        else base_mobility,
                        "recomputed_mobility_m2_per_Vs": ""
                        if new_mobility is None
                        else new_mobility,
                        "sentaurus_exported_mobility_m2_per_Vs": sent_mobility,
                        "mobility_abs_log10_error_dex": ""
                        if mobility_error is None
                        else mobility_error,
                        "baseline_sg_flux_per_m2_s": base_flux,
                        "recomputed_sg_flux_per_m2_s": new_flux,
                        "sentaurus_endpoint_current_proxy_flux_per_m2_s": sent_flux,
                        "current_symmetric_relative_residual": symmetric_relative_residual(
                            new_flux, sent_flux
                        ),
                        "current_abs_log10_error_dex": ""
                        if current_error is None
                        else current_error,
                        "current_sign_agreement": ""
                        if sign_agreement is None
                        else sign_agreement,
                        "baseline_alpha_per_m": base_alpha,
                        "recomputed_alpha_per_m": new_alpha,
                        "sentaurus_endpoint_mean_alpha_per_m": sent_alpha,
                        "alpha_abs_log10_error_dex": ""
                        if alpha_error is None
                        else alpha_error,
                        "baseline_source_integral_per_m_s": base_source,
                        "recomputed_source_integral_per_m_s": new_source,
                    }
                )

        sent_generation = 0.0
        for triangle in replaced_triangles:
            area = abs(float(triangle["signed_double_area_m2"])) * 0.5
            nodes = [int(triangle[f"node{i}"]) for i in range(3)]
            mean_generation = sum(
                observations[
                    (
                        "sentaurus",
                        topology,
                        bias,
                        node,
                        "ImpactIonization",
                        "component0",
                    )
                ]
                for node in nodes
            ) / 3.0
            sent_generation += mean_generation * area
        base_error = absolute_log10_error(baseline_total, sent_generation)
        new_error = absolute_log10_error(replaced_total, sent_generation)
        if base_error is not None:
            generation_errors["baseline"].append(base_error)
        if new_error is not None:
            generation_errors["self_consistent"].append(new_error)
        generation_rows.append(
            {
                "topology": topology,
                "bias_V": bias,
                "baseline_vela_source_integral_per_m_s": baseline_total,
                "self_consistent_vela_source_integral_per_m_s": replaced_total,
                "sentaurus_native_impact_integral_per_m_s": sent_generation,
                "baseline_abs_log10_error_dex": ""
                if base_error is None
                else base_error,
                "self_consistent_abs_log10_error_dex": ""
                if new_error is None
                else new_error,
            }
        )

    summaries: list[dict[str, object]] = []
    zero_potential_errors = [0.0] * (len(states) * 6)
    for field in ("psi", "phin", "phip"):
        summaries.append(
            _metric_row(
                stage=STAGE_ORDER[0],
                metric=f"{field}_absolute_error",
                carrier="combined",
                support="node",
                unit="V",
                values=zero_potential_errors,
            )
        )
    for carrier in ("electron", "hole"):
        summaries.append(
            _metric_row(
                stage=STAGE_ORDER[1],
                metric="density_abs_log10_error",
                carrier=carrier,
                support="all_nodes",
                unit="dex",
                values=density_errors[f"{carrier}_all"],
            )
        )
        summaries.append(
            _metric_row(
                stage=STAGE_ORDER[1],
                metric="density_abs_log10_error",
                carrier=carrier,
                support="internal_nodes_1_5",
                unit="dex",
                values=density_errors[f"{carrier}_internal"],
            )
        )
        summaries.append(
            _metric_row(
                stage=STAGE_ORDER[2],
                metric="mobility_abs_log10_error",
                carrier=carrier,
                support="identifiable_edges",
                unit="dex",
                values=mobility_errors[carrier],
            )
        )
        summaries.append(
            _metric_row(
                stage=STAGE_ORDER[3],
                metric="sg_current_abs_log10_error",
                carrier=carrier,
                support="endpoint_current_proxy_edges",
                unit="dex",
                values=current_errors[carrier],
            )
        )
        summaries.append(
            _metric_row(
                stage=STAGE_ORDER[3],
                metric="sg_current_sign_agreement",
                carrier=carrier,
                support="nonzero_endpoint_current_proxy_edges",
                unit="fraction",
                values=current_sign[carrier],
            )
        )
        summaries.append(
            _metric_row(
                stage=STAGE_ORDER[4],
                metric="alpha_abs_log10_error",
                carrier=carrier,
                support="endpoint_mean_node_alpha_edges",
                unit="dex",
                values=alpha_errors[carrier],
            )
        )
    for branch in ("baseline", "self_consistent"):
        summaries.append(
            _metric_row(
                stage=STAGE_ORDER[4],
                metric=f"{branch}_integrated_generation_abs_log10_error",
                carrier="combined",
                support="state",
                unit="dex",
                values=generation_errors[branch],
            )
        )

    node_path = output / "self_consistent_node_samples.csv"
    edge_path = output / "self_consistent_edge_samples.csv"
    generation_path = output / "self_consistent_generation_samples.csv"
    summary_path = output / "self_consistent_summary.csv"
    _write_csv(node_path, node_rows)
    _write_csv(edge_path, edge_rows)
    _write_csv(generation_path, generation_rows)
    _write_csv(summary_path, summaries)
    manifest: dict[str, object] = {
        "schema_version": 1,
        "status": "valid",
        "experiment": "minimal6_self_consistent_potential_replacement",
        "state_contract": {
            "state_count": len(states),
            "topologies": ["mirror", "sketch"],
            "biases_V": list(range(-20, 0)),
            "node_count": 6,
            "edge_count": 9,
            "cell_count": 4,
        },
        "stage_order": list(STAGE_ORDER),
        "replacement_fields": [
            "Sentaurus psi at all nodes",
            "Sentaurus electron QFP at all nodes",
            "Sentaurus hole QFP at all nodes",
        ],
        "derived_fields": [
            "n and p from Vela effective-ni/BGN and replaced potentials",
            "mobility from the Vela production model",
            "variable-ni SG current from the Vela production operator",
            "Van Overstraeten alpha and Genius-truncated avalanche source",
        ],
        "maximum_recovered_ni_gap_dex": maximum_ni_gap,
        "reference_limitations": [
            "Sentaurus edge current is endpoint-mean node current projected onto the edge",
            "Sentaurus edge alpha is an endpoint-mean node alpha",
            "Sentaurus generation anchor is P1 integration of native node ImpactIonization",
            "the experiment is fixed-state and does not run a nonlinear solve",
        ],
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
            "node_samples_csv": node_path.name,
            "node_samples_sha256": _sha256(node_path),
            "edge_samples_csv": edge_path.name,
            "edge_samples_sha256": _sha256(edge_path),
            "generation_samples_csv": generation_path.name,
            "generation_samples_sha256": _sha256(generation_path),
            "summary_csv": summary_path.name,
            "summary_sha256": _sha256(summary_path),
        },
    }
    report_path = output / "report.md"
    report_path.write_text(_markdown(manifest, summaries), encoding="utf-8")
    manifest["outputs"]["report_md"] = report_path.name
    manifest["outputs"]["report_sha256"] = _sha256(report_path)
    manifest_path = output / "manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest
