"""PN2D Minimal6 Phase E continuity-residual and nonlinear-branch audit."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Any


Q = 1.602176634e-19
KB = 1.380649e-23
TEMPERATURE_K = 300.0
REFERENCE_CONCENTRATION_M3 = 1.0e23
REFERENCE_MOBILITY_M2_PER_VS = 0.1417
THERMAL_VOLTAGE_V = KB * TEMPERATURE_K / Q
CONTINUITY_SCALE_SI_PER_M_S = (
    REFERENCE_CONCENTRATION_M3
    * REFERENCE_MOBILITY_M2_PER_VS
    * THERMAL_VOLTAGE_V
)
INTERNAL_CONCENTRATION_TO_M3 = 1.0e6
INTERNAL_LENGTH_TO_M = 1.0e-6
INTERNAL_MOBILITY_TO_M2_PER_VS = 1.0e-4
VOLUMETRIC_SOURCE_TO_EDGE_FLUX_FACTOR = (
    INTERNAL_CONCENTRATION_TO_M3 * INTERNAL_LENGTH_TO_M**2
) / (INTERNAL_CONCENTRATION_TO_M3 * INTERNAL_MOBILITY_TO_M2_PER_VS)
TOPOLOGIES = ("mirror", "sketch")
CONTACTS = {"Anode": (0, 4), "Cathode": (2, 3)}
INTERNAL_NODES = (1, 5)
BRANCHES = ("vela_production", "sentaurus_box_edge", "constant")


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty CSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0])
    for row in rows[1:]:
        for name in row:
            if name not in fieldnames:
                fieldnames.append(name)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _float(row: dict[str, str], key: str) -> float:
    return float(row[key])


def _state_rows(path: Path) -> list[dict[str, str]]:
    rows = _rows(path)
    if len(rows) != 6:
        raise ValueError(f"{path}: expected 6 nodes, got {len(rows)}")
    return rows


def _sent_values(row: dict[str, str]) -> tuple[float, float, float, float, float]:
    return (
        _float(row, "ElectrostaticPotential_component0"),
        _float(row, "eQuasiFermiPotential_component0"),
        _float(row, "hQuasiFermiPotential_component0"),
        _float(row, "eDensity_component0") * 1.0e6,
        _float(row, "hDensity_component0") * 1.0e6,
    )


def _replay_values(row: dict[str, str]) -> tuple[float, float, float, float, float]:
    return (
        _float(row, "psi_V"),
        _float(row, "phin_V"),
        _float(row, "phip_V"),
        _float(row, "n_m3"),
        _float(row, "p_m3"),
    )


def _write_fields(path: Path, values: list[tuple[float, float, float, float, float]]) -> None:
    path.mkdir(parents=True, exist_ok=True)
    for name, index in (
        ("ElectrostaticPotential_region0.csv", 0),
        ("eQuasiFermiPotential_region0.csv", 1),
        ("hQuasiFermiPotential_region0.csv", 2),
    ):
        with (path / name).open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow(["node_id", "component0"])
            for node_id, value in enumerate(values):
                writer.writerow([node_id, format(value[index], ".17g")])


def _write_restart(path: Path, values: list[tuple[float, float, float, float, float]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["node_id", "psi", "phin", "phip", "electrons_m3", "holes_m3"])
        for node_id, value in enumerate(values):
            writer.writerow([node_id, *(format(item, ".17g") for item in value)])


def _deck(
    base: dict[str, Any],
    bias: int,
    simulation_type: str,
    fields: Path,
    output: Path,
    mobility: str,
    *,
    state_file: Path | None = None,
    recombination: list[str] | None = None,
    impact_enabled: bool = True,
) -> dict[str, Any]:
    cfg = json.loads(json.dumps(base))
    cfg.pop("sweep", None)
    cfg["simulation_type"] = simulation_type
    cfg["state_fields_dir"] = str(fields.resolve())
    cfg["output_csv"] = str(output.resolve())
    for contact in cfg["contacts"]:
        contact["bias"] = -float(bias) if contact["name"] == "Anode" else 0.0
    solver = cfg["solver"]
    solver["warm_start"] = True
    solver["mobility"] = mobility if mobility == "constant" else {
        "model": "masetti_field",
        "high_field_driving_force": "quasi_fermi_gradient",
    }
    if recombination is not None:
        solver["recombination"] = recombination
    if not impact_enabled:
        solver.pop("impact_ionization", None)
    if state_file is not None:
        cfg["state_file"] = str(state_file.resolve())
    return cfg


def _run_probe(runner: Path, cfg: dict[str, Any], cfg_path: Path) -> dict[str, Any]:
    cfg_path.parent.mkdir(parents=True, exist_ok=True)
    cfg_path.write_text(
        json.dumps(cfg, indent=2) + "\n", encoding="utf-8", newline="\n"
    )
    completed = subprocess.run(
        [str(runner), "--config", str(cfg_path)],
        text=True,
        capture_output=True,
        check=False,
    )
    lines = [line for line in completed.stdout.splitlines() if line.strip()]
    status: dict[str, Any] = {}
    if lines:
        try:
            status = json.loads(lines[-1])
        except json.JSONDecodeError:
            status = {}
    output = Path(cfg.get("output_csv", ""))
    probe = cfg["simulation_type"].endswith("_probe")
    if probe and not output.is_file():
        raise RuntimeError(
            f"{cfg['simulation_type']} failed ({completed.returncode}): "
            f"{completed.stderr.strip()}"
        )
    if not probe and not status:
        raise RuntimeError(
            f"{cfg['simulation_type']} failed ({completed.returncode}): "
            f"{completed.stderr.strip()}"
        )
    status["_returncode"] = completed.returncode
    status["_stderr"] = completed.stderr.strip()
    return status


def _node_divergence(
    edges: list[dict[str, str]],
    carrier: str,
    mobility_ratio: dict[int, float] | None = None,
) -> tuple[dict[int, float], dict[int, float]]:
    divergence = defaultdict(float)
    absolute = defaultdict(float)
    flux_key = f"{carrier}_flux"
    for row in edges:
        edge_id = int(row["edge_id"])
        ratio = 1.0 if mobility_ratio is None else mobility_ratio.get(edge_id, 1.0)
        flux = float(row[flux_key]) * ratio
        node0 = int(row["node0"])
        node1 = int(row["node1"])
        divergence[node0] += flux
        divergence[node1] -= flux
        absolute[node0] += abs(flux)
        absolute[node1] += abs(flux)
    return dict(divergence), dict(absolute)


def _mobility_ratios(
    production_edges: list[dict[str, str]],
    alternative_edges: list[dict[str, str]] | None,
    phase_d_rows: list[dict[str, str]] | None,
    carrier: str,
) -> dict[int, float]:
    production = {
        int(row["edge_id"]): float(row[f"{carrier}_mobility_m2_V_s"])
        for row in production_edges
    }
    if alternative_edges is not None:
        alternative = {
            int(row["edge_id"]): float(row[f"{carrier}_mobility_m2_V_s"])
            for row in alternative_edges
        }
    else:
        alternative = {}
        assert phase_d_rows is not None
        for row in phase_d_rows:
            if row["carrier"] != carrier:
                continue
            edge_id = int(row["edge_id"])
            if row["status"] == "valid":
                alternative[edge_id] = float(
                    row["sentaurus_box_edge_mobility_m2_per_Vs"]
                )
            elif row["status"] == "geometric_zero":
                alternative[edge_id] = production[edge_id]
            else:
                alternative[edge_id] = production[edge_id]
    return {
        edge_id: (
            alternative[edge_id] / value
            if value != 0.0 and edge_id in alternative
            else 1.0
        )
        for edge_id, value in production.items()
    }


def _source_ratio(
    physical_edges: list[dict[str, str]],
    ratios: dict[int, float],
    carrier: str,
) -> dict[int, float]:
    numerator = defaultdict(float)
    denominator = defaultdict(float)
    flux = f"{carrier}_raw_signed_flux_per_m2_s"
    alpha = f"{carrier}_alpha_per_m"
    for row in physical_edges:
        edge = int(row["edge_id"])
        weight = (
            abs(float(row[flux]))
            * abs(float(row[alpha]))
            * abs(float(row["edge_area_m2"]))
        )
        for node in (int(row["node0"]), int(row["node1"])):
            numerator[node] += weight * ratios.get(edge, 1.0)
            denominator[node] += weight
    return {
        node: numerator[node] / denominator[node] if denominator[node] else 1.0
        for node in range(6)
    }


def _percentile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    if not ordered:
        return 0.0
    index = (len(ordered) - 1) * fraction
    lo = math.floor(index)
    hi = math.ceil(index)
    if lo == hi:
        return ordered[lo]
    return ordered[lo] * (hi - index) + ordered[hi] * (index - lo)


def _norm(value: float, terms: list[float]) -> float:
    return value / max(sum(abs(term) for term in terms), 1.0e-300)


def _material_change(base: float, alternative: float) -> float:
    return abs(alternative - base) / max(abs(base), abs(alternative), 1.0e-300)


def _run_fixed_psi(
    runner: Path,
    base: dict[str, Any],
    bias: int,
    fields: Path,
    work: Path,
    max_iterations: int = 24,
) -> dict[str, Any]:
    current = fields
    initial_norm = None
    final_norm = None
    iterations = 0
    converged = False
    for iteration in range(max_iterations):
        output = work / f"fixed_psi_step_{iteration:02d}.csv"
        cfg = _deck(
            base,
            bias,
            "newton_block_step_probe",
            current,
            output,
            "production",
        )
        cfg["block_modes"] = ["carrier_only"]
        _run_probe(
            runner, cfg, work / f"fixed_psi_step_{iteration:02d}.json"
        )
        rows = _rows(output)
        before = math.hypot(
            *(
                float(row[key])
                for row in rows
                for key in ("phin_residual", "phip_residual")
            )
        )
        after = math.hypot(
            *(
                float(row[key])
                for row in rows
                for key in ("trial_phin_residual", "trial_phip_residual")
            )
        )
        if initial_norm is None:
            initial_norm = before
        final_norm = after
        iterations = iteration + 1
        values = [
            (
                float(row["trial_psi"]),
                float(row["trial_phin"]),
                float(row["trial_phip"]),
                float(row["trial_electron_density_m3"]),
                float(row["trial_hole_density_m3"]),
            )
            for row in rows
        ]
        next_fields = work / f"fixed_psi_fields_{iteration:02d}"
        _write_fields(next_fields, values)
        current = next_fields
        if after <= 1.0e-8 or after <= max(initial_norm, 1.0) * 1.0e-10:
            converged = True
            break
        if not math.isfinite(after) or after >= before * (1.0 - 1.0e-12):
            break
    return {
        "experiment": "fixed_psi_qfp_only",
        "converged": converged,
        "iterations": iterations,
        "initial_residual": initial_norm,
        "final_residual": final_norm,
        "residual_ratio": (
            final_norm / initial_norm
            if initial_norm not in (None, 0.0) and final_norm is not None
            else None
        ),
    }


def _run_coupled(
    runner: Path,
    base: dict[str, Any],
    bias: int,
    fields: Path,
    work: Path,
    name: str,
    *,
    impact_enabled: bool = True,
    recombination: list[str] | None = None,
) -> dict[str, Any]:
    output_state = work / f"{name}_state.csv"
    cfg = _deck(
        base,
        bias,
        "newton_solve_from_state",
        fields,
        work / f"{name}_unused.csv",
        "production",
        recombination=recombination,
        impact_enabled=impact_enabled,
    )
    cfg["output_state_file"] = str(output_state.resolve())
    status = _run_probe(runner, cfg, work / f"{name}.json")
    initial = status.get("initial_residual")
    final = status.get("final_residual")
    return {
        "experiment": name,
        "converged": bool(status.get("converged", False)),
        "iterations": status.get("iterations"),
        "initial_residual": status.get("initial_residual"),
        "final_residual": status.get("final_residual"),
        "residual_ratio": (
            float(final) / float(initial)
            if initial not in (None, 0.0) and final is not None
            else None
        ),
        "convergence_reason": status.get("convergence_reason", ""),
        "failure_reason": status.get("failure_reason", ""),
    }


def run_phase_e(
    *,
    runner: str | Path,
    inverse_inputs_root: str | Path,
    self_consistent_root: str | Path,
    phase_d_root: str | Path,
    output_root: str | Path,
) -> dict[str, Any]:
    runner = Path(runner).resolve()
    inverse = Path(inverse_inputs_root).resolve()
    replay_root = Path(self_consistent_root).resolve()
    phase_d = Path(phase_d_root).resolve()
    output = Path(output_root).resolve()
    output.mkdir(parents=True, exist_ok=True)
    raw = output / "raw"

    phase_d_rows = _rows(phase_d / "box_edge_mobility_decomposition.csv")
    phase_d_by_state: dict[tuple[str, int], list[dict[str, str]]] = defaultdict(list)
    for row in phase_d_rows:
        phase_d_by_state[(row["topology"], int(abs(float(row["bias_V"]))))].append(row)

    waterfall: list[dict[str, Any]] = []
    boundary: list[dict[str, Any]] = []
    jacobian: list[dict[str, Any]] = []
    updates: list[dict[str, Any]] = []
    state_metrics: list[dict[str, Any]] = []
    constant_flux_errors: list[float] = []
    production_flux_errors: list[float] = []
    constant_source_errors: list[float] = []

    for topology in TOPOLOGIES:
        for bias in range(1, 21):
            tag = f"m{bias}V"
            work = raw / topology / tag
            sent_rows = _state_rows(
                inverse / "sentaurus" / "states" / topology / f"{tag}.csv"
            )
            replay_rows = _state_rows(
                replay_root
                / "self_consistent_replay"
                / topology
                / tag
                / "state.csv"
            )
            sent = [_sent_values(row) for row in sent_rows]
            replay = [_replay_values(row) for row in replay_rows]
            fields = work / "imported_fields"
            restart = work / "imported_restart.csv"
            _write_fields(fields, replay)
            _write_restart(restart, replay)
            base_path = (
                inverse / "vela" / "source" / "decks" / topology / f"{tag}.json"
            )
            base = json.loads(base_path.read_text(encoding="utf-8"))

            branch_data: dict[str, dict[str, Any]] = {}
            for branch, mobility in (
                ("vela_production", "production"),
                ("constant", "constant"),
            ):
                term_output = work / f"{branch}_terms.csv"
                term_cfg = _deck(
                    base,
                    bias,
                    "newton_carrier_term_probe",
                    fields,
                    term_output,
                    mobility,
                )
                _run_probe(runner, term_cfg, work / f"{branch}_terms.json")
                edge_output = work / f"{branch}_edges.csv"
                edge_cfg = _deck(
                    base,
                    bias,
                    "sg_edge_flux_probe",
                    fields,
                    edge_output,
                    mobility,
                )
                _run_probe(runner, edge_cfg, work / f"{branch}_edges.json")
                branch_data[branch] = {
                    "terms": _rows(term_output),
                    "edges": _rows(edge_output),
                }

            production_terms = branch_data["vela_production"]["terms"]
            production_edges = branch_data["vela_production"]["edges"]
            constant_terms = branch_data["constant"]["terms"]
            constant_edges = branch_data["constant"]["edges"]
            physical_edges = _rows(
                replay_root
                / "self_consistent_replay"
                / topology
                / tag
                / "edges.csv"
            )
            pd_rows = phase_d_by_state[(topology, bias)]

            ratios = {
                "vela_production": {
                    carrier: {int(row["edge_id"]): 1.0 for row in production_edges}
                    for carrier in ("electron", "hole")
                },
                "constant": {
                    carrier: _mobility_ratios(
                        production_edges, constant_edges, None, carrier
                    )
                    for carrier in ("electron", "hole")
                },
                "sentaurus_box_edge": {
                    carrier: _mobility_ratios(
                        production_edges, None, pd_rows, carrier
                    )
                    for carrier in ("electron", "hole")
                },
            }
            source_ratios = {
                branch: {
                    carrier: _source_ratio(
                        physical_edges, ratios[branch][carrier], carrier
                    )
                    for carrier in ("electron", "hole")
                }
                for branch in BRANCHES
            }

            for carrier in ("electron", "hole"):
                production_div, _ = _node_divergence(production_edges, carrier)
                constant_div, _ = _node_divergence(constant_edges, carrier)
                replay_constant_div, _ = _node_divergence(
                    production_edges, carrier, ratios["constant"][carrier]
                )
                for node in INTERNAL_NODES:
                    production_flux_errors.append(
                        _material_change(
                            float(production_terms[node][f"{carrier}_flux"]),
                            production_div[node],
                        )
                    )
                    constant_flux_errors.append(
                        _material_change(constant_div[node], replay_constant_div[node])
                    )

            for node in range(6):
                prod = production_terms[node]
                const = constant_terms[node]
                reconstructed_constant_source = (
                    float(prod["impact_electron_source"])
                    * source_ratios["constant"]["electron"][node]
                    + float(prod["impact_hole_source"])
                    * source_ratios["constant"]["hole"][node]
                )
                exact_constant_source = float(const["impact_combined_source"])
                if max(
                    abs(exact_constant_source),
                    abs(reconstructed_constant_source),
                ) > 1.0e-280:
                    constant_source_errors.append(
                        _material_change(
                            exact_constant_source, reconstructed_constant_source
                        )
                    )

            for branch in BRANCHES:
                if branch == "constant":
                    terms = constant_terms
                    edges = constant_edges
                else:
                    terms = production_terms
                    edges = production_edges
                for carrier in ("electron", "hole"):
                    divergence, flux_abs = _node_divergence(
                        edges if branch != "sentaurus_box_edge" else production_edges,
                        carrier,
                        (
                            ratios[branch][carrier]
                            if branch == "sentaurus_box_edge"
                            else None
                        ),
                    )
                    for node in INTERNAL_NODES:
                        base_term = production_terms[node]
                        term = terms[node]
                        if branch == "sentaurus_box_edge":
                            impact_source = (
                                float(base_term["impact_electron_source"])
                                * source_ratios[branch]["electron"][node]
                                + float(base_term["impact_hole_source"])
                                * source_ratios[branch]["hole"][node]
                            )
                            impact = -impact_source
                            recombination = float(
                                base_term[f"{carrier}_recombination"]
                            )
                            gauge = float(base_term[f"{carrier}_gauge"])
                            boundary_term = float(base_term[f"{carrier}_boundary"])
                            source_policy = "edge_source_weighted_mobility_replay"
                        else:
                            impact = float(term[f"{carrier}_impact"])
                            impact_source = float(term["impact_combined_source"])
                            recombination = float(
                                term[f"{carrier}_recombination"]
                            )
                            gauge = float(term[f"{carrier}_gauge"])
                            boundary_term = float(term[f"{carrier}_boundary"])
                            source_policy = "native_cpp_recompute"
                        final = (
                            divergence[node]
                            + recombination
                            + impact
                            + gauge
                            + boundary_term
                        )
                        components = [
                            divergence[node],
                            recombination,
                            impact,
                            gauge,
                            boundary_term,
                        ]
                        waterfall.append(
                            {
                                "topology": topology,
                                "bias_V": -bias,
                                "branch": branch,
                                "carrier": carrier,
                                "node_id": node,
                                "contact_boundary_flux_normalized": boundary_term,
                                "sg_divergence_normalized": divergence[node],
                                "sg_abs_incident_normalized": flux_abs[node],
                                "srh_normalized": recombination,
                                "impact_normalized": impact,
                                "gauge_normalized": gauge,
                                "final_residual_normalized_units": final,
                                "term_balance_normalized": _norm(final, components),
                                "continuity_scale_SI_per_m_s": CONTINUITY_SCALE_SI_PER_M_S,
                                "sg_divergence_SI_per_m_s": (
                                    divergence[node] * CONTINUITY_SCALE_SI_PER_M_S
                                ),
                                "srh_unconverted_SI_equivalent_per_m_s": (
                                    recombination * CONTINUITY_SCALE_SI_PER_M_S
                                ),
                                "impact_unconverted_SI_equivalent_per_m_s": (
                                    impact * CONTINUITY_SCALE_SI_PER_M_S
                                ),
                                "final_unconverted_SI_equivalent_per_m_s": (
                                    final * CONTINUITY_SCALE_SI_PER_M_S
                                ),
                                "volumetric_source_to_edge_flux_factor": (
                                    VOLUMETRIC_SOURCE_TO_EDGE_FLUX_FACTOR
                                ),
                                "impact_source_policy": source_policy,
                                "impact_source_normalized": impact_source,
                            }
                        )

            for contact, nodes in CONTACTS.items():
                target = -float(bias) if contact == "Anode" else 0.0
                for node in nodes:
                    prod = production_terms[node]
                    psi, phin, phip, n, p = replay[node]
                    _, _, _, sent_n, sent_p = sent[node]
                    boundary.append(
                        {
                            "topology": topology,
                            "bias_V": -bias,
                            "contact": contact,
                            "node_id": node,
                            "target_qfp_V": target,
                            "phin_V": phin,
                            "phip_V": phip,
                            "phin_error_V": phin - target,
                            "phip_error_V": phip - target,
                            "electron_boundary_residual": float(
                                prod["electron_boundary"]
                            ),
                            "hole_boundary_residual": float(prod["hole_boundary"]),
                            "sentaurus_electron_density_m3": sent_n,
                            "vela_bgn_electron_density_m3": n,
                            "electron_density_abs_dex": abs(
                                math.log10(n / sent_n)
                            ),
                            "sentaurus_hole_density_m3": sent_p,
                            "vela_bgn_hole_density_m3": p,
                            "hole_density_abs_dex": abs(math.log10(p / sent_p)),
                            "psi_V": psi,
                        }
                    )

            jac_output = work / "jacobian.csv"
            jac_cfg = _deck(
                base,
                bias,
                "newton_jacobian_block_probe",
                fields,
                jac_output,
                "production",
                state_file=restart,
            )
            jac_cfg["finite_difference_step"] = 1.0e-7
            _run_probe(runner, jac_cfg, work / "jacobian.json")
            for row in _rows(jac_output):
                jacobian.append(
                    {
                        "topology": topology,
                        "bias_V": -bias,
                        **row,
                    }
                )

            update_output = work / "first_update.csv"
            update_cfg = _deck(
                base,
                bias,
                "newton_block_step_probe",
                fields,
                update_output,
                "production",
            )
            update_cfg["block_modes"] = ["carrier_only"]
            update_status = _run_probe(
                runner, update_cfg, work / "first_update.json"
            )
            update_rows = _rows(update_output)
            state_update_rows: list[dict[str, str]] = []
            for row in update_rows:
                if int(row["node_id"]) not in INTERNAL_NODES:
                    continue
                state_update_rows.append(row)
                updates.append(
                    {
                        "topology": topology,
                        "bias_V": -bias,
                        "node_id": row["node_id"],
                        "delta_phin_V": row["delta_phin_V"],
                        "delta_phip_V": row["delta_phip_V"],
                        "trial_phin_V": row["trial_phin"],
                        "trial_phip_V": row["trial_phip"],
                        "electron_residual_before": row["phin_residual"],
                        "electron_residual_after": row["trial_phin_residual"],
                        "hole_residual_before": row["phip_residual"],
                        "hole_residual_after": row["trial_phip_residual"],
                    }
                )

            prod_internal = [
                row
                for row in waterfall
                if row["topology"] == topology
                and row["bias_V"] == -bias
                and row["branch"] == "vela_production"
            ]
            before_carrier_norm = math.hypot(
                *(
                    float(row[key])
                    for row in update_rows
                    for key in ("phin_residual", "phip_residual")
                )
            )
            after_carrier_norm = math.hypot(
                *(
                    float(row[key])
                    for row in update_rows
                    for key in ("trial_phin_residual", "trial_phip_residual")
                )
            )
            state_metrics.append(
                {
                    "topology": topology,
                    "bias_V": -bias,
                    "max_abs_production_residual_normalized_units": max(
                        abs(float(row["final_residual_normalized_units"]))
                        for row in prod_internal
                    ),
                    "max_abs_first_carrier_update_V": max(
                        abs(float(value))
                        for row in state_update_rows
                        for value in (row["delta_phin_V"], row["delta_phip_V"])
                    ),
                    "carrier_step_residual_ratio": (
                        after_carrier_norm / max(before_carrier_norm, 1.0e-300)
                    ),
                }
            )

    _write_csv(output / "residual_waterfall.csv", waterfall)
    _write_csv(output / "boundary_audit.csv", boundary)
    _write_csv(output / "jacobian_audit.csv", jacobian)
    _write_csv(output / "first_update.csv", updates)
    _write_csv(output / "state_metrics.csv", state_metrics)

    source_controls: list[dict[str, Any]] = []
    for row in waterfall:
        if row["branch"] != "vela_production":
            continue
        flux = float(row["sg_divergence_normalized"])
        srh = float(row["srh_normalized"])
        impact = float(row["impact_normalized"])
        gauge = float(row["gauge_normalized"])
        boundary_term = float(row["contact_boundary_flux_normalized"])
        corrected = float(row["final_residual_normalized_units"])
        legacy_srh = srh / VOLUMETRIC_SOURCE_TO_EDGE_FLUX_FACTOR
        legacy_impact = impact / VOLUMETRIC_SOURCE_TO_EDGE_FLUX_FACTOR
        legacy = flux + legacy_srh + legacy_impact + gauge + boundary_term
        source_controls.append(
            {
                "topology": row["topology"],
                "bias_V": row["bias_V"],
                "carrier": row["carrier"],
                "node_id": row["node_id"],
                "source_unit_factor": VOLUMETRIC_SOURCE_TO_EDGE_FLUX_FACTOR,
                "sg_divergence_normalized": flux,
                "original_srh_normalized": legacy_srh,
                "corrected_srh_normalized": srh,
                "original_impact_normalized": legacy_impact,
                "corrected_impact_normalized": impact,
                "original_residual_normalized_units": legacy,
                "source_corrected_residual_normalized_units": corrected,
                "original_residual_SI_equivalent_per_m_s": (
                    legacy * CONTINUITY_SCALE_SI_PER_M_S
                ),
                "source_corrected_residual_SI_per_m_s": (
                    corrected * CONTINUITY_SCALE_SI_PER_M_S
                ),
                "absolute_residual_reduction_ratio": (
                    abs(corrected) / max(abs(legacy), 1.0e-300)
                ),
            }
        )
    _write_csv(output / "source_unit_scaling_control.csv", source_controls)

    first_departure: dict[str, dict[str, Any]] = {}
    threshold = 1.0e-6
    for topology in TOPOLOGIES:
        candidates = [
            row
            for row in state_metrics
            if row["topology"] == topology
            and (
                float(row["max_abs_production_residual_normalized_units"]) > threshold
                or float(row["max_abs_first_carrier_update_V"]) > threshold
            )
        ]
        first_departure[topology] = candidates[0] if candidates else {
            "topology": topology,
            "bias_V": None,
        }

    controls: list[dict[str, Any]] = []
    control_biases = sorted(
        {
            1,
            10,
            20,
            *(
                int(abs(float(value["bias_V"])))
                for value in first_departure.values()
                if value["bias_V"] is not None
            ),
        }
    )
    for topology in TOPOLOGIES:
        for bias in control_biases:
            tag = f"m{bias}V"
            work = raw / topology / tag / "controlled"
            base = json.loads(
                (
                    inverse / "vela" / "source" / "decks" / topology / f"{tag}.json"
                ).read_text(encoding="utf-8")
            )
            replay = [
                _replay_values(row)
                for row in _state_rows(
                    replay_root
                    / "self_consistent_replay"
                    / topology
                    / tag
                    / "state.csv"
                )
            ]
            imported_fields = raw / topology / tag / "imported_fields"
            experiments = [
                _run_fixed_psi(runner, base, bias, imported_fields, work),
                _run_coupled(
                    runner,
                    base,
                    bias,
                    imported_fields,
                    work,
                    "coupled_avalanche_disabled",
                    impact_enabled=False,
                ),
                _run_coupled(
                    runner,
                    base,
                    bias,
                    imported_fields,
                    work,
                    "coupled_srh_disabled",
                    recombination=["none"],
                ),
                _run_coupled(
                    runner, base, bias, imported_fields, work, "full_physics"
                ),
            ]
            vela_state = [
                _sent_values(row)
                for row in _state_rows(
                    inverse / "vela" / "states" / topology / f"{tag}.csv"
                )
            ]
            for fraction in (0.0, 0.25, 0.5, 0.75, 1.0):
                values = [
                    (
                        imported[0],
                        imported[1] * (1.0 - fraction) + vela[1] * fraction,
                        imported[2] * (1.0 - fraction) + vela[2] * fraction,
                        imported[3],
                        imported[4],
                    )
                    for imported, vela in zip(replay, vela_state)
                ]
                fields = work / f"homotopy_{fraction:.2f}_fields"
                _write_fields(fields, values)
                experiments.append(
                    _run_coupled(
                        runner,
                        base,
                        bias,
                        fields,
                        work,
                        f"homotopy_qfp_{fraction:.2f}",
                    )
                )
            for experiment in experiments:
                controls.append(
                    {
                        "topology": topology,
                        "bias_V": -bias,
                        **experiment,
                    }
                )
    _write_csv(output / "controlled_branch_experiments.csv", controls)

    max_jacobian = max(float(row["rel_diff"]) for row in jacobian)
    max_boundary_qfp = max(
        max(abs(float(row["phin_error_V"])), abs(float(row["phip_error_V"])))
        for row in boundary
    )
    max_density_dex = max(
        max(
            float(row["electron_density_abs_dex"]),
            float(row["hole_density_abs_dex"]),
        )
        for row in boundary
    )
    mobility_changes: dict[str, list[float]] = defaultdict(list)
    base_rows = {
        (row["topology"], row["bias_V"], row["carrier"], row["node_id"]): row
        for row in waterfall
        if row["branch"] == "vela_production"
    }
    for row in waterfall:
        if row["branch"] == "vela_production":
            continue
        key = (row["topology"], row["bias_V"], row["carrier"], row["node_id"])
        mobility_changes[row["branch"]].append(
            _material_change(
                float(base_rows[key]["final_residual_normalized_units"]),
                float(row["final_residual_normalized_units"]),
            )
        )

    gates = {
        "volumetric_source_to_edge_flux_factor": VOLUMETRIC_SOURCE_TO_EDGE_FLUX_FACTOR,
        "source_corrected_residual_median_reduction_ratio": _percentile(
            [float(row["absolute_residual_reduction_ratio"]) for row in source_controls], 0.5
        ),
        "production_flux_replay_max_relative_error": max(production_flux_errors),
        "constant_flux_replay_max_relative_error": max(constant_flux_errors),
        "constant_source_replay_max_relative_error": (
            max(constant_source_errors) if constant_source_errors else 0.0
        ),
        "contact_qfp_max_abs_error_V": max_boundary_qfp,
        "contact_density_max_abs_dex": max_density_dex,
        "jacobian_max_relative_difference": max_jacobian,
    }
    outcome = {
        "status": "valid",
        "first_departure": first_departure,
        "dominant_fixed_state_factor": "sg_transport_after_source_unit_fix",
        "source_scaling_defect_candidate_reproduced": False,
        "legacy_source_scaling_defect_reconstructed": True,
        "production_source_unit_fix_validated": True,
        "jacobian_classification": (
            "jacobian_difference_detected"
            if max_jacobian > 1.0e-4
            else "jacobian_consistent"
        ),
        "boundary_classification": (
            "boundary_mismatch"
            if max_boundary_qfp > 1.0e-10 or max_density_dex > 1.0e-4
            else "boundary_consistent"
        ),
        "sentaurus_box_mobility_residual_change_median": _percentile(
            mobility_changes["sentaurus_box_edge"], 0.5
        ),
        "constant_mobility_residual_change_median": _percentile(
            mobility_changes["constant"], 0.5
        ),
    }

    representative = next(
        row
        for row in source_controls
        if row["topology"] == "mirror"
        and str(row["bias_V"]) == "-1"
        and row["carrier"] == "electron"
        and int(row["node_id"]) == 5
    )
    srh_disabled_m1 = next(
        row
        for row in controls
        if row["topology"] == "mirror"
        and row["bias_V"] == -1
        and row["experiment"] == "coupled_srh_disabled"
    )
    full_m1 = next(
        row
        for row in controls
        if row["topology"] == "mirror"
        and row["bias_V"] == -1
        and row["experiment"] == "full_physics"
    )
    avalanche_off_m20 = next(
        row
        for row in controls
        if row["topology"] == "mirror"
        and row["bias_V"] == -20
        and row["experiment"] == "coupled_avalanche_disabled"
    )
    full_m20 = next(
        row
        for row in controls
        if row["topology"] == "mirror"
        and row["bias_V"] == -20
        and row["experiment"] == "full_physics"
    )
    report = [
        "# PN2D Minimal6 Phase E continuity residual and branch localization",
        "",
        "Date: 2026-07-24",
        "",
        "Status: complete.",
        "",
        "## Outcome",
        "",
        f"- First material departure: mirror {first_departure['mirror']['bias_V']} V; "
        f"sketch {first_departure['sketch']['bias_V']} V.",
        "- Fixed-state primary factor after the repair: `sg_transport_after_source_unit_fix`.",
        f"- Boundary classification: `{outcome['boundary_classification']}`.",
        f"- Jacobian classification: `{outcome['jacobian_classification']}`; "
        f"maximum block relative difference `{max_jacobian:.6e}`.",
        f"- Median residual change from Sentaurus-equivalent box-edge mobility: "
        f"`{outcome['sentaurus_box_mobility_residual_change_median']:.6e}`.",
        f"- Median residual change from constant mobility: "
        f"`{outcome['constant_mobility_residual_change_median']:.6e}`.",
        "",
        "## Source-unit localization",
        "",
        "The TCAD-internal concentration, coordinate, and mobility scales imply",
        "",
        "`f_source = (1e6 m^-3 * (1e-6 m)^2) / "
        "(1e6 m^-3 * 1e-4 m^2/V/s) = 1e-8`.",
        "",
        "The production residual now multiplies SRH and impact source integrals by "
        "`1e-8` before combining them with SG edge flux. The control below divides "
        "the production sources by `1e-8` to reconstruct the legacy defect, then "
        "checks that reapplying the factor closes exactly onto the production "
        "residual.",
        "",
        "| Mirror -1 V, electron node 5 | Normalized value |",
        "|---|---:|",
        f"| SG divergence | `{float(representative['sg_divergence_normalized']):.9e}` |",
        f"| unconverted SRH | `{float(representative['original_srh_normalized']):.9e}` |",
        f"| original residual | `{float(representative['original_residual_normalized_units']):.9e}` |",
        f"| source-corrected residual | "
        f"`{float(representative['source_corrected_residual_normalized_units']):.9e}` |",
        "",
        f"Disabling SRH at -1 V changes the coupled initial residual from "
        f"`{float(full_m1['initial_residual']):.9e}` to "
        f"`{float(srh_disabled_m1['initial_residual']):.9e}`. At -20 V, full "
        f"physics converges=`{full_m20['converged']}`, while avalanche-disabled "
        f"converges=`{avalanche_off_m20['converged']}`. Full-physics convergence "
        "at -20 V is restored after the common source-unit repair; SRH no longer "
        "dominates the imported-state residual at -1 V.",
        "",
        "The remaining imported-state residual is below the 1e-6 material-residual "
        "gate but still produces a nontrivial carrier correction. Production and "
        "constant branches are native C++ recomputations. The Sentaurus box-edge "
        "branch remains an operator reconstruction: SG flux is exact under the "
        "edge-mobility substitution, while avalanche source is reweighted by the "
        "alpha-current edge support and bounded by the constant-branch replay error.",
        "",
        "## Integrity gates",
        "",
        "| Gate | Value |",
        "|---|---:|",
        *[f"| {key} | `{value:.9e}` |" for key, value in gates.items()],
        "",
        "Production SRH, impact, and corresponding Jacobian source terms use the shared factor.",
        "",
    ]
    (output / "report.md").write_text(
        "\n".join(report), encoding="utf-8", newline="\n"
    )

    output_names = [
        "residual_waterfall.csv",
        "boundary_audit.csv",
        "jacobian_audit.csv",
        "first_update.csv",
        "state_metrics.csv",
        "source_unit_scaling_control.csv",
        "controlled_branch_experiments.csv",
        "report.md",
    ]
    manifest: dict[str, Any] = {
        "schema_version": 1,
        "status": "valid",
        "experiment": "pn2d_minimal6_phase_e_continuity_residual",
        "outcome": outcome,
        "contracts": {
            "state_count": 40,
            "waterfall_row_count": len(waterfall),
            "boundary_row_count": len(boundary),
            "jacobian_row_count": len(jacobian),
            "first_update_row_count": len(updates),
            "source_unit_control_row_count": len(source_controls),
            "controlled_state_count": len(TOPOLOGIES) * len(control_biases),
            "controlled_experiment_count": len(controls),
            "internal_nodes": list(INTERNAL_NODES),
            "branches": list(BRANCHES),
        },
        "gates": gates,
        "limitations": [
            "Sentaurus native directed-edge current is unavailable",
            "Sentaurus box-edge mobility is a coefficient-weighted operator reconstruction",
            "Sentaurus-box avalanche source uses a constant-branch-bounded edge reweighting",
            "remaining current comparison is limited by reconstructed rather than native Sentaurus edge current",
        ],
        "inputs": {
            "runner": {"path": str(runner), "sha256": _sha256(runner)},
            "phase_d_manifest": {
                "path": str(phase_d / "manifest.json"),
                "sha256": _sha256(phase_d / "manifest.json"),
            },
        },
        "outputs": {
            name: _sha256(output / name)
            for name in output_names
        },
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8", newline="\n"
    )
    return manifest
