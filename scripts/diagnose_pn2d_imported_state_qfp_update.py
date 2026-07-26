#!/usr/bin/env python3
"""Audit production-triangle versus opt-in element-edge first QFP updates."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
from pathlib import Path
from typing import Any

try:
    from diagnose_pn2d_general_tri3_element_edge_avalanche import parse_log
except ModuleNotFoundError:
    from scripts.diagnose_pn2d_general_tri3_element_edge_avalanche import parse_log


BIASES = (1, 10, 20)
CONTINUITY_SCALE = 1.0e23 * 0.1417 * (1.380649e-23 * 300.0 / 1.602176634e-19)


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, values: list[dict[str, Any]]) -> None:
    if not values:
        raise ValueError(f"refusing to write empty CSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    names = list(values[0])
    for value in values[1:]:
        for name in value:
            if name not in names:
                names.append(name)
    with path.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=names)
        writer.writeheader()
        writer.writerows(values)


def write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n", encoding="ascii", newline="\n")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run_probe(runner: Path, cfg: dict[str, Any], cfg_path: Path) -> dict[str, Any]:
    write_json(cfg_path, cfg)
    completed = subprocess.run(
        [str(runner), "--config", str(cfg_path)], text=True,
        capture_output=True, check=False,
    )
    lines = [line for line in completed.stdout.splitlines() if line.strip()]
    status: dict[str, Any] = {}
    if lines:
        try:
            status = json.loads(lines[-1])
        except json.JSONDecodeError:
            pass
    output = Path(cfg["output_csv"])
    if completed.returncode != 0 or not output.is_file():
        raise RuntimeError(
            f"{cfg['simulation_type']} failed ({completed.returncode}): "
            f"{completed.stderr.strip()}"
        )
    status["returncode"] = completed.returncode
    status["stderr"] = completed.stderr.strip()
    return status


def write_fields(path: Path, values: dict[int, tuple[float, float, float]]) -> None:
    path.mkdir(parents=True, exist_ok=True)
    for filename, index in (
        ("ElectrostaticPotential_region0.csv", 0),
        ("eQuasiFermiPotential_region0.csv", 1),
        ("hQuasiFermiPotential_region0.csv", 2),
    ):
        with (path / filename).open("w", newline="", encoding="ascii") as stream:
            writer = csv.writer(stream)
            writer.writerow(("node_id", "component0"))
            for node_id in sorted(values):
                writer.writerow((node_id, format(values[node_id][index], ".17g")))


def write_restart(path: Path, fields: Path, term_rows: list[dict[str, str]]) -> None:
    values: dict[int, list[float]] = {}
    names = (
        "ElectrostaticPotential_region0.csv",
        "eQuasiFermiPotential_region0.csv",
        "hQuasiFermiPotential_region0.csv",
    )
    for index, name in enumerate(names):
        for row in rows(fields / name):
            values.setdefault(int(row["node_id"]), [0.0, 0.0, 0.0])[index] = float(row["component0"])
    with path.open("w", newline="", encoding="ascii") as stream:
        writer = csv.writer(stream)
        writer.writerow(("node_id", "psi", "phin", "phip", "electrons_m3", "holes_m3"))
        for row in term_rows:
            node = int(row["node_id"])
            # The block audit packs psi/QFP only; finite positive density fields keep
            # the restart parser schema explicit and deterministic.
            writer.writerow((node, *(format(v, ".17g") for v in values[node]), "1", "1"))


def impact_config(variant: str) -> dict[str, Any]:
    common: dict[str, Any] = {
        "model": "van_overstraeten",
        "driving_force": "quasi_fermi_gradient",
        "generation": "current_density",
        "quasi_fermi_gradient_discretization": "cell_gradient",
    }
    if variant == "production_triangle":
        return {
            **common,
            "current_approximation": "cell_reconstructed",
            "current_magnitude_mode": "edge_scalar_abs",
            "cell_reconstructed_midpoint_density": "gss_logistic",
            "source_volume_policy": "genius_truncated",
            "source_volume_factor": 0.0,
            "source_geometry_scale": 1.0,
            "edge_source_partition": "symmetric",
            "driving_force_interpolation": "none",
            "quasi_fermi_carrier_truncation": 0.0,
            "minimum_field_V_m": 0.0,
            "electron_driving_force_ref_density_m3": 0.0,
            "hole_driving_force_ref_density_m3": 0.0,
            "source_mapping_mode": "triangle_gss_gradqf_truncated",
        }
    if variant == "element_edge_opt_in":
        return {
            **common,
            "current_approximation": "element_edge_sg_gss_laux",
            "source_mapping_mode": "element_vertex_box_measure",
        }
    raise ValueError(variant)


def deck(base: dict[str, Any], bias: int, fields: Path, output: Path,
         simulation_type: str, variant: str, control: str = "full") -> dict[str, Any]:
    cfg = json.loads(json.dumps(base))
    cfg.pop("sweep", None)
    cfg["simulation_type"] = simulation_type
    cfg["state_fields_dir"] = str(fields.resolve())
    cfg["output_csv"] = str(output.resolve())
    for contact in cfg["contacts"]:
        contact["bias"] = -float(bias) if contact["name"] == "Anode" else 0.0
    solver = cfg.setdefault("solver", {})
    solver.update({"method": "gummel_newton", "warm_start": True,
                   "reltol": 1.0e-8, "abstol": 1.0e-9})
    solver["mobility"] = {
        "model": "masetti_field",
        "high_field_driving_force": "quasi_fermi_gradient",
    }
    solver["recombination"] = [] if control == "srh_off" else ["srh"]
    if control == "avalanche_off":
        solver.pop("impact_ionization", None)
    else:
        solver["impact_ionization"] = impact_config(variant)
    return cfg


def coarse_inputs(log: Path, mesh_root: Path, output: Path, materials_file: Path) -> tuple[dict[str, Any], dict[int, Path], list[dict[str, Any]]]:
    groups = parse_log(log)
    mesh_path = mesh_root / "mesh.json"
    mesh = json.loads(mesh_path.read_text(encoding="utf-8"))
    mesh_nodes = {int(row["id"]): (float(row["x"]), float(row["y"])) for row in mesh["nodes"]}
    topology: list[dict[str, Any]] = []
    fields_by_bias: dict[int, Path] = {}
    for bias in BIASES:
        selected = {
            int(row["vertex"]): row for row in groups["vertices"]
            if float(row["bias_V"]) == -float(bias) and int(row["vertex"]) in mesh_nodes
        }
        if set(selected) != set(mesh_nodes):
            raise ValueError(f"coarse m{bias}V node mismatch")
        maximum = max(
            max(abs(float(selected[node][key]) - mesh_nodes[node][axis]) for key, axis in (("x_um", 0), ("y_um", 1)))
            for node in mesh_nodes
        )
        values = {
            node: (float(row["psi_V"]), float(row["eQFP_V"]), float(row["hQFP_V"]))
            for node, row in selected.items()
        }
        fields = output / "coarse_fields" / f"m{bias}V"
        write_fields(fields, values)
        fields_by_bias[bias] = fields
        topology.append({"topology": "coarse7x3", "bias_V": -bias,
                         "nodes": len(values), "max_coordinate_error_um": maximum})
    measures = [row for row in groups["measures"] if float(row["bias_V"]) == -1.0]
    permutations: dict[int, list[tuple[int, int]]] = {}
    for row in measures:
        permutations.setdefault(int(row["element"]), []).append((int(row["local_vertex"]), int(row["vertex"])))
    probe_cells = {key: [node for _, node in sorted(value)] for key, value in permutations.items()}
    mesh_cells = {int(row["id"]): [int(node) for node in row["node_ids"]] for row in mesh["triangles"]}
    if probe_cells != mesh_cells:
        raise ValueError("coarse element/local-vertex topology mismatch")
    base = {
        "simulation_type": "newton_carrier_term_probe",
        "mesh_file": str(mesh_path.resolve()),
        "materials_file": str(materials_file.resolve()),
        "node_doping_file": str((mesh_root / "doping.csv").resolve()),
        "scaling": {"mode": "unit_scaling"},
        "contacts": [{"name": "Cathode", "bias": 0.0}, {"name": "Anode", "bias": 0.0}],
        "solver": {"bandgap_narrowing": "old_slotboom", "contact_boundary_minority_electron_relaxation": False},
    }
    return base, fields_by_bias, topology


def internal_nodes(base: dict[str, Any]) -> set[int]:
    mesh = json.loads(Path(base["mesh_file"]).read_text(encoding="utf-8"))
    contacts = {int(node) for contact in mesh["contacts"] for node in contact["node_ids"]}
    return {int(row["id"]) for row in mesh["nodes"]} - contacts


def execute_state(runner: Path, topology: str, bias: int, base: dict[str, Any],
                  fields: Path, work: Path, residuals: list[dict[str, Any]],
                  updates: list[dict[str, Any]], jacobians: list[dict[str, Any]],
                  controls: list[dict[str, Any]]) -> None:
    internal = internal_nodes(base)
    edge_output = work / "native_mobility_edges.csv"
    edge_cfg = deck(base, bias, fields, edge_output, "sg_edge_flux_probe", "production_triangle")
    run_probe(runner, edge_cfg, work / "native_mobility_edges.json")
    edge_rows = rows(edge_output)
    incident = {carrier: {node: 0.0 for node in internal} for carrier in ("electron", "hole")}
    for edge in edge_rows:
        for carrier in incident:
            value = abs(float(edge[f"{carrier}_flux"]))
            for node in (int(edge["node0"]), int(edge["node1"])):
                if node in internal:
                    incident[carrier][node] += value

    for variant in ("production_triangle", "element_edge_opt_in"):
        term_output = work / f"{variant}_terms.csv"
        term_cfg = deck(base, bias, fields, term_output, "newton_carrier_term_probe", variant)
        run_probe(runner, term_cfg, work / f"{variant}_terms.json")
        term_rows = rows(term_output)
        for row in term_rows:
            node = int(row["node_id"])
            for carrier in ("electron", "hole"):
                terms = [float(row[f"{carrier}_{name}"]) for name in ("flux", "recombination", "impact", "gauge", "boundary")]
                final = sum(terms)
                assembler = float(row[f"{carrier}_residual"])
                scale = max(incident[carrier].get(node, 0.0), abs(terms[1]), abs(terms[2]))
                residuals.append({
                    "topology": topology, "bias_V": -bias, "variant": variant,
                    "carrier": carrier, "node_id": node, "is_boundary": int(node not in internal),
                    "sg_divergence_normalized": terms[0], "srh_normalized": terms[1],
                    "avalanche_normalized": terms[2], "gauge_normalized": terms[3],
                    "boundary_normalized": terms[4], "row_scale_normalized": scale,
                    "final_residual_normalized": final,
                    "assembled_residual_normalized": assembler,
                    "closure_relative": abs(final - assembler) / max(abs(final), abs(assembler), 1.0e-300),
                    "final_residual_physical_per_m_s": final * CONTINUITY_SCALE,
                    "row_scale_physical_per_m_s": scale * CONTINUITY_SCALE,
                })

        update_specs = (
            ("carrier_only", "newton_block_step_probe"),
            ("coupled", "newton_step_probe"),
        )
        for mode, simulation_type in update_specs:
            update_output = work / f"{variant}_{mode}_first_update.csv"
            update_cfg = deck(base, bias, fields, update_output, simulation_type, variant)
            if mode == "carrier_only":
                update_cfg["block_modes"] = ["carrier_only"]
            run_probe(runner, update_cfg, work / f"{variant}_{mode}_first_update.json")
            for row in rows(update_output):
                node = int(row["node_id"])
                if node not in internal:
                    continue
                for carrier, qfp, delta, before, after in (
                    ("electron", "phin", "delta_phin_V", "phin_residual", "trial_phin_residual"),
                    ("hole", "phip", "delta_phip_V", "phip_residual", "trial_phip_residual"),
                ):
                    updates.append({
                        "topology": topology, "bias_V": -bias, "variant": variant,
                        "mode": mode, "carrier": carrier, "node_id": node,
                        "initial_qfp_V": row[qfp], "delta_qfp_V": row[delta],
                        "sentaurus_error_before_V": 0.0,
                        "sentaurus_error_after_V": abs(float(row[delta])),
                        "residual_before": row[before], "residual_after": row[after],
                    })

        restart = work / f"{variant}_restart.csv"
        write_restart(restart, fields, term_rows)
        jac_output = work / f"{variant}_jacobian.csv"
        jac_cfg = deck(base, bias, fields, jac_output, "newton_jacobian_block_probe", variant)
        jac_cfg["state_file"] = str(restart.resolve())
        jac_cfg["finite_difference_step"] = 1.0e-7
        run_probe(runner, jac_cfg, work / f"{variant}_jacobian.json")
        for row in rows(jac_output):
            jacobians.append({"topology": topology, "bias_V": -bias, "variant": variant, **row})

    for control in ("avalanche_off", "srh_off"):
        output = work / f"element_edge_opt_in_{control}_terms.csv"
        cfg = deck(base, bias, fields, output, "newton_carrier_term_probe", "element_edge_opt_in", control)
        run_probe(runner, cfg, work / f"element_edge_opt_in_{control}_terms.json")
        for row in rows(output):
            for carrier in ("electron", "hole"):
                controls.append({
                    "topology": topology, "bias_V": -bias, "control": control,
                    "carrier": carrier, "node_id": row["node_id"],
                    "flux": row[f"{carrier}_flux"], "srh": row[f"{carrier}_recombination"],
                    "avalanche": row[f"{carrier}_impact"], "boundary": row[f"{carrier}_boundary"],
                    "residual": row[f"{carrier}_residual"],
                })


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--minimal-phase-e-root", type=Path, required=True)
    parser.add_argument("--coarse-log", type=Path, required=True)
    parser.add_argument("--coarse-mesh-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()
    output = args.output_root.resolve()
    output.mkdir(parents=True, exist_ok=True)
    minimal_seed = json.loads((args.minimal_phase_e_root.resolve() / "raw" / "mirror" / "m1V" / "vela_production_terms.json").read_text(encoding="utf-8"))
    coarse_base, coarse_fields, topology_gate = coarse_inputs(args.coarse_log.resolve(), args.coarse_mesh_root.resolve(), output, Path(minimal_seed["materials_file"]))
    cases: list[tuple[str, int, dict[str, Any], Path]] = []
    for topology in ("minimal6_mirror", "minimal6_sketch"):
        source_name = topology.removeprefix("minimal6_")
        for bias in BIASES:
            source = args.minimal_phase_e_root.resolve() / "raw" / source_name / f"m{bias}V"
            base = json.loads((source / "vela_production_terms.json").read_text(encoding="utf-8"))
            cases.append((topology, bias, base, source / "imported_fields"))
            topology_gate.append({"topology": topology, "bias_V": -bias, "nodes": 6, "max_coordinate_error_um": 0.0})
    for bias in BIASES:
        cases.append(("coarse7x3", bias, coarse_base, coarse_fields[bias]))

    residuals: list[dict[str, Any]] = []
    updates: list[dict[str, Any]] = []
    jacobians: list[dict[str, Any]] = []
    controls: list[dict[str, Any]] = []
    for topology, bias, base, fields in cases:
        execute_state(args.runner.resolve(), topology, bias, base, fields,
                      output / "raw" / topology / f"m{bias}V",
                      residuals, updates, jacobians, controls)
    write_csv(output / "topology_gate.csv", topology_gate)
    write_csv(output / "residual_decomposition.csv", residuals)
    write_csv(output / "first_qfp_updates.csv", updates)
    write_csv(output / "jacobian_blocks.csv", jacobians)
    write_csv(output / "disabled_controls.csv", controls)

    internal_residuals = [row for row in residuals if not row["is_boundary"]]
    boundaries = [row for row in residuals if row["is_boundary"]]
    max_closure = max(float(row["closure_relative"]) for row in residuals)
    max_jac = max(float(row["rel_diff"]) for row in jacobians if row["block"] == "sg_avalanche")
    boundary_pairs: dict[tuple[Any, ...], dict[str, float]] = {}
    for row in boundaries:
        key = (row["topology"], row["bias_V"], row["carrier"], row["node_id"])
        boundary_pairs.setdefault(key, {})[str(row["variant"])] = float(row["final_residual_normalized"])
    max_boundary_difference = max(abs(pair["production_triangle"] - pair["element_edge_opt_in"]) for pair in boundary_pairs.values())
    update_groups: dict[tuple[Any, ...], dict[str, float]] = {}
    for row in updates:
        key = (row["topology"], row["bias_V"], row["mode"], row["carrier"], row["node_id"])
        update_groups.setdefault(key, {})[str(row["variant"])] = abs(float(row["delta_qfp_V"]))
    topology_metrics: dict[str, dict[str, Any]] = {}
    for topology in sorted({row["topology"] for row in residuals}):
        residual_improved = []
        update_improved = []
        for bias in BIASES:
            for carrier in ("electron", "hole"):
                base_norm = math.sqrt(sum(float(row["final_residual_normalized"]) ** 2 for row in internal_residuals if row["topology"] == topology and row["bias_V"] == -bias and row["carrier"] == carrier and row["variant"] == "production_triangle"))
                candidate_norm = math.sqrt(sum(float(row["final_residual_normalized"]) ** 2 for row in internal_residuals if row["topology"] == topology and row["bias_V"] == -bias and row["carrier"] == carrier and row["variant"] == "element_edge_opt_in"))
                residual_improved.append(candidate_norm < base_norm)
            relevant = [pair for key, pair in update_groups.items() if key[0] == topology and key[1] == -bias]
            update_improved.append(sum(pair["element_edge_opt_in"] for pair in relevant) < sum(pair["production_triangle"] for pair in relevant))
        topology_metrics[topology] = {"residual_all_bias_carriers_improved": all(residual_improved), "first_update_all_bias_improved": all(update_improved)}
    minimal_authorized = all(topology_metrics[name]["residual_all_bias_carriers_improved"] and topology_metrics[name]["first_update_all_bias_improved"] for name in ("minimal6_mirror", "minimal6_sketch"))
    coarse_authorized = topology_metrics["coarse7x3"]["residual_all_bias_carriers_improved"] and topology_metrics["coarse7x3"]["first_update_all_bias_improved"]
    authorized = minimal_authorized and coarse_authorized
    outcome = "source_support_causes_qfp_update" if authorized else "operator_improvement_without_qfp_causality"
    outputs = ["topology_gate.csv", "residual_decomposition.csv", "first_qfp_updates.csv", "jacobian_blocks.csv", "disabled_controls.csv"]
    manifest = {
        "schema": "pn2d_imported_state_qfp_update_v1",
        "biases_V": [-value for value in BIASES],
        "variants": ["production_triangle", "element_edge_opt_in"],
        "global_defaults_unchanged": {"impact_model": "van_overstraeten", "driving_force": "quasi_fermi_gradient", "mobility": "masetti_field_qfp_gradient"},
        "gates": {"residual_closure_max_relative": max_closure, "boundary_max_abs_difference": max_boundary_difference, "sg_avalanche_jacobian_max_relative": max_jac, "minimal6_authorized": minimal_authorized, "coarse7x3_authorized": coarse_authorized, "task8_authorized": authorized},
        "topology_metrics": topology_metrics,
        "typed_outcome": outcome,
        "hashes": {name: sha256(output / name) for name in outputs},
    }
    write_json(output / "manifest.json", manifest)
    print(json.dumps({"typed_outcome": outcome, "task8_authorized": authorized}, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
