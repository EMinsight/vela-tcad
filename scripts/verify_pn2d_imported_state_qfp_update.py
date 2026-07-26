#!/usr/bin/env python3
"""Independently verify PN2D imported-state QFP-update evidence."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any

FILES = (
    "topology_gate.csv",
    "residual_decomposition.csv",
    "first_qfp_updates.csv",
    "jacobian_blocks.csv",
    "disabled_controls.csv",
    "causality_groups.csv",
)
EXPECTED_COUNTS = {
    "topology_gate.csv": 9,
    "residual_decomposition.csv": 468,
    "first_qfp_updates.csv": 600,
    "jacobian_blocks.csv": 90,
    "disabled_controls.csv": 936,
    "causality_groups.csv": 54,
}
BIASES = (1, 10, 20)
VARIANTS = ("production_triangle", "element_edge_opt_in")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="ascii") as stream:
        return list(csv.DictReader(stream))


def direction(baseline: float, candidate: float) -> str:
    change = abs(candidate - baseline) / max(abs(baseline), abs(candidate), 1.0e-300)
    if change <= 1.0e-12:
        return "equal"
    return "improved" if candidate < baseline else "worsened"


def parse_probe_records(path: Path, prefix: str) -> list[dict[str, str]]:
    output: list[dict[str, str]] = []
    with path.open(encoding="ascii") as stream:
        for line in stream:
            if not line.startswith(prefix + " "):
                continue
            output.append(dict(token.split("=", 1) for token in line.split()[1:] if "=" in token))
    return output


FIELD_CONTRACT = (
    ("ElectrostaticPotential_region0.csv", "psi_V"),
    ("eQuasiFermiPotential_region0.csv", "eQFP_V"),
    ("hQuasiFermiPotential_region0.csv", "hQFP_V"),
)


def expected_impact(variant: str) -> dict[str, Any]:
    common: dict[str, Any] = {
        "model": "van_overstraeten",
        "driving_force": "quasi_fermi_gradient",
        "generation": "current_density",
        "quasi_fermi_gradient_discretization": "cell_gradient",
    }
    if variant == "element_edge_opt_in":
        return {**common, "current_approximation": "element_edge_sg_gss_laux",
                "source_mapping_mode": "element_vertex_box_measure"}
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


def expected_input_hash_keys() -> set[str]:
    keys = {"runner", "coarse_log", "coarse_mesh", "coarse_doping", "materials"}
    for topology in ("mirror", "sketch"):
        for bias in BIASES:
            prefix = f"minimal6_{topology}_m{bias}V"
            keys.add(prefix + "_config")
            keys.add(prefix + "_mesh")
            keys.add(prefix + "_doping")
            for filename, _ in FIELD_CONTRACT:
                keys.add(prefix + "_" + filename)
    return keys


def independently_check_config_lattice(root: Path, errors: list[str]) -> None:
    manifest = json.loads((root / "manifest.json").read_text(encoding="ascii"))
    input_hashes = manifest["input_hashes"]
    sealed_materials = Path(input_hashes["materials"]["path"]).resolve()
    specs = {
        "terms": ("newton_carrier_term_probe", "full"),
        "carrier_only_first_update": ("newton_block_step_probe", "full"),
        "coupled_first_update": ("newton_step_probe", "full"),
        "jacobian": ("newton_jacobian_block_probe", "full"),
        "avalanche_off_terms": ("newton_carrier_term_probe", "avalanche_off"),
        "srh_off_terms": ("newton_carrier_term_probe", "srh_off"),
    }
    for topology in ("minimal6_mirror", "minimal6_sketch", "coarse7x3"):
        for bias in BIASES:
            work = root / "raw" / topology / f"m{bias}V"
            if topology == "coarse7x3":
                expected_mesh = Path(input_hashes["coarse_mesh"]["path"]).resolve()
                expected_doping = Path(input_hashes["coarse_doping"]["path"]).resolve()
                expected_materials = sealed_materials
                expected_fields = (root / "coarse_fields" / f"m{bias}V").resolve()
            else:
                source_topology = topology.removeprefix("minimal6_")
                source_key = f"minimal6_{source_topology}_m{bias}V_config"
                source_path = Path(input_hashes[source_key]["path"]).resolve()
                source_cfg = json.loads(source_path.read_text(encoding="utf-8"))
                expected_mesh = Path(source_cfg["mesh_file"]).resolve()
                expected_doping = Path(source_cfg["node_doping_file"]).resolve()
                expected_materials = Path(source_cfg["materials_file"]).resolve()
                expected_fields = (source_path.parent / "imported_fields").resolve()
            expected_json = {"native_mobility_edges.json"} | {
                f"{variant}_{suffix}.json" for variant in VARIANTS for suffix in specs
            }
            actual_json = {path.name for path in work.glob("*.json")}
            if actual_json != expected_json:
                errors.append(f"{topology} m{bias}V generated configuration lattice mismatch")
                continue
            configs: dict[tuple[str, str], dict[str, Any]] = {}
            for variant in VARIANTS:
                for suffix, (simulation_type, control) in specs.items():
                    cfg = json.loads((work / f"{variant}_{suffix}.json").read_text(encoding="ascii"))
                    configs[(variant, suffix)] = cfg
                    expected_bindings = {
                        "mesh_file": expected_mesh,
                        "node_doping_file": expected_doping,
                        "materials_file": expected_materials,
                        "state_fields_dir": expected_fields,
                    }
                    for key, expected_path in expected_bindings.items():
                        if Path(cfg.get(key, "")).resolve() != expected_path:
                            errors.append(f"{topology} m{bias}V {variant} {suffix} {key} binding drift")
                    expected_output = (work / f"{variant}_{suffix}.csv").resolve()
                    if Path(cfg.get("output_csv", "")).resolve() != expected_output:
                        errors.append(f"{topology} m{bias}V {variant} {suffix} output binding drift")
                    solver = cfg.get("solver", {})
                    expected_recombination = [] if control == "srh_off" else ["srh"]
                    if cfg.get("simulation_type") != simulation_type:
                        errors.append(f"{topology} m{bias}V {variant} {suffix} simulation type drift")
                    if sorted(cfg.get("contacts", []), key=lambda item: item["name"]) != [
                        {"name": "Anode", "bias": -float(bias)},
                        {"name": "Cathode", "bias": 0.0},
                    ]:
                        errors.append(f"{topology} m{bias}V {variant} {suffix} contact drift")
                    if solver.get("method") != "gummel_newton" or solver.get("warm_start") is not True:
                        errors.append(f"{topology} m{bias}V {variant} {suffix} solver method drift")
                    if solver.get("reltol") != 1.0e-8 or solver.get("abstol") != 1.0e-9:
                        errors.append(f"{topology} m{bias}V {variant} {suffix} tolerance drift")
                    if solver.get("mobility") != {"model": "masetti_field", "high_field_driving_force": "quasi_fermi_gradient"}:
                        errors.append(f"{topology} m{bias}V {variant} {suffix} mobility drift")
                    if solver.get("recombination") != expected_recombination:
                        errors.append(f"{topology} m{bias}V {variant} {suffix} recombination control drift")
                    impact = solver.get("impact_ionization")
                    if control == "avalanche_off":
                        if impact is not None:
                            errors.append(f"{topology} m{bias}V {variant} avalanche-off drift")
                    elif impact != expected_impact(variant):
                        errors.append(f"{topology} m{bias}V {variant} {suffix} impact configuration drift")
                    if suffix == "carrier_only_first_update" and cfg.get("block_modes") != ["carrier_only"]:
                        errors.append(f"{topology} m{bias}V {variant} carrier-only mode drift")
                    if suffix != "carrier_only_first_update" and "block_modes" in cfg:
                        errors.append(f"{topology} m{bias}V {variant} {suffix} unexpected block mode")
                    if suffix == "jacobian":
                        if cfg.get("finite_difference_step") != 1.0e-7 or not str(cfg.get("state_file", "")).endswith(f"{variant}_restart.csv"):
                            errors.append(f"{topology} m{bias}V {variant} Jacobian contract drift")
                    elif "finite_difference_step" in cfg or "state_file" in cfg:
                        errors.append(f"{topology} m{bias}V {variant} {suffix} unexpected Jacobian setting")
            for suffix in specs:
                left = json.loads(json.dumps(configs[(VARIANTS[0], suffix)]))
                right = json.loads(json.dumps(configs[(VARIANTS[1], suffix)]))
                for cfg in (left, right):
                    cfg.pop("output_csv", None)
                    cfg.pop("state_file", None)
                    cfg.get("solver", {}).pop("impact_ionization", None)
                if left != right:
                    errors.append(f"{topology} m{bias}V {suffix} non-source variant drift")
            native = json.loads((work / "native_mobility_edges.json").read_text(encoding="ascii"))
            native_bindings = {
                "mesh_file": expected_mesh,
                "node_doping_file": expected_doping,
                "materials_file": expected_materials,
                "state_fields_dir": expected_fields,
                "output_csv": (work / "native_mobility_edges.csv").resolve(),
            }
            if any(Path(native.get(key, "")).resolve() != expected for key, expected in native_bindings.items()):
                errors.append(f"{topology} m{bias}V native path binding drift")
            if native.get("simulation_type") != "sg_edge_flux_probe" or native.get("solver", {}).get("impact_ionization") != expected_impact("production_triangle"):
                errors.append(f"{topology} m{bias}V native mobility probe drift")

def independently_check_topology(root: Path, topology_rows: list[dict[str, str]], errors: list[str]) -> None:
    by_key = {(row["topology"], int(float(row["bias_V"]))): row for row in topology_rows}
    for topology in ("minimal6_mirror", "minimal6_sketch", "coarse7x3"):
        for bias in BIASES:
            cfg_path = root / "raw" / topology / f"m{bias}V" / "production_triangle_terms.json"
            cfg = json.loads(cfg_path.read_text(encoding="ascii"))
            mesh = json.loads(Path(cfg["mesh_file"]).read_text(encoding="utf-8"))
            node_ids = {int(row["id"]) for row in mesh["nodes"]}
            contact_ids = sorted({int(node) for contact in mesh["contacts"] for node in contact["node_ids"]})
            summary = by_key[(topology, -bias)]
            if int(summary["nodes"]) != len(node_ids):
                errors.append(f"{topology} m{bias}V node summary mismatch")
            if summary["contact_nodes"] != ";".join(str(node) for node in contact_ids):
                errors.append(f"{topology} m{bias}V contact summary mismatch")
            fields = Path(cfg["state_fields_dir"])
            for name in (
                "ElectrostaticPotential_region0.csv",
                "eQuasiFermiPotential_region0.csv",
                "hQuasiFermiPotential_region0.csv",
            ):
                field_ids = {int(row["node_id"]) for row in rows(fields / name)}
                if field_ids != node_ids:
                    errors.append(f"{topology} m{bias}V field-node mismatch: {name}")

    coarse_cfg = json.loads((root / "raw" / "coarse7x3" / "m1V" / "production_triangle_terms.json").read_text(encoding="ascii"))
    coarse_mesh = json.loads(Path(coarse_cfg["mesh_file"]).read_text(encoding="utf-8"))
    log_path = Path(json.loads((root / "manifest.json").read_text(encoding="ascii"))["input_hashes"]["coarse_log"]["path"])
    vertex_records = parse_probe_records(log_path, "AVAL_PROBE_VERTEX")
    coarse_node_ids = {int(node["id"]) for node in coarse_mesh["nodes"]}
    for bias in BIASES:
        expected_values = {
            int(row["vertex"]): row for row in vertex_records
            if float(row["bias_V"]) == -float(bias) and int(row["vertex"]) in coarse_node_ids
        }
        if set(expected_values) != coarse_node_ids:
            errors.append(f"coarse m{bias}V frozen-log field coverage mismatch")
            continue
        fields_dir = root / "coarse_fields" / f"m{bias}V"
        for filename, probe_key in FIELD_CONTRACT:
            actual_values = {int(row["node_id"]): float(row["component0"]) for row in rows(fields_dir / filename)}
            expected_field = {node: float(record[probe_key]) for node, record in expected_values.items()}
            if actual_values != expected_field:
                errors.append(f"coarse m{bias}V exact field mismatch: {filename}")
    vertices = {
        int(row["vertex"]): (float(row["x_um"]), float(row["y_um"]))
        for row in vertex_records if float(row["bias_V"]) == -1.0
    }
    for node in coarse_mesh["nodes"]:
        node_id = int(node["id"])
        if node_id not in vertices or vertices[node_id] != (float(node["x"]), float(node["y"])):
            errors.append(f"coarse node-coordinate mismatch at {node_id}")
    measures: dict[int, list[tuple[int, int]]] = {}
    for row in parse_probe_records(log_path, "AVAL_PROBE_MEASURE"):
        if float(row["bias_V"]) != -1.0:
            continue
        measures.setdefault(int(row["element"]), []).append((int(row["local_vertex"]), int(row["vertex"])))
    probe_cells = {cell: [node for _, node in sorted(values)] for cell, values in measures.items()}
    mesh_cells = {int(cell["id"]): [int(node) for node in cell["node_ids"]] for cell in coarse_mesh["triangles"]}
    if probe_cells != mesh_cells:
        errors.append("coarse element/local-vertex mapping mismatch")


def derive_causality(residuals: list[dict[str, str]], updates: list[dict[str, str]]) -> tuple[list[dict[str, Any]], bool]:
    internal = [row for row in residuals if int(row["is_boundary"]) == 0]
    groups: list[dict[str, Any]] = []
    topologies = sorted({row["topology"] for row in residuals})
    for topology in topologies:
        for bias in BIASES:
            for carrier in ("electron", "hole"):
                norms = {
                    variant: math.sqrt(sum(
                        float(row["final_residual_normalized"]) ** 2
                        for row in internal
                        if row["topology"] == topology and float(row["bias_V"]) == -bias
                        and row["carrier"] == carrier and row["variant"] == variant
                    )) for variant in VARIANTS
                }
                groups.append({"topology": topology, "bias_V": -bias, "quantity": "residual", "mode": "assembled", "carrier": carrier, "production_norm": norms[VARIANTS[0]], "candidate_norm": norms[VARIANTS[1]], "direction": direction(norms[VARIANTS[0]], norms[VARIANTS[1]])})
            for mode in ("carrier_only", "coupled"):
                for carrier in ("electron", "hole"):
                    norms = {
                        variant: math.sqrt(sum(
                            float(row["delta_qfp_V"]) ** 2
                            for row in updates
                            if row["topology"] == topology and float(row["bias_V"]) == -bias
                            and row["mode"] == mode and row["carrier"] == carrier
                            and row["variant"] == variant
                        )) for variant in VARIANTS
                    }
                    groups.append({"topology": topology, "bias_V": -bias, "quantity": "first_qfp_update", "mode": mode, "carrier": carrier, "production_norm": norms[VARIANTS[0]], "candidate_norm": norms[VARIANTS[1]], "direction": direction(norms[VARIANTS[0]], norms[VARIANTS[1]])})
    first_biases: dict[str, int | None] = {}
    all_improved: dict[str, bool] = {}
    for topology in topologies:
        active = [bias for bias in BIASES if any(row["topology"] == topology and row["bias_V"] == -bias and row["direction"] != "equal" for row in groups)]
        first = active[0] if active else None
        first_biases[topology] = first
        first_rows = [row for row in groups if row["topology"] == topology and first is not None and row["bias_V"] == -first]
        all_improved[topology] = bool(first_rows) and all(row["direction"] == "improved" for row in first_rows)
    common_first = len({value for value in first_biases.values() if value is not None}) == 1
    return groups, common_first and all(all_improved.values())


def verify(root_a: Path, root_b: Path) -> dict[str, Any]:
    errors: list[str] = []
    manifests = [json.loads((root / "manifest.json").read_text(encoding="ascii")) for root in (root_a, root_b)]
    for index, (root, manifest) in enumerate(zip((root_a, root_b), manifests), start=1):
        if manifest.get("schema") != "pn2d_imported_state_qfp_update_v2":
            errors.append(f"root {index}: schema mismatch")
        for name in FILES:
            if sha256(root / name) != manifest.get("hashes", {}).get(name):
                errors.append(f"root {index}: manifest hash mismatch for {name}")
            if len(rows(root / name)) != EXPECTED_COUNTS[name]:
                errors.append(f"root {index}: row-count mismatch for {name}")
        input_hashes = manifest.get("input_hashes", {})
        if set(input_hashes) != expected_input_hash_keys():
            errors.append(f"root {index}: input hash key set mismatch")
        for name, record in input_hashes.items():
            path = Path(record["path"])
            if not path.is_file() or sha256(path) != record["sha256"]:
                errors.append(f"root {index}: input hash mismatch for {name}")
    for name in FILES:
        if sha256(root_a / name) != sha256(root_b / name):
            errors.append(f"A/B deterministic hash mismatch for {name}")
    if {name: value["sha256"] for name, value in manifests[0]["input_hashes"].items()} != {name: value["sha256"] for name, value in manifests[1]["input_hashes"].items()}:
        errors.append("A/B input provenance mismatch")

    residuals = rows(root_a / "residual_decomposition.csv")
    updates = rows(root_a / "first_qfp_updates.csv")
    jacobians = rows(root_a / "jacobian_blocks.csv")
    controls = rows(root_a / "disabled_controls.csv")
    topology = rows(root_a / "topology_gate.csv")
    for root in (root_a, root_b):
        independently_check_config_lattice(root, errors)
        independently_check_topology(root, rows(root / "topology_gate.csv"), errors)
    if {float(row["bias_V"]) for row in updates} != {-1.0, -10.0, -20.0}:
        errors.append("first-update bias lattice mismatch")
    if {row["mode"] for row in updates} != {"carrier_only", "coupled"}:
        errors.append("first-update mode lattice mismatch")
    if {row["variant"] for row in updates} != set(VARIANTS):
        errors.append("first-update variant lattice mismatch")
    if {row["update_direction_class"] for row in updates} != {"undefined_zero_reference"}:
        errors.append("exact-state direction classification mismatch")
    if not all(math.isfinite(float(row["delta_qfp_V"])) for row in updates):
        errors.append("non-finite first update")

    max_closure = max(float(row["closure_relative"]) for row in residuals)
    if max_closure > 1.0e-12:
        errors.append(f"residual decomposition closure {max_closure}")
    boundary: dict[tuple[Any, ...], dict[str, float]] = {}
    for row in residuals:
        if int(row["is_boundary"]) != 1:
            continue
        key = (row["topology"], row["bias_V"], row["carrier"], row["node_id"])
        boundary.setdefault(key, {})[row["variant"]] = float(row["final_residual_normalized"])
    max_boundary = max(abs(pair[VARIANTS[0]] - pair[VARIANTS[1]]) for pair in boundary.values())
    if max_boundary != 0.0:
        errors.append(f"boundary difference {max_boundary}")

    control_pairs: dict[tuple[Any, ...], dict[str, dict[str, float]]] = {}
    for row in controls:
        key = (row["topology"], row["bias_V"], row["control"], row["carrier"], row["node_id"])
        control_pairs.setdefault(key, {})[row["variant"]] = {name: float(row[name]) for name in ("flux", "srh", "avalanche", "boundary", "residual")}
    if any(set(pair) != set(VARIANTS) for pair in control_pairs.values()):
        errors.append("unpaired disabled control")
    avalanche_off_difference = 0.0
    for key, pair in control_pairs.items():
        if key[2] == "avalanche_off":
            if pair[VARIANTS[0]]["avalanche"] != 0.0 or pair[VARIANTS[1]]["avalanche"] != 0.0:
                errors.append("avalanche-off control contains avalanche source")
            avalanche_off_difference = max(avalanche_off_difference, abs(pair[VARIANTS[0]]["residual"] - pair[VARIANTS[1]]["residual"]))
        if key[2] == "srh_off" and (pair[VARIANTS[0]]["srh"] != 0.0 or pair[VARIANTS[1]]["srh"] != 0.0):
            errors.append("SRH-off control contains SRH source")
    if avalanche_off_difference != 0.0:
        errors.append(f"avalanche-off branch difference {avalanche_off_difference}")

    derived_groups, cross_topology = derive_causality(residuals, updates)
    sealed_groups = rows(root_a / "causality_groups.csv")
    sealed_by_key = {(row["topology"], int(float(row["bias_V"])), row["quantity"], row["mode"], row["carrier"]): row for row in sealed_groups}
    for row in derived_groups:
        key = (row["topology"], row["bias_V"], row["quantity"], row["mode"], row["carrier"])
        sealed = sealed_by_key.get(key)
        if sealed is None or sealed["direction"] != row["direction"]:
            errors.append(f"causality group mismatch for {key}")

    source_rows = [row for row in jacobians if row["block"] == "sg_avalanche"]
    nonzero = [row for row in source_rows if max(abs(float(row["analytic_norm"])), abs(float(row["fd_norm"]))) > 1.0e-12]
    near_zero = [row for row in source_rows if row not in nonzero]
    max_jac_relative = max((float(row["rel_diff"]) for row in nonzero), default=0.0)
    max_jac_absolute = max((float(row["diff_norm"]) for row in near_zero), default=0.0)
    jacobian_gate = max_jac_relative <= 1.0e-8 and max_jac_absolute <= 1.0e-12
    if not jacobian_gate:
        errors.append("source-specific Jacobian gate failed")
    authorized = cross_topology and jacobian_gate
    typed_outcome = "source_support_causes_qfp_update" if authorized else "operator_improvement_without_qfp_causality"
    declared = manifests[0]
    if bool(declared["gates"]["cross_topology_first_material_causality"]) != cross_topology:
        errors.append("declared cross-topology causality mismatch")
    if bool(declared["gates"]["sg_avalanche_jacobian_gate_1e_8"]) != jacobian_gate:
        errors.append("declared Jacobian gate mismatch")
    if bool(declared["gates"]["task8_authorized"]) != authorized:
        errors.append("declared Task 8 authorization mismatch")
    if declared["typed_outcome"] != typed_outcome:
        errors.append("declared typed outcome mismatch")

    return {
        "schema": "pn2d_imported_state_qfp_update_independent_verification_v2",
        "pass": not errors,
        "errors": errors,
        "sealed_hashes_match": all(sha256(root_a / name) == sha256(root_b / name) for name in FILES),
        "input_hashes_verified": not any("input hash" in error for error in errors),
        "residual_closure_max_relative": max_closure,
        "boundary_max_abs_difference": max_boundary,
        "avalanche_off_branch_max_abs_difference": avalanche_off_difference,
        "sg_avalanche_jacobian_max_nonzero_relative": max_jac_relative,
        "sg_avalanche_jacobian_max_near_zero_absolute": max_jac_absolute,
        "sg_avalanche_jacobian_gate_1e_8": jacobian_gate,
        "cross_topology_first_material_causality": cross_topology,
        "task8_authorized": authorized,
        "typed_outcome": typed_outcome,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root-a", type=Path, required=True)
    parser.add_argument("--root-b", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = verify(args.root_a.resolve(), args.root_b.resolve())
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="ascii", newline="\n")
    print(json.dumps({"pass": result["pass"], "typed_outcome": result["typed_outcome"]}, sort_keys=True))
    return 0 if result["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
