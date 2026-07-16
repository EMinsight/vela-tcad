import csv
import json
import hashlib
import math
import subprocess
import sys
import shutil
import tempfile
import unittest
from unittest import mock
from pathlib import Path
from scripts.export_pn2d_minimal6_states import collect_member_hashes, validate_member_hashes
from scripts.pn2d_minimal6_diagnostics.counterfactual import validate_formula_input, evaluate_counterfactual_paths, native_source_anchor, integrate_native_nodal_per_unit_depth, integrate_vela_reconstructed_per_unit_depth, sentaurus_alpha_current_nodal, source_log_gap, validate_dependency_dag, interaction_dex, assert_counterfactual_closure, build_adjacent_interactions, symmetric_contributions, score_dominance, validate_field_units, validate_source_anchor_kind
from scripts.pn2d_minimal6_diagnostics.counterfactual import DependencyCounterfactualEngine, FACTOR_DEPENDENCIES, make_formula_operator_engine, evaluate_formula_counterfactual
from scripts.diagnose_pn2d_minimal6_formula_difference import _node_state_rows, validate_audit_binding
from scripts.pn2d_minimal6_diagnostics.schemas import validate_formula_difference_v1
from tests.regression.test_pn2d_minimal6_diagnostic_contracts import schema_document, validate_schema_document
import scripts.audit_pn2d_minimal6_fixed_state as fixed_audit
import scripts.diagnose_pn2d_minimal6_formula_difference as formula_cli
from scripts.pn2d_minimal6_diagnostics.plots import render_formula_difference_figures

def _prepare_formula_fixture(temp: str) -> tuple[Path, Path]:
    source = Path(__file__).parents[1] / "fixtures" / "pn2d_minimal6_synthetic"
    state_root = Path(temp) / "state"
    shutil.copytree(source, state_root)
    manifest_path = state_root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    installed = {
        "ImpactIonization": "cm^-3*s^-1",
        "eVelocity": "cm*s^-1",
        "hVelocity": "cm*s^-1",
        "eIonIntegral": "1",
        "hIonIntegral": "1",
        "MeanIonIntegral": "1",
        "LatticeTemperature": "K",
    }
    for state_index, state in enumerate(manifest["states"]):
        export = state_root / state["export_dir"]
        bundle = export.parent / "source"
        bundle.mkdir(exist_ok=True)
        shutil.copy2(
            Path(__file__).parents[2] / "reference_tcad" / "pn2d_sentaurus2018_minimal6" / "source" / "models.par",
            bundle / "models.par",
        )
        state["bundle_dir"] = bundle.relative_to(state_root).as_posix()
        fields = export / "fields"
        node_ids = sorted(_read_node_ids(fields / "eDensity_region0.csv"))
        bias_scale = {0.0: 1.0, -12.0: 1.0e4, -19.0: 1.0e7}[float(state["requested_bias_V"])]
        topology_scale = 1.0 if state["topology_id"] == "sketch" else 1.25
        values = {
            "ImpactIonization": [bias_scale * topology_scale * (index + 1) * 1.0e8 for index in range(6)],
            "eVelocity": [1.0e5 + 100.0 * index for index in range(6)],
            "hVelocity": [0.5e5 + 80.0 * index for index in range(6)],
            "eIonIntegral": [0.01 * bias_scale * (index + 1) for index in range(6)],
            "hIonIntegral": [0.005 * bias_scale * (index + 1) for index in range(6)],
            "MeanIonIntegral": [0.0075 * bias_scale * (index + 1) for index in range(6)],
            "LatticeTemperature": [300.0 for _ in range(6)],
        }
        for name, field_values in values.items():
            rows = "node_id,component0\n" + "".join(
                f"{node_id},{value:.17g}\n" for node_id, value in zip(node_ids, field_values)
            )
            (fields / f"{name}_region0.csv").write_text(rows, encoding="utf-8")
        field_manifest_path = export / "field_manifest.json"
        field_manifest = json.loads(field_manifest_path.read_text(encoding="utf-8"))
        names = {row["name"] for row in field_manifest["fields"]}
        for name, unit in installed.items():
            if name not in names:
                field_manifest["fields"].append({
                    "name": name, "region": 0, "components": 1, "unit": unit,
                    "mapping_status": "complete",
                    "global_node_mapping": "global_vertex_order",
                })
        field_manifest_path.write_text(
            json.dumps(field_manifest, sort_keys=True, indent=2) + "\n", encoding="utf-8"
        )
    manifest_path.write_text(json.dumps(manifest, sort_keys=True, indent=2) + "\n", encoding="utf-8")

    report = fixed_audit.build_report(state_root)
    audit_root = Path(temp) / "audit"
    fixed_audit.write_report(report, audit_root)
    return state_root, audit_root


def _read_node_ids(path: Path) -> list[int]:
    with path.open(newline="", encoding="utf-8") as handle:
        return [int(row["node_id"]) for row in csv.DictReader(handle)]

class FormulaDifferenceTest(unittest.TestCase):
    def test_requires_exact_six_state_matrix_and_emits_named_residual(self):
        states = [{"topology_id":t,"requested_bias_V":b,"actual_bias_V":b,"status":"passed"}
                  for t in ("sketch","mirror") for b in (0.0,-12.0,-19.0)]
        report = validate_formula_input({"outputs_complete":True,"states":states})
        self.assertEqual(report["row_counts"], {"node":36,"edge":54,"triangle":24})
        self.assertIn("sentaurus_internal_semantics_residual", report)
    def test_forward_reverse_paths_and_residual_close_exactly(self):
        factors = {"ni_eff/BGN": 1.0, "mobility": 10.0, "alpha_law": 100.0}
        dependencies = {"ni_eff/BGN": (), "mobility": ("ni_eff/BGN",), "alpha_law": ("mobility",)}
        result = evaluate_counterfactual_paths(native=1000.0, baseline=1.0, factors=factors, dependencies=dependencies)
        self.assertEqual(result["forward"]["order"], list(factors))
        self.assertEqual(result["reverse"]["order"], list(reversed(factors)))
        self.assertAlmostEqual(result["residual_dex"], 0.0)
    def test_dependency_dag_interaction_and_closure_contracts(self):
        order = validate_dependency_dag({"ni_eff/BGN": (), "mobility": ("ni_eff/BGN",), "alpha_law": ("mobility",)})
        self.assertEqual(order, ["ni_eff/BGN", "mobility", "alpha_law"])
        with self.assertRaises(ValueError):
            validate_dependency_dag({"mobility": ("missing",)})
        self.assertAlmostEqual(interaction_dex(baseline=1., a_only=10., b_only=100., both=1000.), 0.0)
        assert_counterfactual_closure(native_gap_dex=2., contributions_dex=[0.5, 1.5], residual_dex=0.)
        with self.assertRaises(ValueError):
            assert_counterfactual_closure(native_gap_dex=2., contributions_dex=[0.5], residual_dex=0.)

    def test_dependency_engine_recomputes_only_declared_downstream_operators(self):
        dependencies = {"a": (), "b": ("a",), "c": ("b",)}
        operators = {
            "a": lambda inputs, raw: raw,
            "b": lambda inputs, raw: inputs["a"] * raw,
            "c": lambda inputs, raw: inputs["b"] * raw,
        }
        engine = DependencyCounterfactualEngine(
            dependencies=dependencies,
            baseline_values={"a": 1.0, "b": 1.0, "c": 1.0},
            replacement_values={"a": 2.0, "b": 3.0, "c": 5.0},
            operators=operators,
            output_factor="c",
        )
        paths = engine.evaluate_paths(native=30.0)
        self.assertEqual(
            [row["recomputed"] for row in paths["forward"]["contributions"]],
            [["a", "b", "c"], ["b", "c"], ["c"]],
        )
        self.assertEqual(
            [row["recomputed"] for row in paths["reverse"]["contributions"]],
            [["c"], ["b", "c"], ["a", "b", "c"]],
        )
        self.assertAlmostEqual(paths["residual_dex"], 0.0, places=14)

        invalid = dict(operators)
        invalid["b"] = lambda inputs, raw: inputs["c"] * raw
        with self.assertRaisesRegex(ValueError, "undeclared dependency"):
            DependencyCounterfactualEngine(
                dependencies=dependencies,
                baseline_values={"a": 1., "b": 1., "c": 1.},
                replacement_values={"a": 2., "b": 3., "c": 5.},
                operators=invalid,
                output_factor="c",
            ).evaluate_paths(native=30.)
    def test_all_eight_formula_operators_independently_change_source(self):
        baseline = {
            "ni_eff/BGN": (1.0, 2.0),
            "gradient_recovery": (2.0, 3.0),
            "mobility": (3.0, 4.0),
            "current_semantics": None,
            "impact_driving_field": (2.0, 2.0),
            "alpha_law": ((2.0, 1.0), (2.0, 1.0)),
            "partial_volume": (1.0, 1.0),
            "source_to_node_mapping": (1.0, 1.0),
        }
        replacements = {
            "ni_eff/BGN": (2.0, 4.0),
            "gradient_recovery": (4.0, 6.0),
            "mobility": (6.0, 8.0),
            "current_semantics": (20.0, 30.0),
            "impact_driving_field": (4.0, 4.0),
            "alpha_law": ((3.0, 1.0), (3.0, 1.0)),
            "partial_volume": (2.0, 2.0),
            "source_to_node_mapping": (2.0, 2.0),
        }
        engine = make_formula_operator_engine(
            baseline_values=baseline, replacement_values=replacements
        )
        base = engine.evaluate_replacements(set())
        for factor in FACTOR_DEPENDENCIES:
            with self.subTest(factor=factor):
                self.assertNotEqual(engine.evaluate_replacements({factor}), base)

    def test_unavailable_formula_input_stays_in_named_residual(self):
        baseline = {
            "ni_eff/BGN": (1.0,), "gradient_recovery": (2.0,),
            "mobility": (3.0,), "current_semantics": None,
            "impact_driving_field": (2.0,), "alpha_law": ((2.0, 1.0),),
            "partial_volume": (1.0,), "source_to_node_mapping": (1.0,),
        }
        replacements = {
            "ni_eff/BGN": (2.0,), "gradient_recovery": (4.0,),
            "mobility": (6.0,), "current_semantics": (20.0,),
            "impact_driving_field": (4.0,), "alpha_law": ((3.0, 1.0),),
            "partial_volume": (2.0,), "source_to_node_mapping": (2.0,),
        }
        native = make_formula_operator_engine(
            baseline_values=baseline, replacement_values=replacements
        ).evaluate_replacements(set(FACTOR_DEPENDENCIES))
        missing = dict(replacements)
        del missing["alpha_law"]
        result = evaluate_formula_counterfactual(
            native=native, baseline_values=baseline, replacement_values=missing,
            unavailable_reasons={"alpha_law": "missing Sentaurus coefficient provenance"},
        )
        alpha = next(row for row in result["factor_availability"] if row["factor"] == "alpha_law")
        self.assertEqual(alpha["status"], "unavailable")
        self.assertIn("missing Sentaurus", alpha["reason"])
        self.assertNotEqual(result["residual_dex"], 0.0)
        for identity in ("forward", "reverse"):
            closure = sum(row["contribution_dex"] for row in result[identity]["contributions"]) + result["residual_dex"]
            self.assertLessEqual(abs(result["native_gap_dex"] - closure), 1.0e-10)

    def test_interaction_and_dominance_gate_require_complete_matrix(self):
        forward = [{"factor":"gradient_recovery", "contribution_dex":0.8}, {"factor":"mobility", "contribution_dex":0.1}]
        reverse = [{"factor":"mobility", "contribution_dex":0.1}, {"factor":"gradient_recovery", "contribution_dex":0.2}]
        source = {frozenset():1., frozenset({"gradient_recovery"}):10., frozenset({"mobility"}):2., frozenset({"gradient_recovery", "mobility"}):30.}
        interactions = build_adjacent_interactions(forward, reverse, lambda replaced: source[frozenset(replaced)])
        self.assertEqual(len(interactions), 2)
        self.assertAlmostEqual(interactions[0]["interaction_dex"], math.log10(1.5))
        symmetric = symmetric_contributions(forward, reverse)
        self.assertAlmostEqual(symmetric["gradient_recovery"], 0.5)
        states = [
            {"topology":topology, "bias_V":bias, "native_gap_dex":2., "residual_dex":0.1,
             "symmetric_contributions":{"gradient_recovery":1.2, "mobility":0.2},
             "factor_availability":[{"factor":"gradient_recovery","status":"available"},
                                    {"factor":"mobility","status":"available"}]}
            for topology in ("sketch", "mirror") for bias in (-12., -19.)
        ]
        score = score_dominance(states)
        self.assertEqual(score["status"], "available")
        self.assertEqual(score["dominant_factor"], "gradient_recovery")
        states[0]["residual_dex"] = 0.6
        self.assertEqual(score_dominance(states)["status"], "insufficient_data")
        states[0]["residual_dex"] = 0.1
        states[0]["factor_availability"][0] = {
            "factor": "gradient_recovery", "status": "unavailable",
            "reason": "missing provenance",
        }
        unavailable = score_dominance(states)
        self.assertEqual(unavailable["status"], "insufficient_data")
        self.assertNotIn("dominant_factor", unavailable)

    def test_adjacent_interactions_include_forward_and_reverse_path_identities(self):
        forward = [
            {"factor": "gradient_recovery", "contribution_dex": 0.8},
            {"factor": "mobility", "contribution_dex": 0.1},
        ]
        reverse = [
            {"factor": "mobility", "contribution_dex": 0.1},
            {"factor": "gradient_recovery", "contribution_dex": 0.2},
        ]
        sources = {
            frozenset(): 1.0,
            frozenset({"gradient_recovery"}): 10.0,
            frozenset({"mobility"}): 2.0,
            frozenset({"gradient_recovery", "mobility"}): 30.0,
        }
        interactions = build_adjacent_interactions(
            forward, reverse, lambda replaced: sources[frozenset(replaced)]
        )
        self.assertEqual(
            {row["path_identity"] for row in interactions},
            {"forward_adjacent", "reverse_adjacent"},
        )
    def test_hash_mutation_is_rejected(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            member = root / "immutable.csv"
            member.write_text("original\n", encoding="utf-8")
            hashes = collect_member_hashes(root)
            member.write_text("mutated\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                validate_member_hashes(root, hashes)
    def test_audit_binding_rejects_fabrication_state_mismatch_and_artifact_mutation(self):
        with tempfile.TemporaryDirectory() as temp:
            state_root, audit_root = _prepare_formula_fixture(temp)
            validate_audit_binding(state_root, audit_root)

            fabricated = Path(temp) / "fabricated"
            shutil.copytree(audit_root, fabricated)
            fabricated_manifest = json.loads((fabricated / "manifest.json").read_text(encoding="utf-8"))
            del fabricated_manifest["input_sha256"]
            (fabricated / "manifest.json").write_text(json.dumps(fabricated_manifest), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "input_sha256"):
                validate_audit_binding(state_root, fabricated)

            state_manifest_path = state_root / "manifest.json"
            state_manifest = json.loads(state_manifest_path.read_text(encoding="utf-8"))
            state_manifest["run_id"] = "unrelated-run"
            state_manifest_path.write_text(json.dumps(state_manifest), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "state.*hash|run"):
                validate_audit_binding(state_root, audit_root)
            state_manifest["run_id"] = "committed_task5_cpp_replay"
            state_manifest_path.write_text(json.dumps(state_manifest, sort_keys=True, indent=2) + "\n", encoding="utf-8")

            mutated = Path(temp) / "mutated"
            shutil.copytree(audit_root, mutated)
            edge_path = mutated / "edge_audit.csv"
            edge_path.write_bytes(edge_path.read_bytes() + b"\n")
            with self.assertRaisesRegex(ValueError, "edge_audit.*mutation|artifact"):
                validate_audit_binding(state_root, mutated)

    def test_adversarial_unit_and_source_kind_contracts(self):
        fields = [{"name":"ImpactIonization", "unit":"cm^-3*s^-1"}, {"name":"eAlphaAvalanche", "unit":"cm^-1"}]
        validate_field_units(fields, {"ImpactIonization":"cm^-3*s^-1", "eAlphaAvalanche":"cm^-1"})
        with self.assertRaises(ValueError):
            validate_field_units(fields, {"ImpactIonization":"m^-3*s^-1"})
        validate_source_anchor_kind(
            "sentaurus_native_avalanche_generation", "sentaurus", native=True
        )
        validate_source_anchor_kind(
            "sentaurus_alpha_current_reconstruction", "derived", native=False
        )
        with self.assertRaises(ValueError):
            validate_source_anchor_kind(
                "sentaurus_alpha_current_reconstruction", "derived", native=True
            )
        with self.assertRaises(ValueError):
            validate_source_anchor_kind(
                "sentaurus_native_avalanche_generation", "derived", native=True
            )
        with self.assertRaises(ValueError):
            validate_source_anchor_kind(
                "sentaurus_alpha_current_reconstruction", "sentaurus", native=False
            )
        mesh = {"nodes":[{"id":0,"x":0.,"y":0.},{"id":0,"x":1.,"y":0.},{"id":2,"x":0.,"y":1.}], "triangles":[{"node_ids":[0,2,0]}]}
        with self.assertRaises(ValueError):
            integrate_native_nodal_per_unit_depth(mesh, {0:1., 2:1.})
    def test_native_nodal_anchor_rejects_reversed_topology(self):
        mesh = {"nodes":[{"id":0,"x":0.,"y":0.},{"id":1,"x":1.e-6,"y":0.},{"id":2,"x":0.,"y":1.e-6}], "triangles":[{"node_ids":[0,2,1]}]}
        with self.assertRaises(ValueError):
            integrate_native_nodal_per_unit_depth(mesh, {0:1., 1:1., 2:1.})
    def test_node_ledger_rejects_missing_raw_field(self):
        with tempfile.TemporaryDirectory() as temp:
            fields = Path(temp) / "fields"
            fields.mkdir()
            (fields / "ElectrostaticPotential_region0.csv").write_text("node_id,component0\n0,0\n", encoding="utf-8")
            with self.assertRaises(FileNotFoundError):
                list(_node_state_rows({"topology":"sketch", "bias_V":0., "export_dir":temp}))
    def test_native_anchor_refuses_missing_volume(self):
        result = native_source_anchor([1.0, 2.0], volume_m3=None)
        self.assertEqual(result["status"], "insufficient_data")
        self.assertIsNone(result["value"])
    def test_native_nodal_anchor_uses_explicit_unit_depth(self):
        mesh = {"nodes":[{"id":0,"x":0.,"y":0.},{"id":1,"x":1.e-6,"y":0.},{"id":2,"x":0.,"y":1.e-6}],"triangles":[{"node_ids":[0,1,2]}]}
        result = integrate_native_nodal_per_unit_depth(mesh, {0:3.,1:6.,2:9.})
        self.assertEqual(result["status"], "available")
        self.assertAlmostEqual(result["value_s_inv_per_unit_depth"], 3.0e-8)
        self.assertEqual(result["depth_convention"], "unit_out_of_plane_length_cm")
    def test_vela_reconstructed_source_uses_unit_depth_conversion(self):
        rows = [{"local_edge0_electron_source_integral_per_m_s":"2", "local_edge0_hole_source_integral_per_m_s":"3"}]
        self.assertAlmostEqual(integrate_vela_reconstructed_per_unit_depth(rows), 0.05)
    def test_sentaurus_alpha_current_reconstruction_is_explicit(self):
        value = sentaurus_alpha_current_nodal({0:2.}, {0:(3.,4.)}, {0:0.}, {0:(0.,0.)}, elementary_charge_C=1.)
        self.assertEqual(value, {0:10.})
    def test_source_log_gap_classifies_zero_and_reports_dex(self):
        self.assertAlmostEqual(source_log_gap(100., 1.)["dex"], 2.0)
        self.assertEqual(source_log_gap(0., 0.)["classification"], "geometric_zero")
    def test_real_mapping_is_conservative_and_never_calibrated_to_target(self):
        with tempfile.TemporaryDirectory() as temp:
            state_root, audit_root = _prepare_formula_fixture(temp)
            manifest = json.loads((state_root / "manifest.json").read_text(encoding="utf-8"))
            audit = formula_cli.validate_audit_binding(state_root, audit_root)
            state = next(
                row for row in manifest["states"]
                if row["topology_id"] == "sketch" and row["requested_bias_V"] == -19.0
            )
            resolved = formula_cli._resolved_state(state, state_root)
            audit_state = {
                kind: [
                    row for row in audit[kind]
                    if row["topology_id"] == "sketch" and float(row["bias_V"]) == -19.0
                ]
                for kind in ("node", "edge", "triangle")
            }
            record = formula_cli._state_sources(resolved, audit_state)
            changed_target = dict(record)
            target_name = "sentaurus_alpha_current_reconstruction_s_inv_per_unit_depth"
            changed_target[target_name] *= 17.0
            baseline, replacements, unavailable = formula_cli._formula_operator_inputs(record)
            changed_baseline, changed_replacements, _ = formula_cli._formula_operator_inputs(changed_target)
            factor = "source_to_node_mapping"
            self.assertEqual(baseline[factor], changed_baseline[factor])
            self.assertEqual(replacements.get(factor), changed_replacements.get(factor))
            self.assertNotIn(factor, replacements)
            self.assertIn("conservative", unavailable[factor])
            for triangle in audit_state["triangle"]:
                for source in ("vela", "python"):
                    for carrier in ("electron", "hole"):
                        self.assertAlmostEqual(
                            formula_cli._mapping_scale(triangle, source, carrier), 1.0
                        )
    def test_cli_writes_closed_deterministic_exact_identity_artifacts(self):
        with tempfile.TemporaryDirectory() as temp:
            state_root, audit_root = _prepare_formula_fixture(temp)
            outputs = [Path(temp) / "out-a", Path(temp) / "out-b"]
            command_prefix = [
                sys.executable,
                str(Path(__file__).parents[2] / "scripts" / "diagnose_pn2d_minimal6_formula_difference.py"),
                "--state-root", str(state_root),
                "--audit-root", str(audit_root),
            ]
            for out in outputs:
                completed = subprocess.run(
                    command_prefix + ["--out-dir", str(out), "--qa-status", "reviewed"],
                    capture_output=True, text=True,
                )
                self.assertEqual(completed.returncode, 0, completed.stderr)

            artifacts = (
                "quantity_ledger.csv", "factor_waterfall.csv",
                "root_cause_summary.json", "root_cause_summary.md",
            )
            first_hashes = {
                name: hashlib.sha256((outputs[0] / name).read_bytes()).hexdigest()
                for name in artifacts
            }
            second_hashes = {
                name: hashlib.sha256((outputs[1] / name).read_bytes()).hexdigest()
                for name in artifacts
            }
            self.assertEqual(first_hashes, second_hashes)

            with (outputs[0] / "quantity_ledger.csv").open(newline="", encoding="utf-8") as handle:
                ledger = list(csv.DictReader(handle))
            with (outputs[0] / "factor_waterfall.csv").open(newline="", encoding="utf-8") as handle:
                waterfall = list(csv.DictReader(handle))
            report = json.loads((outputs[0] / "root_cause_summary.json").read_text(encoding="utf-8"))
            validate_formula_difference_v1(report)
            validate_schema_document(report, schema_document("vela.pn2d_minimal6_formula_difference.v1"))

            node_ids = {(row["topology"], row["bias_V"], row["node_id"])
                        for row in ledger if row["record_kind"] == "node_state"}
            edge_ids = {(row["topology"], row["bias_V"], row["edge_id"])
                        for row in ledger if row["record_kind"] == "edge_raw"}
            cell_ids = {(row["topology"], row["bias_V"], row["cell_id"])
                        for row in ledger if row["record_kind"] == "cell_replay"}
            self.assertEqual((len(node_ids), len(edge_ids), len(cell_ids)), (36, 54, 24))

            source_rows = [row for row in ledger if row["record_kind"] == "source_integral"]
            self.assertEqual(len(source_rows), 18)
            kinds = {row["source"]: row["source_kind"] for row in source_rows}
            self.assertEqual(kinds["sentaurus_native_avalanche_generation"], "sentaurus")
            self.assertEqual(kinds["sentaurus_alpha_current_reconstruction"], "derived")
            self.assertEqual(kinds["vela_alpha_flux_partial_volume_reconstruction"], "derived")
            self.assertTrue(any(row["quantity"] == "eIonIntegral" for row in ledger))
            self.assertTrue(any(row["quantity"] == "ni_eff_relative_residual" for row in ledger))

            direction_rows = [row for row in ledger if row["component"] == "direction_rad"]
            self.assertTrue(direction_rows)
            self.assertEqual({row["unit"] for row in direction_rows}, {"rad"})

            unavailable_factors = {
                "ni_eff/BGN": "independently inferred ni_eff/BGN",
                "impact_driving_field": "coefficient provenance",
                "alpha_law": "coefficient provenance",
                "source_to_node_mapping": "conservative",
            }
            for path in report["waterfall_paths"]:
                availability = {row["factor"]: row for row in path["factor_availability"]}
                for factor, reason in unavailable_factors.items():
                    self.assertEqual(availability[factor]["status"], "unavailable")
                    self.assertIn(reason, availability[factor]["reason"])
                for identity in ("forward", "reverse"):
                    contributions = {
                        row["factor"]: row["contribution_dex"]
                        for row in path[identity]["contributions"]
                    }
                    for factor in unavailable_factors:
                        self.assertAlmostEqual(contributions[factor], 0.0)
            self.assertEqual(report["dominance_rules"]["status"], "insufficient_data")
            self.assertNotIn("dominant_factor", report["dominance_rules"])

            residual_by_state = {
                (row["topology"], row["bias_V"]): row["dex"]
                for row in report["sentaurus_internal_semantics_residual"]
            }
            factor_order = list(FACTOR_DEPENDENCIES)
            for path in report["waterfall_paths"]:
                self.assertEqual(path["forward"]["order"], factor_order)
                self.assertEqual(path["reverse"]["order"], list(reversed(factor_order)))
                self.assertLessEqual(abs(path["residual_dex"] - residual_by_state[(path["topology"], path["bias_V"])]), 1.0e-10)
                if path["bias_V"] in (-12.0, -19.0):
                    for identity in ("forward", "reverse"):
                        closure = sum(
                            row["contribution_dex"] for row in path[identity]["contributions"]
                        ) + path["residual_dex"]
                        self.assertLessEqual(abs(path["native_gap_dex"] - closure), 1.0e-10)

            groups = {}
            for row in waterfall:
                if row["path_identity"] not in ("forward", "reverse"):
                    continue
                key = (row["topology"], row["bias_V"], row["path_identity"])
                groups.setdefault(key, []).append(row)
            self.assertEqual(len(groups), 12)
            for (topology, bias, identity), rows in groups.items():
                rows.sort(key=lambda row: int(row["order_index"]))
                expected = factor_order if identity == "forward" else list(reversed(factor_order))
                self.assertEqual([row["factor"] for row in rows], expected)

            self.assertEqual(len(report["sentaurus_internal_semantics_residual"]), 6)
            self.assertTrue(all(
                row["name"] == "sentaurus_internal_semantics_residual"
                for row in report["sentaurus_internal_semantics_residual"]
            ))
            if report["dominance_rules"]["status"] == "insufficient_data":
                self.assertNotIn("dominant_factor", report["dominance_rules"])
    def test_triggered_interaction_renders_first_second_factor_contract(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            ledger = root / "quantity_ledger.csv"
            waterfall = root / "factor_waterfall.csv"
            report = root / "root_cause_summary.json"
            ledger.write_text("record_kind\n", encoding="utf-8")
            waterfall.write_text("path_identity\n", encoding="utf-8")
            report.write_text(json.dumps({
                "interactions": [{
                    "topology": "sketch", "bias_V": -19.0,
                    "first_factor": "gradient_recovery",
                    "second_factor": "mobility",
                    "path_identity": "forward_adjacent",
                    "baseline": 1.0, "a_only": 10.0, "b_only": 2.0,
                    "both": 50.0, "interaction_dex": math.log10(2.5),
                }]
            }), encoding="utf-8")
            manifest = render_formula_difference_figures(
                ledger_path=ledger, waterfall_path=waterfall,
                report_path=report, out_dir=root, qa_status="reviewed",
            )
            self.assertGreater(abs(math.log10(2.5)), 0.3)
            self.assertTrue((root / "interaction.png").is_file())
            self.assertTrue((root / "interaction.pdf").is_file())
            self.assertEqual(
                next(row for row in manifest["figures"] if row["stem"] == "interaction")["unit"],
                "dex",
            )

    def test_cli_adversarial_file_mutations_fail_before_artifacts(self):
        with tempfile.TemporaryDirectory() as temp:
            state_root, audit_root = _prepare_formula_fixture(temp)
            cases = {}

            def missing_field(root, _audit):
                state = json.loads((root / "manifest.json").read_text(encoding="utf-8"))["states"][0]
                (root / state["export_dir"] / "fields" / "eVelocity_region0.csv").unlink()

            def wrong_unit(root, _audit):
                state = json.loads((root / "manifest.json").read_text(encoding="utf-8"))["states"][0]
                path = root / state["export_dir"] / "field_manifest.json"
                fields = json.loads(path.read_text(encoding="utf-8"))
                next(row for row in fields["fields"] if row["name"] == "ImpactIonization")["unit"] = "m^-3*s^-1"
                path.write_text(json.dumps(fields), encoding="utf-8")

            def reversed_topology(root, _audit):
                state = json.loads((root / "manifest.json").read_text(encoding="utf-8"))["states"][0]
                path = root / state["export_dir"] / "mesh.json"
                mesh = json.loads(path.read_text(encoding="utf-8"))
                mesh["triangles"][0]["node_ids"].reverse()
                path.write_text(json.dumps(mesh), encoding="utf-8")

            def duplicate_node(root, _audit):
                state = json.loads((root / "manifest.json").read_text(encoding="utf-8"))["states"][0]
                path = root / state["export_dir"] / "mesh.json"
                mesh = json.loads(path.read_text(encoding="utf-8"))
                mesh["nodes"][1]["id"] = mesh["nodes"][0]["id"]
                path.write_text(json.dumps(mesh), encoding="utf-8")

            def state_hash(root, _audit):
                path = root / "manifest.json"
                manifest = json.loads(path.read_text(encoding="utf-8"))
                manifest["run_id"] = "mutated"
                path.write_text(json.dumps(manifest), encoding="utf-8")

            def audit_hash(_root, audit):
                path = audit / "edge_audit.csv"
                path.write_bytes(path.read_bytes() + b"\n")

            cases.update({
                "missing_field": missing_field,
                "wrong_unit": wrong_unit,
                "reversed_topology": reversed_topology,
                "duplicate_node": duplicate_node,
                "state_hash": state_hash,
                "audit_hash": audit_hash,
            })
            for name, mutate in cases.items():
                with self.subTest(name=name):
                    case_root = Path(temp) / f"state-{name}"
                    case_audit = Path(temp) / f"audit-{name}"
                    shutil.copytree(state_root, case_root)
                    shutil.copytree(audit_root, case_audit)
                    mutate(case_root, case_audit)
                    out = Path(temp) / f"out-{name}"
                    completed = subprocess.run([
                        sys.executable,
                        str(Path(__file__).parents[2] / "scripts" / "diagnose_pn2d_minimal6_formula_difference.py"),
                        "--state-root", str(case_root),
                        "--audit-root", str(case_audit),
                        "--out-dir", str(out),
                    ], capture_output=True, text=True)
                    self.assertNotEqual(completed.returncode, 0)
                    self.assertFalse(out.exists())

    def test_cli_rejects_undeclared_dependency_and_false_native_kind_before_artifacts(self):
        with tempfile.TemporaryDirectory() as temp:
            state_root, audit_root = _prepare_formula_fixture(temp)
            argv = [
                "diagnose_pn2d_minimal6_formula_difference.py",
                "--state-root", str(state_root),
                "--audit-root", str(audit_root),
            ]
            out_dependency = Path(temp) / "out-dependency"
            with mock.patch.object(
                formula_cli, "FACTOR_DEPENDENCIES", {"mobility": ("missing",)}
            ), mock.patch.object(
                sys, "argv", argv + ["--out-dir", str(out_dependency)]
            ):
                with self.assertRaisesRegex(ValueError, "undeclared"):
                    formula_cli.main()
            self.assertFalse(out_dependency.exists())

            out_kind = Path(temp) / "out-kind"
            false_kinds = dict(formula_cli.SOURCE_FAMILIES)
            false_kinds["sentaurus_native_avalanche_generation"] = "derived"
            with mock.patch.object(
                formula_cli, "SOURCE_FAMILIES", false_kinds
            ), mock.patch.object(
                sys, "argv", argv + ["--out-dir", str(out_kind)]
            ):
                with self.assertRaisesRegex(ValueError, "native|SourceKind"):
                    formula_cli.main()
            self.assertFalse(out_kind.exists())
    def test_rejects_inexact_bias(self):
        states = [{"topology_id":t,"requested_bias_V":b,"actual_bias_V":b,"status":"passed"}
                  for t in ("sketch","mirror") for b in (0.0,-12.0,-19.0)]
        states[-1]["actual_bias_V"] += 2e-12
        with self.assertRaises(ValueError): validate_formula_input({"outputs_complete":True,"states":states})

if __name__ == '__main__': unittest.main()
