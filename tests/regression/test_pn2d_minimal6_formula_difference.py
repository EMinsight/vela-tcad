import csv
import json
import hashlib
import math
import subprocess
import sys
import shutil
import tempfile
import unittest
from pathlib import Path
from scripts.export_pn2d_minimal6_states import collect_member_hashes, validate_member_hashes
from scripts.pn2d_minimal6_diagnostics.counterfactual import validate_formula_input, evaluate_counterfactual_paths, native_source_anchor, integrate_native_nodal_per_unit_depth, integrate_vela_reconstructed_per_unit_depth, sentaurus_alpha_current_nodal, source_log_gap, validate_dependency_dag, interaction_dex, assert_counterfactual_closure, build_adjacent_interactions, symmetric_contributions, score_dominance, validate_field_units, validate_source_anchor_kind
from scripts.pn2d_minimal6_diagnostics.counterfactual import DependencyCounterfactualEngine, FACTOR_DEPENDENCIES
from scripts.diagnose_pn2d_minimal6_formula_difference import _node_state_rows
from scripts.pn2d_minimal6_diagnostics.schemas import validate_formula_difference_v1
from tests.regression.test_pn2d_minimal6_diagnostic_contracts import schema_document, validate_schema_document
import scripts.audit_pn2d_minimal6_fixed_state as fixed_audit

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
    audit_root.mkdir()
    fixed_audit.write_csv(audit_root / "node_state.csv", report.node_rows)
    fixed_audit.write_csv(audit_root / "edge_audit.csv", report.edge_rows)
    fixed_audit.write_csv(audit_root / "triangle_audit.csv", report.triangle_rows)
    summary = dict(report.summary)
    summary["status"] = "PASS"
    (audit_root / "summary.json").write_text(
        json.dumps(summary, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    (audit_root / "manifest.json").write_text(json.dumps({
        "schema": "vela.pn2d_minimal6_fixed_state_audit.v1",
        "row_counts": summary["row_counts"],
    }, sort_keys=True, indent=2) + "\n", encoding="utf-8")
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
             "symmetric_contributions":{"gradient_recovery":1.2, "mobility":0.2}}
            for topology in ("sketch", "mirror") for bias in (-12., -19.)
        ]
        score = score_dominance(states)
        self.assertEqual(score["status"], "available")
        self.assertEqual(score["dominant_factor"], "gradient_recovery")
        states[0]["residual_dex"] = 0.6
        self.assertEqual(score_dominance(states)["status"], "insufficient_data")

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
    def test_adversarial_unit_and_source_kind_contracts(self):
        fields = [{"name":"ImpactIonization", "unit":"cm^-3*s^-1"}, {"name":"eAlphaAvalanche", "unit":"cm^-1"}]
        validate_field_units(fields, {"ImpactIonization":"cm^-3*s^-1", "eAlphaAvalanche":"cm^-1"})
        with self.assertRaises(ValueError):
            validate_field_units(fields, {"ImpactIonization":"m^-3*s^-1"})
        with self.assertRaises(ValueError):
            validate_source_anchor_kind("sentaurus_alpha_current_reconstruction", native=True)
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
    def test_rejects_inexact_bias(self):
        states = [{"topology_id":t,"requested_bias_V":b,"actual_bias_V":b,"status":"passed"}
                  for t in ("sketch","mirror") for b in (0.0,-12.0,-19.0)]
        states[-1]["actual_bias_V"] += 2e-12
        with self.assertRaises(ValueError): validate_formula_input({"outputs_complete":True,"states":states})

if __name__ == '__main__': unittest.main()
