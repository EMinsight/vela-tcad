import csv
import json
import math
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path
from scripts.export_pn2d_minimal6_states import collect_member_hashes, validate_member_hashes
from scripts.pn2d_minimal6_diagnostics.counterfactual import validate_formula_input, evaluate_counterfactual_paths, native_source_anchor, integrate_native_nodal_per_unit_depth, integrate_vela_reconstructed_per_unit_depth, sentaurus_alpha_current_nodal, source_log_gap, validate_dependency_dag, interaction_dex, assert_counterfactual_closure, build_adjacent_interactions, symmetric_contributions, score_dominance, validate_field_units, validate_source_anchor_kind
from scripts.diagnose_pn2d_minimal6_formula_difference import _node_state_rows

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
    def test_interaction_and_dominance_gate_require_complete_matrix(self):
        forward = [{"factor":"gradient_recovery", "contribution_dex":0.8}, {"factor":"mobility", "contribution_dex":0.1}]
        reverse = [{"factor":"mobility", "contribution_dex":0.1}, {"factor":"gradient_recovery", "contribution_dex":0.2}]
        source = {frozenset():1., frozenset({"gradient_recovery"}):10., frozenset({"mobility"}):2., frozenset({"gradient_recovery", "mobility"}):30.}
        interactions = build_adjacent_interactions(forward, reverse, lambda replaced: source[frozenset(replaced)])
        self.assertEqual(len(interactions), 1)
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
    def test_cli_writes_source_family_ledger_and_named_residual(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp) / "state"
            root.mkdir()
            states = []
            for topology in ("sketch", "mirror"):
                for bias in (0.0, -12.0, -19.0):
                    export = root / f"{topology}_{bias}"
                    fields = export / "fields"
                    fields.mkdir(parents=True)
                    mesh = {"nodes":[{"id":0,"x":0.,"y":0.},{"id":1,"x":1.e-6,"y":0.},{"id":2,"x":0.,"y":1.e-6}], "triangles":[{"node_ids":[0,1,2]}]}
                    (export / "mesh.json").write_text(json.dumps(mesh), encoding="utf-8")
                    (fields / "ImpactIonization_region0.csv").write_text("node_id,component0\n0,2\n1,2\n2,2\n", encoding="utf-8")
                    field_rows = {
                        "eAlphaAvalanche":"node_id,component0\n0,1\n1,1\n2,1\n",
                        "hAlphaAvalanche":"node_id,component0\n0,0\n1,0\n2,0\n",
                        "eCurrentDensity":"node_id,component0,component1\n0,1,0\n1,1,0\n2,1,0\n",
                        "hCurrentDensity":"node_id,component0,component1\n0,0,0\n1,0,0\n2,0,0\n",
                        "ElectrostaticPotential":"node_id,component0\n0,0.1\n1,0.2\n2,0.3\n",
                        "eDensity":"node_id,component0\n0,1e10\n1,2e10\n2,3e10\n",
                        "hDensity":"node_id,component0\n0,4e10\n1,5e10\n2,6e10\n",
                        "eQuasiFermiPotential":"node_id,component0\n0,0.4\n1,0.5\n2,0.6\n",
                        "hQuasiFermiPotential":"node_id,component0\n0,0.7\n1,0.8\n2,0.9\n",
                        "eMobility":"node_id,component0\n0,100\n1,101\n2,102\n",
                        "hMobility":"node_id,component0\n0,50\n1,51\n2,52\n",
                        "eVelocity":"node_id,component0\n0,1000\n1,1001\n2,1002\n",
                        "hVelocity":"node_id,component0\n0,2000\n1,2001\n2,2002\n",
                        "LatticeTemperature":"node_id,component0\n0,300\n1,300\n2,300\n",
                    }
                    for name, content in field_rows.items():
                        (fields / f"{name}_region0.csv").write_text(content, encoding="utf-8")
                    audit_header = "cell_id,node0,node1,node2,grad_psi_x_V_per_m,grad_psi_y_V_per_m,grad_phin_x_V_per_m,grad_phin_y_V_per_m,grad_phip_x_V_per_m,grad_phip_y_V_per_m,local_edge0_edge_id,local_edge0_node0,local_edge0_node1,local_edge0_electron_cell_qf_field_V_per_m,local_edge0_hole_cell_qf_field_V_per_m,local_edge0_electron_midpoint_density_m3,local_edge0_hole_midpoint_density_m3,local_edge0_electron_mobility_m2_per_V_s,local_edge0_hole_mobility_m2_per_V_s,local_edge0_electron_alpha_per_m,local_edge0_hole_alpha_per_m,local_edge0_electron_flux_proxy_per_m2_s,local_edge0_hole_flux_proxy_per_m2_s,local_edge0_electron_source_integral_per_m_s,local_edge0_hole_source_integral_per_m_s"
                    audit_row = "0,0,1,2,1,2,3,4,5,6,9,0,1,7,8,9,10,11,12,13,14,15,16,2,3"
                    (export / "vela_triangle_audit.csv").write_text(audit_header + "\n" + audit_row + "\n", encoding="utf-8")
                    states.append({"topology_id":topology,"requested_bias_V":bias,"actual_bias_V":bias,"status":"passed","export_dir":str(export)})
            model_source = root / "source"
            model_source.mkdir()
            (model_source / "models.par").write_text(
                "vanOverstraetendeMan * Impact Ionization {\n"
                "  a(low) = 1e6, 2e6\n  a(high) = 3e6, 4e6\n"
                "  b(low) = 2e5, 3e5\n  b(high) = 4e5, 5e5\n"
                "  E0 = 4e5, 4e5\n  hbarOmega = 0.063, 0.063\n}\n",
                encoding="utf-8")
            (root / "manifest.json").write_text(json.dumps({"outputs_complete":True,"states":states}), encoding="utf-8")
            out = Path(temp) / "out"
            command = [sys.executable, str(Path(__file__).parents[2] / "scripts" / "diagnose_pn2d_minimal6_formula_difference.py"), "--state-root", str(root), "--audit-root", temp, "--out-dir", str(out), "--qa-status", "reviewed"]
            completed = subprocess.run(command, capture_output=True, text=True)
            self.assertEqual(completed.returncode, 0, completed.stderr)
            with (out / "quantity_ledger.csv").open(newline="", encoding="utf-8") as handle:
                ledger = list(csv.DictReader(handle))
            with (out / "factor_waterfall.csv").open(newline="", encoding="utf-8") as handle:
                waterfall = list(csv.DictReader(handle))
            report = json.loads((out / "root_cause_summary.json").read_text(encoding="utf-8"))
            source_rows = [row for row in ledger if row["record_kind"] == "source_integral"]
            node_rows = [row for row in ledger if row["record_kind"] == "node_state"]
            self.assertEqual(len(source_rows), 18)
            self.assertEqual({row["source"] for row in source_rows}, {"sentaurus_native_avalanche_generation", "sentaurus_alpha_current_reconstruction", "vela_alpha_flux_partial_volume_reconstruction"})
            self.assertEqual(len(node_rows), 270)
            self.assertEqual({row["unit"] for row in node_rows if row["quantity"] == "eMobility"}, {"cm^2*V^-1*s^-1"})
            self.assertTrue(any(row["quantity"] == "ni_eff_electron" for row in node_rows))
            self.assertTrue(any(row["quantity"] == "ni_eff_relative_residual" for row in node_rows))
            self.assertTrue(any(row["record_kind"] == "cell_replay" for row in ledger))
            self.assertTrue(any(row["record_kind"] == "edge_replay" for row in ledger))
            self.assertTrue(any(row["quantity"] == "sentaurus_dem_electron_alpha_recomputed" and row["source"] == "models_par_reference" for row in ledger))
            self.assertEqual(len(waterfall), 54)
            self.assertEqual(waterfall[0]["factor"], "ni_eff/BGN")
            self.assertEqual(waterfall[-1]["factor"], "unattributed_residual")
            self.assertEqual(len(report["sentaurus_internal_semantics_residual"]), 6)
            self.assertEqual(len(report["waterfall_paths"]), 6)
            self.assertEqual(report["vela_parameter_agreement"][0]["status"], "unavailable")
            self.assertEqual(report["waterfall_paths"][0]["factor_availability"][0]["status"], "unavailable")
            self.assertEqual(report["dominance_rules"]["status"], "insufficient_data")
            self.assertIn("vela_native_minus_reconstruction", report["records"][0])
            summary_markdown = (out / "root_cause_summary.md").read_text(encoding="utf-8")
            self.assertIn("reference_tcad/pn2d_sentaurus2018_minimal6/source/models.par", summary_markdown)
            self.assertIn("scripts/pn2d_minimal6_diagnostics/physics.py", summary_markdown)
            self.assertIn("src/physics/ImpactIonizationModel.cpp", summary_markdown)
            figure_manifest_path = out / "figure_manifest.json"
            self.assertTrue(figure_manifest_path.is_file(), figure_manifest_path)
            figure_manifest = json.loads(figure_manifest_path.read_text(encoding="utf-8"))
            self.assertEqual(
                figure_manifest["diagnostic_disclaimer"],
                "minimal6 diagnostic sweep; not a physical BV curve",
            )
            self.assertEqual(figure_manifest["manual_qa"]["status"], "reviewed")
            self.assertEqual(
                [item["stem"] for item in figure_manifest["figures"]],
                ["gradient", "current_alpha", "source_waterfall", "interaction", "topology_symmetry"],
            )
            for item in figure_manifest["figures"]:
                self.assertIn("unit", item)
                for relative_path in item["artifacts"]:
                    path = out / relative_path
                    self.assertTrue(path.is_file(), path)
                    payload = path.read_bytes()
                    if path.suffix == ".png":
                        self.assertTrue(payload.startswith(b"\x89PNG\r\n\x1a\n"))
                        from PIL import Image
                        with Image.open(path) as image:
                            self.assertGreaterEqual(image.width, 640)
                            self.assertGreaterEqual(image.height, 360)
                            self.assertGreater(len(image.getcolors(maxcolors=image.width * image.height)), 1)
                    else:
                        self.assertTrue(payload.startswith(b"%PDF-"))
    def test_rejects_inexact_bias(self):
        states = [{"topology_id":t,"requested_bias_V":b,"actual_bias_V":b,"status":"passed"}
                  for t in ("sketch","mirror") for b in (0.0,-12.0,-19.0)]
        states[-1]["actual_bias_V"] += 2e-12
        with self.assertRaises(ValueError): validate_formula_input({"outputs_complete":True,"states":states})

if __name__ == '__main__': unittest.main()
