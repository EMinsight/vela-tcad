from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path
import shutil
import tempfile
import unittest

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "audit_pn2d_minimal6_fixed_state.py"
FIXTURE = REPO / "tests" / "fixtures" / "pn2d_minimal6_synthetic"
spec = importlib.util.spec_from_file_location("pn2d_minimal6_fixed_state_audit", SCRIPT)
if spec is None or spec.loader is None:
    raise RuntimeError(f"cannot load audit module from {SCRIPT}")
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


class FormulaReferenceTests(unittest.TestCase):
    def test_linear_triangle_gradient_matches_closed_form(self):
        self.assertEqual(audit.triangle_gradient([(0, 0), (1, 0), (0, 1)], [2, 5, 7]), (3.0, 5.0))

    def test_stable_bernoulli_limits(self):
        self.assertEqual(audit.bernoulli(0.0), 1.0)
        self.assertAlmostEqual(audit.bernoulli(1.0e-10), 1.0 - 0.5e-10)
        self.assertEqual(audit.bernoulli(-1000.0), 1000.0)
        self.assertEqual(audit.bernoulli(1000.0), 0.0)

    def test_electron_and_hole_sg_flux_signs_match_production(self):
        electron = audit.sg_electron_flux(2.0, 5.0, 0.0, 0.025, 3.0, 2.0)
        hole = audit.sg_hole_flux(2.0, 5.0, 0.0, 0.025, 3.0, 2.0)
        self.assertAlmostEqual(electron, 0.1125)
        self.assertAlmostEqual(hole, 0.1125)

    def test_gss_logistic_midpoint_uses_carrier_orientation(self):
        vt = 0.025852
        self.assertGreater(audit.gss_logistic_midpoint(1e12, 8e17, -0.2, 0.1, vt, "electron"), 7.9e17)
        self.assertGreater(audit.gss_logistic_midpoint(3e17, 2e11, -0.2, 0.1, vt, "hole"), 2.9e17)

    def test_canonical_projection_uses_endpoint_direction(self):
        self.assertAlmostEqual(audit.canonical_projection((3.0, 4.0), (0.0, 0.0), (2.0, 0.0)), 3.0)
        self.assertAlmostEqual(audit.canonical_projection((3.0, 4.0), (2.0, 0.0), (0.0, 0.0)), -3.0)

    def test_genius_truncated_partial_volume_for_right_triangle(self):
        points = [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)]
        self.assertAlmostEqual(audit.genius_truncated_partial_volume(points, 0), 0.25)
        self.assertAlmostEqual(audit.genius_truncated_partial_volume(points, 1), 0.0)
        self.assertAlmostEqual(audit.genius_truncated_partial_volume(points, 2), 0.25)

    def test_hybrid_error_and_both_zero_classification(self):
        self.assertEqual(audit.hybrid_error(0.0, 0.0), 0.0)
        self.assertEqual(audit.hybrid_error(2.0, 1.0), 0.5)
        ratio = audit.classify_orientation_pair(0.0, 0.0)
        self.assertEqual(ratio["zero_classification"], "both_zero")
        self.assertIsNone(ratio["absolute_log10_ratio"])


class ReviewContractTests(unittest.TestCase):
    def _copy(self):
        temp = tempfile.TemporaryDirectory()
        root = Path(temp.name) / "fixture"
        shutil.copytree(FIXTURE, root)
        self.addCleanup(temp.cleanup)
        return root

    @staticmethod
    def _manifest(root):
        return json.loads((root / "manifest.json").read_text(encoding="utf-8"))

    @staticmethod
    def _write(path, value):
        path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")

    @staticmethod
    def _export(root, state):
        path = Path(state["export_dir"])
        return path if path.is_absolute() else root / path

    @staticmethod
    def _mutate_csv(path, column, value, row=0):
        with path.open(encoding="utf-8", newline="") as source:
            rows = list(csv.DictReader(source))
        rows[row][column] = value
        with path.open("w", encoding="utf-8", newline="") as output:
            writer = csv.DictWriter(output, fieldnames=rows[0])
            writer.writeheader(); writer.writerows(rows)

    def test_exact_variable_ni_sg_both_carriers(self):
        vt = 1.380649e-23 * 300.0 / 1.602176634e-19
        coef = 0.135 * vt / 1.25e-6
        args = (1.45e16, 2.2e16, -0.03, 0.07)
        self.assertAlmostEqual(audit.sg_electron_variable_ni_flux(*args, -0.021, 0.013, vt, coef), 9.085130342836316e19, delta=1e5)
        self.assertAlmostEqual(audit.sg_hole_variable_ni_flux(*args, 0.025, 0.082, vt, coef), -3.0985907840586436e20, delta=1e6)

    def test_actual_task3_root_and_field_manifest_list(self):
        manifest = self._manifest(FIXTURE)
        self.assertEqual(manifest["schema"], "vela.pn2d_minimal6_states.v1")
        self.assertEqual(len(manifest["states"]), 6)
        export = self._export(FIXTURE, manifest["states"][0])
        fields = json.loads((export / "field_manifest.json").read_text(encoding="utf-8"))["fields"]
        self.assertIsInstance(fields, list)
        report = audit.build_report(FIXTURE)
        self.assertEqual((len(report.node_rows), len(report.edge_rows), len(report.triangle_rows)), (36,54,24))

    def test_requested_actual_and_optional_field_bias_mismatch(self):
        root = self._copy(); manifest = self._manifest(root)
        manifest["states"][0]["actual_bias_V"] = 1e-15
        self._write(root / "manifest.json", manifest)
        with self.assertRaisesRegex(audit.ContractError, "requested.*actual bias"):
            audit.build_report(root)
        root = self._copy(); manifest = self._manifest(root); state = manifest["states"][0]
        field_path = self._export(root, state) / "field_manifest.json"
        field = json.loads(field_path.read_text(encoding="utf-8")); field["bias_V"] = -1.0
        self._write(field_path, field)
        with self.assertRaisesRegex(audit.ContractError, "field manifest bias"):
            audit.build_report(root)

    def test_noncompensated_doping_and_local_edge_mapping_fail(self):
        root = self._copy(); state = self._manifest(root)["states"][0]; export = self._export(root,state)
        self._mutate_csv(export / "doping.csv", "donors_cm3", "1e17", 0)
        with self.assertRaisesRegex(audit.ContractError, "doping semantics"):
            audit.build_report(root)
        root = self._copy(); state = self._manifest(root)["states"][0]; export = self._export(root,state)
        self._mutate_csv(export / "vela_triangle_audit.csv", "local_edge0_node0", "5", 0)
        with self.assertRaisesRegex(audit.ContractError, "local edge"):
            audit.build_report(root)

    def test_wide_schema_orientation_and_abs_log_ratio(self):
        self.assertAlmostEqual(audit.classify_orientation_pair(10.0,1.0)["absolute_log10_ratio"], 1.0)
        report = audit.build_report(FIXTURE)
        for key in ("raw_eDensity_cm3","sentaurus_n_m3","vela_n_m3","abs_error_n_m3","hybrid_error_n_m3"):
            self.assertIn(key, report.node_rows[0])
        for key in ("node0_psi_V","electron_mobility_m2_per_V_s","python_electron_flux_per_m2_s","vela_electron_current_A_per_m2","electron_formula_hybrid_error"):
            self.assertIn(key, report.edge_rows[0])
        for key in ("shape_grad_N0_x_per_m","vela_grad_psi_x_V_per_m","reconstructed_electron_current_x_A_per_m2","sentaurus_vs_vela_total_source_diagnostic"):
            self.assertIn(key, report.triangle_rows[0])
        self.assertGreater(len({x["quantity"] for x in report.summary["orientation_sensitivity"]}), 51)

    def test_cpp_provenance_replay_and_no_python_vela_producer(self):
        manifest = self._manifest(FIXTURE); provenance = manifest["task4_provenance"]
        self.assertEqual(provenance["producer"], "build-release/pn2d_minimal6_operator_audit.exe")
        self.assertEqual(provenance["task4_source_commit"], "37a95459dc5f360bb24b9afa00439301935e98de")
        self.assertEqual(len(provenance["replays"]), 6)
        if (REPO / "build-release" / "pn2d_minimal6_operator_audit.exe").exists():
            self.assertEqual(audit.verify_task4_replay(FIXTURE, REPO / "build-release" / "pn2d_minimal6_operator_audit.exe"), [])
        source = SCRIPT.read_text(encoding="utf-8")
        self.assertNotIn("make_synthetic_fixture", source)
        self.assertNotIn("--make-synthetic-fixture", source)

    def test_fourteen_figures_and_numeric_summary(self):
        report = audit.build_report(FIXTURE)
        with tempfile.TemporaryDirectory() as directory:
            out = Path(directory); audit.write_report(report,out)
            figures = list((out / "figures").iterdir())
            self.assertEqual(len(figures),14); self.assertTrue(all(x.stat().st_size for x in figures))
            md=(out/"summary.md").read_text(encoding="utf-8")
            self.assertIn("Maximum state parity hybrid error:",md)
            self.assertIn("Maximum C++/Python formula hybrid error:",md)
    def test_complete_matrix_unique_keys_and_formula_gate(self):
        report = audit.build_report(FIXTURE)
        self.assertEqual((len(report.node_rows),len(report.edge_rows),len(report.triangle_rows)),(36,54,24))
        audit.require_unique_keys(report.edge_rows,("topology_id","bias_V","node0","node1"))
        self.assertLess(report.summary["gates"]["max_cpp_python_formula_hybrid_error"],5e-12)
        self.assertIsNone(report.summary["gates"]["sentaurus_vs_vela_current_source_threshold"])

    def test_partial_state_fails_closed(self):
        root=self._copy(); manifest=self._manifest(root); state=manifest["states"][0]
        path=root/state["state_csv"]; lines=path.read_text(encoding="utf-8").splitlines(); path.write_text("\n".join(lines[:-1])+"\n",encoding="utf-8")
        with self.assertRaisesRegex(audit.ContractError,"partial state"): audit.build_report(root)

    def test_missing_required_field_fails_closed(self):
        root=self._copy(); state=self._manifest(root)["states"][0]; path=self._export(root,state)/"field_manifest.json"
        value=json.loads(path.read_text(encoding="utf-8")); value["fields"]=[x for x in value["fields"] if x["name"]!="eMobility"]; self._write(path,value)
        with self.assertRaisesRegex(audit.ContractError,"missing required field"): audit.build_report(root)

    def test_duplicate_edge_key_fails_closed(self):
        root=self._copy(); state=self._manifest(root)["states"][0]; path=self._export(root,state)/"vela_edge_audit.csv"
        lines=path.read_text(encoding="utf-8").splitlines(); path.write_text("\n".join(lines+[lines[1]])+"\n",encoding="utf-8")
        with self.assertRaisesRegex(audit.ContractError,"row count|duplicate"): audit.build_report(root)

    def test_inexact_requested_bias_fails_closed(self):
        root=self._copy(); manifest=self._manifest(root); manifest["states"][0]["requested_bias_V"]=-11.999999999; manifest["states"][0]["actual_bias_V"]=-11.999999999; self._write(root/"manifest.json",manifest)
        with self.assertRaisesRegex(audit.ContractError,"unexpected topology, bias"): audit.build_report(root)

    def test_wrong_vector_components_and_units_fail_closed(self):
        for key,value,pattern in (("components",1,"missing required field"),("unit","V/cm","wrong unit")):
            root=self._copy(); state=self._manifest(root)["states"][0]; path=self._export(root,state)/"field_manifest.json"; manifest=json.loads(path.read_text(encoding="utf-8"))
            field=next(x for x in manifest["fields"] if x["name"]=="ElectricField" and x["components"]==2); field[key]=value; self._write(path,manifest)
            with self.subTest(key=key),self.assertRaisesRegex(audit.ContractError,pattern): audit.build_report(root)

    def test_changed_topology_fails_closed(self):
        root=self._copy(); state=self._manifest(root)["states"][0]; path=self._export(root,state)/"elements.csv"; self._mutate_csv(path,"node2","5",0)
        with self.assertRaisesRegex(audit.ContractError,"topology|triangle orientation"): audit.build_report(root)

    def test_nonfinite_state_fails_closed(self):
        root=self._copy(); state=self._manifest(root)["states"][0]; self._mutate_csv(root/state["state_csv"],"psi_V","nan",0)
        with self.assertRaisesRegex(audit.ContractError,"non-finite"): audit.build_report(root)

    def test_formula_error_above_threshold_fails_closed(self):
        root=self._copy(); state=self._manifest(root)["states"][0]; path=self._export(root,state)/"vela_edge_audit.csv"
        with path.open(encoding="utf-8",newline="") as source: rows=list(csv.DictReader(source))
        rows[0]["electron_raw_signed_flux_per_m2_s"]=format(float(rows[0]["electron_raw_signed_flux_per_m2_s"])*(1+5.1e-12),".17g")
        with path.open("w",encoding="utf-8",newline="") as output: writer=csv.DictWriter(output,fieldnames=rows[0]); writer.writeheader(); writer.writerows(rows)
        with self.assertRaisesRegex(audit.ContractError,"formula error"): audit.build_report(root)

    def test_state_error_above_threshold_fails_closed(self):
        root=self._copy(); state=self._manifest(root)["states"][0]; path=self._export(root,state)/"vela_node_state.csv"
        with path.open(encoding="utf-8",newline="") as source: rows=list(csv.DictReader(source))
        rows[0]["n_m3"]=format(float(rows[0]["n_m3"])*(1+1.1e-12),".17g")
        with path.open("w",encoding="utf-8",newline="") as output: writer=csv.DictWriter(output,fieldnames=rows[0]); writer.writeheader(); writer.writerows(rows)
        with self.assertRaisesRegex(audit.ContractError,"state parity"): audit.build_report(root)

    def test_raw_unit_conversions_and_normalized_zero_are_visible(self):
        report=audit.build_report(FIXTURE); node=report.node_rows[0]
        self.assertEqual(node["sentaurus_n_m3"],node["raw_eDensity_cm3"]*1e6)
        self.assertEqual(node["sentaurus_electric_field_x_V_per_m"],node["raw_ElectricField_x_V_per_cm"]*100)
        self.assertTrue(any(row["normalized_geometric_zero"] for row in report.triangle_rows))


if __name__ == "__main__":
    unittest.main()