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


class FullMatrixAndContractTests(unittest.TestCase):
    def _copy_fixture(self):
        temp = tempfile.TemporaryDirectory()
        destination = Path(temp.name) / "fixture"
        shutil.copytree(FIXTURE, destination)
        return temp, destination

    @staticmethod
    def _manifest(path):
        return json.loads((path / "manifest.json").read_text(encoding="utf-8"))

    @staticmethod
    def _write_manifest(path, manifest):
        (path / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    @staticmethod
    def _first_state(manifest):
        return manifest["topologies"][0]["states"][0]

    def test_complete_matrix_has_exact_row_counts_and_unique_keys(self):
        report = audit.build_report(FIXTURE)
        self.assertEqual(len(report.node_rows), 36)
        self.assertEqual(len(report.edge_rows), 54)
        self.assertEqual(len(report.triangle_rows), 24)
        audit.require_unique_keys(report.edge_rows, ("topology_id", "bias_V", "node0", "node1"))
        self.assertEqual(report.summary["schema"], "vela.pn2d_minimal6_fixed_state_audit.v1")
        self.assertTrue(report.summary["gates"]["passed"])
        self.assertIsNone(report.summary["gates"]["sentaurus_vs_vela_current_source_threshold"])

    def _mutate_csv(self, root, key, column, transform):
        manifest = self._manifest(root)
        csv_path = root / self._first_state(manifest)[key]
        with csv_path.open(encoding="utf-8", newline="") as source:
            rows = list(csv.DictReader(source))
        rows[0][column] = transform(rows[0][column])
        with csv_path.open("w", encoding="utf-8", newline="") as output:
            writer = csv.DictWriter(output, fieldnames=rows[0].keys())
            writer.writeheader(); writer.writerows(rows)

    def test_partial_state_fails_closed(self):
        temp, root = self._copy_fixture(); self.addCleanup(temp.cleanup)
        manifest = self._manifest(root); csv_path = root / self._first_state(manifest)["sentaurus_state_csv"]
        lines = csv_path.read_text(encoding="utf-8").splitlines()
        csv_path.write_text("\n".join(lines[:-1]) + "\n", encoding="utf-8")
        with self.assertRaisesRegex(audit.ContractError, "partial state"): audit.build_report(root)
    def test_missing_required_field_fails_closed(self):
        temp, root = self._copy_fixture(); self.addCleanup(temp.cleanup)
        manifest = self._manifest(root)
        field_path = root / self._first_state(manifest)["field_manifest"]
        fields = json.loads(field_path.read_text(encoding="utf-8")); del fields["fields"]["eMobility"]
        field_path.write_text(json.dumps(fields, indent=2), encoding="utf-8")
        with self.assertRaisesRegex(audit.ContractError, "missing required field"): audit.build_report(root)

    def test_duplicate_key_fails_closed(self):
        temp, root = self._copy_fixture(); self.addCleanup(temp.cleanup)
        manifest = self._manifest(root); csv_path = root / self._first_state(manifest)["vela_edge_csv"]
        lines = csv_path.read_text(encoding="utf-8").splitlines(); csv_path.write_text("\n".join(lines + [lines[1]]) + "\n", encoding="utf-8")
        with self.assertRaisesRegex(audit.ContractError, "duplicate"): audit.build_report(root)

    def test_inexact_bias_fails_closed(self):
        temp, root = self._copy_fixture(); self.addCleanup(temp.cleanup)
        manifest = self._manifest(root); self._first_state(manifest)["bias_V"] = -11.999999999; self._write_manifest(root, manifest)
        with self.assertRaisesRegex(audit.ContractError, "exact biases"): audit.build_report(root)

    def test_wrong_vector_component_count_fails_closed(self):
        temp, root = self._copy_fixture(); self.addCleanup(temp.cleanup)
        manifest = self._manifest(root); field_path = root / self._first_state(manifest)["field_manifest"]
        fields = json.loads(field_path.read_text(encoding="utf-8")); fields["fields"]["ElectricField"]["components"] = 1
        field_path.write_text(json.dumps(fields, indent=2), encoding="utf-8")
        with self.assertRaisesRegex(audit.ContractError, "component count"): audit.build_report(root)
    def test_wrong_unit_fails_closed(self):
        temp, root = self._copy_fixture(); self.addCleanup(temp.cleanup)
        manifest = self._manifest(root); field_path = root / self._first_state(manifest)["field_manifest"]
        fields = json.loads(field_path.read_text(encoding="utf-8")); fields["fields"]["ElectricField"]["unit"] = "V/cm"
        field_path.write_text(json.dumps(fields, indent=2), encoding="utf-8")
        with self.assertRaisesRegex(audit.ContractError, "unit"): audit.build_report(root)

    def test_wrong_topology_fails_closed(self):
        temp, root = self._copy_fixture(); self.addCleanup(temp.cleanup)
        manifest = self._manifest(root); topology_path = root / manifest["topologies"][0]["topology_file"]
        topology = json.loads(topology_path.read_text(encoding="utf-8")); topology["triangles"][0] = [1, 5, 6]
        topology_path.write_text(json.dumps(topology, indent=2), encoding="utf-8")
        with self.assertRaisesRegex(audit.ContractError, "topology"): audit.build_report(root)

    def test_nonfinite_value_fails_closed(self):
        temp, root = self._copy_fixture(); self.addCleanup(temp.cleanup)
        self._mutate_csv(root, "sentaurus_state_csv", "psi_V", lambda _: "nan")
        with self.assertRaisesRegex(audit.ContractError, "non-finite"): audit.build_report(root)

    def test_formula_error_at_threshold_fails_closed(self):
        temp, root = self._copy_fixture(); self.addCleanup(temp.cleanup)
        self._mutate_csv(root, "vela_edge_csv", "electron_raw_signed_flux_per_m2_s", lambda x: format(float(x) * (1 + 5e-12), ".17g"))
        with self.assertRaisesRegex(audit.ContractError, "formula error"): audit.build_report(root)

    def test_state_error_at_threshold_fails_closed(self):
        temp, root = self._copy_fixture(); self.addCleanup(temp.cleanup)
        self._mutate_csv(root, "vela_node_csv", "n_m3", lambda x: format(float(x) * (1 + 1e-12), ".17g"))
        with self.assertRaisesRegex(audit.ContractError, "state parity"): audit.build_report(root)


if __name__ == "__main__":
    unittest.main()