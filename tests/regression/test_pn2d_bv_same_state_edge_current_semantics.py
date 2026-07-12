from __future__ import annotations

import csv
import importlib.util
import json
import math
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT_PATH = REPO_ROOT / "scripts" / "audit_pn2d_bv_same_state_edge_current_semantics.py"


def _load_audit():
    if not SCRIPT_PATH.is_file():
        return None
    spec = importlib.util.spec_from_file_location("same_state_audit", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


audit = _load_audit()


class TestAuditModuleExists(unittest.TestCase):
    def test_audit_script_exists(self):
        self.assertTrue(SCRIPT_PATH.is_file(), f"missing audit script: {SCRIPT_PATH}")


@unittest.skipIf(audit is None, "audit script not implemented yet")
class TestSameStateCurrentSemantics(unittest.TestCase):
    def test_adjacent_cell_aggregation_and_area_gate(self):
        triangle = [
            {
                "bias_V": "-12", "edge_id": "50", "node0": "2", "node1": "1",
                "truncated_partial_volume_m2": "2", "electron_flux_proxy": "10",
                "hole_flux_proxy": "4", "electron_alpha_m_inv": "3",
                "hole_alpha_m_inv": "5", "electron_cell_qf_field_V_per_m": "7",
                "hole_cell_qf_field_V_per_m": "9", "electron_edge_qf_field_V_per_m": "11",
                "hole_edge_qf_field_V_per_m": "13", "electron_midpoint_density_m3": "17",
                "hole_midpoint_density_m3": "19", "electron_mobility_m2_V_s": "0.1",
                "hole_mobility_m2_V_s": "0.2",
            },
            {
                "bias_V": "-12", "edge_id": "50", "node0": "2", "node1": "1",
                "truncated_partial_volume_m2": "1", "electron_flux_proxy": "4",
                "hole_flux_proxy": "1", "electron_alpha_m_inv": "6",
                "hole_alpha_m_inv": "8", "electron_cell_qf_field_V_per_m": "10",
                "hole_cell_qf_field_V_per_m": "12", "electron_edge_qf_field_V_per_m": "11",
                "hole_edge_qf_field_V_per_m": "13", "electron_midpoint_density_m3": "23",
                "hole_midpoint_density_m3": "29", "electron_mobility_m2_V_s": "0.3",
                "hole_mobility_m2_V_s": "0.4",
            },
        ]
        aggregate = audit.aggregate_triangle_rows(triangle)[(-12.0, 50)]
        self.assertAlmostEqual(aggregate["partial_volume_sum_m2"], 3.0)
        self.assertAlmostEqual(aggregate["electron_pdf_gradqf_flux_m2_s"], 8.0)
        self.assertAlmostEqual(aggregate["hole_pdf_gradqf_flux_m2_s"], 3.0)
        self.assertAlmostEqual(aggregate["electron_qf_source_proxy"], 84.0)
        self.assertAlmostEqual(aggregate["hole_qf_source_proxy"], 48.0)
        gate = audit.enforce_shared_edge_area_gate(
            {(-12.0, 50): aggregate},
            [{"bias_V": "-12", "edge_id": "50", "edge_area_proxy_m2": "3"}],
        )
        self.assertEqual(gate["common_edge_count"], 1)
        self.assertLess(gate["max_relative_error"], 1.0e-12)

    def test_area_mismatch_fails(self):
        aggregate = {(-12.0, 50): {"partial_volume_sum_m2": 3.0}}
        with self.assertRaisesRegex(ValueError, "shared edge area mismatch"):
            audit.enforce_shared_edge_area_gate(
                aggregate,
                [{"bias_V": "-12", "edge_id": "50", "edge_area_proxy_m2": "3.01"}],
            )

    def test_canonical_reversal_and_electron_continuity_sign(self):
        forward = audit.project_endpoint_vector(
            (0.0, 0.0), (1.0, 0.0), (2.0, 3.0), (4.0, 5.0),
            electron_continuity=True,
        )
        reverse = audit.project_endpoint_vector(
            (1.0, 0.0), (0.0, 0.0), (4.0, 5.0), (2.0, 3.0),
            electron_continuity=True,
        )
        self.assertEqual(forward["canonical_projection"], -3.0)
        self.assertEqual(reverse["canonical_projection"], -3.0)
        self.assertAlmostEqual(forward["magnitude"], math.hypot(3.0, 4.0))

    def test_legacy_vtk_vector_parser(self):
        with tempfile.TemporaryDirectory() as tmp:
            vtk = Path(tmp) / "sample.vtk"
            vtk.write_text(
                "# vtk DataFile Version 3.0\nfixture\nASCII\nDATASET UNSTRUCTURED_GRID\n"
                "POINTS 3 double\n0 0 0\n1 0 0\n0 1 0\n"
                "CELLS 1 4\n3 0 1 2\nCELL_TYPES 1\n5\n"
                "CELL_DATA 1\nSCALARS region_id int 1\nLOOKUP_TABLE default\n7\n"
                "POINT_DATA 3\nSCALARS Potential double 1\nLOOKUP_TABLE default\n1\n2\n3\n"
                "VECTORS ElectronCurrentDensityVector double\n1 2 0\n3 4 0\n5 6 0\n",
                encoding="ascii",
            )
            parsed = audit.parse_legacy_ascii_vtk(vtk)
        self.assertEqual(parsed["points"][1], (1.0, 0.0, 0.0))
        self.assertEqual(parsed["cells"], [(0, 1, 2)])
        self.assertEqual(parsed["cell_data"]["region_id"], [7.0])
        self.assertEqual(
            parsed["point_data"]["ElectronCurrentDensityVector"][2],
            (5.0, 6.0, 0.0),
        )

    def test_nearest_sentaurus_bias_is_not_in_exact_summary(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "sentaurus_-19v").mkdir()
            (root / "sentaurus_-20v").mkdir()
            selected = audit.select_sentaurus_export(root, -19.4)
        self.assertEqual(selected["selected_bias_V"], -19.0)
        self.assertFalse(selected["exact_match"])
        rows = [
            {"exact_match": True, "sentaurus_vector_magnitude_flux_m2_s": 10.0,
             "pdf_log10_abs_error": 1.0},
            {"exact_match": True, "sentaurus_vector_magnitude_flux_m2_s": 20.0,
             "pdf_log10_abs_error": 2.0},
            {"exact_match": False, "sentaurus_vector_magnitude_flux_m2_s": 1000.0,
             "pdf_log10_abs_error": 9.0},
        ]
        summary = audit.summarize_active_support(rows, ["pdf"])
        self.assertEqual(summary["exact_row_count"], 2)
        self.assertAlmostEqual(summary["positive_p80_threshold"], 18.0)
        self.assertEqual(summary["candidates"]["pdf"]["median_abs_log10_error"], 2.0)

    def test_active_support_median_and_p95(self):
        rows = []
        for magnitude, error in zip([1, 2, 3, 4, 5, 6, 7, 8, 9, 10], range(10)):
            rows.append({
                "exact_match": True,
                "sentaurus_vector_magnitude_flux_m2_s": float(magnitude),
                "pdf_log10_abs_error": float(error),
            })
        summary = audit.summarize_active_support(rows, ["pdf"])
        self.assertAlmostEqual(summary["positive_p80_threshold"], 8.2)
        self.assertEqual(summary["active_row_count"], 2)
        self.assertAlmostEqual(summary["candidates"]["pdf"]["median_abs_log10_error"], 8.5)
        self.assertAlmostEqual(summary["candidates"]["pdf"]["p95_abs_log10_error"], 8.95)

    def test_contact_gate_true_and_false(self):
        favorable = [
            {"exact_match": True, "edge_class": "contact_edge", "active_support": True,
             "qf_log10_abs_error": 1.0, "fallback_log10_abs_error": 0.5},
            {"exact_match": True, "edge_class": "contact_edge", "active_support": True,
             "qf_log10_abs_error": 0.9, "fallback_log10_abs_error": 0.4},
            {"exact_match": True, "edge_class": "interior_bulk", "active_support": True,
             "qf_log10_abs_error": 0.2, "fallback_log10_abs_error": 0.4},
        ]
        result = audit.evaluate_contact_policy_gate(favorable)
        self.assertEqual(result["improvement_coverage"], 1.0)
        self.assertGreater(result["interior_median_worsening_dex"], 0.0)
        self.assertTrue(result["recommend_explicit_contact_policy"])
        unfavorable = [dict(row) for row in favorable]
        unfavorable[0]["fallback_log10_abs_error"] = 0.9
        unfavorable[1]["fallback_log10_abs_error"] = 0.8
        self.assertFalse(
            audit.evaluate_contact_policy_gate(unfavorable)["recommend_explicit_contact_policy"]
        )

    def test_cli_writes_csv_json_and_markdown(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            triangle = root / "triangle.csv"
            sg = root / "sg.csv"
            vtk_root = root / "vtk"
            sent_root = root / "sentaurus"
            out = root / "out"
            vtk_root.mkdir()
            self._write_triangle(triangle)
            self._write_sg(sg)
            self._write_vtk(vtk_root / "dc_sweep_0001_-12V.vtk")
            self._write_sentaurus(sent_root / "sentaurus_-12v")
            completed = subprocess.run(
                [sys.executable, str(SCRIPT_PATH), "--triangle-csv", str(triangle),
                 "--sg-edge-csv", str(sg), "--vtk-root", str(vtk_root),
                 "--sentaurus-root", str(sent_root), "--out-dir", str(out),
                 "--biases=-12", "--focus-edge", "50", "--top-n", "3"],
                text=True, capture_output=True, check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)
            paths = [
                out / "same_state_edge_current_semantics.csv",
                out / "same_state_edge_current_semantics.json",
                out / "same_state_edge_current_semantics.md",
            ]
            self.assertTrue(all(path.is_file() for path in paths))
            payload = json.loads(paths[1].read_text(encoding="utf-8"))
            self.assertEqual(payload["row_count"], 1)
            with paths[0].open(newline="", encoding="utf-8") as handle:
                row = next(csv.DictReader(handle))
            required = {
                "bias_V", "edge_id", "edge_class", "contact", "source_rank",
                "endpoint_n0_m3", "endpoint_p0_m3", "endpoint_psi0_V",
                "endpoint_phin0_V", "endpoint_phip0_V", "electron_gss_midpoint_density_m3",
                "electron_mobility_m2_V_s", "electron_cell_qf_field_V_per_m",
                "electron_edge_qf_field_V_per_m", "electric_field_V_per_m",
                "electron_qf_alpha_m_inv", "electron_electric_fallback_alpha_m_inv",
                "electron_pdf_gradqf_flux_m2_s", "electron_genius_sg_flux_m2_s",
                "electron_vela_vector_projection_flux_m2_s",
                "electron_sentaurus_vector_projection_flux_m2_s",
                "electron_pdf_log10_abs_error", "partial_volume_sum_m2",
                "edge_area_proxy_m2", "area_relative_error",
                "electron_pdf_qf_source_proxy", "electron_fallback_source_proxy",
            }
            self.assertTrue(required.issubset(row), sorted(required - set(row)))
            self.assertNotEqual(row["electron_sentaurus_source_proxy"], "")
            self.assertIn("Data Contract Issues", paths[2].read_text(encoding="utf-8"))

    @staticmethod
    def _write_triangle(path: Path):
        fields = [
            "bias_V", "edge_id", "node0", "node1", "x0_um", "y0_um", "x1_um", "y1_um",
            "truncated_partial_volume_m2", "electron_flux_proxy", "hole_flux_proxy",
            "electron_alpha_m_inv", "hole_alpha_m_inv", "electron_cell_qf_field_V_per_m",
            "hole_cell_qf_field_V_per_m", "electron_edge_qf_field_V_per_m",
            "hole_edge_qf_field_V_per_m", "electron_midpoint_density_m3",
            "hole_midpoint_density_m3", "electron_mobility_m2_V_s", "hole_mobility_m2_V_s",
        ]
        rows = []
        for volume in (1.0, 2.0):
            rows.append(dict(zip(fields, [
                -12, 50, 0, 1, 0, 0, 1, 0, volume, 10, 5, 2, 3, 4, 5, 6, 7,
                8, 9, 0.1, 0.2,
            ])))
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)

    @staticmethod
    def _write_sg(path: Path):
        row = {
            "bias_V": -12, "edge_id": 50, "node0": 0, "node1": 1,
            "x0_um": 0, "y0_um": 0, "x1_um": 1, "y1_um": 0,
            "edge_area_proxy_m2": 3, "edge_class": "contact_edge",
            "electric_field_V_per_m": 1.0e8,
            "electron_sg_n0": 11, "electron_sg_n1": 12,
            "electron_sg_psi0": 1, "electron_sg_psi1": 2,
            "electron_sg_phin0": 3, "electron_sg_phin1": 4,
            "electron_sg_production_abs_continuity_particle_flux_m2_s": 20,
            "hole_raw_flux_proxy": 7,
        }
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(row))
            writer.writeheader()
            writer.writerow(row)

    @staticmethod
    def _write_vtk(path: Path):
        path.write_text(
            "# vtk DataFile Version 3.0\nfixture\nASCII\nDATASET UNSTRUCTURED_GRID\n"
            "POINTS 2 double\n0 0 0\n1 0 0\nCELLS 0 0\nCELL_TYPES 0\nPOINT_DATA 2\n"
            "VECTORS ElectronCurrentDensityVector double\n-1 0 0\n-1 0 0\n"
            "VECTORS HoleCurrentDensityVector double\n2 0 0\n2 0 0\n",
            encoding="ascii",
        )

    @staticmethod
    def _write_sentaurus(path: Path):
        fields = path / "fields"
        fields.mkdir(parents=True)
        (path / "nodes.csv").write_text("id,x_um,y_um\n0,0,0\n1,1,0\n", encoding="utf-8")
        (fields / "eCurrentDensity_region0.csv").write_text(
            "node_id,component0,component1\n0,-1,0\n1,-1,0\n", encoding="utf-8"
        )
        (fields / "hCurrentDensity_region0.csv").write_text(
            "node_id,component0,component1\n0,2,0\n1,2,0\n", encoding="utf-8"
        )
        (fields / "eAlphaAvalanche_region0.csv").write_text(
            "node_id,component0\n0,3\n1,3\n", encoding="utf-8"
        )
        (fields / "hAlphaAvalanche_region0.csv").write_text(
            "node_id,component0\n0,4\n1,4\n", encoding="utf-8"
        )


if __name__ == "__main__":
    unittest.main()
