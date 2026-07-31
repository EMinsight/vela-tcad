import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = (
    Path(__file__).resolve().parents[2]
    / "scripts"
    / "prepare_pn2d_bv_off_srh_mesh_matrix.py"
)
SPEC = importlib.util.spec_from_file_location("pn2d_srh_mesh_matrix", SCRIPT)
matrix = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = matrix
SPEC.loader.exec_module(matrix)

RUNNER_SCRIPT = (
    Path(__file__).resolve().parents[2]
    / "scripts"
    / "run_pn2d_bv_off_srh_mesh_matrix.py"
)
RUNNER_SPEC = importlib.util.spec_from_file_location(
    "pn2d_srh_mesh_matrix_runner", RUNNER_SCRIPT
)
runner = importlib.util.module_from_spec(RUNNER_SPEC)
assert RUNNER_SPEC.loader is not None
sys.modules[RUNNER_SPEC.name] = runner
RUNNER_SPEC.loader.exec_module(runner)


class TestMeshMatrix(unittest.TestCase):
    def test_nested_junction_spacing(self):
        self.assertEqual(matrix.LEVELS["M0"]["junction_x"], 1.0 / 3.0)
        self.assertEqual(
            matrix.LEVELS["M0"]["junction_x"]
            / matrix.LEVELS["M1"]["junction_x"],
            2.0,
        )
        self.assertEqual(
            matrix.LEVELS["M0"]["junction_x"]
            / matrix.LEVELS["M2"]["junction_x"],
            4.0,
        )

    def test_refinement_replacement_is_scoped(self):
        source = '(sdedr:define-refinement-size "Junction.Mesh" 1 2 3 4)'
        updated = matrix.replace_size(source, "Junction.Mesh", 0.125)
        self.assertIn("0.125 0.125 0.125 0.125", updated)

    def test_all_levels_use_balanced_half_species_junction(self):
        source = (
            "(define XJ 1.0)      ; PN junction position\n"
            "(sdedr:define-refeval-window\n"
            '  "P.Window"\n'
            '  "Rectangle"\n'
            "  (position 0.0 0.0 0.0)\n"
            "  (position XJ H 0.0)\n"
            ")\n"
            "(sdedr:define-refeval-window\n"
            '  "N.Window"\n'
            '  "Rectangle"\n'
            "  (position XJ 0.0 0.0)\n"
            "  (position L H 0.0)\n"
            ")\n"
            '(sdedr:define-constant-profile\n'
            '  "N.Doping"\n'
            '  "PhosphorusActiveConcentration"\n'
            "  1e17\n"
            ")\n"
            '(sdedr:define-constant-profile-placement\n'
            '  "N.Place"\n'
            '  "N.Doping"\n'
            '  "N.Window"\n'
            ")\n"
        )
        updated = matrix.make_dose_preserving_junction(source)
        self.assertIn("(define XJ_LEFT (- XJ 0.001))", updated)
        self.assertIn("(define XJ_RIGHT (+ XJ 0.001))", updated)
        self.assertIn("(position XJ_LEFT H 0.0)", updated)
        self.assertIn("(position XJ_RIGHT 0.0 0.0)", updated)
        self.assertEqual(updated.count("(position XJ H 0.0)"), 0)
        self.assertEqual(updated.count("(position XJ 0.0 0.0)"), 0)
        self.assertIn('"Junction.Doping.Window"', updated)
        self.assertIn('"Junction.P.Half"', updated)
        self.assertIn('"Junction.N.Half"', updated)
        self.assertEqual(updated.count("  5e16"), 2)

    def test_refined_levels_add_window_without_changing_profile_contract(self):
        source = "; header\n;----------------------------------------------------------\n; Build mesh\n"
        refined = matrix.add_junction_refinement(source, 1.0 / 6.0, 0.125)
        self.assertIn('"Junction.Window"', refined)
        self.assertIn('"Junction.Mesh"', refined)
        self.assertTrue(refined.endswith("; Build mesh\n"))

    def test_mesh_metrics_preserve_area_contacts_and_dose(self):
        mesh = {
            "nodes": [
                {"id": 0, "x": 0.0, "y": 0.0},
                {"id": 1, "x": 2.0, "y": 0.0},
                {"id": 2, "x": 0.0, "y": 1.0},
            ],
            "triangles": [{"id": 0, "region_id": 0, "node_ids": [0, 1, 2]}],
            "contacts": [
                {"name": "Anode", "node_ids": [0, 2]},
                {"name": "Cathode", "node_ids": [1]},
            ],
        }
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            mesh_path = root / "mesh.json"
            doping_path = root / "doping.csv"
            mesh_path.write_text(json.dumps(mesh), encoding="utf-8")
            doping_path.write_text(
                "node_id,donors_cm3,acceptors_cm3\n"
                "0,1e17,0\n1,1e17,0\n2,1e17,0\n",
                encoding="utf-8",
            )
            actual = runner.mesh_metrics(mesh_path, doping_path)
        self.assertEqual(actual["bounds_um"], [0.0, 0.0, 2.0, 1.0])
        self.assertAlmostEqual(actual["area_um2"], 1.0)
        self.assertAlmostEqual(actual["total_impurity_dose_cm3_um2"], 1.0e17)


if __name__ == "__main__":
    unittest.main()
