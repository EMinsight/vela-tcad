from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
MODULE_PATH = REPO / "scripts" / "analyze_slot_ldmos_bvds_reference.py"
SPEC = importlib.util.spec_from_file_location(
    "analyze_slot_ldmos_bvds_reference", MODULE_PATH
)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class SlotLdmosBvdsReferenceAnalysisTests(unittest.TestCase):
    def test_load_line_and_threshold_interpolation(self) -> None:
        row = {
            "time": 1.0,
            "outer_voltage_V": 1.0,
            "inner_voltage_V": 0.9994,
            "drain_total_current_A_per_um": 6.0e-16,
        }
        self.assertAlmostEqual(MODULE.load_line(row)["residual_V"], 0.0)
        crossing = MODULE.interpolate_crossing(
            [
                {**row, "outer_voltage_V": 40.0,
                 "inner_voltage_V": 38.4,
                 "drain_total_current_A_per_um": 0.8e-7},
                {**row, "outer_voltage_V": 60.0,
                 "inner_voltage_V": 38.6,
                 "drain_total_current_A_per_um": 1.2e-7},
            ],
            1.0e-7,
        )
        self.assertIsNotNone(crossing)
        assert crossing is not None
        self.assertAlmostEqual(crossing["outer_voltage_V"], 50.0)
        self.assertAlmostEqual(crossing["inner_voltage_V"], 38.5)

    def test_parses_sentaurus_mesh_and_iic_log(self) -> None:
        text = """
 Total  3.8195124e+01 3.8195130e+01 1.5e-05 27275 2 ( 0.01 %) 5.3797532e-05 ( 0.0001)
 Maximum electric field: 2.3515e+06 V/cm at (6.0533e-03,4.42555) um
 Ionization-Integrals:
 Electron: 0.12874
 Hole: 0.0983688
"""
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "stage.log"
            path.write_text(text)
            mesh = MODULE.parse_mesh_log(path)
            iic = MODULE.parse_iic_log(path)
        self.assertEqual(mesh["non_delaunay_cell_count"], 2)
        self.assertEqual(mesh["cell_count"], 27275)
        self.assertAlmostEqual(iic["electron_ionization_integral"], 0.12874)


if __name__ == "__main__":
    unittest.main()
