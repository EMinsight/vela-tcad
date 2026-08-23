from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from scripts.analyze_slot_ldmos_vela_ialmob_ablation import (
    read_vtk_node_scalar_peak,
    relative_delta,
)


class SlotLdmosVelaIALMobAnalysisTest(unittest.TestCase):
    def test_reads_node_scalar_peak_coordinate(self) -> None:
        vtk = """# vtk DataFile Version 3.0
test
ASCII
DATASET UNSTRUCTURED_GRID
POINTS 3 double
0 0 0
1 0 0
0 2 0
POINT_DATA 3
SCALARS AvalancheGeneration double 1
LOOKUP_TABLE default
1
9
3
"""
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "state.vtk"
            path.write_text(vtk, encoding="utf-8")
            peak = read_vtk_node_scalar_peak(path, "AvalancheGeneration")
        self.assertEqual(peak["node_index"], 1)
        self.assertEqual(peak["x_um"], 1.0)
        self.assertEqual(peak["value"], 9.0)

    def test_relative_delta_uses_off_denominator(self) -> None:
        self.assertAlmostEqual(relative_delta(9.0, 10.0), -0.1)


if __name__ == "__main__":
    unittest.main()
