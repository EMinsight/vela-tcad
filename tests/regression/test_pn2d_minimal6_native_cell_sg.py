import math
import unittest

from scripts.pn2d_minimal6_diagnostics.native_cell_sg_experiment import (
    _edge_cells,
    _native_edge_current,
)


class Minimal6NativeCellSgTest(unittest.TestCase):
    def test_edge_cell_adjacency_is_exact(self):
        triangles = ((0, 4, 5), (0, 5, 1), (1, 5, 2), (5, 3, 2))
        adjacent = _edge_cells(triangles)
        self.assertEqual(adjacent[(0, 4)], (0,))
        self.assertEqual(adjacent[(0, 5)], (0, 1))
        self.assertEqual(adjacent[(1, 5)], (1, 2))
        self.assertEqual(adjacent[(2, 5)], (2, 3))
        self.assertEqual(len(adjacent), 9)

    def test_native_internal_edge_uses_both_cell_projections(self):
        currents = {
            ("mirror", -1.0, "electron", 0): (2.0, 0.0),
            ("mirror", -1.0, "electron", 1): (6.0, 4.0),
        }
        tangent = (1.0 / math.sqrt(2.0), 1.0 / math.sqrt(2.0))
        mean, projections = _native_edge_current(
            currents,
            "mirror",
            -1.0,
            "electron",
            (0, 1),
            tangent,
        )
        self.assertAlmostEqual(projections[0], math.sqrt(2.0))
        self.assertAlmostEqual(projections[1], 5.0 * math.sqrt(2.0))
        self.assertAlmostEqual(mean, 3.0 * math.sqrt(2.0))


if __name__ == "__main__":
    unittest.main()
