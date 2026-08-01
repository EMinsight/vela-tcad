from __future__ import annotations

import unittest

from scripts.create_obtuse_dfise_probe import replace_vertex, triangle_angles


class CreateObtuseDfiseProbeTests(unittest.TestCase):
    def test_replaces_one_vertex_and_preserves_header(self) -> None:
        source = """DF-ISE text
Data {
  Vertices (3) {
 0.0 0.0
 1.0 0.0
 0.0 1.0
  }
}
"""
        result, old = replace_vertex(source, 1, 0.1, 0.0)
        self.assertEqual(old, (1.0, 0.0))
        self.assertIn("Vertices (3)", result)
        self.assertIn("1.000000000000000e-01 0.000000000000000e+00", result)

    def test_expected_probe_triangle_is_obtuse(self) -> None:
        angles = triangle_angles([(0.25, 0.25), (0.0, 0.0), (0.05, 0.0)])
        self.assertGreater(max(angles), 90.0)


if __name__ == "__main__":
    unittest.main()
