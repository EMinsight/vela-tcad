import importlib.util
import math
import sys
import unittest
from pathlib import Path


SCRIPT = (
    Path(__file__).resolve().parents[2]
    / "scripts"
    / "build_pn2d_bv_off_srh_spatial_report.py"
)
SPEC = importlib.util.spec_from_file_location(
    "build_pn2d_bv_off_srh_spatial_report", SCRIPT
)
report = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = report
SPEC.loader.exec_module(report)


class TestTriangleSourceIntegration(unittest.TestCase):
    def test_constant_source_is_exact(self):
        points = ((0.0, 0.0), (2.0, 0.0), (0.0, 3.0))
        value = 7.25
        integral = report.integrate_linear_triangle(points, (value, value, value))
        expected = 3.0 * value
        self.assertLessEqual(abs(integral - expected) / abs(expected), 1.0e-12)

    def test_linear_source_is_exact(self):
        points = ((-0.5, 0.25), (1.75, -0.25), (0.25, 2.0))

        def field(x, y):
            return 2.5 - 0.75 * x + 1.25 * y

        values = tuple(field(*point) for point in points)
        area = report.triangle_area(points)
        centroid = (
            sum(point[0] for point in points) / 3.0,
            sum(point[1] for point in points) / 3.0,
        )
        expected = area * field(*centroid)
        integral = report.integrate_linear_triangle(points, values)
        scale = max(abs(expected), 1.0)
        self.assertLessEqual(abs(integral - expected) / scale, 1.0e-12)

    def test_triangle_area_ignores_orientation(self):
        points = ((0.0, 0.0), (0.0, 3.0), (2.0, 0.0))
        self.assertTrue(math.isclose(report.triangle_area(points), 3.0))


if __name__ == "__main__":
    unittest.main()
