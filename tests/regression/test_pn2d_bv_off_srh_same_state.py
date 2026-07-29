import importlib.util
import itertools
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = (
    Path(__file__).resolve().parents[2]
    / "scripts"
    / "decompose_pn2d_bv_off_srh_same_state.py"
)
SPEC = importlib.util.spec_from_file_location("pn2d_bv_off_srh_same_state", SCRIPT)
audit = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = audit
SPEC.loader.exec_module(audit)


class TestSameStateElementIntegration(unittest.TestCase):
    def test_exact_current_accepts_native_plt_export_column(self):
        with tempfile.TemporaryDirectory() as directory:
            curve = Path(directory) / "curve.csv"
            curve.write_text(
                "bias_V,current_total_A_per_um,current_total\n"
                "-1,,1.25e-17\n"
                "-2,2.5e-17,9.9e-17\n",
                encoding="utf-8",
            )
            actual = audit.exact_current(
                curve, ("current_total_A_per_um", "current_total")
            )
        self.assertEqual(actual, {1: 1.25e-17, 2: 2.5e-17})

    def test_permutation_invariance(self):
        points = ((0.0, 0.0), (2.0, 0.0), (0.0, 3.0))
        values = (1.0, 4.0, -2.0)
        expected = audit.integrate_element(points, values)
        for permutation in itertools.permutations(range(3)):
            actual = audit.integrate_element(
                tuple(points[index] for index in permutation),
                tuple(values[index] for index in permutation),
            )
            self.assertAlmostEqual(actual, expected, places=14)

    def test_barycentric_conservation(self):
        triangle = ((0.0, 0.0), (2.0, 0.0), (0.0, 2.0))
        weights = audit.barycentric((0.5, 0.75), triangle)
        self.assertIsNotNone(weights)
        self.assertAlmostEqual(sum(weights), 1.0, places=15)
        self.assertGreaterEqual(min(weights), 0.0)


if __name__ == "__main__":
    unittest.main()
