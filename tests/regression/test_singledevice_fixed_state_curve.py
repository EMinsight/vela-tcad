import importlib.util
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "singledevice_fixed", REPO / "scripts" / "run_singledevice_fixed_state_curve.py")
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class SingleDeviceFixedStateCurveTest(unittest.TestCase):
    def test_bias_grid_matches_sentaurus_current_plot(self) -> None:
        self.assertAlmostEqual(MODULE.bias_for_index(0), -0.5)
        self.assertAlmostEqual(MODULE.bias_for_index(10), 0.85)
        self.assertAlmostEqual(MODULE.bias_for_index(20), 2.2)


if __name__ == "__main__":
    unittest.main()
