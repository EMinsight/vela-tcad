import math
import unittest

from scripts.pn2d_minimal6_diagnostics.transport_element_closure import (
    ELEMENTARY_CHARGE_C,
    _node_density_models,
    effective_density,
)


class Minimal6TransportElementClosureTest(unittest.TestCase):
    def test_recovers_density_from_collinear_transport_fields(self):
        density = 2.5e21
        mobility = 0.04
        grad = (300.0, -400.0)
        current = tuple(
            ELEMENTARY_CHARGE_C * mobility * density * value
            for value in grad
        )
        result = effective_density(current, mobility, grad)
        self.assertEqual(result["status"], "valid")
        self.assertAlmostEqual(result["density_m3"], density)
        self.assertAlmostEqual(result["angle_deg"], 0.0)
        self.assertLess(result["orthogonal_residual"], 1.0e-15)

    def test_detects_sign_incompatible_current(self):
        result = effective_density((-1.0, 0.0), 0.1, (1.0, 0.0))
        self.assertEqual(result["status"], "sign_incompatible")
        self.assertIsNone(result["density_m3"])
        self.assertAlmostEqual(result["angle_deg"], 180.0)

    def test_density_controls_are_exact(self):
        models = _node_density_models([1.0, 4.0, 16.0])
        self.assertEqual(models["arithmetic"], 7.0)
        self.assertAlmostEqual(models["geometric"], 4.0)
        self.assertAlmostEqual(models["harmonic"], 16.0 / 7.0)
        self.assertTrue(math.isfinite(models["harmonic"]))


if __name__ == "__main__":
    unittest.main()
