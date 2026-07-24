import unittest

from scripts.pn2d_minimal6_diagnostics.phase_f_self_consistent import (
    classify_log_error,
    directed_current_A_per_um,
    first_failed_metric,
    percentile,
)


class Minimal6PhaseFSelfConsistentTest(unittest.TestCase):
    def test_percentile_uses_frozen_linear_interpolation(self):
        self.assertEqual(percentile([0.0, 10.0], 0.5), 5.0)
        self.assertAlmostEqual(percentile([0.0, 10.0], 0.95), 9.5)

    def test_log_error_preserves_exact_zero_and_zero_mismatch(self):
        self.assertEqual(classify_log_error(0.0, 0.0), ("exact_zero", None, None))
        self.assertEqual(
            classify_log_error(0.0, 1.0),
            ("zero_reference_mismatch", None, None),
        )
        status, error, sign = classify_log_error(-2.0, -1.0)
        self.assertEqual(status, "valid")
        self.assertAlmostEqual(error, 0.3010299956639812)
        self.assertEqual(sign, 1.0)

    def test_directed_current_respects_reference_orientation(self):
        value = directed_current_A_per_um(
            current_density_A_per_m2=4.0,
            dual_length_m=0.25,
            candidate_node0=5,
            candidate_node1=1,
            reference_node0=1,
            reference_node1=5,
        )
        self.assertEqual(value, -1.0e-6)

    def test_first_failed_metric_follows_dependency_order(self):
        gates = {
            "psi": True,
            "electron_qfp": False,
            "hole_qfp": False,
            "electron_density": True,
        }
        self.assertEqual(first_failed_metric(gates), "electron_qfp")


if __name__ == "__main__":
    unittest.main()
