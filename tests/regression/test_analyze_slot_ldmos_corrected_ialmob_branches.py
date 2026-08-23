import math
import unittest

from scripts.analyze_slot_ldmos_corrected_ialmob_branches import locate_bvds


class CorrectedIalmobBranchAnalysisTest(unittest.TestCase):
    def test_log_interpolates_threshold(self) -> None:
        result = locate_bvds(
            [
                {"voltage_V": 10.0, "current_A_per_um": 1.0e-8},
                {"voltage_V": 12.0, "current_A_per_um": 1.0e-6},
            ],
            1.0e-7,
        )
        self.assertEqual(result["status"], "located")
        self.assertEqual(result["interpolation"], "log_current")
        self.assertAlmostEqual(result["bvds_V"], 11.0)

    def test_reports_unreached_criterion(self) -> None:
        result = locate_bvds(
            [{"voltage_V": 12.0, "current_A_per_um": 2.0e-9}], 1.0e-7
        )
        self.assertEqual(result["status"], "criterion_not_reached")
        self.assertTrue(math.isclose(result["maximum_current_A_per_um"], 2.0e-9))


if __name__ == "__main__":
    unittest.main()
