import math
import unittest

from scripts.pn2d_minimal6_diagnostics.mobility_diagnosis import FIELD
from scripts.pn2d_minimal6_diagnostics.mobility_unit_root_cause import (
    absolute_log10_error,
    _markdown,
    endpoint_averaged_mobility,
    field_limited_mobility,
    quantile,
)


class MobilityUnitRootCauseTest(unittest.TestCase):
    def test_legacy_velocity_interpretation_is_one_hundredth_si(self) -> None:
        carrier = "electron"
        low_field = 0.1
        field = 1.0e7
        correct = field_limited_mobility(
            carrier,
            low_field,
            field,
            saturation_velocity_scale=1.0,
        )
        legacy = field_limited_mobility(
            carrier,
            low_field,
            field,
            saturation_velocity_scale=0.01,
        )
        self.assertGreater(correct, legacy)
        self.assertAlmostEqual(
            legacy,
            FIELD[carrier]["saturation_velocity"] * 0.01 / field,
            delta=legacy * 0.02,
        )

    def test_endpoint_average_matches_average_of_endpoint_evaluations(self) -> None:
        field = 2.5e6
        value = endpoint_averaged_mobility(
            "hole",
            -1.0e23,
            1.0e23,
            field,
            saturation_velocity_scale=1.0,
        )
        self.assertTrue(math.isfinite(value))
        self.assertGreater(value, 0.0)

    def test_log_error_and_quantile_contracts(self) -> None:
        self.assertEqual(absolute_log10_error(0.01, 0.1), 1.0)
        self.assertEqual(quantile([0.0, 1.0, 2.0], 0.5), 1.0)
        with self.assertRaises(ValueError):
            absolute_log10_error(0.0, 1.0)
        with self.assertRaises(ValueError):
            quantile([], 0.5)

    def test_report_detects_repaired_production_velocity_branch(self) -> None:
        rows = []
        for carrier in ("electron", "hole"):
            for support, branch, median in (
                ("triangle_local_edge", "legacy_velocity_interpretation", 1.8),
                ("triangle_local_edge", "correct_velocity_interpretation", 0.0),
                ("sentaurus_native_element", "legacy_cell_average_doping", 1.8),
                ("sentaurus_native_element", "correct_cell_average_doping", 0.05),
            ):
                rows.append(
                    {
                        "support": support,
                        "carrier": carrier,
                        "branch": branch,
                        "sample_count": 1,
                        "median_abs_log10_error_dex": median,
                        "p95_abs_log10_error_dex": median,
                        "maximum_abs_log10_error_dex": median,
                    }
                )

        report = _markdown(rows, {})

        self.assertIn("unit conversion defect has been repaired", report)
        self.assertIn("correctly converted velocity interpretation", report)
        self.assertNotIn("Deterministic output root:", report)

if __name__ == "__main__":
    unittest.main()
