#!/usr/bin/env python3
"""Tests for the Sentaurus avalanche driving-force comparison."""

from __future__ import annotations

import math
import unittest

from scripts.compare_pn2d_minimal6_sentaurus_avalanche_drive_controls import (
    abs_dex,
    compare_group,
    symmetric_relative_error,
    vector_angle_error_deg,
    vector_relative_error,
)


class SentaurusAvalancheDriveComparisonTest(unittest.TestCase):
    def test_scalar_metrics_distinguish_exact_ratio_and_sign(self) -> None:
        self.assertEqual(abs_dex(100.0, 10.0), 1.0)
        self.assertEqual(abs_dex(-100.0, 10.0), 1.0)
        self.assertIsNone(abs_dex(0.0, 10.0))
        self.assertEqual(symmetric_relative_error(5.0, 5.0), 0.0)
        self.assertEqual(symmetric_relative_error(-5.0, 5.0), 2.0)

    def test_vector_metrics_capture_magnitude_and_direction(self) -> None:
        self.assertAlmostEqual(
            vector_relative_error((2.0, 0.0), (1.0, 0.0)),
            0.5,
        )
        self.assertAlmostEqual(
            vector_angle_error_deg((0.0, 1.0), (1.0, 0.0)),
            90.0,
        )
        self.assertEqual(
            vector_angle_error_deg((3.0, 4.0), (3.0, 4.0)), 0.0
        )
        self.assertIsNone(vector_angle_error_deg((0.0, 0.0), (1.0, 0.0)))

    def test_group_comparison_requires_same_entities_and_marks_exact(self) -> None:
        reference = [
            {
                "bias_V": -10,
                "vertex": 0,
                "psi_V": -1.0,
                "alpha_n_cm_inv": 10.0,
            }
        ]
        exact = compare_group(
            topology="mirror",
            candidate_variant="implicit_default",
            group="vertices",
            reference_rows=reference,
            candidate_rows=[dict(reference[0])],
        )
        self.assertTrue(all(row["exact"] == 1 for row in exact))

        changed = [dict(reference[0], alpha_n_cm_inv=100.0)]
        compared = compare_group(
            topology="mirror",
            candidate_variant="explicit_electric_field",
            group="vertices",
            reference_rows=reference,
            candidate_rows=changed,
        )
        alpha = next(
            row for row in compared if row["field"] == "alpha_n_cm_inv"
        )
        self.assertEqual(alpha["exact"], 0)
        self.assertTrue(
            math.isclose(alpha["absolute_log10_ratio_dex"], 1.0)
        )

        with self.assertRaisesRegex(ValueError, "identity mismatch"):
            compare_group(
                topology="mirror",
                candidate_variant="implicit_default",
                group="vertices",
                reference_rows=reference,
                candidate_rows=[dict(reference[0], vertex=1)],
            )


if __name__ == "__main__":
    unittest.main()
