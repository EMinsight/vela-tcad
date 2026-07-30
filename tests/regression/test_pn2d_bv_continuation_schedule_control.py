from __future__ import annotations

import csv
import tempfile
import unittest
from pathlib import Path

from scripts.analyze_pn2d_bv_continuation_schedule_control import (
    classify,
    curve_metrics,
    reversal_metric,
    state_metrics,
)


class ContinuationScheduleControlTests(unittest.TestCase):
    def test_classification_requires_deterministic_controls(self) -> None:
        self.assertEqual(
            classify(
                controls_valid=False,
                state_invariant=True,
                baseline_reversal=True,
                refined_reversal=False,
                refined_resolved=True,
            ),
            "insufficient_or_nondeterministic_control",
        )

    def test_classification_separates_resolution_from_invariance(self) -> None:
        self.assertEqual(
            classify(
                controls_valid=True,
                state_invariant=False,
                baseline_reversal=True,
                refined_reversal=False,
                refined_resolved=True,
            ),
            "continuation_schedule_resolves_cross_block_reversal",
        )
        self.assertEqual(
            classify(
                controls_valid=True,
                state_invariant=True,
                baseline_reversal=True,
                refined_reversal=True,
                refined_resolved=False,
            ),
            "continuation_invariant_cross_block_reversal",
        )

    def test_reversal_metric_requires_carrier_improvement_and_full_reversal(self) -> None:
        metric = reversal_metric(
            {
                "carrier_only_qfp_error_improvement_fraction": 0.13,
                "carrier_only_update_direction_cosine": 0.64,
                "carrier_only_no_carrier_worsening": True,
                "qfp_error_improvement_fraction": -0.07,
                "update_direction_cosine": -0.2,
                "no_carrier_worsening": False,
            }
        )
        self.assertTrue(metric["cross_block_reversal"])
        self.assertFalse(metric["cross_block_reversal_resolved"])

    def test_state_metrics_compare_qfp_and_log_density(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            left = root / "left.csv"
            right = root / "right.csv"
            fields = [
                "node_id",
                "psi",
                "phin",
                "phip",
                "electrons_m3",
                "holes_m3",
            ]
            rows = [
                {
                    "node_id": 0,
                    "psi": 1.0,
                    "phin": 0.5,
                    "phip": -0.5,
                    "electrons_m3": 1.0e20,
                    "holes_m3": 1.0e21,
                }
            ]
            for path, delta in ((left, 0.0), (right, 1.0e-9)):
                with path.open("w", newline="", encoding="utf-8") as handle:
                    writer = csv.DictWriter(handle, fieldnames=fields)
                    writer.writeheader()
                    row = dict(rows[0])
                    row["psi"] += delta
                    row["phin"] += delta
                    row["phip"] += delta
                    writer.writerow(row)
            metrics = state_metrics(left, right)
            self.assertAlmostEqual(metrics["psi_max_abs_V"], 1.0e-9)
            self.assertAlmostEqual(metrics["qfp_max_abs_V"], 1.0e-9)
            self.assertEqual(metrics["density_log_max_abs_dex"], 0.0)

    def test_curve_metrics_use_exact_bias_and_log_current(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            fields = ["node_id", "bias_V", "current_total_A_per_um"]
            for path, scale in (
                (root / "left.csv", 1.0),
                (root / "right.csv", 10.0),
            ):
                with path.open("w", newline="", encoding="utf-8") as handle:
                    writer = csv.DictWriter(handle, fieldnames=fields)
                    writer.writeheader()
                    writer.writerow(
                        {
                            "node_id": 0,
                            "bias_V": -19.7,
                            "current_total_A_per_um": 1.0e-12 * scale,
                        }
                    )
            metrics = curve_metrics(root / "left.csv", root / "right.csv")
            self.assertEqual(metrics["point_count"], 1)
            self.assertAlmostEqual(metrics["log_current_max_abs_dex"], 1.0)


if __name__ == "__main__":
    unittest.main()
