from __future__ import annotations

import unittest

from scripts.analyze_pn2d_bv_poisson_qfp_cross_block import (
    classify_bias,
    mode_metrics,
)


class PoissonQfpCrossBlockTests(unittest.TestCase):
    def test_mode_metrics_use_interior_qfp_target(self) -> None:
        rows = [
            {
                "is_contact": "0",
                "baseline_phin_V": "0",
                "baseline_phip_V": "0",
                "replacement_phin_V": "1",
                "replacement_phip_V": "-1",
                "target_delta_phin_V": "1",
                "target_delta_phip_V": "-1",
                "independent_delta_phin_V": "0.5",
                "independent_delta_phip_V": "-0.5",
            }
        ]
        metric = mode_metrics(rows, "independent")
        self.assertAlmostEqual(metric["update_direction_cosine"], 1.0)
        self.assertEqual(metric["qfp_error_improvement_fraction"], 0.5)

    def test_classification_requires_both_cross_blocks_for_reversal(self) -> None:
        def metric(improvement: float, direction: float) -> dict:
            return {
                "qfp_error_improvement_fraction": improvement,
                "update_direction_cosine": direction,
            }

        metrics = {
            "independent": metric(0.13, 0.64),
            "no_psi_qfp": metric(0.10, 0.50),
            "no_qfp_psi": metric(0.13, 0.64),
            "full_raw": metric(-0.07, -0.20),
        }
        self.assertEqual(
            classify_bias(metrics),
            "bidirectional_poisson_qfp_closed_loop_cause",
        )

    def test_direction_evidence_is_not_rejected_by_oversized_unit_step(
        self,
    ) -> None:
        def metric(improvement: float, direction: float) -> dict:
            return {
                "qfp_error_improvement_fraction": improvement,
                "update_direction_cosine": direction,
            }

        metrics = {
            "independent": metric(-1000.0, 0.75),
            "no_psi_qfp": metric(-1000.0, 0.75),
            "no_qfp_psi": metric(-1000.0, 0.75),
            "full_raw": metric(-1.5, -0.44),
        }
        self.assertEqual(
            classify_bias(metrics),
            "bidirectional_poisson_qfp_closed_loop_cause",
        )

    def test_direct_qfp_from_psi_feed_is_separate(self) -> None:
        def metric(improvement: float, direction: float) -> dict:
            return {
                "qfp_error_improvement_fraction": improvement,
                "update_direction_cosine": direction,
            }

        metrics = {
            "independent": metric(0.13, 0.64),
            "no_psi_qfp": metric(-0.10, -0.20),
            "no_qfp_psi": metric(0.13, 0.64),
            "full_raw": metric(-0.07, -0.20),
        }
        self.assertEqual(
            classify_bias(metrics),
            "J_qfp_psi_direct_feed_cause",
        )


if __name__ == "__main__":
    unittest.main()
