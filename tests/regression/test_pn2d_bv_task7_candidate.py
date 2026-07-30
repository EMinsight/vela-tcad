from __future__ import annotations

import unittest

from scripts.analyze_pn2d_bv_task7_candidate import (
    estimator_error,
    estimator_improved,
    knee_rmse,
    paired_metric,
    vela_on_reverse_intervals,
)


class Task7CandidateTests(unittest.TestCase):
    def test_curve_metrics_are_exact(self) -> None:
        acceptance = {
            "knee_error_rows": [
                {"absolute_log_error_dex": 3.0},
                {"absolute_log_error_dex": 4.0},
            ],
            "knee_estimators": {
                "sentaurus": {"V_break": -20.0, "V_slope": -19.8},
                "vela": {"V_break": -19.5, "V_slope": None},
            },
        }
        self.assertAlmostEqual(knee_rmse(acceptance), 5.0 / (2.0**0.5))
        self.assertEqual(estimator_error(acceptance, "V_break"), 0.5)
        self.assertIsNone(estimator_error(acceptance, "V_slope"))
        self.assertTrue(estimator_improved(None, 0.02))
        self.assertTrue(estimator_improved(0.2, 0.02))
        self.assertFalse(estimator_improved(0.2, None))

    def test_reverse_interval_gate_tracks_candidate_on_curve(self) -> None:
        acceptance = {
            "curve_rows": [
                {
                    "bias_V": 0.0,
                    "vela_on_A_per_um": 1.0,
                },
                {
                    "bias_V": -1.0,
                    "vela_on_A_per_um": -2.0,
                },
                {
                    "bias_V": -2.0,
                    "vela_on_A_per_um": -1.5,
                },
            ]
        }
        self.assertEqual(
            vela_on_reverse_intervals(acceptance),
            [
                {
                    "left_bias_V": -1.0,
                    "right_bias_V": -2.0,
                    "left_abs_current_A_per_um": 2.0,
                    "right_abs_current_A_per_um": 1.5,
                }
            ],
        )

    def test_internal_metric_matches_exact_support(self) -> None:
        def payload(offset: float) -> dict:
            return {
                "records": [
                    {
                        "bias_V": -19.7,
                        "branch": "avalanche_on",
                        "stage": "state",
                        "quantity": "quasi_fermi",
                        "carrier": carrier,
                        "support_key": f"node:{node}",
                        "values": [float(node) + offset],
                    }
                    for carrier in ("electron", "hole")
                    for node in range(2)
                ]
            }

        self.assertEqual(
            paired_metric(
                payload(0.0),
                payload(0.25),
                stage="state",
                quantity="quasi_fermi",
                logarithmic=False,
            ),
            0.25,
        )

    def test_density_metric_uses_log_space(self) -> None:
        sentaurus = {
            "records": [
                {
                    "bias_V": -19.7,
                    "branch": "avalanche_on",
                    "stage": "density",
                    "quantity": "density",
                    "carrier": "hole",
                    "support_key": "node:14",
                    "values": [100.0],
                }
            ]
        }
        vela = {
            "records": [
                {
                    **sentaurus["records"][0],
                    "values": [1000.0],
                }
            ]
        }
        self.assertEqual(
            paired_metric(
                sentaurus,
                vela,
                stage="density",
                quantity="density",
                logarithmic=True,
            ),
            1.0,
        )


if __name__ == "__main__":
    unittest.main()
