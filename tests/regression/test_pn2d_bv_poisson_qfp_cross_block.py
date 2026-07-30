from __future__ import annotations

import unittest

from scripts.analyze_pn2d_bv_poisson_qfp_cross_block import (
    classify_bias,
    leave_out_classification,
    mode_metrics,
    schur_loop_decomposition,
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

    def test_leave_out_classification_names_unique_model(self) -> None:
        metrics = {
            "leave_out_transport_boundary": {
                "update_direction_cosine": 0.2,
            },
            "leave_out_srh_auger": {
                "update_direction_cosine": -0.3,
            },
            "leave_out_sg_avalanche": {
                "update_direction_cosine": -0.4,
            },
            "only_transport_boundary": {
                "update_direction_cosine": 0.1,
            },
            "only_srh_auger": {
                "update_direction_cosine": 0.2,
            },
            "only_sg_avalanche": {
                "update_direction_cosine": 0.3,
            },
        }
        result = leave_out_classification(metrics)
        self.assertEqual(
            result["classification"],
            "transport_boundary_necessary_for_reversal",
        )

    def test_isolated_transport_and_avalanche_can_both_be_adverse(
        self,
    ) -> None:
        metrics = {}
        for component in (
            "transport_boundary",
            "srh_auger",
            "sg_avalanche",
        ):
            metrics[f"leave_out_{component}"] = {
                "update_direction_cosine": -0.2,
            }
            metrics[f"only_{component}"] = {
                "update_direction_cosine":
                    0.3 if component == "srh_auger" else -0.3,
            }
        result = leave_out_classification(metrics)
        self.assertEqual(
            result["classification"],
            "transport_and_avalanche_independently_sustain_reversal",
        )

    def test_schur_loop_decomposition_splits_carrier_signs(self) -> None:
        rows = []
        for component, value in (
            ("transport_boundary", 2.0),
            ("srh_auger", -3.0),
            ("sg_avalanche", 4.0),
        ):
            rows.append(
                {
                    "matrix": "C_Ainv_B_component",
                    "component": component,
                    "row_carrier": "electron",
                    "row_node": "1",
                    "row_x": "1",
                    "row_y": "0",
                    "col_carrier": "hole",
                    "col_node": "2",
                    "col_x": "2",
                    "col_y": "0",
                    "value": str(value),
                }
            )
        rows.append(
            {
                "matrix": "C_Ainv_B",
                "component": "all",
                "row_carrier": "electron",
                "row_node": "1",
                "row_x": "1",
                "row_y": "0",
                "col_carrier": "hole",
                "col_node": "2",
                "col_x": "2",
                "col_y": "0",
                "value": "3",
            }
        )
        result = schur_loop_decomposition(rows)
        pair = result["components"]["srh_auger"]["carrier_pairs"][
            "electron_from_hole"
        ]
        self.assertEqual(pair["negative_l1"], 3.0)
        self.assertEqual(result["total_l2_norm"], 3.0)


if __name__ == "__main__":
    unittest.main()
