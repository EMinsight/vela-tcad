from __future__ import annotations

import unittest

from scripts.analyze_pn2d_bv_feedback_state_substitution import (
    VARIANTS,
    boundary_identity,
    classify,
    variant_metrics,
)
from scripts.run_pn2d_bv_feedback_state_substitution import (
    bias_token,
    sentaurus_node_records,
)


def row(variant: str, node: int, *, contact: bool, trial: float) -> dict[str, str]:
    values = {
        "variant": variant,
        "node_id": str(node),
        "is_contact": "1" if contact else "0",
        "baseline_phin_V": "0",
        "baseline_phip_V": "0",
        "replacement_phin_V": "1",
        "replacement_phip_V": "1",
        "electron_residual": "1",
        "hole_residual": "1",
        "desired_electron_residual": "1",
        "desired_hole_residual": "1",
        "delta_phin_V": str(trial),
        "delta_phip_V": str(trial),
        "trial_phin_V": str(trial),
        "trial_phip_V": str(trial),
        "carrier_only_delta_phin_V": str(trial),
        "carrier_only_delta_phip_V": str(trial),
        "carrier_only_trial_phin_V": str(trial),
        "carrier_only_trial_phip_V": str(trial),
        "electron_closure_error": "0",
        "hole_closure_error": "0",
        "psi_residual": "0",
    }
    if contact:
        values.update(
            {
                "electron_residual": "0",
                "hole_residual": "0",
                "desired_electron_residual": "0",
                "desired_hole_residual": "0",
                "delta_phin_V": "0",
                "delta_phip_V": "0",
                "trial_phin_V": "0",
                "trial_phip_V": "0",
                "carrier_only_delta_phin_V": "0",
                "carrier_only_delta_phip_V": "0",
                "carrier_only_trial_phin_V": "0",
                "carrier_only_trial_phip_V": "0",
            }
        )
    return values


class FeedbackStateSubstitutionTest(unittest.TestCase):
    def test_material_toward_target_update_passes_causal_gate(self) -> None:
        baseline = [
            row("baseline", 0, contact=False, trial=0.0),
            row("baseline", 1, contact=True, trial=0.0),
        ]
        candidate = [
            row("density_only", 0, contact=False, trial=0.5),
            row("density_only", 1, contact=True, trial=0.0),
        ]
        metrics = variant_metrics(candidate, baseline, "density_only")
        self.assertTrue(metrics["causal_gate_passed"])
        self.assertAlmostEqual(metrics["qfp_error_improvement_fraction"], 0.5)
        self.assertAlmostEqual(metrics["residual_direction_cosine"], 1.0)
        self.assertAlmostEqual(metrics["update_direction_cosine"], 1.0)

    def test_boundary_identity_is_fail_closed(self) -> None:
        matrix = {
            variant: [
                row(variant, 0, contact=False, trial=0.0),
                row(variant, 1, contact=True, trial=0.0),
            ]
            for variant in VARIANTS
        }
        self.assertTrue(boundary_identity(matrix)["passed"])
        matrix["qfp_only"][1]["electron_residual"] = "1e-30"
        result = boundary_identity(matrix)
        self.assertFalse(result["passed"])
        self.assertEqual(result["mismatch_count"], 1)

    def test_cross_bias_classification_is_typed(self) -> None:
        self.assertEqual(
            classify(
                {
                    "density_only": True,
                    "qfp_only": False,
                    "density_qfp": True,
                }
            ),
            "density_feedback_cause",
        )
        self.assertEqual(
            classify(
                {
                    "density_only": False,
                    "qfp_only": False,
                    "density_qfp": False,
                }
            ),
            "no_cross_bias_causal_substitution",
        )

    def test_sentaurus_contact_duplicates_are_explicitly_excluded(self) -> None:
        records = []
        for node, coordinate in ((0, [0.0, 0.0]), (1, [1.0, 0.0]), (2, [0.0, 0.0])):
            records.append(
                {
                    "branch": "avalanche_on",
                    "bias_V": -19.7,
                    "quantity": "density",
                    "carrier": "electron",
                    "support_kind": "physical_node",
                    "provenance": "native",
                    "unit": "cm^-3",
                    "support_key": f"node:{node}",
                    "coordinates_um": coordinate,
                    "values": [1.0],
                }
            )
        selected = sentaurus_node_records(
            {"records": records},
            "avalanche_on",
            -19.7,
            "density",
            "electron",
            "cm^-3",
            2,
        )
        self.assertEqual(sorted(selected), [0, 1])
        self.assertEqual(bias_token(-19.7), "m19p700000")


if __name__ == "__main__":
    unittest.main()
