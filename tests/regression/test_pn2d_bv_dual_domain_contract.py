from __future__ import annotations

import copy
import json
import unittest
from pathlib import Path

from scripts.review_pn2d_bv_dual_domain_contract import evaluate


ROOT = Path(__file__).resolve().parents[2]
CONTRACT = json.loads(
    (
        ROOT
        / "docs"
        / "validation"
        / "contracts"
        / "pn2d_bv_dual_domain_acceptance_v1.json"
    ).read_text(encoding="utf-8")
)


def fixtures() -> tuple[dict, dict, dict, dict]:
    rows = []
    all_biases = (
        CONTRACT["domains"]["bv_model_consistency"]["exact_biases_V"]
        + CONTRACT["domains"]["low_current_solver_precision"]["exact_biases_V"]
    )
    for bias in all_biases:
        magnitude = abs(float(bias))
        sentaurus = 3.0e-17 * (1.0 + max(magnitude - 15.0, 0.0) * 0.2)
        vela = sentaurus * 1.01
        rows.append(
            {
                "bias_V": float(bias),
                "vela_on_A_per_um": vela,
                "vela_off_A_per_um": 3.0e-17,
                "sentaurus_on_A_per_um": sentaurus,
                "sentaurus_off_A_per_um": 3.0e-17,
                "vela_gain": vela / 3.0e-17,
                "sentaurus_gain": sentaurus / 3.0e-17,
            }
        )
    curve = {
        "curve_rows": rows,
        "knee_metrics": {
            "median_absolute_log_error_dex": 0.01,
            "maximum_absolute_log_error_dex": 0.02,
        },
        "knee_estimators": {
            "vela": {"V_break": -19.65, "V_slope": -19.81},
            "sentaurus": {"V_break": -19.66, "V_slope": -19.82},
        },
        "adjacent_slope_rmse_dex_per_V": 0.03,
    }
    candidate = {
        "outcome": "tradeoff_without_parity",
        "gates": {
            name: True
            for name in CONTRACT["domains"]["bv_model_consistency"][
                "required_self_consistent_gates"
            ]
        },
    }
    frozen = {
        "fixed_state_gates": {
            name: True
            for name in CONTRACT["domains"]["bv_model_consistency"][
                "required_fixed_state_gates"
            ]
        }
    }
    low = {
        "observation_only": True,
        "production_default_changed": False,
        "outcome": CONTRACT["domains"]["low_current_solver_precision"][
            "required_typed_outcome"
        ],
        "gates": {
            name: True
            for name in CONTRACT["domains"]["low_current_solver_precision"][
                "required_evidence_gates"
            ]
        },
    }
    return curve, candidate, frozen, low


class DualDomainContractTest(unittest.TestCase):
    def test_separates_bv_pass_from_classified_precision_floor(self) -> None:
        result = evaluate(CONTRACT, *fixtures())
        self.assertTrue(result["bv_model_consistency"]["passed"])
        self.assertTrue(
            result["low_current_solver_precision"][
                "classified_as_precision_floor"
            ]
        )
        self.assertEqual(
            result["decision"],
            "bv_model_consistent_low_current_precision_floor_open",
        )
        self.assertEqual(
            result["next_stage_authorization"],
            "opt_in_bv_model_validation_only",
        )
        self.assertFalse(result["production_default_change_authorized"])

    def test_high_current_runaway_cannot_use_low_current_waiver(self) -> None:
        curve, candidate, frozen, low = fixtures()
        row = next(row for row in curve["curve_rows"] if row["bias_V"] == -3.0)
        row["vela_on_A_per_um"] = 1.0e-6
        result = evaluate(CONTRACT, curve, candidate, frozen, low)
        self.assertTrue(result["bv_model_consistency"]["passed"])
        self.assertFalse(
            result["low_current_solver_precision"][
                "classified_as_precision_floor"
            ]
        )
        self.assertEqual(
            result["decision"],
            "bv_model_consistent_low_current_unclassified_fail_closed",
        )

    def test_failed_intermediate_quantity_gate_fails_bv_domain(self) -> None:
        curve, candidate, frozen, low = fixtures()
        frozen = copy.deepcopy(frozen)
        frozen["fixed_state_gates"][
            "integrated_source_within_0p02_relative"
        ] = False
        result = evaluate(CONTRACT, curve, candidate, frozen, low)
        self.assertFalse(result["bv_model_consistency"]["passed"])
        self.assertEqual(result["decision"], "bv_model_consistency_failed")
        self.assertEqual(result["next_stage_authorization"], "none")

    def test_domains_are_disjoint_and_historical_score_is_preserved(self) -> None:
        result = evaluate(CONTRACT, *fixtures())
        checks = result["independent_review"]["checks"]
        self.assertTrue(checks["bias_domains_disjoint"])
        self.assertTrue(checks["historical_score_preserved"])
        self.assertTrue(
            checks["low_current_raw_monotonicity_excluded_from_bv_gate"]
        )


if __name__ == "__main__":
    unittest.main()
