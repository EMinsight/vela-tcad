import math
import unittest

from scripts.pn2d_minimal6_diagnostics.inverse_contracts import (
    AcceptanceThresholds,
    Identifiability,
    SampleStatus,
    SupportKind,
)
from scripts.pn2d_minimal6_diagnostics.inverse_replacements import (
    INVERSE_DEPENDENCIES,
    classify_candidate,
    metric_summary,
    rank_candidates,
    run_replacement_matrix,
    run_state_localization_control,
)


FACTORS = tuple(INVERSE_DEPENDENCIES)


def operand(factor, value, **overrides):
    row = {
        "factor": factor,
        "value": value,
        "status": SampleStatus.VALID,
        "support_kind": SupportKind.NODE,
        "support_id": 7,
        "unit_si": "1",
        "carrier": "electron",
        "topology": "sketch",
        "bias_V": -12.0,
    }
    row.update(overrides)
    return row


def evidence(candidate, split, error, *, prediction, status=SampleStatus.VALID,
             metric="gradient_abs_dex", bias=-12.0, **overrides):
    row = {
        "record_kind": "formula_candidate",
        "candidate": candidate,
        "factor": "gradient_recovery",
        "split": split,
        "topology": "sketch" if split == "discovery" else "mirror",
        "bias_V": bias,
        "support_kind": SupportKind.NODE,
        "support_id": 7,
        "carrier": "electron",
        "metric": metric,
        "error": error,
        "prediction_dex": prediction,
        "status": status,
    }
    row.update(overrides)
    return row


class InverseReplacementsTest(unittest.TestCase):
    def test_dependency_order_is_the_declared_physical_chain(self):
        expected = (
            "gradient_recovery", "mobility", "current_semantics",
            "impact_driving_field", "alpha_law", "geometric_integration",
            "source_to_node_mapping",
        )
        self.assertEqual(tuple(INVERSE_DEPENDENCIES), expected)

    def test_one_factor_staged_and_reverse_replacements_close(self):
        contributions = dict(zip(FACTORS, (0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7)))
        baseline = {factor: operand(factor, 1.0) for factor in FACTORS}
        replacement = {
            factor: operand(factor, 10.0 ** contribution)
            for factor, contribution in contributions.items()
        }
        target = 10.0 ** sum(contributions.values())

        result = run_replacement_matrix(
            baseline, replacement, direct_target=target
        )

        for row in result["one_factor"]:
            factor = row["factor"]
            index = FACTORS.index(factor)
            self.assertAlmostEqual(row["delta_dex"], contributions[factor], places=12)
            self.assertEqual(row["changed_stages"], list(FACTORS[index:]))
            for upstream in FACTORS[:index]:
                self.assertEqual(row["stage_values"][upstream], 1.0)

        self.assertEqual(
            [row["factor"] for row in result["forward"]], list(FACTORS)
        )
        for row in result["forward"]:
            self.assertAlmostEqual(
                row["incremental_dex"], contributions[row["factor"]], places=12
            )
        self.assertEqual(
            [row["factor"] for row in result["reverse"]], list(reversed(FACTORS))
        )
        for row in result["reverse"]:
            self.assertAlmostEqual(
                row["incremental_dex"], -contributions[row["factor"]], places=12
            )
        self.assertAlmostEqual(result["full_replacement"], target, places=12)
        self.assertAlmostEqual(result["reverse"][-1]["value"], result["baseline"], places=12)
        self.assertLessEqual(result["closure"]["forward_abs_dex"], 1.0e-10)
        self.assertLessEqual(result["closure"]["reverse_abs_dex"], 1.0e-10)
        self.assertLessEqual(result["closure"]["direct_abs_dex"], 1.0e-10)
        self.assertTrue(all(abs(row["interaction_dex"]) <= 1.0e-12
                            for row in result["adjacent_interactions"]))

    def test_support_carrier_state_and_unit_mismatches_are_rejected(self):
        baseline = {factor: operand(factor, 1.0) for factor in FACTORS}
        for field, value in (
            ("support_kind", SupportKind.EDGE),
            ("carrier", "hole"),
            ("bias_V", -19.0),
            ("unit_si", "V/m"),
        ):
            with self.subTest(field=field):
                replacement = {
                    factor: operand(factor, 2.0) for factor in FACTORS
                }
                replacement["mobility"] = operand("mobility", 2.0, **{field: value})
                with self.assertRaisesRegex(ValueError, field):
                    run_replacement_matrix(baseline, replacement)

    def test_invalid_replacement_values_remain_typed_and_are_not_coerced(self):
        baseline = {factor: operand(factor, 1.0) for factor in FACTORS}
        replacement = {factor: operand(factor, 2.0) for factor in FACTORS}
        replacement["alpha_law"] = operand(
            "alpha_law", None, status=SampleStatus.MISSING_FIELD
        )
        result = run_replacement_matrix(baseline, replacement)
        self.assertEqual(result["status"], SampleStatus.MISSING_FIELD.value)
        self.assertIsNone(result["full_replacement"])
        self.assertEqual(result["unavailable_factor"], "alpha_law")

        replacement["alpha_law"] = operand(
            "alpha_law", math.inf, status=SampleStatus.VALID
        )
        with self.assertRaisesRegex(ValueError, "finite"):
            run_replacement_matrix(baseline, replacement)

    def test_metric_summary_uses_linear_p95_and_typed_valid_samples(self):
        rows = [
            {"status": SampleStatus.VALID, "error": value}
            for value in (0.0, 1.0, 2.0)
        ] + [
            {"status": SampleStatus.GEOMETRIC_ZERO, "error": 0.0},
            {"status": SampleStatus.MISSING_FIELD, "error": None},
            {"status": SampleStatus.NONFINITE, "error": None},
        ]
        summary = metric_summary(rows)
        self.assertEqual(summary["valid_count"], 3)
        self.assertEqual(summary["status_counts"]["geometric_zero"], 1)
        self.assertAlmostEqual(summary["median_abs_error"], 1.0)
        self.assertAlmostEqual(summary["p95_abs_error"], 1.9)
        with self.assertRaisesRegex(ValueError, "valid.*finite"):
            metric_summary([{"status": SampleStatus.VALID, "error": math.nan}])

    def test_identified_confounded_insufficient_and_holdout_rejected(self):
        passing = [
            evidence("passing", "discovery", 0.02, prediction=1.0),
            evidence("passing", "holdout", 0.03, prediction=1.1),
        ]
        self.assertIs(
            classify_candidate("passing", passing), Identifiability.IDENTIFIED
        )

        confounded = [
            evidence("combined", "discovery", 0.02, prediction=1.0,
                     missing_independent_factors=("mobility",)),
            evidence("combined", "holdout", 0.03, prediction=1.1,
                     missing_independent_factors=("mobility",)),
        ]
        self.assertIs(
            classify_candidate("combined", confounded), Identifiability.CONFOUNDED
        )

        missing = [
            evidence("missing", "discovery", None, prediction=None,
                     status=SampleStatus.MISSING_FIELD),
            evidence("missing", "holdout", None, prediction=None,
                     status=SampleStatus.MISSING_FIELD),
        ]
        self.assertIs(
            classify_candidate("missing", missing),
            Identifiability.INSUFFICIENT_DATA,
        )

        holdout_failure = [
            evidence("overfit", "discovery", 0.02, prediction=1.0),
            evidence("overfit", "holdout", 0.4, prediction=1.1),
        ]
        self.assertIs(
            classify_candidate("overfit", holdout_failure),
            Identifiability.REJECTED,
        )

    def test_indistinguishable_passing_candidates_are_consistent_nonunique(self):
        rows = []
        for split, prediction in (("discovery", 1.0), ("holdout", 1.2)):
            rows.extend([
                evidence("a", split, 0.02, prediction=prediction),
                evidence("b", split, 0.03, prediction=prediction + 5.0e-11),
            ])
        self.assertIs(
            classify_candidate("a", rows), Identifiability.CONSISTENT_NONUNIQUE
        )
        ranked = rank_candidates(rows)
        self.assertEqual([row["candidate"] for row in ranked], ["a", "b"])
        self.assertTrue(all(row["classification"] is Identifiability.CONSISTENT_NONUNIQUE
                            for row in ranked))

    def test_ranking_uses_discovery_only_and_rejects_local_fit_leakage(self):
        rows = [
            evidence("discovery_best", "discovery", 0.01, prediction=1.0),
            evidence("discovery_best", "holdout", 0.09, prediction=1.2),
            evidence("holdout_best", "discovery", 0.02, prediction=2.0),
            evidence("holdout_best", "holdout", 0.01, prediction=2.2),
            evidence("bias_fit", "discovery", 0.0, prediction=3.0,
                     fit_dimensions=("bias",)),
            evidence("bias_fit", "holdout", 0.0, prediction=3.2,
                     fit_dimensions=("bias",)),
        ]
        ranked = rank_candidates(rows)
        self.assertEqual(
            [row["candidate"] for row in ranked],
            ["discovery_best", "holdout_best", "bias_fit"],
        )
        self.assertIs(ranked[-1]["classification"], Identifiability.REJECTED)
        self.assertIn("fit", ranked[-1]["reason"])

    def test_closure_gate_is_fixed_even_if_threshold_object_is_relaxed(self):
        rows = [
            evidence("loose", "discovery", 5.0e-10, prediction=1.0,
                     metric="replacement_closure_abs_dex"),
            evidence("loose", "holdout", 5.0e-10, prediction=1.2,
                     metric="replacement_closure_abs_dex"),
        ]
        relaxed = AcceptanceThresholds(replacement_closure_abs_dex=1.0)
        self.assertIs(
            classify_candidate("loose", rows, thresholds=relaxed),
            Identifiability.REJECTED,
        )

    def test_whole_state_replay_is_localization_only_and_excluded(self):
        baseline = {
            "potential": operand("potential", 1.0, carrier=None, unit_si="V"),
            "carrier_state": operand("carrier_state", 2.0, unit_si="m^-3"),
            "quasi_fermi_state": operand("quasi_fermi_state", 3.0, unit_si="V"),
        }
        replacement = {
            "potential": operand("potential", 2.0, carrier=None, unit_si="V"),
            "carrier_state": operand("carrier_state", 4.0, unit_si="m^-3"),
            "quasi_fermi_state": operand("quasi_fermi_state", 5.0, unit_si="V"),
        }
        control = run_state_localization_control(
            baseline, replacement,
            evaluate=lambda state: state["potential"] * state["carrier_state"]
                                   * state["quasi_fermi_state"],
        )
        self.assertEqual(control["classification"], "localization_control")
        self.assertFalse(control["eligible_for_candidate_ranking"])
        self.assertFalse(control["eligible_for_formula_classification"])
        self.assertEqual(control["layers_replaced"], [
            "potential", "carrier_state", "quasi_fermi_state"
        ])
        self.assertEqual(rank_candidates([control]), [])
        with self.assertRaisesRegex(ValueError, "localization"):
            classify_candidate(control["candidate"], [control])


if __name__ == "__main__":
    unittest.main()
