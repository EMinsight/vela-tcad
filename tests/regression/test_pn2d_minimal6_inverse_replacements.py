import math
import unittest

from scripts.pn2d_minimal6_diagnostics.inverse_contracts import (
    AcceptanceThresholds,
    Identifiability,
    SampleStatus,
    SupportKind,
)
from scripts.pn2d_minimal6_diagnostics.inverse_avalanche import (
    AvalancheCandidateResult,
    AvalancheCandidateSample,
    GenerationError,
    GenerationMetricSummary,
)
from scripts.pn2d_minimal6_diagnostics.inverse_fields import (
    FieldCandidateResult,
    FieldCandidateSample,
    VectorErrorResult,
)
from scripts.pn2d_minimal6_diagnostics.inverse_replacements import (
    INVERSE_DEPENDENCIES,
    classify_candidate,
    metric_summary,
    rank_candidates,
    run_replacement_matrix,
    run_state_localization_control,
)
from scripts.pn2d_minimal6_diagnostics.inverse_transport import (
    TransportCandidateResult,
    TransportCandidateSample,
    TransportConfoundingRecord,
    TransportVectorError,
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
                    run_replacement_matrix(baseline, replacement, direct_target=128.0)

    def test_invalid_replacement_values_remain_typed_and_are_not_coerced(self):
        baseline = {factor: operand(factor, 1.0) for factor in FACTORS}
        replacement = {factor: operand(factor, 2.0) for factor in FACTORS}
        replacement["alpha_law"] = operand(
            "alpha_law", None, status=SampleStatus.MISSING_FIELD
        )
        result = run_replacement_matrix(baseline, replacement, direct_target=128.0)
        self.assertEqual(result["status"], SampleStatus.MISSING_FIELD.value)
        self.assertIsNone(result["full_replacement"])
        self.assertEqual(result["unavailable_factor"], "alpha_law")

        replacement["alpha_law"] = operand(
            "alpha_law", math.inf, status=SampleStatus.VALID
        )
        with self.assertRaisesRegex(ValueError, "finite"):
            run_replacement_matrix(baseline, replacement, direct_target=128.0)

    def test_direct_target_is_independent_required_and_exact(self):
        baseline = {factor: operand(factor, 1.0) for factor in FACTORS}
        replacement = {factor: operand(factor, 2.0) for factor in FACTORS}
        with self.assertRaisesRegex(ValueError, "direct_target.*required"):
            run_replacement_matrix(baseline, replacement)
        with self.assertRaisesRegex(ValueError, "closure"):
            run_replacement_matrix(baseline, replacement, direct_target=127.0)

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
            evidence("combined", "discovery", None, prediction=None,
                     status=SampleStatus.MISSING_FIELD,
                     missing_independent_factors=("mobility",)),
            evidence("combined", "holdout", None, prediction=None,
                     status=SampleStatus.MISSING_FIELD,
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
            evidence("bias_fit", "discovery", 0.03, prediction=3.0,
                     fit_dimensions=("bias",)),
            evidence("bias_fit", "holdout", 0.03, prediction=3.2,
                     fit_dimensions=("bias",)),
        ]
        ranked = rank_candidates(rows)
        self.assertEqual(
            [row["candidate"] for row in ranked],
            ["discovery_best", "holdout_best", "bias_fit"],
        )
        self.assertIs(ranked[-1]["classification"], Identifiability.REJECTED)
        self.assertIn("fit", ranked[-1]["reason"])

    def test_holdout_rejection_does_not_reorder_discovery_ranking(self):
        rows = [
            evidence("discovery_best", "discovery", 0.01, prediction=1.0),
            evidence("discovery_best", "holdout", 0.4, prediction=1.2),
            evidence("holdout_pass", "discovery", 0.02, prediction=2.0),
            evidence("holdout_pass", "holdout", 0.02, prediction=2.2),
        ]
        ranked = rank_candidates(rows)
        self.assertEqual(
            [row["candidate"] for row in ranked],
            ["discovery_best", "holdout_pass"],
        )
        self.assertIs(ranked[0]["classification"], Identifiability.REJECTED)
        self.assertIs(ranked[1]["classification"], Identifiability.IDENTIFIED)

    def test_replacement_matrix_requires_one_global_context_but_not_one_unit(self):
        baseline = {
            factor: operand(factor, 1.0, unit_si=f"unit:{factor}")
            for factor in FACTORS
        }
        replacement = {
            factor: operand(factor, 2.0, unit_si=f"unit:{factor}")
            for factor in FACTORS
        }
        self.assertEqual(
            run_replacement_matrix(
                baseline, replacement, direct_target=128.0
            )["status"],
            SampleStatus.VALID.value,
        )
        for field, value in (
            ("topology", "mirror"), ("bias_V", -19.0),
            ("carrier", "hole"), ("support_kind", SupportKind.EDGE),
            ("support_id", "edge-7"),
        ):
            with self.subTest(field=field):
                mixed_baseline = dict(baseline)
                mixed_replacement = dict(replacement)
                mixed_baseline["mobility"] = operand(
                    "mobility", 1.0, unit_si="unit:mobility", **{field: value}
                )
                mixed_replacement["mobility"] = operand(
                    "mobility", 2.0, unit_si="unit:mobility", **{field: value}
                )
                with self.assertRaisesRegex(ValueError, f"global.*{field}"):
                    run_replacement_matrix(
                        mixed_baseline, mixed_replacement, direct_target=128.0
                    )

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

        for field, value in (("bias_V", -19.0), ("support_id", 8),
                             ("carrier", "hole")):
            with self.subTest(field=field):
                mixed_baseline = dict(baseline)
                mixed_replacement = dict(replacement)
                mixed_baseline["carrier_state"] = operand(
                    "carrier_state", 2.0, unit_si="m^-3", **{field: value}
                )
                mixed_replacement["carrier_state"] = operand(
                    "carrier_state", 4.0, unit_si="m^-3", **{field: value}
                )
                with self.assertRaisesRegex(ValueError, f"global.*{field}"):
                    run_state_localization_control(
                        mixed_baseline, mixed_replacement,
                        evaluate=lambda state: sum(state.values()),
                    )


        unavailable_baseline = dict(baseline)
        unavailable_replacement = dict(replacement)
        unavailable_baseline["carrier_state"] = operand(
            "carrier_state", None, unit_si="m^-3", bias_V=-19.0,
            status=SampleStatus.MISSING_FIELD,
        )
        unavailable_replacement["carrier_state"] = operand(
            "carrier_state", None, unit_si="m^-3", bias_V=-19.0,
            status=SampleStatus.MISSING_FIELD,
        )
        with self.assertRaisesRegex(ValueError, "global.*bias_V"):
            run_state_localization_control(
                unavailable_baseline, unavailable_replacement,
                evaluate=lambda state: sum(state.values()),
            )

    def test_real_task4_task5_and_task6_results_are_flattened(self):
        def field_result(topology):
            error = VectorErrorResult(
                SampleStatus.VALID, SampleStatus.VALID, 0.01, 0.0
            )
            sample = FieldCandidateSample(
                "triangle_minus_grad_psi", "sentaurus", topology, -12.0,
                SupportKind.NODE, 7, (10.0, 0.0), (10.0, 0.0), "V/m",
                "sentaurus_xy", "global_vector", "triangle_to_node",
                "native_node", error,
            )
            return FieldCandidateResult(
                sample.candidate, sample.solver, sample.topology, sample.bias_V,
                sample.support_kind, (sample,), 1, 1, 0.01, 0.0,
                Identifiability.IDENTIFIED,
            )

        def transport_result(split, topology):
            error = TransportVectorError(
                SampleStatus.VALID, SampleStatus.VALID, 0.01, 0.0
            )
            sample = TransportCandidateSample(
                "node_area_weighted_qf_gradient_current", "sentaurus",
                topology, -12.0, "electron", split, SupportKind.NODE, 7,
                (10.0, 0.0), (10.0, 0.0), "A/m^2", "sentaurus_xy",
                "global_vector", "cell_to_node", "native_node", error,
            )
            return TransportCandidateResult(
                sample.candidate, sample.solver, sample.topology, sample.bias_V,
                sample.carrier, sample.split, sample.support_kind, sample.unit_si,
                (sample,), 1, 1, 0.01, 0.01, 0.0,
                Identifiability.IDENTIFIED, (),
            )

        def avalanche_result(split, topology):
            error = GenerationError(SampleStatus.VALID, 0.01)
            sample = AvalancheCandidateSample(
                candidate="electric_field_magnitude", solver="sentaurus",
                topology=topology, bias_V=-12.0, split=split,
                support_kind=SupportKind.NODE, support_id=7,
                electron_driver_V_m=1.0, hole_driver_V_m=1.0,
                electron_reference_driver_V_m=1.0,
                hole_reference_driver_V_m=1.0,
                electron_alpha_m_inv=1.0, hole_alpha_m_inv=1.0,
                candidate_generation_m3_s=10.0,
                reference_generation_m3_s=10.0, error=error,
            )
            summary = GenerationMetricSummary(
                (error,), error, 1, 0.01, 0.01, 0.01,
                Identifiability.IDENTIFIED,
            )
            return AvalancheCandidateResult(
                sample.candidate, sample.solver, sample.topology, sample.bias_V,
                sample.split, sample.support_kind, (sample,), None, summary, ()
            )

        cases = (
            ("triangle_minus_grad_psi",
             [field_result("sketch"), field_result("mirror")]),
            ("node_area_weighted_qf_gradient_current",
             [transport_result("discovery", "sketch"),
              transport_result("holdout", "mirror")]),
            ("electric_field_magnitude",
             [avalanche_result("discovery", "sketch"),
              avalanche_result("holdout", "mirror")]),
        )
        for candidate, records in cases:
            with self.subTest(candidate=candidate):
                self.assertGreater(metric_summary(records)["valid_count"], 0)
                self.assertIs(
                    classify_candidate(candidate, records),
                    Identifiability.IDENTIFIED,
                )
                self.assertEqual(rank_candidates(records)[0]["candidate"], candidate)

    def test_real_transport_missing_mobility_is_confounded_before_support_count(self):
        def result(split, topology):
            confounding = TransportConfoundingRecord(
                "current_inverted_qf_gradient", "electron", topology, -12.0,
                SupportKind.NODE, 7, SampleStatus.MISSING_FIELD,
                Identifiability.CONFOUNDED, ("mobility",),
                "mu_times_gradient", (1.0, 0.0), "m*s^-1",
            )
            return TransportCandidateResult(
                confounding.candidate, "sentaurus", topology, -12.0,
                confounding.carrier, split, confounding.support_kind, "V/m",
                (), 0, 0, None, None, None, Identifiability.CONFOUNDED,
                (confounding,),
            )

        records = [result("discovery", "sketch"),
                   result("holdout", "mirror")]
        self.assertIs(
            classify_candidate("current_inverted_qf_gradient", records),
            Identifiability.CONFOUNDED,
        )

if __name__ == "__main__":
    unittest.main()
