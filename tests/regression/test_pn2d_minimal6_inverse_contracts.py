import math
from dataclasses import FrozenInstanceError
import unittest

from scripts.pn2d_minimal6_diagnostics.inverse_contracts import (
    AcceptanceThresholds,
    Identifiability,
    Observation,
    SampleStatus,
    SupportKind,
    classify_numeric_sample,
    validate_inverse_report_v1,
)


class InverseContractsTest(unittest.TestCase):
    def make_valid_report(self):
        return {
            "schema": "vela.pn2d_minimal6_physics_inverse_audit.v1",
            "diagnostic_only": True,
            "phase_base": "a5524cf",
            "payload": {
                "input_manifest_sha256": "0" * 64,
                "discovery_keys": [["sketch", -1.0]],
                "holdout_keys": [["mirror", -1.0]],
                "thresholds": AcceptanceThresholds().__dict__,
                "field_inventory": {},
                "sample_status_counts": {"valid": 1},
                "candidate_metrics": [],
                "classifications": [
                    {
                        "candidate": "triangle_gradient",
                        "classification": Identifiability.IDENTIFIED.value,
                    }
                ],
                "replacement_closure": [],
                "localization_control": {},
                "sentaurus_version": "O-2018.06-SP2",
                "production_cpp_changed": False,
            },
        }

    def test_numeric_statuses_do_not_turn_missing_into_zero(self):
        self.assertEqual(
            classify_numeric_sample(None, floor=1e-30),
            SampleStatus.MISSING_FIELD,
        )
        self.assertEqual(
            classify_numeric_sample(0.0, floor=1e-30, geometric_zero=True),
            SampleStatus.GEOMETRIC_ZERO,
        )
        self.assertEqual(
            classify_numeric_sample(1e-40, floor=1e-30),
            SampleStatus.BELOW_FLOOR,
        )
        self.assertEqual(
            classify_numeric_sample(math.inf, floor=1e-30),
            SampleStatus.NONFINITE,
        )
        self.assertEqual(
            classify_numeric_sample(2.0, floor=1e-30),
            SampleStatus.VALID,
        )

    def test_observation_key_is_complete_and_immutable(self):
        row = Observation(
            "sentaurus",
            "sketch",
            -12.0,
            SupportKind.NODE,
            1,
            "electric_field",
            "x",
            -64036.5,
            "V*cm^-1",
            -6.40365e6,
            "V/m",
            "sentaurus_xy",
            "global_vector",
            "multiply_by_100",
            SampleStatus.VALID,
            "field.csv",
            "0" * 64,
        )
        self.assertEqual(
            row.key,
            ("sentaurus", "sketch", -12.0, "node", 1, "electric_field", "x"),
        )
        with self.assertRaises(FrozenInstanceError):
            row.value_si = 0.0

    def test_thresholds_and_report_schema_are_exact(self):
        limits = AcceptanceThresholds()
        self.assertEqual(limits.gradient_median_abs_dex, 0.1)
        self.assertEqual(limits.gradient_p95_abs_dex, 0.3)
        report = {
            "schema": "vela.pn2d_minimal6_physics_inverse_audit.v1",
            "diagnostic_only": True,
            "phase_base": "a5524cf",
            "payload": {
                "input_manifest_sha256": "0" * 64,
                "discovery_keys": [["sketch", -1.0]],
                "holdout_keys": [["mirror", -1.0]],
                "thresholds": limits.__dict__,
                "field_inventory": {},
                "sample_status_counts": {"valid": 1},
                "candidate_metrics": [],
                "classifications": [
                    {
                        "candidate": "triangle_gradient",
                        "classification": Identifiability.IDENTIFIED.value,
                    }
                ],
                "replacement_closure": [],
                "localization_control": {},
                "sentaurus_version": "O-2018.06-SP2",
                "production_cpp_changed": False,
            },
        }
        validate_inverse_report_v1(report)
    def test_report_rejects_extra_or_missing_contract_keys(self):
        cases = []
        extra_top_level = self.make_valid_report()
        extra_top_level["unexpected"] = True
        cases.append(extra_top_level)
        missing_top_level = self.make_valid_report()
        del missing_top_level["schema"]
        cases.append(missing_top_level)
        extra_payload = self.make_valid_report()
        extra_payload["payload"]["unexpected"] = True
        cases.append(extra_payload)
        missing_payload = self.make_valid_report()
        del missing_payload["payload"]["field_inventory"]
        cases.append(missing_payload)

        for report in cases:
            with self.subTest(report=report):
                with self.assertRaises(ValueError):
                    validate_inverse_report_v1(report)

    def test_report_rejects_wrong_provenance(self):
        non_diagnostic = self.make_valid_report()
        non_diagnostic["diagnostic_only"] = False
        wrong_phase = self.make_valid_report()
        wrong_phase["phase_base"] = "not-a5524cf"

        for report in (non_diagnostic, wrong_phase):
            with self.subTest(report=report):
                with self.assertRaises(ValueError):
                    validate_inverse_report_v1(report)

    def test_report_rejects_unknown_classifications_in_all_row_types(self):
        classifications_row = self.make_valid_report()
        classifications_row["payload"]["classifications"][0]["classification"] = "unknown"
        candidate_metric_row = self.make_valid_report()
        candidate_metric_row["payload"]["candidate_metrics"] = [
            {"candidate": "triangle_gradient", "classification": "unknown"}
        ]

        for report in (classifications_row, candidate_metric_row):
            with self.subTest(report=report):
                with self.assertRaises(ValueError):
                    validate_inverse_report_v1(report)

    def test_report_rejects_nonfinite_candidate_metrics(self):
        for value in (math.nan, math.inf, -math.inf):
            with self.subTest(value=value):
                report = self.make_valid_report()
                report["payload"]["candidate_metrics"] = [
                    {"candidate": "triangle_gradient", "median_abs_error": value}
                ]
                with self.assertRaises(ValueError):
                    validate_inverse_report_v1(report)
