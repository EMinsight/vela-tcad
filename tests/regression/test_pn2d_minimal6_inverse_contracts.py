import math
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
        with self.assertRaises(Exception):
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
