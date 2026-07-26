#!/usr/bin/env python3
"""Contracts for PN2D high-bias process and source-Jacobian evidence."""

from __future__ import annotations

import copy
import json
import unittest
from pathlib import Path

from scripts.pn2d_high_bias_process_contract import (
    EXACT_HIGH_BIAS_V,
    REQUIRED_DERIVATIVE_STAGES,
    SCHEMA_ID,
    true_relative_error,
    validate_derivative_lattice,
    validate_manifest_pair,
    validate_process_record,
)


REPO = Path(__file__).resolve().parents[2]
SCHEMA_PATH = (
    REPO / "schemas" / "vela.pn2d_high_bias_process_jacobian.v1.schema.json"
)
HASH = "a" * 64


def valid_record(stage: str = REQUIRED_DERIVATIVE_STAGES[0]) -> dict:
    return {
        "schema": SCHEMA_ID,
        "topology": "coarse7x3",
        "bias_V": -19.95,
        "carrier": "electron",
        "residual_variable": "electron_continuity",
        "state_variable": "electron_qfp",
        "process_stage": stage,
        "analytic_contribution": 0.1,
        "fd_contribution": 0.1,
        "analytic_total": 1.0,
        "fd_total": 1.0,
        "derivative_status": "finite",
        "observation_label": "operator_replay",
        "observation_provenance": "reconstructed",
        "support_status": "valid",
        "error_dex": None,
        "unit": "A/V",
        "fd_steps_V": [1.0e-6, 3.0e-7, 1.0e-7],
        "residual_config_sha256": HASH,
        "jacobian_config_sha256": HASH,
        "state_sha256": HASH,
        "deck_sha256": HASH,
        "tdr_sha256": HASH,
        "mesh_sha256": HASH,
        "parameters_sha256": HASH,
        "sentaurus_release": "O-2018.06-SP2",
    }


def valid_lattice() -> list[dict]:
    return [valid_record(stage) for stage in REQUIRED_DERIVATIVE_STAGES]


def valid_manifest() -> dict:
    return {
        "schema": SCHEMA_ID,
        "status": "red_contract_frozen",
        "sentaurus_release": "O-2018.06-SP2",
        "exact_biases_V": list(EXACT_HIGH_BIAS_V),
        "input_hashes": {
            "deck": HASH,
            "tdr": HASH,
            "mesh": HASH,
            "parameters": HASH,
            "state": HASH,
        },
        "residual_config_sha256": HASH,
        "jacobian_config_sha256": HASH,
    }


class HighBiasProcessContractTest(unittest.TestCase):
    def test_versioned_schema_document_is_frozen(self) -> None:
        document = json.loads(SCHEMA_PATH.read_text(encoding="ascii"))
        self.assertEqual(document["title"], SCHEMA_ID)
        self.assertEqual(document["properties"]["schema"]["const"], SCHEMA_ID)
        self.assertIn("records", document["required"])
        self.assertFalse(document["additionalProperties"])

    def test_valid_record_lattice_and_manifest_pair_pass(self) -> None:
        for record in valid_lattice():
            validate_process_record(record)
        validate_derivative_lattice(valid_lattice())
        validate_manifest_pair(valid_manifest(), valid_manifest())

    def test_missing_dependency_contribution_is_rejected(self) -> None:
        records = valid_lattice()
        records.pop()
        with self.assertRaisesRegex(ValueError, "missing derivative contribution"):
            validate_derivative_lattice(records)

    def test_residual_and_jacobian_config_mismatch_is_rejected(self) -> None:
        record = valid_record()
        record["jacobian_config_sha256"] = "b" * 64
        with self.assertRaisesRegex(ValueError, "residual/Jacobian configuration"):
            validate_process_record(record)

    def test_reconstructed_current_cannot_be_labeled_native(self) -> None:
        record = valid_record()
        record["observation_label"] = "native_directed_edge"
        with self.assertRaisesRegex(ValueError, "reconstructed current mislabeled as native"):
            validate_process_record(record)

    def test_exact_zero_cannot_be_converted_to_finite_dex(self) -> None:
        record = valid_record()
        record["support_status"] = "zero"
        record["error_dex"] = 0.0
        with self.assertRaisesRegex(ValueError, "zero converted to finite dex"):
            validate_process_record(record)

    def test_true_relative_denominator_does_not_use_absolute_floor_one(self) -> None:
        self.assertEqual(true_relative_error(1.0e-15, 2.0e-15), 0.5)
        self.assertEqual(true_relative_error(0.0, 0.0), 0.0)

    def test_nonexact_bias_is_rejected(self) -> None:
        record = valid_record()
        record["bias_V"] = -19.94
        with self.assertRaisesRegex(ValueError, "non-exact bias"):
            validate_process_record(record)

    def test_manifest_provenance_hash_mismatch_is_rejected(self) -> None:
        changed = copy.deepcopy(valid_manifest())
        changed["input_hashes"]["tdr"] = "b" * 64
        with self.assertRaisesRegex(ValueError, "input provenance mismatch"):
            validate_manifest_pair(valid_manifest(), changed)


if __name__ == "__main__":
    unittest.main()
