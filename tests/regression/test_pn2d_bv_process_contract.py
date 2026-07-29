#!/usr/bin/env python3

from __future__ import annotations

import copy
import hashlib
import json
import math
import tempfile
import unittest
from pathlib import Path

from scripts.pn2d_bv_process_contract import (
    SCHEMA_ID,
    ProcessContractError,
    normalized_records_sha256,
    validate_process_run,
)


ZERO_HASH = hashlib.sha256(b"").hexdigest()


def source(index: int = 0) -> dict[str, object]:
    return {"file": "raw/run.out", "dataset": "probe", "index": index}


def field(
    *,
    support_kind: str = "physical_node",
    support_key: str = "node:0",
    centering: str = "vertex",
    provenance: str = "native",
    carrier: str = "electron",
    quantity: str = "density",
    unit: str = "cm^-3",
    value: float = 1.0,
) -> dict[str, object]:
    return {
        "branch": "avalanche_off",
        "requested_bias_V": -19.0,
        "actual_bias_V": -19.0,
        "support_kind": support_kind,
        "support_key": support_key,
        "centering": centering,
        "provenance": provenance,
        "carrier": carrier,
        "quantity": quantity,
        "components": ["scalar"],
        "unit": unit,
        "values": [value],
        "coordinates_um": [0.0, 0.0],
        "source": source(),
    }


def valid_manifest() -> dict[str, object]:
    artifact = {"path": "raw/empty", "sha256": ZERO_HASH}
    return {
        "schema": SCHEMA_ID,
        "status": "passed",
        "outcome": "process_contract_verified",
        "run_id": "fixture",
        "simulator": "sentaurus",
        "release": "O-2018.06-SP2",
        "missing_value_policy": "reject",
        "input_hashes": {"input.dat": ZERO_HASH},
        "normalized_output_hashes": {"normalized/aggregate.csv": ZERO_HASH},
        "branch_records": [
            {
                "branch": "avalanche_off",
                "requested_biases_V": [-19.0],
                "bias_records": [
                    {
                        "requested_bias_V": -19.0,
                        "actual_bias_V": -19.0,
                        "snapshot_tdr": copy.deepcopy(artifact),
                        "currentplot": copy.deepcopy(artifact),
                        "process_record": copy.deepcopy(artifact),
                    }
                ],
            }
        ],
        "field_records": [field()],
        "aggregate_records": [
            {
                "branch": "avalanche_off",
                "requested_bias_V": -19.0,
                "actual_bias_V": -19.0,
                "carrier": "total",
                "quantity": "terminal_current",
                "unit": "A/um",
                "value": 1.0e-12,
                "provenance": "native",
                "source": source(),
            }
        ],
        "newton_attempt_records": [
            {
                "branch": "avalanche_off",
                "attempt_id": "segment-0-attempt-0",
                "requested_bias_V": -19.0,
                "actual_bias_V": -19.0,
                "status": "accepted",
                "reason": "converged",
                "source": source(),
            }
        ],
    }


class BVProcessContractTest(unittest.TestCase):
    def assert_reason(self, manifest: dict[str, object], reason: str) -> None:
        with self.assertRaises(ProcessContractError) as captured:
            validate_process_run(manifest)
        self.assertEqual(captured.exception.reason, reason)

    def test_schema_matches_runtime_contract(self) -> None:
        schema = json.loads(
            (
                Path(__file__).resolve().parents[2]
                / "schemas"
                / "vela.pn2d_bv_process_run.v1.schema.json"
            ).read_text(encoding="ascii")
        )
        self.assertEqual(schema["title"], SCHEMA_ID)
        self.assertEqual(schema["properties"]["schema"]["const"], SCHEMA_ID)
        self.assertEqual(schema["properties"]["missing_value_policy"]["const"], "reject")

    def test_minimal_valid_paired_run(self) -> None:
        manifest = valid_manifest()
        paired = copy.deepcopy(manifest["branch_records"][0])
        paired["branch"] = "iic_postprocess"
        manifest["branch_records"].append(paired)
        validate_process_run(manifest)

    def test_contact_support_duplicate_is_distinct_from_physical_node(self) -> None:
        manifest = valid_manifest()
        contact = field(
            support_kind="contact_support_vertex",
            support_key="node:0",
            centering="vertex",
        )
        manifest["field_records"].append(contact)
        validate_process_run(manifest)

    def test_duplicate_support_key_fails(self) -> None:
        manifest = valid_manifest()
        manifest["field_records"].append(copy.deepcopy(manifest["field_records"][0]))
        self.assert_reason(manifest, "duplicate_support_key")

    def test_duplicate_and_missing_bias_rows_fail(self) -> None:
        duplicate = valid_manifest()
        duplicate["branch_records"][0]["bias_records"].append(
            copy.deepcopy(duplicate["branch_records"][0]["bias_records"][0])
        )
        self.assert_reason(duplicate, "duplicate_bias_row")

        missing = valid_manifest()
        missing["branch_records"][0]["requested_biases_V"].append(-20.0)
        self.assert_reason(missing, "missing_bias_row")

    def test_wrong_centering_and_unit_fail(self) -> None:
        wrong_centering = valid_manifest()
        wrong_centering["field_records"][0]["centering"] = "cell"
        self.assert_reason(wrong_centering, "wrong_centering")

        wrong_unit = valid_manifest()
        wrong_unit["field_records"][0]["unit"] = "V"
        self.assert_reason(wrong_unit, "wrong_unit")

        unknown_unit = valid_manifest()
        unknown_unit["field_records"][0]["unit"] = "arbitrary"
        self.assert_reason(unknown_unit, "unknown_unit")

    def test_support_remapping_data_is_required(self) -> None:
        missing_coordinates = valid_manifest()
        missing_coordinates["field_records"][0].pop("coordinates_um")
        self.assert_reason(missing_coordinates, "missing_support_coordinates")

        missing_connectivity = valid_manifest()
        missing_connectivity["field_records"] = [
            field(
                support_kind="cell",
                support_key="cell:0",
                centering="cell",
            )
        ]
        self.assert_reason(missing_connectivity, "missing_support_connectivity")

    def test_unsupported_native_edge_current_claim_fails(self) -> None:
        manifest = valid_manifest()
        manifest["field_records"] = [
            field(
                support_kind="element_local_edge",
                support_key="cell:0/edge:0",
                centering="element_edge",
                quantity="current_density",
                unit="A/cm^2",
            )
        ]
        manifest["field_records"][0]["connectivity"] = [0, 1]
        self.assert_reason(manifest, "unsupported_native_edge_claim")

    def test_unknown_nested_fields_fail_closed(self) -> None:
        manifest = valid_manifest()
        manifest["field_records"][0]["undeclared"] = "value"
        self.assert_reason(manifest, "unexpected_record_field")

    def test_nonfinite_missing_carrier_and_implicit_fill_fail(self) -> None:
        nonfinite = valid_manifest()
        nonfinite["field_records"][0]["values"] = [math.nan]
        self.assert_reason(nonfinite, "nonfinite_value")

        missing_carrier = valid_manifest()
        missing_carrier["field_records"][0]["carrier"] = ""
        self.assert_reason(missing_carrier, "missing_or_unknown_carrier")

        implicit_fill = valid_manifest()
        implicit_fill["field_records"][0]["values"] = []
        self.assert_reason(implicit_fill, "implicit_zero_fill")

    def test_nearest_bias_substitution_fails(self) -> None:
        manifest = valid_manifest()
        manifest["branch_records"][0]["bias_records"][0]["actual_bias_V"] = -18.999
        self.assert_reason(manifest, "nearest_bias_substitution")

    def test_hash_drift_fails(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / "raw").mkdir()
            (root / "normalized").mkdir()
            for relative in ("input.dat", "raw/empty", "normalized/aggregate.csv"):
                (root / relative).write_bytes(b"")
            manifest = valid_manifest()
            validate_process_run(manifest, base_dir=root)
            (root / "input.dat").write_bytes(b"drift")
            self.assert_reason_with_base(manifest, "hash_drift", root)

    def assert_reason_with_base(
        self,
        manifest: dict[str, object],
        reason: str,
        base_dir: Path,
    ) -> None:
        with self.assertRaises(ProcessContractError) as captured:
            validate_process_run(manifest, base_dir=base_dir)
        self.assertEqual(captured.exception.reason, reason)

    def test_reordered_rows_have_same_normalized_hash(self) -> None:
        rows = [
            {"support": "node:1", "value": 2.0},
            {"support": "node:0", "value": 1.0},
        ]
        self.assertEqual(
            normalized_records_sha256(rows),
            normalized_records_sha256(list(reversed(rows))),
        )


if __name__ == "__main__":
    unittest.main()
