#!/usr/bin/env python3

from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from scripts.pn2d_high_bias_process_contract import EXACT_HIGH_BIAS_V
from scripts.pn2d_sentaurus_process_run_contract import (
    SCHEMA_ID,
    build_run_manifest,
    validate_case,
    validate_run_manifest,
)
from scripts.run_pn2d_high_bias_oracle_matrix import passed


class SentaurusProcessRunContractTest(unittest.TestCase):
    def test_versioned_schema_matches_runtime_contract(self) -> None:
        schema = json.loads(
            (
                Path(__file__).resolve().parents[2]
                / "schemas"
                / "vela.pn2d_sentaurus_process_run.v1.schema.json"
            ).read_text(encoding="ascii")
        )
        self.assertEqual(schema["title"], SCHEMA_ID)
        self.assertEqual(schema["properties"]["schema"]["const"], SCHEMA_ID)
        self.assertIn("passed", schema["properties"]["status"]["enum"])

    def test_case_validation_closes_manifest_file_hashes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            case = Path(tmp)
            bundle = case / "bundle"
            fetched = case / "fetched"
            bundle.mkdir()
            fetched.mkdir()
            (bundle / "deck.cmd").write_text("deck\n", encoding="ascii")
            (fetched / "run.out").write_text("run\n", encoding="ascii")
            manifest = build_run_manifest(
                status="passed",
                experiment="pn2d_high_bias_process_probe",
                variant="implicit_default",
                exact_biases=(-19.95,),
                observed_biases=(-19.95,),
                remote_root="/remote/task13",
                bundle=bundle,
                fetched=fetched,
            )
            (case / "manifest.json").write_text(
                json.dumps(manifest),
                encoding="ascii",
            )
            validate_case(
                case,
                experiment="pn2d_high_bias_process_probe",
                variant="implicit_default",
                exact_biases=(-19.95,),
            )
            (fetched / "run.out").write_text("tampered\n", encoding="ascii")
            with self.assertRaisesRegex(ValueError, "hash closure mismatch"):
                validate_case(
                    case,
                    experiment="pn2d_high_bias_process_probe",
                    variant="implicit_default",
                    exact_biases=(-19.95,),
                )

    def test_jacobian_schema_cannot_masquerade_as_run_manifest(self) -> None:
        manifest = {
            "schema": "vela.pn2d_high_bias_process_jacobian.v1",
            "status": "passed",
        }
        with self.assertRaisesRegex(ValueError, "schema mismatch"):
            validate_run_manifest(
                manifest,
                experiment="pn2d_high_bias_process_probe",
                variant="implicit_default",
                exact_biases=(-19.95,),
            )

    def test_matrix_resume_rejects_stale_or_tampered_case(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            case = Path(tmp)
            bundle = case / "bundle"
            fetched = case / "fetched"
            bundle.mkdir()
            fetched.mkdir()
            (bundle / "deck.cmd").write_text("deck\n", encoding="ascii")
            (fetched / "run.out").write_text("run\n", encoding="ascii")
            manifest = build_run_manifest(
                status="passed",
                experiment="pn2d_exact_high_bias_oracle_variant",
                variant="implicit_default",
                exact_biases=EXACT_HIGH_BIAS_V,
                observed_biases=EXACT_HIGH_BIAS_V,
                remote_root="/remote/task14",
                bundle=bundle,
                fetched=fetched,
            )
            (case / "manifest.json").write_text(
                json.dumps(manifest),
                encoding="ascii",
            )
            self.assertTrue(passed(case, "implicit_default"))
            (fetched / "run.out").write_text("stale\n", encoding="ascii")
            self.assertFalse(passed(case, "implicit_default"))


if __name__ == "__main__":
    unittest.main()
