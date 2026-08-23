from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT / "scripts") not in sys.path:
    sys.path.insert(0, str(ROOT / "scripts"))

from compare_bvmethods_release_summaries import compare  # noqa: E402


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class BvSchottkyFrozenLedgersTest(unittest.TestCase):
    def test_bvmethods_release_comparator_detects_drift_and_missing_methods(self) -> None:
        baseline = {"results": [
            {"method": "resistor", "bv_V": 6.0, "rows": 10, "criterion": "I"},
            {"method": "continuation", "bv_V": 6.1, "rows": 11, "criterion": "I"},
        ]}
        candidate = {"results": [
            {"method": "resistor", "bv_V": 6.03, "rows": 12, "criterion": "I"},
        ]}
        result = compare(baseline, candidate, max_relative_drift=0.01)
        self.assertEqual("fail", result["status"])
        methods = {row["method"]: row for row in result["methods"]}
        self.assertEqual("pass", methods["resistor"]["status"])
        self.assertEqual("missing", methods["continuation"]["status"])

    def test_bvmethods_nontransient_ledger_and_curves_are_sealed(self) -> None:
        fixture = ROOT / "reference_tcad" / "bvmethods_sentaurus2018"
        ledger = json.loads((
            fixture / "bvmethods_nontransient_validation_20260817.json"
        ).read_text(encoding="utf-8"))
        self.assertEqual(
            ["ABA_poisson", "ABA_coupled", "resistor", "voltage2current", "continuation"],
            ledger["scope"],
        )
        self.assertEqual(["transient"], ledger["excluded"])
        self.assertEqual("pass", ledger["sentaurus_reference"]["status"])
        for method in ledger["sentaurus_reference"]["methods"]:
            self.assertTrue(method["reference_matches"], method["method"])
            curve = fixture / method["curve"]
            self.assertTrue(curve.is_file())
            self.assertEqual(
                ledger["sha256"][method["curve"]], sha256(curve))
        self.assertIn(
            ledger["status"], {"pass", "pass_with_continuation_pending"})
        continuation = ledger["vela_acceptance"]["continuation"]
        if continuation["status"] == "pending":
            diagnostic = fixture / continuation["diagnostic"]
            self.assertTrue(diagnostic.is_file())
            self.assertEqual(
                ledger["sha256"][continuation["diagnostic"]],
                sha256(diagnostic),
            )
            evidence = json.loads(diagnostic.read_text(encoding="utf-8"))
            self.assertEqual("none", evidence["physics_changes"])

    def test_simple_schottky_clean_workflow_is_passed(self) -> None:
        ledger = json.loads((
            ROOT / "reference_tcad" / "schottky_charon_sentaurus2018"
            / "schottky_workflow_validation_20260817.json"
        ).read_text(encoding="utf-8"))
        self.assertEqual("pass", ledger["status"])
        self.assertEqual(24, ledger["curve_acceptance"]["points_compared"])
        self.assertLessEqual(
            ledger["curve_acceptance"]["maximum_log10_current_error_dex"], 0.5)
        self.assertEqual("pass", ledger["stage_b"]["status"])
        self.assertEqual(
            [0.0, 0.4, 1.0],
            [point["requested_bias_V"]
             for point in ledger["three_bias_boundary_kcl"]["points"]],
        )
        self.assertTrue(all(
            point["status"] == "pass"
            for point in ledger["three_bias_boundary_kcl"]["points"]))


if __name__ == "__main__":
    unittest.main()
