from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "singledevice_closeout", ROOT / "scripts" / "close_singledevice_validation.py")
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class SingleDeviceCloseoutTest(unittest.TestCase):
    def rows(self, scale: float = 1.0) -> list[dict[str, float]]:
        return [
            {"gate_voltage_V": -0.5, "sentaurus_current_A_per_um": 1.0e-14,
             "vela_current_A_per_um": 0.88e-14 * scale},
            {"gate_voltage_V": -0.2, "sentaurus_current_A_per_um": 1.0e-10,
             "vela_current_A_per_um": 0.96e-10 * scale},
            {"gate_voltage_V": 0.1, "sentaurus_current_A_per_um": 1.0e-7,
             "vela_current_A_per_um": 0.98e-7 * scale},
            {"gate_voltage_V": 2.2, "sentaurus_current_A_per_um": 1.0e-3,
             "vela_current_A_per_um": 1.01e-3 * scale},
        ]

    def test_low_current_hybrid_is_narrow_and_fail_closed(self) -> None:
        row = {
            "gate_voltage_V": -0.5,
            "sentaurus_current_A_per_um": 1.0e-14,
            "vela_current_A_per_um": 0.88e-14,
            "absolute_error_A_per_um": 1.2e-15,
            "relative_error": 0.12,
            "log_error_dex": 0.0555,
        }
        accepted = MODULE.point_gate(row, MODULE.DEFAULT_POLICY)
        self.assertEqual("pass", accepted["status"])
        self.assertEqual("low_current_hybrid", accepted["acceptance_path"])

        row["sentaurus_current_A_per_um"] = 2.0e-13
        rejected = MODULE.point_gate(row, MODULE.DEFAULT_POLICY)
        self.assertEqual("fail", rejected["status"])

    def test_absolute_guard_still_applies_below_current_boundary(self) -> None:
        row = {
            "gate_voltage_V": -0.5,
            "sentaurus_current_A_per_um": 9.0e-14,
            "vela_current_A_per_um": 1.05e-13,
            "absolute_error_A_per_um": 1.5e-14,
            "relative_error": 1.5e-14 / 9.0e-14,
            "log_error_dex": 0.067,
        }
        self.assertEqual("fail", MODULE.point_gate(row, MODULE.DEFAULT_POLICY)["status"])

    def test_real_frozen_comparisons_close(self) -> None:
        fixture = ROOT / "reference_tcad" / "singledevice_sentaurus2018"
        runtime = (ROOT.parents[1] / "build-release" / "reference_tcad" /
                   "singledevice_sentaurus2018" / "reports" /
                   "self_consistent_idvg_20260815")
        if not runtime.exists():
            self.skipTest("ignored SingleDevice runtime artifacts are unavailable")
        reference = json.loads(
            (fixture / "singledevice_sentaurus2018_reference.json").read_text())
        policy = dict(MODULE.DEFAULT_POLICY)
        policy.update(reference.get("low_current_acceptance", {}))
        report = MODULE.closeout(
            MODULE.load_comparison(runtime / "lin_compare.csv"),
            MODULE.load_comparison(runtime / "sat_compare.csv"),
            reference, policy)
        self.assertEqual("pass", report["status"])
        self.assertEqual(
            "low_current_hybrid",
            report["branches"]["linear"]["ioff"]["acceptance_path"])
        self.assertTrue(all(report["checks"].values()))


if __name__ == "__main__":
    unittest.main()
