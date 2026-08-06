#!/usr/bin/env python3
"""Unit coverage for the Sentaurus BVmethods result extractor."""

from __future__ import annotations

import importlib.util
import sys
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
MODULE_PATH = REPO / "scripts" / "analyze_sentaurus_bvmethods.py"
if str(MODULE_PATH.parent) not in sys.path:
    sys.path.insert(0, str(MODULE_PATH.parent))
SPEC = importlib.util.spec_from_file_location("analyze_sentaurus_bvmethods", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class SentaurusBVMethodsAnalysisTest(unittest.TestCase):
    def test_first_upward_crossing_interpolates(self) -> None:
        value = MODULE.first_upward_crossing(
            [0.0, 2.0, 4.0], [0.0, 0.5, 1.5], 1.0)
        self.assertAlmostEqual(value, 3.0)

    def test_first_upward_crossing_honors_minimum_x(self) -> None:
        value = MODULE.first_upward_crossing(
            [0.0, 0.5, 2.0, 4.0],
            [-1.0, 1.0, -1.0, 1.0],
            0.0,
            minimum_x=1.0,
        )
        self.assertAlmostEqual(value, 3.0)

    def test_aba_poisson_uses_smaller_carrier_crossing(self) -> None:
        datasets = [
            "drain InnerVoltage", "drain OuterVoltage", "drain TotalCurrent",
            "PhiElectron", "PhiHole",
        ]
        rows = [
            [0.0, 0.0, 0.0, 0.1, 0.2],
            [4.0, 4.0, 0.0, 0.9, 0.8],
            [6.0, 6.0, 0.0, 1.1, 1.2],
        ]
        result = MODULE.analyze_method(
            "ABA_poisson", datasets, rows, 1.0e-4, 1.05)
        self.assertAlmostEqual(result["electron_bv_V"], 5.5)
        self.assertAlmostEqual(result["hole_bv_V"], 5.25)
        self.assertAlmostEqual(result["bv_V"], 5.25)


if __name__ == "__main__":
    unittest.main()
