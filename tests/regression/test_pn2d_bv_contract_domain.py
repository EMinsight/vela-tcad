#!/usr/bin/env python3
"""Tests for the predeclared PN2D BV contract-domain analyzer."""

from __future__ import annotations

import sys
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
if str(REPO / "scripts") not in sys.path:
    sys.path.insert(0, str(REPO / "scripts"))

from analyze_pn2d_avalanche_on_bv_parity import CurvePoint, KNEE_BIASES_V  # noqa: E402
from analyze_pn2d_bv_contract_domain import analyze  # noqa: E402


class Pn2dBvContractDomainTest(unittest.TestCase):
    def test_zero_current_outside_predeclared_domain_is_not_floored(self) -> None:
        biases = (-15.0, -16.0, -17.0, *KNEE_BIASES_V)
        unique_biases = tuple(dict.fromkeys(biases))
        curves = {}
        for name in ("vela_on", "vela_off", "sentaurus_on", "sentaurus_off"):
            points = [CurvePoint(0.0, 0.0)]
            for bias in unique_biases:
                magnitude = 10.0 ** (abs(bias) - 25.0)
                if name.endswith("_off"):
                    magnitude *= 0.5
                points.append(CurvePoint(bias, magnitude))
            curves[name] = points
        result = analyze(curves, unique_biases)
        self.assertEqual(result["outcome"], "contract_domain_metrics_complete")
        self.assertEqual(
            result["effective_metrics"]["maximum_absolute_log_error_dex"], 0.0
        )
        self.assertTrue(
            result["zero_current_outside_contract_domain_is_not_evaluated"]
        )


if __name__ == "__main__":
    unittest.main()
