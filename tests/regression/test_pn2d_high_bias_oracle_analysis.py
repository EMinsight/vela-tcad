#!/usr/bin/env python3

from __future__ import annotations

import math
import tempfile
import unittest
from pathlib import Path

from scripts.analyze_pn2d_high_bias_oracle import adjacent, records


class HighBiasOracleAnalysisTest(unittest.TestCase):
    def test_runtime_record_parser_keeps_decimal_bias_and_components(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "run.out"
            path.write_text(
                "ignored\n"
                "AVAL_PROBE_PROCESS bias_V=-19.95 vertex=2 "
                "total_current_x_A_cm2=1.25e-3 "
                "total_current_y_A_cm2=-2.5e-4\n",
                encoding="ascii",
            )
            parsed = records(path)
        self.assertEqual(len(parsed), 1)
        self.assertEqual(parsed[0]["kind"], "process")
        self.assertEqual(parsed[0]["bias_V"], -19.95)
        self.assertEqual(parsed[0]["vertex"], 2.0)
        self.assertEqual(parsed[0]["total_current_y_A_cm2"], -2.5e-4)

    def test_adjacent_growth_uses_actual_voltage_spacing(self) -> None:
        rows = [
            {"variant": "implicit_default", "bias_V": -19.9, "current": 2.0},
            {"variant": "implicit_default", "bias_V": -19.95, "current": 3.0},
        ]
        result = adjacent(rows)
        self.assertEqual(len(result), 1)
        self.assertAlmostEqual(result[0]["ratio"], 1.5)
        self.assertAlmostEqual(
            result[0]["log_slope_per_V"],
            math.log(1.5) / 0.05,
        )


if __name__ == "__main__":
    unittest.main()
