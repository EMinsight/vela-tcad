from __future__ import annotations

import unittest

from scripts.audit_pn2d_task7_low_current_nonmonotonicity import (
    reverse_intervals,
)


class Task7LowCurrentNonmonotonicityTests(unittest.TestCase):
    def test_reverse_intervals_follow_increasing_reverse_bias(self) -> None:
        rows = [
            {"bias_V": -5.0, "current": -3.2e-17},
            {"bias_V": -3.0, "current": -3.1e-17},
            {"bias_V": -4.0, "current": -3.0e-17},
        ]

        self.assertEqual(
            reverse_intervals(rows, "current"),
            [{
                "left_bias_V": -3.0,
                "right_bias_V": -4.0,
                "left_abs_current_A_per_um": 3.1e-17,
                "right_abs_current_A_per_um": 3.0e-17,
            }],
        )


if __name__ == "__main__":
    unittest.main()
