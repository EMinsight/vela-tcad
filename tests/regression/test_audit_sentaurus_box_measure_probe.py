from __future__ import annotations

import math
import unittest

from scripts.audit_sentaurus_box_measure_probe import (
    geometry,
    raw_circumcentric_shares,
    vela_mixed_shares,
)


class AuditSentaurusBoxMeasureProbeTests(unittest.TestCase):
    def test_raw_circumcentric_has_negative_share_for_obtuse_triangle(self) -> None:
        points = [(0.25, 0.25), (0.0, 0.0), (0.05, 0.0)]
        area, angles = geometry(points)
        raw = raw_circumcentric_shares(points)
        mixed = vela_mixed_shares(points, area, angles)

        self.assertGreater(max(angles), 90.0)
        self.assertTrue(any(value < 0.0 for value in raw))
        self.assertTrue(all(value > 0.0 for value in mixed))
        self.assertTrue(math.isclose(sum(mixed), area))


if __name__ == "__main__":
    unittest.main()
