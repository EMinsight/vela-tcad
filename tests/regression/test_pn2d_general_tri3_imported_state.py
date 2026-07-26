#!/usr/bin/env python3
"""Regression tests for general-Tri3 imported-state replay formulas."""

from __future__ import annotations

import unittest

from scripts.diagnose_pn2d_general_tri3_imported_state import (
    carrier_density,
    effective_ni,
    field_limit,
    masetti,
    p1_gradient,
)


class GeneralTri3ImportedStateTest(unittest.TestCase):
    def test_old_slotboom_effective_ni_is_frozen(self) -> None:
        ni_eff, delta_ev = effective_ni(1.0e17)
        self.assertAlmostEqual(delta_ev, 0.006363961030678928, places=16)
        self.assertAlmostEqual(ni_eff, 16556319846.864452, places=5)

    def test_masetti_defaults_match_sentaurus_2018_parameter_set(self) -> None:
        self.assertAlmostEqual(
            masetti(1.0e17, "electron"),
            727.0544030121931,
            places=11,
        )
        self.assertAlmostEqual(
            masetti(-1.0e17, "hole"),
            319.0980929489702,
            places=11,
        )
        self.assertAlmostEqual(
            field_limit(727.0544030121931, 1.0e5, "electron"),
            96.6502588118563,
            places=12,
        )

    def test_carrier_statistics_use_imported_qfp_signs(self) -> None:
        ni_eff, _ = effective_ni(1.0e17)
        self.assertGreater(
            carrier_density(ni_eff, 0.2, 0.0, "electron"),
            ni_eff,
        )
        self.assertGreater(
            carrier_density(ni_eff, 0.0, 0.2, "hole"),
            ni_eff,
        )

    def test_p1_gradient_preserves_original_cell_order(self) -> None:
        points = [(0.0, 0.0), (0.2, 0.8), (1.0, 0.0)]
        values = [0.0, 1.8, 2.0]
        gx, gy = p1_gradient(points, values)
        self.assertAlmostEqual(gx, 2.0)
        self.assertAlmostEqual(gy, 1.75)


if __name__ == "__main__":
    unittest.main()
