#!/usr/bin/env python3
"""Unit tests for the general Tri3 SG-current diagnostic helpers."""

from __future__ import annotations

import math
import unittest

from scripts.diagnose_pn2d_general_tri3_element_edge_current import (
    bernoulli,
    geometric_partial_volumes,
    gss_laux_vector,
    variable_ni_flux,
)
from scripts.diagnose_pn2d_general_tri3_imported_state import VT_300K_V


class GeneralTri3ElementEdgeCurrentTest(unittest.TestCase):
    def test_bernoulli_identity(self) -> None:
        for value in (
            -60.0,
            -1.0,
            -1.0e-10,
            0.0,
            1.0e-10,
            1.0,
            60.0,
        ):
            self.assertAlmostEqual(
                bernoulli(-value) - bernoulli(value),
                value,
                places=12,
            )

    def test_variable_ni_flux_is_zero_for_flat_qfp(self) -> None:
        for carrier in ("electron", "hole"):
            self.assertEqual(
                variable_ni_flux(
                    1.0e10,
                    2.0e10,
                    -0.3,
                    0.4,
                    0.2,
                    0.2,
                    500.0,
                    1.0e-4,
                    carrier,
                ),
                0.0,
            )

    def test_hole_tcl_current_is_negative_continuity_flux(self) -> None:
        ni0 = 1.1e10
        ni1 = 1.7e10
        psi0 = -0.1
        psi1 = 0.2
        phip0 = 0.03
        phip1 = -0.04
        mobility = 350.0
        length_cm = 2.0e-4
        flux = variable_ni_flux(
            ni0,
            ni1,
            psi0,
            psi1,
            phip0,
            phip1,
            mobility,
            length_cm,
            "hole",
        )
        p0 = ni0 * math.exp((phip0 - psi0) / VT_300K_V)
        p1 = ni1 * math.exp((phip1 - psi1) / VT_300K_V)
        xp = math.log(p0 / p1) - (phip0 - phip1) / VT_300K_V
        tcl_current_without_q = (
            mobility
            * VT_300K_V
            / length_cm
            * (p1 * bernoulli(-xp) - p0 * bernoulli(xp))
        )
        self.assertTrue(
            math.isclose(
                tcl_current_without_q,
                -flux,
                rel_tol=3.0e-15,
                abs_tol=0.0,
            )
        )

    def test_right_triangle_hypotenuse_partial_is_zero(self) -> None:
        partials = geometric_partial_volumes(
            [(0.0, 0.0), (2.0, 0.0), (0.0, 1.0)]
        )
        self.assertAlmostEqual(partials[0], 0.5)
        self.assertAlmostEqual(partials[1], 0.0)
        self.assertAlmostEqual(partials[2], 0.5)

    def test_obtuse_opposite_support_is_truncated(self) -> None:
        partials = geometric_partial_volumes(
            [(0.0, 0.0), (2.0, 0.0), (0.2, 0.1)]
        )
        self.assertEqual(partials[0], 0.0)
        self.assertGreater(partials[1], 0.0)
        self.assertGreater(partials[2], 0.0)

    def test_gss_laux_recovers_constant_vector(self) -> None:
        points = [(0.0, 0.0), (1.7, 0.2), (0.3, 1.1)]
        vector = (3.0, -2.0)
        edge_values = []
        for index in range(3):
            start = points[index]
            end = points[(index + 1) % 3]
            delta = (end[0] - start[0], end[1] - start[1])
            length = math.hypot(*delta)
            tangent = (delta[0] / length, delta[1] / length)
            edge_values.append(
                vector[0] * tangent[0] + vector[1] * tangent[1]
            )
        reconstructed = gss_laux_vector(
            points,
            edge_values,
            geometric_partial_volumes(points),
        )
        self.assertAlmostEqual(reconstructed[0], vector[0], places=12)
        self.assertAlmostEqual(reconstructed[1], vector[1], places=12)


if __name__ == "__main__":
    unittest.main()
