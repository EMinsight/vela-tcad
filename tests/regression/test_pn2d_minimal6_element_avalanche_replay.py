#!/usr/bin/env python3
"""Unit tests for the Minimal6 element avalanche replay geometry."""

from __future__ import annotations

import math
import unittest

from scripts.diagnose_pn2d_minimal6_element_avalanche_replay import (
    charon_whitney_vector,
    gss_laux_vector,
    least_squares,
    solve_pair,
)


class ElementAvalancheReplayTest(unittest.TestCase):
    def test_pair_and_least_squares_recover_constant_vector(self) -> None:
        vector = (3.0, -2.0)
        tangents = (
            (1.0, 0.0),
            (0.0, 1.0),
            (math.sqrt(0.5), math.sqrt(0.5)),
        )
        values = [
            tangent[0] * vector[0] + tangent[1] * vector[1]
            for tangent in tangents
        ]
        self.assertEqual(
            solve_pair(tangents[0], values[0], tangents[1], values[1]),
            vector,
        )
        recovered = least_squares(
            (tangent, value, 1.0)
            for tangent, value in zip(tangents, values, strict=True)
        )
        self.assertAlmostEqual(recovered[0], vector[0], places=14)
        self.assertAlmostEqual(recovered[1], vector[1], places=14)

    def test_gss_laux_recovers_constant_vector_with_zero_diagonal(self) -> None:
        vector = (4.0, -5.0)
        raw = [
            ((0.0, -1.0), 1.0, 0.5),
            ((math.sqrt(0.8), -math.sqrt(0.2)), 0.0, math.sqrt(1.25)),
            ((1.0, 0.0), 0.25, 1.0),
        ]
        edges = []
        for tangent, kappa, length in raw:
            edges.append(
                {
                    "tangent_x": tangent[0],
                    "tangent_y": tangent[1],
                    "kappa": kappa,
                    "length_um": length,
                    "sg_jn_A_cm2": (
                        tangent[0] * vector[0] + tangent[1] * vector[1]
                    ),
                }
            )
        recovered = gss_laux_vector(edges, "sg_jn_A_cm2")
        self.assertAlmostEqual(recovered[0], vector[0], places=14)
        self.assertAlmostEqual(recovered[1], vector[1], places=14)

    def test_charon_whitney_recovers_constant_vector(self) -> None:
        vector = (2.5, -1.25)
        vertices = {
            0: {"x_um": 0.0, "y_um": 0.0},
            1: {"x_um": 1.0, "y_um": 0.0},
            2: {"x_um": 0.0, "y_um": 0.5},
        }
        edges = []
        for start, end in ((0, 1), (1, 2), (2, 0)):
            dx = vertices[end]["x_um"] - vertices[start]["x_um"]
            dy = vertices[end]["y_um"] - vertices[start]["y_um"]
            length = math.hypot(dx, dy)
            tangent = (dx / length, dy / length)
            edges.append(
                {
                    "start": start,
                    "end": end,
                    "length_um": length,
                    "sg_jn_A_cm2": (
                        tangent[0] * vector[0] + tangent[1] * vector[1]
                    ),
                }
            )
        recovered = charon_whitney_vector(
            edges, [0, 1, 2], vertices, "sg_jn_A_cm2"
        )
        self.assertAlmostEqual(recovered[0], vector[0], places=13)
        self.assertAlmostEqual(recovered[1], vector[1], places=13)


if __name__ == "__main__":
    unittest.main()
