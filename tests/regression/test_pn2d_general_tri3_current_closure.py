#!/usr/bin/env python3
"""Unit tests for general-Tri3 current-vector closure helpers."""

from __future__ import annotations

import math
import unittest

from scripts.close_pn2d_general_tri3_element_edge_current import (
    candidate_vectors,
    charon_whitney_vector,
    least_squares,
    signal_status,
)


def edge_projections(
    points: list[tuple[float, float]],
    vector: tuple[float, float],
) -> list[float]:
    result = []
    for index in range(3):
        start = points[index]
        end = points[(index + 1) % 3]
        delta = (end[0] - start[0], end[1] - start[1])
        length = math.hypot(*delta)
        result.append(
            vector[0] * delta[0] / length
            + vector[1] * delta[1] / length
        )
    return result


class GeneralTri3CurrentClosureTest(unittest.TestCase):
    def test_least_squares_recovers_constant_vector(self) -> None:
        points = [(0.0, 0.0), (1.7, 0.2), (0.3, 1.1)]
        expected = (3.0, -2.0)
        values = edge_projections(points, expected)
        rows = []
        for index, value in enumerate(values):
            start = points[index]
            end = points[(index + 1) % 3]
            delta = (end[0] - start[0], end[1] - start[1])
            length = math.hypot(*delta)
            rows.append(
                ((delta[0] / length, delta[1] / length), value, 1.0)
            )
        actual = least_squares(rows)
        self.assertAlmostEqual(actual[0], expected[0], places=12)
        self.assertAlmostEqual(actual[1], expected[1], places=12)

    def test_whitney_recovers_constant_vector(self) -> None:
        points = [(0.0, 0.0), (1.7, 0.2), (0.3, 1.1)]
        expected = (3.0, -2.0)
        actual = charon_whitney_vector(
            points,
            edge_projections(points, expected),
        )
        self.assertAlmostEqual(actual[0], expected[0], places=12)
        self.assertAlmostEqual(actual[1], expected[1], places=12)

    def test_all_candidates_recover_constant_vector(self) -> None:
        points = [(0.0, 0.0), (1.7, 0.2), (0.3, 1.1)]
        expected = (3.0, -2.0)
        values = edge_projections(points, expected)
        for method, actual in candidate_vectors(
            points,
            values,
            [0.4, 0.3, 0.2],
        ).items():
            with self.subTest(method=method):
                self.assertAlmostEqual(actual[0], expected[0], places=12)
                self.assertAlmostEqual(actual[1], expected[1], places=12)

    def test_right_triangle_active_edge_exact_ignores_zero_support(self) -> None:
        points = [(0.0, 0.0), (2.0, 0.0), (0.0, 1.0)]
        expected = (2.5, -0.75)
        actual = candidate_vectors(
            points,
            edge_projections(points, expected),
            [0.5, 0.0, 0.5],
        )["box_active_edge_exact"]
        self.assertAlmostEqual(actual[0], expected[0], places=12)
        self.assertAlmostEqual(actual[1], expected[1], places=12)

    def test_signal_status_is_relative_to_each_state(self) -> None:
        self.assertEqual(signal_status(0.0, 10.0), "zero")
        self.assertEqual(
            signal_status(1.0e-13, 1.0),
            "below_state_relative_floor",
        )
        self.assertEqual(signal_status(1.0e-11, 1.0), "valid")


if __name__ == "__main__":
    unittest.main()
