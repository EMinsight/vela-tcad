import math
import unittest

from scripts.pn2d_minimal6_diagnostics.mobility_diagnosis import (
    CONSTANT_MOBILITY_M2_PER_VS,
    cell_inverted_gradient,
    field_limited_mobility,
    masetti_low_field_mobility,
    unique_edges,
    vela_masetti_edge_mobility,
)


class Minimal6MobilityDiagnosisTest(unittest.TestCase):
    def test_masetti_zero_doping_and_high_field_limits_match_production_defaults(self):
        self.assertEqual(masetti_low_field_mobility("electron", 0.0), 0.14170)
        self.assertEqual(masetti_low_field_mobility("hole", 0.0), 0.04705)

        electron_low = masetti_low_field_mobility("electron", 1.0e23)
        hole_low = masetti_low_field_mobility("hole", -1.0e23)
        self.assertAlmostEqual(electron_low, 0.0727054403012193)
        self.assertAlmostEqual(hole_low, 0.03190980929489702)
        self.assertLess(field_limited_mobility("electron", electron_low, 2.0e6), electron_low)
        self.assertLess(field_limited_mobility("hole", hole_low, 2.0e6), hole_low)

    def test_edge_mobility_uses_average_net_doping_and_qf_difference(self):
        value = vela_masetti_edge_mobility(
            "electron",
            net_doping0_m3=-1.0e23,
            net_doping1_m3=1.0e23,
            qf0_V=-1.0,
            qf1_V=0.0,
            length_m=2.0e-6,
        )
        expected = field_limited_mobility("electron", 0.14170, 5.0e5)
        self.assertAlmostEqual(value, expected)
        self.assertAlmostEqual(CONSTANT_MOBILITY_M2_PER_VS["electron"], 0.14170)
        self.assertAlmostEqual(CONSTANT_MOBILITY_M2_PER_VS["hole"], 0.04705)

    def test_unique_edges_are_canonical_and_order_independent(self):
        triangles = {"1": ("2", "1", "0"), "0": ("0", "1", "3")}
        self.assertEqual(
            unique_edges(triangles),
            (("0", "1"), ("0", "2"), ("0", "3"), ("1", "2"), ("1", "3")),
        )

    def test_scalar_mobility_changes_magnitude_but_cannot_rotate_cell_inversion(self):
        current = (-24.0, 30.0)
        first = cell_inverted_gradient("electron", 2.0, 3.0, current, q=1.0)
        second = cell_inverted_gradient("electron", 2.0, 6.0, current, q=1.0)
        self.assertEqual(first, (4.0, -5.0))
        self.assertEqual(second, (2.0, -2.5))
        cross = first[0] * second[1] - first[1] * second[0]
        self.assertAlmostEqual(cross, 0.0)
        self.assertGreater(first[0] * second[0] + first[1] * second[1], 0.0)
        self.assertEqual(
            cell_inverted_gradient("hole", 2.0, 3.0, current, q=1.0),
            (4.0, -5.0),
        )

    def test_invalid_denominators_fail_closed(self):
        for density, mobility in ((0.0, 1.0), (1.0, 0.0), (math.nan, 1.0)):
            with self.subTest(density=density, mobility=mobility):
                with self.assertRaises(ValueError):
                    cell_inverted_gradient("hole", density, mobility, (1.0, 0.0))
        self.assertEqual(
            cell_inverted_gradient("hole", 2.0, 3.0, (-24.0, 30.0), q=1.0),
            (4.0, -5.0),
        )


if __name__ == "__main__":
    unittest.main()
