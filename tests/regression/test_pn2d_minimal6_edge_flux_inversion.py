import math
import unittest

from scripts.pn2d_minimal6_diagnostics.edge_flux_inversion import (
    canonical_edges,
    continuity_flux_from_current,
    edge_current_supports,
    required_positive_mobility,
    staged_sg_flux,
)


class Minimal6EdgeFluxInversionTest(unittest.TestCase):
    def test_canonical_edges_are_unique_and_sorted(self):
        triangles = ((0, 4, 5), (0, 5, 1), (1, 5, 2), (5, 3, 2))
        self.assertEqual(
            canonical_edges(triangles),
            (
                (0, 1),
                (0, 4),
                (0, 5),
                (1, 2),
                (1, 5),
                (2, 3),
                (2, 5),
                (3, 5),
                (4, 5),
            ),
        )

    def test_support_reconstructions_are_explicit_and_p1_line_mean_is_not_new(self):
        coordinates = {
            0: (0.0, 0.0),
            1: (1.0, 0.0),
            2: (0.0, 1.0),
            3: (1.0, 1.0),
        }
        triangles = ((0, 1, 2), (1, 3, 2))
        vectors = {
            0: (0.0, 0.0),
            1: (2.0, 0.0),
            2: (8.0, 0.0),
            3: (20.0, 0.0),
        }
        supports = edge_current_supports(coordinates, triangles, vectors)

        shared = supports[(1, 2)]
        expected_endpoint = -((2.0 + 8.0) * 0.5) / math.sqrt(2.0)
        cell0 = (0.0 + 2.0 + 8.0) / 3.0
        cell1 = (2.0 + 20.0 + 8.0) / 3.0
        expected_cell = -((cell0 + cell1) * 0.5) / math.sqrt(2.0)
        self.assertAlmostEqual(shared["endpoint_mean_tangent"], expected_endpoint)
        self.assertAlmostEqual(shared["p1_line_mean_tangent"], expected_endpoint)
        self.assertAlmostEqual(shared["adjacent_cell_mean_tangent"], expected_cell)
        self.assertNotEqual(
            shared["endpoint_mean_tangent"],
            shared["adjacent_cell_mean_tangent"],
        )

    def test_current_conversion_uses_continuity_sign_conventions(self):
        q = 1.602176634e-19
        self.assertEqual(continuity_flux_from_current("electron", q), -1.0)
        self.assertEqual(continuity_flux_from_current("hole", q), 1.0)

    def test_required_mobility_is_typed_and_never_negative(self):
        self.assertEqual(
            required_positive_mobility(reference_flux=12.0, unit_mobility_flux=3.0),
            {"classification": "available", "mobility_m2_per_Vs": 4.0},
        )
        self.assertEqual(
            required_positive_mobility(reference_flux=-12.0, unit_mobility_flux=3.0),
            {
                "classification": "sign_incompatible",
                "mobility_m2_per_Vs": None,
            },
        )
        self.assertEqual(
            required_positive_mobility(reference_flux=1.0, unit_mobility_flux=0.0),
            {"classification": "zero_operator", "mobility_m2_per_Vs": None},
        )

    def test_staged_dependencies_make_irrelevant_replacements_exact_no_ops(self):
        state = {
            "psi0_V": -0.2,
            "psi1_V": 0.3,
            "qf0_V": -0.4,
            "qf1_V": -0.1,
            "density0_m3": 2.0e19,
            "density1_m3": 4.0e19,
            "ni0_m3": 1.0e16,
            "ni1_m3": 1.0e16,
            "mobility_m2_per_Vs": 0.08,
            "length_m": 1.0e-6,
            "thermal_voltage_V": 0.025851999786435535,
        }
        qf_baseline = staged_sg_flux("qf_sg", "electron", state)
        density_changed = dict(state, density0_m3=1.0, density1_m3=3.0)
        self.assertEqual(
            staged_sg_flux("qf_sg", "electron", density_changed),
            qf_baseline,
        )

        density_baseline = staged_sg_flux("density_sg", "electron", state)
        qf_changed = dict(state, qf0_V=20.0, qf1_V=-30.0)
        self.assertEqual(
            staged_sg_flux("density_sg", "electron", qf_changed),
            density_baseline,
        )


if __name__ == "__main__":
    unittest.main()
