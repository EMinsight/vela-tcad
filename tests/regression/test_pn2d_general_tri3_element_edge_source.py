#!/usr/bin/env python3
"""Unit tests for general-Tri3 element-edge avalanche source helpers."""

from __future__ import annotations

import unittest

from scripts.diagnose_pn2d_general_tri3_element_edge_source import (
    CURRENT_METHODS,
    integrated_qg_A_um,
    matching_driver,
)


class GeneralTri3ElementEdgeSourceTest(unittest.TestCase):
    def test_driver_contract_uses_qfp_away_from_contacts(self) -> None:
        self.assertEqual(
            matching_driver("interior"),
            "quasi_fermi_gradient",
        )

    def test_driver_contract_keeps_contact_fallback_explicit(self) -> None:
        self.assertEqual(matching_driver("contact"), "electric_field")

    def test_integrated_source_unit_conversion(self) -> None:
        self.assertAlmostEqual(
            integrated_qg_A_um(2.0e4, 3.0e2, 5.0),
            3.0e-5,
            places=18,
        )

    def test_all_current_vector_candidates_are_frozen(self) -> None:
        self.assertEqual(
            CURRENT_METHODS,
            (
                "gss_laux_truncated_support",
                "charon_whitney_hcurl_cell_average",
                "genius_least_squares_tangent",
                "box_active_edge_exact",
            ),
        )


if __name__ == "__main__":
    unittest.main()
