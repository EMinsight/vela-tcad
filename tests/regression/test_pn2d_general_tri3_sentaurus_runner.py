#!/usr/bin/env python3
"""Tests for general-Tri3 Sentaurus avalanche control generation."""

from __future__ import annotations

import unittest

from scripts.diagnose_pn2d_general_tri3_element_edge_avalanche import (
    geometry_rows,
)
from scripts.run_pn2d_general_tri3_sentaurus_avalanche_controls_vm import (
    VARIANTS,
    make_variant_deck,
    validate_case_name,
)


BASE_DECK = """File {
  Grid = "pn2d_minimal6.tdr"
  Plot = "runtime_element_avalanche_probe_default.tdr"
  Current = "runtime_element_avalanche_probe_default.plt"
  Output = "runtime_element_avalanche_probe_default"
}
Physics {
  Mobility(
    DopingDependence
    HighFieldSaturation
  )
  Recombination(SRH Avalanche(VanOverstraeten))
}
Math {
  Extrapolate
}
Solve {
  Coupled { Poisson Electron Hole }
}
"""


class GeneralTri3SentaurusRunnerTest(unittest.TestCase):
    def test_case_name_rejects_remote_path_syntax(self) -> None:
        self.assertEqual(validate_case_name("skewed_tri3"), "skewed_tri3")
        for value in ("../bad", "has space", "case;touch_bad", ""):
            with self.subTest(value=value):
                with self.assertRaises(ValueError):
                    validate_case_name(value)

    def test_variant_matrix_is_exact(self) -> None:
        self.assertEqual(
            tuple(VARIANTS),
            (
                "implicit_default",
                "explicit_grad_qf",
                "explicit_electric_field",
                "grad_qf_use_qf_contacts",
                "electric_field_use_qf_contacts",
                "grad_qf_aval_dens_grad_qf",
                "lowfield_mobility_avalanche_electric_field",
                "lowfield_mobility_avalanche_grad_qf",
            ),
        )

    def test_lowfield_pair_differs_only_by_avalanche_drive(self) -> None:
        electric = make_variant_deck(
            BASE_DECK,
            "lowfield_mobility_avalanche_electric_field",
            (-1, -10, -20),
        )
        grad_qf = make_variant_deck(
            BASE_DECK,
            "lowfield_mobility_avalanche_grad_qf",
            (-1, -10, -20),
        )
        for deck in (electric, grad_qf):
            self.assertNotIn("HighFieldSaturation", deck)
            self.assertEqual(
                deck.count(
                    "ComputeGradQuasiFermiAtContacts=UseQuasiFermi"
                ),
                1,
            )
            self.assertIn('Grid = "pn2d_msh.tdr"', deck)
        normalized = electric.replace("ElectricField", "GradQuasiFermi")
        normalized = normalized.replace(
            "lowfield_mobility_avalanche_electric_field",
            "lowfield_mobility_avalanche_grad_qf",
        )
        self.assertEqual(normalized, grad_qf)

    def test_geometry_uses_original_element_vertex_permutation(self) -> None:
        groups = {
            "vertices": [
                {"bias_V": -1.0, "vertex": 0, "x_um": 0.0, "y_um": 0.0},
                {"bias_V": -1.0, "vertex": 1, "x_um": 1.0, "y_um": 0.0},
                {"bias_V": -1.0, "vertex": 2, "x_um": 0.2, "y_um": 0.8},
            ],
            "measures": [
                {
                    "bias_V": -1.0,
                    "element": 0,
                    "local_vertex": 0,
                    "vertex": 0,
                },
                {
                    "bias_V": -1.0,
                    "element": 0,
                    "local_vertex": 1,
                    "vertex": 2,
                },
                {
                    "bias_V": -1.0,
                    "element": 0,
                    "local_vertex": 2,
                    "vertex": 1,
                },
            ],
            "edges": [
                {
                    "bias_V": -1.0,
                    "element": 0,
                    "kappa": 0.1,
                    "start_x_um": 0.0,
                    "start_y_um": 0.0,
                    "end_x_um": 1.0,
                    "end_y_um": 0.0,
                },
                {
                    "bias_V": -1.0,
                    "element": 0,
                    "kappa": 0.2,
                    "start_x_um": 1.0,
                    "start_y_um": 0.0,
                    "end_x_um": 0.2,
                    "end_y_um": 0.8,
                },
                {
                    "bias_V": -1.0,
                    "element": 0,
                    "kappa": 0.3,
                    "start_x_um": 0.2,
                    "start_y_um": 0.8,
                    "end_x_um": 0.0,
                    "end_y_um": 0.0,
                },
            ],
        }
        row = geometry_rows(groups, -1.0)[0]
        self.assertEqual(row["cell_vertex_permutation"], "0;2;1")
        self.assertEqual(row["orientation"], "cw")
        self.assertAlmostEqual(row["signed_area_um2"], -0.4)


if __name__ == "__main__":
    unittest.main()
