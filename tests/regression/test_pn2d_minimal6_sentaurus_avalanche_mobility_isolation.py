#!/usr/bin/env python3
"""Tests for low-field mobility avalanche-drive isolation decks."""

from __future__ import annotations

import unittest
from unittest import mock

from scripts import (
    run_pn2d_minimal6_sentaurus_avalanche_mobility_isolation_vm as controls,
)


BASE_DECK = """File {
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


class SentaurusAvalancheMobilityIsolationTest(unittest.TestCase):
    def setUp(self) -> None:
        self.original = dict(controls.base.VARIANTS)
        controls.base.VARIANTS.clear()
        controls.base.VARIANTS.update(controls.ISOLATION_VARIANTS)

    def tearDown(self) -> None:
        controls.base.VARIANTS.clear()
        controls.base.VARIANTS.update(self.original)

    def test_isolation_removes_hfs_and_forces_true_qf_at_contacts(self) -> None:
        for variant in controls.ISOLATION_VARIANTS:
            deck = controls.make_isolation_variant_deck(
                BASE_DECK,
                variant,
                (-1, -10, -20),
            )
            self.assertNotIn("HighFieldSaturation", deck)
            self.assertEqual(deck.count("DopingDependence"), 1)
            self.assertEqual(
                deck.count(
                    "ComputeGradQuasiFermiAtContacts=UseQuasiFermi"
                ),
                1,
            )

    def test_only_avalanche_drive_differs_between_isolation_decks(self) -> None:
        electric = controls.make_isolation_variant_deck(
            BASE_DECK,
            "lowfield_mobility_avalanche_electric_field",
            (-20,),
        )
        grad_qf = controls.make_isolation_variant_deck(
            BASE_DECK,
            "lowfield_mobility_avalanche_grad_qf",
            (-20,),
        )
        self.assertIn(
            "Avalanche(VanOverstraeten ElectricField)",
            electric,
        )
        self.assertIn(
            "Avalanche(VanOverstraeten GradQuasiFermi)",
            grad_qf,
        )
        normalized = electric.replace("ElectricField", "GradQuasiFermi")
        normalized = normalized.replace(
            "lowfield_mobility_avalanche_electric_field",
            "lowfield_mobility_avalanche_grad_qf",
        )
        self.assertEqual(normalized, grad_qf)

    def test_isolation_main_restores_base_module_on_failure(self) -> None:
        original_variants = dict(controls.base.VARIANTS)
        original_builder = controls.base.make_variant_deck
        with mock.patch.object(
            controls.base,
            "main",
            side_effect=RuntimeError("synthetic failure"),
        ):
            with self.assertRaisesRegex(RuntimeError, "synthetic failure"):
                controls.main()
        self.assertEqual(controls.base.VARIANTS, original_variants)
        self.assertIs(controls.base.make_variant_deck, original_builder)


if __name__ == "__main__":
    unittest.main()
