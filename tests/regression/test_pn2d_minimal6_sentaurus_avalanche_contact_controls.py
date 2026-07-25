#!/usr/bin/env python3
"""Tests for contact-forced Sentaurus avalanche controls."""

from __future__ import annotations

import unittest
from unittest import mock

from scripts import (
    run_pn2d_minimal6_sentaurus_avalanche_contact_controls_vm as controls,
)


BASE_DECK = """File {
  Plot = "runtime_element_avalanche_probe_default.tdr"
  Current = "runtime_element_avalanche_probe_default.plt"
  Output = "runtime_element_avalanche_probe_default"
}
Physics {
  Recombination(SRH Avalanche(VanOverstraeten))
}
Math {
  Extrapolate
}
Solve {
  Coupled { Poisson Electron Hole }
}
"""


class SentaurusAvalancheContactControlsTest(unittest.TestCase):
    def setUp(self) -> None:
        self.original = dict(controls.base.VARIANTS)
        controls.base.VARIANTS.clear()
        controls.base.VARIANTS.update(controls.CONTACT_VARIANTS)

    def tearDown(self) -> None:
        controls.base.VARIANTS.clear()
        controls.base.VARIANTS.update(self.original)

    def test_contact_variants_force_qf_gradient_at_all_elements(self) -> None:
        for variant in controls.CONTACT_VARIANTS:
            deck = controls.make_contact_variant_deck(
                BASE_DECK,
                variant,
                (-1, -10, -20),
            )
            self.assertEqual(
                deck.count(
                    "ComputeGradQuasiFermiAtContacts=UseQuasiFermi"
                ),
                1,
            )
            self.assertIn(
                "Avalanche(VanOverstraeten GradQuasiFermi)",
                deck,
            )

    def test_aval_dens_variant_is_orthogonal(self) -> None:
        baseline = controls.make_contact_variant_deck(
            BASE_DECK,
            "grad_qf_use_qf_contacts",
            (-1,),
        )
        candidate = controls.make_contact_variant_deck(
            BASE_DECK,
            "grad_qf_use_qf_contacts_aval_dens_grad_qf",
            (-1,),
        )
        self.assertNotIn("AvalDensGradQF", baseline)
        self.assertEqual(candidate.count("AvalDensGradQF"), 1)

    def test_contact_main_restores_base_module_on_failure(self) -> None:
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
