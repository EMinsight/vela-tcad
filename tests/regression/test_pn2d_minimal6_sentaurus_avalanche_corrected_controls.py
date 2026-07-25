#!/usr/bin/env python3
"""Schema tests for the corrected Sentaurus avalanche control matrix."""

from __future__ import annotations

import unittest

from scripts.compare_pn2d_minimal6_sentaurus_avalanche_corrected_controls import (
    PAIR_SPECS,
)


class SentaurusAvalancheCorrectedControlsTest(unittest.TestCase):
    def test_pair_matrix_separates_drive_contact_fallback_and_current(self) -> None:
        self.assertEqual(
            PAIR_SPECS,
            (
                (
                    "implicit_default",
                    "base",
                    "explicit_grad_qf",
                    "base",
                ),
                (
                    "explicit_electric_field",
                    "base",
                    "explicit_grad_qf",
                    "base",
                ),
                (
                    "grad_qf_aval_dens_grad_qf",
                    "base",
                    "explicit_grad_qf",
                    "base",
                ),
                (
                    "grad_qf_use_qf_contacts",
                    "contact",
                    "explicit_grad_qf",
                    "base",
                ),
                (
                    "grad_qf_use_qf_contacts_aval_dens_grad_qf",
                    "contact",
                    "grad_qf_use_qf_contacts",
                    "contact",
                ),
            ),
        )


if __name__ == "__main__":
    unittest.main()
