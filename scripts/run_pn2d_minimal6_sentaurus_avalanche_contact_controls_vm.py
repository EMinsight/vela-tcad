#!/usr/bin/env python3
"""Run Minimal6 contact-forced GradQuasiFermi avalanche controls."""

from __future__ import annotations

import re
import sys
from pathlib import Path

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts import (
    run_pn2d_minimal6_sentaurus_avalanche_drive_controls_vm as base,
)


CONTACT_VARIANTS = {
    "grad_qf_use_qf_contacts": {
        "avalanche": "Avalanche(VanOverstraeten GradQuasiFermi)",
        "aval_dens_grad_qf": False,
    },
    "grad_qf_use_qf_contacts_aval_dens_grad_qf": {
        "avalanche": "Avalanche(VanOverstraeten GradQuasiFermi)",
        "aval_dens_grad_qf": True,
    },
}
BASE_MAKE_VARIANT_DECK = base.make_variant_deck


def make_contact_variant_deck(
    template: str,
    variant: str,
    biases: tuple[int, ...],
) -> str:
    result = BASE_MAKE_VARIANT_DECK(template, variant, biases)
    result, count = re.subn(
        r"Math\s*\{",
        "Math {\n"
        "  ComputeGradQuasiFermiAtContacts=UseQuasiFermi",
        result,
        count=1,
    )
    if count != 1:
        raise ValueError("Math block was not found exactly once")
    return result


def main() -> int:
    original_variants = dict(base.VARIANTS)
    original_builder = base.make_variant_deck
    try:
        base.VARIANTS.clear()
        base.VARIANTS.update(CONTACT_VARIANTS)
        base.make_variant_deck = make_contact_variant_deck
        return base.main()
    finally:
        base.VARIANTS.clear()
        base.VARIANTS.update(original_variants)
        base.make_variant_deck = original_builder


if __name__ == "__main__":
    raise SystemExit(main())
