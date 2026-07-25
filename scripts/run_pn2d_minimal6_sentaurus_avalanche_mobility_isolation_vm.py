#!/usr/bin/env python3
"""Run low-field-mobility controls isolating the avalanche driving force."""

from __future__ import annotations

import re
import sys
from pathlib import Path

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts import (
    run_pn2d_minimal6_sentaurus_avalanche_drive_controls_vm as base,
)


ISOLATION_VARIANTS = {
    "lowfield_mobility_avalanche_electric_field": {
        "avalanche": "Avalanche(VanOverstraeten ElectricField)",
        "aval_dens_grad_qf": False,
    },
    "lowfield_mobility_avalanche_grad_qf": {
        "avalanche": "Avalanche(VanOverstraeten GradQuasiFermi)",
        "aval_dens_grad_qf": False,
    },
}
BASE_MAKE_VARIANT_DECK = base.make_variant_deck


def make_isolation_variant_deck(
    template: str,
    variant: str,
    biases: tuple[int, ...],
) -> str:
    result = BASE_MAKE_VARIANT_DECK(template, variant, biases)
    result, mobility_count = re.subn(
        r"Mobility\s*\(\s*DopingDependence\s*"
        r"HighFieldSaturation\s*\)",
        "Mobility(\n    DopingDependence\n  )",
        result,
        count=1,
        flags=re.S,
    )
    if mobility_count != 1:
        raise ValueError("HighFieldSaturation block was not removed once")
    result, math_count = re.subn(
        r"Math\s*\{",
        "Math {\n"
        "  ComputeGradQuasiFermiAtContacts=UseQuasiFermi",
        result,
        count=1,
    )
    if math_count != 1:
        raise ValueError("Math block was not found exactly once")
    return result


def main() -> int:
    original_variants = dict(base.VARIANTS)
    original_builder = base.make_variant_deck
    try:
        base.VARIANTS.clear()
        base.VARIANTS.update(ISOLATION_VARIANTS)
        base.make_variant_deck = make_isolation_variant_deck
        return base.main()
    finally:
        base.VARIANTS.clear()
        base.VARIANTS.update(original_variants)
        base.make_variant_deck = original_builder


if __name__ == "__main__":
    raise SystemExit(main())
