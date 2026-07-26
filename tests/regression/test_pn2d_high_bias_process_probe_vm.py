#!/usr/bin/env python3

from __future__ import annotations

import unittest

from scripts.run_pn2d_high_bias_process_probe_vm import (
    BIAS_V,
    PROCESS_FIELDS,
    process_deck,
    process_tcl,
)


DECK = """\
File {
  Grid = "pn2d_minimal6.tdr"
  Plot = "runtime_element_avalanche_probe_default.tdr"
  Current = "runtime_element_avalanche_probe_default.plt"
  Output = "runtime_element_avalanche_probe_default"
}
Physics {
  Mobility(DopingDependence HighFieldSaturation)
  Recombination(Avalanche(VanOverstraeten))
}
Plot {
  hAlphaAvalanche
}
Math {
}
Solve {
}
"""

TCL = """\
proc tcl_cp_Compute_Plot_Values {} {
    set data [$tcl_cp_adr Data]
    set eqfp [$data ReadScalar $::des_data_vertex "eQuasiFermiPotential"]
    set qfp0 [tcl_cp_get_double $eqfp 0]

    set target ""
    foreach candidate {-1 -10 -20} {
        if {[expr {abs($qfp0-$candidate)}] < 1.0e-8} {
            set target $candidate
        }
    }
    puts [format "AVAL_PROBE_BEGIN bias_V=%d" $target]
}
"""


class HighBiasProcessProbeVmTest(unittest.TestCase):
    def test_exact_decimal_bias_and_process_fields_are_frozen(self) -> None:
        deck = process_deck(DECK)
        self.assertIn(f'Voltage={BIAS_V}', deck)
        for field in PROCESS_FIELDS:
            self.assertIn(field, deck)

    def test_tcl_logs_decimal_bias_without_integer_truncation(self) -> None:
        tcl = process_tcl(TCL)
        self.assertIn("foreach candidate {-19.95}", tcl)
        self.assertIn("bias_V=%.17g", tcl)
        self.assertNotIn("bias_V=%d", tcl)


if __name__ == "__main__":
    unittest.main()
