#!/usr/bin/env python3

from __future__ import annotations

import unittest

from scripts.pn2d_high_bias_process_contract import EXACT_HIGH_BIAS_V
from scripts.run_pn2d_high_bias_oracle_variant_vm import (
    AVALANCHE_DISABLED,
    oracle_deck,
    oracle_tcl,
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
  Recombination(
    Avalanche(VanOverstraeten)
  )
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
    set generation_total [$data ReadScalar $::des_data_vertex "AvalancheGeneration"]
    puts [format "AVAL_PROBE_VERTEX bias_V=%d generation_total_cm3_s=%.17g" \
            $target \
            [tcl_cp_get_double $generation_total $vertex_index]]
    }
}
"""


class HighBiasOracleVmTest(unittest.TestCase):
    def test_avalanche_disabled_branch_removes_only_avalanche_selector(self) -> None:
        deck = oracle_deck(DECK, AVALANCHE_DISABLED, EXACT_HIGH_BIAS_V)
        self.assertNotIn("Avalanche(VanOverstraeten)", deck)
        self.assertIn("Mobility(DopingDependence HighFieldSaturation)", deck)
        for bias in EXACT_HIGH_BIAS_V:
            self.assertIn(f"Voltage={bias}", deck)

    def test_every_exact_state_emits_extended_process_record(self) -> None:
        tcl = oracle_tcl(TCL, EXACT_HIGH_BIAS_V)
        self.assertIn("AVAL_PROBE_PROCESS bias_V=%.17g", tcl)
        self.assertIn("TotalCurrentDensity", tcl)
        self.assertIn("eIonIntegral", tcl)
        self.assertIn("srhRecombination", tcl)
        self.assertIn("-19.95", tcl)


if __name__ == "__main__":
    unittest.main()
