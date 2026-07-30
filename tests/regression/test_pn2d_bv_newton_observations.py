from __future__ import annotations

import math
import unittest

from scripts.build_pn2d_bv_newton_chain_inputs import normalized
from scripts.run_pn2d_bv_sentaurus_newton_probe_vm import (
    make_probe_deck,
    source_prefix_from_artifact,
    transition_sources,
)


TEMPLATE = """\
File {
  Grid = "mesh.tdr"
  Plot = "old.tdr"
  Current = "old.plt"
  Output = "old"
}
Electrode {
  { Name="Anode" Voltage=0 }
  { Name="Cathode" Voltage=0 }
}
Physics { Recombination(SRH) }
Math {
  CNormPrint
  NewtonPlot(Error MinError Residual)
  Extrapolate
  Iterations=80
}
Solve {
  Coupled { Poisson Electron Hole }
}
"""


class Pn2dBvNewtonObservationTest(unittest.TestCase):
    def test_fixed_transition_deck_writes_residual_and_raw_update(self) -> None:
        deck = make_probe_deck(
            TEMPLATE,
            source_prefix="/remote/source/snapshot_023",
            target_bias=-19.8,
            stem="probe",
        )
        self.assertIn(
            'Load(FilePrefix="/remote/source/snapshot_023")',
            deck,
        )
        self.assertIn("NewtonPlot(Error Residual Update)", deck)
        self.assertNotIn("MinError", deck)
        self.assertNotIn("Extrapolate", deck)
        self.assertIn("NewtonPlotStep=2", deck)
        self.assertIn("Coupled(Iterations=1)", deck)
        self.assertIn('Goal { Name="Anode" Voltage=-19.800000000000001 }', deck)
        self.assertEqual(deck.count("Quasistationary("), 1)

    def test_transition_source_is_previous_exact_snapshot(self) -> None:
        manifest = {
            "branch_records": [
                {
                    "branch": "avalanche_off",
                    "bias_records": [
                        {
                            "requested_bias_V": -19.5,
                            "snapshot_tdr": {
                                "path": "avalanche_off/raw/previous_des.tdr"
                            },
                        },
                        {
                            "requested_bias_V": -19.7,
                            "snapshot_tdr": {
                                "path": "avalanche_off/raw/target_des.tdr"
                            },
                        },
                    ],
                }
            ]
        }
        result = transition_sources(
            manifest, "avalanche_off", (-19.7,)
        )
        self.assertEqual(result[0][0], -19.7)
        self.assertEqual(result[0][1]["requested_bias_V"], -19.5)
        self.assertEqual(
            source_prefix_from_artifact(
                "/remote/run", result[0][1]["snapshot_tdr"]
            ),
            "/remote/run/avalanche_off/raw/previous",
        )

    def test_spatial_signature_has_unit_l2_norm_and_preserves_sign(self) -> None:
        values = normalized([3.0, -4.0, 0.0])
        self.assertEqual(values, [0.6, -0.8, 0.0])
        self.assertAlmostEqual(
            math.sqrt(sum(value * value for value in values)),
            1.0,
        )
        self.assertEqual(normalized([0.0, 0.0]), [0.0, 0.0])


if __name__ == "__main__":
    unittest.main()
