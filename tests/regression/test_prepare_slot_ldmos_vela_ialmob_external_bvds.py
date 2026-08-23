from __future__ import annotations

import copy
import unittest

from scripts.prepare_slot_ldmos_vela_ialmob_external_bvds import (
    build_case,
    normalized_external_physics,
)
from tests.regression.test_prepare_slot_ldmos_vela_ialmob_ablation import BASE


EXTERNAL_BASE = copy.deepcopy(BASE)
EXTERNAL_BASE["sweep"].update(
    {
        "external_circuit": {
            "mode": "series_resistor",
            "initial_inner_voltage_V": 38.5,
            "max_inner_voltage_step_V": 0.01,
        },
        "boundary_control": {
            "resume": False,
            "predictor_max_step_factor": 2.0,
            "evaluation_csv": (
                "outputs/stages/06_bvds_external_resistor_final/"
                "boundary_control_evaluations.csv"
            ),
        },
        "diagnostics": {"release_bv_config_audit": {"enabled": True}},
    }
)
EXTERNAL_BASE["output_csv"] = (
    "outputs/stages/06_bvds_external_resistor_final/iv.csv"
)


class SlotLdmosVelaIALMobExternalBVDSTest(unittest.TestCase):
    def test_pair_preserves_external_resistor_and_caps_actual_step(self) -> None:
        off = build_case(copy.deepcopy(EXTERNAL_BASE), "ialmob_off", 0.733)
        on = build_case(copy.deepcopy(EXTERNAL_BASE), "ialmob_on", 0.807)
        self.assertEqual(
            off["sweep"]["external_circuit"]["max_inner_voltage_step_V"], 1.0
        )
        self.assertEqual(
            off["sweep"]["external_circuit"]["solver"], "coupled_newton"
        )
        self.assertEqual(
            off["sweep"]["external_circuit"]["coupled_initial_outer_step_V"],
            25.0,
        )
        self.assertEqual(
            off["sweep"]["external_circuit"]["coupled_min_outer_step_V"],
            0.1,
        )
        self.assertEqual(
            off["sweep"]["external_circuit"]["coupled_max_outer_step_V"],
            5000.0,
        )
        self.assertEqual(
            off["sweep"]["external_circuit"]["coupled_max_step_retries"], 16
        )
        self.assertFalse(
            off["sweep"]["external_circuit"][
                "coupled_apply_device_update_limit"
            ]
        )
        self.assertEqual(
            off["sweep"]["external_circuit"]["coupled_line_search_mode"],
            "residual_filter",
        )
        self.assertEqual(
            off["sweep"]["external_circuit"][
                "coupled_filter_envelope_factor"
            ],
            1.25,
        )
        self.assertEqual(
            off["solver"]["impact_ionization"]["source_jacobian"],
            "local_ad",
        )
        self.assertEqual(
            on["sweep"]["external_circuit"]["max_iterations"], 80
        )
        self.assertEqual(
            on["sweep"]["boundary_control"]["predictor_max_step_factor"], 1.0
        )
        self.assertTrue(off["sweep"]["boundary_control"]["resume"])
        self.assertTrue(
            off["sweep"]["boundary_control"]["adaptive_device_continuation"]
        )
        self.assertTrue(on["sweep"]["diagnostics"]["newton_history"]["enabled"])
        self.assertEqual(off["solver"]["handoff"]["gummel_max_iter"], 50)
        self.assertEqual(
            normalized_external_physics(off), normalized_external_physics(on)
        )


if __name__ == "__main__":
    unittest.main()
