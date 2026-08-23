from __future__ import annotations

import copy
import unittest

from scripts.prepare_slot_ldmos_vela_ialmob_ablation import normalized_physics
from scripts.prepare_slot_ldmos_vela_ialmob_bvds import build_bvds_case
from tests.regression.test_prepare_slot_ldmos_vela_ialmob_ablation import BASE


class SlotLdmosVelaIALMobBVDSTest(unittest.TestCase):
    def test_pair_has_shared_direct_voltage_continuation(self) -> None:
        off = build_bvds_case(copy.deepcopy(BASE), "ialmob_off", 3)
        on = build_bvds_case(copy.deepcopy(BASE), "ialmob_on", 3)
        self.assertEqual(
            off["sweep"]["bias_points"],
            [0.85, 0.90, 0.95, 1.0, 1.5, 2.0, 2.5, 3.0],
        )
        self.assertNotIn("external_circuit", off["sweep"])
        self.assertNotIn("boundary_control", on["sweep"])
        self.assertNotIn("diagnostics", on["sweep"])
        self.assertEqual(off["solver"]["handoff"]["gummel_max_iter"], 0)
        self.assertEqual(normalized_physics(off), normalized_physics(on))


if __name__ == "__main__":
    unittest.main()
