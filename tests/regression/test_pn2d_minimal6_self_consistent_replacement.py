import math
import unittest

from scripts.pn2d_minimal6_diagnostics.self_consistent_replacement import (
    THERMAL_VOLTAGE_300K_V,
    _current_edge_key,
    infer_edge_mobility_m2_per_Vs,
    replace_potentials_and_recompute_carriers,
)
from scripts.pn2d_minimal6_diagnostics.qfp_sg_replacement import qf_sg_flux


class Minimal6SelfConsistentReplacementTest(unittest.TestCase):
    def test_current_edge_reference_key_uses_unordered_node_pair(self):
        forward = _current_edge_key("diag_a", -20.0, "electron", 1, 4)
        reverse = _current_edge_key("diag_a", -20.0, "electron", 4, 1)

        self.assertEqual(forward, reverse)
        self.assertEqual(forward, ("diag_a", -20.0, "electron", 1, 4))

    def test_all_potentials_are_replaced_and_carriers_are_recomputed(self):
        vt = THERMAL_VOLTAGE_300K_V
        ni = 1.8e16
        vela = {}
        sentaurus = {}
        for node in range(6):
            psi = -0.2 + 0.01 * node
            phin = psi - 0.05
            phip = psi + 0.04
            vela[node] = {
                "psi_V": psi,
                "phin_V": phin,
                "phip_V": phip,
                "n_m3": ni * math.exp((psi - phin) / vt),
                "p_m3": ni * math.exp((phip - psi) / vt),
            }
            sentaurus[node] = {
                "psi_V": psi + 0.03,
                "phin_V": phin - 0.07,
                "phip_V": phip + 0.08,
                "n_m3": 1.0,
                "p_m3": 2.0,
            }

        replaced, effective_ni, maximum_ni_gap_dex = (
            replace_potentials_and_recompute_carriers(vela, sentaurus)
        )

        self.assertLess(maximum_ni_gap_dex, 1.0e-12)
        for node in range(6):
            self.assertEqual(replaced[node]["psi_V"], sentaurus[node]["psi_V"])
            self.assertEqual(replaced[node]["phin_V"], sentaurus[node]["phin_V"])
            self.assertEqual(replaced[node]["phip_V"], sentaurus[node]["phip_V"])
            self.assertNotEqual(replaced[node]["n_m3"], sentaurus[node]["n_m3"])
            self.assertNotEqual(replaced[node]["p_m3"], sentaurus[node]["p_m3"])
            self.assertAlmostEqual(effective_ni[node], ni, delta=ni * 2.0e-14)
            self.assertAlmostEqual(
                replaced[node]["n_m3"],
                ni
                * math.exp(
                    (
                        sentaurus[node]["psi_V"]
                        - sentaurus[node]["phin_V"]
                    )
                    / vt
                ),
                delta=replaced[node]["n_m3"] * 2.0e-14,
            )
            self.assertAlmostEqual(
                replaced[node]["p_m3"],
                ni
                * math.exp(
                    (
                        sentaurus[node]["phip_V"]
                        - sentaurus[node]["psi_V"]
                    )
                    / vt
                ),
                delta=replaced[node]["p_m3"] * 2.0e-14,
            )

    def test_edge_mobility_is_recovered_from_recomputed_sg_flux(self):
        vt = THERMAL_VOLTAGE_300K_V
        ni = {0: 1.7e16, 1: 1.9e16}
        state = {
            0: {"psi_V": -0.2, "phin_V": -0.3, "phip_V": -0.1},
            1: {"psi_V": 0.4, "phin_V": 0.1, "phip_V": 0.5},
        }
        length = 1.2e-6
        expected_mobility = 0.047
        electron_flux = qf_sg_flux(
            "electron",
            ni[0],
            ni[1],
            state[0]["psi_V"],
            state[1]["psi_V"],
            state[0]["phin_V"],
            state[1]["phin_V"],
            vt,
            expected_mobility * vt / length,
        )

        recovered = infer_edge_mobility_m2_per_Vs(
            carrier="electron",
            node0=0,
            node1=1,
            length_m=length,
            signed_flux_per_m2_s=electron_flux,
            state=state,
            effective_ni_m3=ni,
        )

        self.assertIsNotNone(recovered)
        self.assertAlmostEqual(
            recovered, expected_mobility, delta=expected_mobility * 2.0e-14
        )


if __name__ == "__main__":
    unittest.main()
