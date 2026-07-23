import unittest

from scripts.pn2d_minimal6_diagnostics.qfp_sg_replacement import (
    continuity_flux_from_current_proxy,
    density_sg_flux,
    qf_sg_flux,
    replace_internal_qfp,
)


class Minimal6QfpSgReplacementTest(unittest.TestCase):
    def test_qf_sg_matches_density_sg_for_consistent_endpoint_densities(self):
        vt = 0.025851999786435535
        ni = 1.0e16
        psi0, psi1 = -0.2, 0.35
        phin0, phin1 = -0.42, -0.11
        phip0, phip1 = 0.07, 0.43
        coef = 0.04 * vt / 2.0e-6

        n0 = ni * __import__("math").exp((psi0 - phin0) / vt)
        n1 = ni * __import__("math").exp((psi1 - phin1) / vt)
        p0 = ni * __import__("math").exp((phip0 - psi0) / vt)
        p1 = ni * __import__("math").exp((phip1 - psi1) / vt)

        self.assertAlmostEqual(
            qf_sg_flux("electron", ni, ni, psi0, psi1, phin0, phin1, vt, coef),
            density_sg_flux("electron", n0, n1, psi0, psi1, vt, coef),
            delta=abs(density_sg_flux("electron", n0, n1, psi0, psi1, vt, coef))
            * 2.0e-14,
        )
        self.assertAlmostEqual(
            qf_sg_flux("hole", ni, ni, psi0, psi1, phip0, phip1, vt, coef),
            density_sg_flux("hole", p0, p1, psi0, psi1, vt, coef),
            delta=abs(density_sg_flux("hole", p0, p1, psi0, psi1, vt, coef))
            * 2.0e-14,
        )

    def test_replacement_changes_only_requested_internal_qfp(self):
        state = {
            node: {
                "psi_V": 0.1 * node,
                "phin_V": -0.2 * node,
                "phip_V": 0.3 * node,
                "n_m3": 1.0e20 + node,
                "p_m3": 2.0e20 + node,
            }
            for node in range(6)
        }
        sentaurus = {
            node: {"phin_V": -1.0 - node, "phip_V": 2.0 + node}
            for node in range(6)
        }

        replaced = replace_internal_qfp(
            state, sentaurus, replace_electron=True, replace_hole=False
        )
        for node in range(6):
            self.assertEqual(replaced[node]["psi_V"], state[node]["psi_V"])
            self.assertEqual(replaced[node]["phip_V"], state[node]["phip_V"])
            self.assertEqual(replaced[node]["n_m3"], state[node]["n_m3"])
            self.assertEqual(replaced[node]["p_m3"], state[node]["p_m3"])
            expected = sentaurus[node]["phin_V"] if node in (1, 5) else state[node]["phin_V"]
            self.assertEqual(replaced[node]["phin_V"], expected)

    def test_strict_frozen_density_sg_is_qfp_independent_negative_control(self):
        args = ("electron", 2.0e19, 4.0e19, -0.2, 0.4, 0.025, 700.0)
        baseline = density_sg_flux(*args)
        state = {
            node: {
                "psi_V": float(node),
                "phin_V": 0.0,
                "phip_V": 0.0,
                "n_m3": 2.0e19,
                "p_m3": 3.0e19,
            }
            for node in range(6)
        }
        sentaurus = {
            node: {"phin_V": 100.0 + node, "phip_V": -100.0 - node}
            for node in range(6)
        }
        replaced = replace_internal_qfp(
            state, sentaurus, replace_electron=True, replace_hole=True
        )
        after = density_sg_flux(*args)
        self.assertEqual(after, baseline)
        self.assertNotEqual(replaced[1]["phin_V"], state[1]["phin_V"])

    def test_current_proxy_conversion_matches_continuity_conventions(self):
        q = 1.602176634e-19
        self.assertEqual(continuity_flux_from_current_proxy("electron", q), -1.0)
        self.assertEqual(continuity_flux_from_current_proxy("hole", q), 1.0)


if __name__ == "__main__":
    unittest.main()
