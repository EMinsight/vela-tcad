from __future__ import annotations

import unittest
from pathlib import Path

from scripts.build_pn2d_bv_exact_lattice_manifest import (
    bias_token,
    triangle_gradient,
)
from scripts.run_pn2d_bv_exact_lattice_process import (
    branch_config,
    exact_bias_lattice,
    parse_branch_list,
    qualify_rows,
)


class ExactLatticeProcessTests(unittest.TestCase):
    def setUp(self) -> None:
        self.biases = [0.0, -1.0, -19.7, -20.0]
        self.manifest = {
            "branch_records": [
                {"branch": branch, "requested_biases_V": self.biases}
                for branch in (
                    "avalanche_off",
                    "iic_postprocess",
                    "avalanche_on",
                )
            ]
        }
        self.base = {
            "solver": {
                "max_iter": 40,
                "impact_ionization": {
                    "model": "van_overstraeten",
                    "driving_force": "quasi_fermi_gradient",
                },
                "handoff": {"newton_max_iter": 40},
            },
            "sweep": {
                "start": 0.0,
                "stop": -20.0,
                "step": -0.05,
                "contact": "Anode",
            },
        }

    def test_extracts_only_a_common_exact_lattice(self) -> None:
        self.assertEqual(exact_bias_lattice(self.manifest), self.biases)
        self.manifest["branch_records"][1]["requested_biases_V"] = [0.0, -1.0]
        with self.assertRaisesRegex(ValueError, "lattice differs"):
            exact_bias_lattice(self.manifest)

    def test_branch_configs_are_opt_in_and_isolated(self) -> None:
        root = Path("qualification")
        off = branch_config(self.base, "avalanche_off", self.biases, root / "off", 80)
        iic = branch_config(
            self.base, "iic_postprocess", self.biases, root / "iic", 80
        )
        on = branch_config(self.base, "avalanche_on", self.biases, root / "on", 80)
        candidate = branch_config(
            self.base,
            "avalanche_on",
            self.biases,
            root / "candidate",
            80,
            1.0e-2,
        )

        self.assertEqual(self.base["solver"]["max_iter"], 40)
        self.assertEqual(off["solver"]["impact_ionization"], {"model": "none"})
        self.assertNotIn(
            "bv_process_probe", off["sweep"]["diagnostics"]
        )
        self.assertEqual(
            iic["solver"]["impact_ionization"]["coupling_mode"],
            "postprocess_only",
        )
        self.assertEqual(
            on["solver"]["impact_ionization"]["coupling_mode"],
            "self_consistent",
        )
        self.assertNotIn(
            "quasi_fermi_carrier_truncation",
            self.base["solver"]["impact_ionization"],
        )
        self.assertEqual(
            candidate["solver"]["impact_ionization"][
                "quasi_fermi_carrier_truncation"
            ],
            1.0e-2,
        )
        for config in (off, iic, on):
            self.assertEqual(config["solver"]["max_iter"], 80)
            self.assertEqual(
                config["solver"]["handoff"]["newton_max_iter"], 80
            )
            self.assertEqual(config["sweep"]["bias_points"], self.biases)
            self.assertTrue(config["sweep"]["stop_on_failure"])

    def test_qualification_rejects_missing_or_inexact_rows(self) -> None:
        rows = [
            {"bias_V": str(bias), "converged": "1"}
            for bias in self.biases
        ]
        self.assertTrue(
            qualify_rows(rows, self.biases)["complete_exact_lattice"]
        )
        self.assertFalse(
            qualify_rows(rows[:-1], self.biases)["complete_exact_lattice"]
        )
        rows[2]["bias_V"] = "-19.69"
        self.assertFalse(
            qualify_rows(rows, self.biases)["complete_exact_lattice"]
        )

    def test_branch_selection_rejects_unknown_or_duplicate(self) -> None:
        self.assertEqual(
            parse_branch_list("avalanche_on,iic_postprocess"),
            ("avalanche_on", "iic_postprocess"),
        )
        with self.assertRaisesRegex(ValueError, "unknown"):
            parse_branch_list("invalid")
        with self.assertRaisesRegex(ValueError, "duplicate"):
            parse_branch_list("avalanche_on,avalanche_on")

    def test_manifest_helpers_preserve_bias_names_and_linear_gradients(self) -> None:
        self.assertEqual(bias_token(0.0), "0p000000")
        self.assertEqual(bias_token(-19.7), "m19p700000")
        nodes = {0: (0.0, 0.0), 1: (1.0, 0.0), 2: (0.0, 1.0)}
        values = {node: 2.0 * x - 3.0 * y + 4.0 for node, (x, y) in nodes.items()}
        self.assertEqual(
            triangle_gradient([0, 1, 2], nodes, values),
            (2.0, -3.0),
        )


if __name__ == "__main__":
    unittest.main()
