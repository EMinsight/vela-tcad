#!/usr/bin/env python3
"""Regression coverage for deep-off TransportModels convergence gates."""

from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "run_transportmodels_idvg_srh_strict_sweeps.py"
SPEC = importlib.util.spec_from_file_location("transportmodels_srh_strict", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
WORKFLOW = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(WORKFLOW)


class TransportModelsStrictConfigTest(unittest.TestCase):
    def test_deep_off_checks_are_hard_gates(self) -> None:
        config = {
            "solver": {
                "electron_quantum_potential": {"enabled": True},
            }
        }

        WORKFLOW.strict_solver(config, quantum=True)

        solver = config["solver"]
        self.assertEqual("contact_majority", solver["quasi_fermi_reference"])
        self.assertEqual(2.0e-11, solver["stall_residual_floor"])
        self.assertTrue(solver["line_search"])
        self.assertEqual(1.0, solver["damping_factor"])
        self.assertEqual(2.5e-2, solver["quasi_fermi_update_limit_V"])
        self.assertTrue(solver["carrier_row_qualified_stall_acceptance"])
        rows = solver["carrier_row_convergence"]
        self.assertEqual("enforce", rows["mode"])
        self.assertEqual(1.0e-3, rows["eps_row"])
        self.assertEqual(0.0, rows["min_source_scale_fraction"])
        self.assertEqual(1.0e-18, rows["min_source_scale"])
        closure = solver["global_continuity_closure"]
        self.assertEqual("enforce", closure["mode"])
        self.assertEqual(0.1, closure["tolerance"])
        self.assertEqual(1.0e-18, closure["source_floor"])

    def test_dd_still_disables_quantum_potential(self) -> None:
        config = {
            "solver": {
                "electron_quantum_potential": {"enabled": True},
            }
        }

        WORKFLOW.strict_solver(config, quantum=False)

        self.assertFalse(
            config["solver"]["electron_quantum_potential"]["enabled"]
        )


if __name__ == "__main__":
    unittest.main()
