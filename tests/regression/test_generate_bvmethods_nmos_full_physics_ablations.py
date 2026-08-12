import tempfile
import unittest
from pathlib import Path

from scripts.generate_bvmethods_nmos_full_physics_ablations import (
    build_deck,
    build_vela_external_resistor_config,
    build_vela_config,
)


TEMPLATE = '''
File{
  Plot      = "n6_des.tdr"
  Parameter = "pp6_des.par"
  Current   = "n6_des.plt"
  Output    = "n6_des.log"
}
Physics{
  Mobility(DopingDep HighFieldsaturation(GradQuasiFermi) Enormal)
  Recombination(SRH(DopingDep) Band2Band(E2) Avalanche(Eparallel))
}
Plot{ eDensity hDensity }
Math { BreakCriteria{ Current(Contact="drain" AbsVal=1.443e-3) } }
Solve { Goal { name=drain current=1.443e-3 } }
'''


class FullPhysicsAblationDeckTests(unittest.TestCase):
    def test_constant_no_enormal_is_explicit(self):
        deck = build_deck(TEMPLATE, "a_constant_no_enormal", False, False)
        self.assertIn('Parameter = "constant_srh.par"', deck)
        self.assertNotIn("Enormal", deck)
        self.assertNotIn("SRH(DopingDep)", deck)
        self.assertIn("current=1e-4", deck.lower())
        self.assertIn("eMobility hMobility", deck)

    def test_full_model_keeps_original_physics(self):
        deck = build_deck(TEMPLATE, "d_doping_enormal", True, True)
        self.assertIn('Parameter = "pp6_des.par"', deck)
        self.assertIn("Enormal", deck)
        self.assertIn("SRH(DopingDep)", deck)
        self.assertIn("AbsVal=1.2e-4", deck)

    def test_vela_ablation_maps_srh_and_lombardi_independently(self):
        base = {
            "solver": {
                "mobility": {"model": "masetti_field"},
                "carrier_row_convergence": {
                    "diagnostic_csv": "old_rows.csv",
                    "trace_csv": "old_trace.csv",
                },
            },
            "sweep": {
                "vtk_prefix": "old/vtk",
                "diagnostics": {
                    "qf_bounds": {"csv_file": "old_qf.csv"},
                    "newton_history": {
                        "csv_file": "old_newton.csv",
                        "attempts_csv_file": "old_attempts.csv",
                        "iterations_csv_file": "old_iterations.csv",
                    },
                },
                "boundary_control": {
                    "evaluation_csv": "old_boundary.csv",
                    "checkpoint_directory": "old_checkpoints",
                },
            },
            "output_csv": "old.csv",
        }
        output = Path("out")
        baseline = build_vela_config(base, output, "a", False, False)
        full = build_vela_config(base, output, "d", True, True)
        self.assertEqual(
            baseline["solver"]["mobility"]["model"], "masetti_field"
        )
        self.assertNotIn("srh_doping_dependence", baseline["solver"])
        self.assertEqual(
            full["solver"]["mobility"]["model"],
            "masetti_field_lombardi",
        )
        self.assertTrue(
            full["solver"]["srh_doping_dependence"]["enabled"]
        )
        self.assertEqual(
            full["solver"]["mobility"]["surface"]["surface_interface"],
            ["R.Substrate", "R.Gateox"],
        )
        case_dir = (output / "d").resolve()
        self.assertEqual(
            Path(full["solver"]["carrier_row_convergence"]["diagnostic_csv"]),
            case_dir / "carrier_row_convergence.csv",
        )
        self.assertEqual(
            Path(full["sweep"]["diagnostics"]["qf_bounds"]["csv_file"]),
            case_dir / "qf_bounds.csv",
        )
        self.assertEqual(
            Path(
                full["sweep"]["diagnostics"]["newton_history"][
                    "iterations_csv_file"
                ]
            ),
            case_dir / "newton_iterations.csv",
        )
        self.assertFalse(full["sweep"]["boundary_control"]["resume"])

    def test_external_resistor_cross_check_replaces_current_control(self):
        full = {
            "output_csv": "old.csv",
            "solver": {
                "carrier_row_convergence": {
                    "diagnostic_csv": "old.csv",
                    "trace_csv": "old.csv",
                }
            },
            "sweep": {
                "voltage_to_current": {"switch_voltage_V": 6.0},
                "continuation": {"predictor": {"mode": "secant"}},
                "boundary_control": {"resume": True},
                "diagnostics": {
                    "qf_bounds": {"csv_file": "old.csv"},
                    "newton_history": {
                        "csv_file": "old.csv",
                        "attempts_csv_file": "old.csv",
                        "iterations_csv_file": "old.csv",
                    },
                },
            },
        }
        case_dir = Path("out/external")
        config = build_vela_external_resistor_config(
            full, case_dir, Path("out/state.csv")
        )
        self.assertNotIn("voltage_to_current", config["sweep"])
        self.assertNotIn("continuation", config["sweep"])
        self.assertEqual(config["sweep"]["bias_points"], [1006.0])
        self.assertEqual(
            config["sweep"]["external_circuit"]["resistance_ohm_um"], 1.0e7
        )
        self.assertEqual(
            config["sweep"]["external_circuit"]["initial_inner_voltage_V"], 6.4069
        )
        self.assertEqual(
            config["sweep"]["external_circuit"]["max_inner_voltage_step_V"],
            0.005,
        )
        self.assertEqual(
            config["sweep"]["external_circuit"]["residual_tolerance_V"], 0.1
        )
        self.assertFalse(config["sweep"]["boundary_control"]["resume"])


if __name__ == "__main__":
    unittest.main()
