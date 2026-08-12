import csv
import tempfile
import unittest
from pathlib import Path

from scripts.summarize_bvmethods_nmos_full_physics_ablations import (
    VARIANTS,
    summarize,
)


class FullPhysicsAblationSummaryTests(unittest.TestCase):
    def test_extracts_target_rows_and_computes_independent_increments(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            voltages = [6.38, 6.382, 6.383, 6.386]
            fields = [
                "boundary_control_mode",
                "target_current_A_per_um",
                "converged",
                "inner_voltage_V",
                "current_total_A_per_um",
                "current_boundary_residual_A_per_um",
                "boundary_control_evaluations",
                "newton_iterations",
                "global_electron_continuity_closure_ratio",
                "global_hole_continuity_closure_ratio",
            ]
            for name, voltage in zip(VARIANTS, voltages):
                case = root / name
                case.mkdir()
                with (case / "sweep.csv").open("w", newline="", encoding="utf-8") as out:
                    writer = csv.DictWriter(out, fieldnames=fields)
                    writer.writeheader()
                    writer.writerow(
                        {
                            "boundary_control_mode": "current",
                            "target_current_A_per_um": "0.0001",
                            "converged": "1",
                            "inner_voltage_V": str(voltage),
                            "current_total_A_per_um": "0.0001",
                            "current_boundary_residual_A_per_um": "0",
                            "boundary_control_evaluations": "3",
                            "newton_iterations": "2",
                            "global_electron_continuity_closure_ratio": "1e-4",
                            "global_hole_continuity_closure_ratio": "2e-4",
                        }
                    )
            result = summarize(root)
            self.assertEqual(result["status"], "PASS")
            increments = result["voltage_increments_V"]
            self.assertAlmostEqual(increments["srh_doping_dependence_B_minus_A"], 0.002)
            self.assertAlmostEqual(increments["enormal_C_minus_A"], 0.003)
            self.assertAlmostEqual(increments["combined_D_minus_A"], 0.006)
            self.assertAlmostEqual(
                increments["interaction_D_minus_B_minus_C_plus_A"], 0.001
            )
            acceptance = result["sentaurus_full_model_acceptance"]
            self.assertEqual(acceptance["status"], "PASS")
            self.assertLess(acceptance["relative_error"], 0.02)

    def test_external_resistor_cross_check_is_included(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fields = [
                "boundary_control_mode",
                "target_current_A_per_um",
                "converged",
                "inner_voltage_V",
                "current_total_A_per_um",
                "current_boundary_residual_A_per_um",
                "boundary_control_evaluations",
                "newton_iterations",
                "global_electron_continuity_closure_ratio",
                "global_hole_continuity_closure_ratio",
                "load_line_residual_V",
            ]
            for name in VARIANTS:
                case = root / name
                case.mkdir()
                with (case / "sweep.csv").open("w", newline="", encoding="utf-8") as out:
                    writer = csv.DictWriter(out, fieldnames=fields)
                    writer.writeheader()
                    writer.writerow(
                        {
                            "boundary_control_mode": "current",
                            "target_current_A_per_um": "0.0001",
                            "converged": "1",
                            "inner_voltage_V": "6.4",
                            "current_total_A_per_um": "0.0001",
                            "current_boundary_residual_A_per_um": "0",
                            "boundary_control_evaluations": "2",
                            "newton_iterations": "2",
                            "global_electron_continuity_closure_ratio": "0",
                            "global_hole_continuity_closure_ratio": "0",
                        }
                    )
            external = root / "external.csv"
            with external.open("w", newline="", encoding="utf-8") as out:
                writer = csv.DictWriter(out, fieldnames=fields)
                writer.writeheader()
                writer.writerow(
                    {
                        "boundary_control_mode": "external_resistor",
                        "converged": "1",
                        "inner_voltage_V": "6.4001",
                        "current_total_A_per_um": "9.999e-5",
                        "load_line_residual_V": "0.05",
                    }
                )
            result = summarize(root, external_resistor_csv=external)
            cross = result["external_resistor_cross_check"]
            self.assertEqual(cross["status"], "PASS")
            self.assertAlmostEqual(
                cross["absolute_difference_vs_voltage_to_current_V"], 0.0001
            )


if __name__ == "__main__":
    unittest.main()
