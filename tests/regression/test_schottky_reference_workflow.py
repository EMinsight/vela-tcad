from __future__ import annotations

import importlib.util
from pathlib import Path
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "schottky_reference_workflow",
    ROOT / "scripts" / "run_schottky_reference_workflow.py")
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class SchottkyReferenceWorkflowTest(unittest.TestCase):
    def test_checked_in_curve_passes_declared_gate(self) -> None:
        rows = [
            {"bias_V": float(row["bias_V"]),
             "current_total_A_per_um": float(row["current_total_A_per_um"])}
            for row in MODULE.read_csv(
                MODULE.FIXTURE / "vela_schottky_iv_combined.csv")
        ]
        result = MODULE.compare_curves(rows)
        self.assertEqual("pass", result["status"])
        self.assertEqual(24, result["points_compared"])
        self.assertLessEqual(result["maximum_log10_current_error_dex"], 0.5)

    def test_materialized_second_stage_reuses_exact_0p82_state(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            first, second = MODULE.materialize(Path(directory))
            self.assertTrue(first.is_file())
            config = __import__("json").loads(second.read_text(encoding="utf-8"))
            self.assertTrue(config["sweep"]["continuation"]["arclength"]["enabled"])
            self.assertTrue(config["sweep"]["initial_state_file"].endswith(
                "schottky_bias_0p820000.csv"))
            self.assertEqual(
                "report", config["solver"]["global_continuity_closure"]["mode"])


if __name__ == "__main__":
    unittest.main()
