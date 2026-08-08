from __future__ import annotations

import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[2] / "scripts/run_bv_performance_benchmarks.py"
SPEC = importlib.util.spec_from_file_location("bv_performance_benchmarks", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
BENCHMARKS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(BENCHMARKS)


class BVPerformanceBenchmarksTest(unittest.TestCase):
    def test_voltage_to_current_seed_stops_before_measured_evaluations(self) -> None:
        self.assertTrue(
            BENCHMARKS.keep_seed_evaluation("voltage_to_current_final", 6e-5, 7)
        )
        self.assertTrue(
            BENCHMARKS.keep_seed_evaluation("voltage_to_current_final", 1e-4, 6)
        )
        self.assertFalse(
            BENCHMARKS.keep_seed_evaluation("voltage_to_current_final", 1e-4, 7)
        )

    def test_external_resistor_seed_stops_before_1206_target(self) -> None:
        self.assertTrue(
            BENCHMARKS.keep_seed_evaluation("external_resistor_1206", 1006.0, 7)
        )
        self.assertFalse(
            BENCHMARKS.keep_seed_evaluation("external_resistor_1206", 1206.0, 1)
        )

    def test_extract_bv_interpolates_at_threshold(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sweep.csv"
            with path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=[
                        "converged",
                        "inner_voltage_V",
                        "current_total_A_per_um",
                    ],
                )
                writer.writeheader()
                writer.writerows(
                    [
                        {
                            "converged": "1",
                            "inner_voltage_V": "6.3",
                            "current_total_A_per_um": "8e-5",
                        },
                        {
                            "converged": "1",
                            "inner_voltage_V": "6.5",
                            "current_total_A_per_um": "1.2e-4",
                        },
                    ]
                )
            self.assertAlmostEqual(BENCHMARKS.extract_bv(path), 6.4)

    def test_parse_gprof_flat_writes_machine_readable_rows(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            flat = root / "flat.txt"
            flat.write_text(
                " 50.00  1.00  1.00  10  0.10  0.20  solver::hot()\n"
                " 25.00  1.50  0.50  startup()\n",
                encoding="utf-8",
            )
            output = root / "hotspots.csv"
            self.assertEqual(BENCHMARKS.parse_gprof_flat(flat, output), 2)
            with output.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(rows[0]["function"], "solver::hot()")
            self.assertEqual(rows[0]["calls"], "10")
            self.assertEqual(rows[1]["function"], "startup()")
            self.assertEqual(rows[1]["calls"], "")


if __name__ == "__main__":
    unittest.main()
