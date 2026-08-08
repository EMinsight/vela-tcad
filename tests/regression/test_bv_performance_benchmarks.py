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

    def test_parse_gprof_callgraph_writes_cumulative_rows(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            callgraph = root / "callgraph.txt"
            callgraph.write_text(
                "[1]  75.0  1.50  3.00  12  solver::solve() [1]\n"
                "[2]  25.0  1.00  0.00      startup() [2]\n"
                "    0.25  0.50  3/12       child() [3]\n",
                encoding="utf-8",
            )
            output = root / "callgraph.csv"
            self.assertEqual(
                BENCHMARKS.parse_gprof_callgraph(callgraph, output), 2
            )
            with output.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(rows[0]["function"], "solver::solve()")
            self.assertEqual(rows[0]["calls"], "12")
            self.assertAlmostEqual(float(rows[0]["cumulative_seconds"]), 4.5)
            self.assertEqual(rows[1]["function"], "startup()")
            self.assertEqual(rows[1]["calls"], "")

    def test_summarize_gprof_applies_candidate_thresholds(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "gmon.out").write_bytes(b"profile")
            for name in ("gprof_flat.txt", "gprof_callgraph.txt"):
                (root / name).write_text("report\n", encoding="utf-8")
            with (root / "gprof_hotspots.csv").open(
                "w", newline="", encoding="utf-8"
            ) as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=[
                        "percent_time", "cumulative_seconds", "self_seconds",
                        "calls", "self_seconds_per_call",
                        "total_seconds_per_call", "function",
                    ],
                )
                writer.writeheader()
                writer.writerows([
                    {"percent_time": "20", "function": "_mcount_private"},
                    {"percent_time": "6", "calls": "50", "function": "hot()"},
                    {"percent_time": "4", "calls": "100", "function": "frequent()"},
                ])
            with (root / "gprof_callgraph_hotspots.csv").open(
                "w", newline="", encoding="utf-8"
            ) as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=[
                        "index", "percent_time", "self_seconds",
                        "children_seconds", "cumulative_seconds", "calls",
                        "function",
                    ],
                )
                writer.writeheader()
                writer.writerows([
                    {"index": "1", "percent_time": "12", "function": "hot()"},
                    {"index": "2", "percent_time": "8", "function": "other()"},
                ])
            summary = BENCHMARKS.summarize_gprof(root)
            self.assertIsNotNone(summary)
            assert summary is not None
            self.assertEqual(
                summary["self_time_candidates_over_5_percent"][0]["function"],
                "hot()",
            )
            self.assertEqual(
                summary["cumulative_time_candidates_over_10_percent"][0]["function"],
                "hot()",
            )
            self.assertEqual(summary["top_call_counts"][0]["function"], "frequent()")
            self.assertEqual(summary["profiler_runtime_percent"], 20.0)


if __name__ == "__main__":
    unittest.main()
