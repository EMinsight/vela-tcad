import csv
import json
import tempfile
import unittest
from pathlib import Path

from scripts.compare_pn2d_minimal6_diagnostic_sweeps import (
    compare_sweeps,
    ratio_record,
    verify_comparison_artifacts,
    write_comparison_package,
)
from scripts.pn2d_minimal6_diagnostics.schemas import validate_bv_comparison_v1


def checkpoint(solver, topology, bias, anode, cathode, field, native, reconstructed):
    return {
        "solver": solver, "topology": topology, "start_bias_V": bias + 1.0,
        "target_bias_V": bias, "actual_bias_V": bias, "status": "accepted",
        "state_sha256": f"{solver}-{topology}-{bias}",
        "observables": {
            "anode_current_A_per_um": anode,
            "cathode_current_A_per_um": cathode,
            "max_field_V_per_m": field,
            "native_source_integral_s_inv_per_cm": native,
            "reconstructed_source_integral_s_inv_per_cm": reconstructed,
        },
    }


def manifest(solver, accepted, failures=()):
    return {
        "schema": "vela.pn2d_minimal6_sweep_manifest.v1",
        "diagnostic_disclaimer": "minimal6 diagnostic sweep; not a physical BV curve",
        "template": {"sha256": f"{solver}-template"},
        "topology_input_sha256": {"sketch": {"mesh.json": f"{solver}-mesh"}},
        "accepted_checkpoints": accepted,
        "failed_transitions": list(failures),
        "failed_transition": list(failures)[0] if failures else None,
        "segments": [], "sentaurus_segments": [], "targets_V": [0.0, -1.0, -2.0],
        "interpolation": "forbidden",
    }


class SweepComparisonTest(unittest.TestCase):
    def test_ratio_record_aligns_terminal_current_sign_and_one_volt_growth(self):
        vela = manifest("vela", [
            checkpoint("vela", "sketch", -1.0, -2.0, 2.0, 4.0, 5.0, 6.0),
            checkpoint("vela", "sketch", -2.0, -8.0, 8.0, 8.0, 10.0, 12.0),
        ])
        sentaurus = manifest("sentaurus", [
            checkpoint("sentaurus", "sketch", -1.0, -1.0, 1.0, 2.0, 2.5, 3.0),
            checkpoint("sentaurus", "sketch", -2.0, -4.0, 4.0, 4.0, 5.0, 6.0),
        ])
        report = compare_sweeps(vela, sentaurus, fixed_state_report={})
        row = report["checkpoints"][0]
        self.assertEqual(row["terminal_current_ratio"]["classification"], "available")
        self.assertEqual(row["terminal_current_ratio"]["value"], 2.0)
        self.assertEqual(row["vela_one_volt_current_growth"]["value"], 4.0)
        self.assertEqual(row["sentaurus_one_volt_current_growth"]["value"], 4.0)
        self.assertEqual(row["maximum_field_ratio"]["value"], 2.0)
        self.assertEqual(row["native_source_ratio"]["value"], 2.0)

    def test_comparison_uses_only_exact_common_biases_and_preserves_missing_tail(self):
        vela = manifest("vela", [checkpoint("vela", "sketch", -1.0, -1.0, 1.0, 2.0, 3.0, 4.0)])
        sentaurus = manifest("sentaurus", [
            checkpoint("sentaurus", "sketch", -1.0, -1.0, 1.0, 2.0, 3.0, 4.0),
            checkpoint("sentaurus", "sketch", -2.0, -2.0, 2.0, 3.0, 4.0, 5.0),
        ])
        report = compare_sweeps(vela, sentaurus, fixed_state_report={})
        self.assertEqual([row["bias_V"] for row in report["checkpoints"]], [-1.0])
        self.assertEqual(report["deepest_common_bias_V"]["value"], -1.0)
        self.assertEqual(report["missing_tails"][0]["solver"], "sentaurus")
        self.assertEqual(report["missing_tails"][0]["biases_V"], [-2.0])

    def test_zero_ratio_is_typed_and_never_divided(self):
        self.assertEqual(ratio_record(1.0, 0.0)["classification"], "zero_denominator")
        self.assertEqual(ratio_record(0.0, 1.0)["classification"], "zero_numerator")
        self.assertEqual(ratio_record(None, 1.0)["classification"], "unavailable")

    def test_empty_common_prefix_reports_failures_and_unidentifiable_fixed_state_recheck(self):
        failure = {"solver": "vela", "topology": "sketch", "start_bias_V": 0.0,
                   "target_bias_V": -1.0, "status": "rejected", "observables": None,
                   "incomplete_reason": "native source unavailable"}
        vela = manifest("vela", [], [failure])
        sentaurus = manifest("sentaurus", [checkpoint("sentaurus", "sketch", -1.0, -1.0, 1.0, 2.0, 3.0, 4.0)])
        fixed = {"root_cause_status": "insufficient_data", "root_cause_reason": "no substitutions"}
        report = compare_sweeps(vela, sentaurus, fixed_state_report=fixed)
        self.assertEqual(report["deepest_common_bias_V"]["classification"], "unavailable")
        self.assertEqual(report["failure_transitions"], [failure])
        self.assertEqual(len(report["fixed_state_recheck"]), 3)
        self.assertTrue(all(row["status"] == "unidentifiable" for row in report["fixed_state_recheck"]))

    def test_topology_sensitivity_uses_common_exact_topology_points(self):
        vela = manifest("vela", [
            checkpoint("vela", "sketch", -1.0, -2.0, 2.0, 4.0, 6.0, 8.0),
            checkpoint("vela", "mirror", -1.0, -1.0, 1.0, 2.0, 3.0, 4.0),
        ])
        sentaurus = manifest("sentaurus", [
            checkpoint("sentaurus", "sketch", -1.0, -1.0, 1.0, 2.0, 3.0, 4.0),
            checkpoint("sentaurus", "mirror", -1.0, -0.5, 0.5, 1.0, 1.5, 2.0),
        ])
        report = compare_sweeps(vela, sentaurus, fixed_state_report={})
        sensitivity = {row["solver"]: row for row in report["topology_sensitivity"]}
        self.assertEqual(sensitivity["vela"]["terminal_current_sketch_over_mirror"]["value"], 2.0)
        self.assertEqual(sensitivity["sentaurus"]["terminal_current_sketch_over_mirror"]["value"], 2.0)

    def test_written_package_has_schema_closure_hashes_and_diagnostic_figures(self):
        vela = manifest("vela", [])
        sentaurus = manifest("sentaurus", [checkpoint("sentaurus", "sketch", -1.0, -1.0, 1.0, 2.0, 3.0, 4.0)])
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            result = write_comparison_package(root, vela, sentaurus, fixed_state_report={})
            validate_bv_comparison_v1(result)
            self.assertTrue((root / "sweep_comparison.csv").is_file())
            self.assertTrue((root / "sweep_comparison.md").is_file())
            self.assertTrue((root / "terminal_current.png").is_file())
            self.assertTrue((root / "topology.png").is_file())
            with (root / "sweep_comparison.csv").open(newline="", encoding="utf-8") as handle:
                self.assertEqual(next(csv.DictReader(handle))["classification"], "side_only")
            payload = json.loads((root / "sweep_comparison.json").read_text(encoding="utf-8"))
            self.assertTrue(payload["artifact_hashes"])
            self.assertEqual(payload["closure"]["status"], "closed")

    def test_package_records_and_verifies_input_and_generated_artifact_hashes(self):
        vela = manifest("vela", [])
        sentaurus = manifest("sentaurus", [checkpoint("sentaurus", "sketch", -1.0, -1.0, 1.0, 2.0, 3.0, 4.0)])
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            input_path = root / "vela_manifest.json"
            input_path.write_text("vela evidence\n", encoding="utf-8")
            write_comparison_package(root, vela, sentaurus, fixed_state_report={}, input_artifacts={"vela_manifest": input_path})
            report_path = root / "sweep_comparison.json"
            self.assertTrue(verify_comparison_artifacts(report_path))
            (root / "sweep_comparison.csv").write_text("tampered\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                verify_comparison_artifacts(report_path)

    def test_solver_configuration_exposes_both_solver_deck_hashes(self):
        vela = manifest("vela", [])
        sentaurus = manifest("sentaurus", [])
        vela["segments"] = [{"deck": "vela/sketch.json", "deck_sha256": "vela-deck"}]
        sentaurus["sentaurus_segments"] = [{"deck": "sentaurus/sketch.cmd", "deck_sha256": "sentaurus-deck"}]
        report = compare_sweeps(vela, sentaurus, fixed_state_report={})
        self.assertEqual(report["solver_configurations"]["vela"]["deck_sha256"], ["vela-deck"])
        self.assertEqual(report["solver_configurations"]["sentaurus"]["deck_sha256"], ["sentaurus-deck"])

if __name__ == "__main__":
    unittest.main()
