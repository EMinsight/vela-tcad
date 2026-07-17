import csv
import hashlib
import json
import tempfile
import unittest
from unittest import mock
from pathlib import Path

from PIL import Image, ImageDraw

import scripts.compare_pn2d_minimal6_diagnostic_sweeps as comparison_module

from scripts.compare_pn2d_minimal6_diagnostic_sweeps import (
    compare_sweeps,
    ratio_record,
    verify_comparison_artifacts,
    write_comparison_package,
)
from scripts.pn2d_minimal6_diagnostics.schemas import validate_bv_comparison_v1
from scripts.run_pn2d_minimal6_diagnostic_sweep import BRANCH_THRESHOLD_VERSION


def digest(label):
    return hashlib.sha256(label.encode("utf-8")).hexdigest()


def canonical_digest(value):
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def checkpoint(solver, topology, bias, anode, cathode, field, native, reconstructed):
    state_path = f"states/{solver}/{topology}/{abs(int(bias)):02d}.json"
    return {
        "solver": solver, "topology": topology, "start_bias_V": bias + 1.0,
        "target_bias_V": bias, "actual_bias_V": bias, "status": "accepted",
        "state_path": state_path,
        "state_sha256": digest(f"{solver}-{topology}-{bias}"),
        "observables": {
            "anode_current_A_per_um": anode,
            "cathode_current_A_per_um": cathode,
            "max_field_V_per_m": field,
            "native_source_integral_s_inv_per_cm": native,
            "reconstructed_source_integral_s_inv_per_cm": reconstructed,
        },
        "branch_classification": "unidentified",
        "branch_threshold_version": BRANCH_THRESHOLD_VERSION,
        "convergence_metadata": {"converged": True, "exit_code": 0},
    }


def failed_transition(solver, topology, start_bias, target_bias, reason):
    return {
        "solver": solver,
        "topology": topology,
        "start_bias_V": start_bias,
        "target_bias_V": target_bias,
        "actual_bias_V": None,
        "exit_code": 7,
        "status": "rejected",
        "state_path": None,
        "state_sha256": None,
        "observables": None,
        "branch_classification": "unidentified",
        "branch_threshold_version": BRANCH_THRESHOLD_VERSION,
        "convergence_metadata": {"converged": False, "exit_code": 7},
        "incomplete_reason": reason,
        "stdout": "",
        "stderr": "synthetic rejected transition",
    }


def manifest(solver, accepted, failures=()):
    return {
        "schema": "vela.pn2d_minimal6_sweep_manifest.v1",
        "diagnostic_disclaimer": "minimal6 diagnostic sweep; not a physical BV curve",
        "template": {"sha256": digest(f"{solver}-template")},
        "topology_input_sha256": {
            "sketch": {"mesh.json": digest(f"{solver}-sketch-mesh")},
            "mirror": {"mesh.json": digest(f"{solver}-mirror-mesh")},
        },
        "accepted_checkpoints": accepted,
        "failed_transitions": list(failures),
        "failed_transition": list(failures)[0] if failures else None,
        "segments": [], "sentaurus_segments": [], "targets_V": [0.0, -1.0, -2.0],
        "interpolation": "forbidden",
        "branch_threshold_version": BRANCH_THRESHOLD_VERSION,
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
        self.assertEqual(row["branch_classification"], "multiplication_like")
        self.assertEqual(row["branch_threshold_version"], BRANCH_THRESHOLD_VERSION)
        self.assertEqual(report["branch_threshold_version"], row["branch_threshold_version"])
        validate_bv_comparison_v1(report)
        tampered = json.loads(json.dumps(report))
        del tampered["checkpoints"][0]["branch_classification"]
        with self.assertRaisesRegex(ValueError, "typed branch classification"):
            validate_bv_comparison_v1(tampered)
        sentaurus["branch_threshold_version"] = "v2"
        with self.assertRaisesRegex(ValueError, "same non-empty branch threshold version"):
            compare_sweeps(vela, sentaurus, fixed_state_report={})

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
        failure = failed_transition("vela", "sketch", 0.0, -1.0, "native source unavailable")
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
            input_path.write_text(json.dumps(vela), encoding="utf-8")
            write_comparison_package(root, vela, sentaurus, fixed_state_report={}, input_artifacts={"vela_manifest": input_path})
            report_path = root / "sweep_comparison.json"
            self.assertTrue(verify_comparison_artifacts(report_path))
            (root / "sweep_comparison.csv").write_text("tampered\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                verify_comparison_artifacts(report_path)

    def test_solver_configuration_exposes_both_solver_deck_hashes(self):
        vela = manifest("vela", [])
        sentaurus = manifest("sentaurus", [])
        vela["segments"] = [{"deck": "vela/sketch.json", "deck_sha256": digest("vela-deck")}]
        sentaurus["sentaurus_segments"] = [{"deck": "sentaurus/sketch.cmd", "deck_sha256": digest("sentaurus-deck")}]
        report = compare_sweeps(vela, sentaurus, fixed_state_report={})
        self.assertEqual(report["solver_configurations"]["vela"]["deck_sha256"], [digest("vela-deck")])
        self.assertEqual(report["solver_configurations"]["sentaurus"]["deck_sha256"], [digest("sentaurus-deck")])


class SweepComparisonRedBatchOneTest(unittest.TestCase):
    def test_manifest_rows_require_typed_branch_and_convergence_evidence(self):
        row = checkpoint("vela", "sketch", -1.0, -1.0, 1.0, 2.0, 3.0, 4.0)
        vela = manifest("vela", [row])
        sentaurus = manifest("sentaurus", [
            checkpoint("sentaurus", "sketch", -1.0, -1.0, 1.0, 2.0, 3.0, 4.0)
        ])
        for missing in ("branch_classification", "branch_threshold_version", "convergence_metadata"):
            with self.subTest(missing=missing):
                tampered = json.loads(json.dumps(vela))
                del tampered["accepted_checkpoints"][0][missing]
                with self.assertRaises(ValueError):
                    compare_sweeps(tampered, sentaurus, fixed_state_report={})

    def test_common_checkpoint_serializes_ratio_evidence_and_geometric_zero(self):
        vela = manifest("vela", [
            checkpoint("vela", "sketch", -1.0, -2.0, 2.0, 4.0, 5.0, 6.0)
        ])
        sentaurus = manifest("sentaurus", [
            checkpoint("sentaurus", "sketch", -1.0, -1.0, 1.0, 2.0, 2.5, 3.0)
        ])
        row = compare_sweeps(vela, sentaurus, fixed_state_report={})["checkpoints"][0]
        self.assertEqual(row["branch_ratio_evidence"], {
            "vela_anode_current_A_per_um": -2.0,
            "sentaurus_anode_current_A_per_um": -1.0,
            "absolute_vela_over_sentaurus": 2.0,
            "geometric_zero": False,
            "threshold_version": BRANCH_THRESHOLD_VERSION,
        })

        zero_vela = manifest("vela", [
            checkpoint("vela", "sketch", -1.0, -2.0, 2.0, 4.0, 0.0, 0.0)
        ])
        zero_row = compare_sweeps(zero_vela, sentaurus, fixed_state_report={})["checkpoints"][0]
        self.assertEqual(zero_row["branch_classification"], "unidentified")
        self.assertTrue(zero_row["branch_ratio_evidence"]["geometric_zero"])
        self.assertEqual(zero_row["terminal_current_ratio"]["classification"], "geometric_zero")
        self.assertIsNone(zero_row["terminal_current_ratio"]["value"])

    def test_bias_identity_uses_tolerance_but_never_interpolates_targets(self):
        vela_row = checkpoint("vela", "sketch", -1.0, -1.0, 1.0, 2.0, 3.0, 4.0)
        vela_row["actual_bias_V"] = -1.0 + 1.0e-12
        sentaurus_row = checkpoint("sentaurus", "sketch", -1.0, -1.0, 1.0, 2.0, 3.0, 4.0)
        self.assertEqual(
            len(compare_sweeps(manifest("vela", [vela_row]), manifest("sentaurus", [sentaurus_row]), fixed_state_report={})["checkpoints"]),
            1,
        )
        vela_row["actual_bias_V"] = -1.0 + 1.1e-12
        with self.assertRaisesRegex(ValueError, "exact target bias"):
            compare_sweeps(manifest("vela", [vela_row]), manifest("sentaurus", [sentaurus_row]), fixed_state_report={})

        near_target = checkpoint("vela", "sketch", -1.0 + 5.0e-13, -1.0, 1.0, 2.0, 3.0, 4.0)
        report = compare_sweeps(
            manifest("vela", [near_target]),
            manifest("sentaurus", [sentaurus_row]),
            fixed_state_report={},
        )
        self.assertEqual(report["checkpoints"], [])
        self.assertEqual(report["deepest_common_bias_V"]["classification"], "unavailable")

    def test_ratio_record_rejects_boolean_inputs(self):
        for numerator, denominator in ((True, 1.0), (1.0, False)):
            with self.subTest(numerator=numerator, denominator=denominator):
                with self.assertRaises(ValueError):
                    ratio_record(numerator, denominator)


class SweepComparisonRedBatchTwoTest(unittest.TestCase):
    def test_each_eligible_sweep_log_gap_has_named_residual_closure(self):
        vela = manifest("vela", [
            checkpoint("vela", "sketch", -1.0, -2.0, 2.0, 4.0, 5.0, 6.0)
        ])
        sentaurus = manifest("sentaurus", [
            checkpoint("sentaurus", "sketch", -1.0, -1.0, 1.0, 2.0, 2.5, 3.0)
        ])
        closure = compare_sweeps(vela, sentaurus, fixed_state_report={})["checkpoints"][0]["gap_closure"]
        self.assertEqual(closure["status"], "closed")
        self.assertEqual(closure["tolerance_dex"], 1.0e-10)
        self.assertEqual(
            {gap["quantity"] for gap in closure["gaps"]},
            {"terminal_current", "maximum_field", "native_source", "reconstructed_source"},
        )
        for gap in closure["gaps"]:
            contribution_sum = sum(item["contribution_dex"] for item in gap["named_contributions"])
            reconstructed = contribution_sum + gap["residual"]["value_dex"]
            self.assertLessEqual(abs(gap["log_gap_dex"] - reconstructed), 1.0e-10)
            self.assertLessEqual(abs(gap["closure_error_dex"]), 1.0e-10)
            self.assertEqual(gap["residual"]["name"], "cross_solver_semantics_residual")

    def test_fixed_state_recheck_binds_exact_self_consistent_states_but_is_honestly_unidentifiable_without_raw_ledger(self):
        biases = (0.0, -12.0, -19.0)
        vela_rows = [
            checkpoint("vela", topology, bias, -2.0, 2.0, 4.0, 5.0, 6.0)
            for topology in ("sketch", "mirror")
            for bias in biases
        ]
        sentaurus_rows = [
            checkpoint("sentaurus", topology, bias, -1.0, 1.0, 2.0, 2.5, 3.0)
            for topology in ("sketch", "mirror")
            for bias in biases
        ]
        fixed = {"root_cause_status": "available", "dominant_factor": "gradient_recovery"}
        rechecks = compare_sweeps(
            manifest("vela", vela_rows),
            manifest("sentaurus", sentaurus_rows),
            fixed_state_report=fixed,
        )["fixed_state_recheck"]
        self.assertEqual([row["bias_V"] for row in rechecks], [0.0, -12.0, -19.0])
        for row in rechecks:
            self.assertEqual(row["status"], "unidentifiable")
            self.assertEqual(row["recheck_basis"], "self_consistent_exact_checkpoints")
            self.assertIn("raw quantity-ledger inputs are unavailable", row["reason"])
            self.assertEqual(len(row["self_consistent_states"]), 4)
            self.assertEqual(
                {(state["solver"], state["topology"]) for state in row["self_consistent_states"]},
                {(solver, topology) for solver in ("vela", "sentaurus") for topology in ("sketch", "mirror")},
            )
            for state in row["self_consistent_states"]:
                self.assertEqual(state["bias_V"], row["bias_V"])
                self.assertTrue(state["state_path"])
                self.assertEqual(len(state["state_sha256"]), 64)
                self.assertEqual(state["state_binding_status"], "manifest_hash_addressed")

    def test_no_common_checkpoint_is_a_typed_stop_with_evidence(self):
        vela = manifest("vela", [])
        sentaurus = manifest("sentaurus", [
            checkpoint("sentaurus", "sketch", -1.0, -1.0, 1.0, 2.0, 3.0, 4.0)
        ])
        report = compare_sweeps(vela, sentaurus, fixed_state_report={})
        self.assertEqual(report["comparison_status"], "stopped_with_evidence")
        self.assertEqual(report["validation_failure"]["code"], "no_exact_common_checkpoint")
        self.assertEqual(report["closure"]["eligible_gaps"], 0)


class SweepComparisonRedBatchTwoAdditionalTest(unittest.TestCase):
    def test_zero_terminal_current_is_unidentified_with_typed_sign_and_ratio(self):
        sentaurus = manifest("sentaurus", [
            checkpoint("sentaurus", "sketch", -1.0, -1.0, 1.0, 2.0, 3.0, 4.0)
        ])
        for vela_current, ratio_classification in ((0.0, "zero_numerator"), (-2.0, "zero_denominator")):
            with self.subTest(vela_current=vela_current):
                sentaurus_row = sentaurus["accepted_checkpoints"][0]
                sentaurus_row["observables"]["anode_current_A_per_um"] = 0.0 if ratio_classification == "zero_denominator" else -1.0
                vela = manifest("vela", [
                    checkpoint("vela", "sketch", -1.0, vela_current, -vela_current, 4.0, 5.0, 6.0)
                ])
                row = compare_sweeps(vela, sentaurus, fixed_state_report={})["checkpoints"][0]
                self.assertEqual(row["branch_classification"], "unidentified")
                self.assertEqual(row["terminal_current_ratio"]["classification"], ratio_classification)
                self.assertEqual(row["terminal_current_sign_alignment"]["classification"], "zero_current")
                self.assertIsNone(row["terminal_current_sign_alignment"]["value"])

    def test_one_volt_growth_requires_two_exact_nonzero_nongeometric_rows(self):
        sentaurus_rows = [
            checkpoint("sentaurus", "sketch", -1.0, -1.0, 1.0, 2.0, 3.0, 4.0),
            checkpoint("sentaurus", "sketch", -2.0, -2.0, 2.0, 3.0, 4.0, 5.0),
        ]
        zero_rows = [
            checkpoint("vela", "sketch", -1.0, 0.0, 0.0, 2.0, 3.0, 4.0),
            checkpoint("vela", "sketch", -2.0, -2.0, 2.0, 3.0, 4.0, 5.0),
        ]
        zero_growth = compare_sweeps(
            manifest("vela", zero_rows), manifest("sentaurus", sentaurus_rows), fixed_state_report={}
        )["checkpoints"][0]["vela_one_volt_current_growth"]
        self.assertEqual(zero_growth["classification"], "zero_current")
        self.assertIsNone(zero_growth["value"])

        geometric_rows = [
            checkpoint("vela", "sketch", -1.0, -1.0, 1.0, 2.0, 0.0, 0.0),
            checkpoint("vela", "sketch", -2.0, -2.0, 2.0, 3.0, 0.0, 0.0),
        ]
        geometric_growth = compare_sweeps(
            manifest("vela", geometric_rows), manifest("sentaurus", sentaurus_rows), fixed_state_report={}
        )["checkpoints"][0]["vela_one_volt_current_growth"]
        self.assertEqual(geometric_growth["classification"], "geometric_zero")
        self.assertIsNone(geometric_growth["value"])

    def test_supplied_self_consistent_ledger_results_validate_ranking_status(self):
        bias = -12.0
        vela_rows = []
        sentaurus_rows = []
        for solver, rows, current in (("vela", vela_rows, -2.0), ("sentaurus", sentaurus_rows, -1.0)):
            for topology in ("sketch", "mirror"):
                row = checkpoint(solver, topology, bias, current, -current, 4.0, 5.0, 6.0)
                row["quantity_ledger_result"] = {
                    "state_sha256": row["state_sha256"],
                    "status": "available",
                    "dominant_factor": "gradient_recovery",
                    "ranking": ["gradient_recovery", "mobility"],
                    "closure": {"status": "closed", "tolerance_dex": 1.0e-10, "closure_error_dex": 0.0},
                }
                rows.append(row)
        fixed = {"root_cause_status": "available", "dominant_factor": "gradient_recovery"}
        recheck = compare_sweeps(
            manifest("vela", vela_rows), manifest("sentaurus", sentaurus_rows), fixed_state_report=fixed
        )["fixed_state_recheck"][1]
        self.assertEqual(recheck["status"], "available")
        self.assertEqual(recheck["ranking_status"], "remains_dominant")
        self.assertEqual(recheck["dominant_factor"], "gradient_recovery")

        vela_rows[0]["quantity_ledger_result"]["state_sha256"] = digest("wrong-state")
        with self.assertRaisesRegex(ValueError, "ledger result state hash"):
            compare_sweeps(
                manifest("vela", vela_rows), manifest("sentaurus", sentaurus_rows), fixed_state_report=fixed
            )

    def test_failure_transition_requires_branch_convergence_and_reason_evidence(self):
        failure = failed_transition("vela", "sketch", 0.0, -1.0, "solver rejected endpoint")
        sentaurus = manifest("sentaurus", [])
        for missing in ("branch_classification", "branch_threshold_version", "convergence_metadata", "incomplete_reason"):
            with self.subTest(missing=missing):
                tampered = json.loads(json.dumps(failure))
                del tampered[missing]
                with self.assertRaises(ValueError):
                    compare_sweeps(manifest("vela", [], [tampered]), sentaurus, fixed_state_report={})

    def test_markdown_is_utf8_clean_and_names_termination_without_extrapolation(self):
        failure = failed_transition("vela", "sketch", 0.0, -1.0, "solver rejected endpoint")
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            write_comparison_package(
                root,
                manifest("vela", [], [failure]),
                manifest("sentaurus", []),
                fixed_state_report={},
            )
            markdown = (root / "sweep_comparison.md").read_text(encoding="utf-8")
            self.assertIn(" V: unidentifiable - ", markdown)
            self.assertNotIn("—", markdown)
            self.assertNotIn("бк", markdown)
            self.assertNotIn("鈥?", markdown)
            self.assertIn("solver rejected endpoint", markdown)
            self.assertIn("not a physical BV curve", markdown)
            self.assertIn("not extrapolated", markdown)


class SweepComparisonRedBatchThreeATest(unittest.TestCase):
    def test_accepted_observables_reject_bool_and_nonfinite(self):
        sentaurus = manifest("sentaurus", [
            checkpoint("sentaurus", "sketch", -1.0, -1.0, 1.0, 2.0, 3.0, 4.0)
        ])
        for invalid in (True, float("inf"), float("nan")):
            with self.subTest(invalid=invalid):
                row = checkpoint("vela", "sketch", -1.0, -1.0, 1.0, 2.0, 3.0, 4.0)
                row["observables"]["anode_current_A_per_um"] = invalid
                with self.assertRaises(ValueError):
                    compare_sweeps(manifest("vela", [row]), sentaurus, fixed_state_report={})

    def test_accepted_checkpoint_requires_state_identity_and_canonical_branch_evidence(self):
        sentaurus = manifest("sentaurus", [
            checkpoint("sentaurus", "sketch", -1.0, -1.0, 1.0, 2.0, 3.0, 4.0)
        ])
        mutations = (
            lambda row: row.update(state_path=""),
            lambda row: row.update(state_sha256="not-a-sha"),
            lambda row: row.update(branch_classification="multiplication_like"),
            lambda row: row.update(branch_threshold_version="noncanonical"),
        )
        for mutate in mutations:
            with self.subTest(mutate=mutate):
                row = checkpoint("vela", "sketch", -1.0, -1.0, 1.0, 2.0, 3.0, 4.0)
                mutate(row)
                with self.assertRaises(ValueError):
                    compare_sweeps(manifest("vela", [row]), sentaurus, fixed_state_report={})

    def test_manifest_configuration_provenance_is_complete_and_hash_addressed(self):
        vela = manifest("vela", [])
        sentaurus = manifest("sentaurus", [])
        for mutate in (
            lambda value: value["template"].update(sha256="bad"),
            lambda value: value["topology_input_sha256"].pop("mirror"),
            lambda value: value["topology_input_sha256"]["sketch"].update({"mesh.json": "bad"}),
        ):
            with self.subTest(mutate=mutate):
                damaged = json.loads(json.dumps(vela))
                mutate(damaged)
                with self.assertRaises(ValueError):
                    compare_sweeps(damaged, sentaurus, fixed_state_report={})
        report = compare_sweeps(vela, sentaurus, fixed_state_report={})
        for solver in ("vela", "sentaurus"):
            config = report["solver_configurations"][solver]
            supplied = config.pop("configuration_sha256")
            self.assertEqual(supplied, canonical_digest(config))
            self.assertEqual(len(supplied), 64)

    def test_common_checkpoint_serializes_and_validates_signed_contact_conservation(self):
        vela = manifest("vela", [
            checkpoint("vela", "sketch", -1.0, -2.0, 2.0, 4.0, 5.0, 6.0)
        ])
        sentaurus = manifest("sentaurus", [
            checkpoint("sentaurus", "sketch", -1.0, -1.0, 1.0, 2.0, 2.5, 3.0)
        ])
        report = compare_sweeps(vela, sentaurus, fixed_state_report={})
        conservation = report["checkpoints"][0]["contact_current_conservation"]
        self.assertEqual(conservation["unit"], "A/um")
        for solver, scale in (("vela", 2.0), ("sentaurus", 1.0)):
            evidence = conservation[solver]
            self.assertEqual(evidence["signed_residual_A_per_um"], 0.0)
            self.assertEqual(evidence["tolerance_A_per_um"], scale * 1.0e-9)
            self.assertEqual(
                evidence["tolerance_formula"],
                "max(1e-18 A/um, 1e-9 * max(|Anode|, |Cathode|))",
            )
            self.assertEqual(evidence["classification"], "conserved")
        validate_bv_comparison_v1(report)
        for mutate in (
            lambda value: value["checkpoints"][0]["contact_current_conservation"]["vela"].update(signed_residual_A_per_um=1.0),
            lambda value: value["checkpoints"][0]["contact_current_conservation"].pop("unit"),
        ):
            with self.subTest(mutate=mutate):
                tampered = json.loads(json.dumps(report))
                mutate(tampered)
                with self.assertRaises(ValueError):
                    validate_bv_comparison_v1(tampered)

class SweepComparisonRedBatchThreeBBranchGapTest(unittest.TestCase):
    def test_validator_independently_recomputes_branch_ratio_evidence(self):
        vela = manifest("vela", [
            checkpoint("vela", "sketch", -1.0, -2.0, 2.0, 4.0, 5.0, 6.0)
        ])
        sentaurus = manifest("sentaurus", [
            checkpoint("sentaurus", "sketch", -1.0, -1.0, 1.0, 2.0, 2.5, 3.0)
        ])
        report = compare_sweeps(vela, sentaurus, fixed_state_report={})
        self.assertIsNone(validate_bv_comparison_v1(report))

        def tamper_nested_current(value):
            row = value["checkpoints"][0]
            row["vela"]["observables"]["anode_current_A_per_um"] = -3.0
            row["vela"]["observables"]["cathode_current_A_per_um"] = 3.0
            conservation = row["contact_current_conservation"]["vela"]
            conservation.update(
                anode_current_A_per_um=-3.0,
                cathode_current_A_per_um=3.0,
                signed_residual_A_per_um=0.0,
                tolerance_A_per_um=3.0e-9,
                classification="conserved",
            )

        mutations = (
            tamper_nested_current,
            lambda value: value["checkpoints"][0]["branch_ratio_evidence"].update(vela_anode_current_A_per_um=-3.0),
            lambda value: value["checkpoints"][0]["branch_ratio_evidence"].update(absolute_vela_over_sentaurus=3.0),
            lambda value: value["checkpoints"][0]["branch_ratio_evidence"].update(geometric_zero=True),
            lambda value: value["checkpoints"][0]["branch_ratio_evidence"].update(threshold_version="tampered"),
            lambda value: value["checkpoints"][0].update(branch_classification="leakage_like"),
        )
        for mutate in mutations:
            with self.subTest(mutate=mutate):
                tampered = json.loads(json.dumps(report))
                mutate(tampered)
                with self.assertRaises(ValueError):
                    validate_bv_comparison_v1(tampered)

    def test_zero_current_branch_evidence_is_typed_and_validated(self):
        cases = (
            (0.0, -1.0, 0.0),
            (-1.0, 0.0, None),
        )
        for vela_current, sentaurus_current, expected_ratio in cases:
            with self.subTest(vela_current=vela_current, sentaurus_current=sentaurus_current):
                vela = manifest("vela", [
                    checkpoint("vela", "sketch", -1.0, vela_current, -vela_current, 4.0, 5.0, 6.0)
                ])
                sentaurus = manifest("sentaurus", [
                    checkpoint("sentaurus", "sketch", -1.0, sentaurus_current, -sentaurus_current, 2.0, 2.5, 3.0)
                ])
                report = compare_sweeps(vela, sentaurus, fixed_state_report={})
                row = report["checkpoints"][0]
                self.assertEqual(row["branch_classification"], "unidentified")
                self.assertEqual(row["branch_ratio_evidence"]["absolute_vela_over_sentaurus"], expected_ratio)
                self.assertIsNone(validate_bv_comparison_v1(report))

    def test_validator_recomputes_exact_gap_set_and_top_level_eligible_count(self):
        vela = manifest("vela", [
            checkpoint("vela", "sketch", -1.0, -2.0, 2.0, 4.0, 5.0, 6.0)
        ])
        sentaurus = manifest("sentaurus", [
            checkpoint("sentaurus", "sketch", -1.0, -1.0, 1.0, 2.0, 2.5, 3.0)
        ])
        report = compare_sweeps(vela, sentaurus, fixed_state_report={})
        self.assertEqual(report["closure"]["eligible_gaps"], 4)
        self.assertIsNone(validate_bv_comparison_v1(report))
        mutations = (
            lambda value: value["checkpoints"][0]["gap_closure"]["gaps"].pop(),
            lambda value: value["checkpoints"][0]["gap_closure"]["gaps"][-1].update(quantity="terminal_current"),
            lambda value: value["checkpoints"][0]["gap_closure"]["gaps"][0]["named_contributions"][0].update(contribution_dex=7.0),
            lambda value: value["checkpoints"][0]["gap_closure"]["gaps"][0]["residual"].update(name="unnamed"),
            lambda value: value["checkpoints"][0]["gap_closure"]["gaps"][0].update(closure_error_dex=1.0e-5),
            lambda value: value["checkpoints"][0]["terminal_current_ratio"].update(classification="zero_numerator", value=0.0),
            lambda value: value["closure"].update(eligible_gaps=3),
        )
        for mutate in mutations:
            with self.subTest(mutate=mutate):
                tampered = json.loads(json.dumps(report))
                mutate(tampered)
                with self.assertRaises(ValueError):
                    validate_bv_comparison_v1(tampered)

class SweepComparisonRedBatchThreeBLedgerInputTest(unittest.TestCase):
    def test_quantity_ledger_result_rejects_inconsistent_ranking_closure_status_and_hash(self):
        bias = -12.0
        vela_row = checkpoint("vela", "sketch", bias, -2.0, 2.0, 4.0, 5.0, 6.0)
        sentaurus_row = checkpoint("sentaurus", "sketch", bias, -1.0, 1.0, 2.0, 2.5, 3.0)
        for row in (vela_row, sentaurus_row):
            row["quantity_ledger_result"] = {
                "state_sha256": row["state_sha256"],
                "status": "available",
                "dominant_factor": "gradient_recovery",
                "ranking": ["gradient_recovery", "mobility"],
                "closure": {
                    "status": "closed",
                    "tolerance_dex": 1.0e-10,
                    "closure_error_dex": 0.0,
                },
            }
        fixed = {"root_cause_status": "available", "dominant_factor": "gradient_recovery"}
        vela_manifest = manifest("vela", [vela_row])
        sentaurus_manifest = manifest("sentaurus", [sentaurus_row])
        report = compare_sweeps(vela_manifest, sentaurus_manifest, fixed_state_report=fixed)
        self.assertEqual(report["fixed_state_recheck"][1]["ranking_status"], "remains_dominant")
        self.assertIsNone(validate_bv_comparison_v1(report))

        mutations = (
            lambda result: result.update(dominant_factor="mobility"),
            lambda result: result["closure"].update(closure_error_dex="not-numeric"),
            lambda result: result["closure"].update(closure_error_dex=1.0e-5),
            lambda result: result.update(status="unavailable"),
            lambda result: result.update(state_sha256=digest("wrong-state")),
        )
        for mutate in mutations:
            with self.subTest(mutate=mutate):
                damaged_vela = json.loads(json.dumps(vela_manifest))
                mutate(damaged_vela["accepted_checkpoints"][0]["quantity_ledger_result"])
                with self.assertRaises(ValueError):
                    compare_sweeps(damaged_vela, sentaurus_manifest, fixed_state_report=fixed)

        tampered_report = json.loads(json.dumps(report))
        tampered_report["accepted_transitions"]["vela"][0]["quantity_ledger_result"]["dominant_factor"] = "mobility"
        with self.assertRaises(ValueError):
            validate_bv_comparison_v1(tampered_report)

    def test_package_preflights_canonical_json_input_equality(self):
        vela = manifest("vela", [])
        sentaurus = manifest("sentaurus", [])
        fixed = {"root_cause_status": "insufficient_data"}
        expected = {
            "vela_manifest": vela,
            "sentaurus_manifest": sentaurus,
            "fixed_state_report": fixed,
        }
        for damaged_name, damaged_payload in (
            ("vela_manifest", sentaurus),
            ("fixed_state_report", "not-json"),
        ):
            with self.subTest(damaged_name=damaged_name), tempfile.TemporaryDirectory() as temp:
                root = Path(temp)
                paths = {}
                for name, payload in expected.items():
                    path = root / f"{name}.json"
                    if name == damaged_name and damaged_payload == "not-json":
                        path.write_text("not-json\n", encoding="utf-8")
                    else:
                        serialized = damaged_payload if name == damaged_name else payload
                        path.write_text(json.dumps(serialized, indent=2) + "\n", encoding="utf-8")
                    paths[name] = path
                out_dir = root / "package"
                with self.assertRaises(ValueError):
                    write_comparison_package(
                        out_dir,
                        vela,
                        sentaurus,
                        fixed_state_report=fixed,
                        input_artifacts=paths,
                    )
                self.assertFalse(out_dir.exists(), "input mismatch must fail before package publication")

        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            paths = {}
            for name, payload in expected.items():
                path = root / f"{name}.json"
                path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
                paths[name] = path
            result = write_comparison_package(
                root / "package",
                vela,
                sentaurus,
                fixed_state_report=fixed,
                input_artifacts=paths,
            )
            self.assertEqual(set(result["input_artifacts"]), set(expected))
            self.assertTrue(verify_comparison_artifacts(root / "package" / "sweep_comparison.json"))

class SweepComparisonRedBatchThreeCFigureContractTest(unittest.TestCase):
    FIGURE_NAMES = {
        "terminal_current.png",
        "one_volt_growth.png",
        "maximum_field.png",
        "source_integrals.png",
        "topology.png",
    }
    SEMANTIC_METADATA = {
        "DiagnosticDisclaimer": "minimal6 diagnostic sweep; not a physical BV curve",
        "SolverTermination": "Every recorded solver failure transition is explicitly marked.",
        "BVExtrapolation": "No physical breakdown voltage (BV) is extrapolated.",
    }

    @staticmethod
    def _full_manifests():
        vela_rows = [
            checkpoint("vela", topology, bias, current, -current, field, native, reconstructed)
            for topology in ("sketch", "mirror")
            for bias, current, field, native, reconstructed in (
                (-1.0, -2.0, 4.0, 5.0, 6.0),
                (-2.0, -8.0, 8.0, 10.0, 12.0),
            )
        ]
        sentaurus_rows = [
            checkpoint("sentaurus", topology, bias, current, -current, field, native, reconstructed)
            for topology in ("sketch", "mirror")
            for bias, current, field, native, reconstructed in (
                (-1.0, -1.0, 2.0, 2.5, 3.0),
                (-2.0, -4.0, 4.0, 5.0, 6.0),
            )
        ]
        vela_failure = failed_transition("vela", "sketch", -2.0, -3.0, "vela terminated")
        sentaurus_failure = failed_transition("sentaurus", "mirror", -2.0, -3.0, "sentaurus terminated")
        return (
            manifest("vela", vela_rows, [vela_failure]),
            manifest("sentaurus", sentaurus_rows, [sentaurus_failure]),
        )

    def test_topology_figure_omits_solver_without_available_ratio(self):
        vela = manifest("vela", [
            checkpoint("vela", "sketch", -1.0, -2.0, 2.0, 4.0, 5.0, 6.0),
            checkpoint("vela", "mirror", -1.0, -1.0, 1.0, 2.0, 2.5, 3.0),
        ])
        sentaurus = manifest("sentaurus", [
            checkpoint("sentaurus", "sketch", -1.0, -1.0, 1.0, 2.0, 2.5, 3.0),
        ])
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            report = write_comparison_package(root, vela, sentaurus, fixed_state_report={})
            self.assertEqual(
                report["figure_contract"]["figures"]["topology.png"]["series_identities"],
                [{"solver": "vela", "quantity": "terminal_current_sketch_over_mirror"}],
            )
            self.assertTrue(verify_comparison_artifacts(root / "sweep_comparison.json"))
    def test_fixed_figure_manifest_and_verifier_rejects_rehashed_replacement(self):
        vela, sentaurus = self._full_manifests()
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            result = write_comparison_package(root, vela, sentaurus, fixed_state_report={})
            contract = result["figure_contract"]
            self.assertEqual(set(contract), {"schema", "figures"})
            self.assertEqual(contract["schema"], "vela.pn2d_minimal6_figure_contract.v1")
            self.assertEqual(set(contract["figures"]), self.FIGURE_NAMES)
            for name, entry in contract["figures"].items():
                with self.subTest(name=name):
                    self.assertEqual(
                        set(entry),
                        {
                            "sha256", "width_px", "height_px", "metadata",
                            "series_identities", "failure_transition_markers",
                        },
                    )
                    self.assertEqual(entry["sha256"], result["artifact_hashes"][name])
                    self.assertEqual((entry["width_px"], entry["height_px"]), (900, 504))
                    self.assertGreaterEqual(entry["width_px"], 640)
                    self.assertGreaterEqual(entry["height_px"], 360)
                    self.assertEqual(entry["metadata"], self.SEMANTIC_METADATA)
                    path = root / name
                    self.assertEqual(path.read_bytes()[:8], b"\x89PNG\r\n\x1a\n")
                    with Image.open(path) as image:
                        image.verify()
                    with Image.open(path) as image:
                        image.load()
                        self.assertEqual(image.format, "PNG")
                        self.assertEqual(image.size, (900, 504))
                        self.assertTrue(
                            any(low != high for low, high in image.convert("RGB").getextrema()),
                            f"{name} must contain nonblank plot pixels",
                        )
                        for key, value in self.SEMANTIC_METADATA.items():
                            self.assertEqual(image.info[key], value)

            attacked_name = "terminal_current.png"
            attacked_path = root / attacked_name
            replacement = Image.new("RGB", (900, 504), "white")
            draw = ImageDraw.Draw(replacement)
            draw.line((0, 0, 899, 503), fill="black", width=5)
            replacement.save(attacked_path, format="PNG")
            attacked_hash = hashlib.sha256(attacked_path.read_bytes()).hexdigest()
            report_path = root / "sweep_comparison.json"
            payload = json.loads(report_path.read_text(encoding="utf-8"))
            payload["artifact_hashes"][attacked_name] = attacked_hash
            payload["figure_contract"]["figures"][attacked_name]["sha256"] = attacked_hash
            report_path.write_text(json.dumps(payload, sort_keys=True, indent=2) + "\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "metadata"):
                verify_comparison_artifacts(report_path)

    def test_all_figures_mark_every_failure_and_sources_keep_both_identities(self):
        self.assertTrue(
            hasattr(comparison_module, "_mark_failure_transitions"),
            "all five plots require one shared failure-transition marker helper",
        )
        vela, sentaurus = self._full_manifests()
        marker_helper = comparison_module._mark_failure_transitions
        with tempfile.TemporaryDirectory() as temp, mock.patch.object(
            comparison_module,
            "_mark_failure_transitions",
            wraps=marker_helper,
        ) as marked:
            root = Path(temp)
            result = write_comparison_package(root, vela, sentaurus, fixed_state_report={})
            self.assertEqual(marked.call_count, 5)
            expected_markers = [
                {
                    "solver": row["solver"],
                    "topology": row["topology"],
                    "start_bias_V": row["start_bias_V"],
                    "target_bias_V": row["target_bias_V"],
                }
                for row in result["failure_transitions"]
            ]
            for name, entry in result["figure_contract"]["figures"].items():
                with self.subTest(name=name):
                    self.assertEqual(entry["failure_transition_markers"], expected_markers)
                    with Image.open(root / name) as image:
                        self.assertEqual(
                            json.loads(image.info["FailureTransitionMarkers"]),
                            expected_markers,
                        )
                        self.assertEqual(
                            json.loads(image.info["SeriesIdentities"]),
                            entry["series_identities"],
                        )

            source_series = result["figure_contract"]["figures"]["source_integrals.png"]["series_identities"]
            self.assertEqual(
                {row["quantity"] for row in source_series},
                {"native_source", "reconstructed_source"},
            )
            for solver in ("vela", "sentaurus"):
                for topology in ("sketch", "mirror"):
                    self.assertIn(
                        {"solver": solver, "topology": topology, "quantity": "native_source"},
                        source_series,
                    )
                    self.assertIn(
                        {"solver": solver, "topology": topology, "quantity": "reconstructed_source"},
                        source_series,
                    )

    def test_one_volt_growth_availability_is_independent_per_solver(self):
        for available_solver in ("vela", "sentaurus"):
            with self.subTest(available_solver=available_solver), tempfile.TemporaryDirectory() as temp:
                vela_rows = [checkpoint("vela", "sketch", -1.0, -2.0, 2.0, 4.0, 5.0, 6.0)]
                sentaurus_rows = [checkpoint("sentaurus", "sketch", -1.0, -1.0, 1.0, 2.0, 2.5, 3.0)]
                if available_solver == "vela":
                    vela_rows.append(checkpoint("vela", "sketch", -2.0, -8.0, 8.0, 8.0, 10.0, 12.0))
                else:
                    sentaurus_rows.append(checkpoint("sentaurus", "sketch", -2.0, -4.0, 4.0, 4.0, 5.0, 6.0))
                result = write_comparison_package(
                    Path(temp),
                    manifest("vela", vela_rows),
                    manifest("sentaurus", sentaurus_rows),
                    fixed_state_report={},
                )
                growth_entry = result["figure_contract"]["figures"]["one_volt_growth.png"]
                self.assertEqual(
                    growth_entry["series_identities"],
                    [{"solver": available_solver, "quantity": "one_volt_growth"}],
                )
                with Image.open(Path(temp) / "one_volt_growth.png") as image:
                    self.assertEqual(
                        json.loads(image.info["SeriesIdentities"]),
                        growth_entry["series_identities"],
                    )

    def test_two_path_independent_writes_are_byte_identical_without_self_hash(self):
        vela, sentaurus = self._full_manifests()
        with tempfile.TemporaryDirectory() as temp:
            first = Path(temp) / "first"
            second = Path(temp) / "second"
            first_report = write_comparison_package(first, vela, sentaurus, fixed_state_report={})
            second_report = write_comparison_package(second, vela, sentaurus, fixed_state_report={})
            self.assertEqual(first_report, second_report)
            self.assertEqual(first_report["artifact_hashes"], second_report["artifact_hashes"])
            self.assertNotIn("sweep_comparison.json", first_report["artifact_hashes"])
            for name in (
                "sweep_comparison.csv",
                "sweep_comparison.md",
                "sweep_comparison.json",
                *sorted(self.FIGURE_NAMES),
            ):
                with self.subTest(name=name):
                    self.assertEqual((first / name).read_bytes(), (second / name).read_bytes())


    def test_runtime_validator_recomputes_fixed_figure_contract(self):
        vela, sentaurus = self._full_manifests()
        with tempfile.TemporaryDirectory() as temp:
            report = write_comparison_package(Path(temp), vela, sentaurus, fixed_state_report={})
            mutations = (
                lambda value: value.pop("figure_contract"),
                lambda value: value["figure_contract"]["figures"]["source_integrals.png"]["series_identities"].pop(),
                lambda value: value["figure_contract"]["figures"]["terminal_current.png"].pop("failure_transition_markers"),
                lambda value: value["figure_contract"]["figures"]["terminal_current.png"]["failure_transition_markers"].append({"solver": "vela", "topology": "sketch", "start_bias_V": 0.0, "target_bias_V": -1.0}),
            )
            for mutate in mutations:
                with self.subTest(mutate=mutate):
                    tampered = json.loads(json.dumps(report))
                    mutate(tampered)
                    with self.assertRaises(ValueError):
                        validate_bv_comparison_v1(tampered)
if __name__ == "__main__":
    unittest.main()
