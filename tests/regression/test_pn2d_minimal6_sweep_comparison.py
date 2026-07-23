import csv
import hashlib
import json
import os
import tempfile
import unittest
from unittest import mock
from pathlib import Path

from PIL import Image, ImageDraw, PngImagePlugin

import scripts.compare_pn2d_minimal6_diagnostic_sweeps as comparison_module

from scripts.compare_pn2d_minimal6_diagnostic_sweeps import (
    compare_sweeps,
    ratio_record,
    verify_comparison_artifacts,
    write_comparison_package as _write_comparison_package,
)
from scripts.pn2d_minimal6_diagnostics.schemas import (
    validate_bv_comparison_v1,
    validate_sweep_manifest_v1,
)
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
    deepest = min(
        [0, *(int(float(row["target_bias_V"])) for row in accepted),
         *(int(float(row["target_bias_V"])) for row in failures)]
    )
    targets = [float(-index) for index in range(abs(deepest) + 1)]
    deck_row = {
        "deck": f"decks/{solver}.deck",
        "deck_sha256": digest(f"{solver}-deck"),
    }
    deck_rows = [deck_row]
    if failures:
        deck_rows = []
        for index, row in enumerate(failures):
            label = (
                f"{solver}-{row['topology']}-{float(row['start_bias_V'])}-"
                f"{float(row['target_bias_V'])}-deck"
            )
            deck_rows.append({
                "solver": solver, "topology": row["topology"],
                "start_bias_V": row["start_bias_V"], "target_bias_V": row["target_bias_V"],
                "deck": f"decks/{solver}-failure-{index}.deck",
                "deck_sha256": digest(label), "status": "rejected",
            })
    return {
        "schema": "vela.pn2d_minimal6_sweep_manifest.v1",
        "diagnostic_disclaimer": "minimal6 diagnostic sweep; not a physical BV curve",
        "template": {"path": "template.json", "sha256": digest(f"{solver}-template")},
        "topology_input_sha256": {
            "sketch": {"mesh.json": digest(f"{solver}-sketch-mesh")},
            "mirror": {"mesh.json": digest(f"{solver}-mirror-mesh")},
        },
        "accepted_checkpoints": accepted,
        "failed_transitions": list(failures),
        "failed_transition": list(failures)[0] if failures else None,
        "segments": deck_rows if solver == "vela" else [],
        "sentaurus_segments": deck_rows if solver == "sentaurus" else [],
        "targets_V": targets,
        "interpolation": "forbidden",
        "branch_threshold_version": BRANCH_THRESHOLD_VERSION,
    }


def phase_a_report():
    states = [
        {
            "topology_id": topology, "requested_bias_V": bias,
            "actual_bias_V": bias, "status": "passed",
        }
        for topology in ("sketch", "mirror") for bias in (0.0, -12.0, -19.0)
    ]
    paths = []
    residuals = []
    for state in states:
        topology, bias = state["topology_id"], state["requested_bias_V"]
        paths.append({
            "topology": topology, "bias_V": bias,
            "dependency_order": ["mobility"],
            "forward": {"order": ["mobility"], "contributions": []},
            "reverse": {"order": ["mobility"], "contributions": []},
            "interactions": [], "native_gap_dex": 0.0, "residual_dex": 0.0,
            "status": "insufficient_data",
        })
        residuals.append({
            "topology": topology, "bias_V": bias,
            "name": "sentaurus_internal_semantics_residual",
            "classification": "available", "dex": 0.0,
        })
    return {
        "schema": "vela.pn2d_minimal6_formula_difference.v1",
        "diagnostic_disclaimer": "minimal6 diagnostic sweep; not a physical BV curve",
        "input_provenance": {"state_manifest": "fixture/manifest.json"},
        "audit_provenance": {"audit_root": "fixture/audit"},
        "state_matrix": states,
        "row_counts": {"node": 36, "edge": 54, "triangle": 24},
        "waterfall_paths": paths, "interactions": [],
        "dominance_rules": {"status": "insufficient_data"},
        "sentaurus_internal_semantics_residual": residuals,
        "vela_parameter_agreement": [],
        "artifact_hashes": {"state_manifest_sha256": digest("phase-a-state-manifest")},
        "records": [],
    }


def _write_bytes(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(payload)


def _payload_for_hash(expected, labels):
    for label in labels:
        payload = label.encode("utf-8")
        if hashlib.sha256(payload).hexdigest() == expected:
            return payload
    raise AssertionError(f"fixture has no payload for declared hash {expected}")


def materialize_bound_manifest(root, value, solver):
    root.mkdir(parents=True, exist_ok=True)
    template = value["template"]
    _write_bytes(
        root / template["path"],
        _payload_for_hash(template["sha256"], [f"{solver}-template"]),
    )
    for topology, entries in value["topology_input_sha256"].items():
        for name, expected in entries.items():
            _write_bytes(
                root / "inputs" / topology / name,
                _payload_for_hash(expected, [f"{solver}-{topology}-mesh"]),
            )
    deck_key = "segments" if solver == "vela" else "sentaurus_segments"
    for row in value[deck_key]:
        identity_label = (
            f"{solver}-{row.get('topology')}-{float(row.get('start_bias_V'))}-"
            f"{float(row.get('target_bias_V'))}-deck"
            if "target_bias_V" in row and "start_bias_V" in row
            else f"{solver}-deck"
        )
        _write_bytes(
            root / row["deck"],
            _payload_for_hash(row["deck_sha256"], [identity_label]),
        )
    for row in value["accepted_checkpoints"]:
        _write_bytes(
            root / row["state_path"],
            _payload_for_hash(
                row["state_sha256"],
                [f"{row['solver']}-{row['topology']}-{float(row['target_bias_V'])}"],
            ),
        )
    manifest_path = root / "sweep_manifest.json"
    manifest_path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")
    return manifest_path


def bound_input_artifacts(root, vela, sentaurus, fixed):
    root.mkdir(parents=True, exist_ok=True)
    vela_path = materialize_bound_manifest(root / "vela", vela, "vela")
    sentaurus_path = materialize_bound_manifest(
        root / "sentaurus", sentaurus, "sentaurus"
    )
    fixed_path = root / "fixed_state_report.json"
    fixed_path.write_text(json.dumps(fixed, indent=2) + "\n", encoding="utf-8")
    return {
        "vela_manifest": vela_path,
        "sentaurus_manifest": sentaurus_path,
        "fixed_state_report": fixed_path,
    }


def write_comparison_package(
    out_dir, vela, sentaurus, *, fixed_state_report, input_artifacts=None,
):
    fixed = fixed_state_report or phase_a_report()
    artifacts = input_artifacts
    if artifacts is None:
        artifacts = bound_input_artifacts(
            Path(out_dir) / "_input_evidence", vela, sentaurus, fixed,
        )
    return _write_comparison_package(
        Path(out_dir),
        vela,
        sentaurus,
        fixed_state_report=fixed,
        input_artifacts=artifacts,
    )


class SweepComparisonTest(unittest.TestCase):
    def test_comparison_prevalidates_both_manifests_and_rejects_noncanonical_targets(self):
        vela = manifest("vela", [
            checkpoint("vela", "sketch", -1.0, -1.0, 1.0, 2.0, 3.0, 4.0),
        ])
        sentaurus = manifest("sentaurus", [
            checkpoint("sentaurus", "sketch", -1.0, -1.0, 1.0, 2.0, 3.0, 4.0),
        ])
        with mock.patch.object(
            comparison_module,
            "validate_sweep_manifest_v1",
            create=True,
            side_effect=validate_sweep_manifest_v1,
        ) as validator:
            compare_sweeps(vela, sentaurus, fixed_state_report={})
        self.assertEqual(
            [call.args[0] for call in validator.call_args_list],
            [vela, sentaurus],
        )

        malformed_prefix = json.loads(json.dumps(vela))
        malformed_prefix["targets_V"] = [0.0, -2.0]
        with self.assertRaisesRegex(ValueError, "exact integer prefix"):
            compare_sweeps(malformed_prefix, sentaurus, fixed_state_report={})

        missing_target = json.loads(json.dumps(vela))
        missing_target["targets_V"] = [0.0]
        with self.assertRaisesRegex(ValueError, "accepted checkpoint target"):
            compare_sweeps(missing_target, sentaurus, fixed_state_report={})

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
        with self.assertRaisesRegex(ValueError, "noncanonical branch threshold version"):
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
            self.assertEqual(payload["closure"]["status"], "not_applicable")

    def test_package_records_and_verifies_input_and_generated_artifact_hashes(self):
        vela = manifest("vela", [])
        sentaurus = manifest("sentaurus", [checkpoint("sentaurus", "sketch", -1.0, -1.0, 1.0, 2.0, 3.0, 4.0)])
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            fixed = phase_a_report()
            artifacts = bound_input_artifacts(root / "inputs", vela, sentaurus, fixed)
            write_comparison_package(
                root, vela, sentaurus, fixed_state_report=fixed,
                input_artifacts=artifacts,
            )
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
        self.assertEqual(zero_row["branch_classification"], "multiplication_like")
        self.assertFalse(zero_row["branch_ratio_evidence"]["geometric_zero"])
        self.assertEqual(zero_row["terminal_current_ratio"]["classification"], "available")
        self.assertEqual(zero_row["terminal_current_ratio"]["value"], 2.0)

    def test_bias_identity_uses_actual_tolerance_but_requires_exact_target_membership(self):
        vela_row = checkpoint("vela", "sketch", -1.0, -1.0, 1.0, 2.0, 3.0, 4.0)
        vela_row["actual_bias_V"] = -1.0 + 1.0e-12
        sentaurus_row = checkpoint("sentaurus", "sketch", -1.0, -1.0, 1.0, 2.0, 3.0, 4.0)
        self.assertEqual(
            len(compare_sweeps(
                manifest("vela", [vela_row]),
                manifest("sentaurus", [sentaurus_row]),
                fixed_state_report={},
            )["checkpoints"]),
            1,
        )
        vela_row["actual_bias_V"] = -1.0 + 1.1e-12
        with self.assertRaisesRegex(ValueError, "not an exact checkpoint"):
            compare_sweeps(
                manifest("vela", [vela_row]),
                manifest("sentaurus", [sentaurus_row]),
                fixed_state_report={},
            )

        near_target = checkpoint("vela", "sketch", -1.0 + 5.0e-13, -1.0, 1.0, 2.0, 3.0, 4.0)
        with self.assertRaisesRegex(ValueError, "target bias"):
            compare_sweeps(
                manifest("vela", [near_target]),
                manifest("sentaurus", [sentaurus_row]),
                fixed_state_report={},
            )

    def test_ratio_record_rejects_boolean_inputs(self):
        for numerator, denominator in ((True, 1.0), (1.0, False)):
            with self.subTest(numerator=numerator, denominator=denominator):
                with self.assertRaises(ValueError):
                    ratio_record(numerator, denominator)


class SweepComparisonRedBatchTwoTest(unittest.TestCase):
    def test_each_eligible_sweep_log_gap_is_retained_without_decomposition(self):
        vela = manifest("vela", [
            checkpoint("vela", "sketch", -1.0, -2.0, 2.0, 4.0, 5.0, 6.0)
        ])
        sentaurus = manifest("sentaurus", [
            checkpoint("sentaurus", "sketch", -1.0, -1.0, 1.0, 2.0, 2.5, 3.0)
        ])
        closure = compare_sweeps(vela, sentaurus, fixed_state_report={})["checkpoints"][0]["gap_closure"]
        self.assertEqual(closure["status"], "unidentifiable")
        self.assertEqual(closure["tolerance_dex"], 1.0e-10)
        self.assertEqual(
            {gap["quantity"] for gap in closure["gaps"]},
            {"terminal_current", "maximum_field", "native_source", "reconstructed_source"},
        )
        for gap in closure["gaps"]:
            self.assertEqual(gap["decomposition_status"], "unidentifiable")
            self.assertEqual(gap["named_contributions"], [])
            self.assertEqual(gap["residual"]["classification"], "unidentifiable")
            self.assertIsNone(gap["residual"]["value_dex"])
            self.assertIsNone(gap["closure_error_dex"])
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
        fixed = phase_a_report()
        rechecks = compare_sweeps(
            manifest("vela", vela_rows),
            manifest("sentaurus", sentaurus_rows),
            fixed_state_report=fixed,
        )["fixed_state_recheck"]
        self.assertEqual([row["bias_V"] for row in rechecks], [0.0, -12.0, -19.0])
        for row in rechecks:
            self.assertEqual(set(row), {
                "bias_V", "status", "ranking_status", "reason_code", "reason",
                "missing_inputs", "fixed_state_status", "fixed_state_dominant_factor",
                "recheck_basis", "topologies", "self_consistent_states",
            })
            self.assertEqual(row["status"], "unidentifiable")
            self.assertEqual(row["ranking_status"], "unidentifiable")
            self.assertEqual(row["reason_code"], "missing_verified_nonlinear_ledger_input_bundle")
            self.assertEqual(row["missing_inputs"], ["verified_nonlinear_ledger_input_bundle"])
            self.assertEqual(row["fixed_state_status"], "insufficient_data")
            self.assertIsNone(row["fixed_state_dominant_factor"])
            self.assertEqual(
                row["recheck_basis"],
                "hash_addressed_self_consistent_states_without_verified_nonlinear_ledger_bundle",
            )
            self.assertNotIn("dominant_factor", row)
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
    def test_zero_or_nonconserved_terminal_current_has_typed_unavailable_claims(self):
        sentaurus = manifest("sentaurus", [
            checkpoint("sentaurus", "sketch", -1.0, -1.0, 1.0, 2.0, 3.0, 4.0)
        ])
        vela_zero = manifest("vela", [
            checkpoint("vela", "sketch", -1.0, 0.0, 0.0, 4.0, 5.0, 6.0)
        ])
        zero_row = compare_sweeps(vela_zero, sentaurus, fixed_state_report={})["checkpoints"][0]
        self.assertEqual(zero_row["branch_classification"], "unidentified")
        self.assertEqual(zero_row["terminal_current_ratio"]["classification"], "zero_numerator")
        self.assertEqual(zero_row["terminal_current_sign_alignment"]["classification"], "zero_current")

        sentaurus["accepted_checkpoints"][0]["observables"]["anode_current_A_per_um"] = 0.0
        vela = manifest("vela", [
            checkpoint("vela", "sketch", -1.0, -2.0, 2.0, 4.0, 5.0, 6.0)
        ])
        gated_row = compare_sweeps(vela, sentaurus, fixed_state_report={})["checkpoints"][0]
        unavailable = {"classification": "unavailable", "value": None, "reason": "contact_current_not_conserved"}
        self.assertEqual(gated_row["branch_classification"], "unidentified")
        self.assertEqual(gated_row["terminal_current_ratio"], unavailable)
        self.assertEqual(gated_row["terminal_current_sign_alignment"], unavailable)
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
        self.assertEqual(geometric_growth["classification"], "available")
        self.assertEqual(geometric_growth["value"], 2.0)

    def test_embedded_self_consistent_ledger_summary_is_rejected_as_unverified(self):
        bias = -12.0
        vela_rows = []
        sentaurus_rows = []
        for solver, rows, current in (("vela", vela_rows, -2.0), ("sentaurus", sentaurus_rows, -1.0)):
            for topology in ("sketch", "mirror"):
                row = checkpoint(solver, topology, bias, current, -current, 4.0, 5.0, 6.0)
                row["quantity_ledger_result"] = {
                    "state_sha256": row["state_sha256"],
                    "status": "available", "dominant_factor": "gradient_recovery",
                    "ranking": ["gradient_recovery", "mobility"],
                    "closure": {"status": "closed", "tolerance_dex": 1.0e-10, "closure_error_dex": 0.0},
                }
                rows.append(row)
        with self.assertRaisesRegex(ValueError, "unverified embedded scientific summary"):
            compare_sweeps(
                manifest("vela", vela_rows), manifest("sentaurus", sentaurus_rows),
                fixed_state_report=phase_a_report(),
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
            self.assertNotIn("\u2014", markdown)
            self.assertNotIn("\u0431\u043a", markdown)
            self.assertNotIn("\u9225?", markdown)
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
            lambda value: value["checkpoints"][0]["gap_closure"]["gaps"][0].update(
                named_contributions=[{"name": "fabricated", "contribution_dex": 7.0}]
            ),
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
    def test_any_embedded_quantity_ledger_result_is_rejected_before_comparison(self):
        vela_row = checkpoint("vela", "sketch", -12.0, -2.0, 2.0, 4.0, 5.0, 6.0)
        sentaurus_row = checkpoint("sentaurus", "sketch", -12.0, -1.0, 1.0, 2.0, 2.5, 3.0)
        vela_row["quantity_ledger_result"] = {
            "state_sha256": vela_row["state_sha256"],
            "status": "available", "dominant_factor": "gradient_recovery",
            "ranking": ["gradient_recovery"],
            "closure": {"status": "closed", "tolerance_dex": 1.0e-10, "closure_error_dex": 0.0},
        }
        with self.assertRaisesRegex(ValueError, "unverified embedded scientific summary"):
            compare_sweeps(
                manifest("vela", [vela_row]), manifest("sentaurus", [sentaurus_row]),
                fixed_state_report=phase_a_report(),
            )
    def test_package_preflights_canonical_json_input_equality(self):
        vela = manifest("vela", [])
        sentaurus = manifest("sentaurus", [])
        fixed = phase_a_report()
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
            paths = bound_input_artifacts(root / "inputs", vela, sentaurus, fixed)
            result = write_comparison_package(
                root / "package",
                vela,
                sentaurus,
                fixed_state_report=fixed,
                input_artifacts=paths,
            )
            self.assertEqual(set(result["input_artifacts"]), set(expected))
            self.assertTrue(verify_comparison_artifacts(root / "package" / "sweep_comparison.json"))

    def test_package_and_preflight_reject_invalid_nonempty_phase_a_report(self):
        vela = manifest("vela", [])
        sentaurus = manifest("sentaurus", [])
        invalid = phase_a_report()
        del invalid["state_matrix"]
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            direct_out = root / "direct"
            with self.assertRaises(ValueError):
                _write_comparison_package(
                    direct_out, vela, sentaurus, fixed_state_report=invalid,
                )
            self.assertFalse(direct_out.exists())

            fixed_path = root / "fixed_state_report.json"
            fixed_path.write_text(json.dumps(invalid, indent=2) + "\n", encoding="utf-8")
            preflight_out = root / "preflight"
            with self.assertRaises(ValueError):
                _write_comparison_package(
                    preflight_out,
                    vela,
                    sentaurus,
                    fixed_state_report=invalid,
                    input_artifacts={"fixed_state_report": fixed_path},
                )
            self.assertFalse(preflight_out.exists())

class SweepComparisonReviewBatchTwoTest(unittest.TestCase):
    def test_rejected_transition_must_match_one_canonical_segment_and_first_failure(self):
        failure = failed_transition("vela", "sketch", 0.0, -1.0, "solver rejected endpoint")
        sentaurus = manifest("sentaurus", [])
        for mutate in (
            lambda value: value["failed_transitions"][0].update(start_bias_V=123.0),
            lambda value: value["failed_transitions"][0].update(target_bias_V=-0.5),
        ):
            damaged = manifest("vela", [], [failure])
            mutate(damaged)
            damaged["failed_transition"] = damaged["failed_transitions"][0]
            with self.subTest(failure=damaged["failed_transition"]), self.assertRaises(ValueError):
                compare_sweeps(damaged, sentaurus, fixed_state_report={})

        first = failed_transition("vela", "sketch", 0.0, -1.0, "first")
        second = failed_transition("vela", "mirror", 0.0, -1.0, "second")
        damaged = manifest("vela", [], [first, second])
        damaged["failed_transition"] = damaged["failed_transitions"][1]
        with self.assertRaises(ValueError):
            compare_sweeps(damaged, sentaurus, fixed_state_report={})

    def test_topology_native_source_ratio_uses_quantity_specific_geometric_zero(self):
        tiny = 1.0e-286
        vela = manifest("vela", [
            checkpoint("vela", "sketch", -1.0, -2.0, 2.0, 4.0, tiny, 6.0),
            checkpoint("vela", "mirror", -1.0, -1.0, 1.0, 2.0, tiny, 3.0),
        ])
        report = compare_sweeps(vela, manifest("sentaurus", []), fixed_state_report={})
        sensitivity = next(row for row in report["topology_sensitivity"] if row["solver"] == "vela")
        self.assertEqual(
            sensitivity["native_source_sketch_over_mirror"],
            {"classification": "geometric_zero", "value": None},
        )
        validate_bv_comparison_v1(report)
    def test_source_geometric_zero_is_per_quantity_and_does_not_gate_current_field_or_growth(self):
        tiny = 1.0e-286
        vela = manifest("vela", [
            checkpoint("vela", "sketch", -1.0, -2.0, 2.0, 4.0, tiny, 6.0),
            checkpoint("vela", "sketch", -2.0, -8.0, 8.0, 8.0, tiny, tiny),
        ])
        sentaurus = manifest("sentaurus", [
            checkpoint("sentaurus", "sketch", -1.0, -1.0, 1.0, 2.0, tiny, 3.0),
            checkpoint("sentaurus", "sketch", -2.0, -4.0, 4.0, 4.0, tiny, tiny),
        ])
        report = compare_sweeps(vela, sentaurus, fixed_state_report={})
        row = next(item for item in report["checkpoints"] if item["bias_V"] == -1.0)
        self.assertEqual(row["terminal_current_ratio"], {"classification": "available", "value": 2.0})
        self.assertEqual(row["maximum_field_ratio"], {"classification": "available", "value": 2.0})
        self.assertEqual(row["native_source_ratio"], {"classification": "geometric_zero", "value": None})
        self.assertEqual(row["reconstructed_source_ratio"], {"classification": "available", "value": 2.0})
        self.assertFalse(row["branch_ratio_evidence"]["geometric_zero"])
        self.assertEqual(row["branch_classification"], "multiplication_like")
        self.assertEqual(row["vela_one_volt_current_growth"], {"classification": "available", "value": 4.0})
        self.assertEqual(row["sentaurus_one_volt_current_growth"], {"classification": "available", "value": 4.0})

    def test_nonconserved_contacts_gate_only_current_claims(self):
        vela = manifest("vela", [
            checkpoint("vela", "sketch", -1.0, -2.0, 0.0, 4.0, 5.0, 6.0),
            checkpoint("vela", "sketch", -2.0, -8.0, 8.0, 8.0, 10.0, 12.0),
            checkpoint("vela", "mirror", -1.0, -1.0, 1.0, 2.0, 2.5, 3.0),
        ])
        sentaurus = manifest("sentaurus", [
            checkpoint("sentaurus", "sketch", -1.0, -1.0, 1.0, 2.0, 2.5, 3.0),
            checkpoint("sentaurus", "sketch", -2.0, -4.0, 4.0, 4.0, 5.0, 6.0),
            checkpoint("sentaurus", "mirror", -1.0, -1.0, 1.0, 2.0, 2.5, 3.0),
        ])
        report = compare_sweeps(vela, sentaurus, fixed_state_report={})
        row = next(item for item in report["checkpoints"] if item["topology"] == "sketch" and item["bias_V"] == -1.0)
        unavailable = {"classification": "unavailable", "value": None, "reason": "contact_current_not_conserved"}
        self.assertEqual(row["contact_current_conservation"]["vela"]["classification"], "not_conserved")
        self.assertEqual(row["branch_classification"], "unidentified")
        self.assertEqual(row["terminal_current_ratio"], unavailable)
        self.assertEqual(row["terminal_current_sign_alignment"], unavailable)
        self.assertEqual(row["vela_one_volt_current_growth"], unavailable)
        self.assertEqual(row["maximum_field_ratio"]["classification"], "available")
        self.assertEqual(row["native_source_ratio"]["classification"], "available")
        vela_topology = next(item for item in report["topology_sensitivity"] if item["solver"] == "vela")
        self.assertEqual(vela_topology["terminal_current_sketch_over_mirror"], unavailable)
        self.assertEqual(vela_topology["maximum_field_sketch_over_mirror"]["classification"], "available")
        self.assertEqual(vela_topology["native_source_sketch_over_mirror"]["classification"], "available")
        validate_bv_comparison_v1(report)

    def test_observed_gaps_are_unidentifiable_not_tautologically_closed(self):
        vela = manifest("vela", [checkpoint("vela", "sketch", -1.0, -2.0, 2.0, 4.0, 5.0, 6.0)])
        sentaurus = manifest("sentaurus", [checkpoint("sentaurus", "sketch", -1.0, -1.0, 1.0, 2.0, 2.5, 3.0)])
        report = compare_sweeps(vela, sentaurus, fixed_state_report={})
        closure = report["checkpoints"][0]["gap_closure"]
        self.assertEqual(closure["status"], "unidentifiable")
        for gap in closure["gaps"]:
            self.assertEqual(gap["decomposition_status"], "unidentifiable")
            self.assertEqual(gap["named_contributions"], [])
            self.assertEqual(gap["residual"]["classification"], "unidentifiable")
            self.assertIsNone(gap["residual"]["value_dex"])
            self.assertIsNone(gap["closure_error_dex"])
        self.assertEqual(report["closure"], {
            "status": "unidentifiable", "eligible_gaps": 4,
            "decomposed_gaps": 0, "unidentifiable_gaps": 4,
            "rule": "observed positive log gaps are retained without fabricated decomposition",
        })
        validate_bv_comparison_v1(report)
        for mutate in (
            lambda value: value["checkpoints"][0]["gap_closure"]["gaps"][0].update(
                named_contributions=[{"name": "terminal_current_difference", "contribution_dex": 0.3010299956639812}]
            ),
            lambda value: value["checkpoints"][0]["gap_closure"]["gaps"][0]["residual"].update(
                classification="available", value_dex=0.0
            ),
            lambda value: value["closure"].update(unidentifiable_gaps=3),
        ):
            damaged = json.loads(json.dumps(report))
            mutate(damaged)
            with self.assertRaises(ValueError):
                validate_bv_comparison_v1(damaged)

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
            with Image.open(attacked_path) as original:
                original_info = dict(original.info)
            pnginfo = PngImagePlugin.PngInfo()
            for key, value in original_info.items():
                if isinstance(value, str):
                    pnginfo.add_text(key, value)
            replacement = Image.new("RGB", (900, 504), "white")
            draw = ImageDraw.Draw(replacement)
            draw.line((0, 0, 899, 503), fill="black", width=5)
            replacement.save(attacked_path, format="PNG", pnginfo=pnginfo)
            attacked_hash = hashlib.sha256(attacked_path.read_bytes()).hexdigest()
            report_path = root / "sweep_comparison.json"
            payload = json.loads(report_path.read_text(encoding="utf-8"))
            payload["artifact_hashes"][attacked_name] = attacked_hash
            payload["figure_contract"]["figures"][attacked_name]["sha256"] = attacked_hash
            report_path.write_text(json.dumps(payload, sort_keys=True, indent=2) + "\n", encoding="utf-8")
            with self.assertRaises(ValueError):
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
            self.assertEqual(marked.call_count, 10)
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
            fixed = phase_a_report()
            artifacts = bound_input_artifacts(Path(temp) / "shared-inputs", vela, sentaurus, fixed)
            first_report = _write_comparison_package(
                first, vela, sentaurus, fixed_state_report=fixed, input_artifacts=artifacts)
            second_report = _write_comparison_package(
                second, vela, sentaurus, fixed_state_report=fixed, input_artifacts=artifacts)
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
class FixedStateIntegrityReviewTest(unittest.TestCase):
    def test_runtime_validator_recomputes_typed_fixed_state_rows_and_state_identities(self):
        biases = (0.0, -12.0, -19.0)
        vela_rows = [
            checkpoint("vela", topology, bias, -2.0, 2.0, 4.0, 5.0, 6.0)
            for topology in ("sketch", "mirror") for bias in biases
        ]
        sentaurus_rows = [
            checkpoint("sentaurus", topology, bias, -1.0, 1.0, 2.0, 2.5, 3.0)
            for topology in ("sketch", "mirror") for bias in biases
        ]
        report = compare_sweeps(
            manifest("vela", vela_rows),
            manifest("sentaurus", sentaurus_rows),
            fixed_state_report=phase_a_report(),
        )
        validate_bv_comparison_v1(report)
        mutations = (
            lambda value: value["fixed_state_recheck"][0].update(status="available"),
            lambda value: value["fixed_state_recheck"][0].update(
                ranking_status="remains_dominant"
            ),
            lambda value: value["fixed_state_recheck"][0]["self_consistent_states"][0].update(
                state_sha256=digest("forged-state")
            ),
            lambda value: value["fixed_state_recheck"][0].update(missing_inputs=[]),
            lambda value: value["fixed_state_recheck"][0].update(
                reason="forged scientific conclusion"
            ),
        )
        for mutate in mutations:
            with self.subTest(mutate=mutate):
                tampered = json.loads(json.dumps(report))
                mutate(tampered)
                with self.assertRaises(ValueError):
                    validate_bv_comparison_v1(tampered)


class PublishedInputBindingReviewTest(unittest.TestCase):
    def _baseline(self):
        vela = manifest("vela", [
            checkpoint("vela", "sketch", -1.0, -2.0, 2.0, 4.0, 5.0, 6.0)
        ])
        sentaurus = manifest("sentaurus", [
            checkpoint("sentaurus", "sketch", -1.0, -1.0, 1.0, 2.0, 2.5, 3.0)
        ])
        return vela, sentaurus, phase_a_report()

    def test_published_package_requires_exact_three_input_artifact_keys(self):
        vela, sentaurus, fixed = self._baseline()
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            complete = bound_input_artifacts(root / "inputs", vela, sentaurus, fixed)
            variants = (
                {key: path for key, path in complete.items() if key != "fixed_state_report"},
                {**complete, "unexpected": complete["fixed_state_report"]},
            )
            for index, artifacts in enumerate(variants):
                out = root / f"out-{index}"
                with self.subTest(keys=set(artifacts)), self.assertRaises(ValueError):
                    _write_comparison_package(
                        out,
                        vela,
                        sentaurus,
                        fixed_state_report=fixed,
                        input_artifacts=artifacts,
                    )
                self.assertFalse(out.exists())

    def test_configuration_rejects_empty_malformed_and_duplicate_deck_rows(self):
        vela, sentaurus, fixed = self._baseline()
        mutations = (
            lambda value: value.update(segments=[]),
            lambda value: value.update(segments=[None]),
            lambda value: value.update(segments=[value["segments"][0], dict(value["segments"][0])]),
        )
        for mutate in mutations:
            damaged = json.loads(json.dumps(vela))
            mutate(damaged)
            with self.subTest(segments=damaged["segments"]), self.assertRaises(ValueError):
                compare_sweeps(damaged, sentaurus, fixed_state_report=fixed)

    def test_published_package_rechecks_every_live_input_byte(self):
        mutations = (
            ("vela-template", lambda paths, vela, sentaurus: Path(paths["vela_manifest"]).parent / vela["template"]["path"]),
            ("vela-sketch-input", lambda paths, vela, sentaurus: Path(paths["vela_manifest"]).parent / "inputs" / "sketch" / "mesh.json"),
            ("sentaurus-mirror-input", lambda paths, vela, sentaurus: Path(paths["sentaurus_manifest"]).parent / "inputs" / "mirror" / "mesh.json"),
            ("vela-deck", lambda paths, vela, sentaurus: Path(paths["vela_manifest"]).parent / vela["segments"][0]["deck"]),
            ("sentaurus-deck", lambda paths, vela, sentaurus: Path(paths["sentaurus_manifest"]).parent / sentaurus["sentaurus_segments"][0]["deck"]),
            ("accepted-state", lambda paths, vela, sentaurus: Path(paths["vela_manifest"]).parent / vela["accepted_checkpoints"][0]["state_path"]),
        )
        for label, locate in mutations:
            with self.subTest(label=label), tempfile.TemporaryDirectory() as temp:
                vela, sentaurus, fixed = self._baseline()
                root = Path(temp)
                paths = bound_input_artifacts(root / "inputs", vela, sentaurus, fixed)
                locate(paths, vela, sentaurus).write_bytes(b"tampered-live-byte")
                out = root / "out"
                with self.assertRaises(ValueError):
                    _write_comparison_package(
                        out,
                        vela,
                        sentaurus,
                        fixed_state_report=fixed,
                        input_artifacts=paths,
                    )
                self.assertFalse(out.exists())

    def test_published_state_path_rejects_relative_traversal_outside_manifest_root(self):
        vela, sentaurus, fixed = self._baseline()
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            paths = bound_input_artifacts(root / "inputs", vela, sentaurus, fixed)
            outside = Path(paths["vela_manifest"]).parent.parent / "outside-state.json"
            outside.write_bytes(b"outside-state")
            vela["accepted_checkpoints"][0]["state_path"] = "../outside-state.json"
            vela["accepted_checkpoints"][0]["state_sha256"] = hashlib.sha256(
                outside.read_bytes()
            ).hexdigest()
            Path(paths["vela_manifest"]).write_text(
                json.dumps(vela, indent=2) + "\n", encoding="utf-8"
            )
            with self.assertRaises(ValueError):
                _write_comparison_package(
                    root / "out", vela, sentaurus, fixed_state_report=fixed,
                    input_artifacts=paths,
                )

    def test_full_canonical_targets_invoke_existing_strict_task6_validator(self):
        vela, sentaurus, fixed = self._baseline()
        vela["targets_V"] = [float(-index) for index in range(21)]
        sentaurus["targets_V"] = [float(-index) for index in range(21)]
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            paths = bound_input_artifacts(root / "inputs", vela, sentaurus, fixed)
            with mock.patch.object(
                comparison_module, "validate_sweep_manifest", create=True,
            ) as strict:
                _write_comparison_package(
                    root / "out", vela, sentaurus, fixed_state_report=fixed,
                    input_artifacts=paths,
                )
            self.assertEqual(strict.call_count, 4)
            self.assertEqual(
                [call.kwargs["package_root"] for call in strict.call_args_list],
                [
                    Path(paths["vela_manifest"]).parent,
                    Path(paths["sentaurus_manifest"]).parent,
                ] * 2,
            )

    def test_solver_configuration_binds_full_canonical_manifest_hash(self):
        vela, sentaurus, fixed = self._baseline()
        report = compare_sweeps(vela, sentaurus, fixed_state_report=fixed)
        for solver, value in (("vela", vela), ("sentaurus", sentaurus)):
            self.assertEqual(
                report["solver_configurations"][solver]["sweep_manifest_sha256"],
                canonical_digest(value),
            )


class RuntimeAliasIntegrityReviewTest(unittest.TestCase):
    def test_runtime_recomputes_checkpoint_aliases_and_convergence_counts(self):
        vela = manifest("vela", [
            checkpoint("vela", "sketch", -1.0, -2.0, 2.0, 4.0, 5.0, 6.0)
        ])
        sentaurus = manifest("sentaurus", [
            checkpoint("sentaurus", "sketch", -1.0, -1.0, 1.0, 2.0, 2.5, 3.0)
        ])
        report = compare_sweeps(vela, sentaurus, fixed_state_report=phase_a_report())
        for name in ("records", "terminal_currents", "maximum_fields", "source_integrals"):
            with self.subTest(alias=name):
                tampered = json.loads(json.dumps(report))
                tampered[name] = []
                with self.assertRaises(ValueError):
                    validate_bv_comparison_v1(tampered)
        for field in ("vela_accepted", "sentaurus_accepted", "common_exact"):
            with self.subTest(convergence=field):
                tampered = json.loads(json.dumps(report))
                tampered["convergence_metadata"][field] += 1
                with self.assertRaises(ValueError):
                    validate_bv_comparison_v1(tampered)


class StandaloneVerifierIntegrityReviewTest(unittest.TestCase):
    def _package(self, root):
        vela = manifest("vela", [
            checkpoint("vela", "sketch", -1.0, -2.0, 2.0, 4.0, 5.0, 6.0)
        ])
        sentaurus = manifest("sentaurus", [
            checkpoint("sentaurus", "sketch", -1.0, -1.0, 1.0, 2.0, 2.5, 3.0)
        ])
        write_comparison_package(
            root, vela, sentaurus, fixed_state_report=phase_a_report()
        )
        path = root / "sweep_comparison.json"
        return path, json.loads(path.read_text(encoding="utf-8"))

    @staticmethod
    def _write(path, report):
        path.write_text(json.dumps(report, sort_keys=True, indent=2) + "\n", encoding="utf-8")

    def test_verifier_requires_exact_report_filename(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            path, report = self._package(root)
            renamed = root / "renamed-comparison.json"
            self._write(renamed, report)
            with self.assertRaises(ValueError):
                verify_comparison_artifacts(renamed)
            self.assertTrue(path.is_file())

    def test_verifier_rejects_extra_missing_traversal_and_absolute_artifact_keys(self):
        mutations = []
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            path, baseline = self._package(root)
            extra = root / "extra.bin"
            extra.write_bytes(b"extra")
            escape = root.parent / "batch5-escape.bin"
            escape.write_bytes(b"escape")
            absolute = root / "absolute.bin"
            absolute.write_bytes(b"absolute")
            mutations.extend((
                lambda value: value["artifact_hashes"].update(
                    {"extra.bin": hashlib.sha256(extra.read_bytes()).hexdigest()}
                ),
                lambda value: value["artifact_hashes"].pop("sweep_comparison.csv"),
                lambda value: value["artifact_hashes"].update(
                    {"../batch5-escape.bin": hashlib.sha256(escape.read_bytes()).hexdigest()}
                ),
                lambda value: value["artifact_hashes"].update(
                    {str(absolute.resolve()): hashlib.sha256(absolute.read_bytes()).hexdigest()}
                ),
            ))
            try:
                for mutate in mutations:
                    tampered = json.loads(json.dumps(baseline))
                    mutate(tampered)
                    self._write(path, tampered)
                    with self.subTest(mutate=mutate), self.assertRaises(ValueError):
                        verify_comparison_artifacts(path)
            finally:
                escape.unlink(missing_ok=True)

    def test_verifier_requires_exact_absolute_input_records_and_hashes(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            path, baseline = self._package(root)
            variants = []
            relative = json.loads(json.dumps(baseline))
            original = Path(relative["input_artifacts"]["vela_manifest"]["path"])
            relative_copy = Path.cwd() / "relative-vela-manifest-batch5.json"
            relative_copy.write_bytes(original.read_bytes())
            self.addCleanup(relative_copy.unlink, missing_ok=True)
            relative["input_artifacts"]["vela_manifest"]["path"] = relative_copy.name
            variants.append(relative)
            bad_hash = json.loads(json.dumps(baseline))
            bad_hash["input_artifacts"]["vela_manifest"]["sha256"] = digest("forged")
            variants.append(bad_hash)
            missing = json.loads(json.dumps(baseline))
            missing["input_artifacts"].pop("fixed_state_report")
            missing["artifact_hashes"].pop("input:fixed_state_report")
            variants.append(missing)
            extra = json.loads(json.dumps(baseline))
            extra["input_artifacts"]["extra"] = dict(
                extra["input_artifacts"]["fixed_state_report"]
            )
            extra["artifact_hashes"]["input:extra"] = extra["input_artifacts"]["extra"]["sha256"]
            variants.append(extra)
            for tampered in variants:
                self._write(path, tampered)
                with self.subTest(keys=set(tampered["input_artifacts"])), self.assertRaises(ValueError):
                    verify_comparison_artifacts(path)

    def test_verifier_rederives_semantics_from_inputs_after_synchronized_report_rewrite(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            path, report = self._package(root)
            report.update(
                comparison_status="stopped_with_evidence",
                validation_failure={
                    "code": "no_exact_common_checkpoint",
                    "message": "no accepted identity-verified exact checkpoint is common to both solvers",
                },
                deepest_common_bias_V={
                    "classification": "unavailable", "value": None,
                    "reason": "no accepted exact common checkpoint",
                },
                checkpoints=[], records=[], terminal_currents=[],
                maximum_fields=[], source_integrals=[], missing_tails=[],
                side_only_checkpoints=[],
                convergence_metadata={
                    "vela_accepted": 1, "sentaurus_accepted": 1, "common_exact": 0,
                },
                closure={
                    "status": "not_applicable", "eligible_gaps": 0,
                    "decomposed_gaps": 0, "unidentifiable_gaps": 0,
                    "rule": "observed positive log gaps are retained without fabricated decomposition",
                },
            )
            self._write(path, report)
            with self.assertRaises(ValueError):
                verify_comparison_artifacts(path)

    def test_verifier_rejects_rebound_manifest_path_and_hash_with_stale_semantics(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            path, report = self._package(root)
            record = report["input_artifacts"]["vela_manifest"]
            source = Path(record["path"])
            manifest_value = json.loads(source.read_text(encoding="utf-8"))
            manifest_value["template"]["path"] = str(
                (source.parent / manifest_value["template"]["path"]).resolve()
            )
            rebound = root / "rebound-vela-manifest.json"
            rebound.write_text(json.dumps(manifest_value, indent=2) + "\n", encoding="utf-8")
            rebound_sha = hashlib.sha256(rebound.read_bytes()).hexdigest()
            record.update(path=str(rebound.resolve()), sha256=rebound_sha)
            report["artifact_hashes"]["input:vela_manifest"] = rebound_sha
            self._write(path, report)
            with self.assertRaises(ValueError):
                verify_comparison_artifacts(path)


if __name__ == "__main__":
    unittest.main()
