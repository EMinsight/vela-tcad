import copy
import sys
import tempfile
from unittest.mock import patch
from pathlib import Path
import unittest

from scripts.run_pn2d_minimal6_diagnostic_sweep import (
    integer_targets,
    classify_branch,
    validate_segment_deck,
    record_transition,
    validate_sweep_manifest,
    copy_topology_inputs,
    execute_segments,
    read_vela_endpoint,
    write_sentaurus_decks,
    main,
)


class DiagnosticSweepTest(unittest.TestCase):
    def test_cli_resolves_relative_output_root_before_running_subprocess(self):
        with tempfile.TemporaryDirectory() as temp:
            with patch("scripts.run_pn2d_minimal6_diagnostic_sweep.initialise_package", return_value=Path(temp) / "sweep_manifest.json") as initialise:
                with patch.object(sys, "argv", ["diagnostic-sweep", "--out-dir", "relative-evidence"]):
                    self.assertEqual(main(), 0)
            self.assertTrue(initialise.call_args.args[0].is_absolute())
    def test_integer_targets_are_exact_and_complete(self):
        self.assertEqual(integer_targets(), tuple(float(-value) for value in range(21)))

    def test_branch_classification_refuses_zero_and_uses_declared_thresholds(self):
        self.assertEqual(classify_branch(1.0, 5.0), "multiplication_like")
        self.assertEqual(classify_branch(1.0, 1.0e-4), "leakage_like")
        self.assertEqual(classify_branch(1.0, 1.0e-2), "unidentified")
        self.assertEqual(classify_branch(0.0, 1.0), "unidentified")

    def test_input_copy_is_separate_from_authoritative_state_root(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            source = root / "authoritative" / "states" / "sketch" / "p0V" / "export"
            source.mkdir(parents=True)
            (source / "mesh.json").write_text("mesh", encoding="utf-8")
            (source / "doping.csv").write_text("doping", encoding="utf-8")
            destination = root / "sweep"
            hashes = copy_topology_inputs(root / "authoritative", destination, ("sketch",))
            self.assertEqual((destination / "inputs" / "sketch" / "mesh.json").read_text(encoding="utf-8"), "mesh")
            self.assertIn("sketch", hashes)
    def test_rejected_transition_keeps_no_fabricated_observables(self):
        manifest = {"schema": "vela.pn2d_minimal6_sweep_manifest.v1", "targets_V": [0.0, -1.0], "segments": [], "accepted_checkpoints": [], "failed_transition": None}
        row = record_transition(manifest, solver="vela", topology="sketch", start_bias_V=0.0, target_bias_V=-1.0, exit_code=9, actual_bias_V=None, state_path=None, observables={"anode_current_A_per_um": 2.0})
        self.assertEqual(row["status"], "rejected")
        self.assertIsNone(row["observables"])
        self.assertEqual(manifest["failed_transition"], row)

    def test_rejected_transition_retains_incomplete_reason_without_partial_observables(self):
        manifest = {"schema": "vela.pn2d_minimal6_sweep_manifest.v1", "targets_V": [0.0, -1.0], "segments": [], "accepted_checkpoints": [], "failed_transition": None}
        row = record_transition(manifest, solver="vela", topology="sketch", start_bias_V=0.0, target_bias_V=-1.0, exit_code=0, actual_bias_V=-1.0, state_path=None, observables=None, incomplete_reason="native source unavailable")
        self.assertEqual(row["incomplete_reason"], "native source unavailable")
        self.assertIsNone(row["observables"])
    def test_sentaurus_decks_are_exact_and_checkpoint_paths_are_unique(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            template = root / "template.cmd"
            template.write_text("Plot=__BIAS_TAG__\nGoal=__TARGET_BIAS_V__\n", encoding="utf-8")
            rows = write_sentaurus_decks(root, template)
            self.assertEqual(len(rows), 42)
            self.assertEqual(len({row["checkpoint_tdr"] for row in rows}), 42)
            final = [row for row in rows if row["topology"] == "mirror" and row["target_bias_V"] == -20.0][0]
            self.assertIn("-20.0", (root / final["deck"]).read_text(encoding="utf-8"))
    def test_vela_endpoint_reader_requires_exact_bias_and_complete_diagnostics(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            curve = root / "curve.csv"
            terminals = root / "terminals.csv"
            curve.write_text("bias_V,max_electric_field_V_per_m\n-1,2\n", encoding="utf-8")
            terminals.write_text("bias_V,contact,I_sgflux_A_per_um,sg_avalanche_source_integral_total\n-1,Anode,-3,5\n-1,Cathode,3,5\n", encoding="utf-8")
            result = read_vela_endpoint(curve, terminals, -1.0)
            self.assertEqual(result["anode_current_A_per_um"], -3.0)
            self.assertEqual(result["reconstructed_source_integral_s_inv_per_cm"], 5.0)
            with self.assertRaises(ValueError): read_vela_endpoint(curve, terminals, -2.0)
    def test_fake_runner_preserves_first_failure_and_stops_later_segments(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            manifest = {"schema": "vela.pn2d_minimal6_sweep_manifest.v1", "targets_V": [0.0, -1.0, -2.0], "segments": [
                {"solver": "vela", "topology": "sketch", "start_bias_V": 0.0, "target_bias_V": -1.0, "status": "pending"},
                {"solver": "vela", "topology": "sketch", "start_bias_V": -1.0, "target_bias_V": -2.0, "status": "pending"}], "accepted_checkpoints": [], "failed_transition": None}
            calls = []
            def runner(segment):
                calls.append(segment["target_bias_V"])
                if len(calls) == 1:
                    state = root / "state_1.csv"
                    state.write_text("state\n", encoding="utf-8")
                    return {"exit_code": 0, "actual_bias_V": -1.0, "state_path": state, "observables": {"anode_current_A_per_um": 1.0, "cathode_current_A_per_um": -1.0, "max_field_V_per_m": 2.0, "native_source_integral_s_inv_per_cm": 3.0, "reconstructed_source_integral_s_inv_per_cm": 4.0}, "stdout": "fake", "stderr": ""}
                return {"exit_code": 7, "actual_bias_V": None, "state_path": None, "observables": None, "stdout": "fake", "stderr": "failure"}
            execute_segments(manifest, root, runner)
            self.assertEqual(calls, [-1.0, -2.0])
            self.assertEqual(manifest["segments"][0]["status"], "accepted")
            self.assertEqual(manifest["failed_transition"]["exit_code"], 7)
            self.assertIsNone(manifest["failed_transition"]["observables"])

    def test_fake_runner_rejects_inexact_bias_and_missing_state(self):
        for label, result in (
            ("inexact", {"exit_code": 0, "actual_bias_V": -0.9, "state_path": None, "observables": None}),
            ("missing_state", {"exit_code": 0, "actual_bias_V": -1.0, "state_path": Path("missing.csv"), "observables": None}),
        ):
            with self.subTest(label=label), tempfile.TemporaryDirectory() as temp:
                root = Path(temp)
                manifest = {"schema": "vela.pn2d_minimal6_sweep_manifest.v1", "targets_V": [0.0, -1.0], "segments": [{"solver": "vela", "topology": "sketch", "start_bias_V": 0.0, "target_bias_V": -1.0, "status": "pending"}], "accepted_checkpoints": [], "failed_transition": None}
                execute_segments(manifest, root, lambda _: result)
                self.assertEqual(manifest["failed_transition"]["status"], "rejected")
                self.assertIsNone(manifest["failed_transition"]["observables"])

    def test_fake_runner_full_success_produces_complete_manifest(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            manifest = {"schema": "vela.pn2d_minimal6_sweep_manifest.v1", "targets_V": [0.0, -1.0], "segments": [
                {"solver": "vela", "topology": "sketch", "start_bias_V": 0.0, "target_bias_V": -1.0, "status": "pending"}], "accepted_checkpoints": [], "failed_transition": None, "interpolation": "forbidden"}
            def runner(segment):
                state = root / f"state_{abs(int(segment['target_bias_V']))}.csv"
                state.write_text("state\n", encoding="utf-8")
                return {"exit_code": 0, "actual_bias_V": segment["target_bias_V"], "state_path": state, "observables": {"anode_current_A_per_um": 1.0, "cathode_current_A_per_um": -1.0, "max_field_V_per_m": 2.0, "native_source_integral_s_inv_per_cm": 3.0, "reconstructed_source_integral_s_inv_per_cm": 4.0}}
            execute_segments(manifest, root, runner)
            validate_sweep_manifest(manifest)
            self.assertEqual(len(manifest["accepted_checkpoints"]), 1)
            self.assertIsNone(manifest["failed_transition"])

    def test_fake_runner_records_sentaurus_early_failure(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            manifest = {"schema": "vela.pn2d_minimal6_sweep_manifest.v1", "targets_V": [0.0, -1.0], "segments": [{"solver": "sentaurus", "topology": "mirror", "start_bias_V": 0.0, "target_bias_V": -1.0, "status": "pending"}], "accepted_checkpoints": [], "failed_transition": None}
            execute_segments(manifest, root, lambda _: {"exit_code": 23, "actual_bias_V": None, "state_path": None, "observables": None, "stderr": "sdevice failure"})
            self.assertEqual(manifest["failed_transition"]["solver"], "sentaurus")
            self.assertEqual(manifest["failed_transition"]["exit_code"], 23)

    def test_manifest_refuses_interpolation(self):
        manifest = {"schema": "vela.pn2d_minimal6_sweep_manifest.v1", "targets_V": [0.0, -1.0], "segments": [], "accepted_checkpoints": [], "failed_transition": None, "interpolation": "linear"}
        with self.assertRaises(ValueError):
            validate_sweep_manifest(manifest)
    def test_manifest_refuses_inexact_or_tampered_accepted_checkpoint(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            state = root / "state.csv"
            state.write_text("state\n", encoding="utf-8")
            manifest = {"schema": "vela.pn2d_minimal6_sweep_manifest.v1", "targets_V": [0.0, -1.0], "segments": [], "accepted_checkpoints": [], "failed_transition": None}
            record_transition(manifest, solver="vela", topology="sketch", start_bias_V=0.0, target_bias_V=-1.0, exit_code=0, actual_bias_V=-1.0, state_path=state, observables={"anode_current_A_per_um": 1.0, "cathode_current_A_per_um": -1.0, "max_field_V_per_m": 2.0, "native_source_integral_s_inv_per_cm": 3.0, "reconstructed_source_integral_s_inv_per_cm": 4.0})
            validate_sweep_manifest(manifest)
            state.write_text("tampered\n", encoding="utf-8")
            with self.assertRaises(ValueError):
                validate_sweep_manifest(manifest)
    def test_deck_contract_allows_only_segment_changes(self):
        template = {"mesh_file": "mesh.json", "node_doping_file": "doping.csv", "solver": {"method": "newton"}, "sweep": {"start": 0.0, "stop": -1.0, "write_state_every_point_prefix": "states/base"}}
        generated = copy.deepcopy(template)
        generated["mesh_file"] = "inputs/sketch/mesh.json"
        generated["node_doping_file"] = "inputs/sketch/doping.csv"
        generated["sweep"].update({"start": -1.0, "stop": -2.0, "initial_state_file": "states/restart.csv", "write_state_every_point_prefix": "states/segment_001"})
        validate_segment_deck(template, generated)
        generated["solver"]["method"] = "gummel"
        with self.assertRaises(ValueError):
            validate_segment_deck(template, generated)


if __name__ == '__main__':
    unittest.main()