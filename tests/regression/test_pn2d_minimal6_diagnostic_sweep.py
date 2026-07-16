import copy
import hashlib
import inspect
import json
import shutil
import sys
import tempfile
from unittest.mock import patch
from pathlib import Path
from types import SimpleNamespace
import unittest

from scripts.run_pn2d_minimal6_diagnostic_sweep import (
    integer_targets,
    classify_branch,
    validate_segment_deck,
    record_transition,
    record_sentaurus_checkpoint,
    record_sentaurus_results,
    validate_sweep_manifest,
    copy_topology_inputs,
    execute_segments,
    run_vela_subprocess_segment,
    segment_state_path,
    read_vela_endpoint,
    read_vela_convergence_metadata,
    read_sentaurus_endpoint,
    initialise_package,
    _required_field_region,
    write_sentaurus_decks,
    main,
)


class DiagnosticSweepTest(unittest.TestCase):
    @staticmethod
    def _sha256(path):
        return hashlib.sha256(path.read_bytes()).hexdigest()

    @staticmethod
    def _validate_package_compat(manifest, root):
        if "package_root" in inspect.signature(validate_sweep_manifest).parameters:
            validate_sweep_manifest(manifest, package_root=root)
        else:
            validate_sweep_manifest(manifest)

    def _complete_package(self, root):
        repo = Path(__file__).resolve().parents[2]
        source_template = repo / "reference_tcad" / "pn2d_sentaurus2018_minimal6" / "vela" / "pn2d_minimal6_sweep_template.json"
        template = root / "immutable_template.json"
        template.write_text(source_template.read_text(encoding="utf-8"), encoding="utf-8")
        fixture = repo / "tests" / "fixtures" / "pn2d_minimal6_synthetic"
        manifest_path = initialise_package(root, template, authoritative_state_root=fixture)
        return json.loads(manifest_path.read_text(encoding="utf-8"))

    def _complete_package_with_common_pair(
        self,
        root,
        *,
        actual_offset_V=0.0,
        sentaurus_current=1.0,
    ):
        manifest = self._complete_package(root)
        target = -1.0
        actual = target + actual_offset_V
        sentaurus_segment = next(
            segment
            for segment in manifest["sentaurus_segments"]
            if segment["topology"] == "sketch"
            and segment["target_bias_V"] == target
        )
        sentaurus_state = root / sentaurus_segment["checkpoint_tdr"]
        sentaurus_state.parent.mkdir(parents=True, exist_ok=True)
        sentaurus_state.write_text("sentaurus state\n", encoding="utf-8")
        vela_state = segment_state_path(root, "sketch", 0.0, target)
        vela_state.parent.mkdir(parents=True, exist_ok=True)
        vela_state.write_text("vela state\n", encoding="utf-8")

        def observables(current):
            return {
                "anode_current_A_per_um": current,
                "cathode_current_A_per_um": -current,
                "max_field_V_per_m": 2.0,
                "native_source_integral_s_inv_per_cm": 3.0,
                "reconstructed_source_integral_s_inv_per_cm": 4.0,
            }

        sentaurus = record_transition(
            manifest,
            solver="sentaurus",
            topology="sketch",
            start_bias_V=0.0,
            target_bias_V=target,
            exit_code=0,
            actual_bias_V=actual,
            state_path=sentaurus_state,
            observables=observables(sentaurus_current),
        )
        vela = record_transition(
            manifest,
            solver="vela",
            topology="sketch",
            start_bias_V=0.0,
            target_bias_V=target,
            exit_code=0,
            actual_bias_V=actual,
            state_path=vela_state,
            observables=observables(5.0),
        )
        sentaurus["state_path"] = str(sentaurus_state.relative_to(root))
        vela["state_path"] = str(vela_state.relative_to(root))
        return manifest, sentaurus, vela

    def test_strict_package_recomputes_common_branch_evidence(self):
        def delete_evidence(index):
            def mutate(manifest, rows):
                del rows[index]["branch_ratio_evidence"]
            return mutate

        def alter_evidence(field, value, *, both=True):
            def mutate(manifest, rows):
                targets = rows if both else rows[:1]
                for row in targets:
                    row["branch_ratio_evidence"][field] = value
            return mutate

        def rows_differ(manifest, rows):
            rows[0]["branch_ratio_evidence"]["absolute_vela_over_sentaurus"] = 6.0

        def wrong_classification(manifest, rows):
            for row in rows:
                row["branch_classification"] = "leakage_like"

        def noncanonical_threshold(manifest, rows):
            manifest["branch_threshold_version"] = "v2-noncanonical"
            for row in rows:
                row["branch_threshold_version"] = "v2-noncanonical"
                row["branch_ratio_evidence"]["threshold_version"] = "v2-noncanonical"

        mutations = (
            ("missing_sentaurus_evidence", delete_evidence(0)),
            ("missing_vela_evidence", delete_evidence(1)),
            ("ratio", alter_evidence("absolute_vela_over_sentaurus", 6.0)),
            ("vela_current", alter_evidence("vela_anode_current_A_per_um", 6.0)),
            ("sentaurus_current", alter_evidence("sentaurus_anode_current_A_per_um", 2.0)),
            ("geometric_zero", alter_evidence("geometric_zero", True)),
            ("ratio_bool", alter_evidence("absolute_vela_over_sentaurus", True)),
            ("sentaurus_current_bool", alter_evidence("sentaurus_anode_current_A_per_um", True)),
            ("geometric_zero_int", alter_evidence("geometric_zero", 0)),
            ("rows_differ", rows_differ),
            ("wrong_allowed_classification", wrong_classification),
            ("noncanonical_threshold", noncanonical_threshold),
        )
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            baseline, _, _ = self._complete_package_with_common_pair(root)
            validate_sweep_manifest(baseline, package_root=root)
        for label, mutate in mutations:
            with self.subTest(label=label), tempfile.TemporaryDirectory() as temp:
                root = Path(temp)
                manifest, sentaurus, vela = self._complete_package_with_common_pair(root)
                mutate(manifest, [sentaurus, vela])
                with self.assertRaises(ValueError):
                    validate_sweep_manifest(manifest, package_root=root)

    def test_strict_common_pair_uses_none_ratio_only_for_zero_reference(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            manifest, sentaurus, vela = self._complete_package_with_common_pair(
                root, sentaurus_current=0.0
            )
            validate_sweep_manifest(manifest, package_root=root)
            for row in (sentaurus, vela):
                self.assertEqual(row["branch_classification"], "unidentified")
                self.assertIsNone(
                    row["branch_ratio_evidence"]["absolute_vela_over_sentaurus"]
                )

    def test_strict_side_only_row_stays_unidentified_without_ratio_evidence(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            manifest = self._complete_package(root)
            state = segment_state_path(root, "sketch", 0.0, -1.0)
            state.parent.mkdir(parents=True, exist_ok=True)
            state.write_text("vela state\n", encoding="utf-8")
            row = record_transition(
                manifest,
                solver="vela",
                topology="sketch",
                start_bias_V=0.0,
                target_bias_V=-1.0,
                exit_code=0,
                actual_bias_V=-1.0,
                state_path=state,
                observables={
                    "anode_current_A_per_um": 5.0,
                    "cathode_current_A_per_um": -5.0,
                    "max_field_V_per_m": 2.0,
                    "native_source_integral_s_inv_per_cm": 3.0,
                    "reconstructed_source_integral_s_inv_per_cm": 4.0,
                },
            )
            row["state_path"] = str(state.relative_to(root))
            validate_sweep_manifest(manifest, package_root=root)
            self.assertEqual(row["branch_classification"], "unidentified")
            self.assertNotIn("branch_ratio_evidence", row)
            row["branch_ratio_evidence"] = {
                "vela_anode_current_A_per_um": 5.0,
                "sentaurus_anode_current_A_per_um": 1.0,
                "absolute_vela_over_sentaurus": 5.0,
                "geometric_zero": False,
                "threshold_version": manifest["branch_threshold_version"],
            }
            with self.assertRaises(ValueError):
                validate_sweep_manifest(manifest, package_root=root)

    def test_common_pair_refresh_uses_exact_bias_tolerance(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            manifest, sentaurus, vela = self._complete_package_with_common_pair(
                root, actual_offset_V=5.0e-13
            )
            validate_sweep_manifest(manifest, package_root=root)
            for row in (sentaurus, vela):
                self.assertEqual(row["branch_classification"], "multiplication_like")
                self.assertEqual(
                    row["branch_ratio_evidence"]["absolute_vela_over_sentaurus"],
                    5.0,
                )

    def test_complete_package_rejects_config_drift_and_identity_tampering(self):
        def template_tamper(root, manifest):
            path = Path(manifest["template"]["path"])
            path.write_text(path.read_text(encoding="utf-8") + "\n", encoding="utf-8")

        def input_tamper(name):
            def mutate(root, manifest):
                path = root / "inputs" / "sketch" / name
                path.write_text(path.read_text(encoding="utf-8") + "\n", encoding="utf-8")
            return mutate

        def missing_topology_hash(root, manifest):
            del manifest["topology_input_sha256"]["mirror"]

        def vela_topology_mismatch(root, manifest):
            row = manifest["segments"][0]
            path = root / row["deck"]
            deck = json.loads(path.read_text(encoding="utf-8"))
            deck["mesh_file"] = str(root / "inputs" / "mirror" / "mesh.json")
            path.write_text(json.dumps(deck, sort_keys=True, indent=2) + "\n", encoding="utf-8")
            row["deck_sha256"] = self._sha256(path)

        def duplicate_state_prefix(root, manifest):
            first_path = root / manifest["segments"][0]["deck"]
            row = manifest["segments"][1]
            path = root / row["deck"]
            first = json.loads(first_path.read_text(encoding="utf-8"))
            deck = json.loads(path.read_text(encoding="utf-8"))
            deck["sweep"]["write_state_every_point_prefix"] = first["sweep"]["write_state_every_point_prefix"]
            path.write_text(json.dumps(deck, sort_keys=True, indent=2) + "\n", encoding="utf-8")
            row["deck_sha256"] = self._sha256(path)

        def sentaurus_stale_hash(root, manifest):
            path = root / manifest["sentaurus_segments"][0]["deck"]
            path.write_text(path.read_text(encoding="utf-8") + "\n", encoding="utf-8")

        def sentaurus_semantic_tamper(root, manifest):
            row = manifest["sentaurus_segments"][0]
            path = root / row["deck"]
            payload = path.read_text(encoding="utf-8")
            path.write_text(payload.replace("Voltage=0.0", "Voltage=0.5", 1), encoding="utf-8")
            row["deck_sha256"] = self._sha256(path)

        def remove_vela(root, manifest):
            manifest["segments"].pop()

        def duplicate_vela(root, manifest):
            manifest["segments"].append(copy.deepcopy(manifest["segments"][0]))

        def noncanonical_vela(root, manifest):
            manifest["segments"][0]["target_bias_V"] = -1.5

        def remove_sentaurus(root, manifest):
            manifest["sentaurus_segments"].pop()

        def duplicate_sentaurus(root, manifest):
            manifest["sentaurus_segments"].append(copy.deepcopy(manifest["sentaurus_segments"][0]))

        def noncanonical_sentaurus(root, manifest):
            manifest["sentaurus_segments"][0]["target_bias_V"] = -0.5

        mutations = (
            ("template_hash", template_tamper),
            ("mesh_hash", input_tamper("mesh.json")),
            ("doping_hash", input_tamper("doping.csv")),
            ("missing_topology_hash", missing_topology_hash),
            ("vela_topology_input", vela_topology_mismatch),
            ("duplicate_state_prefix", duplicate_state_prefix),
            ("sentaurus_stale_hash", sentaurus_stale_hash),
            ("sentaurus_semantics", sentaurus_semantic_tamper),
            ("missing_vela_segment", remove_vela),
            ("duplicate_vela_segment", duplicate_vela),
            ("noncanonical_vela_segment", noncanonical_vela),
            ("missing_sentaurus_segment", remove_sentaurus),
            ("duplicate_sentaurus_segment", duplicate_sentaurus),
            ("noncanonical_sentaurus_segment", noncanonical_sentaurus),
        )
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            self._validate_package_compat(self._complete_package(root), root)
        for label, mutate in mutations:
            with self.subTest(label=label), tempfile.TemporaryDirectory() as temp:
                root = Path(temp)
                manifest = self._complete_package(root)
                mutate(root, manifest)
                with self.assertRaises(ValueError):
                    self._validate_package_compat(manifest, root)

    def test_package_records_and_rechecks_authoritative_state_root(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            manifest = self._complete_package(root)
            fixture = (
                Path(__file__).resolve().parents[2]
                / "tests"
                / "fixtures"
                / "pn2d_minimal6_synthetic"
            ).resolve()
            expected = {
                "resolved_path": str(fixture),
                "manifest_json_sha256": self._sha256(fixture / "manifest.json"),
                "topology_sources": {
                    topology: {
                        "mesh_path": f"states/{topology}/p0V/export/mesh.json",
                        "mesh_sha256": self._sha256(
                            fixture / "states" / topology / "p0V" / "export" / "mesh.json"
                        ),
                        "doping_path": f"states/{topology}/p0V/export/doping.csv",
                        "doping_sha256": self._sha256(
                            fixture / "states" / topology / "p0V" / "export" / "doping.csv"
                        ),
                    }
                    for topology in ("sketch", "mirror")
                },
            }
            self.assertEqual(manifest["authoritative_state_root"], expected)
            missing = copy.deepcopy(manifest)
            del missing["authoritative_state_root"]
            with self.assertRaisesRegex(ValueError, "authoritative state root"):
                validate_sweep_manifest(missing, package_root=root)
            tampered = copy.deepcopy(manifest)
            tampered["authoritative_state_root"]["topology_sources"]["sketch"][
                "mesh_sha256"
            ] = tampered["authoritative_state_root"]["topology_sources"]["mirror"][
                "mesh_sha256"
            ]
            with self.assertRaisesRegex(ValueError, "authoritative"):
                validate_sweep_manifest(tampered, package_root=root)

    def test_package_rejects_synchronized_cross_topology_input_substitution(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            manifest = self._complete_package(root)
            for name in ("mesh.json", "doping.csv"):
                source = root / "inputs" / "mirror" / name
                destination = root / "inputs" / "sketch" / name
                shutil.copyfile(source, destination)
                manifest["topology_input_sha256"]["sketch"][name] = self._sha256(
                    destination
                )
            with self.assertRaisesRegex(ValueError, "authoritative"):
                validate_sweep_manifest(manifest, package_root=root)

    def test_generated_package_rejects_deck_hash_tampering(self):
        signature = inspect.signature(validate_sweep_manifest)
        self.assertIn(
            "package_root",
            signature.parameters,
            "manifest validation must verify relative package artifacts on disk",
        )
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            manifest = self._complete_package(root)
            validate_sweep_manifest(manifest, package_root=root)
            deck = root / manifest["segments"][0]["deck"]
            payload = json.loads(deck.read_text(encoding="utf-8"))
            payload["solver"]["max_iter"] += 1
            deck.write_text(json.dumps(payload), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "hash-tampered"):
                validate_sweep_manifest(manifest, package_root=root)

    def test_sentaurus_endpoint_recovers_complete_observables_from_export(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            export = root / "export"
            fields = export / "fields"
            fields.mkdir(parents=True)
            mesh = root / "mesh.json"
            mesh.write_text('{"nodes":[{"id":0,"x":0,"y":0},{"id":1,"x":1e-6,"y":0},{"id":2,"x":0,"y":1e-6}],"triangles":[{"node_ids":[0,1,2]}]}', encoding="utf-8")
            (export / "field_manifest.json").write_text('{"fields":[{"name":"ElectricField","components":2,"unit":"V*cm^-1","region":0,"region_name":"R.Si"},{"name":"ImpactIonization","components":1,"unit":"cm^-3*s^-1","region":0,"region_name":"R.Si"},{"name":"eAlphaAvalanche","components":1,"unit":"cm^-1","region":0,"region_name":"R.Si"},{"name":"hAlphaAvalanche","components":1,"unit":"cm^-1","region":0,"region_name":"R.Si"},{"name":"eCurrentDensity","components":2,"unit":"A*cm^-2","region":0,"region_name":"R.Si"},{"name":"hCurrentDensity","components":2,"unit":"A*cm^-2","region":0,"region_name":"R.Si"},{"name":"ContactCurrentFlux","components":1,"unit":"A","region":1,"region_name":"Cathode"},{"name":"ContactCurrentFlux","components":1,"unit":"A","region":2,"region_name":"Anode"},{"name":"ContactExternalVoltage","components":1,"unit":"V","region":1,"region_name":"Cathode"},{"name":"ContactExternalVoltage","components":1,"unit":"V","region":2,"region_name":"Anode"}]}', encoding="utf-8")
            def scalar(name, values):
                (fields / f"{name}_region0.csv").write_text("node_id,component0\n" + "\n".join(f"{node},{value}" for node, value in values.items()) + "\n", encoding="utf-8")
            def vector(name, values):
                (fields / f"{name}_region0.csv").write_text("node_id,component0,component1\n" + "\n".join(f"{node},{x},{y}" for node, (x, y) in values.items()) + "\n", encoding="utf-8")
            scalar("ImpactIonization", {0: 2.0, 1: 2.0, 2: 2.0})
            scalar("eAlphaAvalanche", {0: 1.0, 1: 1.0, 2: 1.0})
            scalar("hAlphaAvalanche", {0: 0.0, 1: 0.0, 2: 0.0})
            vector("ElectricField", {0: (3.0, 4.0), 1: (3.0, 4.0), 2: (3.0, 4.0)})
            vector("eCurrentDensity", {0: (1.602176634e-19, 0.0), 1: (1.602176634e-19, 0.0), 2: (1.602176634e-19, 0.0)})
            vector("hCurrentDensity", {0: (0.0, 0.0), 1: (0.0, 0.0), 2: (0.0, 0.0)})
            (fields / "ContactCurrentFlux_region1.csv").write_text("node_id,component0\n2,2\n", encoding="utf-8")
            (fields / "ContactCurrentFlux_region2.csv").write_text("node_id,component0\n0,-2\n", encoding="utf-8")
            (fields / "ContactExternalVoltage_region1.csv").write_text("node_id,component0\n2,0\n", encoding="utf-8")
            (fields / "ContactExternalVoltage_region2.csv").write_text("node_id,component0\n0,-1\n", encoding="utf-8")
            endpoint = read_sentaurus_endpoint(export, mesh, -1.0)
            self.assertEqual(endpoint["actual_bias_V"], -1.0)
            self.assertEqual(endpoint["observables"]["anode_current_A_per_um"], -2.0)
            self.assertEqual(endpoint["observables"]["cathode_current_A_per_um"], 2.0)
            self.assertEqual(endpoint["observables"]["max_field_V_per_m"], 500.0)
            self.assertAlmostEqual(endpoint["observables"]["native_source_integral_s_inv_per_cm"], 1.0e-8)
            self.assertAlmostEqual(endpoint["observables"]["reconstructed_source_integral_s_inv_per_cm"], 5.0e-9)

    def test_required_field_region_selects_exactly_one_contract_match(self):
        fields = [
            {"name": "ElectricField", "components": 1, "unit": "V*cm^-1", "region": 0, "region_name": "R.Si"},
            {"name": "ElectricField", "components": 2, "unit": "V*cm^-1", "region": 0, "region_name": "R.Si"},
        ]
        self.assertEqual(
            _required_field_region(fields, "ElectricField", components=2, unit="V*cm^-1"),
            0,
        )
        fields.append(
            {"name": "ElectricField", "components": 2, "unit": "V*cm^-1", "region": 0, "region_name": "R.Si"}
        )
        with self.assertRaisesRegex(ValueError, "exactly one contract-compatible"):
            _required_field_region(fields, "ElectricField", components=2, unit="V*cm^-1")
    def test_cli_resume_loads_declared_sentaurus_checkpoint_results(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            (root / "sweep_manifest.json").write_text('{"schema":"vela.pn2d_minimal6_sweep_manifest.v1","targets_V":[0.0,-1.0],"segments":[],"sentaurus_segments":[],"accepted_checkpoints":[],"failed_transition":null,"interpolation":"forbidden"}\n', encoding="utf-8")
            results = root / "results.json"
            results.write_text('[{"topology":"sketch","start_bias_V":0.0,"target_bias_V":-1.0,"state_path":"checkpoint.tdr","export_dir":"export","mesh_path":"mesh.json"}]\n', encoding="utf-8-sig")
            with patch("scripts.run_pn2d_minimal6_diagnostic_sweep.record_sentaurus_checkpoint") as record:
                with patch("scripts.run_pn2d_minimal6_diagnostic_sweep.validate_sweep_manifest") as validate:
                    with patch.object(sys, "argv", ["diagnostic-sweep", "--out-dir", str(root), "--resume", "--sentaurus-results-json", str(results)]):
                        self.assertEqual(main(), 0)
            self.assertEqual(record.call_count, 1)
            self.assertEqual(record.call_args.kwargs["topology"], "sketch")
            validate.assert_called_once()
            self.assertEqual(validate.call_args.args[0]["schema"], "vela.pn2d_minimal6_sweep_manifest.v1")
            self.assertEqual(validate.call_args.kwargs["package_root"], root.resolve())
    def test_cli_resolves_relative_output_root_before_running_subprocess(self):
        with tempfile.TemporaryDirectory() as temp:
            with patch("scripts.run_pn2d_minimal6_diagnostic_sweep.initialise_package", return_value=Path(temp) / "sweep_manifest.json") as initialise:
                with patch.object(sys, "argv", ["diagnostic-sweep", "--out-dir", "relative-evidence"]):
                    self.assertEqual(main(), 0)
            self.assertTrue(initialise.call_args.args[0].is_absolute())
    def test_integer_targets_are_exact_and_complete(self):
        self.assertEqual(integer_targets(), tuple(float(-value) for value in range(21)))

    def test_branch_classification_refuses_zero_and_uses_declared_thresholds(self):
        self.assertIn("geometric_zero", inspect.signature(classify_branch).parameters)
        self.assertEqual(classify_branch(1.0, 5.0), "multiplication_like")
        self.assertEqual(classify_branch(1.0, 1.0e-4), "leakage_like")
        self.assertEqual(classify_branch(1.0, 1.0e-2), "unidentified")
        self.assertEqual(classify_branch(0.0, 1.0), "unidentified")
        self.assertEqual(
            classify_branch(1.0, 1.0, geometric_zero=True), "unidentified"
        )

    def test_accepted_common_exact_rows_refresh_ratio_branch_evidence(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            manifest = {
                "schema": "vela.pn2d_minimal6_sweep_manifest.v1",
                "targets_V": [0.0, -1.0],
                "accepted_checkpoints": [],
                "failed_transition": None,
                "failed_transitions": [],
                "interpolation": "forbidden",
                "branch_threshold_version": "v1: multiplication=[0.1,10], leakage<=1e-3",
            }
            sentaurus_state = root / "sentaurus.tdr"
            vela_state = root / "vela.csv"
            sentaurus_state.write_text("state\n", encoding="utf-8")
            vela_state.write_text("state\n", encoding="utf-8")
            sentaurus_observables = {
                "anode_current_A_per_um": 1.0,
                "cathode_current_A_per_um": -1.0,
                "max_field_V_per_m": 2.0,
                "native_source_integral_s_inv_per_cm": 3.0,
                "reconstructed_source_integral_s_inv_per_cm": 4.0,
            }
            vela_observables = dict(sentaurus_observables)
            vela_observables["anode_current_A_per_um"] = 5.0
            record_transition(
                manifest,
                solver="sentaurus",
                topology="sketch",
                start_bias_V=0.0,
                target_bias_V=-1.0,
                exit_code=0,
                actual_bias_V=-1.0,
                state_path=sentaurus_state,
                observables=sentaurus_observables,
            )
            record_transition(
                manifest,
                solver="vela",
                topology="sketch",
                start_bias_V=0.0,
                target_bias_V=-1.0,
                exit_code=0,
                actual_bias_V=-1.0,
                state_path=vela_state,
                observables=vela_observables,
            )
            sentaurus, vela = manifest["accepted_checkpoints"]
            expected_evidence = {
                "vela_anode_current_A_per_um": 5.0,
                "sentaurus_anode_current_A_per_um": 1.0,
                "absolute_vela_over_sentaurus": 5.0,
                "geometric_zero": False,
                "threshold_version": manifest["branch_threshold_version"],
            }
            for row in (sentaurus, vela):
                self.assertEqual(row["branch_classification"], "multiplication_like")
                self.assertEqual(row["branch_ratio_evidence"], expected_evidence)

    def test_branch_refresh_is_exact_only_and_geometric_zero_fails_closed(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            manifest = {
                "schema": "vela.pn2d_minimal6_sweep_manifest.v1",
                "targets_V": [0.0, -1.0, -2.0],
                "accepted_checkpoints": [],
                "failed_transition": None,
                "failed_transitions": [],
                "interpolation": "forbidden",
                "branch_threshold_version": "v1: multiplication=[0.1,10], leakage<=1e-3",
            }
            def accepted(solver, bias, current, native, reconstructed):
                state = root / f"{solver}_{abs(int(bias))}.state"
                state.write_text("state\n", encoding="utf-8")
                return record_transition(
                    manifest,
                    solver=solver,
                    topology="mirror",
                    start_bias_V=bias + 1.0,
                    target_bias_V=bias,
                    exit_code=0,
                    actual_bias_V=bias,
                    state_path=state,
                    observables={
                        "anode_current_A_per_um": current,
                        "cathode_current_A_per_um": -current,
                        "max_field_V_per_m": 2.0,
                        "native_source_integral_s_inv_per_cm": native,
                        "reconstructed_source_integral_s_inv_per_cm": reconstructed,
                    },
                )
            sentaurus = accepted("sentaurus", -1.0, 1.0, 0.0, 1.0e-286)
            vela_other_bias = accepted("vela", -2.0, 5.0, 3.0, 4.0)
            self.assertEqual(sentaurus["branch_classification"], "unidentified")
            self.assertNotIn("branch_ratio_evidence", sentaurus)
            self.assertEqual(vela_other_bias["branch_classification"], "unidentified")
            self.assertNotIn("branch_ratio_evidence", vela_other_bias)
            vela = accepted("vela", -1.0, 5.0, 3.0, 4.0)
            for row in (sentaurus, vela):
                self.assertEqual(row["branch_classification"], "unidentified")
                self.assertTrue(row["branch_ratio_evidence"]["geometric_zero"])
                self.assertEqual(
                    row["branch_ratio_evidence"]["absolute_vela_over_sentaurus"],
                    5.0,
                )

    def test_input_copy_is_separate_from_authoritative_state_root(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            source = root / "authoritative" / "states" / "sketch" / "p0V" / "export"
            source.mkdir(parents=True)
            source_mesh = {
                "nodes": [
                    {"id": 0, "x": 0.0, "y": 0.5e-6},
                    {"id": 1, "x": 1.0e-6, "y": 0.5e-6},
                ],
                "triangles": [],
                "regions": [],
                "contacts": [],
            }
            (source / "mesh.json").write_text(
                json.dumps(source_mesh), encoding="utf-8"
            )
            (source / "doping.csv").write_text("doping", encoding="utf-8")
            destination = root / "sweep"
            hashes = copy_topology_inputs(root / "authoritative", destination, ("sketch",))
            copied_mesh = json.loads(
                (destination / "inputs" / "sketch" / "mesh.json").read_text(
                    encoding="utf-8"
                )
            )
            self.assertEqual(copied_mesh["nodes"][1]["x"], 1.0)
            self.assertEqual(copied_mesh["nodes"][0]["y"], 0.5)
            self.assertEqual(
                json.loads((source / "mesh.json").read_text(encoding="utf-8")),
                source_mesh,
            )
            self.assertEqual(
                (destination / "inputs" / "sketch" / "doping.csv").read_text(
                    encoding="utf-8"
                ),
                "doping",
            )
            self.assertIn("sketch", hashes)
    def test_package_local_sentaurus_rejects_checkpoint_path_escape_before_write(self):
        with tempfile.TemporaryDirectory() as temp:
            base = Path(temp)
            package = base / "package"
            package.mkdir()
            external = base / "external"
            external.mkdir()
            (external / "checkpoint.tdr").write_text("tdr\n", encoding="utf-8")
            (external / "export").mkdir()
            (external / "mesh.json").write_text("{}\n", encoding="utf-8")
            results = external / "results.json"
            results.write_text(
                json.dumps([
                    {
                        "topology": "sketch",
                        "start_bias_V": 0.0,
                        "target_bias_V": -1.0,
                        "state_path": "checkpoint.tdr",
                        "export_dir": "export",
                        "mesh_path": "mesh.json",
                    }
                ]),
                encoding="utf-8",
            )
            manifest = {
                "sentaurus_segments": [
                    {
                        "solver": "sentaurus",
                        "topology": "sketch",
                        "target_bias_V": -1.0,
                        "checkpoint_tdr": "../escape.tdr",
                        "status": "pending",
                    }
                ],
                "accepted_checkpoints": [],
                "failed_transition": None,
                "failed_transitions": [],
            }
            with self.assertRaisesRegex(ValueError, "canonical"):
                record_sentaurus_results(
                    manifest, results, package_root=package
                )
            self.assertFalse((base / "escape.tdr").exists())

    def test_package_local_sentaurus_missing_external_state_records_rejection(self):
        with tempfile.TemporaryDirectory() as temp:
            base = Path(temp)
            package = base / "package"
            package.mkdir()
            external = base / "external"
            external.mkdir()
            export = external / "export"
            export.mkdir()
            (export / "field_manifest.json").write_text("{}\n", encoding="utf-8")
            (external / "mesh.json").write_text("{}\n", encoding="utf-8")
            results = external / "results.json"
            results.write_text(
                json.dumps([
                    {
                        "topology": "sketch",
                        "start_bias_V": 0.0,
                        "target_bias_V": -1.0,
                        "state_path": "missing.tdr",
                        "export_dir": "export",
                        "mesh_path": "mesh.json",
                    }
                ]),
                encoding="utf-8",
            )
            checkpoint = (
                Path("sentaurus")
                / "sketch"
                / "checkpoints"
                / "sketch_m1p000000.tdr"
            )
            manifest = {
                "sentaurus_segments": [
                    {
                        "solver": "sentaurus",
                        "topology": "sketch",
                        "target_bias_V": -1.0,
                        "checkpoint_tdr": str(checkpoint),
                        "status": "pending",
                    }
                ],
                "accepted_checkpoints": [],
                "failed_transition": None,
                "failed_transitions": [],
            }
            rows = record_sentaurus_results(
                manifest, results, package_root=package
            )
            self.assertEqual(rows[0]["status"], "rejected")
            self.assertEqual(manifest["failed_transition"], rows[0])
            self.assertIn("missing", rows[0]["incomplete_reason"].lower())
            self.assertFalse((package / checkpoint).exists())

    def test_package_local_state_paths_are_canonical_and_external_import_is_preserved(self):
        with tempfile.TemporaryDirectory() as temp:
            base = Path(temp)
            package = base / "package"
            package.mkdir()
            external = base / "external"
            external.mkdir()
            state = external / "checkpoint.tdr"
            state.write_text("external tdr\n", encoding="utf-8")
            export = external / "export"
            export.mkdir()
            (export / "field_manifest.json").write_text("{}\n", encoding="utf-8")
            mesh = external / "mesh.json"
            mesh.write_text("{}\n", encoding="utf-8")
            results = external / "results.json"
            results.write_text(
                json.dumps([
                    {
                        "topology": "sketch",
                        "start_bias_V": 0.0,
                        "target_bias_V": -1.0,
                        "state_path": "checkpoint.tdr",
                        "export_dir": "export",
                        "mesh_path": "mesh.json",
                    }
                ]),
                encoding="utf-8",
            )
            checkpoint = (
                Path("sentaurus")
                / "sketch"
                / "checkpoints"
                / "sketch_m1p000000.tdr"
            )
            manifest = {
                "schema": "vela.pn2d_minimal6_sweep_manifest.v1",
                "targets_V": [0.0, -1.0],
                "segments": [],
                "sentaurus_segments": [
                    {
                        "solver": "sentaurus",
                        "topology": "sketch",
                        "target_bias_V": -1.0,
                        "checkpoint_tdr": str(checkpoint),
                        "status": "pending",
                    }
                ],
                "accepted_checkpoints": [],
                "failed_transition": None,
            }
            endpoint = {
                "actual_bias_V": -1.0,
                "depth_convention": "unit_out_of_plane_length_cm",
                "current_conversion": "fixture",
                "observables": {
                    "anode_current_A_per_um": -1.0,
                    "cathode_current_A_per_um": 1.0,
                    "max_field_V_per_m": 2.0,
                    "native_source_integral_s_inv_per_cm": 3.0,
                    "reconstructed_source_integral_s_inv_per_cm": 4.0,
                },
            }
            with patch(
                "scripts.run_pn2d_minimal6_diagnostic_sweep.read_sentaurus_endpoint",
                return_value=endpoint,
            ):
                rows = record_sentaurus_results(
                    manifest, results, package_root=package
                )
            copied = package / checkpoint
            self.assertEqual(copied.read_text(encoding="utf-8"), "external tdr\n")
            self.assertTrue(state.is_file())
            self.assertEqual(rows[0]["state_path"], str(checkpoint))
            self.assertEqual(rows[0]["state_sha256"], self._sha256(copied))

    def test_execute_segments_serializes_package_relative_vela_state_path(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            state = segment_state_path(root, "sketch", 0.0, -1.0)
            state.parent.mkdir(parents=True)
            state.write_text("state\n", encoding="utf-8")
            manifest = {
                "schema": "vela.pn2d_minimal6_sweep_manifest.v1",
                "targets_V": [0.0, -1.0],
                "segments": [
                    {
                        "solver": "vela",
                        "topology": "sketch",
                        "start_bias_V": 0.0,
                        "target_bias_V": -1.0,
                        "status": "pending",
                    }
                ],
                "accepted_checkpoints": [],
                "failed_transition": None,
            }
            execute_segments(
                manifest,
                root,
                lambda _: {
                    "exit_code": 0,
                    "actual_bias_V": -1.0,
                    "state_path": state,
                    "observables": {
                        "anode_current_A_per_um": 1.0,
                        "cathode_current_A_per_um": -1.0,
                        "max_field_V_per_m": 2.0,
                        "native_source_integral_s_inv_per_cm": 3.0,
                        "reconstructed_source_integral_s_inv_per_cm": 4.0,
                    },
                },
                package_local_state_paths=True,
            )
            self.assertEqual(
                manifest["accepted_checkpoints"][0]["state_path"],
                str(state.relative_to(root)),
            )

    def test_strict_package_validator_resolves_canonical_relative_vela_state(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            manifest = self._complete_package(root)
            state = segment_state_path(root, "sketch", 0.0, -1.0)
            state.parent.mkdir(parents=True, exist_ok=True)
            state.write_text("state\n", encoding="utf-8")
            row = record_transition(
                manifest,
                solver="vela",
                topology="sketch",
                start_bias_V=0.0,
                target_bias_V=-1.0,
                exit_code=0,
                actual_bias_V=-1.0,
                state_path=state,
                observables={
                    "anode_current_A_per_um": 1.0,
                    "cathode_current_A_per_um": -1.0,
                    "max_field_V_per_m": 2.0,
                    "native_source_integral_s_inv_per_cm": 3.0,
                    "reconstructed_source_integral_s_inv_per_cm": 4.0,
                },
            )
            row["state_path"] = str(state.relative_to(root))
            validate_sweep_manifest(manifest, package_root=root)

    def test_record_sentaurus_checkpoint_promotes_complete_export_to_manifest(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            state = root / "checkpoint.tdr"
            state.write_text("tdr\n", encoding="utf-8")
            export = root / "export"
            export.mkdir()
            (export / "field_manifest.json").write_text("{}\n", encoding="utf-8")
            manifest = {"schema": "vela.pn2d_minimal6_sweep_manifest.v1", "targets_V": [0.0, -1.0], "segments": [], "sentaurus_segments": [{"solver": "sentaurus", "topology": "sketch", "target_bias_V": -1.0, "status": "pending"}], "accepted_checkpoints": [], "failed_transition": None}
            endpoint = {"actual_bias_V": -1.0, "depth_convention": "unit_out_of_plane_length_cm", "current_conversion": "Sentaurus 2-D ContactCurrentFlux A compared numerically with Vela A/um", "observables": {"anode_current_A_per_um": -1.0, "cathode_current_A_per_um": 1.0, "max_field_V_per_m": 2.0, "native_source_integral_s_inv_per_cm": 3.0, "reconstructed_source_integral_s_inv_per_cm": 4.0}}
            with patch("scripts.run_pn2d_minimal6_diagnostic_sweep.read_sentaurus_endpoint", return_value=endpoint):
                row = record_sentaurus_checkpoint(manifest, topology="sketch", start_bias_V=0.0, target_bias_V=-1.0, state_path=state, export_dir=export, mesh_path=root / "mesh.json")
            self.assertEqual(row["status"], "accepted")
            self.assertEqual(manifest["sentaurus_segments"][0]["status"], "accepted")
            self.assertEqual(manifest["accepted_checkpoints"], [row])
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
    def test_sentaurus_sweep_source_requests_vector_current_densities(self):
        source = Path(__file__).resolve().parents[2] / "reference_tcad" / "pn2d_sentaurus2018_minimal6" / "source" / "pn2d_minimal6_sweep_sdevice.cmd"
        deck = source.read_text(encoding="utf-8")
        self.assertIn("eCurrentDensity/Vector", deck)
        self.assertIn("hCurrentDensity/Vector", deck)

    def test_vela_sweep_template_uses_sentaurus_equivalent_ohmic_qf_pinning(self):
        source = (
            Path(__file__).resolve().parents[2]
            / "reference_tcad"
            / "pn2d_sentaurus2018_minimal6"
            / "vela"
            / "pn2d_minimal6_sweep_template.json"
        )
        deck = json.loads(source.read_text(encoding="utf-8"))
        self.assertIs(
            deck["solver"]["contact_boundary_minority_electron_relaxation"],
            False,
        )

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
            edges = root / "edges.csv"
            curve.write_text("bias_V,max_electric_field_V_per_m\n-1,2\n", encoding="utf-8")
            terminals.write_text("bias_V,contact,I_sgflux_A_per_um,sg_avalanche_source_integral_total\n-1,Anode,-3,5\n-1,Cathode,3,5\n", encoding="utf-8")
            edges.write_text(
                "bias_V,edge_id,edge_source_integral\n"
                "-1,0,1.25\n"
                "-1,1,3.75\n",
                encoding="utf-8",
            )
            result = read_vela_endpoint(curve, terminals, edges, -1.0)
            self.assertEqual(result["anode_current_A_per_um"], -3.0)
            self.assertAlmostEqual(result["native_source_integral_s_inv_per_cm"], 5.0e-8)
            self.assertAlmostEqual(result["reconstructed_source_integral_s_inv_per_cm"], 5.0e-8)
            with self.assertRaises(ValueError): read_vela_endpoint(curve, terminals, edges, -2.0)
            edges.write_text("bias_V,edge_id,edge_source_integral\n-1,0,4.0\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "source closure"):
                read_vela_endpoint(curve, terminals, edges, -1.0)
    def _prepare_vela_segment_artifacts(self, root, topology):
        state = segment_state_path(root, topology, 0.0, -1.0)
        state.parent.mkdir(parents=True, exist_ok=True)
        state.write_text("state\n", encoding="utf-8")
        diagnostics = root / "vela" / topology / "diagnostics"
        diagnostics.mkdir(parents=True, exist_ok=True)
        (diagnostics / "segment_00_terminal_current_method_compare.csv").write_text(
            "placeholder\n", encoding="utf-8"
        )
        (diagnostics / "segment_00_sg_avalanche_edges.csv").write_text(
            "placeholder\n", encoding="utf-8"
        )
        curve = root / "vela" / topology / "segment_00.csv"
        curve.parent.mkdir(parents=True, exist_ok=True)
        curve.write_text("placeholder\n", encoding="utf-8")
        return {
            "solver": "vela",
            "topology": topology,
            "start_bias_V": 0.0,
            "target_bias_V": -1.0,
            "deck": f"vela/{topology}/decks/segment_00.json",
            "status": "pending",
        }

    def test_vela_curve_parser_serializes_exact_convergence_and_failure_fields(self):
        with tempfile.TemporaryDirectory() as temp:
            curve = Path(temp) / "curve.csv"
            curve.write_text(
                "bias_V,max_electric_field_V_per_m,converged,iterations,newton_iterations,step_diagnostics,validation_diagnostics,failure_reason,newton_failure_class,newton_failure_diagnostics_json\n"
                "0,1,1,1,1,accepted_step=0,validation=ok,,,\n"
                "-1,2,0,7,4,attempted_step=-1;retry_count=2,carrier=invalid,minimum step reached,line_search_non_decrease,failure.json\n",
                encoding="utf-8",
            )
            self.assertEqual(
                read_vela_convergence_metadata(curve, -1.0),
                {
                    "converged": False,
                    "iterations": 7,
                    "newton_iterations": 4,
                    "step_diagnostics": "attempted_step=-1;retry_count=2",
                    "validation_diagnostics": "carrier=invalid",
                    "failure_reason": "minimum step reached",
                    "newton_failure_class": "line_search_non_decrease",
                    "newton_failure_diagnostics_json": "failure.json",
                },
            )
            inexact = Path(temp) / "inexact.csv"
            inexact.write_text(
                "bias_V,converged,iterations,newton_iterations,step_diagnostics\n"
                "-0.9,1,2,2,accepted_step=-0.9\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "exactly one target-bias"):
                read_vela_convergence_metadata(inexact, -1.0)
            nonfinite = Path(temp) / "nonfinite.csv"
            nonfinite.write_text(
                "bias_V,converged,iterations,newton_iterations,step_diagnostics\n"
                "-1,1,nan,2,accepted_step=-1\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(ValueError, "finite integer"):
                read_vela_convergence_metadata(nonfinite, -1.0)

    def test_real_vela_adapter_returns_curve_convergence_on_success_and_failure(self):
        observables = {
            "anode_current_A_per_um": 1.0,
            "cathode_current_A_per_um": -1.0,
            "max_field_V_per_m": 2.0,
            "native_source_integral_s_inv_per_cm": 3.0,
            "reconstructed_source_integral_s_inv_per_cm": 4.0,
        }
        success_row = "-1,2,1,4,4,accepted_step=-1;retry_count=0,validation=ok,\n"
        failure_row = "-1,2,0,7,4,attempted_step=-1;retry_count=2,carrier=invalid,minimum step reached\n"
        success_metadata = {
            "converged": True,
            "iterations": 4,
            "newton_iterations": 4,
            "step_diagnostics": "accepted_step=-1;retry_count=0",
            "validation_diagnostics": "validation=ok",
        }
        failure_metadata = {
            "converged": False,
            "iterations": 7,
            "newton_iterations": 4,
            "step_diagnostics": "attempted_step=-1;retry_count=2",
            "validation_diagnostics": "carrier=invalid",
            "failure_reason": "minimum step reached",
        }
        cases = (
            ("success", 0, observables, success_row, success_metadata, True),
            ("endpoint_failure", 0, ValueError("source closure failed"), failure_row, failure_metadata, False),
            ("solver_failure", 9, None, failure_row, failure_metadata, False),
            ("nonconverged_exit_zero", 0, observables, failure_row, failure_metadata, False),
        )
        for label, returncode, endpoint_side_effect, curve_row, expected_metadata, accepted in cases:
            with self.subTest(label=label), tempfile.TemporaryDirectory() as temp:
                root = Path(temp)
                segment = self._prepare_vela_segment_artifacts(root, "sketch")
                curve = root / "vela" / "sketch" / "segment_00.csv"
                curve.write_text(
                    "bias_V,max_electric_field_V_per_m,converged,iterations,newton_iterations,step_diagnostics,validation_diagnostics,failure_reason\n"
                    + curve_row,
                    encoding="utf-8",
                )
                completed = SimpleNamespace(
                    returncode=returncode,
                    stdout="solver stdout",
                    stderr="solver stderr",
                )
                with patch(
                    "scripts.run_pn2d_minimal6_diagnostic_sweep.subprocess.run",
                    return_value=completed,
                ):
                    with patch(
                        "scripts.run_pn2d_minimal6_diagnostic_sweep.read_vela_endpoint",
                        side_effect=endpoint_side_effect,
                    ):
                        result = run_vela_subprocess_segment(
                            root, Path("vela.exe"), segment
                        )
                self.assertEqual(result["convergence_metadata"], expected_metadata)
                if accepted:
                    self.assertIsNotNone(result["observables"])
                else:
                    self.assertIsNone(result["observables"])
                if label == "nonconverged_exit_zero":
                    self.assertIn("not converged", result["incomplete_reason"])

    def test_real_vela_adapter_converts_endpoint_errors_to_rejected_payload(self):
        errors = (
            "Vela curve lacks exactly one target-bias endpoint",
            "could not convert string to float: malformed",
            "Vela native/reconstructed source closure failed",
        )
        for message in errors:
            with self.subTest(message=message), tempfile.TemporaryDirectory() as temp:
                root = Path(temp)
                segment = self._prepare_vela_segment_artifacts(root, "sketch")
                completed = SimpleNamespace(
                    returncode=0, stdout="solver stdout", stderr="solver stderr"
                )
                with patch(
                    "scripts.run_pn2d_minimal6_diagnostic_sweep.subprocess.run",
                    return_value=completed,
                ):
                    with patch(
                        "scripts.run_pn2d_minimal6_diagnostic_sweep.read_vela_endpoint",
                        side_effect=ValueError(message),
                    ):
                        result = run_vela_subprocess_segment(
                            root, Path("vela.exe"), segment
                        )
                self.assertEqual(result["exit_code"], 0)
                self.assertIsNone(result["actual_bias_V"])
                self.assertIsNone(result["state_path"])
                self.assertIsNone(result["observables"])
                self.assertEqual(result["stdout"], "solver stdout")
                self.assertEqual(result["stderr"], "solver stderr")
                self.assertIn(message, result["incomplete_reason"])

    def test_real_vela_adapter_failure_is_retained_and_other_topology_runs(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            sketch = self._prepare_vela_segment_artifacts(root, "sketch")
            mirror = self._prepare_vela_segment_artifacts(root, "mirror")
            manifest = {
                "schema": "vela.pn2d_minimal6_sweep_manifest.v1",
                "targets_V": [0.0, -1.0],
                "segments": [sketch, mirror],
                "accepted_checkpoints": [],
                "failed_transition": None,
                "failed_transitions": [],
            }
            completed = SimpleNamespace(
                returncode=0, stdout="solver stdout", stderr="solver stderr"
            )
            observables = {
                "anode_current_A_per_um": 1.0,
                "cathode_current_A_per_um": -1.0,
                "max_field_V_per_m": 2.0,
                "native_source_integral_s_inv_per_cm": 3.0,
                "reconstructed_source_integral_s_inv_per_cm": 4.0,
            }
            with patch(
                "scripts.run_pn2d_minimal6_diagnostic_sweep.subprocess.run",
                return_value=completed,
            ) as subprocess_run:
                with patch(
                    "scripts.run_pn2d_minimal6_diagnostic_sweep.read_vela_endpoint",
                    side_effect=[
                        ValueError("Vela curve lacks exactly one target-bias endpoint"),
                        observables,
                    ],
                ):
                    execute_segments(
                        manifest,
                        root,
                        lambda segment: run_vela_subprocess_segment(
                            root, Path("vela.exe"), segment
                        ),
                    )
            self.assertEqual(subprocess_run.call_count, 2)
            self.assertEqual(manifest["failed_transition"]["topology"], "sketch")
            self.assertIn(
                "exactly one target-bias endpoint",
                manifest["failed_transition"]["incomplete_reason"],
            )
            mirror_rows = [
                row
                for row in manifest["accepted_checkpoints"]
                if row["topology"] == "mirror"
            ]
            self.assertEqual(len(mirror_rows), 1)

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

    def test_fake_runner_preserves_accepted_and_rejected_convergence_metadata(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            manifest = {
                "schema": "vela.pn2d_minimal6_sweep_manifest.v1",
                "targets_V": [0.0, -1.0, -2.0],
                "segments": [
                    {"solver": "vela", "topology": "sketch", "start_bias_V": 0.0, "target_bias_V": -1.0, "status": "pending"},
                    {"solver": "vela", "topology": "sketch", "start_bias_V": -1.0, "target_bias_V": -2.0, "status": "pending"},
                ],
                "accepted_checkpoints": [],
                "failed_transition": None,
            }
            accepted_metadata = {
                "newton_iterations": 4,
                "continuation_retries": 1,
                "residual_norm": 1.0e-10,
            }
            rejected_metadata = {
                "newton_iterations": 40,
                "continuation_retries": 29,
                "failure_reason": "minimum step reached",
            }

            def runner(segment):
                if segment["target_bias_V"] == -1.0:
                    state = root / "state_m1.csv"
                    state.write_text("state\n", encoding="utf-8")
                    return {
                        "exit_code": 0,
                        "actual_bias_V": -1.0,
                        "state_path": state,
                        "observables": {
                            "anode_current_A_per_um": 1.0,
                            "cathode_current_A_per_um": -1.0,
                            "max_field_V_per_m": 2.0,
                            "native_source_integral_s_inv_per_cm": 3.0,
                            "reconstructed_source_integral_s_inv_per_cm": 4.0,
                        },
                        "convergence_metadata": accepted_metadata,
                    }
                return {
                    "exit_code": 7,
                    "actual_bias_V": None,
                    "state_path": None,
                    "observables": None,
                    "convergence_metadata": rejected_metadata,
                }

            execute_segments(manifest, root, runner)
            self.assertEqual(
                manifest["accepted_checkpoints"][0]["convergence_metadata"],
                accepted_metadata,
            )
            self.assertEqual(
                manifest["failed_transition"]["convergence_metadata"],
                rejected_metadata,
            )

    def test_failed_topology_does_not_skip_other_topology_first_segment(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            manifest = {"schema": "vela.pn2d_minimal6_sweep_manifest.v1", "targets_V": [0.0, -1.0], "segments": [
                {"solver": "vela", "topology": "sketch", "start_bias_V": 0.0, "target_bias_V": -1.0, "status": "pending"},
                {"solver": "vela", "topology": "mirror", "start_bias_V": 0.0, "target_bias_V": -1.0, "status": "pending"}], "accepted_checkpoints": [], "failed_transition": None}
            calls = []
            def runner(segment):
                calls.append(segment["topology"])
                return {"exit_code": 0, "actual_bias_V": -1.0, "state_path": None, "observables": None, "incomplete_reason": "native source unavailable"}
            execute_segments(manifest, root, runner)
            self.assertEqual(calls, ["sketch", "mirror"])
            self.assertEqual(len(manifest["failed_transitions"]), 2)
            self.assertTrue(all(segment["status"] == "rejected" for segment in manifest["segments"]))
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
