#!/usr/bin/env python3
"""Regression tests for the compensated-junction SG replay orchestrator."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from unittest import mock
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]


def load_module(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


compare = load_module(
    REPO / "scripts" / "run_pn2d_coarse7x3_previous_full20_compare.py",
    "run_pn2d_coarse7x3_previous_full20_compare_orchestrator_test",
)
orchestrator = load_module(
    REPO / "scripts" / "reproduce_pn2d_bv_compensated_sg_replay.py",
    "reproduce_pn2d_bv_compensated_sg_replay_test",
)

def manifest_field(name: str, components: int, unit: str) -> dict[str, object]:
    return {
        "name": name,
        "components": components,
        "unit": unit,
        "region": 0,
        "mapping_status": "complete",
        "global_node_mapping": "global_vertex_order",
    }



class CompensatedSgReplayOrchestratorTest(unittest.TestCase):
    EXPECTED_IMPACT_BASE = {
        "model": "van_overstraeten",
        "driving_force": "quasi_fermi_gradient",
        "generation": "current_density",
        "current_magnitude_mode": "edge_scalar_abs",
        "cell_reconstructed_midpoint_density": "bernoulli",
        "quasi_fermi_gradient_discretization": "edge_difference",
        "source_volume_policy": "genius_truncated",
        "source_volume_factor": 0.0,
        "source_geometry_scale": 1.0,
        "edge_source_partition": "symmetric",
    }

    def write_minimal_base(self, root: Path) -> Path:
        (root / "mesh.json").write_text(
            json.dumps({"nodes": [], "triangles": [], "contacts": []}),
            encoding="utf-8",
        )
        (root / "doping.csv").write_text(
            "node_id,donors_cm3,acceptors_cm3\n",
            encoding="utf-8",
        )
        (root / "pn2d_sentaurus2018_iv_materials.json").write_text(
            json.dumps({"materials": [{"name": "Si", "ni": 1.0e10}]}),
            encoding="utf-8",
        )
        base = root / "simulation_bv.json"
        base.write_text(
            json.dumps({
                "mesh_file": "mesh.json",
                "node_doping_file": "doping.csv",
                "solver": {},
                "sweep": {},
            }),
            encoding="utf-8",
        )
        return base

    def test_previous_full20_current_approximation_default_remains_legacy_value(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_compensated_default_") as td:
            root = Path(td)
            deck = compare.write_previous_full20_config(
                base_config=self.write_minimal_base(root),
                out_dir=root / "run",
            )

            data = json.loads(deck.read_text(encoding="utf-8"))
            self.assertEqual(
                data["solver"]["impact_ionization"]["current_approximation"],
                "cell_reconstructed",
            )
            self.assertEqual(
                {
                    key: data["solver"]["impact_ionization"][key]
                    for key in self.EXPECTED_IMPACT_BASE
                },
                self.EXPECTED_IMPACT_BASE,
            )

    def test_previous_full20_current_approximation_can_use_density_gradient(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_compensated_override_") as td:
            root = Path(td)
            deck = compare.write_previous_full20_config(
                base_config=self.write_minimal_base(root),
                out_dir=root / "run",
                current_approximation="density_gradient",
            )

            data = json.loads(deck.read_text(encoding="utf-8"))
            self.assertEqual(
                data["solver"]["impact_ionization"]["current_approximation"],
                "density_gradient",
            )

    def test_previous_full20_cli_current_approximation_default_and_override(self) -> None:
        parser = compare.build_parser()

        self.assertEqual(
            parser.parse_args([]).current_approximation,
            "cell_reconstructed",
        )
        self.assertEqual(
            parser.parse_args([
                "--current-approximation",
                "density_gradient",
            ]).current_approximation,
            "density_gradient",
        )

    def test_multibias_indices_and_missing_tdr_gate_are_strict(self) -> None:
        self.assertEqual(
            {
                bias: orchestrator.multibias_index_for_bias(bias)
                for bias in (-12.0, -19.0, -20.0)
            },
            {-12.0: 240, -19.0: 380, -20.0: 400},
        )
        with self.assertRaisesRegex(ValueError, "unsupported replay bias"):
            orchestrator.multibias_index_for_bias(-18.0)

        with tempfile.TemporaryDirectory(prefix="vela_compensated_tdr_") as td:
            source = Path(td)
            expected = source / "pn2d_bv_multibias_0240_des.tdr"
            expected.write_bytes(b"tdr")
            self.assertEqual(
                orchestrator.required_tdr_path(source, -12.0),
                expected,
            )
            with self.assertRaisesRegex(FileNotFoundError, "0380"):
                orchestrator.required_tdr_path(source, -19.0)

    def test_vector_manifest_allows_scalar_sibling_but_requires_one_vector(self) -> None:
        required = [
            manifest_field("ElectrostaticPotential", 1, "V"),
            manifest_field("eQuasiFermiPotential", 1, "V"),
            manifest_field("eDensity", 1, "cm^-3"),
            manifest_field("eMobility", 1, "cm^2*V^-1*s^-1"),
            manifest_field("eAlphaAvalanche", 1, "cm^-1"),
        ]
        scalar = manifest_field("eCurrentDensity", 1, "A*cm^-2")
        vector = manifest_field("eCurrentDensity", 2, "A*cm^-2")

        selected = orchestrator.validate_vector_field_manifest({
            "fields": [*required, scalar, vector],
        })
        self.assertEqual(selected["components"], 2)

        with self.assertRaisesRegex(ValueError, "components=2"):
            orchestrator.validate_vector_field_manifest({
                "fields": [*required, scalar],
            })
        with self.assertRaisesRegex(ValueError, "exactly one"):
            orchestrator.validate_vector_field_manifest({
                "fields": [*required, scalar, vector, dict(vector)],
            })
        wrong_scalar = list(required)
        wrong_scalar[3] = manifest_field("eMobility", 2, "cm^2*V^-1*s^-1")
        with self.assertRaisesRegex(ValueError, "eMobility components=1"):
            orchestrator.validate_vector_field_manifest({
                "fields": [*wrong_scalar, scalar, vector],
            })

    def test_hash_and_artifact_manifest_are_deterministic_and_exclude_self(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_compensated_manifest_") as td:
            root = Path(td)
            a_path = root / "a.txt"
            b_path = root / "b.txt"
            self_path = root / "artifact_manifest.json"
            a_path.write_bytes(b"alpha")
            b_path.write_bytes(b"beta")
            self_path.write_text("old manifest", encoding="utf-8")

            first = orchestrator.build_artifact_manifest(
                out_dir=root,
                git_head="abc123",
                dirty=True,
                parameters={"z": 2, "a": 1},
                commands=[{"argv": ["tool", "--flag"], "cwd": "work", "returncode": 0}],
                artifact_paths=[b_path, self_path, a_path],
            )
            second = orchestrator.build_artifact_manifest(
                out_dir=root,
                git_head="abc123",
                dirty=True,
                parameters={"a": 1, "z": 2},
                commands=[{"argv": ["tool", "--flag"], "cwd": "work", "returncode": 0}],
                artifact_paths=[a_path, b_path, self_path],
            )

            self.assertEqual(first, second)
            self.assertEqual(
                first["schema"],
                "vela.pn2d_compensated_sg_replay.artifact_manifest.v2",
            )
            self.assertEqual(first["schema_version"], 2)
            self.assertEqual(
                [item["path"] for item in first["artifacts"]],
                ["a.txt", "b.txt"],
            )
            self.assertEqual(first["artifacts"][0]["size_bytes"], 5)
            self.assertEqual(
                orchestrator.sha256_file(a_path),
                hashlib.sha256(b"alpha").hexdigest(),
            )
            self.assertNotIn(
                "artifact_manifest.json",
                [item["path"] for item in first["artifacts"]],
            )
            self.assertIn("edge_mapping", first)
    def test_v1_manifest_is_never_reused(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_compensated_v1_cache_") as td:
            root = Path(td)
            manifest_path = root / orchestrator.ARTIFACT_MANIFEST_NAME
            manifest_path.write_text(json.dumps({
                "schema": "vela.pn2d_compensated_sg_replay.artifact_manifest.v1",
                "git_head": "deadbeef",
                "dirty": False,
            }), encoding="utf-8")

            self.assertFalse(orchestrator.reuse_manifest_matches(
                manifest_path,
                out_dir=root,
                git_head="deadbeef",
                dirty=False,
                tdrs=[],
                signature={},
            ))


    def test_variant_specs_form_explicit_two_by_two_matrix(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_compensated_variants_") as td:
            out = Path(td)
            specs = orchestrator.variant_specs(out)

            self.assertEqual(
                list(specs),
                [
                    "legacy_density_gradient",
                    "legacy_gss_midpoint",
                    "reported_density_gradient",
                    "reported_gss_midpoint",
                ],
            )
            self.assertEqual(
                specs["legacy_density_gradient"]["compensated_doping_policy"],
                "dominant_signed_region",
            )
            self.assertEqual(
                specs["reported_gss_midpoint"]["compensated_doping_policy"],
                "reported",
            )
            self.assertEqual(
                specs["legacy_density_gradient"]["imported_dir"],
                out / "variants" / "legacy_density_gradient" / "imported",
            )
            self.assertEqual(
                specs["reported_gss_midpoint"]["run_dir"],
                out / "variants" / "reported_gss_midpoint" / "run",
            )
            self.assertEqual(
                [spec["current_variant"] for spec in specs.values()],
                ["density_gradient", "gss_midpoint", "density_gradient", "gss_midpoint"],
            )
            self.assertEqual(
                [spec["current_approximation"] for spec in specs.values()],
                ["density_gradient", "cell_reconstructed", "density_gradient", "cell_reconstructed"],
            )

    def test_legacy_doping_replay_and_matrix_gate(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_compensated_matrix_") as td:
            specs = orchestrator.variant_specs(Path(td))
            records = {}
            for spec in specs.values():
                imported = spec["imported_dir"]
                vela = imported / "vela"
                vela.mkdir(parents=True)
                (vela / "mesh.json").write_text(
                    json.dumps({"nodes": [{"id": 0}]}),
                    encoding="utf-8",
                )
                (vela / "doping.csv").write_text(
                    "node_id,donors_cm3,acceptors_cm3\n0,1e17,1e17\n",
                    encoding="utf-8",
                )
                resolution = (
                    "signed_aggregate_zero"
                    if spec["doping_strategy"] == "legacy"
                    else "reported"
                )
                (imported / "doping_metadata.json").write_text(
                    json.dumps({
                        "compensated_nodes": {
                            "nodes": [{
                                "node_id": 0,
                                "resolved": False,
                                "resolution_source": resolution,
                            }],
                        },
                    }),
                    encoding="utf-8",
                )
                records[spec["name"]] = (
                    orchestrator.apply_variant_doping_strategy(spec)
                )

            matrix = orchestrator.validate_doping_strategy_matrix(specs)
            self.assertNotEqual(
                matrix["doping_sha256_by_strategy"]["legacy"],
                matrix["doping_sha256_by_strategy"]["reported"],
            )
            self.assertEqual(
                records["legacy_density_gradient"]["transformed_node_ids"],
                [0],
            )
            self.assertEqual(
                records["reported_density_gradient"]["transformed_node_ids"],
                [],
            )

            reported = (
                specs["reported_density_gradient"]["imported_dir"]
                / "vela"
                / "doping.csv"
            ).read_bytes()
            for name in ("legacy_density_gradient", "legacy_gss_midpoint"):
                (
                    specs[name]["imported_dir"] / "vela" / "doping.csv"
                ).write_bytes(reported)
            with self.assertRaisesRegex(ValueError, "collapsed"):
                orchestrator.validate_doping_strategy_matrix(specs)

    def test_command_builders_and_recording_use_argv_arrays(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_compensated_commands_") as td:
            root = Path(td)
            spec = orchestrator.variant_specs(root)["legacy_density_gradient"]
            import_argv = orchestrator.reference_import_command(
                import_script=root / "sentaurus_import.py",
                reference_config=spec["reference_config"],
                source_dir=root / "source",
                imported_dir=spec["imported_dir"],
                importer=root / "sentaurus_import.exe",
                runner=root / "vela_example_runner.exe",
            )
            self.assertEqual(import_argv[0], sys.executable)
            self.assertEqual(import_argv[1:3], [
                str(root / "sentaurus_import.py"),
                "reference",
            ])
            self.assertIn("--skip-vela-run", import_argv)

            diagnostic_argv = orchestrator.diagnostic_command(
                diagnostic_script=root / "diagnose.py",
                out_dir=root,
            )
            self.assertEqual(diagnostic_argv, [
                sys.executable,
                str(root / "diagnose.py"),
                "--variants-root",
                str(root / "variants"),
                "--sentaurus-root",
                str(root / "sentaurus_exports"),
                "--out-dir",
                str(root / "report"),
            ])

            argv = [str(root / "tool.exe"), "--flag", "value"]
            completed = subprocess.CompletedProcess(argv, 0)
            commands: list[dict[str, object]] = []
            with mock.patch.object(
                orchestrator.subprocess,
                "run",
                return_value=completed,
            ) as run:
                returncode = orchestrator.run_recorded_command(
                    argv,
                    root,
                    commands,
                )

            self.assertEqual(returncode, 0)
            run.assert_called_once_with(argv, cwd=root, check=False)
            self.assertEqual(commands, [{
                "argv": argv,
                "cwd": str(root),
                "returncode": 0,
            }])
    def test_prepare_only_writes_configs_and_manifest_without_running_commands(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_compensated_prepare_") as td:
            root = Path(td)
            source = root / "source"
            source.mkdir()
            for bias, index in [(-12.0, 240), (-19.0, 380), (-20.0, 400)]:
                (source / f"pn2d_bv_multibias_{index:04d}_des.tdr").write_bytes(
                    f"tdr {bias:g}".encode("ascii")
                )
            reference = root / "coarse_reference.json"
            reference.write_text(json.dumps({
                "case": "pn2d_sentaurus2018_coarse7x3",
                "device": "pn_diode",
                "tdr_doping": {"compensated_node_policy": "reported"},
                "simulations": [],
            }), encoding="utf-8")
            out = root / "out"
            args = orchestrator.parse_args([
                "--source-dir", str(source),
                "--reference-config", str(reference),
                "--out-dir", str(out),
                "--runner", str(root / "runner.exe"),
                "--importer", str(root / "importer.exe"),
                "--import-script", str(root / "sentaurus_import.py"),
                "--diagnostic-script", str(root / "diagnose.py"),
                "--prepare-only",
            ])
            prepare_calls: list[list[str]] = []

            def prepare_runner(argv, *, cwd, check):
                del cwd, check
                argv = list(argv)
                prepare_calls.append(argv)
                self.assertEqual(argv[0], sys.executable)
                self.assertEqual(argv[2], "reference")
                imported = Path(argv[argv.index("--output-dir") + 1])
                vela = imported / "vela"
                vela.mkdir(parents=True, exist_ok=True)
                (vela / "mesh.json").write_text(
                    json.dumps({"nodes": [], "triangles": [], "contacts": []}),
                    encoding="utf-8",
                )
                (vela / "doping.csv").write_text(
                    "node_id,donors_cm3,acceptors_cm3\n0,1e17,1e17\n",
                    encoding="utf-8",
                )
                config_path = Path(argv[argv.index("--config") + 1])
                policy = json.loads(
                    config_path.read_text(encoding="utf-8")
                )["tdr_doping"]["compensated_node_policy"]
                (imported / "doping_metadata.json").write_text(
                    json.dumps({
                        "compensated_nodes": {
                            "policy": policy,
                            "nodes": [{
                                "node_id": 0,
                                "resolved": False,
                                "resolution_source": (
                                    "signed_aggregate_zero"
                                    if policy == "dominant_signed_region"
                                    else "reported"
                                ),
                            }],
                        },
                    }),
                    encoding="utf-8",
                )
                (vela / "pn2d_sentaurus2018_iv_materials.json").write_text(
                    json.dumps({"materials": [{"name": "Si", "ni": 1.0e10}]}),
                    encoding="utf-8",
                )
                (vela / "simulation_bv.json").write_text(json.dumps({
                    "mesh_file": "mesh.json",
                    "node_doping_file": "doping.csv",
                    "materials_file": "pn2d_sentaurus2018_iv_materials.json",
                    "solver": {},
                    "sweep": {},
                }), encoding="utf-8")
                return subprocess.CompletedProcess(argv, 0)

            with mock.patch.object(
                orchestrator,
                "git_state",
                return_value=("deadbeef", False),
            ):
                manifest = orchestrator.run_reproduction(
                    args,
                    command_runner=prepare_runner,
                )

            self.assertEqual(len(prepare_calls), 4)
            self.assertEqual(
                [(item["bias_V"], item["index"]) for item in manifest["tdrs"]],
                [(-20.0, 400), (-19.0, 380), (-12.0, 240)],
            )
            self.assertTrue(all("sha256" in item for item in manifest["tdrs"]))
            self.assertTrue(all("size_bytes" in item for item in manifest["tdrs"]))
            self.assertEqual(manifest["schema_version"], 2)
            self.assertEqual(len(manifest["commands"]), 12)
            self.assertEqual(manifest["parameters"]["implementation"], "current_head")
            self.assertEqual(
                list(manifest["parameters"]["variants"]),
                [
                    "legacy_density_gradient",
                    "legacy_gss_midpoint",
                    "reported_density_gradient",
                    "reported_gss_midpoint",
                ],
            )
            self.assertEqual(
                sum(command["returncode"] == 0 for command in manifest["commands"]),
                4,
            )
            self.assertEqual(
                sum(command["returncode"] is None for command in manifest["commands"]),
                8,
            )
            specs = orchestrator.variant_specs(out)
            for name, spec in specs.items():
                config = json.loads(
                    spec["reference_config"].read_text(encoding="utf-8")
                )
                self.assertEqual(
                    config["tdr_doping"]["compensated_node_policy"],
                    spec["compensated_doping_policy"],
                )
                deck_path = spec["run_dir"] / spec["deck_name"]
                deck = json.loads(deck_path.read_text(encoding="utf-8"))
                self.assertEqual(
                    deck["solver"]["impact_ionization"],
                    spec["impact_ionization"],
                )
            self.assertEqual(len(manifest["generated"]["decks"]), 4)
            manifest_path = out / "artifact_manifest.json"
            self.assertEqual(
                json.loads(manifest_path.read_text(encoding="utf-8")),
                manifest,
            )
            self.assertNotIn(
                "artifact_manifest.json",
                [item["path"] for item in manifest["artifacts"]],
            )
    def test_run_reproduction_executes_mocked_pipeline_and_reuses_outputs(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_compensated_execute_") as td:
            root = Path(td)
            source = root / "source"
            source.mkdir()
            for bias, index in [(-12.0, 240), (-19.0, 380), (-20.0, 400)]:
                (source / f"pn2d_bv_multibias_{index:04d}_des.tdr").write_bytes(
                    f"tdr {bias:g}".encode("ascii")
                )
            reference = root / "coarse_reference.json"
            reference.write_text(json.dumps({
                "case": "pn2d_sentaurus2018_coarse7x3",
                "device": "pn_diode",
                "tdr_doping": {"compensated_node_policy": "reported"},
                "simulations": [],
            }), encoding="utf-8")
            out = root / "out"
            importer = root / "sentaurus_import.exe"
            runner = root / "vela_example_runner.exe"
            import_script = root / "sentaurus_import.py"
            diagnostic_script = root / "diagnose.py"
            for tool in (importer, runner, import_script, diagnostic_script):
                tool.write_bytes(tool.name.encode("utf-8"))

            args = orchestrator.parse_args([
                "--source-dir", str(source),
                "--reference-config", str(reference),
                "--out-dir", str(out),
                "--runner", str(runner),
                "--importer", str(importer),
                "--import-script", str(import_script),
                "--diagnostic-script", str(diagnostic_script),
            ])
            calls: list[list[str]] = []

            def fake_runner(argv, *, cwd, check):
                del check
                argv = list(argv)
                calls.append(argv)
                if argv[0] == str(importer):
                    export_dir = Path(argv[argv.index("--export-dir") + 1])
                    export_dir.mkdir(parents=True, exist_ok=True)
                    (export_dir / "field_manifest.json").write_text(json.dumps({
                        "fields": [
                            manifest_field("ElectrostaticPotential", 1, "V"),
                            manifest_field("eQuasiFermiPotential", 1, "V"),
                            manifest_field("eDensity", 1, "cm^-3"),
                            manifest_field("eMobility", 1, "cm^2*V^-1*s^-1"),
                            manifest_field("eAlphaAvalanche", 1, "cm^-1"),
                            manifest_field("eCurrentDensity", 1, "A*cm^-2"),
                            manifest_field("eCurrentDensity", 2, "A*cm^-2"),
                        ],
                    }), encoding="utf-8")
                elif len(argv) > 2 and argv[0] == sys.executable and argv[2] == "reference":
                    imported = Path(argv[argv.index("--output-dir") + 1])
                    vela = imported / "vela"
                    vela.mkdir(parents=True, exist_ok=True)
                    (vela / "mesh.json").write_text(json.dumps({
                        "nodes": [],
                        "triangles": [],
                        "contacts": [],
                    }), encoding="utf-8")
                    (vela / "doping.csv").write_text(
                        "node_id,donors_cm3,acceptors_cm3\n0,1e17,1e17\n",
                        encoding="utf-8",
                    )
                    config_path = Path(argv[argv.index("--config") + 1])
                    policy = json.loads(
                        config_path.read_text(encoding="utf-8")
                    )["tdr_doping"]["compensated_node_policy"]
                    (imported / "doping_metadata.json").write_text(
                        json.dumps({
                            "compensated_nodes": {
                                "policy": policy,
                                "nodes": [{
                                    "node_id": 0,
                                    "resolved": False,
                                    "resolution_source": (
                                        "signed_aggregate_zero"
                                        if policy == "dominant_signed_region"
                                        else "reported"
                                    ),
                                }],
                            },
                        }),
                        encoding="utf-8",
                    )
                    (vela / "pn2d_sentaurus2018_iv_materials.json").write_text(
                        json.dumps({"materials": [{"name": "Si", "ni": 1.0e10}]}),
                        encoding="utf-8",
                    )
                    (vela / "simulation_bv.json").write_text(json.dumps({
                        "mesh_file": "mesh.json",
                        "node_doping_file": "doping.csv",
                        "materials_file": "pn2d_sentaurus2018_iv_materials.json",
                        "solver": {},
                        "sweep": {},
                    }), encoding="utf-8")
                elif argv[0] == str(runner):
                    deck_path = Path(argv[argv.index("--config") + 1])
                    self.assertTrue(deck_path.is_file())
                    deck = json.loads(deck_path.read_text(encoding="utf-8"))
                    run_dir = Path(cwd)
                    (run_dir / deck["output_csv"]).write_text(
                        "bias_V,current_total_A_per_um\n0,0\n",
                        encoding="utf-8",
                    )
                    sg = Path(
                        deck["sweep"]["diagnostics"]["sg_avalanche_edges"]["csv_file"]
                    )
                    sg.parent.mkdir(parents=True, exist_ok=True)
                    sg.write_text("bias_V,edge_id\n-20,9\n", encoding="utf-8")
                    vtk = run_dir / (
                        deck["sweep"]["vtk_prefix"] + "_0400_-20V.vtk"
                    )
                    vtk.parent.mkdir(parents=True, exist_ok=True)
                    vtk.write_text("# vtk DataFile Version 2.0\n", encoding="utf-8")
                elif len(argv) > 1 and argv[1] == str(diagnostic_script):
                    report = Path(argv[argv.index("--out-dir") + 1])
                    report.mkdir(parents=True, exist_ok=True)
                    (report / "compensated_sg_replay.csv").write_text(
                        "variant,bias_V,y_um,side,edge_id\n"
                        "legacy_density_gradient,-12,0,left,9\n",
                        encoding="utf-8",
                    )
                    (report / "compensated_sg_replay.json").write_text(
                        json.dumps({"rows": 1}),
                        encoding="utf-8",
                    )
                    (report / "compensated_sg_replay_report.md").write_text(
                        "# replay\n", encoding="utf-8"
                    )
                else:
                    self.fail(f"unexpected command: {argv}")
                return subprocess.CompletedProcess(argv, 0)

            with mock.patch.object(
                orchestrator,
                "git_state",
                return_value=("deadbeef", False),
            ):
                manifest = orchestrator.run_reproduction(
                    args,
                    command_runner=fake_runner,
                )

            self.assertEqual(len(calls), 12)
            self.assertEqual(len(manifest["commands"]), 12)
            self.assertEqual(manifest["parameters"]["implementation"], "current_head")
            self.assertEqual(
                list(manifest["parameters"]["variants"]),
                [
                    "legacy_density_gradient",
                    "legacy_gss_midpoint",
                    "reported_density_gradient",
                    "reported_gss_midpoint",
                ],
            )
            self.assertTrue(all(
                command["returncode"] == 0 for command in manifest["commands"]
            ))
            specs = orchestrator.variant_specs(out)
            for spec in specs.values():
                deck_path = spec["run_dir"] / spec["deck_name"]
                deck = json.loads(deck_path.read_text(encoding="utf-8"))
                self.assertEqual(
                    deck["solver"]["impact_ionization"],
                    spec["impact_ionization"],
                )
                self.assertEqual(len(deck["sweep"]["bias_points"]), 401)
                self.assertEqual(deck["sweep"]["bias_points"][0], 0.0)
                self.assertEqual(deck["sweep"]["bias_points"][-1], -20.0)
                self.assertEqual(deck["sweep"]["step"], -0.05)
                self.assertTrue(deck["sweep"]["write_vtk"])
                self.assertTrue(
                    deck["sweep"]["diagnostics"]["sg_avalanche_edges"]["enabled"]
                )
            self.assertEqual(len(manifest["generated"]["decks"]), 4)
            self.assertEqual(len(manifest["generated"]["model_configs"]), 4)
            self.assertEqual(len(manifest["inputs"]["mesh"]), 4)
            self.assertEqual(len(manifest["inputs"]["doping"]), 4)
            self.assertEqual(len(manifest["inputs"]["materials"]), 4)
            self.assertEqual(manifest["edge_mapping"], [{
                "variant": "legacy_density_gradient",
                "bias_V": "-12",
                "y_um": "0",
                "side": "left",
                "edge_id": "9",
            }])

            args.reuse_existing = True
            reuse_runner = mock.Mock(side_effect=AssertionError("must reuse"))
            with mock.patch.object(
                orchestrator,
                "git_state",
                return_value=("deadbeef", False),
            ):
                reused = orchestrator.run_reproduction(
                    args,
                    command_runner=reuse_runner,
                )
            reuse_runner.assert_not_called()
            self.assertEqual(len(reused["commands"]), 12)
            self.assertEqual(reused["invocation"]["commands"], [])
            self.assertTrue(reused["invocation"]["reuse_accepted"])

            diagnostic_script.write_bytes(b"changed diagnostic")
            calls.clear()
            with mock.patch.object(
                orchestrator,
                "git_state",
                return_value=("deadbeef", False),
            ):
                refreshed_tool = orchestrator.run_reproduction(
                    args,
                    command_runner=fake_runner,
                )
            self.assertEqual(len(calls), 12)
            self.assertFalse(refreshed_tool["parameters"]["reuse_accepted"])

            (source / "pn2d_bv_multibias_0240_des.tdr").write_bytes(b"changed")
            calls.clear()
            with mock.patch.object(
                orchestrator,
                "git_state",
                return_value=("deadbeef", False),
            ):
                refreshed = orchestrator.run_reproduction(
                    args,
                    command_runner=fake_runner,
                )
            self.assertEqual(len(calls), 12)
            self.assertFalse(refreshed["parameters"]["reuse_accepted"])
    def test_main_parses_cli_and_delegates_to_reproduction(self) -> None:
        with mock.patch.object(
            orchestrator,
            "run_reproduction",
            return_value={"schema": orchestrator.MANIFEST_SCHEMA},
        ) as run, mock.patch("builtins.print") as output:
            returncode = orchestrator.main([
                "--out-dir", "custom-out",
                "--prepare-only",
            ])

        self.assertEqual(returncode, 0)
        parsed = run.call_args.args[0]
        self.assertEqual(parsed.out_dir, Path("custom-out"))
        self.assertTrue(parsed.prepare_only)
        output.assert_called_once()
    def test_orchestrator_parser_defaults_and_path_overrides(self) -> None:
        defaults = orchestrator.parse_args([])
        self.assertEqual(defaults.biases, [-12.0, -19.0, -20.0])
        self.assertEqual(defaults.source_dir, orchestrator.DEFAULT_SOURCE_DIR)
        self.assertEqual(defaults.reference_config, orchestrator.DEFAULT_REFERENCE_CONFIG)
        self.assertEqual(defaults.runner, orchestrator.DEFAULT_RUNNER)
        self.assertEqual(defaults.importer, orchestrator.DEFAULT_IMPORTER)
        self.assertFalse(defaults.prepare_only)
        self.assertFalse(defaults.skip_export)
        self.assertFalse(defaults.skip_vela_run)
        self.assertFalse(defaults.reuse_existing)

        overridden = orchestrator.parse_args([
            "--source-dir", "raw",
            "--reference-config", "coarse.json",
            "--out-dir", "repro",
            "--runner", "runner.exe",
            "--importer", "importer.exe",
            "--biases=-20,-12",
            "--prepare-only",
            "--skip-export",
            "--skip-vela-run",
            "--reuse-existing",
        ])
        self.assertEqual(overridden.source_dir, Path("raw"))
        self.assertEqual(overridden.reference_config, Path("coarse.json"))
        self.assertEqual(overridden.out_dir, Path("repro"))
        self.assertEqual(overridden.runner, Path("runner.exe"))
        self.assertEqual(overridden.importer, Path("importer.exe"))
        self.assertEqual(overridden.biases, [-20.0, -12.0])
        self.assertTrue(overridden.prepare_only)
        self.assertTrue(overridden.skip_export)
        self.assertTrue(overridden.skip_vela_run)
        self.assertTrue(overridden.reuse_existing)


if __name__ == "__main__":
    unittest.main()
