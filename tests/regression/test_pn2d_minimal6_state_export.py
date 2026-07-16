import csv
import importlib.util
import json
import math
import shutil
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
MODULE_PATH = REPO / "scripts" / "export_pn2d_minimal6_states.py"
SPEC = importlib.util.spec_from_file_location("export_pn2d_minimal6_states", MODULE_PATH)
export = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(export)


SCALAR_FIELDS = {
    "ElectrostaticPotential": "V",
    "eQuasiFermiPotential": "V",
    "hQuasiFermiPotential": "V",
    "eDensity": "cm^-3",
    "hDensity": "cm^-3",
    "eMobility": "cm^2*V^-1*s^-1",
    "hMobility": "cm^2*V^-1*s^-1",
    "eAlphaAvalanche": "cm^-1",
    "hAlphaAvalanche": "cm^-1",
    "ImpactIonization": "cm^-3*s^-1",
    "eVelocity": "cm*s^-1", "hVelocity": "cm*s^-1",
    "eIonIntegral": "1", "hIonIntegral": "1", "MeanIonIntegral": "1",
}
VECTOR_FIELDS = {
    "ElectricField": "V*cm^-1",
    "eCurrentDensity": "A*cm^-2",
    "hCurrentDensity": "A*cm^-2",
}


def valid_field_manifest() -> dict[str, object]:
    fields = []
    for name, unit in SCALAR_FIELDS.items():
        fields.append({
            "name": name,
            "region": 0,
            "components": 1,
            "unit": unit,
            "mapping_status": "complete",
            "global_node_mapping": "global_vertex_order",
        })
    for name, unit in VECTOR_FIELDS.items():
        fields.append({
            "name": name,
            "region": 0,
            "components": 2,
            "unit": unit,
            "mapping_status": "complete",
            "global_node_mapping": "global_vertex_order",
        })
    return {"fields": fields}


def write_csv(path: Path, header: list[str], rows: list[list[object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def write_neutral_state(root: Path, missing: str | None = None) -> None:
    root.mkdir(parents=True, exist_ok=True)
    (root / "field_manifest.json").write_text(
        json.dumps(valid_field_manifest(), indent=2) + "\n", encoding="utf-8"
    )
    # Deliberately reverse importer IDs. Canonical order must come from exact coordinates.
    write_csv(
        root / "nodes.csv",
        ["id", "x_um", "y_um"],
        [[15 - i, x, y] for i, (x, y) in enumerate(
            [(0.0, 0.5), (1.0, 0.5), (2.0, 0.5),
             (2.0, 0.0), (0.0, 0.0), (1.0, 0.0)]
        )],
    )
    values = {
        "ElectrostaticPotential": [0.01 * i for i in range(6)],
        "eQuasiFermiPotential": [0.1 + 0.01 * i for i in range(6)],
        "hQuasiFermiPotential": [-0.2 - 0.01 * i for i in range(6)],
        "eDensity": [1.0e10 * (i + 1) for i in range(6)],
        "hDensity": [2.0e10 * (i + 1) for i in range(6)],
    }
    source_ids = [15 - i for i in range(6)]
    for name, field_values in values.items():
        if name == missing:
            continue
        write_csv(
            root / "fields" / f"{name}_region0.csv",
            ["node_id", "component0"],
            [[node_id, value] for node_id, value in zip(source_ids, field_values)],
        )


def write_recovery_candidate(root: Path) -> tuple[Path, str]:
    topology_contracts = {
        "sketch": [[1, 5, 2], [5, 6, 2], [2, 6, 4], [2, 4, 3]],
        "mirror": [[1, 5, 6], [1, 6, 2], [2, 6, 3], [6, 4, 3]],
    }
    states = []
    replays = []
    for topology in ("sketch", "mirror"):
        for bias in (0.0, -12.0, -19.0):
            tag = export._bias_tag(bias)
            export_dir = root / "states" / topology / tag / "export"
            export_dir.mkdir(parents=True)
            (export_dir / "field_manifest.json").write_text(
                json.dumps(valid_field_manifest(), indent=2) + "\n", encoding="utf-8"
            )
            (export_dir / "raw_member.txt").write_text(
                f"{topology} {bias:g}\n", encoding="utf-8"
            )
            audit_members = {
                "mesh.json": "{}\n",
                "doping.csv": "node_id,doping_cm3\n0,0\n",
                "state.csv": "node_id,potential_V\n0,0\n",
                "audit.json": "{}\n",
                "vela_node_state.csv": "node_id\n0\n",
                "vela_edge_audit.csv": "edge_id\n0\n",
                "vela_triangle_audit.csv": "triangle_id\n0\n",
            }
            for name, contents in audit_members.items():
                (export_dir / name).write_text(contents, encoding="utf-8")
            option_names = (
                ("--mesh", "mesh.json"),
                ("--doping", "doping.csv"),
                ("--state", "state.csv"),
                ("--config", "audit.json"),
                ("--node-out", "vela_node_state.csv"),
                ("--edge-out", "vela_edge_audit.csv"),
                ("--triangle-out", "vela_triangle_audit.csv"),
            )
            arguments = [
                item
                for option, name in option_names
                for item in (option, str((export_dir / name).resolve()))
            ]
            replays.append({
                "topology_id": topology,
                "bias_V": bias,
                "producer": "build-release/pn2d_minimal6_operator_audit.exe",
                "command": " ".join([
                    "build-release/pn2d_minimal6_operator_audit.exe", *arguments
                ]),
                "arguments": arguments,
                "exit_code": 0,
                "stdout": "",
                "stderr": "",
                "input_sha256": {
                    str((export_dir / name).resolve()): export._sha256(export_dir / name)
                    for _, name in option_names[:4]
                },
                "output_sha256": {
                    str((export_dir / name).resolve()): export._sha256(export_dir / name)
                    for _, name in option_names[4:]
                },
            })
            states.append({
                "topology_id": topology,
                "requested_bias_V": bias,
                "actual_bias_V": bias,
                "bias_tag": tag,
                "status": "passed",
                "export_dir": str(export_dir.relative_to(root)),
                "member_sha256": export.collect_member_hashes(export_dir),
                "topology_contract": {
                    "nodes": 6,
                    "triangles": 4,
                    "edges": 9,
                    "contact_edges": {"Anode": [1, 5], "Cathode": [3, 4]},
                    "triangle_connectivity": topology_contracts[topology],
                },
            })
    manifest = {
        "schema": "vela.pn2d_minimal6_states.v1",
        "run_id": "candidate-v1",
        "outputs_complete": True,
        "states": states,
        "task4_provenance": {
            "producer": "build-release/pn2d_minimal6_operator_audit.exe",
            "producer_sha256": export._sha256(
                REPO / "build-release" / "pn2d_minimal6_operator_audit.exe"
            ),
            "task4_source_commit": "1" * 40,
            "replay_environment": "test fixture",
            "replays": replays,
        },
    }
    manifest_path = root / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    members = export.collect_member_hashes(root)
    recovery_validation = {
        "schema": "vela.pn2d_minimal6_recovery_validation.v1",
        "run_id": manifest["run_id"],
        "generated_at": "2026-07-16T00:00:00+08:00",
        "archive_root": str(root.resolve()),
        "archive_manifest_sha256": export._sha256(manifest_path),
        "member_count": len(members),
        "member_sha256": members,
        "validation": {
            "outputs_complete": True,
            "exact_state_matrix": True,
            "field_contract": True,
            "member_hashes_verified": True,
        },
    }
    (root / "recovery_validation.json").write_text(
        json.dumps(recovery_validation, indent=2) + "\n", encoding="utf-8"
    )
    return manifest_path, export._sha256(manifest_path)


class PN2DMinimal6StateExportTest(unittest.TestCase):
    def test_v2_native_avalanche_identity_uses_raw_impact_ionization(self) -> None:
        fields = export.validate_v2_field_manifest(valid_field_manifest())
        native = fields["sentaurus_native_avalanche_generation"]
        self.assertEqual(
            native,
            {
                "normalized_name": "sentaurus_native_avalanche_generation",
                "raw_name": "ImpactIonization",
                "components": 1,
                "unit": "cm^-3*s^-1",
                "semantic_role": "native_avalanche_generation",
            },
        )

    def test_v2_velocity_and_ionization_integrals_are_installed_scalars(self) -> None:
        fields = export.validate_v2_field_manifest(valid_field_manifest())
        expected = {
            "sentaurus_electron_speed": ("eVelocity", "cm*s^-1", "carrier_speed"),
            "sentaurus_hole_speed": ("hVelocity", "cm*s^-1", "carrier_speed"),
            "sentaurus_electron_ionization_integral": (
                "eIonIntegral", "1", "path_ionization_integral"
            ),
            "sentaurus_hole_ionization_integral": (
                "hIonIntegral", "1", "path_ionization_integral"
            ),
            "sentaurus_mean_ionization_integral": (
                "MeanIonIntegral", "1", "path_ionization_integral"
            ),
        }
        for normalized, (raw, unit, role) in expected.items():
            with self.subTest(normalized=normalized):
                self.assertEqual(fields[normalized]["raw_name"], raw)
                self.assertEqual(fields[normalized]["components"], 1)
                self.assertEqual(fields[normalized]["unit"], unit)
                self.assertEqual(fields[normalized]["semantic_role"], role)

    def test_v2_selects_one_exact_vector_match_among_scalar_aliases(self) -> None:
        manifest = valid_field_manifest()
        vector = next(
            field for field in manifest["fields"] if field["name"] == "ElectricField"
        )
        scalar = dict(vector)
        scalar["components"] = 1
        manifest["fields"].append(scalar)
        fields = export.validate_v2_field_manifest(manifest)
        self.assertEqual(fields["sentaurus_electric_field"]["components"], 2)

        manifest["fields"].append(dict(vector))
        with self.assertRaisesRegex(ValueError, "duplicate alias.*ElectricField"):
            export.validate_v2_field_manifest(manifest)

    def test_v2_fields_require_one_complete_region_zero_alias(self) -> None:
        manifest = valid_field_manifest()
        manifest["fields"].append(dict(
            next(field for field in manifest["fields"] if field["name"] == "ImpactIonization")
        ))
        with self.assertRaisesRegex(ValueError, "duplicate alias.*ImpactIonization"):
            export.validate_v2_field_manifest(manifest)

        for key, value in (("region", 1), ("mapping_status", "partial"),
                           ("global_node_mapping", "region_node_order")):
            manifest = valid_field_manifest()
            native = next(
                field for field in manifest["fields"] if field["name"] == "ImpactIonization"
            )
            native[key] = value
            with self.subTest(key=key):
                with self.assertRaisesRegex(ValueError, rf"ImpactIonization.*{key}"):
                    export.validate_v2_field_manifest(manifest)

    def test_v1_contract_remains_unchanged_when_v2_is_added(self) -> None:
        self.assertEqual(export.SCHEMA, "vela.pn2d_minimal6_states.v1")
        with tempfile.TemporaryDirectory() as tmp:
            manifest = export.prepare_exports(
                topology_ids=("sketch", "mirror"),
                biases=(0.0, -12.0, -19.0),
                run_id="v1-compatibility",
                output_dir=Path(tmp),
                ssh_target="sentaurus",
            )
        self.assertEqual(manifest["schema"], "vela.pn2d_minimal6_states.v1")

    def test_v2_seal_preserves_exact_bias_artifact_paths_and_hashes(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "candidate"
            source.mkdir()
            manifest_path, digest = write_recovery_candidate(source)
            original = manifest_path.read_bytes()
            sealed = Path(tmp) / "candidate-sealed-v2"
            manifest = export.seal_recovered_archive(
                source, sealed, digest, run_id="candidate-sealed-v2",
                source_kind="regenerated",
            )
            validated = export.validate_sealed_archive(sealed)

            self.assertEqual(manifest_path.read_bytes(), original)
            self.assertEqual(manifest, validated)
            self.assertEqual(manifest["schema"], "vela.pn2d_minimal6_states.v2")
            self.assertEqual(manifest["source_manifest_sha256"], digest)
            source_validation = manifest["source_validation"]
            self.assertEqual(source_validation["kind"], "regenerated")
            self.assertEqual(
                source_validation["sha256"],
                export._sha256(source / "recovery_validation.json"),
            )
            self.assertTrue((sealed / source_validation["path"]).is_file())
            self.assertTrue(source_validation["validation"]["member_hashes_verified"])
            self.assertEqual(len(manifest["states"]), 6)
            self.assertEqual(
                {(state["topology_id"], state["requested_bias_V"], state["actual_bias_V"])
                 for state in manifest["states"]},
                {(topology, bias, bias) for topology in ("sketch", "mirror")
                 for bias in (0.0, -12.0, -19.0)},
            )
            for state in manifest["states"]:
                self.assertTrue(state["raw_artifacts"])
                self.assertTrue(all(
                    set(artifact) == {"path", "sha256"}
                    and not Path(artifact["path"]).is_absolute()
                    and len(artifact["sha256"]) == 64
                    for artifact in state["raw_artifacts"]
                ))
                self.assertTrue(all(
                    set(field) == {
                        "normalized_name", "raw_name", "components", "unit", "semantic_role"
                    }
                    for field in state["fields"]
                ))

    def test_v2_manifest_and_member_hash_mutation_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "candidate"
            source.mkdir()
            _, digest = write_recovery_candidate(source)
            sealed = Path(tmp) / "candidate-sealed-v2"
            export.seal_recovered_archive(
                source, sealed, digest, run_id="candidate-sealed-v2",
                source_kind="regenerated",
            )

            source_manifest = sealed / "source_manifest.json"
            source_manifest.write_bytes(source_manifest.read_bytes() + b"\n")
            with self.assertRaisesRegex(ValueError, "source manifest hash mismatch"):
                export.validate_sealed_archive(sealed)

            shutil.rmtree(sealed)
            export.seal_recovered_archive(
                source, sealed, digest, run_id="candidate-sealed-v2",
                source_kind="regenerated",
            )
            manifest = json.loads((sealed / "manifest.json").read_text(encoding="utf-8"))
            member = sealed / manifest["states"][0]["raw_artifacts"][0]["path"]
            member.write_bytes(member.read_bytes() + b"mutation")
            with self.assertRaisesRegex(
                ValueError, "raw artifact hash mismatch|sealed source member ledger hash"
            ):
                export.validate_sealed_archive(sealed)

            shutil.rmtree(sealed)
            export.seal_recovered_archive(
                source, sealed, digest, run_id="candidate-sealed-v2",
                source_kind="regenerated",
            )
            sidecar = sealed / "source_recovery_validation.json"
            sidecar.write_bytes(sidecar.read_bytes() + b"mutation")
            with self.assertRaisesRegex(
                ValueError, "sealed source member ledger validation hash"
            ):
                export.validate_sealed_archive(sealed)

            shutil.rmtree(sealed)
            export.seal_recovered_archive(
                source, sealed, digest, run_id="candidate-sealed-v2",
                source_kind="regenerated",
            )
            path = sealed / "manifest.json"
            manifest = json.loads(path.read_text(encoding="utf-8"))
            manifest["states"][0]["requested_bias_V"] = -11.999999999
            manifest["states"][0]["actual_bias_V"] = -11.999999999
            path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "exact six-state matrix mismatch"):
                export.validate_sealed_archive(sealed)

    def test_v2_rejects_noncanonical_state_paths(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "candidate"
            source.mkdir()
            _, digest = write_recovery_candidate(source)
            sealed = Path(tmp) / "candidate-sealed-v2"
            export.seal_recovered_archive(
                source, sealed, digest, run_id=sealed.name, source_kind="regenerated"
            )
            path = sealed / "manifest.json"
            baseline = json.loads(path.read_text(encoding="utf-8"))
            variants = []
            for field in ("export_dir", "field_manifest", "state_csv"):
                cross_state = json.loads(json.dumps(baseline))
                cross_state["states"][0][field] = baseline["states"][1][field]
                variants.append((f"cross-state {field}", cross_state))
                absolute = json.loads(json.dumps(baseline))
                absolute["states"][0][field] = str(
                    (sealed / baseline["states"][0][field]).resolve()
                )
                variants.append((f"absolute {field}", absolute))
            for label, manifest in variants:
                with self.subTest(label=label):
                    path.write_text(
                        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
                    )
                    with self.assertRaisesRegex(ValueError, "canonical state path"):
                        export.validate_sealed_archive(sealed)

    def test_v2_requires_complete_six_state_audit_provenance(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "candidate"
            source.mkdir()
            _, digest = write_recovery_candidate(source)
            sealed = Path(tmp) / "candidate-sealed-v2"
            export.seal_recovered_archive(
                source, sealed, digest, run_id=sealed.name, source_kind="regenerated"
            )
            path = sealed / "manifest.json"
            baseline = json.loads(path.read_text(encoding="utf-8"))
            missing = json.loads(json.dumps(baseline))
            missing.pop("audit_provenance")
            short = json.loads(json.dumps(baseline))
            short["audit_provenance"]["replays"].pop()
            bad_argv = json.loads(json.dumps(baseline))
            bad_argv["audit_provenance"]["replays"][0]["arguments"] = "not-an-array"
            bad_reference = json.loads(json.dumps(baseline))
            replay = bad_reference["audit_provenance"]["replays"][0]
            replay["arguments"][7] = baseline["audit_provenance"]["replays"][1][
                "arguments"
            ][7]
            bad_hash = json.loads(json.dumps(baseline))
            bad_hash["audit_provenance"]["producer_sha256"] = "0" * 64
            variants = (
                ("missing", missing),
                ("short", short),
                ("bad argv", bad_argv),
                ("cross-state audit reference", bad_reference),
                ("bad executable hash", bad_hash),
            )
            for label, manifest in variants:
                with self.subTest(label=label):
                    path.write_text(
                        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
                    )
                    with self.assertRaisesRegex(
                        ValueError, "audit provenance|schema validation"
                    ):
                        export.validate_sealed_archive(sealed)

    def test_local_recovery_is_bound_to_named_run_and_manifest_hash(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "candidate"
            source.mkdir()
            _, digest = write_recovery_candidate(source)
            sealed = Path(tmp) / "candidate-sealed-v2"
            with self.assertRaisesRegex(ValueError, "named local recovery"):
                export.seal_recovered_archive(
                    source, sealed, digest, run_id=sealed.name,
                    source_kind="local_recovery",
                )
            export.seal_recovered_archive(
                source, sealed, digest, run_id=sealed.name,
                source_kind="regenerated",
            )
            path = sealed / "manifest.json"
            manifest = json.loads(path.read_text(encoding="utf-8"))
            manifest["source_kind"] = "local_recovery"
            path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "named local recovery"):
                export.validate_sealed_archive(sealed)


    def test_source_member_ledger_fails_closed_for_all_set_mutations(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            mutations = (
                "missing-sidecar",
                "wrong-sidecar-hash",
                "altered",
                "extra",
                "missing",
                "manifest-binding",
                "false-validation",
            )
            for mutation in mutations:
                with self.subTest(mutation=mutation):
                    source = Path(tmp) / f"candidate-{mutation}"
                    source.mkdir()
                    _, digest = write_recovery_candidate(source)
                    export_dir = (
                        source / "states" / "sketch" / "p0V" / "export"
                    )
                    sidecar_path = source / "recovery_validation.json"
                    if mutation == "missing-sidecar":
                        sidecar_path.unlink()
                    elif mutation == "altered":
                        (export_dir / "raw_member.txt").write_text(
                            "altered\n", encoding="utf-8"
                        )
                    elif mutation == "extra":
                        (export_dir / "unexpected.txt").write_text(
                            "extra\n", encoding="utf-8"
                        )
                    elif mutation == "missing":
                        (export_dir / "raw_member.txt").unlink()
                    elif mutation in ("manifest-binding", "false-validation"):
                        sidecar = json.loads(sidecar_path.read_text(encoding="utf-8"))
                        if mutation == "manifest-binding":
                            sidecar["archive_manifest_sha256"] = "0" * 64
                        else:
                            sidecar["validation"]["member_hashes_verified"] = False
                        sidecar_path.write_text(
                            json.dumps(sidecar, indent=2) + "\n", encoding="utf-8"
                        )
                    sealed = Path(tmp) / f"sealed-{mutation}"
                    with self.assertRaisesRegex(
                        ValueError, "source member ledger|recovery validation"
                    ):
                        export.seal_recovered_archive(
                            source, sealed, digest, run_id=sealed.name,
                            source_kind="regenerated",
                            expected_source_validation_sha256=(
                                "0" * 64 if mutation == "wrong-sidecar-hash" else None
                            ),
                        )

    def test_v2_rejects_duplicate_artifact_paths_before_hash_mapping(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / "candidate"
            source.mkdir()
            _, digest = write_recovery_candidate(source)
            sealed = Path(tmp) / "candidate-sealed-v2"
            export.seal_recovered_archive(
                source, sealed, digest, run_id=sealed.name, source_kind="regenerated"
            )
            path = sealed / "manifest.json"
            manifest = json.loads(path.read_text(encoding="utf-8"))
            artifacts = manifest["states"][0]["raw_artifacts"]
            artifacts.insert(0, {"path": artifacts[0]["path"], "sha256": "0" * 64})
            path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "duplicate artifact path"):
                export.validate_sealed_archive(sealed)


    def test_required_state_matrix_is_exact_and_complete(self) -> None:
        states = [
            {"topology_id": topology, "requested_bias_V": bias, "actual_bias_V": bias,
             "status": "passed"}
            for topology in ("sketch", "mirror")
            for bias in (0.0, -12.0, -19.0)
        ]
        matrix = export.validate_state_matrix(states)
        self.assertEqual(
            set(matrix),
            {(t, v) for t in ("sketch", "mirror") for v in (0.0, -12.0, -19.0)},
        )

    def test_nearest_bias_is_rejected_at_stricter_than_one_picovolt(self) -> None:
        states = [
            {"topology_id": topology, "requested_bias_V": bias,
             "actual_bias_V": bias + (2.0e-12 if topology == "mirror" and bias == -19.0 else 0.0),
             "status": "passed"}
            for topology in ("sketch", "mirror")
            for bias in (0.0, -12.0, -19.0)
        ]
        with self.assertRaisesRegex(ValueError, r"mirror.*-19.*1e-12 V"):
            export.validate_state_matrix(states)

    def test_final_contact_voltage_uses_strict_one_picovolt_boundary(self) -> None:
        self.assertEqual(export.validate_final_bias(0.0, 1.0e-12), 1.0e-12)
        inside = math.nextafter(1.0e-12, 0.0)
        self.assertEqual(export.validate_final_bias(0.0, inside), inside)
        review_example = math.nextafter(-19.0 - 1.0e-12, -math.inf)
        review_error = abs(review_example - (-19.0))
        self.assertEqual(review_error, 1.0018652574217413e-12)
        self.assertGreater(review_error, export.BIAS_TOLERANCE_V)
        with self.assertRaisesRegex(ValueError, r"within 1e-12 V"):
            export.validate_final_bias(-19.0, review_example)

    def test_vector_fields_require_two_components_and_cm_units(self) -> None:
        manifest = valid_field_manifest()
        current = next(f for f in manifest["fields"] if f["name"] == "eCurrentDensity")
        current["unit"] = "A*m^-2"
        with self.assertRaisesRegex(ValueError, r"eCurrentDensity.*A\*cm\^-2"):
            export.validate_field_manifest(manifest)
        current["unit"] = "A*cm^-2"
        current["components"] = 1
        with self.assertRaisesRegex(ValueError, r"eCurrentDensity.*components=2"):
            export.validate_field_manifest(manifest)

    def test_every_field_requires_complete_global_region_zero_mapping(self) -> None:
        for key, bad in (("region", 1), ("mapping_status", "partial"),
                         ("global_node_mapping", "region_node_order")):
            manifest = valid_field_manifest()
            field = next(f for f in manifest["fields"] if f["name"] == "ElectricField")
            field[key] = bad
            with self.subTest(key=key):
                with self.assertRaisesRegex(ValueError, rf"ElectricField.*{key}"):
                    export.validate_field_manifest(manifest)

    def test_scalar_field_units_are_explicit(self) -> None:
        for name, expected in SCALAR_FIELDS.items():
            manifest = valid_field_manifest()
            field = next(f for f in manifest["fields"] if f["name"] == name)
            field["unit"] = "wrong"
            with self.subTest(name=name):
                with self.assertRaisesRegex(ValueError, f"{name}.*{expected.replace('*', r'\*').replace('^', r'\^')}"):
                    export.validate_field_manifest(manifest)

    def test_missing_quasi_fermi_field_is_fatal(self) -> None:
        for field, role in (("eQuasiFermiPotential", "phin"),
                            ("hQuasiFermiPotential", "phip")):
            with self.subTest(field=field), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                write_neutral_state(root, missing=field)
                with self.assertRaisesRegex(ValueError, rf"missing.*{field}.*{role}"):
                    export.write_state_csv(root, export.canonical_minimal6_coordinates())

    def test_state_csv_uses_canonical_order_and_si_values(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            write_neutral_state(root)
            path = export.write_state_csv(root, export.canonical_minimal6_coordinates())
            with path.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
        self.assertEqual(
            list(rows[0]),
            ["node_id", "psi_V", "phin_V", "phip_V", "n_m3", "p_m3"],
        )
        self.assertEqual([int(row["node_id"]) for row in rows], list(range(6)))
        self.assertAlmostEqual(float(rows[3]["psi_V"]), 0.03)
        self.assertAlmostEqual(float(rows[3]["phin_V"]), 0.13)
        self.assertAlmostEqual(float(rows[3]["phip_V"]), -0.23)
        self.assertAlmostEqual(float(rows[3]["n_m3"]), 4.0e16)
        self.assertAlmostEqual(float(rows[3]["p_m3"]), 8.0e16)

    def test_deck_generation_is_exact_bias_and_explicit_grid_only(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            manifest = export.prepare_exports(
                topology_ids=("sketch", "mirror"),
                biases=(0.0, -12.0, -19.0),
                run_id="minimal6_states_test",
                output_dir=Path(tmp),
                ssh_target="sentaurus",
            )
            self.assertEqual(manifest["schema"], "vela.pn2d_minimal6_states.v1")
            self.assertEqual(len(manifest["states"]), 6)
            for state in manifest["states"]:
                command_text = " ".join(state["remote_commands"])
                self.assertIn(
                    "tdx -d pn2d_minimal6.grd pn2d_minimal6.dat pn2d_minimal6.tdr",
                    command_text,
                )
                self.assertNotRegex(command_text.lower(), r"\b(?:sde|snmesh)\b")
                self.assertNotIn("nearest", command_text.lower())
                deck = Path(state["bundle_dir"], state["deck_name"]).read_text(encoding="utf-8")
                self.assertNotIn("__TARGET_BIAS_V__", deck)
                self.assertIn(f'Voltage={state["requested_bias_V"]:.17g}', deck)
                self.assertIn('Grid      = "pn2d_minimal6.tdr"', deck)
                self.assertIn('Doping    = "pn2d_minimal6.tdr"', deck)
                self.assertIn("ElectricField/Vector", deck)
                self.assertIn("eCurrentDensity/Vector", deck)
                self.assertIn("hCurrentDensity/Vector", deck)
                plot_stem = f"pn2d_minimal6_state_{state['bias_tag']}"
                self.assertIn(f'Output    = "{plot_stem}"', deck)
                self.assertNotIn(f'Output    = "{plot_stem}.log"', deck)
                self.assertNotIn(".log_des.log", deck)
                self.assertEqual(state["log_name"], f"{plot_stem}_des.log")

    def test_returned_files_include_exact_remote_stdout_without_globs(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            manifest = export.prepare_exports(
                topology_ids=("sketch", "mirror"),
                biases=(0.0, -12.0, -19.0),
                run_id="minimal6_states_returned_files",
                output_dir=Path(tmp),
                ssh_target="sentaurus",
            )
            state = manifest["states"][0]
        plot_stem = f"pn2d_minimal6_state_{state['bias_tag']}"
        stdout_name = f"run_{plot_stem}.out"
        returned = state.get("returned_files")
        self.assertIsInstance(returned, list)
        self.assertIn(stdout_name, returned)
        self.assertEqual(state.get("stdout_name"), stdout_name)
        self.assertEqual(state["log_name"], f"{plot_stem}_des.log")
        self.assertFalse(any(set(name) & set("*?[]") for name in returned))

    def test_sdevice_failure_recovers_exact_remote_stdout_with_argv_arrays(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = export.prepare_exports(
                topology_ids=("sketch", "mirror"),
                biases=(0.0, -12.0, -19.0),
                run_id="minimal6_states_stdout_recovery",
                output_dir=root,
                ssh_target="sentaurus",
            )
            state = manifest["states"][0]
            importer = root / "sentaurus_import.exe"
            importer.touch()
            calls: list[list[str]] = []
            original_run_checked = export.run_checked

            def fail_sdevice(argv: list[str]) -> None:
                calls.append(list(argv))
                if argv[0] == "ssh-test" and "sdevice" in argv[-1]:
                    raise RuntimeError("mock sdevice failure")

            export.run_checked = fail_sdevice
            try:
                with self.assertRaisesRegex(RuntimeError, "mock sdevice failure"):
                    export._live_executor(
                        state,
                        ssh_bin="ssh-test",
                        scp_bin="scp-test",
                        ssh_target="sentaurus",
                        importer=importer,
                    )
            finally:
                export.run_checked = original_run_checked
        stdout_name = f"run_pn2d_minimal6_state_{state['bias_tag']}.out"
        recovered = [argv for argv in calls if argv and argv[0] == "scp-test"]
        self.assertTrue(all(isinstance(argv, list) for argv in calls))
        self.assertTrue(any(argv[1].endswith(f"/{stdout_name}") for argv in recovered))
    def test_nonfinite_actual_biases_are_fail_closed(self) -> None:
        for actual in (math.nan, math.inf, -math.inf):
            with self.subTest(actual=actual), tempfile.TemporaryDirectory() as tmp:
                manifest = export.prepare_exports(
                    topology_ids=("sketch", "mirror"),
                    biases=(0.0, -12.0, -19.0),
                    run_id="minimal6_states_nonfinite_actual",
                    output_dir=Path(tmp),
                    ssh_target="sentaurus",
                )

                def executor(state: dict[str, object]) -> dict[str, object]:
                    neutral = Path(str(state["export_dir"]))
                    write_neutral_state(neutral)
                    return {"actual_bias_V": actual, "export_dir": str(neutral)}

                with self.assertRaisesRegex(ValueError, "not finite"):
                    export.run_exports(manifest, executor=executor)
                saved = json.loads(Path(manifest["manifest_path"]).read_text(encoding="utf-8"))
                self.assertFalse(saved["outputs_complete"])
                self.assertEqual(saved["states"][0]["status"], "failed")
                self.assertTrue(all(
                    state["status"] == "prepared" for state in saved["states"][1:]
                ))
    def test_missing_actual_bias_is_fail_closed_and_stops_later_states(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            manifest = export.prepare_exports(
                topology_ids=("sketch", "mirror"),
                biases=(0.0, -12.0, -19.0),
                run_id="minimal6_states_missing_actual_bias",
                output_dir=Path(tmp),
                ssh_target="sentaurus",
            )
            calls: list[tuple[str, float]] = []

            def executor(state: dict[str, object]) -> dict[str, object]:
                key = (str(state["topology_id"]), float(state["requested_bias_V"]))
                calls.append(key)
                neutral = Path(str(state["export_dir"]))
                write_neutral_state(neutral)
                return {"export_dir": str(neutral)}

            with self.assertRaisesRegex(ValueError, "missing actual_bias_V"):
                export.run_exports(manifest, executor=executor)
            saved = json.loads(Path(manifest["manifest_path"]).read_text(encoding="utf-8"))
        self.assertFalse(saved["outputs_complete"])
        self.assertEqual(calls, [("sketch", 0.0)])
        self.assertEqual(saved["states"][0]["status"], "failed")
        self.assertTrue(all(state["status"] == "prepared" for state in saved["states"][1:]))
    def test_missing_neutral_outputs_are_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            manifest = export.prepare_exports(
                topology_ids=("sketch", "mirror"),
                biases=(0.0, -12.0, -19.0),
                run_id="minimal6_states_missing_output",
                output_dir=Path(tmp),
                ssh_target="sentaurus",
            )
            with self.assertRaisesRegex(ValueError, "missing neutral export directory"):
                export.run_exports(
                    manifest,
                    executor=lambda state: {"actual_bias_V": state["requested_bias_V"]},
                )
            saved = json.loads(Path(manifest["manifest_path"]).read_text(encoding="utf-8"))
        self.assertFalse(saved["outputs_complete"])
        self.assertEqual(saved["states"][0]["status"], "failed")
    def test_partial_manifest_is_fail_closed_and_stops_later_states(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = export.prepare_exports(
                topology_ids=("sketch", "mirror"),
                biases=(0.0, -12.0, -19.0),
                run_id="minimal6_states_partial",
                output_dir=root,
                ssh_target="sentaurus",
            )
            calls: list[tuple[str, float]] = []

            def executor(state: dict[str, object]) -> dict[str, object]:
                key = (str(state["topology_id"]), float(state["requested_bias_V"]))
                calls.append(key)
                if key == ("sketch", -12.0):
                    raise RuntimeError("mock convergence failure")
                neutral = Path(str(state["export_dir"]))
                write_neutral_state(neutral)
                return {"actual_bias_V": key[1], "export_dir": str(neutral)}

            with self.assertRaisesRegex(RuntimeError, "mock convergence failure"):
                export.run_exports(manifest, executor=executor)
            saved = json.loads(Path(manifest["manifest_path"]).read_text(encoding="utf-8"))
        self.assertFalse(saved["outputs_complete"])
        self.assertEqual(calls, [("sketch", 0.0), ("sketch", -12.0)])
        self.assertEqual(saved["states"][1]["status"], "failed")
        self.assertTrue(all(state["status"] == "prepared" for state in saved["states"][2:]))

    def test_mocked_six_state_export_writes_complete_si_states(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            manifest = export.prepare_exports(
                topology_ids=("sketch", "mirror"),
                biases=(0.0, -12.0, -19.0),
                run_id="minimal6_states_mocked",
                output_dir=Path(tmp),
                ssh_target="sentaurus",
            )

            def executor(state: dict[str, object]) -> dict[str, object]:
                neutral = Path(str(state["export_dir"]))
                write_neutral_state(neutral)
                return {
                    "actual_bias_V": float(state["requested_bias_V"]),
                    "export_dir": str(neutral),
                }

            export.run_exports(manifest, executor=executor)
            saved = json.loads(Path(manifest["manifest_path"]).read_text(encoding="utf-8"))
            matrix = export.validate_state_matrix(saved["states"])
            state_paths = [Path(state["state_csv"]) for state in saved["states"]]
        self.assertTrue(saved["outputs_complete"])
        self.assertEqual(len(matrix), 6)
        self.assertTrue(all(path.name == "state.csv" for path in state_paths))
        self.assertTrue(all("member_sha256" in state for state in saved["states"]))

    def test_member_hashes_fail_closed_after_mutation(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            member = root / "member.txt"
            member.write_text("original", encoding="utf-8")
            hashes = export.collect_member_hashes(root)
            export.validate_member_hashes(root, hashes)
            member.write_text("changed", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "hash mismatch"):
                export.validate_member_hashes(root, hashes)

    def test_recovered_archive_requires_manifest_hash_and_six_states(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = {"outputs_complete": True, "states": [
                {"topology_id":t,"requested_bias_V":b,"actual_bias_V":b,"status":"passed"}
                for t in ("sketch","mirror") for b in (0.0,-12.0,-19.0)]}
            path = root / "manifest.json"; path.write_text(json.dumps(manifest), encoding="utf-8")
            digest = export._sha256(path)
            export.validate_recovered_archive(root, digest)
            with self.assertRaisesRegex(ValueError, "manifest hash mismatch"):
                export.validate_recovered_archive(root, "0" * 64)

if __name__ == "__main__":
    unittest.main()
