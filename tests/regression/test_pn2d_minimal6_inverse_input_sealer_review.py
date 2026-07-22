import csv
import json
import shutil
import tempfile
import unittest
from pathlib import Path

from scripts.pn2d_minimal6_diagnostics import inverse_input_sealer as sealer
from tests.regression.pn2d_minimal6_inverse_sealer_fixture import (
    SOURCE_IDS,
    SUPPLEMENTAL,
    fake_importer,
    make_fixture,
    mesh,
    mutate_json,
    raw_tdr,
    sha256,
    write_json,
)


def _bind_supplemental_exports(root: Path, importer: Path) -> dict:
    manifest_path = root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["importer"] = str(importer.resolve())
    for state in manifest["states"]:
        state_root = root / "states" / state["topology_id"] / state["bias_tag"]
        export = state_root / "export"
        shutil.rmtree(export)
        tdr = state_root / "artifacts" / state["final_tdr_name"]
        fake_importer(tdr, export)
        state["member_sha256"] = {
            path.relative_to(export).as_posix(): sha256(path)
            for path in export.rglob("*") if path.is_file()
        }
    write_json(manifest_path, manifest)
    return manifest


def _imported_fixture(base: Path) -> tuple[Path, dict]:
    tdr = base / "state.tdr"
    export = base / "export"
    mesh_path = base / "mesh.json"
    write_json(tdr, raw_tdr("sketch", -1.0, SUPPLEMENTAL))
    write_json(mesh_path, mesh("sketch"))
    fake_importer(tdr, export)
    return export, sealer._mesh_contract(mesh_path)


def _rewrite_csv(path: Path, mutate) -> None:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
        fieldnames = list(rows[0])
    mutate(rows)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


class PN2DMinimal6InverseInputSealerReviewTest(unittest.TestCase):
    def test_rejects_importer_identity_mismatch_before_import(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            vela, sentaurus, supplemental, importer, runner = make_fixture(base)
            _bind_supplemental_exports(supplemental, importer)
            other = base / "other_sentaurus_import.exe"
            other.write_bytes(b"other importer")
            mutate_json(supplemental / "manifest.json",
                        lambda value: value.update(importer=str(other.resolve())))
            calls = []

            def forbidden_import(tdr: Path, output: Path) -> None:
                calls.append(tdr)
                raise AssertionError("import must not run")

            with self.assertRaisesRegex(ValueError, "supplemental importer identity"):
                sealer.seal_inverse_input_roots(
                    vela, sentaurus, supplemental, base / "sealed", importer=importer,
                    vela_executable=runner, phase_base="a5524cf",
                    importer_runner=forbidden_import)
            self.assertEqual(calls, [])

    def test_rejects_mutated_unhashed_supplemental_tdr(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            vela, sentaurus, supplemental, importer, runner = make_fixture(base)
            _bind_supplemental_exports(supplemental, importer)
            tdr = supplemental / "states/sketch/m1V/artifacts/pn2d_minimal6_state_m1V.tdr"
            mutate_json(tdr, lambda value: value["fields"]["eMobility"]["values"]["5"]
                        .__setitem__(0, 987654.0))
            with self.assertRaisesRegex(ValueError, "supplemental reimport mismatch"):
                sealer.seal_inverse_input_roots(
                    vela, sentaurus, supplemental, base / "sealed", importer=importer,
                    vela_executable=runner, phase_base="a5524cf",
                    importer_runner=fake_importer)
            self.assertFalse((base / "sealed").exists())

    def test_rejects_reledgered_supplemental_export_divergence(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            vela, sentaurus, supplemental, importer, runner = make_fixture(base)
            manifest = _bind_supplemental_exports(supplemental, importer)
            state = next(row for row in manifest["states"]
                         if row["topology_id"] == "sketch" and row["requested_bias_V"] == -1.0)
            export = supplemental / "states/sketch/m1V/export"
            field = export / "fields/eMobility_region0.csv"
            _rewrite_csv(field, lambda rows: rows[0].update(component0="987654"))
            state["member_sha256"]["fields/eMobility_region0.csv"] = sha256(field)
            write_json(supplemental / "manifest.json", manifest)
            with self.assertRaisesRegex(ValueError, "supplemental reimport mismatch"):
                sealer.seal_inverse_input_roots(
                    vela, sentaurus, supplemental, base / "sealed", importer=importer,
                    vela_executable=runner, phase_base="a5524cf",
                    importer_runner=fake_importer)

    def test_seal_binds_declared_importer_sha256(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            vela, sentaurus, supplemental, importer, runner = make_fixture(base)
            _bind_supplemental_exports(supplemental, importer)
            roots = sealer.seal_inverse_input_roots(
                vela, sentaurus, supplemental, base / "sealed", importer=importer,
                vela_executable=runner, phase_base="a5524cf",
                importer_runner=fake_importer)
            manifest = json.loads((roots["supplemental"] / "manifest.json")
                                  .read_text(encoding="utf-8"))
            self.assertEqual(manifest["provenance"]["declared_importer"],
                             str(importer.resolve()))
            self.assertEqual(manifest["provenance"]["importer_sha256"], sha256(importer))

    def test_rejects_contact_voltage_contract_and_role_mutations(self) -> None:
        cases = (
            ("unit", lambda export: mutate_json(
                export / "field_manifest.json",
                lambda value: next(row for row in value["fields"]
                                   if row["name"] == "ContactExternalVoltage")
                .update(unit="mV"))),
            ("components", lambda export: mutate_json(
                export / "field_manifest.json",
                lambda value: next(row for row in value["fields"]
                                   if row["name"] == "ContactExternalVoltage")
                .update(components=2))),
            ("mapping", lambda export: mutate_json(
                export / "field_manifest.json",
                lambda value: next(row for row in value["fields"]
                                   if row["name"] == "ContactExternalVoltage")
                .update(mapping_status="partial"))),
            ("region role", lambda export: mutate_json(
                export / "field_manifest.json",
                lambda value: next(row for row in value["fields"]
                                   if row["name"] == "ContactExternalVoltage"
                                   and row["region"] == 1).update(region_name="Anode"))),
            ("scalar role", lambda export: (
                _rewrite_csv(export / "fields/ContactExternalVoltage_region1.csv",
                             lambda rows: rows[0].update(component0="-1")),
                _rewrite_csv(export / "fields/ContactExternalVoltage_region2.csv",
                             lambda rows: rows[0].update(component0="0")))),
            ("contact region", lambda export: _rewrite_csv(
                export / "contacts.csv", lambda rows: rows[0].update(region="R.Other"))),
        )
        for label, mutate in cases:
            with self.subTest(label=label), tempfile.TemporaryDirectory() as tmp:
                export, mesh_contract = _imported_fixture(Path(tmp))
                mutate(export)
                with self.assertRaisesRegex(ValueError, "contact"):
                    sealer._validate_import(export, SUPPLEMENTAL, mesh_contract, -1.0)

    def test_rejects_duplicate_ids_names_and_connectivity_multiplicity(self) -> None:
        import_cases = (
            ("element id", lambda export: _rewrite_csv(
                export / "elements.csv", lambda rows: rows[1].update(id=rows[0]["id"]))),
            ("element multiplicity", lambda export: _rewrite_csv(
                export / "elements.csv", lambda rows: rows.append(dict(rows[0], id="99")))),
            ("contact name", lambda export: _rewrite_csv(
                export / "contacts.csv", lambda rows: rows.append(dict(rows[0])))),
            ("contact node", lambda export: _rewrite_csv(
                export / "contacts.csv",
                lambda rows: rows[0].update(node_ids=rows[0]["node_ids"] + ";" +
                                            rows[0]["node_ids"].split(";")[0]))),
        )
        for label, mutate in import_cases:
            with self.subTest(source="import", label=label), tempfile.TemporaryDirectory() as tmp:
                export, mesh_contract = _imported_fixture(Path(tmp))
                mutate(export)
                with self.assertRaisesRegex(ValueError, "duplicate|topology"):
                    sealer._validate_import(export, SUPPLEMENTAL, mesh_contract, -1.0)

        mesh_cases = (
            ("element id", lambda value: value["triangles"][1].update(
                id=value["triangles"][0]["id"])),
            ("contact id", lambda value: value["contacts"][1].update(
                id=value["contacts"][0]["id"])),
            ("contact name", lambda value: value["contacts"].append(
                dict(value["contacts"][0], id=99))),
            ("contact node", lambda value: value["contacts"][0]["node_ids"].append(0)),
        )
        for label, mutate in mesh_cases:
            with self.subTest(source="mesh", label=label), tempfile.TemporaryDirectory() as tmp:
                path = Path(tmp) / "mesh.json"
                value = mesh("sketch")
                mutate(value)
                write_json(path, value)
                with self.assertRaisesRegex(ValueError, "duplicate|contact"):
                    sealer._mesh_contract(path)


if __name__ == "__main__":
    unittest.main()
