import hashlib
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(Path(__file__).resolve().parent))
from pn2d_minimal6_inverse_sealer_fixture import (  # noqa: E402
    fake_importer, make_fixture, mutate_json, sha256, write_json,
)
from scripts.pn2d_minimal6_diagnostics.inverse_contracts import SampleStatus  # noqa: E402
from scripts.pn2d_minimal6_diagnostics.inverse_inputs import load_input_bundle  # noqa: E402


MODULE_PATH = REPO / "scripts" / "pn2d_minimal6_diagnostics" / "inverse_input_sealer.py"
if MODULE_PATH.is_file():
    SPEC = importlib.util.spec_from_file_location("inverse_input_sealer", MODULE_PATH)
    sealer = importlib.util.module_from_spec(SPEC)
    assert SPEC.loader is not None
    sys.modules[SPEC.name] = sealer
    SPEC.loader.exec_module(sealer)
else:
    sealer = None


class PN2DMinimal6InverseInputSealerContractTest(unittest.TestCase):
    def test_sealer_module_exists(self) -> None:
        self.assertIsNotNone(sealer, "inverse_input_sealer.py is missing")

    @unittest.skipIf(sealer is None, "sealer not implemented")
    def test_seals_realistic_roots_with_missing_vela_derived_fields(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            vela, sentaurus, supplemental, importer, runner = make_fixture(base)
            roots = sealer.seal_inverse_input_roots(
                vela, sentaurus, supplemental, base / "sealed", importer=importer,
                vela_executable=runner, phase_base="a5524cf", importer_runner=fake_importer,
            )
            bundle = load_input_bundle(roots["vela"], roots["sentaurus"], roots["supplemental"])
            self.assertEqual(len(bundle.common_keys), 40)
            vela_field = next(row for row in bundle.observations
                              if row.solver == "vela" and row.topology == "sketch"
                              and row.bias_V == -1.0 and row.support_id == "0"
                              and row.quantity == "ElectricField")
            self.assertIs(vela_field.status, SampleStatus.MISSING_FIELD)
            self.assertIsNone(vela_field.raw_value)
            density = next(row for row in bundle.observations
                           if row.solver == "vela" and row.topology == "sketch"
                           and row.bias_V == -1.0 and row.support_id == "0"
                           and row.quantity == "eDensity")
            self.assertEqual(density.value_si, 1.0e18)
            sent_field = next(row for row in bundle.observations
                              if row.solver == "sentaurus" and row.topology == "sketch"
                              and row.bias_V == -1.0 and row.support_id == "0"
                              and row.quantity == "ElectricField" and row.component == "component0")
            self.assertEqual(sent_field.raw_value, 6001.0)
            mobility = next(row for row in bundle.observations
                            if row.solver == "sentaurus" and row.topology == "sketch"
                            and row.bias_V == -1.0 and row.support_id == "0"
                            and row.quantity == "eMobility")
            self.assertEqual(mobility.raw_value, 1001.0)
            vela_manifest = json.loads((roots["vela"] / "manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(vela_manifest["provenance"]["execution_binding_status"],
                             "post_hoc_observed")
            sent_manifest = json.loads((roots["sentaurus"] / "manifest.json").read_text(encoding="utf-8"))
            self.assertEqual(sent_manifest["provenance"]["remote_solver_binding_status"],
                             "not_declared_by_source_manifest")

    @unittest.skipIf(sealer is None, "sealer not implemented")
    def test_dual_outputs_are_byte_identical_and_old_exports_are_irrelevant(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            vela, sentaurus, supplemental, importer, runner = make_fixture(base)
            old = sentaurus / "sentaurus/sketch/exports/sketch_m1p000000/ignored.csv"
            old.write_text("first\n", encoding="utf-8")
            first = base / "first"
            sealer.seal_inverse_input_roots(vela, sentaurus, supplemental, first,
                                            importer=importer, vela_executable=runner,
                                            phase_base="a5524cf", importer_runner=fake_importer)
            old.write_text("tampered but untrusted\n", encoding="utf-8")
            second = base / "second"
            sealer.seal_inverse_input_roots(vela, sentaurus, supplemental, second,
                                            importer=importer, vela_executable=runner,
                                            phase_base="a5524cf", importer_runner=fake_importer)

            def members(root: Path) -> dict[str, bytes]:
                return {path.relative_to(root).as_posix(): path.read_bytes()
                        for path in root.rglob("*") if path.is_file()}

            self.assertEqual(members(first), members(second))

    @unittest.skipIf(sealer is None, "sealer not implemented")
    def test_rejects_source_contract_errors_before_output(self) -> None:
        mutations = (
            ("incomplete", lambda v, s, x: mutate_json(
                x / "manifest.json", lambda j: j.update(outputs_complete=False))),
            ("state hash", lambda v, s, x: mutate_json(
                v / "sweep_manifest.json",
                lambda j: j["accepted_checkpoints"][0].update(state_sha256="0" * 64))),
            ("path escape", lambda v, s, x: mutate_json(
                v / "sweep_manifest.json",
                lambda j: j["accepted_checkpoints"][0].update(state_path="../escape.csv"))),
            ("version", lambda v, s, x: mutate_json(
                x / "manifest.json", lambda j: j.update(sentaurus_version="P-2019.03"))),
            ("interpolation", lambda v, s, x: mutate_json(
                v / "sweep_manifest.json", lambda j: j.update(interpolation="linear"))),
        )
        for label, mutate in mutations:
            with self.subTest(label=label), tempfile.TemporaryDirectory() as tmp:
                base = Path(tmp)
                vela, sentaurus, supplemental, importer, runner = make_fixture(base)
                mutate(vela, sentaurus, supplemental)
                out = base / "sealed"
                with self.assertRaises(ValueError):
                    sealer.seal_inverse_input_roots(
                        vela, sentaurus, supplemental, out, importer=importer,
                        vela_executable=runner, phase_base="a5524cf",
                        importer_runner=fake_importer)
                self.assertFalse(out.exists())

    @unittest.skipIf(sealer is None, "sealer not implemented")
    def test_rejects_existing_output_without_modification(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            vela, sentaurus, supplemental, importer, runner = make_fixture(base)
            out = base / "sealed"
            out.mkdir()
            marker = out / "keep.txt"
            marker.write_text("keep", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "already exists"):
                sealer.seal_inverse_input_roots(
                    vela, sentaurus, supplemental, out, importer=importer,
                    vela_executable=runner, phase_base="a5524cf", importer_runner=fake_importer)
            self.assertEqual(marker.read_text(encoding="utf-8"), "keep")

    @unittest.skipIf(sealer is None, "sealer not implemented")
    def test_rejects_imported_unit_mapping_and_geometry(self) -> None:
        cases = (
            (lambda raw: raw["fields"]["ElectricField"].update(unit="mV"), "unit"),
            (lambda raw: raw["fields"]["ElectricField"].update(mapping_status="partial"), "mapping"),
            (lambda raw: raw["nodes"][0].update(x_um=9.0), "coordinate"),
        )
        for change, message in cases:
            with self.subTest(message=message), tempfile.TemporaryDirectory() as tmp:
                base = Path(tmp)
                vela, sentaurus, supplemental, importer, runner = make_fixture(base)
                checkpoint = sentaurus / "sentaurus/sketch/checkpoints/sketch_m1p000000.tdr"
                raw = json.loads(checkpoint.read_text(encoding="utf-8"))
                change(raw)
                write_json(checkpoint, raw)
                manifest_path = sentaurus / "sweep_manifest.json"
                manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
                row = next(row for row in manifest["accepted_checkpoints"]
                           if row["topology"] == "sketch" and row["target_bias_V"] == -1.0)
                row["state_sha256"] = sha256(checkpoint)
                write_json(manifest_path, manifest)
                out = base / "sealed"
                with self.assertRaisesRegex(ValueError, message):
                    sealer.seal_inverse_input_roots(
                        vela, sentaurus, supplemental, out, importer=importer,
                        vela_executable=runner, phase_base="a5524cf",
                        importer_runner=fake_importer)
                self.assertFalse(out.exists())

    @unittest.skipIf(sealer is None, "sealer not implemented")
    def test_rejects_importer_failure_and_tampered_seal(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            vela, sentaurus, supplemental, importer, runner = make_fixture(base)

            def failed_importer(tdr: Path, output: Path) -> None:
                raise RuntimeError("injected importer failure")

            with self.assertRaisesRegex(ValueError, "importer"):
                sealer.seal_inverse_input_roots(
                    vela, sentaurus, supplemental, base / "failed", importer=importer,
                    vela_executable=runner, phase_base="a5524cf",
                    importer_runner=failed_importer)
            self.assertFalse((base / "failed").exists())
            roots = sealer.seal_inverse_input_roots(
                vela, sentaurus, supplemental, base / "sealed", importer=importer,
                vela_executable=runner, phase_base="a5524cf", importer_runner=fake_importer)
            target = roots["vela"] / "states/sketch/m1V.csv"
            target.write_text("tampered", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "state hash mismatch"):
                load_input_bundle(roots["vela"], roots["sentaurus"], roots["supplemental"])


if __name__ == "__main__":
    unittest.main()
