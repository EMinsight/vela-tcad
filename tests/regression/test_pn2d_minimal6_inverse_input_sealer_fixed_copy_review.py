import json
import tempfile
import unittest
from pathlib import Path

from scripts.pn2d_minimal6_diagnostics import inverse_input_sealer as sealer
from tests.regression.pn2d_minimal6_inverse_sealer_fixture import (
    fake_importer,
    make_fixture,
    sha256,
)
from tests.regression.test_pn2d_minimal6_inverse_input_sealer_review import (
    _bind_supplemental_exports,
)


class PN2DMinimal6InverseInputSealerFixedCopyReviewTest(unittest.TestCase):
    def test_rejects_staging_importer_mutation_after_first_import(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            vela, sentaurus, supplemental, importer, runner = make_fixture(base)
            _bind_supplemental_exports(supplemental, importer)
            calls = []

            def mutating_importer(tdr: Path, output: Path) -> None:
                fake_importer(tdr, output)
                calls.append(tdr)
                if len(calls) == 1:
                    fixed = output.parents[3] / "_tools" / importer.name
                    self.assertTrue(fixed.is_file())
                    fixed.write_bytes(b"mutated staging importer")

            target = base / "sealed"
            with self.assertRaisesRegex(ValueError, "importer.*SHA|importer.*changed"):
                sealer.seal_inverse_input_roots(
                    vela, sentaurus, supplemental, target, importer=importer,
                    vela_executable=runner, phase_base="a5524cf",
                    importer_runner=mutating_importer)
            self.assertEqual(len(calls), 1)
            self.assertFalse(target.exists())
            self.assertEqual(list(base.glob(".sealed.staging-*")), [])

    def test_sealed_executable_provenance_equals_member_ledger(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            vela, sentaurus, supplemental, importer, runner = make_fixture(base)
            _bind_supplemental_exports(supplemental, importer)
            expected = sha256(importer)
            roots = sealer.seal_inverse_input_roots(
                vela, sentaurus, supplemental, base / "sealed", importer=importer,
                vela_executable=runner, phase_base="a5524cf",
                importer_runner=fake_importer)
            for name in ("sentaurus", "supplemental"):
                manifest = json.loads((roots[name] / "manifest.json").read_text(encoding="utf-8"))
                member = f"source/executables/{importer.name}"
                self.assertEqual(manifest["provenance"]["executable_sha256"], expected)
                self.assertEqual(manifest["member_sha256"][member], expected)


if __name__ == "__main__":
    unittest.main()
