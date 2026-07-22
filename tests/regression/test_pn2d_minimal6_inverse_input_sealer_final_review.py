import tempfile
import unittest
from pathlib import Path

from scripts.pn2d_minimal6_diagnostics import inverse_input_sealer as sealer
from tests.regression.pn2d_minimal6_inverse_sealer_fixture import (
    fake_importer,
    make_fixture,
    mutate_json,
)
from tests.regression.test_pn2d_minimal6_inverse_input_sealer_review import (
    _bind_supplemental_exports,
)


class PN2DMinimal6InverseInputSealerFinalReviewTest(unittest.TestCase):
    def test_rejects_importer_mutation_after_first_import(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            base = Path(tmp)
            vela, sentaurus, supplemental, importer, runner = make_fixture(base)
            _bind_supplemental_exports(supplemental, importer)
            calls = []

            def mutating_importer(tdr: Path, output: Path) -> None:
                fake_importer(tdr, output)
                calls.append(tdr)
                if len(calls) == 1:
                    importer.write_bytes(b"mutated importer")

            target = base / "sealed"
            with self.assertRaisesRegex(ValueError, "importer.*SHA|importer.*changed"):
                sealer.seal_inverse_input_roots(
                    vela, sentaurus, supplemental, target, importer=importer,
                    vela_executable=runner, phase_base="a5524cf",
                    importer_runner=mutating_importer)
            self.assertEqual(len(calls), 1)
            self.assertFalse(target.exists())
            self.assertEqual(list(base.glob(".sealed.staging-*")), [])

    def test_rejects_noncanonical_expected_matrix_rows(self) -> None:
        cases = (
            ("string bias", lambda rows: rows[0].__setitem__(1, "-1.0")),
            ("integer bias", lambda rows: rows[0].__setitem__(1, -1)),
            ("boolean bias", lambda rows: rows[0].__setitem__(1, True)),
            ("boolean topology", lambda rows: rows[0].__setitem__(0, True)),
            ("extra value", lambda rows: rows[0].append("extra")),
            ("missing value", lambda rows: rows[0].pop()),
            ("mapping row", lambda rows: rows.__setitem__(
                0, {"topology": "sketch", "bias": -1.0})),
        )
        for label, mutate in cases:
            with self.subTest(label=label), tempfile.TemporaryDirectory() as tmp:
                base = Path(tmp)
                _, _, supplemental, _, _ = make_fixture(base)
                mutate_json(supplemental / "manifest.json",
                            lambda value: mutate(value["expected_matrix"]))
                with self.assertRaisesRegex(ValueError, "expected matrix"):
                    sealer._validate_supplemental(supplemental)


if __name__ == "__main__":
    unittest.main()
