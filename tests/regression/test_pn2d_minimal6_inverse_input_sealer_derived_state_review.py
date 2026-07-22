import tempfile
import unittest
from pathlib import Path

from scripts.export_pn2d_minimal6_states import write_state_csv
from scripts.pn2d_minimal6_diagnostics import inverse_input_sealer as sealer
from tests.regression.pn2d_minimal6_inverse_sealer_fixture import (
    REQUIRED,
    SUPPLEMENTAL,
    fake_importer,
    raw_tdr,
    sha256,
    write_json,
)


def _exports_with_derived_state(base: Path) -> tuple[Path, Path, dict[str, str]]:
    tdr = base / "state.tdr"
    integral_fields = {
        "eIonIntegral": (1, "1"),
        "hIonIntegral": (1, "1"),
        "MeanIonIntegral": (1, "1"),
    }
    write_json(
        tdr,
        raw_tdr("sketch", -1.0, {**REQUIRED, **SUPPLEMENTAL, **integral_fields}),
    )
    trusted = base / "trusted"
    fresh = base / "fresh"
    fake_importer(tdr, trusted)
    write_state_csv(trusted)
    fake_importer(tdr, fresh)
    ledger = {path.relative_to(trusted).as_posix(): sha256(path)
              for path in trusted.rglob("*") if path.is_file()}
    return trusted, fresh, ledger


class PN2DMinimal6InverseInputSealerDerivedStateReviewTest(unittest.TestCase):
    def test_regenerates_declared_derived_state_before_exact_ledger_check(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            _, fresh, ledger = _exports_with_derived_state(Path(tmp))
            sealer._require_reimport_matches_ledger(fresh, ledger)
            self.assertEqual(sha256(fresh / "state.csv"), ledger["state.csv"])
            self.assertEqual(
                {path.relative_to(fresh).as_posix() for path in fresh.rglob("*") if path.is_file()},
                set(ledger),
            )

    def test_rejects_reledgered_tampered_derived_state(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            trusted, fresh, ledger = _exports_with_derived_state(Path(tmp))
            state = trusted / "state.csv"
            state.write_text(state.read_text(encoding="utf-8") + "tampered\n", encoding="utf-8")
            ledger["state.csv"] = sha256(state)
            with self.assertRaisesRegex(ValueError, "state.csv"):
                sealer._require_reimport_matches_ledger(fresh, ledger)


if __name__ == "__main__":
    unittest.main()
