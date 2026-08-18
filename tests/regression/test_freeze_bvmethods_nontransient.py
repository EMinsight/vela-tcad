from __future__ import annotations

import importlib.util
from pathlib import Path
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "freeze_bvmethods_nontransient",
    ROOT / "scripts" / "freeze_bvmethods_nontransient.py")
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class FreezeBvmethodsNontransientTest(unittest.TestCase):
    def test_scope_excludes_only_transient(self) -> None:
        self.assertEqual(
            ["ABA_poisson", "ABA_coupled", "resistor", "voltage2current", "continuation"],
            list(MODULE.METHODS),
        )
        self.assertNotIn("transient", MODULE.METHODS)

    def test_continuation_crossing_is_interpolated(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "curve.csv"
            path.write_text(
                "bias_V,current_total_A_per_um,converged\n"
                "6.3,8e-5,1\n6.5,1.2e-4,1\n", encoding="utf-8")
            result = MODULE.extract_continuation_bv(path)
            self.assertAlmostEqual(6.4, result["vela_bv_V"])
            self.assertEqual("pass", result["status"])


if __name__ == "__main__":
    unittest.main()
