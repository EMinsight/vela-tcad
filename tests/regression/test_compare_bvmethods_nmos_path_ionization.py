from __future__ import annotations

import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[2] / "scripts/compare_bvmethods_nmos_path_ionization.py"
SPEC = importlib.util.spec_from_file_location("compare_bvmethods_paths", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class CompareBVMethodsPathIonizationTests(unittest.TestCase):
    def test_final_write_all_paths_use_last_inventory_arithmetic_mean_and_alias_count(
        self,
    ) -> None:
        with tempfile.TemporaryDirectory() as directory:
            log = Path(directory) / "n4_des.log"
            log.write_text(
                """
Path number 0
Maximum Field: 1.0e+06
Electron: 0.2
Hole:     0.4
Best Path
Electron: 0.2
Hole:     0.4
Path number 0
Maximum Field: 1.5e+06
Electron: 1.2
Hole:     1.6
Path number 1
Maximum Field: 1.5e+06
Electron: 1.2
Hole:     1.6
Path number 2
Maximum Field: 1.7e+06
Electron: 1.7
Hole:     1.9
Best Path
Electron: 1.7
Hole:     1.9
""".strip()
                + "\n",
                encoding="utf-8",
            )

            rows = MODULE.final_write_all_path_integrals(log)

        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]["path_number"], 2)
        self.assertAlmostEqual(float(rows[0]["mean_ionization_integral"]), 1.8)
        self.assertEqual(rows[0]["multiplicity"], 1)
        self.assertEqual(rows[1]["path_number"], 0)
        self.assertAlmostEqual(float(rows[1]["mean_ionization_integral"]), 1.4)
        self.assertEqual(rows[1]["multiplicity"], 2)


if __name__ == "__main__":
    unittest.main()
