from __future__ import annotations

import importlib.util
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "analyze_singledevice_closeout_fields.py"
SPEC = importlib.util.spec_from_file_location("singledevice_fields", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class SingleDeviceFieldAuditTest(unittest.TestCase):
    def test_percentile_and_field_summaries(self) -> None:
        self.assertEqual(2.5, MODULE.percentile([1.0, 2.0, 3.0, 4.0], 0.5))
        absolute = MODULE.abs_summary({0: 1.0, 1: 2.0}, {0: 1.1, 1: 1.8})
        self.assertAlmostEqual(0.15, absolute["median_abs_error"])
        logarithmic = MODULE.log_summary({0: 1.0, 1: 10.0}, {0: 2.0, 1: 10.0})
        self.assertEqual(2, logarithmic["count"])
        self.assertAlmostEqual(0.1505149978, logarithmic["median_abs_log10_error_dex"])


if __name__ == "__main__":
    unittest.main()
