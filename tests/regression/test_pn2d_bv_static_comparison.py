import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "plot_pn2d_bv_static_comparison.py"
SPEC = importlib.util.spec_from_file_location("plot_pn2d_bv_static_comparison", SCRIPT)
plot = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = plot
SPEC.loader.exec_module(plot)


class TestCurrentSemanticsFigure(unittest.TestCase):
    def test_loads_only_exact_groups_in_bias_then_carrier_order(self):
        groups = [
            self._group(-19.4, "electron", 0, None),
            self._group(-19.0, "hole", 42, 2.0),
            self._group(-12.0, "electron", 42, 1.0),
        ]
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "audit.json"
            path.write_text(json.dumps({"active_support": {"by_bias_and_carrier": groups}}))
            labels, medians, coverage = plot.load_current_semantics_groups(path)

        self.assertEqual(labels, ["-12 V / e", "-19 V / h"])
        self.assertEqual(medians["pdf"], [1.0, 2.0])
        self.assertEqual(coverage["pdf"], [0.75, 0.75])

    @staticmethod
    def _group(bias, carrier, exact_count, median):
        candidates = {}
        for name in ("pdf", "genius", "vela_projection", "vela_magnitude"):
            candidates[name] = {
                "coverage": 0.75,
                "median_abs_log10_error": median,
                "p95_abs_log10_error": None if median is None else median + 1.0,
            }
        return {
            "bias_V": bias,
            "carrier": carrier,
            "exact_row_count": exact_count,
            "candidates": candidates,
        }


if __name__ == "__main__":
    unittest.main()
