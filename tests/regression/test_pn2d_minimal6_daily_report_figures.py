import hashlib
import tempfile
import unittest
from unittest import mock
from pathlib import Path

from PIL import Image

import scripts.compare_pn2d_minimal6_diagnostic_sweeps as comparison_module
from tests.regression.test_pn2d_minimal6_sweep_comparison import (
    checkpoint,
    manifest,
    write_comparison_package,
)
from scripts.render_pn2d_minimal6_daily_report_figures import (
    DAILY_FIGURE_NAMES,
    render_daily_report_figures,
    verify_daily_report_figures,
)


class DailyReportFigureTest(unittest.TestCase):
    def setUp(self):
        strict_validator = mock.patch.object(
            comparison_module, "validate_sweep_manifest", create=True
        )
        strict_validator.start()
        self.addCleanup(strict_validator.stop)

    def comparison_package(self, root: Path) -> Path:
        biases = (-1.0, -2.0, -20.0)
        rows = {
            solver: [
                checkpoint(
                    solver,
                    topology,
                    bias,
                    1.0 + abs(bias),
                    -(1.0 + abs(bias)),
                    2.0 + abs(bias),
                    3.0 + abs(bias),
                    4.0 + abs(bias),
                )
                for topology in ("sketch", "mirror")
                for bias in biases
            ]
            for solver in ("vela", "sentaurus")
        }
        write_comparison_package(
            root,
            manifest("vela", rows["vela"]),
            manifest("sentaurus", rows["sentaurus"]),
            fixed_state_report={},
        )
        return root / "sweep_comparison.json"

    def test_contract_has_only_daily_figures_and_descending_voltage_axis(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            source = self.comparison_package(root / "source")
            output = root / "daily"

            contract = render_daily_report_figures(source, output)

            self.assertEqual(set(contract["figures"]), set(DAILY_FIGURE_NAMES))
            self.assertEqual(
                set(DAILY_FIGURE_NAMES),
                {"terminal_current.png", "maximum_field.png", "source_integrals.png"},
            )
            for name, entry in contract["figures"].items():
                self.assertEqual(entry["x_quantity"], "applied_bias_V")
                self.assertEqual(
                    entry["x_axis_order"], "decreasing_left_to_right"
                )
                self.assertEqual(entry["x_limits_V"], [-1.0, -20.0])
                self.assertGreater(entry["x_limits_V"][0], entry["x_limits_V"][1])
                with Image.open(output / name) as image:
                    self.assertEqual(image.size, (900, 504))
                    self.assertEqual(
                        image.info["VoltageAxisOrder"],
                        "decreasing_left_to_right",
                    )
            self.assertTrue(
                verify_daily_report_figures(
                    output / "daily_report_figure_manifest.json"
                )
            )

    def test_rerender_is_byte_deterministic_and_source_is_unchanged(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            source = self.comparison_package(root / "source")
            before = hashlib.sha256(source.read_bytes()).hexdigest()

            render_daily_report_figures(source, root / "first")
            render_daily_report_figures(source, root / "second")

            self.assertEqual(before, hashlib.sha256(source.read_bytes()).hexdigest())
            for name in (*DAILY_FIGURE_NAMES, "daily_report_figure_manifest.json"):
                self.assertEqual(
                    (root / "first" / name).read_bytes(),
                    (root / "second" / name).read_bytes(),
                )


if __name__ == "__main__":
    unittest.main()
