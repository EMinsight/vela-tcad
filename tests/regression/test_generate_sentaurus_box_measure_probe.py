from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from scripts.generate_sentaurus_box_measure_probe import main, sha256


class SentaurusBoxMeasureProbeTests(unittest.TestCase):
    def test_generates_explicit_mix_average_probe(self) -> None:
        with tempfile.TemporaryDirectory() as raw_tmp:
            root = Path(raw_tmp)
            source = root / "source"
            output = root / "probe"
            source.mkdir()
            for name in (
                "models.par",
                "pn2d_msh.tdr",
                "runtime_element_avalanche_probe.tcl",
            ):
                (source / name).write_text(name, encoding="utf-8")
            (source / "input.cmd").write_text(
                "Plot {\n  Potential\n}\n"
                "CurrentPlot {\n}\n"
                "Math {\n}\n"
                "Solve {\n  Coupled { Poisson Electron Hole }\n}\n",
                encoding="utf-8",
            )
            argv = [
                "generate_sentaurus_box_measure_probe.py",
                "--source-dir",
                str(source),
                "--output-dir",
                str(output),
                "--source-deck",
                "input.cmd",
                "--box-method",
                "mix-average",
            ]
            with patch("sys.argv", argv):
                self.assertEqual(main(), 0)

            deck = (output / "box_measure_probe.cmd").read_text(encoding="utf-8")
            self.assertIn("MixAverageBoxMethod", deck)
            self.assertIn("BoxMeasureFromFile(GrdNumbering)", deck)
            self.assertIn("BM_CoeffIntersectionNonDelaunayElements", deck)
            self.assertIn('Plot(FilePrefix="box_measure_probe")', deck)
            self.assertNotIn("Electron Hole", deck)
            self.assertEqual(len(sha256(output / "box_measure_probe.cmd")), 64)


if __name__ == "__main__":
    unittest.main()
