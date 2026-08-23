from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

from scripts.prepare_slot_ldmos_ialmob_ablation import (
    CONTROL_DECKS,
    build_no_ialmob_deck,
    prepare,
)
from scripts.analyze_slot_ldmos_ialmob_ablation import relative_delta


TEMPLATE = """File {{
  Output= \"{stem}.log\"
  Current= \"{stem}.plt\"
}}
Physics (Material=\"Silicon\") {{
  Mobility(
    Enormal (IALMob)
    HighFieldSaturation
  )
}}
"""


class SlotLdmosIalMobAblationTest(unittest.TestCase):
    def test_relative_delta_preserves_control_direction(self) -> None:
        self.assertAlmostEqual(relative_delta(38.7, 38.5), 0.2 / 38.5)

    def test_build_removes_only_selector_and_isolates_outputs(self) -> None:
        result = build_no_ialmob_deck(TEMPLATE.format(stem="case"), "case")
        self.assertNotIn("IALMob", result)
        self.assertNotIn("Enormal", result)
        self.assertIn("HighFieldSaturation", result)
        self.assertIn("case_no_ialmob.plt", result)

    def test_prepare_writes_two_controls_and_manifest(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source"
            output = root / "output"
            source.mkdir()
            for stem in CONTROL_DECKS:
                (source / f"{stem}.cmd").write_text(
                    TEMPLATE.format(stem=stem), encoding="utf-8"
                )
            (source / "pp2_des.par").write_text("parameter data\n", encoding="utf-8")

            manifest = prepare(source, output)

            self.assertEqual(len(manifest["cases"]), 2)
            self.assertTrue((output / "ialmob_ablation_manifest.json").is_file())
            for stem in CONTROL_DECKS:
                deck = output / f"{stem}_no_ialmob.cmd"
                self.assertTrue(deck.is_file())
                self.assertNotIn("IALMob", deck.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
