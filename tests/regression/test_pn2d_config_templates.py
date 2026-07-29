#!/usr/bin/env python3
"""Regression tests for versioned PN2D simulation templates."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path, PureWindowsPath
from typing import Any


REPO = Path(__file__).resolve().parents[2]
if str(REPO / "scripts") not in sys.path:
    sys.path.insert(0, str(REPO / "scripts"))

from generate_pn2d_config import (  # noqa: E402
    TemplateError,
    render_named_template,
    write_rendered_config,
)


def all_strings(value: Any):
    if isinstance(value, dict):
        for item in value.values():
            yield from all_strings(item)
    elif isinstance(value, list):
        for item in value:
            yield from all_strings(item)
    elif isinstance(value, str):
        yield value


class Pn2dConfigTemplatesTest(unittest.TestCase):
    def test_iv_defaults_capture_qualified_forward_models(self) -> None:
        config, manifest = render_named_template("pn2d_iv")
        self.assertEqual(config["sweep"]["mode"], "iv")
        self.assertEqual(config["sweep"]["initial_step"], 1.0e-3)
        self.assertEqual(config["sweep"]["min_step"], 1.0e-8)
        self.assertEqual(config["sweep"]["max_step"], 0.025)
        self.assertEqual(config["sweep"]["growth_factor"], 1.3)
        self.assertEqual(config["solver"]["mobility"]["model"], "masetti")
        self.assertEqual(
            config["solver"]["mobility"]["doping_concentration_basis"],
            "cell_reconstructed_total_impurity",
        )
        self.assertEqual(config["solver"]["impact_ionization"]["model"], "none")
        self.assertEqual(manifest["template_version"], 1)

    def test_bv_defaults_capture_sentaurus_aligned_models(self) -> None:
        config, _ = render_named_template("pn2d_bv")
        self.assertEqual(config["sweep"]["mode"], "bv_reverse")
        self.assertEqual(config["sweep"]["initial_step"], 1.0e-4)
        self.assertEqual(config["sweep"]["min_step"], 1.0e-10)
        self.assertEqual(config["sweep"]["max_step"], 0.05)
        self.assertEqual(config["sweep"]["growth_factor"], 1.2)
        self.assertEqual(config["solver"]["max_iter"], 80)
        self.assertEqual(config["solver"]["mobility"]["model"], "masetti_field")
        self.assertEqual(
            config["solver"]["mobility"]["doping_concentration_basis"],
            "net_doping",
        )
        self.assertEqual(
            config["solver"]["impact_ionization"]["model"], "van_overstraeten"
        )
        self.assertEqual(
            config["sweep"]["write_state_every_point_prefix"],
            "pn2d_bv_states/state",
        )
        newton_history = config["sweep"]["diagnostics"]["newton_history"]
        self.assertEqual(
            newton_history["attempts_csv_file"],
            "pn2d_bv_newton_attempts.csv",
        )
        self.assertEqual(
            newton_history["iterations_csv_file"],
            "pn2d_bv_newton_iterations.csv",
        )

    def test_defaults_contain_no_absolute_paths_or_placeholders(self) -> None:
        for template in ("pn2d_iv", "pn2d_bv"):
            config, _ = render_named_template(template)
            for value in all_strings(config):
                self.assertNotIn("${", value)
                self.assertFalse(Path(value).is_absolute(), value)
                self.assertFalse(PureWindowsPath(value).is_absolute(), value)

    def test_unknown_and_invalid_overrides_are_rejected(self) -> None:
        with self.assertRaisesRegex(TemplateError, "unknown"):
            render_named_template("pn2d_iv", {"typo": 1})
        with self.assertRaisesRegex(TemplateError, "type number"):
            render_named_template("pn2d_iv", {"stop_voltage": "20"})
        with self.assertRaisesRegex(TemplateError, "must be relative"):
            render_named_template(
                "pn2d_iv", {"mesh_file": r"D:\external\mesh.json"}
            )
        with self.assertRaisesRegex(TemplateError, "positive sweep.step"):
            render_named_template("pn2d_bv", {"stop_voltage": 20.0})

    def test_rendering_is_byte_deterministic_and_manifest_is_separate(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_pn2d_template_") as td:
            root = Path(td)
            first = root / "first.json"
            second = root / "second.json"
            write_rendered_config("pn2d_iv", first)
            write_rendered_config("pn2d_iv", second)
            self.assertEqual(first.read_bytes(), second.read_bytes())
            first_manifest = json.loads(
                first.with_suffix(".manifest.json").read_text(encoding="utf-8")
            )
            rendered = json.loads(first.read_text(encoding="utf-8"))
            self.assertEqual(first_manifest["template"], "pn2d_iv")
            self.assertEqual(first_manifest["overrides"], {})
            self.assertNotIn("template_schema", rendered)

    def test_cli_writes_valid_bv_config_and_manifest(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_pn2d_template_cli_") as td:
            output = Path(td) / "simulation.json"
            result = subprocess.run(
                [
                    sys.executable,
                    str(REPO / "scripts" / "generate_pn2d_config.py"),
                    "--template",
                    "pn2d_bv",
                    "--output",
                    str(output),
                    "--set",
                    "write_vtk=true",
                ],
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(result.returncode, 0, result.stderr)
            config = json.loads(output.read_text(encoding="utf-8"))
            self.assertTrue(config["sweep"]["write_vtk"])
            self.assertTrue(output.with_suffix(".manifest.json").is_file())

    def test_json_schema_artifact_is_draft_2020_12(self) -> None:
        schema = json.loads(
            (REPO / "configs" / "schema" / "vela-simulation.schema.json")
            .read_text(encoding="utf-8")
        )
        self.assertEqual(
            schema["$schema"], "https://json-schema.org/draft/2020-12/schema"
        )
        self.assertIn("sweep", schema["required"])
        self.assertIn("initial_step", schema["properties"]["sweep"]["required"])
        sweep_properties = schema["properties"]["sweep"]["properties"]
        self.assertIn("write_state_every_point_prefix", sweep_properties)
        newton_history = sweep_properties["diagnostics"]["properties"]["newton_history"]
        self.assertIn("attempts_csv_file", newton_history["properties"])
        self.assertIn("iterations_csv_file", newton_history["properties"])


if __name__ == "__main__":
    unittest.main()
