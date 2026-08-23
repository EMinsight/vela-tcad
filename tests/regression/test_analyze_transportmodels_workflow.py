#!/usr/bin/env python3
"""Regression coverage for TransportModels bias-regime analysis."""

from __future__ import annotations

import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "analyze_transportmodels_workflow.py"
SPEC = importlib.util.spec_from_file_location("transportmodels_analysis", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
ANALYSIS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ANALYSIS)


def write_curve(path: Path, column: str,
                rows: list[tuple[float, float]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["bias_V", column])
        writer.writeheader()
        for bias, current in rows:
            writer.writerow({"bias_V": bias, column: current})


class TransportModelsAnalysisTest(unittest.TestCase):
    def test_analyze_splits_idvg_regimes_and_excludes_zero_idvd_bias(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            generated = root / "generated"
            workflow = root / "workflow"
            workflow.mkdir()
            (workflow / "workflow_manifest.json").write_text(json.dumps({
                "status": "pass",
                "comparison_status": "pass",
                "controlled_delta": {"shared_mesh_doping_materials": True},
                "stages": [{"name": "stage", "status": "pass"}],
            }))
            idvg_biases = [-1.0, -0.84, -0.68, -0.52, -0.36,
                           -0.2, -0.04, 0.12, 0.28]
            idvd_biases = [0.0, 0.1, 0.2]
            for kind, biases in (("idvg", idvg_biases),
                                 ("idvd", idvd_biases)):
                reference = [(bias, float(index + 1))
                             for index, bias in enumerate(biases)]
                candidate = [(bias, 2.0 * current)
                             for bias, current in reference]
                write_curve(
                    generated / "reference_curves" /
                    f"transportmodels_sentaurus2022_dd_{kind}_reference.csv",
                    "current_total", reference)
                write_curve(
                    workflow / f"dd_{kind}_curve_comparison_candidate.csv",
                    "current_total_A_per_um", candidate)

            result = ANALYSIS.analyze(workflow, generated, "dd")
            regions = result["curves"]["idvg"]["regions"]
            self.assertEqual(3, regions["off"]["points"])
            self.assertEqual(5, regions["transition"]["points"])
            self.assertEqual(1, regions["on"]["points"])
            self.assertEqual(
                2, result["curves"]["idvd"]["nonzero_bias"]["points"])
            self.assertAlmostEqual(
                1.0, result["curves"]["idvd"]["nonzero_bias"][
                    "max_relative_error"])
            self.assertAlmostEqual(
                2.0, result["curves"]["idvg"]["endpoint"]["vela_A_per_um"] /
                result["curves"]["idvg"]["endpoint"]["reference_A_per_um"])

    def test_aligned_errors_rejects_bias_lattice_drift(self) -> None:
        with self.assertRaisesRegex(ValueError, "bias lattices differ"):
            ANALYSIS.aligned_errors(
                [{"bias_V": 0.0, "current": 1.0}],
                [{"bias_V": 0.1, "current": 1.0}],
            )

    def test_analyze_accepts_split_candidate_paths(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            generated = root / "generated"
            workflow = root / "workflow"
            workflow.mkdir()
            (workflow / "workflow_manifest.json").write_text(json.dumps({
                "status": "pass",
                "comparison_status": "pass",
                "controlled_delta": {},
                "stages": [],
            }), encoding="utf-8")
            idvg_biases = [-1.0, -0.84, -0.68, -0.52, -0.36,
                           -0.2, -0.04, 0.12, 0.28]
            idvd_biases = [0.0, 0.1, 0.2]
            candidate_paths = {}
            for kind, biases in (("idvg", idvg_biases),
                                 ("idvd", idvd_biases)):
                reference = [(bias, float(index + 1))
                             for index, bias in enumerate(biases)]
                write_curve(
                    generated / "reference_curves" /
                    f"transportmodels_sentaurus2022_dg_{kind}_reference.csv",
                    "current_total", reference)
                candidate_path = root / "split" / f"{kind}.csv"
                write_curve(candidate_path, "current_total_A_per_um", reference)
                candidate_paths[kind] = candidate_path

            result = ANALYSIS.analyze(
                workflow, generated, "dg", candidate_paths)
            self.assertEqual(
                str(candidate_paths["idvg"].resolve()),
                result["curves"]["idvg"]["candidate_path"],
            )
            self.assertEqual(
                str(candidate_paths["idvd"].resolve()),
                result["curves"]["idvd"]["candidate_path"],
            )


if __name__ == "__main__":
    unittest.main()
