from __future__ import annotations

import csv
import json
import tempfile
import unittest
from pathlib import Path

from scripts.run_pn2d_bv_m2_mixed_voronoi_self_consistent_control import (
    ATOMIC_DEFAULT_PROFILE,
    determinism_metrics,
    off_golden_metrics,
    prepare_mixed_config,
    validate_actual_default_render,
)


class MixedVoronoiSelfConsistentControlTests(unittest.TestCase):
    def test_prepare_mixed_config_is_opt_in_and_does_not_mutate_base(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            base = root / "base.json"
            output = root / "mixed.json"
            base.write_text(json.dumps({"mesh_geometry": {"node_volume_policy": "barycentric"}}))
            prepared = prepare_mixed_config(base, output)
            original = json.loads(base.read_text())
        self.assertEqual(original["mesh_geometry"]["node_volume_policy"], "barycentric")
        self.assertEqual(prepared["mesh_geometry"]["node_volume_policy"], "mixed_voronoi")

    def test_off_golden_metrics_uses_nonzero_native_total_current(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            iv = root / "iv.csv"
            with iv.open("w", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=["bias_V", "current_total_A_per_um"])
                writer.writeheader()
                writer.writerows([
                    {"bias_V": 0.0, "current_total_A_per_um": 0.0},
                    {"bias_V": -1.0, "current_total_A_per_um": -1.0e-12},
                    {"bias_V": -2.0, "current_total_A_per_um": -2.0e-12},
                ])
            manifest = root / "sentaurus.json"
            manifest.write_text(json.dumps({"aggregate_records": [
                {"branch": "avalanche_off", "quantity": "terminal_current", "carrier": "total", "provenance": "native", "requested_bias_V": -1.0, "value": -1.0e-12},
                {"branch": "avalanche_off", "quantity": "terminal_current", "carrier": "total", "provenance": "native", "requested_bias_V": -2.0, "value": -2.0e-12},
            ]}))
            metrics = off_golden_metrics(iv, manifest)
        self.assertEqual(metrics["compared_bias_count"], 2)
        self.assertEqual(metrics["log10_current_rmse_dex"], 0.0)

    def test_determinism_requires_iv_and_every_state_hash(self) -> None:
        def execution(iv_hash: str, state_hash: str) -> dict:
            return {"branches": [{
                "branch": "avalanche_off",
                "output_csv_sha256": iv_hash,
                "state_files": {"-1": {"sha256": state_hash}},
            }]}
        equal = determinism_metrics(execution("iv", "state"), execution("iv", "state"), "avalanche_off")
        changed = determinism_metrics(execution("iv", "state"), execution("iv", "other"), "avalanche_off")
        self.assertTrue(equal["iv_sha256_equal"])
        self.assertTrue(equal["state_hashes_equal"])
        self.assertFalse(changed["state_hashes_equal"])

    def test_actual_default_render_requires_atomic_profile_without_override(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            config = root / "simulation.json"
            manifest = root / "simulation.manifest.json"
            config.write_text(json.dumps({
                "solver": {"impact_ionization": ATOMIC_DEFAULT_PROFILE["impact_ionization"]},
                "mesh_geometry": ATOMIC_DEFAULT_PROFILE["mesh_geometry"],
            }))
            manifest.write_text(json.dumps({
                "template": "pn2d_bv",
                "template_version": 3,
                "overrides": {},
                "parameters": {
                    "avalanche_current_support_profile": "element_edge_sg_gss_laux"
                },
                "resolved_profile": ATOMIC_DEFAULT_PROFILE,
            }))
            result = validate_actual_default_render(config, manifest)
            self.assertTrue(all(result["gates"].values()))

            mutated = json.loads(config.read_text())
            mutated["mesh_geometry"]["node_volume_policy"] = "barycentric"
            config.write_text(json.dumps(mutated))
            with self.assertRaisesRegex(ValueError, "atomic binding"):
                validate_actual_default_render(config, manifest)

    def test_actual_default_render_rejects_hidden_profile_override(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            config = root / "simulation.json"
            manifest = root / "simulation.manifest.json"
            config.write_text(json.dumps({
                "solver": {"impact_ionization": ATOMIC_DEFAULT_PROFILE["impact_ionization"]},
                "mesh_geometry": ATOMIC_DEFAULT_PROFILE["mesh_geometry"],
            }))
            manifest.write_text(json.dumps({
                "template": "pn2d_bv",
                "template_version": 3,
                "overrides": {
                    "avalanche_current_support_profile": "element_edge_sg_gss_laux"
                },
                "parameters": {
                    "avalanche_current_support_profile": "element_edge_sg_gss_laux"
                },
                "resolved_profile": ATOMIC_DEFAULT_PROFILE,
            }))
            with self.assertRaisesRegex(ValueError, "profile_override"):
                validate_actual_default_render(config, manifest)


if __name__ == "__main__":
    unittest.main()
