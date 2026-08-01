from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from scripts.generate_pn2d_config import render_named_template
from scripts.review_pn2d_node_volume_policy_default_acceptance import (
    default_render_binding,
    on_gate,
    scope_and_rollback_gates,
    sha256,
)


class NodeVolumePolicyDefaultAcceptanceTests(unittest.TestCase):
    def test_no_crossing_is_a_typed_pair_outcome(self) -> None:
        metrics = {
            "all_nonzero": {"rmse_dex": 0.001, "maximum_dex": 0.002},
            "knee": {"rmse_dex": 0.001, "maximum_dex": 0.002},
            "V_break_V": {"vela": -19.4, "sentaurus": -19.41},
            "V_slope_V": {"vela": None, "sentaurus": None},
            "nonmonotonic_intervals_V": [],
        }
        gate = {
            "maximum_all_nonzero_log10_current_rmse_dex": 0.01,
            "maximum_all_nonzero_log10_current_error_dex": 0.03,
            "maximum_knee_log10_current_rmse_dex": 0.01,
            "maximum_knee_log10_current_error_dex": 0.03,
            "maximum_V_break_absolute_error_V": 0.10,
        }
        self.assertTrue(all(on_gate(metrics, gate).values()))

    def test_default_render_binding_replays_manifest_and_rejects_override(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            template = root / "pn2d_bv.template.json"
            source_template = Path("configs/templates/pn2d_bv.template.json")
            template.write_bytes(source_template.read_bytes())
            config, manifest = render_named_template("pn2d_bv")
            config_path = root / "simulation.json"
            manifest_path = root / "simulation.manifest.json"
            config_path.write_text(json.dumps(config))
            manifest_path.write_text(json.dumps(manifest))
            controls = {
                "candidate_origin": "pn2d_bv_template_default_render",
                "candidate_config": {
                    "path": str(config_path), "sha256": sha256(config_path)
                },
                "base_config_manifest": {
                    "path": str(manifest_path), "sha256": sha256(manifest_path)
                },
                "default_render_binding": {"gates": {"atomic": True}},
            }
            contract = {
                "artifact_bindings": {
                    "pn2d_bv_template_sha256": sha256(template)
                },
                "candidate": {"required_template_version_minimum": 3},
            }
            self.assertTrue(default_render_binding(controls, template, contract)["passed"])

            manifest["overrides"] = {
                "avalanche_current_support_profile": "element_edge_sg_gss_laux"
            }
            manifest_path.write_text(json.dumps(manifest))
            controls["base_config_manifest"]["sha256"] = sha256(manifest_path)
            result = default_render_binding(controls, template, contract)
            self.assertFalse(result["passed"])
            self.assertFalse(result["gates"]["manifest_has_no_profile_override"])

    def test_scope_and_rollback_gates_bind_bv_and_leave_iv_unchanged(self) -> None:
        bv_template = Path("configs/templates/pn2d_bv.template.json")
        iv_template = Path("configs/templates/pn2d_iv.template.json")
        contract = {
            "artifact_bindings": {
                "pn2d_bv_template_sha256": sha256(bv_template),
                "pn2d_iv_template_sha256": sha256(iv_template),
            },
            "candidate": {
                "impact_ionization": {
                    "current_approximation": "element_edge_sg_gss_laux",
                    "source_mapping_mode": "element_vertex_box_measure",
                    "cell_reconstructed_midpoint_density": "bernoulli",
                },
                "mesh_geometry": {
                    "node_volume_policy": "mixed_voronoi",
                    "require_non_obtuse": True,
                },
            },
            "rollback": {
                "impact_ionization": {
                    "current_approximation": "cell_reconstructed",
                    "source_mapping_mode": "triangle_gss_gradqf_truncated",
                    "cell_reconstructed_midpoint_density": "gss_logistic",
                },
                "mesh_geometry": {
                    "node_volume_policy": "barycentric",
                    "require_non_obtuse": False,
                },
            },
        }
        result = scope_and_rollback_gates(bv_template, iv_template, contract)
        self.assertTrue(result["passed"], result["gates"])


if __name__ == "__main__":
    unittest.main()
