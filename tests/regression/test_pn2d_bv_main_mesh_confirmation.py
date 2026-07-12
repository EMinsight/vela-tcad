#!/usr/bin/env python3
"""Regression tests for PN2D BV main-mesh mechanism confirmation."""

from __future__ import annotations

import importlib.util
import json
import math
import tempfile
import unittest
from pathlib import Path
from unittest import mock


REPO = Path(__file__).resolve().parents[2]
MODULE_PATH = REPO / "scripts" / "diagnose_pn2d_bv_main_mesh_confirmation.py"


def load_module():
    if not MODULE_PATH.exists():
        raise AssertionError(f"missing main-mesh confirmation script: {MODULE_PATH}")
    spec = importlib.util.spec_from_file_location(
        "diagnose_pn2d_bv_main_mesh_confirmation_test",
        MODULE_PATH,
    )
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def raw_anchor_edges() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    replay_flux = 10.0 / (10.0 ** 0.05)
    for edge_id in range(101):
        vela_source = 0.0
        sentaurus_source = 0.0
        if edge_id == 0:
            vela_source = 100.0
            sentaurus_source = 100.0
        elif edge_id == 1:
            vela_source = 90.0
            sentaurus_source = 1.0
        elif edge_id == 2:
            vela_source = 1.0
            sentaurus_source = 90.0
        rows.append({
            "edge_id": edge_id,
            "node0": edge_id,
            "node1": edge_id + 1,
            "x0_um": float(edge_id),
            "y0_um": 0.0,
            "x1_um": float(edge_id) + 0.5,
            "y1_um": 0.0,
            "edge_type": "p-n",
            "vela_source_physical_m_inv_s": vela_source,
            "sentaurus_same_area_source_proxy_physical_m_inv_s": sentaurus_source,
            "vela_alpha_m_inv": 100.0,
            "sentaurus_alpha_same_edge_m_inv": 100.0,
            "vela_flux_abs_m2_s": 1.0,
            "sentaurus_vector_flux_abs_m2_s": 10.0,
            "sentaurus_replay_flux_abs_m2_s": replay_flux,
            "production_highprec_relative_error": 1.0e-8,
            "cancellation_condition": 10.0,
            "any_exponent_clamped": False,
        })
    return rows


class MainMeshConfirmationTest(unittest.TestCase):
    def test_p99_union_reports_overlap_false_sets_and_replay_recovery(self) -> None:
        module = load_module()
        rows, summary = module.analyze_anchor_rows(
            -19.0,
            raw_anchor_edges(),
            percentile_value=99.0,
        )

        self.assertEqual(len(rows), 101)
        self.assertEqual(summary["union_count"], 3)
        self.assertEqual(summary["overlap_count"], 1)
        self.assertEqual(summary["false_positive_count"], 1)
        self.assertEqual(summary["false_negative_count"], 1)
        self.assertAlmostEqual(summary["vela_source_p99_m_inv_s"], 90.0)
        self.assertAlmostEqual(summary["sentaurus_same_area_source_proxy_p99_m_inv_s"], 90.0)
        self.assertAlmostEqual(summary["median_gap_recovery"], 0.95, places=12)
        self.assertAlmostEqual(summary["median_alpha_gap_dex"], 0.0, places=12)
        self.assertEqual(
            summary["mechanism"]["classification"],
            "vela_internal_state_branch",
        )
        union_rows = [row for row in rows if row["in_active_union"]]
        self.assertEqual(
            {row["support_class"] for row in union_rows},
            {"overlap", "false_positive", "false_negative"},
        )

    def test_confirmation_gate_requires_four_of_five_and_high_bias_recovery(self) -> None:
        module = load_module()
        anchors = []
        for bias in module.DEFAULT_BIASES:
            _rows, summary = module.analyze_anchor_rows(
                bias,
                raw_anchor_edges(),
                percentile_value=99.0,
            )
            anchors.append(summary)
        anchors[0]["mechanism"] = {
            "classification": "sg_discretization_ni_or_current_semantics",
            "rule": "fixture minority",
            "evidence": {},
        }

        gate = module.evaluate_confirmation_gate(anchors)
        self.assertEqual(gate["status"], "pass")
        self.assertEqual(gate["dominant_mechanism"], "vela_internal_state_branch")
        self.assertEqual(gate["same_mechanism_count"], 4)
        self.assertTrue(gate["support_bidirectional_pass"])
        self.assertGreater(gate["support_false_positive_total"], 0)
        self.assertGreater(gate["support_false_negative_total"], 0)
        self.assertGreaterEqual(gate["high_bias_recovery"]["-19"], 0.8)
        self.assertGreaterEqual(gate["high_bias_recovery"]["-20"], 0.8)
        self.assertEqual(gate["next_target"], "main_mesh_continuation_branch_recovery")
        self.assertEqual(
            gate["minimum_failing_test"],
            "test_pn2d_bv_main_mesh_continuation_recovers_multiplication_current",
        )

        failing = json.loads(json.dumps(anchors))
        next(item for item in failing if item["bias_V"] == -20.0)[
            "median_gap_recovery"
        ] = 0.79
        failed_gate = module.evaluate_confirmation_gate(failing)
        self.assertEqual(failed_gate["status"], "fail")
        self.assertEqual(
            failed_gate["next_target"],
            "high_bias_sent_state_replay_recovery",
        )
        self.assertEqual(
            failed_gate["minimum_failing_test"],
            "test_pn2d_bv_main_mesh_confirmation_high_bias_recovery_gate",
        )

        one_sided = json.loads(json.dumps(anchors))
        for anchor in one_sided:
            anchor["false_positive_count"] = 0
        support_gate = module.evaluate_confirmation_gate(one_sided)
        self.assertEqual(support_gate["status"], "fail")
        self.assertFalse(support_gate["support_bidirectional_pass"])
        self.assertEqual(
            support_gate["next_target"],
            "main_mesh_bidirectional_support_explanation",
        )
        self.assertEqual(
            support_gate["minimum_failing_test"],
            "test_pn2d_bv_main_mesh_confirmation_requires_bidirectional_support",
        )

    def test_alpha_gap_prevents_false_branch_ownership(self) -> None:
        module = load_module()
        rows = raw_anchor_edges()
        for row in rows:
            row["vela_flux_abs_m2_s"] = 10.0
            row["sentaurus_vector_flux_abs_m2_s"] = 10.0
            row["sentaurus_replay_flux_abs_m2_s"] = 10.0 / (10.0 ** 0.15)
            row["sentaurus_alpha_same_edge_m_inv"] = 200.0

        analyzed, summary = module.analyze_anchor_rows(-19.0, rows)

        self.assertAlmostEqual(summary["median_alpha_gap_dex"], math.log10(2.0))
        self.assertEqual(
            summary["mechanism"]["classification"],
            "impact_coefficient_or_source_semantics",
        )
        union = [row for row in analyzed if row["in_active_union"]]
        self.assertTrue(all(row["alpha_gap_dex"] > 0.1 for row in union))

    def test_anchor_mechanism_does_not_combine_extrema_from_different_edges(self) -> None:
        module = load_module()
        rows = raw_anchor_edges()
        for row in rows:
            row["sentaurus_replay_flux_abs_m2_s"] = 10.0 / (10.0 ** 0.15)
        rows[0]["production_highprec_relative_error"] = 1.0e-3
        rows[0]["cancellation_condition"] = 1.0
        rows[1]["production_highprec_relative_error"] = 1.0e-9
        rows[1]["cancellation_condition"] = 1.0e13
        rows[1]["any_exponent_clamped"] = True

        analyzed, summary = module.analyze_anchor_rows(-19.0, rows)

        union = [row for row in analyzed if row["in_active_union"]]
        self.assertNotIn(
            "variable_ni_sg_numerical_stability",
            {row["row_mechanism_classification"] for row in union},
        )
        self.assertEqual(summary["mechanism"]["classification"], "inconclusive")

    def test_vtk_resolver_uses_encoded_bias_not_nominal_step_index(self) -> None:
        module = load_module()
        with tempfile.TemporaryDirectory(prefix="vela_main_mesh_vtk_") as td:
            root = Path(td)
            expected = root / "dc_sweep_0702_-10V.vtk"
            expected.write_text("fixture", encoding="utf-8")
            (root / "dc_sweep_0200_-2.85661V.vtk").write_text("fixture", encoding="utf-8")

            resolved = module.resolve_vtk_for_bias(root, "dc_sweep", -10.0)

            self.assertEqual(resolved, expected)

    def test_missing_vector_current_and_alpha_exports_fail_explicit_gate(self) -> None:
        module = load_module()
        with tempfile.TemporaryDirectory(prefix="vela_main_mesh_artifact_") as td:
            root = Path(td)
            export_dir = root / "sentaurus_-10v"
            export_dir.mkdir()
            (export_dir / "field_manifest.json").write_text(json.dumps({
                "fields": [{
                    "name": "eCurrentDensity",
                    "components": 1,
                    "unit": "A*cm^-2",
                    "region": 0,
                }],
            }), encoding="utf-8")
            out_dir = root / "out"

            return_code = module.main([
                "--sg-csv", str(root / "sg.csv"),
                "--vtk-root", str(root / "vtk"),
                "--imported-doping", str(root / "doping.csv"),
                f"--sentaurus-export=-10={export_dir}",
                "--biases=-10",
                "--out-dir", str(out_dir),
            ])

            self.assertEqual(return_code, 0)
            payload = json.loads(
                (out_dir / "main_mesh_confirmation_summary.json").read_text(encoding="utf-8")
            )
            self.assertEqual(payload["gate"]["status"], "fail")
            self.assertFalse(payload["gate"]["artifact_contract_pass"])
            self.assertEqual(
                payload["gate"]["next_target"],
                "sentaurus_main_mesh_vector_current_alpha_export",
            )
            reasons = " ".join(payload["artifact_checks"][0]["reasons"])
            self.assertIn("eCurrentDensity components=2", reasons)
            self.assertIn("eAlphaAvalanche raw export required", reasons)
            report = (out_dir / "main_mesh_confirmation_report.md").read_text(
                encoding="utf-8"
            )
            self.assertIn("Scalar current magnitude is never substituted", report)

    def test_collector_uses_real_enriched_source_schema(self) -> None:
        module = load_module()
        args = module.parse_args([
            "--sg-csv", "sg.csv",
            "--vtk-root", "vtk",
            "--imported-doping", "doping.csv",
            "--sentaurus-root", "sentaurus",
            "--out-dir", "out",
        ])
        edge = {
            "bias_V": "-19",
            "edge_id": "7",
            "node0": "0",
            "node1": "1",
            "x0_um": "0",
            "y0_um": "0",
            "x1_um": "1",
            "y1_um": "0",
            "electron_sg_node0_exponent_clamped_low": "0",
            "electron_sg_node0_exponent_clamped_high": "0",
            "electron_sg_node1_exponent_clamped_low": "0",
            "electron_sg_node1_exponent_clamped_high": "0",
            "electron_sg_production_vs_high_precision_reference_relative_error": "1e-8",
            "electron_sg_cancellation_condition": "10",
            "electron_alpha_m_inv": "5",
        }
        nodes = [
            {"id": 0, "x_um": 0.0, "y_um": 0.0},
            {"id": 1, "x_um": 1.0, "y_um": 0.0},
        ]
        enriched = {
            "vela_e_source_integral_physical_m_inv_s": 1.0,
            "sentaurus_e_source_on_vela_area_physical_m_inv_s": 2.0,
            "sentaurus_e_alpha_edge_average_m_inv": 10.0,
            "vela_e_sg_production_canonical_signed_flux_m2_s": 3.0,
            "sentaurus_e_continuity_edge_signed_flux_m2_s": 4.0,
            "sentaurus_e_sg_vela_mobility_signed_flux_m2_s": 3.8,
        }
        with (
            mock.patch.object(module.compensated, "load_doping", return_value={
                0: {"type": "p"},
                1: {"type": "n"},
            }),
            mock.patch.object(module, "resolve_vtk_for_bias", return_value=Path("state.vtk")),
            mock.patch.object(module.compensated, "parse_vtk", return_value={
                "points": [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)],
                "scalars": {},
            }),
            mock.patch.object(module, "load_sentaurus_nodes_from_export", return_value=nodes),
            mock.patch.object(module.compensated, "load_sentaurus_electron_state", return_value={}),
            mock.patch.object(module.compensated, "read_csv", return_value=[edge]),
            mock.patch.object(
                module.compensated,
                "enrich_edge_with_sentaurus_replay",
                return_value=enriched,
            ),
        ):
            rows = module.collect_anchor_edges(args, -19.0)
        self.assertEqual(
            rows[0]["sentaurus_same_area_source_proxy_physical_m_inv_s"],
            2.0,
        )
        self.assertEqual(rows[0]["vela_alpha_m_inv"], 5.0)
        self.assertEqual(rows[0]["sentaurus_alpha_same_edge_m_inv"], 10.0)

    def test_cli_writes_current_head_diagnostic_deck_under_out_dir(self) -> None:
        module = load_module()
        with tempfile.TemporaryDirectory(prefix="vela_main_mesh_deck_") as td:
            root = Path(td)
            template_path = root / "simulation_dense_full_floor_light.json"
            out_dir = root / "out"
            template_path.write_text(json.dumps({
                "_comment": "fixture dense main deck",
                "simulation_type": "dc_sweep",
                "output_csv": str(root / "old.csv"),
                "solver": {
                    "impact_ionization": {
                        "generation": "current_density",
                        "current_approximation": "density_gradient",
                    },
                },
                "sweep": {
                    "mode": "bv_reverse",
                    "start": 0.0,
                    "stop": -20.0,
                    "step": -0.1,
                    "bias_points": [0.0, -10.0, -13.2, -19.0, -20.0],
                    "write_vtk": False,
                    "vtk_prefix": str(root / "old_vtk"),
                    "write_state_file": str(root / "old_state.csv"),
                },
            }), encoding="utf-8")

            return_code = module.main([
                "--write-diagnostic-deck-from", str(template_path),
                "--out-dir", str(out_dir),
            ])

            self.assertEqual(return_code, 0)
            deck_path = out_dir / "main_mesh_confirmation_diagnostic_deck.json"
            deck = json.loads(deck_path.read_text(encoding="utf-8"))
            self.assertTrue(set(module.DEFAULT_BIASES).issubset(deck["sweep"]["bias_points"]))
            self.assertTrue(deck["sweep"]["write_vtk"])
            self.assertEqual(
                deck["solver"]["impact_ionization"]["generation"],
                "current_density",
            )
            self.assertEqual(
                deck["solver"]["impact_ionization"]["current_approximation"],
                "density_gradient",
            )
            self.assertEqual(
                Path(deck["output_csv"]).resolve().parent,
                out_dir.resolve(),
            )
            self.assertEqual(
                Path(deck["sweep"]["vtk_prefix"]).resolve().parent,
                (out_dir / "vtk").resolve(),
            )
            self.assertEqual(
                Path(deck["sweep"]["diagnostics"]["sg_avalanche_edges"]["csv_file"]).resolve(),
                (out_dir / "sg_avalanche_edges.csv").resolve(),
            )
            self.assertIn("current HEAD", deck["_comment"])

            bad_deck = json.loads(template_path.read_text(encoding="utf-8"))
            bad_deck["solver"]["impact_ionization"]["current_approximation"] = "edge_scalar"
            bad_path = root / "bad_template.json"
            bad_path.write_text(json.dumps(bad_deck), encoding="utf-8")
            with self.assertRaisesRegex(
                ValueError, "current_density.*density_gradient"
            ):

                module.write_diagnostic_deck(bad_path, root / "bad_out")

    def test_cli_writes_csv_json_md_with_unique_followup(self) -> None:
        module = load_module()
        with tempfile.TemporaryDirectory(prefix="vela_main_mesh_confirmation_") as td:
            root = Path(td)
            out_dir = root / "out"
            with (
                mock.patch.object(
                    module,
                    "collect_anchor_edges",
                    side_effect=lambda _args, _bias: raw_anchor_edges(),
                ),
                mock.patch.object(module, "preflight_sentaurus_exports", return_value=[]),
            ):
                return_code = module.main([
                    "--sg-csv", str(root / "sg.csv"),
                    "--vtk-root", str(root / "vtk"),
                    "--imported-doping", str(root / "doping.csv"),
                    "--sentaurus-root", str(root / "sentaurus"),
                    "--out-dir", str(out_dir),
                ])

            self.assertEqual(return_code, 0)
            csv_path = out_dir / "main_mesh_confirmation_edges.csv"
            json_path = out_dir / "main_mesh_confirmation_summary.json"
            md_path = out_dir / "main_mesh_confirmation_report.md"
            self.assertTrue(csv_path.exists())
            self.assertTrue(json_path.exists())
            self.assertTrue(md_path.exists())
            payload = json.loads(json_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["schema"], "vela.pn2d_bv_main_mesh_confirmation.v1")
            self.assertEqual(len(payload["anchors"]), 5)
            self.assertEqual(payload["gate"]["status"], "pass")
            self.assertEqual(
                payload["gate"]["next_target"],
                "main_mesh_continuation_branch_recovery",
            )
            report = md_path.read_text(encoding="utf-8")
            self.assertIn("Unique next target", report)
            self.assertIn("Minimum failing test", report)
            self.assertIn("same-area proxy", report)
            self.assertIn("Bidirectional support", report)


if __name__ == "__main__":
    unittest.main()

