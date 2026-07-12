#!/usr/bin/env python3
"""Focused regression tests for compensated-junction SG replay diagnostics."""

from __future__ import annotations

import importlib.util
import json
import math
import tempfile
import unittest
from unittest import mock
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
MODULE_PATH = REPO / "scripts" / "diagnose_pn2d_bv_compensated_source_proxy.py"
SPEC = importlib.util.spec_from_file_location("diagnose_pn2d_bv_compensated_source_proxy_core_test", MODULE_PATH)
assert SPEC is not None
assert SPEC.loader is not None
diagnostic = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(diagnostic)


class CompensatedSgReplayCoreTest(unittest.TestCase):
    def test_production_bernoulli_handles_zero_and_large_arguments(self) -> None:
        self.assertEqual(diagnostic.production_bernoulli(0.0), 1.0)
        self.assertAlmostEqual(diagnostic.production_bernoulli(1.0e-10), 1.0 - 0.5e-10)
        self.assertTrue(math.isfinite(diagnostic.production_bernoulli(1000.0)))
        self.assertGreaterEqual(diagnostic.production_bernoulli(1000.0), 0.0)
        self.assertEqual(diagnostic.production_bernoulli(-1000.0), 1000.0)

    def test_production_sg_replay_covers_normal_flat_clamp_and_cancellation(self) -> None:
        vt = 0.02585
        normal = diagnostic.replay_electron_variable_ni_sg(
            ni0=1.0, ni1=1.0, psi0=0.0, psi1=0.0, phin0=0.0, phin1=vt,
            vt=vt, mobility_m2_V_s=2.0, length_m=vt,
        )
        self.assertAlmostEqual(normal["double_flux_m2_s"], 2.0 * (1.0 - math.exp(-1.0)), places=12)
        self.assertLess(normal["double_highprec_relative_error"], 1.0e-12)
        flat = diagnostic.replay_electron_variable_ni_sg(
            ni0=1.0, ni1=2.0, psi0=-1.0, psi1=1.0, phin0=0.125, phin1=0.125,
            vt=vt, mobility_m2_V_s=0.1, length_m=1.0e-6,
        )
        self.assertEqual(flat["double_flux_m2_s"], 0.0)
        self.assertEqual(flat["decimal_flux_m2_s"], 0.0)
        self.assertTrue(flat["flat_qf_short_circuit"])
        clamped = diagnostic.replay_electron_variable_ni_sg(
            ni0=1.0, ni1=1.0, psi0=600.0 * vt, psi1=0.0, phin0=0.0, phin1=vt,
            vt=vt, mobility_m2_V_s=0.1, length_m=1.0e-6,
        )
        self.assertTrue(clamped["node0_exponent_clamped_high"])
        self.assertTrue(math.isfinite(clamped["double_flux_m2_s"]))
        cancellation = diagnostic.replay_electron_variable_ni_sg(
            ni0=1.0e16, ni1=1.0e16, psi0=10.0, psi1=10.0, phin0=0.0, phin1=1.0e-15,
            vt=vt, mobility_m2_V_s=0.1, length_m=1.0e-6,
        )
        self.assertGreater(cancellation["cancellation_condition"], 1.0e12)
        self.assertGreater(cancellation["double_highprec_relative_error"], 1.0e-6)
        self.assertTrue(math.isfinite(cancellation["decimal_flux_m2_s"]))

    def test_production_sg_replay_rejects_nonfinite_or_nonphysical_inputs(self) -> None:
        base = {
            "ni0": 1.0, "ni1": 1.0, "psi0": 0.0, "psi1": 0.0,
            "phin0": 0.0, "phin1": 0.1, "vt": 0.02585,
            "mobility_m2_V_s": 0.1, "length_m": 1.0e-6,
        }
        for key, value in [("ni0", math.nan), ("vt", 0.0), ("length_m", -1.0)]:
            with self.subTest(key=key):
                args = dict(base)
                args[key] = value
                with self.assertRaises(ValueError):
                    diagnostic.replay_electron_variable_ni_sg(**args)

    def test_manifest_requires_two_component_e_current_density(self) -> None:
        manifest = {"fields": [
            {"name": "ElectrostaticPotential", "components": 1},
            {"name": "eCurrentDensity", "components": 1},
            {"name": "eCurrentDensity", "components": 2},
        ]}
        field = diagnostic.validate_manifest_vector_components(
            manifest, "eCurrentDensity", expected_components=2
        )
        self.assertEqual(field["components"], 2)
        with self.assertRaisesRegex(ValueError, "components=2"):
            diagnostic.validate_manifest_vector_components(
                {"fields": [{"name": "eCurrentDensity", "components": 1}]},
                "eCurrentDensity", expected_components=2,
            )
        with self.assertRaisesRegex(ValueError, "missing"):
            diagnostic.validate_manifest_vector_components(
                {"fields": []}, "eCurrentDensity", expected_components=2
            )

    def test_vector_projection_uses_canonical_orientation_and_electron_sign(self) -> None:
        charge = diagnostic.ELEMENTARY_CHARGE_C
        cases = [
            ((0.0, 0.0), (2.0, 0.0), (2.0, 4.0), (6.0, 8.0), 4.0, False),
            ((2.0, 0.0), (0.0, 0.0), (6.0, 8.0), (2.0, 4.0), 4.0, True),
            ((0.0, 2.0), (0.0, 0.0), (8.0, 6.0), (4.0, 2.0), 4.0, True),
        ]
        for point0, point1, current0, current1, expected, reversed_input in cases:
            with self.subTest(point0=point0, point1=point1):
                result = diagnostic.project_endpoint_current_to_canonical_edge(
                    point0=point0, point1=point1,
                    current0_A_cm2=current0, current1_A_cm2=current1,
                )
                self.assertAlmostEqual(result["conventional_current_A_cm2"], expected)
                self.assertEqual(result["input_orientation_reversed"], reversed_input)
                self.assertAlmostEqual(
                    result["electron_continuity_flux_m2_s"], -expected * 1.0e4 / charge
                )

    def test_108_row_key_validator_requires_exact_unique_matrix(self) -> None:
        rows = [
            {"variant": variant, "bias_V": bias, "y_um": y_um, "side": side}
            for variant in diagnostic.REPLAY_VARIANTS
            for bias in diagnostic.BIASES
            for y_um in diagnostic.Y_CUTS
            for side in ("left", "right")
        ]
        self.assertEqual(len(diagnostic.validate_108_row_keys(rows)), 108)
        duplicate = [dict(row) for row in rows]
        duplicate[-1] = dict(duplicate[0])
        with self.assertRaisesRegex(ValueError, "duplicate"):
            diagnostic.validate_108_row_keys(duplicate)
        with self.assertRaisesRegex(ValueError, "expected 108"):
            diagnostic.validate_108_row_keys(rows[:-1])

    def test_required_field_validator_rejects_missing_and_nonfinite_values(self) -> None:
        row = {"a": 1.0, "b": -2.0}
        diagnostic.require_finite_fields(row, ("a", "b"), context="fixture")
        with self.assertRaisesRegex(ValueError, "missing"):
            diagnostic.require_finite_fields(row, ("a", "missing"), context="fixture")
        with self.assertRaisesRegex(ValueError, "non-finite"):
            diagnostic.require_finite_fields({"a": math.inf}, ("a",), context="fixture")

    def test_root_cause_classifier_uses_strict_ordered_thresholds(self) -> None:
        cases = [
            ({"double_highprec_relative_error": 1.000001e-6, "cancellation_condition": 1.0e12},
             "variable_ni_sg_numerical_stability"),
            ({"sent_state_gap_recovery": 0.8, "sent_state_replay_residual_dex": 0.1},
             "vela_internal_state_branch"),
            ({"sent_state_vector_residual_dex": 0.200001},
             "sg_discretization_ni_or_current_semantics"),
            ({"raw_edge_residual_dex": 0.1, "alpha_residual_dex": 0.100001,
              "source_residual_dex": 0.3},
             "impact_coefficient_or_source_semantics"),
            ({"raw_edge_residual_dex": 0.1, "alpha_residual_dex": 0.1,
              "source_residual_dex": 0.200001},
             "ownership_support_mapping"),
            ({"coarse_only": True}, "inconclusive"),
            ({"coarse_only": True, "main_comparison_supports_same_failure": False},
             "coarse_artifact"),
            ({}, "inconclusive"),
        ]
        for evidence, expected in cases:
            with self.subTest(expected=expected):
                result = diagnostic.classify_root_cause(evidence)
                self.assertEqual(result["classification"], expected)
                self.assertIn("rule", result)
                self.assertEqual(result["evidence"], evidence)

    def test_standard_variant_root_resolves_explicit_two_by_three_inputs(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_compensated_variant_root_") as td:
            root = Path(td)
            specs = diagnostic.standard_variant_inputs(root)
            self.assertEqual(list(specs), [
                "legacy_density_gradient",
                "legacy_gss_midpoint",
                "legacy_triangle_gss_gradqf",
                "reported_density_gradient",
                "reported_gss_midpoint",
                "reported_triangle_gss_gradqf",
            ])
            self.assertEqual(
                specs["legacy_density_gradient"]["run_root"],
                root / "legacy_density_gradient" / "run",
            )
            self.assertEqual(
                specs["reported_gss_midpoint"]["doping_csv"],
                root / "reported_gss_midpoint" / "imported" / "vela" / "doping.csv",
            )
            self.assertEqual(
                specs["legacy_gss_midpoint"]["compensated_doping_policy"],
                "dominant_signed_region",
            )
            self.assertEqual(
                specs["reported_density_gradient"]["compensated_doping_policy"], "reported"
            )
            self.assertEqual(
                [spec["current_variant"] for spec in specs.values()],
                [
                    "density_gradient", "gss_midpoint", "triangle_gss_gradqf",
                    "density_gradient", "gss_midpoint", "triangle_gss_gradqf",
                ],
            )
            self.assertEqual(
                [spec["current_approximation"] for spec in specs.values()],
                [
                    "density_gradient", "cell_reconstructed", "cell_reconstructed",
                    "density_gradient", "cell_reconstructed", "cell_reconstructed",
                ],
            )

    def test_gss_midpoint_ratios_pair_with_density_gradient_within_doping(self) -> None:
        rows = []
        for variant, current_variant, value in (
            ("legacy_density_gradient", "density_gradient", 2.0),
            ("legacy_gss_midpoint", "gss_midpoint", 6.0),
        ):
            for side in ("left", "right"):
                rows.append({
                    "variant": variant,
                    "doping_strategy": "legacy",
                    "current_variant": current_variant,
                    "bias_V": -12.0,
                    "y_um": 0.0,
                    "side": side,
                    **{
                        field: value
                        for field in diagnostic.STANDARD_REPLAY_RATIO_FIELDS
                    },
                })
        diagnostic._append_standard_ratios(rows)
        pairs = diagnostic.current_discretization_pair_rows(rows)

        self.assertEqual(len(pairs), 2)
        self.assertTrue(all(
            item["gss_midpoint_over_density_gradient_edge_source_integral"] == 3.0
            for item in pairs
        ))

    def test_variant_run_status_records_last_converged_bias_and_handoff(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_compensated_status_") as td:
            root = Path(td)
            spec = diagnostic.standard_variant_inputs(root)["legacy_density_gradient"]
            spec["run_root"].mkdir(parents=True)
            spec["deck_path"].write_text(json.dumps({
                "output_csv": "pn2d_bv_density_gradient.csv",
            }), encoding="utf-8")
            (spec["run_root"] / "pn2d_bv_density_gradient.csv").write_text(
                "bias_V,converged,handoff_stage\n"
                "0,1,newton\n"
                "-0.05,1,newton\n"
                "-0.1,0,newton_failed\n",
                encoding="utf-8",
            )

            status = diagnostic.variant_run_status(spec)

            self.assertEqual(status["run_status"], "partial")
            self.assertEqual(status["last_converged_bias_V"], -0.05)
            self.assertEqual(status["handoff_stage"], "newton_failed")

    def test_cli_accepts_standardized_variant_root_and_legacy_paths(self) -> None:
        standardized = diagnostic.parse_args([
            "--variants-root", "variants",
            "--sentaurus-root", "sentaurus",
            "--out-dir", "report",
        ])
        self.assertEqual(standardized.variants_root, Path("variants"))
        self.assertIsNone(standardized.baseline_report_root)
        legacy = diagnostic.parse_args([
            "--baseline-report-root", "baseline",
            "--probe-root", "probe",
            "--sentaurus-root", "sentaurus",
            "--out-dir", "report",
        ])
        self.assertEqual(legacy.baseline_report_root, Path("baseline"))
        self.assertEqual(legacy.probe_root, Path("probe"))

    def test_sentaurus_field_loader_reads_endpoint_scalars_and_vector(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_compensated_sentaurus_fields_") as td:
            export_dir = Path(td)
            fields = export_dir / "fields"
            fields.mkdir()
            manifest = {
                "fields": [
                    {"name": "ElectrostaticPotential", "components": 1, "unit": "V", "mapping_status": "complete", "region": 0},
                    {"name": "eQuasiFermiPotential", "components": 1, "unit": "V", "mapping_status": "complete", "region": 0},
                    {"name": "eDensity", "components": 1, "unit": "cm^-3", "mapping_status": "complete", "region": 0},
                    {"name": "eMobility", "components": 1, "unit": "cm^2*V^-1*s^-1", "mapping_status": "complete", "region": 0},
                    {"name": "eAlphaAvalanche", "components": 1, "unit": "cm^-1", "mapping_status": "complete", "region": 0},
                    {"name": "eCurrentDensity", "components": 1, "unit": "A*cm^-2", "mapping_status": "complete", "region": 0},
                    {"name": "eCurrentDensity", "components": 2, "unit": "A*cm^-2", "mapping_status": "complete", "region": 0},
                ],
            }
            for field in manifest["fields"]:
                field["global_node_mapping"] = "global_vertex_order"
            manifest_path = export_dir / "field_manifest.json"
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
            scalar_values = {
                "ElectrostaticPotential": (0.1, 0.2),
                "eQuasiFermiPotential": (-0.1, -0.2),
                "eDensity": (2.0e10, 3.0e10),
                "eMobility": (1000.0, 1200.0),
                "eAlphaAvalanche": (3.0e4, 4.0e4),
            }
            for name, values in scalar_values.items():
                (fields / f"{name}_region0.csv").write_text(
                    f"node_id,component0\n0,{values[0]}\n1,{values[1]}\n",
                    encoding="utf-8",
                )
            (fields / "eCurrentDensity_region0.csv").write_text(
                "node_id,component0,component1\n0,2,4\n1,6,8\n",
                encoding="utf-8",
            )
            state = diagnostic.load_sentaurus_electron_state(export_dir)
            self.assertEqual(state["psi_V"], {0: 0.1, 1: 0.2})
            self.assertEqual(state["density_m3"], {0: 2.0e16, 1: 3.0e16})
            self.assertAlmostEqual(state["mobility_m2_V_s"][0], 0.1)
            self.assertAlmostEqual(state["mobility_m2_V_s"][1], 0.12)
            self.assertEqual(state["alpha_m_inv"], {0: 3.0e6, 1: 4.0e6})
            self.assertEqual(state["current_A_cm2"], {0: (2.0, 4.0), 1: (6.0, 8.0)})

            mutations = [
                ("missing scalar", lambda data: data["fields"].pop(2)),
                ("wrong unit", lambda data: data["fields"][3].update(unit="m2/V/s")),
                ("wrong components", lambda data: data["fields"][0].update(components=2)),
                ("incomplete mapping", lambda data: data["fields"][-1].update(mapping_status="partial")),
                ("wrong region", lambda data: data["fields"][1].update(region=1)),
            ]
            for label, mutate in mutations:
                with self.subTest(label=label):
                    broken = json.loads(json.dumps(manifest))
                    mutate(broken)
                    manifest_path.write_text(json.dumps(broken), encoding="utf-8")
                    with self.assertRaises(ValueError):
                        diagnostic.load_sentaurus_electron_state(export_dir)
            manifest_path.write_text(json.dumps(manifest), encoding="utf-8")

    def test_enriched_edge_fields_append_all_sg_sentaurus_and_replay_columns(self) -> None:
        edge_row = {field: "1.0" for field in diagnostic.ELECTRON_SG_FIELDS}
        edge_row.update({
            "x0_um": "0.75",
            "y0_um": "0.0",
            "x1_um": "1.0",
            "y1_um": "0.0",
            "edge_length_m": "1e-6",
            "edge_area_proxy_m2": "2e-12",
            "electron_mobility_m2_V_s": "0.11",
            "electron_alpha_m_inv": "3e6",
            "electron_source_integral": "6.0",
        })
        sentaurus_state = {
            "psi_V": {0: 0.1, 1: 0.2},
            "phin_V": {0: -0.1, 1: -0.2},
            "density_m3": {0: 2.0e16, 1: 3.0e16},
            "mobility_m2_V_s": {0: 0.1, 1: 0.12},
            "alpha_m_inv": {0: 3.0e6, 1: 4.0e6},
            "current_A_cm2": {0: (2.0, 4.0), 1: (6.0, 8.0)},
        }
        enriched = diagnostic.enrich_edge_with_sentaurus_replay(
            edge_row=edge_row,
            sentaurus_state=sentaurus_state,
            sentaurus_node0={"id": 0, "x_um": 0.75, "y_um": 0.0},
            sentaurus_node1={"id": 1, "x_um": 1.0, "y_um": 0.0},
            temperature_K=300.0,
            unit_system="tcad_internal",
        )
        self.assertTrue(set(diagnostic.ELECTRON_SG_FIELDS).issubset(enriched))
        self.assertEqual(enriched["sentaurus_e_current_edge_signed_A_cm2"], 4.0)
        self.assertLess(enriched["sentaurus_e_continuity_edge_signed_flux_m2_s"], 0.0)
        self.assertIn("sentaurus_e_ni_inferred0_m3", enriched)
        self.assertIn("sentaurus_e_sg_vela_mobility_signed_flux_m2_s", enriched)
        self.assertIn("sentaurus_e_sg_sentaurus_mobility_signed_flux_m2_s", enriched)
        self.assertAlmostEqual(enriched["vela_e_source_integral_physical_m_inv_s"], 6.0e-6)
        self.assertAlmostEqual(enriched["vela_e_source_closure_ratio"], 1.0)
        self.assertAlmostEqual(enriched["sentaurus_edge_length_m"], 0.25e-6)
        self.assertGreater(
            enriched["sentaurus_e_source_on_vela_area_physical_m_inv_s"], 0.0
        )
        self.assertGreater(enriched["vela_e_over_sentaurus_source_abs_ratio"], 0.0)
        self.assertTrue(math.isfinite(enriched["vela_e_over_sentaurus_source_abs_ratio"]))
        diagnostic.require_finite_fields(
            enriched,
            diagnostic.REQUIRED_ENRICHED_FIELDS,
            context="enriched fixture",
        )


    def test_zero_ratios_remain_finite_for_flat_flux_cases(self) -> None:
        cases = ((0.0, 1.0), (1.0, 0.0), (0.0, 0.0))
        for numerator, denominator in cases:
            with self.subTest(numerator=numerator, denominator=denominator):
                ratio = diagnostic._finite_abs_ratio(numerator, denominator)
                self.assertGreater(ratio, 0.0)
                self.assertTrue(math.isfinite(ratio))
                self.assertTrue(math.isfinite(diagnostic._log_gap_from_abs_ratio(ratio)))

    def test_load_sg_edges_rejects_duplicate_bias_edge_key(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_duplicate_sg_edge_") as td:
            path = Path(td) / "sg.csv"
            path.write_text(
                "bias_V,edge_id,value\n-19,13,1\n-19,13,2\n",
                encoding="utf-8",
            )
            with self.assertRaisesRegex(
                ValueError, "duplicate SG edge row for bias=-19, edge_id=13"
            ):
                diagnostic.load_sg_edges(path)

    def test_triangle_source_loader_aggregates_adjacent_cell_records_by_edge(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_triangle_source_edge_") as td:
            path = Path(td) / "triangle.csv"
            path.write_text(
                "bias_V,edge_id,node0,node1,x0_um,y0_um,x1_um,y1_um,"
                "truncated_partial_volume_m2,electron_flux_proxy,hole_flux_proxy,"
                "electron_alpha_m_inv,hole_alpha_m_inv,electron_mobility_m2_V_s,"
                "hole_mobility_m2_V_s,electron_midpoint_density_m3,"
                "hole_midpoint_density_m3,electron_source_integral,"
                "hole_source_integral,edge_source_integral,node0_source_integral,"
                "node1_source_integral\n"
                "-19,13,7,10,0,0,1,0,2,4,3,10,5,100,50,1000,500,80,30,110,55,55\n"
                "-19,13,10,7,1,0,0,0,3,6,2,20,8,200,80,2000,800,360,48,408,204,204\n",
                encoding="utf-8",
            )
            rows = diagnostic.load_triangle_gss_source_edges(path)
            self.assertEqual(list(rows), [(-19.0, 13)])
            row = rows[(-19.0, 13)]
            self.assertAlmostEqual(float(row["edge_area_proxy_m2"]), 5.0)
            self.assertAlmostEqual(float(row["electron_source_integral"]), 440.0)
            self.assertAlmostEqual(float(row["hole_source_integral"]), 78.0)
            self.assertAlmostEqual(float(row["edge_source_integral"]), 518.0)
            self.assertAlmostEqual(float(row["node0_source_integral"]), 259.0)
            self.assertAlmostEqual(float(row["node1_source_integral"]), 259.0)
            self.assertAlmostEqual(float(row["electron_flux_proxy"]), 5.2)
            self.assertAlmostEqual(float(row["electron_alpha_m_inv"]), 440.0 / 26.0)
            self.assertAlmostEqual(float(row["electron_mobility_m2_V_s"]), 160.0)
            self.assertAlmostEqual(float(row["electron_density_mid_m3"]), 1600.0)
            self.assertEqual(row["electron_raw_flux_proxy"], row["electron_flux_proxy"])

    def test_sg_edge_mapping_uses_unique_endpoint_nodes_not_historical_id(self) -> None:
        edges = {
            (-12.0, 22): {"node0": "7", "node1": "10"},
            (-12.0, 29): {"node0": "10", "node1": "13"},
        }
        edge_id, row = diagnostic.unique_sg_edge_for_nodes(edges, -12.0, 10, 7)
        self.assertEqual(edge_id, 22)
        self.assertIs(row, edges[(-12.0, 22)])
        with self.assertRaisesRegex(ValueError, "found 0"):
            diagnostic.unique_sg_edge_for_nodes(edges, -12.0, 7, 13)
        ambiguous = dict(edges)
        ambiguous[(-12.0, 99)] = {"node0": "10", "node1": "7"}
        with self.assertRaisesRegex(ValueError, "found 2"):
            diagnostic.unique_sg_edge_for_nodes(ambiguous, -12.0, 7, 10)



    def test_enrichment_requires_finite_vela_coordinates_and_reverses_sign(self) -> None:
        edge_row = {field: "1.0" for field in diagnostic.ELECTRON_SG_FIELDS}
        edge_row.update({
            "x0_um": "1.0",
            "y0_um": "0.0",
            "x1_um": "0.75",
            "y1_um": "0.0",
            "edge_length_m": "1e-6",
            "edge_area_proxy_m2": "2e-12",
            "electron_mobility_m2_V_s": "0.11",
            "electron_alpha_m_inv": "3e6",
            "electron_source_integral": "6.0",
            "electron_sg_production_signed_continuity_particle_flux_m2_s": "1.0",
        })
        state = {
            "psi_V": {0: 0.1, 1: 0.2},
            "phin_V": {0: -0.1, 1: -0.2},
            "density_m3": {0: 2.0e16, 1: 3.0e16},
            "mobility_m2_V_s": {0: 0.1, 1: 0.12},
            "alpha_m_inv": {0: 3.0e6, 1: 4.0e6},
            "current_A_cm2": {0: (2.0, 4.0), 1: (6.0, 8.0)},
        }
        result = diagnostic.enrich_edge_with_sentaurus_replay(
            edge_row=edge_row,
            sentaurus_state=state,
            sentaurus_node0={"id": 0, "x_um": 0.75, "y_um": 0.0},
            sentaurus_node1={"id": 1, "x_um": 1.0, "y_um": 0.0},
            temperature_K=300.0,
            unit_system="tcad_internal",
        )
        self.assertEqual(result["vela_edge_orientation_sign_to_canonical"], -1.0)
        self.assertEqual(result["vela_e_sg_production_canonical_signed_flux_m2_s"], -1.0)
        missing = dict(edge_row)
        del missing["x0_um"]
        with self.assertRaisesRegex(ValueError, "x0_um"):
            diagnostic.enrich_edge_with_sentaurus_replay(
                edge_row=missing,
                sentaurus_state=state,
                sentaurus_node0={"id": 0, "x_um": 0.75, "y_um": 0.0},
                sentaurus_node1={"id": 1, "x_um": 1.0, "y_um": 0.0},
                temperature_K=300.0,
                unit_system="tcad_internal",
            )
        with self.assertRaisesRegex(ValueError, "unit_system"):
            diagnostic.enrich_edge_with_sentaurus_replay(
                edge_row=edge_row,
                sentaurus_state=state,
                sentaurus_node0={"id": 0, "x_um": 0.75, "y_um": 0.0},
                sentaurus_node1={"id": 1, "x_um": 1.0, "y_um": 0.0},
                temperature_K=300.0,
                unit_system="si",
            )

    def test_standard_build_detail_rows_produces_enriched_classified_matrix(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_compensated_standard_build_") as td:
            root = Path(td)
            variants_root = root / "variants"
            sentaurus_root = root / "sentaurus"
            args = diagnostic.parse_args([
                "--variants-root", str(variants_root),
                "--sentaurus-root", str(sentaurus_root),
                "--out-dir", str(root / "report"),
            ])
            points = [
                (x, y, 0.0)
                for y in diagnostic.Y_CUTS
                for x in diagnostic.VELA_X_COLUMNS
            ]
            count = len(points)
            vela_state = {
                "points": points,
                "scalars": {
                    "Potential": [0.1 + index * 1.0e-3 for index in range(count)],
                    "ElectronQuasiFermi": [-0.1 - index * 1.0e-3 for index in range(count)],
                    "HoleQuasiFermi": [0.2 + index * 1.0e-3 for index in range(count)],
                    "Electrons": [2.0e16 for _ in range(count)],
                    "Holes": [3.0e15 for _ in range(count)],
                },
            }
            doping = {
                index: {
                    "type": "compensated" if index % 3 == 1 else ("p" if index % 3 == 0 else "n"),
                    "net_doping_cm3": 0.0,
                }
                for index in range(count)
            }
            sentaurus_nodes = [
                {"id": index, "x_um": x, "y_um": y}
                for index, (x, y, _z) in enumerate([
                    (x, y, 0.0)
                    for y in diagnostic.Y_CUTS
                    for x in diagnostic.SENTAURUS_X_COLUMNS
                ])
            ]
            sentaurus_state = {
                "psi_V": {node["id"]: 0.1 + node["id"] * 1.0e-3 for node in sentaurus_nodes},
                "phin_V": {node["id"]: -0.1 - node["id"] * 1.0e-3 for node in sentaurus_nodes},
                "density_m3": {node["id"]: 2.0e16 for node in sentaurus_nodes},
                "mobility_m2_V_s": {node["id"]: 0.11 for node in sentaurus_nodes},
                "alpha_m_inv": {node["id"]: 3.0e6 for node in sentaurus_nodes},
                "current_A_cm2": {node["id"]: (-0.2, 0.0) for node in sentaurus_nodes},
            }
            sg_rows = {}
            for bias in diagnostic.BIASES:
                for y_um in diagnostic.Y_CUTS:
                    for side in ("left", "right"):
                        edge_id = diagnostic.EDGE_BY_SIDE[y_um][side]
                        x0, x1 = (
                            (diagnostic.VELA_X_COLUMNS[0], diagnostic.VELA_X_COLUMNS[1])
                            if side == "left"
                            else (diagnostic.VELA_X_COLUMNS[1], diagnostic.VELA_X_COLUMNS[2])
                        )
                        node0 = diagnostic.nearest_node(points, x0, y_um)
                        node1 = diagnostic.nearest_node(points, x1, y_um)
                        edge = {field: "1.0" for field in diagnostic.ELECTRON_SG_FIELDS}
                        edge.update({
                            "bias_V": str(bias),
                            "edge_id": str(edge_id),
                            "node0": str(node0),
                            "node1": str(node1),
                            "x0_um": str(x0),
                            "y0_um": str(y_um),
                            "x1_um": str(x1),
                            "y1_um": str(y_um),
                            "edge_length_m": str((x1 - x0) * 1.0e-6),
                            "edge_couple_m": "1e-6",
                            "edge_area_proxy_m2": "2e-12",
                            "electric_field_V_per_m": "1e6",
                            "electron_impact_field_V_per_m": "1e6",
                            "hole_impact_field_V_per_m": "1e6",
                            "electron_alpha_m_inv": "3e6",
                            "hole_alpha_m_inv": "2e6",
                            "electron_mobility_m2_V_s": "0.11",
                            "hole_mobility_m2_V_s": "0.05",
                            "electron_flux_proxy": "1.0",
                            "hole_flux_proxy": "1.0",
                            "electron_raw_flux_proxy": "1.0",
                            "hole_raw_flux_proxy": "1.0",
                            "electron_reconstructed_flux_proxy": "1.0",
                            "hole_reconstructed_flux_proxy": "1.0",
                            "electron_final_over_raw_flux_proxy": "1.0",
                            "hole_final_over_raw_flux_proxy": "1.0",
                            "electron_source_integral": "6.0",
                            "hole_source_integral": "1.0",
                            "edge_source_integral": "7.0",
                            "electron_sg_production_signed_continuity_particle_flux_m2_s": "1.0",
                            "electron_sg_production_vs_high_precision_reference_relative_error": "1e-8",
                            "electron_sg_cancellation_condition": "10.0",
                            "electron_sg_node0_exponent_clamped_low": "0",
                            "electron_sg_node0_exponent_clamped_high": "0",
                            "electron_sg_node1_exponent_clamped_low": "0",
                            "electron_sg_node1_exponent_clamped_high": "0",
                        })
                        sg_rows[(round(bias, 10), edge_id)] = edge

            with (
                mock.patch.object(diagnostic, "load_sg_edges", return_value=sg_rows) as load_sg,
                mock.patch.object(
                    diagnostic, "load_triangle_gss_source_edges",
                    return_value=sg_rows,
                ) as load_triangle,
                mock.patch.object(diagnostic, "load_doping", return_value=doping),
                mock.patch.object(
                    diagnostic,
                    "load_sentaurus_nodes",
                    return_value=(sentaurus_nodes, doping),
                ),
                mock.patch.object(
                    diagnostic,
                    "load_sentaurus_electron_state",
                    return_value=sentaurus_state,
                ),
                mock.patch.object(diagnostic, "vtk_for_bias", return_value=Path("state.vtk")),
                mock.patch.object(diagnostic, "parse_vtk", return_value=vela_state),
            ):
                rows = diagnostic.build_standard_detail_rows(args)

            self.assertEqual(len(rows), 108)
            self.assertEqual({row["variant"] for row in rows}, set(diagnostic.REPLAY_VARIANTS))
            self.assertTrue(all("root_cause_classification" in row for row in rows))
            self.assertTrue(all("classifier_gap_recovery" in row for row in rows))
            self.assertEqual(
                [call.args[0] for call in load_sg.call_args_list],
                [
                    diagnostic.standard_variant_inputs(variants_root)[variant]["sg_csv"]
                    for variant in diagnostic.REPLAY_VARIANTS
                ],
            )
            self.assertEqual(
                [call.args[0] for call in load_triangle.call_args_list],
                [
                    diagnostic.standard_variant_inputs(variants_root)[variant]["source_csv"]
                    for variant in diagnostic.REPLAY_VARIANTS
                    if "triangle_gss_gradqf" in variant
                ],
            )

    def test_standard_main_writes_generic_structured_outputs(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_compensated_standard_main_") as td:
            root = Path(td)
            rows = [
                {
                    "variant": variant,
                    "bias_V": bias,
                    "y_um": y_um,
                    "side": side,
                    "root_cause_classification": "inconclusive",
                    "root_cause_rule": "fixture",
                    **{field: 1.0 for field in (
                        diagnostic.REQUIRED_ENRICHED_FIELDS
                        + diagnostic.CURRENT_DISCRETIZATION_RATIO_FIELDS
                    )},
                }
                for variant in diagnostic.REPLAY_VARIANTS
                for bias in diagnostic.BIASES
                for y_um in diagnostic.Y_CUTS
                for side in ("left", "right")
            ]
            out_dir = root / "report"
            with mock.patch.object(
                diagnostic,
                "build_standard_detail_rows",
                return_value=rows,
            ):
                diagnostic.main([
                    "--variants-root", str(root / "variants"),
                    "--sentaurus-root", str(root / "sentaurus"),
                    "--out-dir", str(out_dir),
                ])
            csv_path = out_dir / "compensated_sg_replay.csv"
            json_path = out_dir / "compensated_sg_replay.json"
            report_path = out_dir / "compensated_sg_replay_report.md"
            self.assertTrue(csv_path.exists())
            self.assertTrue(json_path.exists())
            self.assertTrue(report_path.exists())
            payload = json.loads(json_path.read_text(encoding="utf-8"))
            self.assertEqual(payload["schema"], "vela.pn2d_bv_compensated_sg_replay.v3")
            self.assertEqual(
                payload["summary"]["schema"],
                "vela.pn2d_bv_compensated_sg_replay.summary.v3",
            )
            self.assertEqual(payload["row_count"], 108)
            self.assertEqual(len(payload["classifications"]), 18)
            self.assertEqual(len(payload["current_discretization_pairs"]), 72)
            self.assertEqual(len(payload["run_statuses"]), 6)
            self.assertEqual(
                {item["pair_family"] for item in payload["current_discretization_pairs"]},
                {"gss_over_density", "triangle_over_gss"},
            )
            self.assertEqual(
                {
                    (item["numerator_current_variant"], item["denominator_current_variant"])
                    for item in payload["current_discretization_pairs"]
                },
                {
                    ("gss_midpoint", "density_gradient"),
                    ("triangle_gss_gradqf", "gss_midpoint"),
                },
            )
            self.assertTrue(all("evidence" in item for item in payload["classifications"]))
            self.assertIn("Structured Root-Cause Classifications", report_path.read_text(encoding="utf-8"))

    def test_enriched_row_matrix_requires_all_required_columns(self) -> None:
        rows = [
            {
                "variant": variant,
                "bias_V": bias,
                "y_um": y_um,
                "side": side,
                **{field: 1.0 for field in (
                    diagnostic.REQUIRED_ENRICHED_FIELDS
                    + diagnostic.CURRENT_DISCRETIZATION_RATIO_FIELDS
                )},
            }
            for variant in diagnostic.REPLAY_VARIANTS
            for bias in diagnostic.BIASES
            for y_um in diagnostic.Y_CUTS
            for side in ("left", "right")
        ]
        self.assertEqual(len(diagnostic.validate_enriched_rows(rows)), 108)
        del rows[0][diagnostic.REQUIRED_ENRICHED_FIELDS[0]]
        with self.assertRaisesRegex(ValueError, "missing"):
            diagnostic.validate_enriched_rows(rows)


if __name__ == "__main__":
    unittest.main()

