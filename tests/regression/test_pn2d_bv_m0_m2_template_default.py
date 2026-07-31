#!/usr/bin/env python3
"""Regression tests for prospective M0/M2 template-default acceptance."""

from __future__ import annotations

import csv
import json
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
if str(REPO / "scripts") not in sys.path:
    sys.path.insert(0, str(REPO / "scripts"))

from generate_pn2d_config import render_named_template  # noqa: E402
from review_pn2d_bv_m0_m2_template_default import (  # noqa: E402
    BRANCHES,
    closure_evidence,
    evaluate_level,
    parity_evidence,
    process_probe_evidence,
    sha256,
)


CONTRACT_PATH = (
    REPO
    / "docs"
    / "validation"
    / "contracts"
    / "pn2d_bv_m0_m2_template_default_acceptance_v1.json"
)
CONTRACT_V2_PATH = CONTRACT_PATH.with_name(
    "pn2d_bv_m0_m2_template_default_acceptance_v2.json"
)


def write_json(path: Path, payload) -> None:
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


class Pn2dBvM0M2TemplateDefaultTest(unittest.TestCase):
    def setUp(self) -> None:
        self.contract = json.loads(CONTRACT_PATH.read_text(encoding="utf-8"))
        self.biases = self.contract["bv_domain"]["exact_biases_V"]

    def test_thresholds_are_copied_from_the_frozen_dual_domain_contract(self) -> None:
        prior = json.loads(
            (
                CONTRACT_PATH.parent / "pn2d_bv_dual_domain_acceptance_v1.json"
            ).read_text(encoding="utf-8")
        )
        self.assertEqual(
            self.contract["bv_domain"]["thresholds"],
            prior["domains"]["bv_model_consistency"]["thresholds"],
        )
        self.assertFalse(
            self.contract["scope"]["global_cpp_default_change_authorized"]
        )
        self.assertTrue(
            self.contract["scope"]["cross_mesh_convergence_is_observation_only"]
        )

    def test_parity_evidence_uses_only_the_predeclared_bv_domain(self) -> None:
        parity = self._parity()
        parity["curve_rows"].append(
            {
                "bias_V": -3.0,
                "vela_on_A_per_um": 1.0e9,
                "sentaurus_on_A_per_um": 1.0e-30,
                "vela_gain": 1.0e9,
                "sentaurus_gain": 1.0e-30,
            }
        )
        result = parity_evidence(
            parity,
            tuple(self.biases),
            self.contract["bv_domain"]["thresholds"],
        )
        self.assertTrue(result["passed"])

    def test_v2_accepts_predeclared_shared_absent_slope_crossing(self) -> None:
        contract = json.loads(CONTRACT_V2_PATH.read_text(encoding="utf-8"))
        parity = self._parity()
        parity["knee_estimators"]["vela"]["V_slope"] = None
        parity["knee_estimators"]["sentaurus"]["V_slope"] = None
        result = parity_evidence(
            parity,
            tuple(self.biases),
            contract["bv_domain"]["thresholds"],
            contract["bv_domain"]["V_slope_policy"],
        )
        self.assertTrue(result["gates"]["V_slope_abs_error_V"])
        self.assertEqual(
            result["V_slope_outcome"],
            "shared_no_slope_crossing_in_frozen_window",
        )

    def test_v2_rejects_one_sided_absent_slope_crossing(self) -> None:
        contract = json.loads(CONTRACT_V2_PATH.read_text(encoding="utf-8"))
        parity = self._parity()
        parity["knee_estimators"]["vela"]["V_slope"] = None
        result = parity_evidence(
            parity,
            tuple(self.biases),
            contract["bv_domain"]["thresholds"],
            contract["bv_domain"]["V_slope_policy"],
        )
        self.assertFalse(result["gates"]["V_slope_abs_error_V"])
        self.assertEqual(
            result["V_slope_outcome"],
            "one_sided_no_slope_crossing_in_frozen_window",
        )

    def test_closure_requires_machine_readable_global_columns(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_m0_m2_closure_") as td:
            root = Path(td)
            iv = root / "iv.csv"
            probe = root / "process_probe.csv"
            self._write_iv(iv, closure=True)
            probe.write_text("bias_V,source_integral\n-20,1\n", encoding="utf-8")
            accepted = closure_evidence(
                iv,
                probe,
                tuple(self.biases),
                self.contract["closure"],
            )
            self.assertTrue(accepted["passed"])

            self._write_iv(iv, closure=False)
            rejected = closure_evidence(
                iv,
                probe,
                tuple(self.biases),
                self.contract["closure"],
            )
            self.assertFalse(rejected["passed"])
            self.assertFalse(rejected["gates"]["columns_present"])

    def test_process_probe_requires_columns_and_every_contract_bias(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_m0_m2_probe_") as td:
            path = Path(td) / "process_probe.csv"
            self._write_probe(path, self.biases)
            accepted = process_probe_evidence(path, tuple(self.biases))
            self.assertTrue(accepted["passed"])
            self._write_probe(path, self.biases[:-1])
            rejected = process_probe_evidence(path, tuple(self.biases))
            self.assertFalse(rejected["gates"]["exact_bias_coverage"])

    def test_complete_level_binds_config_closure_iv_and_state_hashes(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_m0_m2_level_") as td:
            root = Path(td)
            paths = self._level_fixture(root)
            result = evaluate_level("M2", self.contract, **paths)
            self.assertEqual(result["status"], "passed")
            self.assertTrue(result["gates"]["duplicate_determinism"])
            self.assertEqual(
                result["artifact_bindings"]["physical_inputs"]["mesh_file"][
                    "sha256"
                ],
                sha256(root / "mesh.json"),
            )

            state_b = json.loads(paths["state_b_path"].read_text(encoding="utf-8"))
            state_b["branch_records"][0]["bias_records"][0]["snapshot_tdr"][
                "sha256"
            ] = "different"
            write_json(paths["state_b_path"], state_b)
            rejected = evaluate_level("M2", self.contract, **paths)
            self.assertEqual(rejected["status"], "failed")
            self.assertFalse(rejected["gates"]["duplicate_determinism"])

    def test_self_consistent_but_wrong_curve_artifact_is_rejected(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_m0_m2_binding_") as td:
            root = Path(td)
            paths = self._level_fixture(root)
            parity = json.loads(paths["parity_path"].read_text(encoding="utf-8"))
            wrong = root / "wrong.csv"
            self._write_iv(wrong, closure=True)
            parity["inputs"]["vela_on"] = {
                "path": str(wrong),
                "sha256": sha256(wrong),
            }
            write_json(paths["parity_path"], parity)
            rejected = evaluate_level("M2", self.contract, **paths)
            self.assertFalse(rejected["binding_checks"]["curve_input_paths_valid"])
            self.assertEqual(rejected["status"], "failed")

    def test_v2_rejects_explicit_profile_override_as_default_evidence(self) -> None:
        contract = json.loads(CONTRACT_V2_PATH.read_text(encoding="utf-8"))
        with tempfile.TemporaryDirectory(prefix="vela_m0_m2_default_path_") as td:
            paths = self._level_fixture(Path(td))
            result = evaluate_level("M0", contract, **paths)
            self.assertFalse(
                result["configuration"]["gates"][
                    "default_render_has_no_profile_override"
                ]
            )
            self.assertEqual(result["status"], "failed")

    def _parity(self):
        return {
            "curve_rows": [
                {
                    "bias_V": bias,
                    "vela_on_A_per_um": 1.0e-10 * (1.0 + abs(bias)),
                    "sentaurus_on_A_per_um": 1.0e-10 * (1.0 + abs(bias)),
                    "vela_gain": 2.0,
                    "sentaurus_gain": 2.0,
                }
                for bias in self.biases
            ],
            "knee_metrics": {
                "median_absolute_log_error_dex": 0.0,
                "maximum_absolute_log_error_dex": 0.0,
            },
            "knee_estimators": {
                "vela": {"V_break": -19.5, "V_slope": -19.6},
                "sentaurus": {"V_break": -19.5, "V_slope": -19.6},
            },
            "adjacent_slope_rmse_dex_per_V": 0.0,
        }

    def _write_iv(self, path: Path, *, closure: bool) -> None:
        fields = ["bias_V"]
        if closure:
            fields += [
                "global_continuity_closure_satisfied",
                "global_electron_continuity_closure_ratio",
                "global_hole_continuity_closure_ratio",
            ]
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            for bias in self.biases:
                row = {"bias_V": bias}
                if closure:
                    row.update(
                        {
                            "global_continuity_closure_satisfied": 1,
                            "global_electron_continuity_closure_ratio": 1.0e-4,
                            "global_hole_continuity_closure_ratio": 1.0e-4,
                        }
                    )
                writer.writerow(row)

    def _write_probe(self, path: Path, biases) -> None:
        fields = [
            "bias_V",
            "configuration_fingerprint",
            "source_integral",
            "electron_residual_contributions",
            "hole_residual_contributions",
        ]
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            for bias in biases:
                writer.writerow(
                    {
                        "bias_V": bias,
                        "configuration_fingerprint": "fixture",
                        "source_integral": 1.0,
                        "electron_residual_contributions": -1.0,
                        "hole_residual_contributions": -1.0,
                    }
                )

    def _level_fixture(self, root: Path):
        for name, content in (
            ("mesh.json", "{}\n"),
            ("doping.csv", "node_id,donors_cm3,acceptors_cm3\n"),
            ("materials.json", "{}\n"),
        ):
            (root / name).write_text(content, encoding="utf-8")
        config, render_manifest = render_named_template(
            "pn2d_bv",
            {
                "avalanche_current_support_profile": "element_edge_sg_gss_laux",
                "mesh_file": str(root / "mesh.json"),
                "node_doping_file": str(root / "doping.csv"),
                "materials_file": str(root / "materials.json"),
            },
            allow_absolute_paths=True,
        )
        base_config = root / "base.json"
        render_manifest_path = root / "base.manifest.json"
        write_json(base_config, config)
        write_json(render_manifest_path, render_manifest)
        sentaurus_manifest_path = root / "sentaurus_manifest.json"
        write_json(sentaurus_manifest_path, {"status": "passed"})
        sentaurus_hash = sha256(sentaurus_manifest_path)

        branches = []
        for branch in BRANCHES:
            case = root / branch
            case.mkdir()
            branch_config = json.loads(json.dumps(config))
            if branch == "avalanche_off":
                branch_config["solver"]["impact_ionization"] = {"model": "none"}
            else:
                branch_config["solver"]["impact_ionization"]["coupling_mode"] = (
                    "postprocess_only"
                    if branch == "iic_postprocess"
                    else "self_consistent"
                )
            config_path = case / "simulation.json"
            write_json(config_path, branch_config)
            iv_path = case / "iv.csv"
            self._write_iv(iv_path, closure=True)
            if branch == "avalanche_on":
                self._write_probe(case / "process_probe.csv", self.biases)
            branches.append(
                {
                    "branch": branch,
                    "config": str(config_path),
                    "config_sha256": sha256(config_path),
                    "output_csv": str(iv_path),
                    "output_csv_sha256": sha256(iv_path),
                    "physics_config_sha256": f"physics-{branch}",
                    "returncode": 0,
                    "complete_exact_lattice": True,
                }
            )
        execution = {
            "status": "passed",
            "current_support": {
                "origin": "base_config",
                "current_approximation": "element_edge_sg_gss_laux",
                "source_mapping_mode": "element_vertex_box_measure",
                "cell_reconstructed_midpoint_density": "bernoulli",
            },
            "requested_biases_V": self.biases,
            "base_config": str(base_config),
            "base_config_sha256": sha256(base_config),
            "base_config_manifest": {
                "path": str(render_manifest_path),
                "sha256": sha256(render_manifest_path),
            },
            "sentaurus_manifest_sha256": sentaurus_hash,
            "branches": branches,
        }
        execution_a_path = root / "execution_a.json"
        execution_b_path = root / "execution_b.json"
        write_json(execution_a_path, execution)
        write_json(execution_b_path, execution)

        branch_records = []
        for branch in BRANCHES:
            bias_records = []
            states_dir = root / f"{branch}_states"
            states_dir.mkdir()
            for bias in self.biases:
                state_path = states_dir / f"{abs(bias):g}.csv"
                state_path.write_text(f"{branch},{bias}\n", encoding="utf-8")
                bias_records.append(
                    {
                        "requested_bias_V": bias,
                        "snapshot_tdr": {
                            "path": str(state_path.relative_to(root)),
                            "sha256": sha256(state_path),
                        },
                    }
                )
            branch_records.append(
                {
                    "branch": branch,
                    "requested_biases_V": self.biases,
                    "bias_records": bias_records,
                }
            )
        state = {
            "status": "passed",
            "branch_records": branch_records,
        }
        state_a_path = root / "state_a.json"
        state_b_path = root / "state_b.json"
        write_json(state_a_path, state)
        write_json(state_b_path, state)
        parity_path = root / "parity.json"
        parity = self._parity()
        parity["inputs"] = {
            "vela_on": {
                "path": branches[2]["output_csv"],
                "sha256": branches[2]["output_csv_sha256"],
            },
            "vela_off": {
                "path": branches[0]["output_csv"],
                "sha256": branches[0]["output_csv_sha256"],
            },
            "sentaurus_on": {
                "path": branches[2]["output_csv"],
                "sha256": branches[2]["output_csv_sha256"],
            },
            "sentaurus_off": {
                "path": branches[0]["output_csv"],
                "sha256": branches[0]["output_csv_sha256"],
            },
        }
        write_json(parity_path, parity)
        return {
            "execution_a_path": execution_a_path,
            "execution_b_path": execution_b_path,
            "state_a_path": state_a_path,
            "state_b_path": state_b_path,
            "parity_path": parity_path,
            "render_manifest_path": render_manifest_path,
            "sentaurus_manifest_path": sentaurus_manifest_path,
            "sentaurus_on_csv_path": Path(branches[2]["output_csv"]),
            "sentaurus_off_csv_path": Path(branches[0]["output_csv"]),
        }


if __name__ == "__main__":
    unittest.main()
