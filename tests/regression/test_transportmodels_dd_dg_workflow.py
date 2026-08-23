#!/usr/bin/env python3
"""Regression coverage for the TransportModels DD/DG workflow driver."""

from __future__ import annotations

import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "run_transportmodels_dd_dg_workflow.py"
SPEC = importlib.util.spec_from_file_location("transportmodels_workflow", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
WORKFLOW = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(WORKFLOW)


def write_json(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload), encoding="utf-8")


def write_reference(path: Path, biases: list[float]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["bias_V", "current_total"])
        writer.writeheader()
        for index, bias in enumerate(biases):
            writer.writerow({"bias_V": bias, "current_total": 1.0e-12 + index})


class TransportModelsWorkflowTest(unittest.TestCase):
    def make_generated_tree(self, root: Path) -> Path:
        generated = root / "generated"
        common_solver = {
            "statistics": "fermi_dirac",
            "bandgap_narrowing": "old_slotboom",
            "mobility": {
                "model": "masetti_field_lombardi",
                "high_field_driving_force": "quasi_fermi_gradient",
            },
            "srh": {"enabled": True},
        }
        for branch in ("dd", "dg"):
            solver = json.loads(json.dumps(common_solver))
            if branch == "dg":
                solver["electron_quantum_potential"] = {"enabled": True}
            for curve in ("idvg", "idvd"):
                write_json(generated / "vela" / f"simulation_{branch}_{curve}.json", {
                    "mesh_file": "mesh.json",
                    "node_doping_file": "doping.csv",
                    "materials_file": "materials.json",
                    "contacts": [
                        {"name": "source", "bias": 0.0},
                        {"name": "drain", "bias": 0.0},
                        {"name": "gate", "bias": 0.0},
                        {"name": "substrate", "bias": 0.0},
                    ],
                    "solver": solver,
                })
            write_reference(
                generated / "reference_curves" /
                f"transportmodels_sentaurus2022_{branch}_idvg_reference.csv",
                [-1.0 + 0.1 * index for index in range(21)],
            )
            write_reference(
                generated / "reference_curves" /
                f"transportmodels_sentaurus2022_{branch}_idvd_reference.csv",
                [0.05 * index for index in range(21)],
            )
        return generated

    def test_materializes_controlled_delta_and_restart_bias_contract(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            generated = self.make_generated_tree(root)
            manifest = WORKFLOW.materialize(generated, root / "run", ["dd", "dg"])

            self.assertEqual(12, len(manifest["stages"]))
            self.assertTrue(manifest["controlled_delta"][
                "dd_to_dg_only_electron_quantum_potential"])
            stages = {stage["name"]: stage for stage in manifest["stages"]}
            for branch in ("dd", "dg"):
                ramp = stages[f"{branch}_idvg_drain_ramp"]
                self.assertAlmostEqual(0.05, ramp["bias_points"][0])
                self.assertAlmostEqual(1.1, ramp["bias_points"][-1])
                self.assertEqual(22, len(ramp["bias_points"]))
                relax = stages[f"{branch}_idvg_final_bias_relax"]
                self.assertEqual([-1.0], relax["bias_points"])
                self.assertEqual([f"{branch}_idvg_drain_ramp"],
                                 relax["depends_on"])
                curve = stages[f"{branch}_idvg_curve"]
                self.assertAlmostEqual(-0.9, curve["bias_points"][0])
                self.assertEqual(
                    f"{branch}_idvg_final_bias_relax",
                    curve["comparison_seed"]["stage"],
                )
                idvd = stages[f"{branch}_idvd_curve"]
                self.assertGreater(idvd["bias_points"][0], 0.0)
                self.assertEqual(
                    f"{branch}_idvd_equilibrium",
                    idvd["comparison_seed"]["stage"],
                )

            dd = json.loads(Path(stages["dd_idvg_curve"]["config"]).read_text())
            dg = json.loads(Path(stages["dg_idvg_curve"]["config"]).read_text())
            self.assertNotIn("electron_quantum_potential", dd["solver"])
            self.assertTrue(dg["solver"]["electron_quantum_potential"]["enabled"])
            self.assertEqual(
                40,
                dg["solver"]["electron_quantum_potential"][
                    "outer_max_iterations"],
            )
            self.assertEqual(
                "aitken",
                dg["solver"]["electron_quantum_potential"][
                    "outer_acceleration"],
            )
            self.assertEqual(
                "transport_cell_vector",
                dd["solver"]["mobility"]["high_field_gradient_discretization"],
            )
            dd_solver = json.loads(json.dumps(dd["solver"]))
            dg_solver = json.loads(json.dumps(dg["solver"]))
            dg_solver.pop("electron_quantum_potential")
            self.assertEqual(dd_solver, dg_solver)

    def test_rejects_non_quantum_dd_dg_physics_drift(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            generated = self.make_generated_tree(root)
            for curve in ("idvg", "idvd"):
                path = generated / "vela" / f"simulation_dg_{curve}.json"
                deck = json.loads(path.read_text())
                deck["solver"]["statistics"] = "boltzmann"
                write_json(path, deck)
            with self.assertRaisesRegex(ValueError, "beyond electron quantum"):
                WORKFLOW.validate_controlled_delta(generated)

    def test_external_idvg_ramp_state_is_hashed_and_skips_only_completed_prefix(
            self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            generated = self.make_generated_tree(root)
            checkpoint = root / "dg_ramp_1p1.csv"
            checkpoint.write_text("state checkpoint\n", encoding="utf-8")

            manifest = WORKFLOW.materialize(
                generated, root / "run", ["dg"], {"dg": checkpoint})
            stages = {stage["name"]: stage for stage in manifest["stages"]}

            self.assertEqual(4, len(stages))
            self.assertNotIn("dg_idvg_equilibrium", stages)
            self.assertNotIn("dg_idvg_drain_ramp", stages)
            relax = stages["dg_idvg_final_bias_relax"]
            self.assertEqual([], relax["depends_on"])
            self.assertEqual(str(checkpoint.resolve()),
                             relax["initial_state_file"])
            self.assertEqual(
                WORKFLOW.sha256(checkpoint),
                relax["external_initial_state"]["sha256"],
            )
            self.assertEqual(
                ["dg_idvg_final_bias_relax"],
                stages["dg_idvg_curve"]["depends_on"],
            )

    def test_external_idvg_curve_prefix_restarts_remaining_exact_biases(
            self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            generated = self.make_generated_tree(root)
            checkpoint = root / "dg_idvg_m0p8.csv"
            checkpoint.write_text("state checkpoint\n", encoding="utf-8")
            prefix = root / "prefix.csv"
            with prefix.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=["bias_V", "current_total_A_per_um"],
                )
                writer.writeheader()
                for index, bias in enumerate((-1.0, -0.9, -0.8)):
                    writer.writerow({
                        "bias_V": bias,
                        "current_total_A_per_um": 1.0e-12 + index,
                    })

            manifest = WORKFLOW.materialize(
                generated,
                root / "run",
                ["dg"],
                idvg_curve_restarts={
                    "dg": {"bias_V": -0.8, "state": checkpoint,
                           "prefix": prefix},
                },
                quantum_outer_max_iterations=80,
            )
            stages = {stage["name"]: stage for stage in manifest["stages"]}

            self.assertEqual(3, len(stages))
            self.assertNotIn("dg_idvg_final_bias_relax", stages)
            curve = stages["dg_idvg_curve"]
            self.assertAlmostEqual(-0.7, curve["bias_points"][0])
            self.assertEqual([], curve["depends_on"])
            self.assertEqual(str(checkpoint.resolve()),
                             curve["initial_state_file"])
            self.assertEqual(str(prefix.resolve()),
                             curve["comparison_seed"]["external_prefix"])
            deck = json.loads(Path(curve["config"]).read_text())
            self.assertEqual(
                80,
                deck["solver"]["electron_quantum_potential"][
                    "outer_max_iterations"],
            )
            provenance = curve["external_initial_state"]
            self.assertEqual(WORKFLOW.sha256(checkpoint),
                             provenance["sha256"])
            self.assertEqual(WORKFLOW.sha256(prefix),
                             provenance["prefix_sha256"])

            with Path(curve["output_csv"]).open(
                    "w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=["bias_V", "current_total_A_per_um",
                                "converged"],
                )
                writer.writeheader()
                for index, bias in enumerate(curve["bias_points"], start=3):
                    writer.writerow({
                        "bias_V": bias,
                        "current_total_A_per_um": 1.0e-12 + index,
                        "converged": 1,
                    })
            candidate = WORKFLOW.prepare_comparison_candidate(
                curve, stages, root / "run")
            with candidate.open(newline="", encoding="utf-8") as handle:
                candidate_rows = list(csv.DictReader(handle))
            self.assertEqual(21, len(candidate_rows))
            self.assertAlmostEqual(-1.0,
                                   float(candidate_rows[0]["bias_V"]))
            self.assertAlmostEqual(1.0,
                                   float(candidate_rows[-1]["bias_V"]))

    def test_external_idvd_curve_prefix_supports_checkpointed_bridge_bias(
            self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            generated = self.make_generated_tree(root)
            checkpoint = root / "dg_idvd_0p4.csv"
            checkpoint.write_text("state checkpoint\n", encoding="utf-8")
            prefix = root / "idvd_prefix.csv"
            with prefix.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=["bias_V", "current_total_A_per_um"],
                )
                writer.writeheader()
                for index, bias in enumerate(
                        (0.0, 0.05, 0.1, 0.15, 0.2,
                         0.25, 0.3, 0.35, 0.4)):
                    writer.writerow({
                        "bias_V": bias,
                        "current_total_A_per_um": 1.0e-6 + index,
                    })

            manifest = WORKFLOW.materialize(
                generated,
                root / "run",
                ["dg"],
                idvd_curve_restarts={
                    "dg": {"bias_V": 0.4, "state": checkpoint,
                           "prefix": prefix},
                },
                idvd_bridge_biases={"dg": [0.425]},
                quantum_outer_max_iterations=80,
                quantum_outer_acceleration="none",
                quantum_outer_relaxation=1.2,
            )

            self.assertEqual(1, len(manifest["stages"]))
            curve = manifest["stages"][0]
            self.assertEqual("dg_idvd_curve", curve["name"])
            self.assertAlmostEqual(0.425, curve["bias_points"][0])
            self.assertAlmostEqual(0.45, curve["bias_points"][1])
            self.assertAlmostEqual(1.0, curve["bias_points"][-1])
            self.assertEqual([], curve["depends_on"])
            self.assertEqual(str(checkpoint.resolve()),
                             curve["initial_state_file"])
            self.assertEqual(str(prefix.resolve()),
                             curve["comparison_seed"]["external_prefix"])
            deck = json.loads(Path(curve["config"]).read_text())
            self.assertEqual(80, deck["solver"]["electron_quantum_potential"]
                             ["outer_max_iterations"])
            self.assertEqual("none", deck["solver"]["electron_quantum_potential"]
                             ["outer_acceleration"])
            self.assertEqual(1.2, deck["solver"]["electron_quantum_potential"]
                             ["outer_relaxation"])
            provenance = curve["external_initial_state"]
            self.assertEqual("completed_idvd_curve_prefix",
                             provenance["role"])
            self.assertEqual(WORKFLOW.sha256(checkpoint),
                             provenance["sha256"])
            self.assertEqual(WORKFLOW.sha256(prefix),
                             provenance["prefix_sha256"])

            with Path(curve["output_csv"]).open(
                    "w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=["bias_V", "current_total_A_per_um",
                                "converged"],
                )
                writer.writeheader()
                for index, bias in enumerate(curve["bias_points"], start=9):
                    writer.writerow({
                        "bias_V": bias,
                        "current_total_A_per_um": 1.0e-6 + index,
                        "converged": 1,
                    })
            candidate = WORKFLOW.prepare_comparison_candidate(
                curve, {curve["name"]: curve}, root / "run")
            with candidate.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(21, len(rows))
            self.assertNotIn(0.425, [float(row["bias_V"]) for row in rows])

    def test_external_idvd_curve_prefix_rejects_reference_lattice_gap(
            self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            generated = self.make_generated_tree(root)
            checkpoint = root / "dg_idvd_0p4.csv"
            checkpoint.write_text("state checkpoint\n", encoding="utf-8")
            prefix = root / "idvd_prefix_with_gap.csv"
            with prefix.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=["bias_V", "current_total_A_per_um"],
                )
                writer.writeheader()
                for bias in (0.0, 0.05, 0.1, 0.2, 0.25, 0.3, 0.35, 0.4):
                    writer.writerow({
                        "bias_V": bias,
                        "current_total_A_per_um": 1.0e-6,
                    })

            with self.assertRaisesRegex(ValueError, "prefix lattice"):
                WORKFLOW.materialize(
                    generated,
                    root / "run",
                    ["dg"],
                    idvd_curve_restarts={
                        "dg": {"bias_V": 0.4, "state": checkpoint,
                               "prefix": prefix},
                    },
                )


if __name__ == "__main__":
    unittest.main()
