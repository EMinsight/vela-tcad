#!/usr/bin/env python3
"""Regression coverage for the Sentaurus SingleDevice reference contract."""

from __future__ import annotations

import csv
import json
import sys
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
FIXTURE = REPO / "reference_tcad" / "singledevice_sentaurus2018"
if str(REPO / "scripts") not in sys.path:
    sys.path.insert(0, str(REPO / "scripts"))

from sentaurus_execution_ir import build_execution_ir  # noqa: E402
from sentaurus_import import apply_solver_physics, parse_cmd, sentaurus_models  # noqa: E402
from run_singledevice_sentaurus_vm import prepare_bundle, remote_commands  # noqa: E402


class SingleDeviceReferenceFixtureTest(unittest.TestCase):
    def test_vm_runner_preprocesses_a_complete_bundle(self) -> None:
        import tempfile

        with tempfile.TemporaryDirectory() as directory:
            bundle = Path(directory) / "bundle"
            prepare_bundle(FIXTURE / "source", bundle)
            self.assertEqual(
                {"singledevice_sde.cmd", "singledevice_sdevice.cmd", "Silicon.par", "sdevice.par"},
                {path.name for path in bundle.iterdir()},
            )
            self.assertNotIn("@node@", (bundle / "singledevice_sde.cmd").read_text())
            sdevice = (bundle / "singledevice_sdevice.cmd").read_text()
            self.assertNotIn("@tdr@", sdevice)
            self.assertIn('Grid      = "n2_msh.tdr"', sdevice)

        commands = remote_commands("/tmp/singledevice")
        self.assertIn("sde -e -l", commands[0])
        self.assertIn("tdx -mtt", commands[1])
        self.assertIn("sdevice singledevice_sdevice.cmd", commands[2])

    def test_reference_curves_are_two_complete_matching_bias_grids(self) -> None:
        curves = []
        for name in (
                "singledevice_idvg_lin_reference.csv",
                "singledevice_idvg_sat_reference.csv"):
            with (FIXTURE / name).open(newline="") as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(21, len(rows))
            self.assertAlmostEqual(-0.5, float(rows[0]["bias_V"]))
            self.assertAlmostEqual(2.2, float(rows[-1]["bias_V"]))
            self.assertTrue(all(float(row["current_total"]) > 0.0 for row in rows))
            curves.append([float(row["bias_V"]) for row in rows])
        self.assertEqual(curves[0], curves[1])

    def test_manifest_preserves_branch_and_mesh_contract(self) -> None:
        manifest = json.loads(
            (FIXTURE / "singledevice_sentaurus2018_reference.json").read_text())
        self.assertEqual("O-2018.06-SP2", manifest["sentaurus_version"])
        self.assertEqual(3584, manifest["mesh_inventory"]["vertices"])
        self.assertEqual(6972, manifest["mesh_inventory"]["triangles"])
        self.assertEqual([0.1, 1.1], [
            item["drain_voltage_V"] for item in manifest["simulations"]])
        self.assertEqual(
            "same saved equilibrium state as idvg_lin",
            manifest["simulations"][1]["initial_state"],
        )

    def test_import_config_keeps_two_quantum_corrected_branches(self) -> None:
        config = json.loads(
            (FIXTURE / "singledevice_import_config.json").read_text())
        self.assertEqual(["idvg_lin", "idvg_sat"], [
            sim["name"] for sim in config["simulations"]])
        self.assertTrue(
            config["vela_solver"]["electron_quantum_potential"]["enabled"])
        quantum = config["vela_solver"]["electron_quantum_potential"]
        self.assertEqual(1.0618016171622988,
                         quantum["effective_mass_ratio"])
        self.assertFalse(quantum["include_insulators"])
        self.assertEqual(0.42, quantum["insulator_effective_mass_ratio"])
        self.assertEqual("homogeneous_neumann", quantum["interface_boundary"])
        self.assertEqual(0.5, quantum["theta"])
        self.assertEqual(0.5, quantum["conduction_band_narrowing_fraction"])
        drain_biases = []
        for sim in config["simulations"]:
            drain = next(
                item for item in sim["vela_contact_overrides"]
                if item["name"] == "drain")
            drain_biases.append(drain["bias"])
        self.assertEqual([0.1, 1.1], drain_biases)

    def test_custom_materials_use_unit_scaling_native_units(self) -> None:
        materials = json.loads(
            (FIXTURE / "vela" / "materials_sentaurus2018.json").read_text())
        poly = next(item for item in materials["materials"]
                    if item["name"] == "PolySilicon")
        oxide = next(item for item in materials["materials"]
                     if item["name"] == "SiO2")
        self.assertEqual(1.0e10, poly["ni"])
        self.assertEqual(1417.0, poly["mun"])
        self.assertEqual(2.8e19, poly["Nc_m3"])
        self.assertEqual(0.9, oxide["electron_affinity_eV"])

    def test_sdevice_frontend_maps_complete_singledevice_physics(self) -> None:
        cmd = FIXTURE / "source" / "singledevice_sdevice.cmd"
        variables = {
            "previous": "1", "tdr": "n2_msh.tdr", "tdrdat": "singledevice.tdr",
            "parameter": "sdevice.par", "plot": "singledevice.plt",
            "log": "singledevice.log", "node": "2",
        }
        summary = parse_cmd(cmd, variables)
        models = sentaurus_models(summary)
        self.assertIn("eQuantumPotential", models)
        self.assertIn("Enormal", models)
        self.assertIn("eHighFieldsaturation", models)
        self.assertIn("hHighFieldsaturation", models)

        deck = {"solver": {}}
        apply_solver_physics(deck, summary, {"name": "idvg", "kind": "iv"})
        self.assertEqual(
            "masetti_field_lombardi", deck["solver"]["mobility"]["model"])
        self.assertEqual("old_slotboom", deck["solver"]["bandgap_narrowing"])
        self.assertTrue(deck["solver"]["electron_quantum_potential"]["enabled"])
        self.assertAlmostEqual(
            3.6, deck["solver"]["electron_quantum_potential"]["gamma"])
        srh = deck["solver"]["srh_doping_dependence"]
        self.assertTrue(srh["temperature_dependence"])
        self.assertAlmostEqual(3.0e-8, srh["electron"]["tau_max_s"])
        self.assertAlmostEqual(3.0e-6, srh["hole"]["tau_max_s"])

        ir = build_execution_ir(summary, str(cmd), models)
        self.assertEqual([], ir["unsupported"])
        stages = ir["stages"]
        load = next(stage for stage in stages if stage["phase"] == "load")
        save = next(stage for stage in stages if stage["phase"] == "save")
        self.assertEqual([save["index"]], load["depends_on"])
        self.assertEqual(
            [0.1, 2.2, 1.1, 2.2],
            [stage["goal_voltage"] for stage in stages if stage["phase"] == "sweep"],
        )


if __name__ == "__main__":
    unittest.main()
