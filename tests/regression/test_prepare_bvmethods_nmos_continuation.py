from __future__ import annotations

import importlib.util
from pathlib import Path
import tempfile
import unittest


ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "prepare_bvmethods_continuation",
    ROOT / "scripts" / "prepare_bvmethods_nmos_continuation.py")
assert SPEC and SPEC.loader
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class PrepareBvmethodsContinuationTest(unittest.TestCase):
    def base(self) -> dict:
        return {
            "simulation_type": "dc_sweep",
            "mesh_file": "mesh.json",
            "node_doping_file": "doping.csv",
            "materials_file": "materials.json",
            "contacts": [{"name": "drain", "type": "ohmic", "bias": 0.0}],
            "solver": {
                "method": "newton",
                "mobility": {"model": "masetti_field_lombardi"},
                "impact_ionization": {
                    "model": "van_overstraeten",
                    "coupling_mode": "self_consistent",
                },
                "srh_doping_dependence": {"enabled": True},
            },
            "sweep": {
                "mode": "bv_reverse",
                "voltage_to_current": {"switch_voltage_V": 6.0},
            },
        }

    def test_only_sweep_control_changes(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            state = root / "state.csv"
            state.write_text("node_id,psi,phin,phip\n")
            previous = root / "previous.csv"
            previous.write_text("node_id,psi,phin,phip\n")
            base = self.base()
            result = MODULE.prepare_config(
                base, state, previous, root / "output")
            self.assertEqual(base["solver"], result["solver"])
            self.assertEqual(base["contacts"], result["contacts"])
            self.assertNotIn("voltage_to_current", result["sweep"])
            self.assertTrue(
                result["sweep"]["continuation"]["arclength"]["enabled"])
            self.assertEqual(
                0.0,
                result["sweep"]["continuation"]["arclength"]["state_weight"],
            )
            self.assertEqual(
                0.02,
                result["sweep"]["continuation"]["arclength"]["max_parameter_update"],
            )
            self.assertEqual(
                previous.resolve().as_posix(),
                Path(result["sweep"]["continuation"]["arclength"]
                     ["initial_secant_state_file"]).as_posix(),
            )
            self.assertEqual(
                6.0,
                result["sweep"]["continuation"]["arclength"]
                      ["initial_secant_bias_V"],
            )
            self.assertEqual(6.383727168968036,
                             result["_validation_case"]["sentaurus_reference_BV_V"])

    def test_explicit_state_weight_is_preserved(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            state = root / "state.csv"
            state.write_text("node_id,psi,phin,phip\n")
            previous = root / "previous.csv"
            previous.write_text("node_id,psi,phin,phip\n")
            result = MODULE.prepare_config(
                self.base(), state, previous, root / "output",
                state_weight=6.25e-7)
            self.assertEqual(
                6.25e-7,
                result["sweep"]["continuation"]["arclength"]["state_weight"],
            )

    def test_arclength_source_jacobian_does_not_change_start_solver(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            state = root / "state.csv"
            state.write_text("node_id,psi,phin,phip\n")
            previous = root / "previous.csv"
            previous.write_text("node_id,psi,phin,phip\n")
            result = MODULE.prepare_config(
                self.base(), state, previous, root / "output",
                source_jacobian="finite_difference")
            self.assertNotIn(
                "source_jacobian", result["solver"]["impact_ionization"])
            self.assertEqual(
                "finite_difference",
                result["sweep"]["continuation"]["arclength"]
                      ["source_jacobian"],
            )

    def test_requires_self_consistent_avalanche(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            state = Path(directory) / "state.csv"
            state.write_text("node_id,psi,phin,phip\n")
            previous = Path(directory) / "previous.csv"
            previous.write_text("node_id,psi,phin,phip\n")
            base = self.base()
            base["solver"]["impact_ionization"]["coupling_mode"] = "postprocess_only"
            with self.assertRaisesRegex(ValueError, "self-consistent"):
                MODULE.prepare_config(
                    base, state, previous, Path(directory) / "output")


if __name__ == "__main__":
    unittest.main()
