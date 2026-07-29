#!/usr/bin/env python3

from __future__ import annotations

import re
import tempfile
import unittest
from pathlib import Path

from scripts.analyze_pn2d_bv_process_matrix_pair import (
    process_matrix_shape,
    source_closure,
)
from scripts.run_pn2d_bv_process_matrix_vm import (
    BRANCHES,
    CURRENT_PLOT_GENERATION_TO_A_PER_UM,
    bias_tag,
    currentplot_aggregates,
    exact_currentplot_rows,
    make_branch_deck,
    remote_command,
    remote_shell_text,
)


TEMPLATE = """\
File {
  Grid = "pn2d_minimal6.tdr"
  Doping = "pn2d_minimal6.tdr"
  Parameter = "models.par"
  Plot = "runtime_element_avalanche_probe_default.tdr"
  Current = "runtime_element_avalanche_probe_default.plt"
  Output    = "runtime_element_avalanche_probe_default"
}
Electrode {
  { Name="Anode" Voltage=0.0 }
  { Name="Cathode" Voltage=0.0 }
}
Physics {
  Mobility(DopingDependence HighFieldSaturation)
  Recombination(
    SRH
    Avalanche(VanOverstraeten)
  )
  EffectiveIntrinsicDensity(OldSlotboom)
}
Plot {
  Potential
  hAlphaAvalanche
}
CurrentPlot {
  eAvalancheGeneration(
    Integrate(Name="eAvalancheIntegral" Semiconductor)
  )
  hAvalancheGeneration(
    Integrate(Name="hAvalancheIntegral" Semiconductor)
  )
  AvalancheGeneration(
    Integrate(Name="AvalancheIntegral" Semiconductor)
  )
  Tcl (tcl = "source runtime_element_avalanche_probe.tcl")
}
Math {
  Extrapolate
  CurrentPlot(IntegrationUnit=um Digits=12)
}
Solve {
  Coupled { Poisson Electron Hole }
}
"""


class BVProcessMatrixVMTest(unittest.TestCase):
    def test_pair_shape_accepts_any_common_exact_lattice(self) -> None:
        manifest = {
            "branch_records": [
                {
                    "branch": branch,
                    "requested_biases_V": [-18.0, -18.5, -19.0],
                    "bias_records": [{}, {}, {}],
                }
                for branch in BRANCHES
            ]
        }
        shape = process_matrix_shape(manifest)
        self.assertTrue(shape["required_branches_present"])
        self.assertTrue(shape["common_bias_lattice"])
        self.assertEqual(shape["snapshot_count"], 12)
        self.assertEqual(shape["expected_snapshot_count"], 12)

        manifest["branch_records"][1]["requested_biases_V"] = [-18.0, -19.0]
        self.assertFalse(process_matrix_shape(manifest)["common_bias_lattice"])

    def test_source_closure_reports_missing_provenance(self) -> None:
        base = {
            "branch": "avalanche_off",
            "requested_bias_V": -20.0,
            "actual_bias_V": -20.0,
            "carrier": "total",
            "quantity": "integrated_source",
            "unit": "A/um",
            "value": 1.0,
            "source": {"file": "run.out", "dataset": "probe", "index": 0},
        }
        native = dict(base, provenance="native")
        replay = dict(base, provenance="operator_replay")
        complete = source_closure({"aggregate_records": [native, replay]})
        self.assertEqual(complete["compared_records"], 1)
        self.assertEqual(complete["incomplete_records"], 0)

        incomplete = source_closure({"aggregate_records": [native]})
        self.assertEqual(incomplete["compared_records"], 0)
        self.assertEqual(incomplete["incomplete_records"], 1)

    def test_declared_branches_and_iic_controls(self) -> None:
        biases = (-10.0, -19.95)
        decks = {
            branch: make_branch_deck(TEMPLATE, branch, biases)
            for branch in BRANCHES
        }
        self.assertNotIn("Avalanche(VanOverstraeten", decks["avalanche_off"])
        self.assertIn(
            "Avalanche(VanOverstraeten GradQuasiFermi)",
            decks["iic_postprocess"],
        )
        self.assertIn("ComputeIonizationIntegrals", decks["iic_postprocess"])
        self.assertIn("AvalPostProcessing", decks["iic_postprocess"])
        self.assertNotIn("AvalPostProcessing", decks["avalanche_on"])
        self.assertNotIn("AvalDerivatives", decks["avalanche_on"])
        self.assertIn(
            "AvalDerivatives",
            decks["avalanche_on_aval_derivatives"],
        )

    def test_every_target_has_explicit_goal_and_unique_snapshot(self) -> None:
        biases = (-10.0, -19.0, -19.95)
        deck = make_branch_deck(TEMPLATE, "avalanche_on", biases)
        goals = re.findall(r'Goal \{ Name="Anode" Voltage=([^ ]+) \}', deck)
        prefixes = re.findall(r'Plot\(FilePrefix="([^"]+)"', deck)
        self.assertEqual([float(value) for value in goals], list(biases))
        self.assertEqual(
            prefixes,
            [bias_tag(index, bias) for index, bias in enumerate(biases)],
        )
        self.assertEqual(len(set(prefixes)), len(biases))
        self.assertEqual(deck.count("CurrentPlot(Time=(1))"), len(biases))

    def test_grad_qf_currentplot_units_maxima_and_failure_diagnostics(self) -> None:
        deck = make_branch_deck(TEMPLATE, "avalanche_on", (-19.95,))
        self.assertIn(
            "Avalanche(VanOverstraeten GradQuasiFermi)",
            deck,
        )
        self.assertNotIn(
            "Avalanche(VanOverstraeten ElectricField)",
            deck,
        )
        self.assertIn("CurrentPlot(IntegrationUnit=um Digits=12)", deck)
        self.assertIn(
            "ImpactIonization(Maximum(Semiconductor Coordinates))",
            deck,
        )
        self.assertIn(
            "ElectricField(Maximum(Semiconductor Coordinates))",
            deck,
        )
        self.assertIn("CNormPrint", deck)
        self.assertIn("NewtonPlot(Error MinError Residual)", deck)
        self.assertIn('NewtonPlot = "newton_avalanche_on_%d_%d_des.tdr"', deck)

    def test_iic_and_on_differ_only_by_declared_feedback_controls(self) -> None:
        iic = make_branch_deck(TEMPLATE, "iic_postprocess", (-19.95,))
        on = make_branch_deck(TEMPLATE, "avalanche_on", (-19.95,))

        def normalize(text: str) -> str:
            text = re.sub(r"pn2d_bv_process_[A-Za-z_]+", "STEM", text)
            text = re.sub(r'\n\s*NewtonPlot\s*=\s*"[^"]+"', "", text)
            for line in (
                "  ComputeIonizationIntegrals\n",
                "  AvalPostProcessing\n",
                "  CNormPrint\n",
                "  NewtonPlot(Error MinError Residual)\n",
            ):
                text = text.replace(line, "")
            return text

        self.assertEqual(normalize(iic), normalize(on))

    def test_remote_root_is_validated_and_commands_remain_argv(self) -> None:
        command = remote_command(
            "/home/tcad/codex/wp2-a/avalanche_on",
            ["sdevice", "deck.cmd"],
        )
        self.assertEqual(command[:3], ["cd", "/home/tcad/codex/wp2-a/avalanche_on", "&&"])
        self.assertEqual(command[-2:], ["sdevice", "deck.cmd"])
        rendered = remote_shell_text(
            "/home/tcad/codex/wp2-a/avalanche_on",
            ["sdevice", "deck.cmd"],
        )
        self.assertEqual(
            rendered,
            "cd /home/tcad/codex/wp2-a/avalanche_on && sdevice deck.cmd",
        )
        with self.assertRaises(ValueError):
            remote_command("../../unsafe", ["sdevice", "deck.cmd"])

    def test_currentplot_selection_rejects_nearest_bias(self) -> None:
        text = """\
DF-ISE text
Info { datasets = ["time" "Anode OuterVoltage" "Anode TotalCurrent"] }
Data {
0 -19.95 1e-12
}
"""
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "run.plt"
            path.write_text(text, encoding="ascii")
            rows = exact_currentplot_rows(path, (-19.95,))
            self.assertEqual(rows[0]["actual_bias_V"], -19.95)
            with self.assertRaisesRegex(ValueError, "nearest_bias_substitution"):
                exact_currentplot_rows(path, (-19.9,))

    def test_currentplot_generation_is_converted_to_canonical_current(self) -> None:
        rows = [
            {
                "requested_bias_V": -20.0,
                "actual_bias_V": -20.0,
                "IntegrAvalancheIntegral AvalancheGeneration": 2.0,
            }
        ]
        records = currentplot_aggregates("iic_postprocess", rows)
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0]["unit"], "A/um")
        self.assertEqual(
            records[0]["value"],
            2.0 * CURRENT_PLOT_GENERATION_TO_A_PER_UM,
        )


if __name__ == "__main__":
    unittest.main()
