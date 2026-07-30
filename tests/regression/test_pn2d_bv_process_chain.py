from __future__ import annotations

import random
import tempfile
import unittest
from pathlib import Path

from scripts.analyze_pn2d_bv_process_chain import (
    CHAIN_SCHEMA,
    STAGES,
    analyze,
)


def fixture(simulator: str) -> dict:
    records = []
    for branch in ("avalanche_off", "iic_postprocess", "avalanche_on"):
        for bias in (-19.8, -19.9):
            for stage_index, stage in enumerate(STAGES):
                for support_index, support in enumerate(("cell:1", "cell:2")):
                    value = 10.0 + stage_index + support_index
                    records.append(
                        {
                            "simulator": simulator,
                            "branch": branch,
                            "bias_V": bias,
                            "stage": stage,
                            "quantity": stage,
                            "carrier": "electron",
                            "support_kind": "cell",
                            "support_key": support,
                            "values": [value, value * 0.5] if stage in {"drive", "current"} else [value],
                            "unit": "1",
                            "coordinates_um": [float(support_index), 0.0],
                            "provenance": "native" if simulator == "sentaurus" else "solver_used",
                        }
                    )
    return {
        "schema": CHAIN_SCHEMA,
        "simulator": simulator,
        "records": records,
        "closures": [
            {
                "branch": branch,
                "bias_V": bias,
                "source_native": 1.0,
                "source_reintegrated": 1.0,
                "terminal_source": 1.0,
                "terminal_current": 1.0,
            }
            for branch in ("avalanche_off", "iic_postprocess", "avalanche_on")
            for bias in (-19.8, -19.9)
        ],
        "newton_updates": [],
    }


def inject(dataset: dict, stage: str, *, one_bias: bool = False) -> None:
    for row in dataset["records"]:
        if (
            row["branch"] == "avalanche_on"
            and row["stage"] == stage
            and row["support_key"] == "cell:1"
            and (not one_bias or row["bias_V"] == -19.8)
        ):
            row["values"] = [value * 2.0 for value in row["values"]]


class ProcessChainTests(unittest.TestCase):
    def test_recovers_injected_error_at_every_stage(self) -> None:
        for stage in STAGES:
            with self.subTest(stage=stage), tempfile.TemporaryDirectory() as tmp:
                sentaurus = fixture("sentaurus")
                vela = fixture("vela")
                inject(vela, stage)
                report = analyze(sentaurus, vela, Path(tmp))
                self.assertEqual(report["causal_stage"], stage)
                self.assertEqual(report["adjacent_biases_V"], [-19.8, -19.9])
                self.assertNotEqual(report["outcome"], "insufficient_observation")

    def test_row_order_and_unrelated_tail_do_not_change_first_stage(self) -> None:
        with tempfile.TemporaryDirectory() as tmp_a, tempfile.TemporaryDirectory() as tmp_b:
            sentaurus = fixture("sentaurus")
            vela = fixture("vela")
            inject(vela, "mobility")
            first = analyze(sentaurus, vela, Path(tmp_a))
            for dataset in (sentaurus, vela):
                random.Random(17).shuffle(dataset["records"])
                dataset["records"].append(
                    {
                        "simulator": dataset["simulator"],
                        "branch": "avalanche_on",
                        "bias_V": -19.8,
                        "stage": "state",
                        "quantity": "state_tail",
                        "carrier": "electron",
                        "support_kind": "cell",
                        "support_key": "tail",
                        "values": [1.0e-300 if dataset["simulator"] == "sentaurus" else 1.0e-250],
                        "unit": "1",
                        "coordinates_um": [99.0, 99.0],
                        "provenance": "native",
                    }
                )
            second = analyze(sentaurus, vela, Path(tmp_b))
            self.assertEqual(first["causal_stage"], "mobility")
            self.assertEqual(second["causal_stage"], "mobility")

    def test_requires_two_adjacent_biases(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            sentaurus = fixture("sentaurus")
            vela = fixture("vela")
            inject(vela, "generation", one_bias=True)
            report = analyze(sentaurus, vela, Path(tmp))
            self.assertEqual(report["outcome"], "insufficient_observation")
            self.assertIsNone(report["causal_stage"])

    def test_missing_simulator_fails_closed_and_writes_contract(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            report = analyze(fixture("sentaurus"), None, Path(tmp))
            self.assertEqual(report["outcome"], "insufficient_observation")
            for name in (
                "stage_summary.csv",
                "support_summary.csv",
                "hotspot_chain.csv",
                "first_departure.json",
                "source_terminal_closure.csv",
                "newton_first_update.csv",
                "process_chain.svg",
                "hotspot.svg",
                "acceptance.json",
            ):
                self.assertTrue((Path(tmp) / name).is_file(), name)


if __name__ == "__main__":
    unittest.main()
