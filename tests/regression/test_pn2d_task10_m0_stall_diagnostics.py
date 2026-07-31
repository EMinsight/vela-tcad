import csv
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def load_script(module_name: str, script_name: str):
    spec = importlib.util.spec_from_file_location(
        module_name, ROOT / "scripts" / script_name
    )
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


ownership = load_script(
    "pn2d_task10_m0_junction_ownership",
    "diagnose_pn2d_task10_m0_junction_ownership.py",
)
failure = load_script(
    "pn2d_task10_m0_first_failure",
    "diagnose_pn2d_task10_m0_first_failure.py",
)


class TestM0StallDiagnostics(unittest.TestCase):
    def test_balanced_half_preserves_total_and_zero_net_junction_doping(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            mesh = root / "mesh.json"
            mesh.write_text(
                json.dumps(
                    {
                        "nodes": [
                            {"id": 0, "x": 0.0, "y": 0.0},
                            {"id": 1, "x": 1.0, "y": 0.0},
                            {"id": 2, "x": 1.0, "y": 1.0},
                        ],
                        "cells": [{"node_ids": [0, 1, 2]}],
                    }
                ),
                encoding="utf-8",
            )
            source = root / "source.csv"
            with source.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle)
                writer.writerow(["node_id", "donors_cm3", "acceptors_cm3"])
                writer.writerows(
                    [
                        [0, 0.0, 1.0e17],
                        [1, 1.0e17, 0.0],
                        [2, 1.0e17, 0.0],
                    ]
                )
            balanced = root / "balanced.csv"
            junction_ids = ownership.junction_nodes(mesh, 1.0)
            ownership.variant_doping(
                source,
                balanced,
                junction_ids,
                donor_fraction=0.5,
                acceptor_fraction=0.5,
                concentration_cm3=1.0e17,
            )
            metrics = ownership.doping_metrics(
                balanced,
                ownership.nodal_control_areas_um2(mesh),
                junction_ids,
            )
            self.assertEqual(metrics["junction_net_doping_cm3"], [0.0, 0.0])
            self.assertEqual(
                metrics["junction_total_impurity_cm3"], [1.0e17, 1.0e17]
            )

    def test_first_rejected_attempt_requires_all_state_snapshots(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            states = [root / f"{name}.csv" for name in ("parent", "initial", "final")]
            for state in states:
                state.write_text("node_id,psi,phin,phip\n0,0,0,0\n", encoding="utf-8")
            attempts = root / "attempts.csv"
            with attempts.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.DictWriter(
                    handle,
                    fieldnames=[
                        "status",
                        "rejected_parent_state_file",
                        "rejected_initial_state_file",
                        "rejected_final_state_file",
                    ],
                )
                writer.writeheader()
                writer.writerow(
                    {
                        "status": "accepted",
                        "rejected_parent_state_file": "",
                        "rejected_initial_state_file": "",
                        "rejected_final_state_file": "",
                    }
                )
                writer.writerow(
                    {
                        "status": "rejected",
                        "rejected_parent_state_file": str(states[0]),
                        "rejected_initial_state_file": str(states[1]),
                        "rejected_final_state_file": str(states[2]),
                    }
                )
            selected = failure.first_rejected_attempt(attempts)
            self.assertEqual(selected["status"], "rejected")
            self.assertEqual(
                Path(selected["rejected_initial_state_file"]), states[1]
            )

    def test_probe_config_separates_off_iic_and_self_consistent_feedback(self):
        base_path = Path("base.json").resolve()
        base = {
            "mesh_file": "mesh.json",
            "node_doping_file": "doping.csv",
            "materials_file": "materials.json",
            "contacts": [
                {"name": "Anode", "bias": 0.0},
                {"name": "Cathode", "bias": 0.0},
            ],
            "solver": {
                "impact_ionization": {
                    "model": "van_overstraeten",
                    "coupling_mode": "self_consistent",
                }
            },
            "sweep": {},
        }
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            configs = {
                branch: failure.probe_config(
                    base_path,
                    base,
                    branch,
                    "residual_probe",
                    -17.15,
                    root / "fields",
                    root / f"{branch}.csv",
                )
                for branch in failure.BRANCHES
            }
        self.assertEqual(
            configs["avalanche_off"]["solver"]["impact_ionization"]["model"],
            "none",
        )
        self.assertEqual(
            configs["iic_postprocess"]["solver"]["impact_ionization"][
                "coupling_mode"
            ],
            "postprocess_only",
        )
        self.assertEqual(
            configs["avalanche_on"]["solver"]["impact_ionization"][
                "coupling_mode"
            ],
            "self_consistent",
        )
        for config in configs.values():
            anode = next(
                contact
                for contact in config["contacts"]
                if contact["name"] == "Anode"
            )
            self.assertEqual(anode["bias"], -17.15)
            self.assertNotIn("sweep", config)


if __name__ == "__main__":
    unittest.main()
