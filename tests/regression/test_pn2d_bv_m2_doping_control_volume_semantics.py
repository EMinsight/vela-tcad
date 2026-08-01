from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]


def write_csv(path: Path, fields: list[str], rows: list[list[object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(fields)
        writer.writerows(rows)


class M2DopingControlVolumeSemanticsTest(unittest.TestCase):
    def test_exact_nodal_export_passes_and_records_unobserved_volume(self) -> None:
        with tempfile.TemporaryDirectory(prefix="vela_m2_doping_semantics_") as tmp:
            root = Path(tmp)
            mesh = {
                "nodes": [
                    {"id": 0, "x": 0.75, "y": 0.0},
                    {"id": 1, "x": 1.0, "y": 0.0},
                    {"id": 2, "x": 1.25, "y": 0.0},
                    {"id": 3, "x": 1.0, "y": 0.25},
                ],
                "triangles": [
                    {"id": 0, "node_ids": [0, 1, 3]},
                    {"id": 1, "node_ids": [1, 2, 3]},
                ],
            }
            (root / "mesh.json").write_text(json.dumps(mesh) + "\n", encoding="utf-8")
            doping_rows = [
                [0, 0.0, 1.0e17],
                [1, 5.0e16, 5.0e16],
                [2, 1.0e17, 0.0],
                [3, 5.0e16, 5.0e16],
            ]
            write_csv(root / "doping.csv", ["node_id", "donors_cm3", "acceptors_cm3"], doping_rows)
            (root / "config.json").write_text(json.dumps({
                "mesh_file": "mesh.json",
                "node_doping_file": "doping.csv",
                "solver": {"mobility": {"doping_concentration_basis": "net_doping"}},
            }) + "\n", encoding="utf-8")

            export = root / "sentaurus"
            write_csv(export / "nodes.csv", ["id", "x_um", "y_um"], [
                [node["id"], node["x"], node["y"]] for node in mesh["nodes"]
            ])
            write_csv(export / "elements.csv", ["id", "node0", "node1", "node2", "region", "material"], [
                [0, 0, 1, 3, "R.Si", "Si"],
                [1, 1, 2, 3, "R.Si", "Si"],
            ])
            write_csv(export / "doping.csv", ["node_id", "donors_cm3", "acceptors_cm3"], doping_rows)
            write_csv(export / "fields" / "DopingConcentration_region0.csv", ["node_id", "component0"], [
                [0, -1.0e17], [1, 0.0], [2, 1.0e17], [3, 0.0],
            ])
            write_csv(export / "fields" / "PhosphorusActiveConcentration_region0.csv", ["node_id", "component0"], [
                [row[0], row[1]] for row in doping_rows
            ])
            write_csv(export / "fields" / "BoronActiveConcentration_region0.csv", ["node_id", "component0"], [
                [row[0], row[2]] for row in doping_rows
            ])
            write_csv(root / "modes.csv", [
                "bias_V", "variant", "step_energy_rank", "top_right_node",
                "top_left_node", "step_energy_fraction",
            ], [[-20.0, "sent_qfp_only", 1, 0, 2, 0.9]])

            output = root / "output"
            completed = subprocess.run([
                sys.executable,
                str(REPO / "scripts" / "audit_pn2d_bv_m2_doping_control_volume_semantics.py"),
                "--vela-config", str(root / "config.json"),
                "--sentaurus-export", str(export),
                "--soft-modes", str(root / "modes.csv"),
                "--output-root", str(output),
            ], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)

            self.assertEqual(completed.returncode, 0, completed.stderr)
            result = json.loads((output / "result.json").read_text(encoding="utf-8"))
            self.assertEqual(result["status"], "passed")
            self.assertTrue(result["comparison"]["unordered_triangle_connectivity_exact"])
            self.assertEqual(result["comparison"]["maximum_net_relative_error"], 0.0)
            self.assertFalse(result["control_volume"]["sentaurus_control_volume_semantics_directly_observable"])
            self.assertEqual(result["configuration"]["node_volume_policy"], "barycentric")


if __name__ == "__main__":
    unittest.main()
