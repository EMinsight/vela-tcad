#!/usr/bin/env python3
"""Regression coverage for the PN2D minimal6 Sentaurus topology gate."""

from __future__ import annotations

import csv
import importlib.util
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
TOPOLOGY_FIXTURE = (
    REPO
    / "reference_tcad"
    / "pn2d_sentaurus2018_minimal6"
    / "source"
    / "minimal6_topologies.json"
)

TOPOLOGY_SPEC = importlib.util.spec_from_file_location(
    "pn2d_minimal6_topology",
    REPO / "scripts" / "pn2d_minimal6_topology.py",
)
assert TOPOLOGY_SPEC is not None and TOPOLOGY_SPEC.loader is not None
topology_module = importlib.util.module_from_spec(TOPOLOGY_SPEC)
sys.modules[TOPOLOGY_SPEC.name] = topology_module
TOPOLOGY_SPEC.loader.exec_module(topology_module)

GATE_SPEC = importlib.util.spec_from_file_location(
    "run_pn2d_minimal6_sentaurus_gate",
    REPO / "scripts" / "run_pn2d_minimal6_sentaurus_gate.py",
)
assert GATE_SPEC is not None and GATE_SPEC.loader is not None
gate = importlib.util.module_from_spec(GATE_SPEC)
sys.modules[GATE_SPEC.name] = gate
GATE_SPEC.loader.exec_module(gate)


def write_csv(path: Path, fieldnames: list[str], rows: list[list[object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(fieldnames)
        writer.writerows(rows)


def write_neutral_export(root: Path, topology: object) -> None:
    nodes = sorted(topology.nodes.items())
    write_csv(
        root / "nodes.csv",
        ["id", "x_um", "y_um"],
        [[node_id - 1, point[0], point[1]] for node_id, point in nodes],
    )
    write_csv(
        root / "elements.csv",
        ["id", "node0", "node1", "node2", "region", "material"],
        [
            [index, *(node_id - 1 for node_id in triangle), "R.Si", "Si"]
            for index, triangle in enumerate(topology.triangles)
        ],
    )
    write_csv(
        root / "contacts.csv",
        ["name", "node_ids", "region"],
        [
            [name, ";".join(str(node_id - 1) for node_id in edge), "R.Si"]
            for name, edge in topology.contacts.items()
        ],
    )
    write_csv(
        root / "doping.csv",
        ["node_id", "donors_cm3", "acceptors_cm3"],
        [
            [
                node_id - 1,
                topology.donors_cm3[node_id],
                topology.acceptors_cm3[node_id],
            ]
            for node_id, _ in nodes
        ],
    )
    (root / "metadata.json").write_text(
        json.dumps(
            {
                "vertex_count": 6,
                "regions": [
                    {
                        "name": "R.Si",
                        "material": "Si",
                        "type": 0,
                        "triangles": 4,
                        "edges": 9,
                    },
                    {"name": "Cathode", "type": 1, "triangles": 0, "edges": 1},
                    {"name": "Anode", "type": 1, "triangles": 0, "edges": 1},
                ],
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )


class PN2DMinimal6SentaurusGateTest(unittest.TestCase):
    def test_prepare_gate_converts_dfise_to_tdr_without_remeshing(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            manifest = gate.prepare_gate(
                topology_ids=("sketch", "mirror"),
                run_id="minimal6_gate_test",
                output_dir=Path(tmp),
                ssh_target="sentaurus",
            )

            self.assertEqual(manifest["schema"], "vela.pn2d_minimal6_sentaurus_gate.v1")
            self.assertEqual(manifest["doping_relative_tolerance"], 1.0e-15)
            self.assertEqual(
                [run["topology_id"] for run in manifest["runs"]],
                ["sketch", "mirror"],
            )
            for run in manifest["runs"]:
                command_text = " ".join(run["remote_commands"]).lower()
                self.assertNotRegex(command_text, r"(?:^|&&\s*)sde\s")
                self.assertEqual(len(run["remote_commands"]), 2)
                self.assertIn(
                    "tdx -d pn2d_minimal6.grd pn2d_minimal6.dat "
                    "pn2d_minimal6.tdr",
                    run["remote_commands"][0],
                )
                self.assertIn(
                    "sdevice pn2d_minimal6_gate_sdevice.cmd",
                    run["remote_commands"][1],
                )
                self.assertEqual(
                    set(run["staged_files"]),
                    {
                        "pn2d_minimal6.grd",
                        "pn2d_minimal6.dat",
                        "pn2d_minimal6_gate_sdevice.cmd",
                        "models.par",
                    },
                )
                self.assertEqual(set(run["file_sha256"]), set(run["staged_files"]))
                for digest in run["file_sha256"].values():
                    self.assertRegex(digest, r"^[0-9a-f]{64}$")
                bundle = Path(run["bundle_dir"])
                self.assertEqual(
                    {path.name for path in bundle.iterdir() if path.is_file()},
                    set(run["staged_files"]),
                )
                deck_path = bundle / "pn2d_minimal6_gate_sdevice.cmd"
                self.assertFalse(deck_path.read_bytes().startswith(b"\xef\xbb\xbf"))
                deck = deck_path.read_text(encoding="utf-8")
                self.assertIn('Grid      = "pn2d_minimal6.tdr"', deck)
                self.assertIn('Doping    = "pn2d_minimal6.tdr"', deck)

    def test_live_argv_uses_argument_arrays(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            manifest = gate.prepare_gate(
                topology_ids=("sketch",),
                run_id="minimal6_gate_argv",
                output_dir=Path(tmp),
                ssh_target="sentaurus",
            )
            commands = gate.build_live_argv(
                manifest["runs"][0],
                ssh_bin="ssh-test",
                scp_bin="scp-test",
                ssh_target="sentaurus",
            )

        self.assertTrue(commands)
        self.assertTrue(all(isinstance(command, list) for command in commands))
        self.assertEqual(commands[0][:2], ["ssh-test", "sentaurus"])
        self.assertTrue(any(command[0] == "scp-test" for command in commands))
        self.assertFalse(any(len(command) == 1 for command in commands))

    def test_live_failure_attempts_to_recover_all_remote_artifacts(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = gate.prepare_gate(
                topology_ids=("sketch",),
                run_id="minimal6_gate_failure",
                output_dir=root,
                ssh_target="sentaurus",
            )
            importer = root / "sentaurus_import.exe"
            importer.touch()
            calls: list[list[str]] = []
            original_run_checked = gate.run_checked

            def fail_sdevice_and_missing_tdr(argv: list[str]) -> None:
                calls.append(list(argv))
                command = " ".join(argv)
                if "sdevice pn2d_minimal6_gate_sdevice.cmd" in command:
                    raise subprocess.CalledProcessError(1, argv)
                if "pn2d_minimal6_gate_des.tdr" in command:
                    raise subprocess.CalledProcessError(1, argv)

            gate.run_checked = fail_sdevice_and_missing_tdr
            try:
                with self.assertRaises(subprocess.CalledProcessError):
                    gate.run_live(
                        manifest,
                        ssh_bin="ssh-test",
                        scp_bin="scp-test",
                        importer=importer,
                    )
            finally:
                gate.run_checked = original_run_checked

        recovered = " ".join(" ".join(argv) for argv in calls)
        for name in (
            "pn2d_minimal6_gate_des.log",
            "pn2d_minimal6.grd",
            "pn2d_minimal6.dat",
            "run_pn2d_minimal6_gate.out",
        ):
            self.assertIn(name, recovered)
        run = manifest["runs"][0]
        self.assertEqual(run["status"], "failed")
        self.assertIn("artifact_recovery_errors", run)
    def test_returned_tdr_gate_accepts_exact_contract(self) -> None:
        topology = topology_module.load_topology(TOPOLOGY_FIXTURE, "sketch")
        with tempfile.TemporaryDirectory() as tmp:
            export = Path(tmp)
            write_neutral_export(export, topology)
            report = gate.validate_returned_tdr(topology, export)

        self.assertTrue(report["passed"])
        self.assertEqual(
            (report["node_count"], report["triangle_count"], report["edge_count"]),
            (6, 4, 9),
        )
        self.assertEqual(report["contact_edges"], {"Anode": [1, 5], "Cathode": [3, 4]})

    def test_returned_tdr_gate_accepts_binary_rounding_of_doping(self) -> None:
        topology = topology_module.load_topology(TOPOLOGY_FIXTURE, "sketch")
        with tempfile.TemporaryDirectory() as tmp:
            export = Path(tmp)
            write_neutral_export(export, topology)
            path = export / "doping.csv"
            path.write_text(
                path.read_text(encoding="utf-8").replace(
                    "1e+17", "1.0000000000000002e+17"
                ),
                encoding="utf-8",
            )
            report = gate.validate_returned_tdr(topology, export)

        self.assertTrue(report["doping_matches"])
        self.assertEqual(report["doping_relative_tolerance"], 1.0e-15)

    def test_returned_tdr_gate_rejects_added_node(self) -> None:
        topology = topology_module.load_topology(TOPOLOGY_FIXTURE, "sketch")
        with tempfile.TemporaryDirectory() as tmp:
            export = Path(tmp)
            write_neutral_export(export, topology)
            with (export / "nodes.csv").open("a", encoding="utf-8") as handle:
                handle.write("6,9.0,9.0\n")
            with self.assertRaisesRegex(ValueError, "expected 6 nodes"):
                gate.validate_returned_tdr(topology, export)

    def test_returned_tdr_gate_rejects_topology_contact_and_doping_changes(self) -> None:
        topology = topology_module.load_topology(TOPOLOGY_FIXTURE, "sketch")
        corruptions = (
            ("elements.csv", "0,0,4,1", "0,0,5,1", "triangle connectivity"),
            ("contacts.csv", "Anode,0;4", "Anode,0;1", "contact edges"),
            ("doping.csv", "0,0.0,1e+17", "0,1.0,1e+17", "doping"),
        )
        for filename, old, new, message in corruptions:
            with self.subTest(filename=filename):
                with tempfile.TemporaryDirectory() as tmp:
                    export = Path(tmp)
                    write_neutral_export(export, topology)
                    path = export / filename
                    original = path.read_text(encoding="utf-8")
                    changed = original.replace(old, new, 1)
                    self.assertNotEqual(changed, original)
                    path.write_text(changed, encoding="utf-8")
                    with self.assertRaisesRegex(ValueError, message):
                        gate.validate_returned_tdr(topology, export)


if __name__ == "__main__":
    unittest.main()
