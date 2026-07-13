#!/usr/bin/env python3
"""Regression coverage for the PN2D minimal6 canonical topology fixtures."""

from __future__ import annotations

import importlib.util
import json
import re
import subprocess
import sys
import tempfile
import unittest
from dataclasses import replace
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
FIXTURE = REPO / "reference_tcad" / "pn2d_sentaurus2018_minimal6" / "source" / "minimal6_topologies.json"
MODULE_SPEC = importlib.util.spec_from_file_location(
    "pn2d_minimal6_topology",
    REPO / "scripts" / "pn2d_minimal6_topology.py",
)
assert MODULE_SPEC is not None
module = importlib.util.module_from_spec(MODULE_SPEC)
assert MODULE_SPEC.loader is not None
sys.modules[MODULE_SPEC.name] = module
MODULE_SPEC.loader.exec_module(module)


class PN2DMinimal6TopologyTest(unittest.TestCase):
    def test_loaded_topology_nested_collections_are_immutable(self) -> None:
        topology = module.load_topology(FIXTURE, "sketch")

        with self.assertRaises(TypeError):
            topology.nodes[1] = (0.1, 0.5)
        with self.assertRaises(AttributeError):
            topology.triangles.append((1, 2, 3))
        with self.assertRaises(TypeError):
            topology.contacts["Anode"] = (1, 2)
        with self.assertRaises(TypeError):
            topology.acceptors_cm3[1] = 0.0
        with self.assertRaises(TypeError):
            topology.donors_cm3[3] = 0.0


    def test_sketch_and_mirror_have_exact_contract(self) -> None:
        sketch = module.load_topology(FIXTURE, "sketch")
        mirror = module.load_topology(FIXTURE, "mirror")
        self.assertEqual(sketch.triangles, ((1, 5, 2), (5, 6, 2), (2, 6, 4), (2, 4, 3)))
        self.assertEqual(mirror.triangles, ((1, 5, 6), (1, 6, 2), (2, 6, 3), (6, 4, 3)))
        for topology in (sketch, mirror):
            summary = module.validate_topology(topology)
            self.assertEqual((summary.nodes, summary.triangles, summary.edges), (6, 4, 9))
            self.assertEqual(summary.contact_edges, {"Anode": (1, 5), "Cathode": (3, 4)})


    def test_dfise_roundtrip_preserves_exact_contract(self) -> None:
        for topology_id in ("sketch", "mirror"):
            with self.subTest(topology_id=topology_id):
                topology = module.load_topology(FIXTURE, topology_id)
                with tempfile.TemporaryDirectory() as tmp:
                    grd, dat = Path(tmp) / "mesh.grd", Path(tmp) / "mesh.dat"
                    module.write_dfise_grid(topology, grd)
                    module.write_dfise_doping(topology, dat)
                    report = module.validate_dfise_roundtrip(topology, grd, dat)
                self.assertTrue(report["passed"])
                self.assertEqual(report["vertices"], topology.nodes)
                self.assertEqual(report["edge_count"], 9)
                self.assertEqual(report["location_counts"], {"e": 6, "i": 3})
                self.assertEqual(report["triangles"], topology.triangles)
                self.assertEqual(report["silicon_triangle_count"], 4)
                self.assertEqual(
                    report["contact_element_counts"],
                    {"Anode": 1, "Cathode": 1},
                )
                self.assertEqual(report["contact_edges"], topology.contacts)
                self.assertEqual(
                    report["region_ownership"],
                    {"R.Si": [0, 1, 2, 3], "Cathode": [4], "Anode": [5]},
                )


    def test_dfise_roundtrip_rejects_incomplete_grid_metadata(self) -> None:
        topology = module.load_topology(FIXTURE, "sketch")
        corruptions = (
            ("version = 1.1", "version = 9.9"),
            ("type    = grid", "type    = dataset"),
            ("dimension   = 2", "dimension   = 3"),
            ("nb_vertices = 6", "nb_vertices = 7"),
            ("nb_edges = 9", "nb_edges = 10"),
            ("nb_faces = 0", "nb_faces = 1"),
            ("nb_elements = 6", "nb_elements = 7"),
            ("nb_regions = 3", "nb_regions = 2"),
            (
                "translate = [  0.000000000000000e+00 0.000000000000000e+00 0.000000000000000e+00 ]",
                "translate = [  0.000000000000000e+00 0.000000000000000e+00 ]",
            ),
            (
                "0.000000000000000e+00 0.000000000000000e+00 1.000000000000000e+00",
                "0.000000000000000e+00 0.000000000000000e+00 2.000000000000000e+00",
            ),
        )
        for old, new in corruptions:
            with self.subTest(old=old):
                with tempfile.TemporaryDirectory() as tmp:
                    grd, dat = Path(tmp) / "mesh.grd", Path(tmp) / "mesh.dat"
                    module.write_dfise_grid(topology, grd)
                    module.write_dfise_doping(topology, dat)
                    original = grd.read_text()
                    changed = original.replace(old, new, 1)
                    self.assertNotEqual(changed, original)
                    grd.write_text(changed)
                    report = module.validate_dfise_roundtrip(topology, grd, dat)
                self.assertFalse(report["passed"])


    def test_dfise_roundtrip_rejects_incomplete_dataset_metadata(self) -> None:
        topology = module.load_topology(FIXTURE, "sketch")
        corruptions = (
            ("version = 1.0", "version = 9.9"),
            ("type    = dataset", "type    = grid"),
            ("dimension   = 2", "dimension   = 3"),
            ("nb_vertices = 6", "nb_vertices = 7"),
            ("nb_edges    = 9", "nb_edges    = 10"),
            ("nb_faces    = 0", "nb_faces    = 1"),
            ("nb_elements = 6", "nb_elements = 7"),
            ("nb_regions  = 3", "nb_regions  = 2"),
            ("function  = DopingConcentration", "function  = WrongFunction"),
            ("dimension = 1", "dimension = 2"),
            ("type      = scalar", "type      = vector"),
            ("location  = vertex", "location  = edge"),
            ('validity  = [ "R.Si" ]', 'validity  = [ "Cathode" ]'),
        )
        for old, new in corruptions:
            with self.subTest(old=old):
                with tempfile.TemporaryDirectory() as tmp:
                    grd, dat = Path(tmp) / "mesh.grd", Path(tmp) / "mesh.dat"
                    module.write_dfise_grid(topology, grd)
                    module.write_dfise_doping(topology, dat)
                    original = dat.read_text()
                    changed = original.replace(old, new, 1)
                    self.assertNotEqual(changed, original)
                    dat.write_text(changed)
                    report = module.validate_dfise_roundtrip(topology, grd, dat)
                self.assertFalse(report["passed"])


    def test_dfise_roundtrip_rejects_unconsumed_tokens_and_wrong_material(self) -> None:
        topology = module.load_topology(FIXTURE, "sketch")
        corruptions = (
            ("grd", "material = Silicon", "material = Oxide"),
            ("grd", "eeeiiieee", "xeeeiiieee"),
            ("dat", 'validity  = [ "R.Si" ]', 'validity  = [ "R.Si" !!! ]'),
        )
        for target, old, new in corruptions:
            with self.subTest(target=target, old=old):
                with tempfile.TemporaryDirectory() as tmp:
                    grd, dat = Path(tmp) / "mesh.grd", Path(tmp) / "mesh.dat"
                    module.write_dfise_grid(topology, grd)
                    module.write_dfise_doping(topology, dat)
                    path = grd if target == "grd" else dat
                    original = path.read_text()
                    changed = original.replace(old, new, 1)
                    self.assertNotEqual(changed, original)
                    path.write_text(changed)
                    report = module.validate_dfise_roundtrip(topology, grd, dat)
                self.assertFalse(report["passed"])


    def test_dfise_roundtrip_rejects_duplicate_unknown_and_trailing_metadata(self) -> None:
        topology = module.load_topology(FIXTURE, "sketch")
        with tempfile.TemporaryDirectory() as tmp:
            grd, dat = Path(tmp) / "mesh.grd", Path(tmp) / "mesh.dat"
            module.write_dfise_grid(topology, grd)
            module.write_dfise_doping(topology, dat)
            grid_text = grd.read_text()
            dataset_text = dat.read_text()

        info = re.search(r"Info \{.*?\n\}", grid_text, re.S)
        region = re.search(
            r'  Region \("R.Si"\) \{.*?(?=  Region \("Cathode"\))',
            grid_text,
            re.S,
        )
        dataset = re.search(
            r'  Dataset \("DopingConcentration"\) \{.*?'
            r'(?=  Dataset \("PhosphorusActiveConcentration"\))',
            dataset_text,
            re.S,
        )
        self.assertIsNotNone(info)
        self.assertIsNotNone(region)
        self.assertIsNotNone(dataset)
        assert info is not None and region is not None and dataset is not None

        corruptions = (
            ("grd", grid_text.replace("Data {", info.group(0) + "\n\nData {", 1)),
            (
                "grd",
                grid_text.replace(
                    '  Region ("Cathode")',
                    region.group(0) + '  Region ("Cathode")',
                    1,
                ),
            ),
            (
                "grd",
                grid_text.replace(
                    "nb_regions = 3",
                    "nb_regions = 3\n  mystery = 1",
                    1,
                ),
            ),
            (
                "dat",
                dataset_text.replace(
                    '  Dataset ("PhosphorusActiveConcentration")',
                    dataset.group(0)
                    + '  Dataset ("PhosphorusActiveConcentration")',
                    1,
                ),
            ),
            ("dat", dataset_text + "\ntrailing-garbage\n"),
        )
        for target, changed in corruptions:
            with self.subTest(target=target):
                with tempfile.TemporaryDirectory() as tmp:
                    grd, dat = Path(tmp) / "mesh.grd", Path(tmp) / "mesh.dat"
                    module.write_dfise_grid(topology, grd)
                    module.write_dfise_doping(topology, dat)
                    path = grd if target == "grd" else dat
                    path.write_text(changed)
                    report = module.validate_dfise_roundtrip(topology, grd, dat)
                self.assertFalse(report["passed"])


    def test_dfise_roundtrip_rejects_junk_in_numeric_blocks(self) -> None:
        topology = module.load_topology(FIXTURE, "sketch")
        corruptions = (
            ("grd", " 0 1\n", " 0 1 junk\n"),
            ("grd", " 0 1 2 3\n", " 0 1 2 3 junk\n"),
            (
                "dat",
                " -1.000000000000000e+17\n",
                " -1.000000000000000e+17 junk\n",
            ),
        )
        for target, old, new in corruptions:
            with self.subTest(target=target, old=old):
                with tempfile.TemporaryDirectory() as tmp:
                    grd, dat = Path(tmp) / "mesh.grd", Path(tmp) / "mesh.dat"
                    module.write_dfise_grid(topology, grd)
                    module.write_dfise_doping(topology, dat)
                    path = grd if target == "grd" else dat
                    original = path.read_text()
                    changed = original.replace(old, new, 1)
                    self.assertNotEqual(changed, original)
                    path.write_text(changed)
                    report = module.validate_dfise_roundtrip(topology, grd, dat)
                self.assertFalse(report["passed"])


    def test_dfise_roundtrip_rejects_invalid_signed_edge_reference(self) -> None:
        topology = module.load_topology(FIXTURE, "sketch")
        with tempfile.TemporaryDirectory() as tmp:
            grd, dat = Path(tmp) / "mesh.grd", Path(tmp) / "mesh.dat"
            module.write_dfise_grid(topology, grd)
            module.write_dfise_doping(topology, dat)
            original = grd.read_text()
            changed = original.replace(" 2 1 -5 -1", " 2 99 -5 -1", 1)
            self.assertNotEqual(changed, original)
            grd.write_text(changed)
            with self.assertRaisesRegex(ValueError, "edge reference"):
                module.validate_dfise_roundtrip(topology, grd, dat)


    def test_dfise_roundtrip_rejects_disconnected_triangle_edge_loop(self) -> None:
        topology = module.load_topology(FIXTURE, "sketch")
        with tempfile.TemporaryDirectory() as tmp:
            grd, dat = Path(tmp) / "mesh.grd", Path(tmp) / "mesh.dat"
            module.write_dfise_grid(topology, grd)
            module.write_dfise_doping(topology, dat)
            original = grd.read_text()
            changed = original.replace(" 2 1 -5 -1", " 2 1 -5 0", 1)
            self.assertNotEqual(changed, original)
            grd.write_text(changed)
            with self.assertRaisesRegex(ValueError, "edge loop is disconnected"):
                module.validate_dfise_roundtrip(topology, grd, dat)


    def test_dfise_roundtrip_rejects_wrong_locations_and_region_ownership(self) -> None:
        topology = module.load_topology(FIXTURE, "sketch")
        corruptions = (
            ("eeeiiieee", "ieeiiieee"),
            (" 0 1 2 3", " 0 1 2 2"),
        )
        for old, new in corruptions:
            with self.subTest(old=old):
                with tempfile.TemporaryDirectory() as tmp:
                    grd, dat = Path(tmp) / "mesh.grd", Path(tmp) / "mesh.dat"
                    module.write_dfise_grid(topology, grd)
                    module.write_dfise_doping(topology, dat)
                    original = grd.read_text()
                    changed = original.replace(old, new, 1)
                    self.assertNotEqual(changed, original)
                    grd.write_text(changed)
                    report = module.validate_dfise_roundtrip(topology, grd, dat)
                self.assertFalse(report["passed"])


    def test_dfise_roundtrip_rejects_dataset_names_counts_and_values(self) -> None:
        topology = module.load_topology(FIXTURE, "sketch")
        corruptions = (
            ('Dataset ("DopingConcentration")', 'Dataset ("RenamedDoping")'),
            ("Values (6)", "Values (5)"),
            (" -1.000000000000000e+17", " -9.000000000000000e+16"),
        )
        for old, new in corruptions:
            with self.subTest(old=old):
                with tempfile.TemporaryDirectory() as tmp:
                    grd, dat = Path(tmp) / "mesh.grd", Path(tmp) / "mesh.dat"
                    module.write_dfise_grid(topology, grd)
                    module.write_dfise_doping(topology, dat)
                    original = dat.read_text()
                    changed = original.replace(old, new, 1)
                    self.assertNotEqual(changed, original)
                    dat.write_text(changed)
                    report = module.validate_dfise_roundtrip(topology, grd, dat)
                self.assertFalse(report["passed"])


    def test_dfise_roundtrip_rejects_changed_contact_endpoint(self) -> None:
        topology = module.load_topology(FIXTURE, "sketch")
        with tempfile.TemporaryDirectory() as tmp:
            grd, dat = Path(tmp) / "mesh.grd", Path(tmp) / "mesh.dat"
            module.write_dfise_grid(topology, grd)
            module.write_dfise_doping(topology, dat)
            original = grd.read_text()
            changed = original.replace(" 1 2 3\n 1 0 4", " 1 1 3\n 1 0 4")
            self.assertNotEqual(changed, original)
            grd.write_text(changed)
            report = module.validate_dfise_roundtrip(topology, grd, dat)
        self.assertFalse(report["passed"])
        self.assertEqual(report["contact_edges"]["Cathode"], (2, 4))


    def test_dfise_writers_follow_the_local_accepted_forms(self) -> None:
        topology = module.load_topology(FIXTURE, "sketch")
        with tempfile.TemporaryDirectory() as tmp:
            grd, dat = Path(tmp) / "mesh.grd", Path(tmp) / "mesh.dat"
            module.write_dfise_grid(topology, grd)
            module.write_dfise_doping(topology, dat)
            grid_text = grd.read_text()
            doping_text = dat.read_text()

        info, data = grid_text.split("Data {", maxsplit=1)
        self.assertIn("version = 1.1", info)
        for declaration in (
            "nb_vertices = 6",
            "nb_edges = 9",
            "nb_faces = 0",
            "nb_elements = 6",
            "nb_regions = 3",
            'regions = [ "R.Si" "Cathode" "Anode" ]',
            "materials = [ Silicon Contact Contact ]",
        ):
            self.assertIn(declaration, info)
            self.assertNotIn(declaration, data)
        self.assertRegex(
            data,
            r"CoordSystem\s*\{\s*translate\s*=\s*\[\s*0(?:\.0+)?e[+-]00\s+0(?:\.0+)?e[+-]00\s+0(?:\.0+)?e[+-]00\s*\]"
            r"\s*transform\s*=\s*\[\s*1(?:\.0+)?e[+-]00\s+0(?:\.0+)?e[+-]00\s+0(?:\.0+)?e[+-]00",
        )

        dat_info, dat_data = doping_text.split("Data {", maxsplit=1)
        self.assertIn("type    = dataset", dat_info)
        self.assertIn("nb_vertices = 6", dat_info)
        self.assertIn("nb_edges    = 9", dat_info)
        self.assertIn("nb_elements = 6", dat_info)
        self.assertIn("nb_regions  = 3", dat_info)
        self.assertIn(
            'datasets    = [ "DopingConcentration" "PhosphorusActiveConcentration" "BoronActiveConcentration" ]',
            dat_info,
        )
        for name in ("DopingConcentration", "PhosphorusActiveConcentration", "BoronActiveConcentration"):
            block = re.search(rf'Dataset \("{name}"\) \{{(?P<body>.*?)Values \(6\)', dat_data, re.S)
            self.assertIsNotNone(block)
            assert block is not None
            self.assertIn(f"function  = {name}", block.group("body"))
            self.assertIn("location  = vertex", block.group("body"))
            self.assertIn('validity  = [ "R.Si" ]', block.group("body"))


    def test_load_rejects_alias_duplicate_and_extra_json_keys(self) -> None:
        original = FIXTURE.read_text(encoding="utf-8")
        corruptions = (
            (
                '"1": [0.0, 0.5],',
                '"01": [0.0, 0.5], "1": [0.0, 0.5],',
                "nodes keys",
            ),
            (
                '"1": [0.0, 0.5],',
                '"1": [0.0, 0.5], "1": [0.0, 0.5],',
                "duplicate JSON key",
            ),
            (
                '"acceptors_cm3": {"1": 1e17,',
                '"acceptors_cm3": {"01": 1e17, "1": 1e17,',
                "acceptors_cm3 keys",
            ),
            (
                '"donors_cm3": {"1": 0.0,',
                '"donors_cm3": {"extra": 0.0, "1": 0.0,',
                "donors_cm3 keys",
            ),
        )
        for old, new, message in corruptions:
            with self.subTest(message=message):
                changed = original.replace(old, new, 1)
                self.assertNotEqual(changed, original)
                with tempfile.TemporaryDirectory() as tmp:
                    fixture = Path(tmp) / "topologies.json"
                    fixture.write_text(changed, encoding="utf-8")
                    with self.assertRaisesRegex(ValueError, message):
                        module.load_topology(fixture, "sketch")


    def test_validation_rejects_noncanonical_fixed_state(self) -> None:
        topology = module.load_topology(FIXTURE, "sketch")
        invalid_topologies = (
            replace(topology, triangles=[(1, 2, 5), *topology.triangles[1:]]),
            replace(topology, triangles=topology.triangles[:-1]),
            replace(topology, nodes={**topology.nodes, 1: (0.1, 0.5)}),
            replace(topology, contacts={**topology.contacts, "Anode": (1, 2)}),
            replace(topology, donors_cm3={**topology.donors_cm3, 2: 0.0}),
            replace(topology, acceptors_cm3={**topology.acceptors_cm3, 6: 0.0}),
        )
        for invalid in invalid_topologies:
            with self.subTest(topology=invalid):
                with self.assertRaises(ValueError):
                    module.validate_topology(invalid)


    def test_mirror_is_the_required_labelled_vertical_reflection(self) -> None:
        sketch = module.load_topology(FIXTURE, "sketch")
        mirror = module.load_topology(FIXTURE, "mirror")
        reflection = {1: 5, 2: 6, 3: 4, 4: 3, 5: 1, 6: 2}
        reflected_nodes = {
            reflection[node_id]: (x, 0.5 - y)
            for node_id, (x, y) in sketch.nodes.items()
        }
        reflected_triangles = {
            module.canonical_triangle(
                tuple(reversed(tuple(reflection[node_id] for node_id in triangle)))
            )
            for triangle in sketch.triangles
        }
        self.assertEqual(mirror.nodes, reflected_nodes)
        self.assertEqual(
            {module.canonical_triangle(triangle) for triangle in mirror.triangles},
            reflected_triangles,
        )


    def test_cli_default_validation_emits_json(self) -> None:
        completed = subprocess.run(
            [
                sys.executable,
                str(REPO / "scripts" / "pn2d_minimal6_topology.py"),
                "--topology",
                "sketch",
            ],
            cwd=REPO,
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        report = json.loads(completed.stdout)
        self.assertEqual(report["contact_edges"], {"Anode": [1, 5], "Cathode": [3, 4]})


    def test_cli_imports_are_robust_from_repo_root(self) -> None:
        completed = subprocess.run(
            [sys.executable, str(REPO / "scripts" / "pn2d_minimal6_topology.py"), "--help"],
            cwd=REPO,
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(completed.returncode, 0, completed.stderr)
        self.assertIn("topology", completed.stdout.lower())


if __name__ == "__main__":
    unittest.main()
