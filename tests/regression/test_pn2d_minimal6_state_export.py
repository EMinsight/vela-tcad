import csv
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
MODULE_PATH = REPO / "scripts" / "export_pn2d_minimal6_states.py"
SPEC = importlib.util.spec_from_file_location("export_pn2d_minimal6_states", MODULE_PATH)
export = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(export)


SCALAR_FIELDS = {
    "ElectrostaticPotential": "V",
    "eQuasiFermiPotential": "V",
    "hQuasiFermiPotential": "V",
    "eDensity": "cm^-3",
    "hDensity": "cm^-3",
    "eMobility": "cm^2*V^-1*s^-1",
    "hMobility": "cm^2*V^-1*s^-1",
    "eAlphaAvalanche": "cm^-1",
    "hAlphaAvalanche": "cm^-1",
}
VECTOR_FIELDS = {
    "ElectricField": "V*cm^-1",
    "eCurrentDensity": "A*cm^-2",
    "hCurrentDensity": "A*cm^-2",
}


def valid_field_manifest() -> dict[str, object]:
    fields = []
    for name, unit in SCALAR_FIELDS.items():
        fields.append({
            "name": name,
            "region": 0,
            "components": 1,
            "unit": unit,
            "mapping_status": "complete",
            "global_node_mapping": "global_vertex_order",
        })
    for name, unit in VECTOR_FIELDS.items():
        fields.append({
            "name": name,
            "region": 0,
            "components": 2,
            "unit": unit,
            "mapping_status": "complete",
            "global_node_mapping": "global_vertex_order",
        })
    return {"fields": fields}


def write_csv(path: Path, header: list[str], rows: list[list[object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerows(rows)


def write_neutral_state(root: Path, missing: str | None = None) -> None:
    root.mkdir(parents=True, exist_ok=True)
    (root / "field_manifest.json").write_text(
        json.dumps(valid_field_manifest(), indent=2) + "\n", encoding="utf-8"
    )
    # Deliberately reverse importer IDs. Canonical order must come from exact coordinates.
    write_csv(
        root / "nodes.csv",
        ["id", "x_um", "y_um"],
        [[15 - i, x, y] for i, (x, y) in enumerate(
            [(0.0, 0.5), (1.0, 0.5), (2.0, 0.5),
             (2.0, 0.0), (0.0, 0.0), (1.0, 0.0)]
        )],
    )
    values = {
        "ElectrostaticPotential": [0.01 * i for i in range(6)],
        "eQuasiFermiPotential": [0.1 + 0.01 * i for i in range(6)],
        "hQuasiFermiPotential": [-0.2 - 0.01 * i for i in range(6)],
        "eDensity": [1.0e10 * (i + 1) for i in range(6)],
        "hDensity": [2.0e10 * (i + 1) for i in range(6)],
    }
    source_ids = [15 - i for i in range(6)]
    for name, field_values in values.items():
        if name == missing:
            continue
        write_csv(
            root / "fields" / f"{name}_region0.csv",
            ["node_id", "component0"],
            [[node_id, value] for node_id, value in zip(source_ids, field_values)],
        )


class PN2DMinimal6StateExportTest(unittest.TestCase):
    def test_required_state_matrix_is_exact_and_complete(self) -> None:
        states = [
            {"topology_id": topology, "requested_bias_V": bias, "actual_bias_V": bias,
             "status": "passed"}
            for topology in ("sketch", "mirror")
            for bias in (0.0, -12.0, -19.0)
        ]
        matrix = export.validate_state_matrix(states)
        self.assertEqual(
            set(matrix),
            {(t, v) for t in ("sketch", "mirror") for v in (0.0, -12.0, -19.0)},
        )

    def test_nearest_bias_is_rejected_at_stricter_than_one_picovolt(self) -> None:
        states = [
            {"topology_id": topology, "requested_bias_V": bias,
             "actual_bias_V": bias + (2.0e-12 if topology == "mirror" and bias == -19.0 else 0.0),
             "status": "passed"}
            for topology in ("sketch", "mirror")
            for bias in (0.0, -12.0, -19.0)
        ]
        with self.assertRaisesRegex(ValueError, r"mirror.*-19.*1e-12 V"):
            export.validate_state_matrix(states)

    def test_final_contact_voltage_at_one_picovolt_is_accepted(self) -> None:
        self.assertEqual(export.validate_final_bias(-19.0, -19.0 + 1.0e-12), -19.0 + 1.0e-12)

    def test_vector_fields_require_two_components_and_cm_units(self) -> None:
        manifest = valid_field_manifest()
        current = next(f for f in manifest["fields"] if f["name"] == "eCurrentDensity")
        current["unit"] = "A*m^-2"
        with self.assertRaisesRegex(ValueError, r"eCurrentDensity.*A\*cm\^-2"):
            export.validate_field_manifest(manifest)
        current["unit"] = "A*cm^-2"
        current["components"] = 1
        with self.assertRaisesRegex(ValueError, r"eCurrentDensity.*components=2"):
            export.validate_field_manifest(manifest)

    def test_every_field_requires_complete_global_region_zero_mapping(self) -> None:
        for key, bad in (("region", 1), ("mapping_status", "partial"),
                         ("global_node_mapping", "region_node_order")):
            manifest = valid_field_manifest()
            field = next(f for f in manifest["fields"] if f["name"] == "ElectricField")
            field[key] = bad
            with self.subTest(key=key):
                with self.assertRaisesRegex(ValueError, rf"ElectricField.*{key}"):
                    export.validate_field_manifest(manifest)

    def test_scalar_field_units_are_explicit(self) -> None:
        for name, expected in SCALAR_FIELDS.items():
            manifest = valid_field_manifest()
            field = next(f for f in manifest["fields"] if f["name"] == name)
            field["unit"] = "wrong"
            with self.subTest(name=name):
                with self.assertRaisesRegex(ValueError, f"{name}.*{expected.replace('*', r'\*').replace('^', r'\^')}"):
                    export.validate_field_manifest(manifest)

    def test_missing_quasi_fermi_field_is_fatal(self) -> None:
        for field, role in (("eQuasiFermiPotential", "phin"),
                            ("hQuasiFermiPotential", "phip")):
            with self.subTest(field=field), tempfile.TemporaryDirectory() as tmp:
                root = Path(tmp)
                write_neutral_state(root, missing=field)
                with self.assertRaisesRegex(ValueError, rf"missing.*{field}.*{role}"):
                    export.write_state_csv(root, export.canonical_minimal6_coordinates())

    def test_state_csv_uses_canonical_order_and_si_values(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            write_neutral_state(root)
            path = export.write_state_csv(root, export.canonical_minimal6_coordinates())
            with path.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
        self.assertEqual(
            list(rows[0]),
            ["node_id", "psi_V", "phin_V", "phip_V", "n_m3", "p_m3"],
        )
        self.assertEqual([int(row["node_id"]) for row in rows], list(range(6)))
        self.assertAlmostEqual(float(rows[3]["psi_V"]), 0.03)
        self.assertAlmostEqual(float(rows[3]["phin_V"]), 0.13)
        self.assertAlmostEqual(float(rows[3]["phip_V"]), -0.23)
        self.assertAlmostEqual(float(rows[3]["n_m3"]), 4.0e16)
        self.assertAlmostEqual(float(rows[3]["p_m3"]), 8.0e16)

    def test_deck_generation_is_exact_bias_and_explicit_grid_only(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            manifest = export.prepare_exports(
                topology_ids=("sketch", "mirror"),
                biases=(0.0, -12.0, -19.0),
                run_id="minimal6_states_test",
                output_dir=Path(tmp),
                ssh_target="sentaurus",
            )
            self.assertEqual(manifest["schema"], "vela.pn2d_minimal6_states.v1")
            self.assertEqual(len(manifest["states"]), 6)
            for state in manifest["states"]:
                command_text = " ".join(state["remote_commands"])
                self.assertIn(
                    "tdx -d pn2d_minimal6.grd pn2d_minimal6.dat pn2d_minimal6.tdr",
                    command_text,
                )
                self.assertNotRegex(command_text.lower(), r"\b(?:sde|snmesh)\b")
                self.assertNotIn("nearest", command_text.lower())
                deck = Path(state["bundle_dir"], state["deck_name"]).read_text(encoding="utf-8")
                self.assertNotIn("__TARGET_BIAS_V__", deck)
                self.assertIn(f'Voltage={state["requested_bias_V"]:.17g}', deck)
                self.assertIn('Grid      = "pn2d_minimal6.tdr"', deck)
                self.assertIn('Doping    = "pn2d_minimal6.tdr"', deck)
                self.assertIn("ElectricField/Vector", deck)
                self.assertIn("eCurrentDensity/Vector", deck)
                self.assertIn("hCurrentDensity/Vector", deck)

    def test_missing_neutral_outputs_are_fail_closed(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            manifest = export.prepare_exports(
                topology_ids=("sketch", "mirror"),
                biases=(0.0, -12.0, -19.0),
                run_id="minimal6_states_missing_output",
                output_dir=Path(tmp),
                ssh_target="sentaurus",
            )
            with self.assertRaisesRegex(ValueError, "missing neutral export directory"):
                export.run_exports(manifest, executor=lambda state: None)
            saved = json.loads(Path(manifest["manifest_path"]).read_text(encoding="utf-8"))
        self.assertFalse(saved["outputs_complete"])
        self.assertEqual(saved["states"][0]["status"], "failed")
    def test_partial_manifest_is_fail_closed_and_stops_later_states(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            manifest = export.prepare_exports(
                topology_ids=("sketch", "mirror"),
                biases=(0.0, -12.0, -19.0),
                run_id="minimal6_states_partial",
                output_dir=root,
                ssh_target="sentaurus",
            )
            calls: list[tuple[str, float]] = []

            def executor(state: dict[str, object]) -> dict[str, object]:
                key = (str(state["topology_id"]), float(state["requested_bias_V"]))
                calls.append(key)
                if key == ("sketch", -12.0):
                    raise RuntimeError("mock convergence failure")
                neutral = Path(str(state["export_dir"]))
                write_neutral_state(neutral)
                return {"actual_bias_V": key[1], "export_dir": str(neutral)}

            with self.assertRaisesRegex(RuntimeError, "mock convergence failure"):
                export.run_exports(manifest, executor=executor)
            saved = json.loads(Path(manifest["manifest_path"]).read_text(encoding="utf-8"))
        self.assertFalse(saved["outputs_complete"])
        self.assertEqual(calls, [("sketch", 0.0), ("sketch", -12.0)])
        self.assertEqual(saved["states"][1]["status"], "failed")
        self.assertTrue(all(state["status"] == "prepared" for state in saved["states"][2:]))

    def test_mocked_six_state_export_writes_complete_si_states(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            manifest = export.prepare_exports(
                topology_ids=("sketch", "mirror"),
                biases=(0.0, -12.0, -19.0),
                run_id="minimal6_states_mocked",
                output_dir=Path(tmp),
                ssh_target="sentaurus",
            )

            def executor(state: dict[str, object]) -> dict[str, object]:
                neutral = Path(str(state["export_dir"]))
                write_neutral_state(neutral)
                return {
                    "actual_bias_V": float(state["requested_bias_V"]),
                    "export_dir": str(neutral),
                }

            export.run_exports(manifest, executor=executor)
            saved = json.loads(Path(manifest["manifest_path"]).read_text(encoding="utf-8"))
            matrix = export.validate_state_matrix(saved["states"])
            state_paths = [Path(state["state_csv"]) for state in saved["states"]]
        self.assertTrue(saved["outputs_complete"])
        self.assertEqual(len(matrix), 6)
        self.assertTrue(all(path.name == "state.csv" for path in state_paths))

if __name__ == "__main__":
    unittest.main()
