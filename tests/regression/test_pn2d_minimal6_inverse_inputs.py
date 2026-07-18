import csv
import hashlib
import importlib.util
import json
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
MODULE_PATH = REPO / "scripts" / "pn2d_minimal6_diagnostics" / "inverse_inputs.py"
SPEC = importlib.util.spec_from_file_location("inverse_inputs", MODULE_PATH)
inverse_inputs = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = inverse_inputs
SPEC.loader.exec_module(inverse_inputs)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_csv(path: Path, header: list[str], row: list[object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerow(row)


REQUIRED = {
    "ElectrostaticPotential": (1, "V"), "eQuasiFermiPotential": (1, "V"),
    "hQuasiFermiPotential": (1, "V"), "eDensity": (1, "cm^-3"),
    "hDensity": (1, "cm^-3"), "ElectricField": (2, "V*cm^-1"),
    "eCurrentDensity": (2, "A*cm^-2"), "hCurrentDensity": (2, "A*cm^-2"),
    "eAlphaAvalanche": (1, "cm^-1"), "hAlphaAvalanche": (1, "cm^-1"),
    "ImpactIonization": (1, "cm^-3*s^-1"),
}
SUPPLEMENTAL = {
    "eMobility": (1, "cm^2*V^-1*s^-1"), "hMobility": (1, "cm^2*V^-1*s^-1"),
    "eVelocity": (1, "cm*s^-1"), "hVelocity": (1, "cm*s^-1"),
}


def write_root(root: Path, solver: str, fields: dict[str, tuple[int, str]],
               *, omit_field: str | None = None) -> None:
    states = []
    member_sha256 = {}
    for topology in ("sketch", "mirror"):
        for bias in range(-1, -21, -1):
            token = f"m{abs(bias)}p000000"
            relative = f"{solver}/{topology}/states/segment_00_bias_{token}.csv"
            values: list[object] = [7, 1.0, 2.0]
            header = ["canonical_node_id", "x_um", "y_um"]
            for name, (components, _) in fields.items():
                for component in range(components):
                    header.append(f"{name}_component{component}")
                    values.append("" if name == omit_field else float(abs(bias) + component + 1))
            state_path = root / relative
            write_csv(state_path, header, values)
            digest = sha256(state_path)
            member_sha256[relative] = digest
            states.append({
                "topology": topology, "requested_bias_V": float(bias),
                "actual_bias_V": float(bias), "state_path": relative,
                "state_sha256": digest, "support_kind": "node",
                "coordinate_frame": "minimal6_cartesian", "orientation": "+x,+y",
            })
    manifest = {
        "schema": "vela.pn2d_minimal6_inverse_input.v1", "solver": solver,
        "bias_tolerance_V": 1.0e-12,
        "fields": [
            {"name": name, "components": components, "unit": unit,
             "support_kind": "node"}
            for name, (components, unit) in fields.items()
        ],
        "states": states, "member_sha256": member_sha256,
        "provenance": {
            "executable_sha256": "a" * 64,
            "tracked_source_sha256": {"tracked/source.cpp": "b" * 64},
        },
    }
    manifest_path = root / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    (root / "seal.json").write_text(json.dumps({"manifest_sha256": sha256(manifest_path)}) + "\n", encoding="utf-8")


def rewrite_seal(root: Path) -> None:
    path = root / "manifest.json"
    (root / "seal.json").write_text(json.dumps({"manifest_sha256": sha256(path)}) + "\n", encoding="utf-8")


def make_fixture(root: Path, *, omit_field: str | None = None) -> tuple[Path, Path, Path]:
    vela, sentaurus, supplemental = root / "vela-root", root / "sentaurus-root", root / "supplemental-root"
    write_root(vela, "vela", REQUIRED)
    write_root(sentaurus, "sentaurus", REQUIRED)
    write_root(supplemental, "supplemental", SUPPLEMENTAL, omit_field=omit_field)
    return vela, sentaurus, supplemental


class PN2DMinimal6InverseInputsTest(unittest.TestCase):
    def test_loads_exact_matrix_hashes_fields_and_supplement(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            vela, sentaurus, supplemental = make_fixture(Path(tmp))
            bundle = inverse_inputs.load_input_bundle(vela, sentaurus, supplemental)
            self.assertEqual(len(bundle.common_keys), 40)
            self.assertEqual(bundle.discovery_keys, tuple(("sketch", float(v)) for v in (-1, -4, -8, -12, -16, -19, -20)))
            self.assertEqual(len(bundle.holdout_keys), 33)
            self.assertEqual(inverse_inputs.field_inventory(bundle)["eMobility"]["unit"], "cm^2*V^-1*s^-1")
            self.assertEqual(bundle.executable_hashes, ("a" * 64,))
            observations = inverse_inputs.canonical_observations(bundle)
            electric = next(item for item in observations if item.quantity == "ElectricField" and item.component == "component0")
            mobility = next(item for item in observations if item.quantity == "eMobility")
            self.assertEqual(electric.value_si, electric.raw_value * 100.0)
            self.assertEqual(mobility.value_si, mobility.raw_value * 1.0e-4)
            coordinate = next(item for item in observations if item.solver == "vela" and item.quantity == "coordinate" and item.component == "x")
            self.assertEqual((coordinate.raw_value, coordinate.raw_unit), (1.0, "um"))
            self.assertEqual((coordinate.value_si, coordinate.unit_si), (1.0e-6, "m"))
            self.assertEqual(coordinate.conversion, "um_to_m")
            output = Path(tmp) / "out" / "input_manifest.json"
            inverse_inputs.write_input_manifest(bundle, output)
            self.assertTrue(output.is_file())
            self.assertEqual(json.loads(output.read_text(encoding="utf-8"))["executable_sha256"], ["a" * 64])

    def test_rejects_tampered_state_before_output(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            vela, sentaurus, supplemental = make_fixture(Path(tmp))
            target = vela / "vela/sketch/states/segment_00_bias_m1p000000.csv"
            target.write_text("tampered", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "state hash mismatch"):
                inverse_inputs.load_input_bundle(vela, sentaurus, supplemental)
            self.assertFalse((Path(tmp) / "out").exists())

    def test_rejects_duplicate_and_inexact_biases(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            vela, sentaurus, supplemental = make_fixture(Path(tmp))
            path = vela / "manifest.json"
            manifest = json.loads(path.read_text(encoding="utf-8"))
            manifest["states"].append(dict(manifest["states"][0]))
            path.write_text(json.dumps(manifest), encoding="utf-8")
            rewrite_seal(vela)
            with self.assertRaisesRegex(ValueError, "duplicate"):
                inverse_inputs.load_input_bundle(vela, sentaurus, supplemental)
            make_fixture(Path(tmp))
            path = vela / "manifest.json"
            manifest = json.loads(path.read_text(encoding="utf-8"))
            manifest["states"][0]["actual_bias_V"] = -1.0001
            path.write_text(json.dumps(manifest), encoding="utf-8")
            rewrite_seal(vela)
            with self.assertRaisesRegex(ValueError, "exact bias"):
                inverse_inputs.load_input_bundle(vela, sentaurus, supplemental)

    def test_rejects_escape_supplement_mismatch_and_field_contract_mismatch(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            vela, sentaurus, supplemental = make_fixture(Path(tmp))
            path = vela / "manifest.json"
            manifest = json.loads(path.read_text(encoding="utf-8"))
            manifest["states"][0]["state_path"] = "../outside.csv"
            path.write_text(json.dumps(manifest), encoding="utf-8")
            rewrite_seal(vela)
            with self.assertRaisesRegex(ValueError, "escapes root"):
                inverse_inputs.load_input_bundle(vela, sentaurus, supplemental)
            make_fixture(Path(tmp))
            path = supplemental / "manifest.json"
            manifest = json.loads(path.read_text(encoding="utf-8"))
            manifest["states"][0]["topology"] = "other"
            path.write_text(json.dumps(manifest), encoding="utf-8")
            rewrite_seal(supplemental)
            with self.assertRaisesRegex(ValueError, "supplemental.*matrix"):
                inverse_inputs.load_input_bundle(vela, sentaurus, supplemental)
            make_fixture(Path(tmp))
            path = sentaurus / "manifest.json"
            manifest = json.loads(path.read_text(encoding="utf-8"))
            manifest["fields"][0]["unit"] = "mV"
            path.write_text(json.dumps(manifest), encoding="utf-8")
            rewrite_seal(sentauraus := sentaurus)
            with self.assertRaisesRegex(ValueError, "field.*unit"):
                inverse_inputs.load_input_bundle(vela, sentaurus, supplemental)

    def test_missing_supplemental_field_is_typed_missing_not_zero(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            vela, sentaurus, supplemental = make_fixture(Path(tmp), omit_field="eMobility")
            bundle = inverse_inputs.load_input_bundle(vela, sentaurus, supplemental)
            observations = inverse_inputs.canonical_observations(bundle)
            mobility = [item for item in observations if item.quantity == "eMobility"]
            self.assertEqual(len(mobility), 40)
            self.assertTrue(all(item.raw_value is None and item.value_si is None for item in mobility))
            self.assertTrue(all(item.status is inverse_inputs.SampleStatus.MISSING_FIELD for item in mobility))


if __name__ == "__main__":
    unittest.main()
