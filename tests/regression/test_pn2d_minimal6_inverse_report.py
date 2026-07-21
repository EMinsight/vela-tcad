import hashlib
import importlib.util
import json
import shutil
import tempfile
import unittest
from pathlib import Path

from scripts.diagnose_pn2d_minimal6_physics_inverse_audit import build_report_package
from scripts.verify_pn2d_minimal6_physics_inverse_audit import verify_report


REPO = Path(__file__).resolve().parents[2]
INPUT_FIXTURE_SPEC = importlib.util.spec_from_file_location(
    "inverse_input_fixture",
    REPO / "tests" / "regression" / "test_pn2d_minimal6_inverse_inputs.py",
)
INPUT_FIXTURE = importlib.util.module_from_spec(INPUT_FIXTURE_SPEC)
assert INPUT_FIXTURE_SPEC.loader is not None
INPUT_FIXTURE_SPEC.loader.exec_module(INPUT_FIXTURE)

ARTIFACTS = (
    "input_manifest.json", "observations_node.csv", "observations_edge.csv",
    "observations_cell.csv", "observations_contact.csv",
    "observations_integrated.csv", "candidate_metrics.csv",
    "candidate_classifications.json", "replacement_matrix.csv",
    "physics_inverse_audit.json", "physics_inverse_audit.md",
    "figure_manifest.json", "report_manifest.json", "verification.json",
    "package_manifest.json",
)
FIGURES = (
    "potential_field", "qf_gradient", "current_density",
    "alpha_generation", "replacement_matrix",
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


class PN2DMinimal6InverseReportTest(unittest.TestCase):
    def build(self, root: Path, name: str) -> Path:
        vela, sentaurus, supplemental = INPUT_FIXTURE.make_fixture(root / "inputs")
        out = root / name
        result = build_report_package(
            vela_root=vela,
            sentaurus_root=sentaurus,
            supplemental_sentaurus_root=supplemental,
            out_dir=out,
            phase_base="a5524cf",
        )
        self.assertTrue(result["passed"])
        return out

    def test_two_runs_are_byte_identical_and_report_contract_is_complete(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            out_a = self.build(root, "out-a")
            out_b = self.build(root, "out-b")
            self.assertEqual(
                {name: sha256(out_a / name) for name in ARTIFACTS},
                {name: sha256(out_b / name) for name in ARTIFACTS},
            )
            for stem in FIGURES:
                for suffix in (".png", ".pdf"):
                    self.assertEqual(
                        sha256(out_a / "figures" / f"{stem}{suffix}"),
                        sha256(out_b / "figures" / f"{stem}{suffix}"),
                    )
            first_package = (out_a / "package_manifest.json").read_bytes()
            first_verification = (out_a / "verification.json").read_bytes()
            self.assertTrue(verify_report(out_a)["passed"])
            self.assertEqual(first_package, (out_a / "package_manifest.json").read_bytes())
            self.assertEqual(first_verification, (out_a / "verification.json").read_bytes())

            report = json.loads((out_a / "physics_inverse_audit.json").read_text())
            self.assertEqual(report["phase_base"], "a5524cf")
            self.assertEqual(len(report["payload"]["discovery_keys"]), 7)
            self.assertEqual(len(report["payload"]["holdout_keys"]), 33)
            markdown = (out_a / "physics_inverse_audit.md").read_text()
            for heading in (
                "# PN2D Minimal6 Physics Inverse Audit", "## Technical summary",
                "## Scope, data, and metric definitions", "## Methodology",
                "## Limitations, uncertainty, and robustness",
                "## Recommended next steps", "## Further questions",
            ):
                self.assertIn(heading, markdown)
            manifest = json.loads((out_a / "figure_manifest.json").read_text())
            self.assertEqual([item["name"] for item in manifest["figures"]], list(FIGURES))
            for item in manifest["figures"]:
                self.assertEqual(
                    set(item["chart_contract"]),
                    {"question", "takeaway", "family", "variant", "row_grain_sufficiency",
                     "fields", "palette_policy", "output_paths", "qa_surface"},
                )

    def test_independent_verifier_rejects_csv_json_png_and_input_mutations(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pristine = self.build(root, "pristine")
            cases = []
            for label in ("csv", "json", "png", "input"):
                target = root / label
                shutil.copytree(pristine, target)
                cases.append((label, target))

            csv_path = cases[0][1] / "candidate_metrics.csv"
            csv_path.write_bytes(csv_path.read_bytes().replace(b"identified", b"rejected", 1))

            json_path = cases[1][1] / "candidate_classifications.json"
            payload = json.loads(json_path.read_text())
            payload["classifications"][0]["classification"] = "rejected"
            json_path.write_text(json.dumps(payload, sort_keys=True, indent=2) + "\n")

            png_path = cases[2][1] / "figures" / "potential_field.png"
            png = bytearray(png_path.read_bytes())
            png[-20] ^= 1
            png_path.write_bytes(png)

            report_manifest = json.loads((cases[3][1] / "report_manifest.json").read_text())
            raw_path = Path(report_manifest["inputs"][0]["path"])
            raw_path.write_bytes(raw_path.read_bytes() + b"mutation")

            for label, target in cases:
                with self.subTest(label=label):
                    with self.assertRaises(ValueError):
                        verify_report(target)


if __name__ == "__main__":
    unittest.main()
