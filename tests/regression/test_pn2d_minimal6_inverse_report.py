import csv
import hashlib
import importlib.util
import json
import shutil
import tempfile
import unittest
from pathlib import Path

from PIL import Image

from scripts.diagnose_pn2d_minimal6_physics_inverse_audit import build_report_package
from scripts.pn2d_minimal6_diagnostics.inverse_inputs import (
    canonical_observations, load_input_bundle,
)
from scripts.pn2d_minimal6_diagnostics.inverse_report import _qf_series
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
    "observations_cell.csv", "observations_contact.csv", "observations_integrated.csv",
    "candidate_metrics.csv", "candidate_classifications.json", "replacement_matrix.csv",
    "physics_inverse_audit.json", "physics_inverse_audit.md", "figure_manifest.json",
    "report_manifest.json", "verification.json", "package_manifest.json",
)
FIGURES = (
    "potential_field", "qf_gradient", "current_density",
    "alpha_generation", "replacement_matrix",
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path: Path, payload: object) -> None:
    path.write_text(json.dumps(payload, sort_keys=True, indent=2, allow_nan=False) + "\n",
                    encoding="utf-8", newline="\n")


def make_three_node_fixture(root: Path, *, omit_field: str | None = None):
    roots = INPUT_FIXTURE.make_fixture(root, omit_field=omit_field)
    points = (("0", 0.0, 0.0), ("1", 2.0, 0.0), ("2", 0.0, 1.0))
    for input_root in roots:
        manifest = json.loads((input_root / "manifest.json").read_text())
        for state_index, state in enumerate(manifest["states"]):
            path = input_root / state["state_path"]
            with path.open(newline="", encoding="utf-8") as handle:
                original = list(csv.reader(handle))
            header, prototype = original[0], original[1]
            generated = []
            for node, x_um, y_um in points:
                row = list(prototype)
                row[header.index("canonical_node_id")] = node
                row[header.index("x_um")] = format(x_um, ".17g")
                row[header.index("y_um")] = format(y_um, ".17g")
                if "ElectrostaticPotential_component0" in header:
                    row[header.index("ElectrostaticPotential_component0")] = format(
                        3.0 * x_um - 4.0 * y_um + abs(float(state["actual_bias_V"])), ".17g")
                    row[header.index("ElectricField_component0")] = "-30000"
                    row[header.index("ElectricField_component1")] = "40000"
                generated.append(row)
            with path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle, lineterminator="\n")
                writer.writerow(header)
                writer.writerows(generated)
            INPUT_FIXTURE.rebind_state(input_root, state_index)
    return roots


def blank_root_field(root: Path, field: str) -> None:
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    for state_index, state in enumerate(manifest["states"]):
        path = root / state["state_path"]
        with path.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.reader(handle))
        columns = [index for index, name in enumerate(rows[0])
                   if name.startswith(f"{field}_component")]
        for row in rows[1:]:
            for index in columns:
                row[index] = ""
        with path.open("w", newline="", encoding="utf-8") as handle:
            csv.writer(handle, lineterminator="\n").writerows(rows)
        INPUT_FIXTURE.rebind_state(root, state_index)


def scale_root_fields_at_state(root: Path, *, topology: str, bias_V: float,
                               fields: tuple[str, ...], factor: float) -> None:
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    for state_index, state in enumerate(manifest["states"]):
        if state["topology"] != topology or state["actual_bias_V"] != bias_V:
            continue
        path = root / state["state_path"]
        with path.open(newline="", encoding="utf-8") as handle:
            rows = list(csv.reader(handle))
        columns = [index for index, name in enumerate(rows[0])
                   if any(name.startswith(f"{field}_component") for field in fields)]
        for row in rows[1:]:
            for index in columns:
                row[index] = format(float(row[index]) * factor, ".17g")
        with path.open("w", newline="", encoding="utf-8") as handle:
            csv.writer(handle, lineterminator="\n").writerows(rows)
        INPUT_FIXTURE.rebind_state(root, state_index)
        return
    raise AssertionError(f"state not found: {topology}@{bias_V:g}V")


def refresh_integrity(root: Path) -> None:
    report_manifest_path = root / "report_manifest.json"
    report_manifest = json.loads(report_manifest_path.read_text())
    report_manifest["artifacts"] = {
        relative: sha256(root / relative)
        for relative in sorted(report_manifest["artifacts"])
    }
    write_json(report_manifest_path, report_manifest)
    verification_path = root / "verification.json"
    verification = json.loads(verification_path.read_text())
    verification["report_manifest_sha256"] = sha256(report_manifest_path)
    write_json(verification_path, verification)
    package_path = root / "package_manifest.json"
    package = json.loads(package_path.read_text())
    package["artifacts"] = {
        path.relative_to(root).as_posix(): sha256(path)
        for path in sorted(item for item in root.rglob("*") if item.is_file())
        if path.name != "package_manifest.json"
    }
    write_json(package_path, package)


class PN2DMinimal6InverseReportTest(unittest.TestCase):
    def build(self, root: Path, name: str, *, omit_field: str | None = None) -> tuple[Path, tuple[Path, Path, Path]]:
        roots = make_three_node_fixture(root / f"inputs-{name}", omit_field=omit_field)
        out = root / name
        result = build_report_package(
            vela_root=roots[0], sentaurus_root=roots[1],
            supplemental_sentaurus_root=roots[2], out_dir=out, phase_base="a5524cf",
        )
        self.assertTrue(result["passed"])
        return out, roots

    def test_two_runs_are_byte_identical_and_report_contract_is_complete(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            roots = make_three_node_fixture(root / "shared-inputs")
            outputs = []
            for name in ("out-a", "out-b"):
                out = root / name
                build_report_package(vela_root=roots[0], sentaurus_root=roots[1],
                                     supplemental_sentaurus_root=roots[2], out_dir=out,
                                     phase_base="a5524cf")
                outputs.append(out)
            out_a, out_b = outputs
            self.assertEqual({name: sha256(out_a / name) for name in ARTIFACTS},
                             {name: sha256(out_b / name) for name in ARTIFACTS})
            for stem in FIGURES:
                for suffix in (".png", ".pdf"):
                    self.assertEqual(sha256(out_a / "figures" / f"{stem}{suffix}"),
                                     sha256(out_b / "figures" / f"{stem}{suffix}"))
            package_bytes = (out_a / "package_manifest.json").read_bytes()
            verification_bytes = (out_a / "verification.json").read_bytes()
            self.assertTrue(verify_report(out_a)["passed"])
            self.assertEqual(package_bytes, (out_a / "package_manifest.json").read_bytes())
            self.assertEqual(verification_bytes, (out_a / "verification.json").read_bytes())

    def test_authoritative_input_ledger_package_raw_inputs_and_renderer_are_truthful(self):
        with tempfile.TemporaryDirectory() as tmp:
            out, roots = self.build(Path(tmp), "report")
            report = json.loads((out / "physics_inverse_audit.json").read_text())
            provenance = report["payload"]["localization_control"]["input_provenance"]
            self.assertEqual(provenance["input_roots"], {
                "vela_root": str(roots[0].resolve()), "sentaurus_root": str(roots[1].resolve()),
                "supplemental_sentaurus_root": str(roots[2].resolve()),
            })
            self.assertEqual(len(provenance["raw_inputs"]), 126)
            self.assertTrue({"vela:manifest.json", "vela:seal.json",
                             "sentaurus:manifest.json", "sentaurus:seal.json",
                             "supplemental:manifest.json", "supplemental:seal.json"}
                            <= {item["logical_id"] for item in provenance["raw_inputs"]})
            self.assertTrue(all(Path(item["path"]).is_file() for item in provenance["raw_inputs"]))
            package = json.loads((out / "package_manifest.json").read_text())
            self.assertEqual(package["raw_inputs"], provenance["raw_inputs"])
            self.assertEqual(package["exclusions"], ["package_manifest.json"])
            renderer = json.loads((out / "figure_manifest.json").read_text())["renderer"]
            self.assertEqual(renderer["backend"], "Pillow PNG + ReportLab invariant PDF")
            self.assertEqual(renderer["determinism_scope"], "same runtime and library versions")
            self.assertNotIn("Agg", json.dumps(renderer))

    def test_missing_qf_inputs_remain_gaps_not_zero_points(self):
        with tempfile.TemporaryDirectory() as tmp:
            roots = make_three_node_fixture(Path(tmp) / "inputs", omit_field="eMobility")
            bundle = load_input_bundle(*roots)
            _, series = _qf_series(canonical_observations(bundle))
            electron = next(item for item in series if item["label"].startswith("Electron"))
            self.assertTrue(electron["values"])
            self.assertTrue(all(value is None for value in electron["values"]))
            self.assertNotIn(0.0, electron["values"])

    def test_composite_candidate_requires_every_declared_quantity_in_every_split(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            complete_roots = make_three_node_fixture(root / "complete-inputs")
            complete = root / "complete"
            build_report_package(
                vela_root=complete_roots[0], sentaurus_root=complete_roots[1],
                supplemental_sentaurus_root=complete_roots[2], out_dir=complete,
                phase_base="a5524cf",
            )
            complete_report = json.loads(
                (complete / "physics_inverse_audit.json").read_text(encoding="utf-8")
            )["payload"]
            complete_status = next(
                row for row in complete_report["classifications"]
                if row["candidate"] == "potential_field_direct"
            )
            self.assertEqual(complete_status["classification"], "identified")

            partial_roots = make_three_node_fixture(root / "partial-inputs")
            blank_root_field(partial_roots[0], "ElectricField")
            partial = root / "partial"
            build_report_package(
                vela_root=partial_roots[0], sentaurus_root=partial_roots[1],
                supplemental_sentaurus_root=partial_roots[2], out_dir=partial,
                phase_base="a5524cf",
            )
            partial_report = json.loads(
                (partial / "physics_inverse_audit.json").read_text(encoding="utf-8")
            )["payload"]
            classification = next(
                row for row in partial_report["classifications"]
                if row["candidate"] == "potential_field_direct"
            )
            self.assertEqual(classification["classification"], "insufficient_data")
            self.assertIn("declared quantities", classification["reason"])
            metrics = [
                row for row in partial_report["candidate_metrics"]
                if row["candidate"] == "potential_field_direct"
            ]
            self.assertEqual([row["split"] for row in metrics],
                             ["discovery", "holdout", "combined"])
            self.assertTrue(all(row["classification"] == "insufficient_data"
                                for row in metrics))

    def test_discovery_numerical_rejection_propagates_to_final_and_verifier(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            roots = make_three_node_fixture(root / "inputs")
            scale_root_fields_at_state(
                roots[0], topology="sketch", bias_V=-1.0,
                fields=("ElectrostaticPotential", "ElectricField"), factor=1.0e6,
            )
            out = root / "report"
            build_report_package(
                vela_root=roots[0], sentaurus_root=roots[1],
                supplemental_sentaurus_root=roots[2], out_dir=out,
                phase_base="a5524cf",
            )
            report_path = out / "physics_inverse_audit.json"
            report = json.loads(report_path.read_text(encoding="utf-8"))
            payload = report["payload"]
            metrics = {
                row["split"]: row["classification"]
                for row in payload["candidate_metrics"]
                if row["candidate"] == "potential_field_direct"
            }
            self.assertEqual(metrics, {
                "discovery": "rejected", "holdout": "identified",
                "combined": "identified",
            })
            classification = next(
                row for row in payload["classifications"]
                if row["candidate"] == "potential_field_direct"
            )
            self.assertEqual(classification["classification"], "rejected")

            classification["classification"] = "identified"
            classification["reason"] = (
                "discovery, holdout, and combined numerical gates passed without local fitting"
            )
            write_json(report_path, report)
            classifications_path = out / "candidate_classifications.json"
            classifications = json.loads(classifications_path.read_text(encoding="utf-8"))
            persisted = next(
                row for row in classifications["classifications"]
                if row["candidate"] == "potential_field_direct"
            )
            persisted.update(classification)
            write_json(classifications_path, classifications)
            refresh_integrity(out)

            with self.assertRaisesRegex(ValueError, "classification reconstruction mismatch"):
                verify_report(out)

    def test_verifier_rejects_partial_composite_promoted_to_identified(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            roots = make_three_node_fixture(root / "inputs")
            blank_root_field(roots[0], "ElectricField")
            out = root / "report"
            build_report_package(
                vela_root=roots[0], sentaurus_root=roots[1],
                supplemental_sentaurus_root=roots[2], out_dir=out,
                phase_base="a5524cf",
            )

            report_path = out / "physics_inverse_audit.json"
            report = json.loads(report_path.read_text(encoding="utf-8"))
            for row in report["payload"]["candidate_metrics"]:
                if row["candidate"] == "potential_field_direct":
                    row["classification"] = "identified"
            for row in report["payload"]["classifications"]:
                if row["candidate"] == "potential_field_direct":
                    row["classification"] = "identified"
                    row["reason"] = "combined and holdout numerical gates passed without local fitting"
            write_json(report_path, report)

            metrics_path = out / "candidate_metrics.csv"
            with metrics_path.open(newline="", encoding="utf-8") as handle:
                metric_rows = list(csv.reader(handle))
            candidate_column = metric_rows[0].index("candidate")
            classification_column = metric_rows[0].index("classification")
            for row in metric_rows[1:]:
                if row[candidate_column] == "potential_field_direct":
                    row[classification_column] = "identified"
            with metrics_path.open("w", newline="", encoding="utf-8") as handle:
                csv.writer(handle, lineterminator="\n").writerows(metric_rows)

            classifications_path = out / "candidate_classifications.json"
            classifications = json.loads(classifications_path.read_text(encoding="utf-8"))
            for row in classifications["classifications"]:
                if row["candidate"] == "potential_field_direct":
                    row["classification"] = "identified"
                    row["reason"] = "combined and holdout numerical gates passed without local fitting"
            write_json(classifications_path, classifications)
            refresh_integrity(out)

            with self.assertRaisesRegex(ValueError, "candidate metric reconstruction"):
                verify_report(out)

    def test_independent_verifier_rejects_simple_byte_mutations(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pristine, _ = self.build(root, "pristine")
            cases = []
            for label in ("csv", "json", "png", "input"):
                target = root / label
                shutil.copytree(pristine, target)
                cases.append((label, target))
            (cases[0][1] / "candidate_metrics.csv").write_bytes(
                (cases[0][1] / "candidate_metrics.csv").read_bytes() + b"mutation")
            with self.assertRaises(ValueError):
                verify_report(cases[0][1])

            classification_path = cases[1][1] / "candidate_classifications.json"
            classification_path.write_bytes(classification_path.read_bytes() + b" ")
            with self.assertRaises(ValueError):
                verify_report(cases[1][1])

            png_path = cases[2][1] / "figures" / "potential_field.png"
            with Image.open(png_path) as source:
                png = source.convert("RGB")
            original = png.getpixel((0, 0))
            png.putpixel((0, 0), (original[0] ^ 1, original[1], original[2]))
            png.save(png_path, format="PNG", optimize=False, compress_level=9,
                     dpi=(120, 120))
            with self.assertRaises(ValueError):
                verify_report(cases[2][1])

            report_manifest = json.loads((cases[3][1] / "report_manifest.json").read_text())
            raw_path = Path(report_manifest["inputs"][0]["path"])
            raw_path.write_bytes(raw_path.read_bytes() + b"mutation")
            with self.assertRaises(ValueError):
                verify_report(cases[3][1])

    def test_verifier_rejects_self_consistent_semantic_tampering_after_rehash(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pristine, _ = self.build(root, "pristine")
            cases = {}
            for label in ("observation", "candidate", "persisted_status"):
                target = root / label
                shutil.copytree(pristine, target)
                cases[label] = target

            observation = cases["observation"] / "observations_node.csv"
            with observation.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.reader(handle))
            rows[1][rows[0].index("value_si")] = "123456789"
            with observation.open("w", newline="", encoding="utf-8") as handle:
                csv.writer(handle, lineterminator="\n").writerows(rows)
            refresh_integrity(cases["observation"])

            candidate_csv = cases["candidate"] / "candidate_metrics.csv"
            with candidate_csv.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.reader(handle))
            rows[1][rows[0].index("valid_count")] = "999"
            with candidate_csv.open("w", newline="", encoding="utf-8") as handle:
                csv.writer(handle, lineterminator="\n").writerows(rows)
            report_path = cases["candidate"] / "physics_inverse_audit.json"
            report = json.loads(report_path.read_text())
            report["payload"]["candidate_metrics"][0]["valid_count"] = 999
            write_json(report_path, report)
            refresh_integrity(cases["candidate"])

            report_path = cases["persisted_status"] / "physics_inverse_audit.json"
            report = json.loads(report_path.read_text())
            report["payload"]["localization_control"]["semantic_replay"]["triangle_gradient"] = {
                "status": "insufficient_data"
            }
            write_json(report_path, report)
            refresh_integrity(cases["persisted_status"])

            for label, target in cases.items():
                with self.subTest(label=label), self.assertRaises(ValueError):
                    verify_report(target)


if __name__ == "__main__":
    unittest.main()
