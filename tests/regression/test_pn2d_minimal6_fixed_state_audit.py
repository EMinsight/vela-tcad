from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path
import shutil
import struct
import tempfile
import unittest

import fitz
from PIL import Image, ImageStat

REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "audit_pn2d_minimal6_fixed_state.py"
FIXTURE = REPO / "tests" / "fixtures" / "pn2d_minimal6_synthetic"
spec = importlib.util.spec_from_file_location("pn2d_minimal6_fixed_state_audit", SCRIPT)
if spec is None or spec.loader is None:
    raise RuntimeError(f"cannot load audit module from {SCRIPT}")
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


class FormulaReferenceTests(unittest.TestCase):
    def test_linear_triangle_gradient_matches_closed_form(self):
        self.assertEqual(audit.triangle_gradient([(0, 0), (1, 0), (0, 1)], [2, 5, 7]), (3.0, 5.0))

    def test_stable_bernoulli_limits(self):
        self.assertEqual(audit.bernoulli(0.0), 1.0)
        self.assertAlmostEqual(audit.bernoulli(1.0e-10), 1.0 - 0.5e-10)
        self.assertEqual(audit.bernoulli(-1000.0), 1000.0)
        self.assertEqual(audit.bernoulli(1000.0), 0.0)

    def test_electron_and_hole_sg_flux_signs_match_production(self):
        electron = audit.sg_electron_flux(2.0, 5.0, 0.0, 0.025, 3.0, 2.0)
        hole = audit.sg_hole_flux(2.0, 5.0, 0.0, 0.025, 3.0, 2.0)
        self.assertAlmostEqual(electron, 0.1125)
        self.assertAlmostEqual(hole, 0.1125)

    def test_gss_logistic_midpoint_uses_carrier_orientation(self):
        vt = 0.025852
        self.assertGreater(audit.gss_logistic_midpoint(1e12, 8e17, -0.2, 0.1, vt, "electron"), 7.9e17)
        self.assertGreater(audit.gss_logistic_midpoint(3e17, 2e11, -0.2, 0.1, vt, "hole"), 2.9e17)

    def test_canonical_projection_uses_endpoint_direction(self):
        self.assertAlmostEqual(audit.canonical_projection((3.0, 4.0), (0.0, 0.0), (2.0, 0.0)), 3.0)
        self.assertAlmostEqual(audit.canonical_projection((3.0, 4.0), (2.0, 0.0), (0.0, 0.0)), -3.0)

    def test_genius_truncated_partial_volume_for_right_triangle(self):
        points = [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0)]
        self.assertAlmostEqual(audit.genius_truncated_partial_volume(points, 0), 0.25)
        self.assertEqual(audit.genius_truncated_partial_volume(points, 1), 0.0)
        self.assertAlmostEqual(audit.genius_truncated_partial_volume(points, 2), 0.25)

        scaled_points = [(0.0, 0.0), (1.0e-6, 0.0), (1.0e-6, 0.5e-6)]
        self.assertEqual(
            audit.genius_truncated_partial_volume(scaled_points, 2), 0.0)

    def test_hybrid_error_and_both_zero_classification(self):
        self.assertEqual(audit.hybrid_error(0.0, 0.0), 0.0)
        self.assertEqual(audit.hybrid_error(2.0, 1.0), 0.5)
        ratio = audit.classify_orientation_pair(0.0, 0.0)
        self.assertEqual(ratio["zero_classification"], "both_zero")
        self.assertIsNone(ratio["absolute_log10_ratio"])


class ReviewContractTests(unittest.TestCase):
    def _copy(self):
        temp = tempfile.TemporaryDirectory()
        root = Path(temp.name) / "fixture"
        shutil.copytree(FIXTURE, root)
        self.addCleanup(temp.cleanup)
        return root

    @staticmethod
    def _manifest(root):
        return json.loads((root / "manifest.json").read_text(encoding="utf-8"))

    @staticmethod
    def _write(path, value):
        path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")

    @staticmethod
    def _export(root, state):
        path = Path(state["export_dir"])
        return path if path.is_absolute() else root / path

    @staticmethod
    def _mutate_csv(path, column, value, row=0):
        with path.open(encoding="utf-8", newline="") as source:
            rows = list(csv.DictReader(source))
        rows[row][column] = value
        with path.open("w", encoding="utf-8", newline="") as output:
            writer = csv.DictWriter(output, fieldnames=rows[0])
            writer.writeheader(); writer.writerows(rows)

    @staticmethod
    def _rewrite_live_topology_schema(export):
        elements_path = export / "elements.csv"
        with elements_path.open(encoding="utf-8", newline="") as source:
            elements = list(csv.DictReader(source))
        for row in elements:
            row["region"] = "R.Si"
            row["material"] = "Si"
        with elements_path.open("w", encoding="utf-8", newline="") as output:
            writer = csv.DictWriter(output, fieldnames=("id", "node0", "node1", "node2", "region", "material"))
            writer.writeheader(); writer.writerows(elements)

        contacts_path = export / "contacts.csv"
        with contacts_path.open(encoding="utf-8", newline="") as source:
            contacts = list(csv.DictReader(source))
        contacts = [{"name": row["name"], "node_ids": row["node_ids"], "region": "R.Si"} for row in contacts]
        with contacts_path.open("w", encoding="utf-8", newline="") as output:
            writer = csv.DictWriter(output, fieldnames=("name", "node_ids", "region"))
            writer.writeheader(); writer.writerows(contacts)


    @staticmethod
    def _producer():
        producer = REPO / "build-release" / "pn2d_minimal6_operator_audit.exe"
        if not producer.is_file():
            raise AssertionError(f"Task4 producer is required for provenance tests: {producer}")
        return producer

    def test_exact_variable_ni_sg_both_carriers(self):
        vt = 1.380649e-23 * 300.0 / 1.602176634e-19
        coef = 0.135 * vt / 1.25e-6
        args = (1.45e16, 2.2e16, -0.03, 0.07)
        self.assertAlmostEqual(audit.sg_electron_variable_ni_flux(*args, -0.021, 0.013, vt, coef), 9.085130342836316e19, delta=1e5)
        self.assertAlmostEqual(audit.sg_hole_variable_ni_flux(*args, 0.025, 0.082, vt, coef), -3.0985907840586436e20, delta=1e6)

    def test_actual_task3_root_and_field_manifest_list(self):
        manifest = self._manifest(FIXTURE)
        self.assertEqual(manifest["schema"], "vela.pn2d_minimal6_states.v1")
        self.assertEqual(len(manifest["states"]), 6)
        export = self._export(FIXTURE, manifest["states"][0])
        fields = json.loads((export / "field_manifest.json").read_text(encoding="utf-8"))["fields"]
        self.assertIsInstance(fields, list)
        report = audit.build_report(FIXTURE)
        self.assertEqual((len(report.node_rows), len(report.edge_rows), len(report.triangle_rows)), (36,54,24))

    def test_requested_actual_and_optional_field_bias_mismatch(self):
        root = self._copy(); manifest = self._manifest(root)
        manifest["states"][0]["actual_bias_V"] = 1e-15
        self._write(root / "manifest.json", manifest)
        with self.assertRaisesRegex(audit.ContractError, "requested.*actual bias"):
            audit.build_report(root)
        root = self._copy(); manifest = self._manifest(root); state = manifest["states"][0]
        field_path = self._export(root, state) / "field_manifest.json"
        field = json.loads(field_path.read_text(encoding="utf-8")); field["bias_V"] = -1.0
        self._write(field_path, field)
        with self.assertRaisesRegex(audit.ContractError, "field manifest bias"):
            audit.build_report(root)

    def test_noncompensated_doping_and_local_edge_mapping_fail(self):
        root = self._copy(); state = self._manifest(root)["states"][0]; export = self._export(root,state)
        self._mutate_csv(export / "doping.csv", "donors_cm3", "1e17", 0)
        with self.assertRaisesRegex(audit.ContractError, "doping semantics"):
            audit.build_report(root)
        root = self._copy(); state = self._manifest(root)["states"][0]; export = self._export(root,state)
        self._mutate_csv(export / "vela_triangle_audit.csv", "local_edge0_node0", "5", 0)
        with self.assertRaisesRegex(audit.ContractError, "local edge"):
            audit.build_report(root)

    def test_wide_schema_orientation_and_abs_log_ratio(self):
        self.assertAlmostEqual(audit.classify_orientation_pair(10.0,1.0)["absolute_log10_ratio"], 1.0)
        report = audit.build_report(FIXTURE)
        for key in ("raw_eDensity_cm3","sentaurus_n_m3","vela_n_m3","abs_error_n_m3","hybrid_error_n_m3"):
            self.assertIn(key, report.node_rows[0])
        for key in ("node0_psi_V","electron_mobility_m2_per_V_s","python_electron_flux_per_m2_s","vela_electron_current_A_per_m2","electron_formula_hybrid_error"):
            self.assertIn(key, report.edge_rows[0])
        for key in ("shape_grad_N0_x_per_m","vela_grad_psi_x_V_per_m","reconstructed_electron_current_x_A_per_m2","sentaurus_vs_vela_total_source_diagnostic"):
            self.assertIn(key, report.triangle_rows[0])
        self.assertGreater(len({x["quantity"] for x in report.summary["orientation_sensitivity"]}), 51)

    def test_cpp_provenance_replay_and_no_python_vela_producer(self):
        manifest = self._manifest(FIXTURE); provenance = manifest["task4_provenance"]
        self.assertEqual(provenance["producer"], "build-release/pn2d_minimal6_operator_audit.exe")
        self.assertEqual(provenance["task4_source_commit"], "dfe2611975742779e6d27e164d3b695a5d189e44")
        self.assertEqual(len(provenance["replays"]), 6)
        self.assertEqual(audit.verify_task4_replay(FIXTURE, self._producer()), [])
        source = SCRIPT.read_text(encoding="utf-8")
        self.assertNotIn("make_synthetic_fixture", source)
        self.assertNotIn("--make-synthetic-fixture", source)

    def test_cpp_producer_has_reproducible_pe_timestamp(self):
        data = self._producer().read_bytes()
        pe_offset = struct.unpack_from("<I", data, 0x3C)[0]
        self.assertEqual(data[pe_offset : pe_offset + 4], b"PE\0\0")
        timestamp = struct.unpack_from("<I", data, pe_offset + 8)[0]
        self.assertEqual(timestamp, 0)

    def test_fourteen_figures_and_numeric_summary(self):
        report = audit.build_report(FIXTURE)
        with tempfile.TemporaryDirectory() as directory:
            out = Path(directory); audit.write_report(report,out)
            expected = {
                "minimal6-topologies.png", "minimal6-topologies.pdf",
                *(f"minimal6-edge-current-audit-{bias}.{suffix}" for bias in ("0v", "minus12v", "minus19v") for suffix in ("png", "pdf")),
                *(f"minimal6-triangle-source-audit-{bias}.{suffix}" for bias in ("0v", "minus12v", "minus19v") for suffix in ("png", "pdf")),
            }
            figures = {x.name: x for x in (out / "figures").iterdir()}
            self.assertEqual(set(figures), expected)
            for name, path in figures.items():
                self.assertGreater(path.stat().st_size, 0, name)
                if path.suffix == ".pdf":
                    self.assertTrue(path.read_bytes().startswith(b"%PDF-"), name)
                    with fitz.open(filename=str(path)) as document:
                        self.assertEqual(document.page_count, 1, name)
                        pixmap = document[0].get_pixmap(colorspace=fitz.csGRAY, alpha=False)
                        pixels = memoryview(pixmap.samples)
                        self.assertGreater(max(pixels) - min(pixels), 1, name)
                else:
                    self.assertTrue(path.read_bytes().startswith(b"\x89PNG\r\n\x1a\n"), name)
                    with Image.open(path) as image:
                        image.verify()
                    with Image.open(path).convert("L") as image:
                        self.assertGreater(ImageStat.Stat(image).stddev[0], 1.0, name)
            md=(out/"summary.md").read_text(encoding="utf-8")
            self.assertIn("Maximum state parity hybrid error:",md)
            self.assertIn("Maximum C++/Python formula hybrid error:",md)
            self.assertIn("fixed-state operator audit, not a BV curve", md)
            self.assertNotIn("committed state root is synthetic", md)
            self.assertEqual(json.loads((out/"summary.json").read_text(encoding="utf-8"))["status"], "PASS")
    def test_complete_matrix_unique_keys_and_formula_gate(self):
        report = audit.build_report(FIXTURE)
        self.assertEqual((len(report.node_rows),len(report.edge_rows),len(report.triangle_rows)),(36,54,24))
        audit.require_unique_keys(report.edge_rows,("topology_id","bias_V","node0","node1"))
        self.assertLess(report.summary["gates"]["max_cpp_python_formula_hybrid_error"],5e-12)
        self.assertIsNone(report.summary["gates"]["sentaurus_vs_vela_current_source_threshold"])

    def test_partial_state_fails_closed(self):
        root=self._copy(); manifest=self._manifest(root); state=manifest["states"][0]
        path=root/state["state_csv"]; lines=path.read_text(encoding="utf-8").splitlines(); path.write_text("\n".join(lines[:-1])+"\n",encoding="utf-8")
        with self.assertRaisesRegex(audit.ContractError,"partial state"): audit.build_report(root)

    def test_missing_required_field_fails_closed(self):
        root=self._copy(); state=self._manifest(root)["states"][0]; path=self._export(root,state)/"field_manifest.json"
        value=json.loads(path.read_text(encoding="utf-8")); value["fields"]=[x for x in value["fields"] if x["name"]!="eMobility"]; self._write(path,value)
        with self.assertRaisesRegex(audit.ContractError,"missing required field"): audit.build_report(root)

    def test_duplicate_edge_key_fails_closed(self):
        root=self._copy(); state=self._manifest(root)["states"][0]; path=self._export(root,state)/"vela_edge_audit.csv"
        lines=path.read_text(encoding="utf-8").splitlines(); path.write_text("\n".join(lines+[lines[1]])+"\n",encoding="utf-8")
        with self.assertRaisesRegex(audit.ContractError,"row count|duplicate"): audit.build_report(root)

    def test_inexact_requested_bias_fails_closed(self):
        root=self._copy(); manifest=self._manifest(root); manifest["states"][0]["requested_bias_V"]=-11.999999999; manifest["states"][0]["actual_bias_V"]=-11.999999999; self._write(root/"manifest.json",manifest)
        with self.assertRaisesRegex(audit.ContractError,"unexpected topology, bias"): audit.build_report(root)

    def test_wrong_vector_components_and_units_fail_closed(self):
        for key,value,pattern in (("components",1,"missing required field"),("unit","V/cm","wrong unit")):
            root=self._copy(); state=self._manifest(root)["states"][0]; path=self._export(root,state)/"field_manifest.json"; manifest=json.loads(path.read_text(encoding="utf-8"))
            field=next(x for x in manifest["fields"] if x["name"]=="ElectricField" and x["components"]==2); field[key]=value; self._write(path,manifest)
            with self.subTest(key=key),self.assertRaisesRegex(audit.ContractError,pattern): audit.build_report(root)

    def test_changed_topology_fails_closed(self):
        root=self._copy(); state=self._manifest(root)["states"][0]; path=self._export(root,state)/"elements.csv"; self._mutate_csv(path,"node2","5",0)
        with self.assertRaisesRegex(audit.ContractError,"topology|triangle orientation"): audit.build_report(root)

    def test_live_task3_topology_csv_schema_is_strictly_supported(self):
        root = self._copy(); state = self._manifest(root)["states"][0]; export = self._export(root, state)
        self._rewrite_live_topology_schema(export)
        audit.build_report(root)

        root = self._copy(); state = self._manifest(root)["states"][0]; export = self._export(root, state)
        self._rewrite_live_topology_schema(export); self._mutate_csv(export / "elements.csv", "region", "R.Other")
        with self.assertRaisesRegex(audit.ContractError, "element region"):
            audit.build_report(root)

        root = self._copy(); state = self._manifest(root)["states"][0]; export = self._export(root, state)
        self._rewrite_live_topology_schema(export); self._mutate_csv(export / "elements.csv", "material", "Ge")
        with self.assertRaisesRegex(audit.ContractError, "element material"):
            audit.build_report(root)

        root = self._copy(); state = self._manifest(root)["states"][0]; export = self._export(root, state)
        self._rewrite_live_topology_schema(export); self._mutate_csv(export / "contacts.csv", "region", "R.Other")
        with self.assertRaisesRegex(audit.ContractError, "contact region"):
            audit.build_report(root)

    def test_live_task3_doping_roundoff_is_tolerated_but_drift_fails(self):
        root = self._copy(); state = self._manifest(root)["states"][0]; export = self._export(root, state)
        rounded = "1.0000000000000002e17"
        self._mutate_csv(export / "doping.csv", "donors_cm3", rounded, 1)
        self._mutate_csv(export / "doping.csv", "acceptors_cm3", rounded, 1)
        self._mutate_csv(export / "fields" / "DonorConcentration_region0.csv", "component0", rounded, 1)
        self._mutate_csv(export / "fields" / "AcceptorConcentration_region0.csv", "component0", rounded, 1)
        audit.build_report(root)

        root = self._copy(); state = self._manifest(root)["states"][0]; export = self._export(root, state)
        drifted = "1.00000001e17"
        self._mutate_csv(export / "doping.csv", "donors_cm3", drifted, 1)
        self._mutate_csv(export / "fields" / "DonorConcentration_region0.csv", "component0", drifted, 1)
        with self.assertRaisesRegex(audit.ContractError, "wrong doping semantics"):
            audit.build_report(root)

    def test_nonfinite_state_fails_closed(self):
        root=self._copy(); state=self._manifest(root)["states"][0]; self._mutate_csv(root/state["state_csv"],"psi_V","nan",0)
        with self.assertRaisesRegex(audit.ContractError,"non-finite"): audit.build_report(root)

    def test_near_zero_triangle_gradient_uses_absolute_formula_tolerance(self):
        root = self._copy(); state = self._manifest(root)["states"][0]; path = self._export(root, state) / "vela_triangle_audit.csv"
        self._mutate_csv(path, "grad_psi_x_V_per_m", "2e-9", 0)
        audit.build_report(root)

        root = self._copy(); state = self._manifest(root)["states"][0]; path = self._export(root, state) / "vela_triangle_audit.csv"
        self._mutate_csv(path, "grad_psi_x_V_per_m", "6e-9", 0)
        with self.assertRaisesRegex(audit.ContractError, "formula error"):
            audit.build_report(root)

    def test_formula_error_above_threshold_fails_closed(self):
        root=self._copy(); state=self._manifest(root)["states"][0]; path=self._export(root,state)/"vela_edge_audit.csv"
        with path.open(encoding="utf-8",newline="") as source: rows=list(csv.DictReader(source))
        rows[0]["electron_raw_signed_flux_per_m2_s"]=format(float(rows[0]["electron_raw_signed_flux_per_m2_s"])*(1+5.1e-12),".17g")
        with path.open("w",encoding="utf-8",newline="") as output: writer=csv.DictWriter(output,fieldnames=rows[0]); writer.writeheader(); writer.writerows(rows)
        with self.assertRaisesRegex(audit.ContractError,"formula error"): audit.build_report(root)

    def test_state_error_above_threshold_fails_closed(self):
        root=self._copy(); state=self._manifest(root)["states"][0]; path=self._export(root,state)/"vela_node_state.csv"
        with path.open(encoding="utf-8",newline="") as source: rows=list(csv.DictReader(source))
        rows[0]["n_m3"]=format(float(rows[0]["n_m3"])*(1+1.1e-12),".17g")
        with path.open("w",encoding="utf-8",newline="") as output: writer=csv.DictWriter(output,fieldnames=rows[0]); writer.writeheader(); writer.writerows(rows)
        with self.assertRaisesRegex(audit.ContractError,"state parity"): audit.build_report(root)

    def test_raw_unit_conversions_and_normalized_zero_are_visible(self):
        report=audit.build_report(FIXTURE); node=report.node_rows[0]
        self.assertEqual(node["sentaurus_n_m3"],node["raw_eDensity_cm3"]*1e6)
        self.assertEqual(node["sentaurus_electric_field_x_V_per_m"],node["raw_ElectricField_x_V_per_cm"]*100)
        self.assertTrue(any(row["normalized_geometric_zero"] for row in report.triangle_rows))

    def test_provenance_missing_duplicate_and_forged_records_fail_closed(self):
        cases = []
        root = self._copy(); manifest = self._manifest(root)
        manifest["task4_provenance"]["replays"].pop(); self._write(root/"manifest.json", manifest)
        cases.append((root, "six exact replay identities"))
        root = self._copy(); manifest = self._manifest(root)
        manifest["task4_provenance"]["replays"][1] = dict(manifest["task4_provenance"]["replays"][0]); self._write(root/"manifest.json", manifest)
        cases.append((root, "six exact replay identities"))
        root = self._copy(); manifest = self._manifest(root)
        manifest["task4_provenance"]["replays"][0]["producer"] = "forged.exe"; self._write(root/"manifest.json", manifest)
        cases.append((root, "replay producer identity"))
        for root, expected in cases:
            with self.subTest(expected=expected):
                failures = audit.verify_task4_replay(root, self._producer())
                self.assertTrue(any(expected in failure for failure in failures), failures)

    def test_provenance_committed_and_fresh_output_hash_tamper_fail_closed(self):
        root = self._copy(); manifest = self._manifest(root); replay = manifest["task4_provenance"]["replays"][0]
        committed = root / replay["arguments"][replay["arguments"].index("--edge-out") + 1]
        committed.write_bytes(committed.read_bytes() + b"\n")
        failures = audit.verify_task4_replay(root, self._producer())
        self.assertTrue(any("committed output hash mismatch" in failure for failure in failures), failures)

        root = self._copy(); manifest = self._manifest(root); replay = manifest["task4_provenance"]["replays"][0]
        key = replay["arguments"][replay["arguments"].index("--edge-out") + 1]
        replay["output_sha256"][key] = "0" * 64; self._write(root/"manifest.json", manifest)
        failures = audit.verify_task4_replay(root, self._producer())
        self.assertTrue(any("fresh replay output hash mismatch" in failure for failure in failures), failures)

    def test_cli_report_generation_is_blocked_before_artifacts_on_bad_provenance(self):
        root = self._copy(); manifest = self._manifest(root)
        manifest["task4_provenance"]["producer_sha256"] = "0" * 64; self._write(root/"manifest.json", manifest)
        report = audit.build_report(root)
        with tempfile.TemporaryDirectory() as directory:
            out = Path(directory)/"must_not_exist"
            with self.assertRaisesRegex(audit.ContractError, "producer hash mismatch"):
                audit.write_report(report, out)
            self.assertFalse(out.exists())

    def test_vela_and_python_triangle_aggregates_are_separate_and_gated(self):
        report = audit.build_report(FIXTURE)
        for row in report.triangle_rows:
            raw_e = sum(row[f"vela_local_edge{i}_electron_source_integral_per_m_s"] for i in range(3))
            raw_h = sum(row[f"vela_local_edge{i}_hole_source_integral_per_m_s"] for i in range(3))
            py_e = sum(row[f"python_local_edge{i}_electron_source_integral_per_m_s"] for i in range(3))
            py_h = sum(row[f"python_local_edge{i}_hole_source_integral_per_m_s"] for i in range(3))
            self.assertLess(
                audit.hybrid_error(row["vela_electron_source_integral_per_m_s"], raw_e),
                audit.FORMULA_LIMIT,
            )
            self.assertLess(
                audit.hybrid_error(row["vela_hole_source_integral_per_m_s"], raw_h),
                audit.FORMULA_LIMIT,
            )
            self.assertLess(
                audit.hybrid_error(row["python_electron_source_integral_per_m_s"], py_e),
                audit.FORMULA_LIMIT,
            )
            self.assertLess(
                audit.hybrid_error(row["python_hole_source_integral_per_m_s"], py_h),
                audit.FORMULA_LIMIT,
            )
            self.assertLess(row["vela_vs_python_total_source_hybrid_error"], audit.FORMULA_LIMIT)

    def test_elements_accept_approved_ccw_tuple_set_in_any_id_order(self):
        root = self._copy()
        state = self._manifest(root)["states"][0]
        path = self._export(root, state) / "elements.csv"
        with path.open(encoding="utf-8", newline="") as source:
            rows = list(csv.DictReader(source))
        rows[0]["id"], rows[1]["id"] = rows[1]["id"], rows[0]["id"]
        with path.open("w", encoding="utf-8", newline="") as output:
            writer = csv.DictWriter(output, fieldnames=rows[0])
            writer.writeheader()
            writer.writerows(rows)

        report = audit.build_report(root)
        actual = [
            (row["node0"], row["node1"], row["node2"])
            for row in report.triangle_rows
            if row["topology_id"] == "sketch" and row["bias_V"] == 0.0
        ]
        self.assertEqual(actual, list(audit.TRIS["sketch"]))

    def test_nonzero_partial_volume_cannot_normalize_tiny_source_residue(self):
        with self.assertRaisesRegex(audit.ContractError, "formula error"):
            audit.geometric_source_gate(1e-290, 2e-290, 1e-14, 1e-14, audit.FORMULA_LIMIT, "formula error")
        self.assertEqual(audit.geometric_source_gate(1e-290, 2e-290, 0.0, 1e-30, audit.FORMULA_LIMIT, "formula error"), 0.0)

    def test_corrupted_cpp_mobility_and_alpha_fail_independent_gates(self):
        for column, expected in (("local_edge0_electron_mobility_m2_per_V_s", "mobility"), ("local_edge0_electron_alpha_per_m", "alpha")):
            root = self._copy(); state = self._manifest(root)["states"][0]; path = self._export(root,state)/"vela_triangle_audit.csv"
            with path.open(encoding="utf-8",newline="") as source: rows=list(csv.DictReader(source))
            rows[0][column] = format(float(rows[0][column]) * 1.01 + 1e-30, ".17g")
            with path.open("w",encoding="utf-8",newline="") as output: writer=csv.DictWriter(output,fieldnames=rows[0]); writer.writeheader(); writer.writerows(rows)
            with self.subTest(column=column), self.assertRaisesRegex(audit.ContractError, expected):
                audit.build_report(root)

    def test_raw_task3_doping_and_inline_triangle_order_fail_closed(self):
        root = self._copy(); state = self._manifest(root)["states"][0]; export = self._export(root,state)
        self._mutate_csv(export/"fields"/"DonorConcentration_region0.csv", "component0", "9e16", 1)
        with self.assertRaisesRegex(audit.ContractError, "raw Task3 donor"):
            audit.build_report(root)

        root = self._copy(); manifest = self._manifest(root)
        triangle = manifest["states"][0]["topology_contract"]["triangle_connectivity"][0]
        manifest["states"][0]["topology_contract"]["triangle_connectivity"][0] = [triangle[0], triangle[2], triangle[1]]
        self._write(root/"manifest.json", manifest)
        with self.assertRaisesRegex(audit.ContractError, "exact approved CCW topology"):
            audit.build_report(root)

if __name__ == "__main__":
    unittest.main()