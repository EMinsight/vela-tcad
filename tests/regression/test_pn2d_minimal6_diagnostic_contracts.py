import tempfile
import unittest
from scripts.pn2d_minimal6_diagnostics.contracts import BranchKind, QuantityRecord, SourceKind, StateIdentity, SupportKind, classify_pair
from scripts.pn2d_minimal6_diagnostics.ledger import DiagnosticLedger
from scripts.pn2d_minimal6_diagnostics.units import convert_value

class DiagnosticContractsTest(unittest.TestCase):
    def test_converts_supported_units_and_rejects_mismatch(self):
        self.assertEqual(convert_value(1.0, "cm^-3", "m^-3"), 1.0e6)
        self.assertEqual(convert_value(1.0, "A/cm^2", "A/m^2"), 1.0e4)
        self.assertEqual(convert_value(1.0, "cm^2/(V s)", "m^2/(V s)"), 1.0e-4)
        self.assertEqual(convert_value(1.0, "cm^-1", "m^-1"), 1.0e2)
        self.assertEqual(convert_value(1.0, "cm^-3*s^-1", "m^-3*s^-1"), 1.0e6)
        with self.assertRaises(ValueError): convert_value(1.0, "V/cm", "A/cm^2")
    def test_ledger_key_is_unique_and_branch_is_typed(self):
        record = QuantityRecord(StateIdentity("run", "sketch", -19.0), "electron", SupportKind.NODE, "0", "G", SourceKind.VELA, "v1", 1.0, "cm^-3*s^-1", "positive", "fixture", "0000000000000000000000000000000000000000000000000000000000000000")
        ledger = DiagnosticLedger(); ledger.add(record)
        with self.assertRaises(ValueError): ledger.add(record)
        self.assertEqual(classify_pair(1e-12, 2e-12), BranchKind.LEAKAGE_LIKE)
        self.assertEqual(classify_pair(1e-8, 2e-12), BranchKind.MULTIPLICATION_LIKE)

    def test_record_requires_provenance_and_sign_convention(self):
        with self.assertRaises(ValueError):
            QuantityRecord(StateIdentity("run", "sketch", -19.0), "electron", SupportKind.NODE, "0", "G", SourceKind.VELA, "v1", 1.0, "cm^-3*s^-1", "", "", "")
    def test_report_schema_rejects_nonfinite_and_missing_disclaimer(self):
        from scripts.pn2d_minimal6_diagnostics.schemas import validate_formula_difference_v1
        report = {"schema":"vela.pn2d_minimal6_formula_difference.v1", "diagnostic_disclaimer":"minimal6 diagnostic sweep; not a physical BV curve", "records":[], "input_provenance":{}, "audit_provenance":{}, "state_matrix":[], "row_counts":{}, "waterfall_paths":[], "interactions":[], "dominance_rules":{}, "sentaurus_internal_semantics_residual":0.0, "vela_parameter_agreement":[], "artifact_hashes":{}}
        self.assertIsNone(validate_formula_difference_v1(report))
        with self.assertRaises(ValueError): validate_formula_difference_v1({"schema":report["schema"], "records":[{"value":float("nan")} ]})
        with self.assertRaises(ValueError): validate_formula_difference_v1({"schema":report["schema"], "diagnostic_disclaimer":report["diagnostic_disclaimer"]})

    def test_ledger_json_is_deterministic_and_forbids_nan(self):
        from pathlib import Path
        ledger = DiagnosticLedger()
        ledger.add(QuantityRecord(StateIdentity("run", "mirror", -12.0), "hole", SupportKind.NODE, "1", "G", SourceKind.VELA, "v1", 1.0, "cm^-3*s^-1", "positive", "fixture", "0000000000000000000000000000000000000000000000000000000000000000"))
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "ledger.json"
            ledger.write_json(path)
            first = path.read_bytes(); ledger.write_json(path)
            self.assertEqual(first, path.read_bytes())
            ledger.write_csv(path.with_suffix(".csv"))
            self.assertIn("raw_source_sha256", path.with_suffix(".csv").read_text(encoding="utf-8"))
            self.assertIn("raw_source_path", path.read_text(encoding="utf-8"))

if __name__ == "__main__": unittest.main()
