import copy
import csv
import json
import math
import tempfile
import unittest
from dataclasses import FrozenInstanceError
from pathlib import Path

from scripts.pn2d_minimal6_diagnostics import schemas
from scripts.pn2d_minimal6_diagnostics.contracts import (
    BranchKind, QuantityRecord, SourceKind, StateIdentity, SupportKind, classify_pair,
)
from scripts.pn2d_minimal6_diagnostics.ledger import DiagnosticLedger
from scripts.pn2d_minimal6_diagnostics.units import convert_value

DISCLAIMER = "minimal6 diagnostic sweep; not a physical BV curve"
SHA_A, SHA_B = "a" * 64, "b" * 64
THRESHOLD = "v1: multiplication=[0.1,10], leakage<=1e-3"


def state_matrix():
    return [
        {"topology_id": topology, "requested_bias_V": bias,
         "actual_bias_V": bias, "status": "passed"}
        for topology in ("sketch", "mirror") for bias in (0.0, -12.0, -19.0)
    ]


def formula_report():
    paths, residuals = [], []
    for state in state_matrix():
        topology, bias = state["topology_id"], state["requested_bias_V"]
        paths.append({
            "topology": topology, "bias_V": bias, "dependency_order": ["mobility"],
            "forward": {"order": ["mobility"], "contributions": []},
            "reverse": {"order": ["mobility"], "contributions": []},
            "interactions": [], "native_gap_dex": 0.0, "residual_dex": 0.0,
            "status": "insufficient_data",
        })
        residuals.append({"topology": topology, "bias_V": bias,
                          "classification": "available", "dex": 0.0})
    return {
        "schema": "vela.pn2d_minimal6_formula_difference.v1",
        "diagnostic_disclaimer": DISCLAIMER,
        "input_provenance": {"state_manifest": "fixture/manifest.json"},
        "audit_provenance": {"audit_root": "fixture/audit"},
        "state_matrix": state_matrix(),
        "row_counts": {"node": 36, "edge": 54, "triangle": 24},
        "waterfall_paths": paths, "interactions": [],
        "dominance_rules": {"status": "insufficient_data"},
        "sentaurus_internal_semantics_residual": residuals,
        "vela_parameter_agreement": [],
        "artifact_hashes": {"state_manifest_sha256": SHA_A}, "records": [],
    }


def transition(solver, topology="sketch", bias=-1.0):
    return {
        "solver": solver, "topology": topology, "start_bias_V": bias + 1.0,
        "target_bias_V": bias, "actual_bias_V": bias, "exit_code": 0,
        "status": "accepted", "state_path": f"{solver}/{topology}/state.json",
        "state_sha256": SHA_A if solver == "vela" else SHA_B,
        "observables": {
            "anode_current_A_per_um": -1.0e-8,
            "cathode_current_A_per_um": 1.0e-8,
            "max_field_V_per_m": 2.0e7,
            "native_source_integral_s_inv_per_cm": 3.0,
            "reconstructed_source_integral_s_inv_per_cm": 2.5,
        },
        "branch_classification": "multiplication_like",
        "branch_threshold_version": THRESHOLD,
        "convergence_metadata": {"iterations": 4, "residual_norm": 1.0e-10},
        "stdout": "", "stderr": "",
    }


def comparison_report():
    vela, sentaurus = transition("vela"), transition("sentaurus")
    checkpoint = {
        "topology": "sketch", "bias_V": -1.0, "classification": "common_exact",
        "vela": copy.deepcopy(vela), "sentaurus": copy.deepcopy(sentaurus),
        "branch_classification": "multiplication_like",
        "branch_threshold_version": THRESHOLD,
    }
    return {
        "schema": "vela.pn2d_minimal6_bv_comparison.v1",
        "diagnostic_disclaimer": DISCLAIMER, "interpolation": "forbidden",
        "branch_threshold_version": THRESHOLD,
        "solver_configurations": {"vela": {}, "sentaurus": {}},
        "accepted_transitions": {"vela": [vela], "sentaurus": [sentaurus]},
        "failed_transitions": [], "failure_transitions": [],
        "checkpoints": [checkpoint], "records": [checkpoint],
        "terminal_currents": [checkpoint], "maximum_fields": [checkpoint],
        "source_integrals": [checkpoint],
        "convergence_metadata": {"vela_accepted": 1, "sentaurus_accepted": 1,
                                 "common_exact": 1},
        "curve_artifact_hashes": {"vela_manifest": SHA_A, "sentaurus_manifest": SHA_B},
        "deepest_common_bias_V": {"classification": "available", "value": -1.0},
        "missing_tails": [], "topology_sensitivity": [], "fixed_state_recheck": [],
        "artifact_hashes": {}, "input_artifacts": {}, "closure": {"status": "closed"},
    }


def sweep_manifest():
    accepted = transition("vela")
    return {
        "schema": "vela.pn2d_minimal6_sweep_manifest.v1",
        "diagnostic_disclaimer": DISCLAIMER, "targets_V": [0.0, -1.0],
        "template": {"path": "fixture/template.json", "sha256": SHA_A},
        "topology_input_sha256": {"sketch": {"mesh.json": SHA_A, "doping.csv": SHA_B}},
        "segments": [], "sentaurus_segments": [], "accepted_checkpoints": [accepted],
        "failed_transition": None, "failed_transitions": [], "interpolation": "forbidden",
        "branch_threshold_version": THRESHOLD,
    }


def record(*, bias=-2.0, support_id="2", value=1.0):
    return QuantityRecord(
        StateIdentity("run", "sketch", bias), "electron", SupportKind.NODE, support_id,
        "avalanche_generation", SourceKind.VELA, "v1", value, "cm^-3*s^-1",
        "positive means carrier-pair generation", "fixture/raw.csv", SHA_A,
    )


class DiagnosticContractsTest(unittest.TestCase):
    def test_converts_all_supported_si_quantities(self):
        conversions = (
            ("V/cm", "V/m", 1.0e2), ("A/cm^2", "A/m^2", 1.0e4),
            ("cm^-3", "m^-3", 1.0e6), ("cm^2/(V s)", "m^2/(V s)", 1.0e-4),
            ("cm^-1", "m^-1", 1.0e2), ("cm/s", "m/s", 1.0e-2),
            ("cm^-3*s^-1", "m^-3*s^-1", 1.0e6),
        )
        for source, target, expected in conversions:
            with self.subTest(source=source, target=target):
                self.assertEqual(convert_value(1.0, source, target), expected)
                self.assertAlmostEqual(convert_value(expected, target, source), 1.0)

    def test_rejects_mismatch_unknown_units_and_nonfinite_values(self):
        for source, target in (("V/cm", "A/cm^2"), ("widget", "widget")):
            with self.subTest(source=source, target=target), self.assertRaises(ValueError):
                convert_value(1.0, source, target)
        with self.assertRaises(ValueError):
            convert_value(math.inf, "V/cm", "V/m")

    def test_enums_and_branch_zero_classification_are_typed(self):
        self.assertEqual({item.value for item in SourceKind}, {"sentaurus", "vela", "derived"})
        self.assertEqual(SourceKind.SENTAURUS.value, "sentaurus")
        self.assertEqual({item.value for item in SupportKind}, {"node", "edge", "cell", "terminal"})
        self.assertEqual({item.value for item in BranchKind},
                         {"leakage_like", "multiplication_like", "unidentified"})
        self.assertEqual(classify_pair(1.0e-12, 2.0e-12), BranchKind.LEAKAGE_LIKE)
        self.assertEqual(classify_pair(1.0e-8, 2.0e-12), BranchKind.MULTIPLICATION_LIKE)
        self.assertEqual(classify_pair(0.0, 2.0e-12), BranchKind.UNIDENTIFIED)
        self.assertEqual(classify_pair(math.nan, 2.0e-12), BranchKind.UNIDENTIFIED)

    def test_records_are_immutable_and_validate_identity(self):
        item = record()
        with self.assertRaises(FrozenInstanceError):
            item.value = 2.0
        with self.assertRaises(ValueError):
            StateIdentity("run", "sketch", math.nan)
        values = [item.state, item.carrier, "node", item.support_id, item.quantity,
                  item.source, item.formula_version, item.value, item.unit,
                  item.sign_convention, item.raw_source_path, item.raw_source_sha256]
        with self.assertRaises(TypeError):
            QuantityRecord(*values)
        values[2], values[-1] = item.support_kind, "z" * 64
        with self.assertRaises(ValueError):
            QuantityRecord(*values)

    def test_ledger_key_is_unique_and_sort_order_is_canonical(self):
        items = [record(bias=-2.0, support_id="10"), record(bias=-12.0, support_id="10"),
                 record(bias=-2.0, support_id="2")]
        ledger = DiagnosticLedger()
        for item in reversed(items):
            ledger.add(item)
        with self.assertRaises(ValueError):
            ledger.add(items[0])
        self.assertEqual([(row.state.bias_V, row.support_id) for row in ledger.records()],
                         [(-12.0, "10"), (-2.0, "2"), (-2.0, "10")])

    def test_serialization_is_deterministic_and_rejects_nonfinite(self):
        items = [record(bias=-2.0, support_id="10"), record(bias=-12.0, support_id="10"),
                 record(bias=-2.0, support_id="2")]
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            first, second = DiagnosticLedger(), DiagnosticLedger()
            for item in items: first.add(item)
            for item in reversed(items): second.add(item)
            for ledger, stem in ((first, "first"), (second, "second")):
                ledger.write_json(root / f"{stem}.json")
                ledger.write_csv(root / f"{stem}.csv")
            self.assertEqual((root / "first.json").read_bytes(), (root / "second.json").read_bytes())
            self.assertEqual((root / "first.csv").read_bytes(), (root / "second.csv").read_bytes())
            with (root / "first.csv").open(newline="", encoding="utf-8") as handle:
                self.assertEqual([(float(row["bias_V"]), row["support_id"])
                                  for row in csv.DictReader(handle)],
                                 [(-12.0, "10"), (-2.0, "2"), (-2.0, "10")])
            corrupted, ledger = record(), DiagnosticLedger()
            ledger.add(corrupted)
            object.__setattr__(corrupted, "value", math.nan)
            with self.assertRaises(ValueError): ledger.write_json(root / "nan.json")
            with self.assertRaises(ValueError): ledger.write_csv(root / "nan.csv")

    def test_all_three_schema_files_declare_their_contracts(self):
        root = Path(__file__).parents[2] / "schemas"
        names = ("vela.pn2d_minimal6_formula_difference.v1",
                 "vela.pn2d_minimal6_bv_comparison.v1",
                 "vela.pn2d_minimal6_sweep_manifest.v1")
        for name in names:
            with self.subTest(name=name):
                payload = json.loads((root / f"{name}.schema.json").read_text(encoding="utf-8"))
                self.assertEqual(payload["$schema"], "https://json-schema.org/draft/2020-12/schema")
                self.assertEqual(payload["title"], name)
                self.assertIn("diagnostic_disclaimer", payload["required"])

    def test_formula_schema_golden_and_invalid_fixture(self):
        report = formula_report()
        self.assertIsNone(schemas.validate_formula_difference_v1(report))
        invalid = copy.deepcopy(report)
        invalid["state_matrix"][-1] = copy.deepcopy(invalid["state_matrix"][0])
        with self.assertRaises(ValueError): schemas.validate_formula_difference_v1(invalid)

    def test_comparison_schema_golden_and_invalid_fixture(self):
        report = comparison_report()
        self.assertIsNone(schemas.validate_bv_comparison_v1(report))
        invalid = copy.deepcopy(report)
        del invalid["accepted_transitions"]["vela"][0]["observables"]["cathode_current_A_per_um"]
        with self.assertRaises(ValueError): schemas.validate_bv_comparison_v1(invalid)

    def test_sweep_manifest_validator_is_exported(self):
        self.assertTrue(hasattr(schemas, "validate_sweep_manifest_v1"),
                        "validate_sweep_manifest_v1 is missing")

    def test_sweep_schema_golden_and_invalid_fixture(self):
        validator = getattr(schemas, "validate_sweep_manifest_v1", None)
        if validator is None: self.skipTest("validate_sweep_manifest_v1 is missing")
        report = sweep_manifest()
        self.assertIsNone(validator(report))
        invalid = copy.deepcopy(report)
        invalid["raw_states"] = [{"potential": [0.0] * 6}]
        with self.assertRaises(ValueError): validator(invalid)


if __name__ == "__main__":
    unittest.main()
