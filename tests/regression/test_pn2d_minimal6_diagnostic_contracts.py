import copy
import csv
import hashlib
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
CONTACT_TOLERANCE_FORMULA = "max(1e-18 A/um, 1e-9 * max(|Anode|, |Cathode|))"


def canonical_sha(value):
    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def schema_root():
    for parent in Path(__file__).resolve().parents:
        candidate = parent / "schemas" / "vela.pn2d_minimal6_bv_comparison.v1.schema.json"
        if candidate.is_file():
            return parent / "schemas"
    raise FileNotFoundError("repository schemas directory is unavailable")


def solver_configuration(solver):
    primary, secondary = (SHA_A, SHA_B) if solver == "vela" else (SHA_B, SHA_A)
    config = {
        "template": {"path": f"fixture/{solver}/template.json", "sha256": primary},
        "topology_input_sha256": {
            "sketch": {"mesh.json": primary},
            "mirror": {"mesh.json": secondary},
        },
        "deck_sha256": [primary],
    }
    return {**config, "configuration_sha256": canonical_sha(config)}


def contact_conservation(row):
    anode = row["observables"]["anode_current_A_per_um"]
    cathode = row["observables"]["cathode_current_A_per_um"]
    residual = anode + cathode
    tolerance = max(1.0e-18, 1.0e-9 * max(abs(anode), abs(cathode)))
    return {
        "anode_current_A_per_um": anode,
        "cathode_current_A_per_um": cathode,
        "signed_residual_A_per_um": residual,
        "tolerance_A_per_um": tolerance,
        "tolerance_formula": CONTACT_TOLERANCE_FORMULA,
        "classification": "conserved" if abs(residual) <= tolerance else "not_conserved",
    }


def closed_gap(quantity):
    return {
        "quantity": quantity,
        "classification": "available",
        "log_gap_dex": 0.0,
        "named_contributions": [
            {"name": f"{quantity}_difference", "contribution_dex": 0.0}
        ],
        "residual": {
            "name": "cross_solver_semantics_residual",
            "classification": "available",
            "value_dex": 0.0,
        },
        "closure_error_dex": 0.0,
    }


def validate_schema_document(instance, schema, *, root=None, path="$"):
    """Validate the Draft 2020-12 subset used by the three tracked schemas."""
    root = schema if root is None else root
    if "$ref" in schema:
        pointer = schema["$ref"]
        if not isinstance(pointer, str) or not pointer.startswith("#/"):
            raise ValueError(f"{path}: only local JSON pointers are supported")
        target = root
        for token in pointer[2:].split("/"):
            token = token.replace("~1", "/").replace("~0", "~")
            target = target[token]
        validate_schema_document(instance, target, root=root, path=path)

    def accepts(candidate):
        try:
            validate_schema_document(instance, candidate, root=root, path=path)
            return True
        except ValueError:
            return False

    if "allOf" in schema and not all(accepts(item) for item in schema["allOf"]):
        raise ValueError(f"{path}: allOf failed")
    if "anyOf" in schema and not any(accepts(item) for item in schema["anyOf"]):
        raise ValueError(f"{path}: anyOf failed")
    if "oneOf" in schema and sum(accepts(item) for item in schema["oneOf"]) != 1:
        raise ValueError(f"{path}: oneOf failed")
    if "not" in schema and accepts(schema["not"]):
        raise ValueError(f"{path}: forbidden by not")

    expected = schema.get("type")
    expected_types = expected if isinstance(expected, list) else [expected]
    type_matches = {
        "object": isinstance(instance, dict),
        "array": isinstance(instance, list),
        "string": isinstance(instance, str),
        "number": isinstance(instance, (int, float)) and not isinstance(instance, bool),
        "integer": isinstance(instance, int) and not isinstance(instance, bool),
        "boolean": isinstance(instance, bool),
        "null": instance is None,
    }
    if expected is not None and not any(type_matches.get(item, False) for item in expected_types):
        raise ValueError(f"{path}: expected {expected}")
    if isinstance(instance, float) and not math.isfinite(instance):
        raise ValueError(f"{path}: number must be finite")
    if isinstance(instance, (int, float)) and not isinstance(instance, bool):
        if "minimum" in schema and instance < schema["minimum"]:
            raise ValueError(f"{path}: below minimum")
        if "maximum" in schema and instance > schema["maximum"]:
            raise ValueError(f"{path}: above maximum")
    if "const" in schema and instance != schema["const"]:
        raise ValueError(f"{path}: const mismatch")
    if "enum" in schema and instance not in schema["enum"]:
        raise ValueError(f"{path}: enum mismatch")
    if isinstance(instance, str):
        if len(instance) < int(schema.get("minLength", 0)):
            raise ValueError(f"{path}: string too short")
        if "pattern" in schema:
            import re
            if re.search(schema["pattern"], instance) is None:
                raise ValueError(f"{path}: pattern mismatch")
    if isinstance(instance, list):
        if len(instance) < int(schema.get("minItems", 0)):
            raise ValueError(f"{path}: too few items")
        if "maxItems" in schema and len(instance) > int(schema["maxItems"]):
            raise ValueError(f"{path}: too many items")
        if schema.get("uniqueItems"):
            encoded = [json.dumps(item, sort_keys=True, allow_nan=False) for item in instance]
            if len(encoded) != len(set(encoded)):
                raise ValueError(f"{path}: duplicate items")
        prefixes = schema.get("prefixItems", [])
        if not isinstance(prefixes, list):
            raise ValueError(f"{path}: prefixItems must be an array")
        for index, item in enumerate(instance):
            if index < len(prefixes):
                validate_schema_document(item, prefixes[index], root=root, path=f"{path}[{index}]")
            else:
                item_schema = schema.get("items")
                if item_schema is False:
                    raise ValueError(f"{path}: additional array item")
                if isinstance(item_schema, dict):
                    validate_schema_document(item, item_schema, root=root, path=f"{path}[{index}]")
        if not prefixes and isinstance(schema.get("items"), dict):
            for index, item in enumerate(instance):
                validate_schema_document(item, schema["items"], root=root, path=f"{path}[{index}]")
    if isinstance(instance, dict):
        if len(instance) < int(schema.get("minProperties", 0)):
            raise ValueError(f"{path}: too few properties")
        missing = [name for name in schema.get("required", []) if name not in instance]
        if missing:
            raise ValueError(f"{path}: missing {missing}")
        properties = schema.get("properties", {})
        for name, value in instance.items():
            child = properties.get(name)
            if isinstance(child, dict):
                validate_schema_document(value, child, root=root, path=f"{path}.{name}")
            elif schema.get("additionalProperties") is False:
                raise ValueError(f"{path}: unexpected property {name}")
            elif isinstance(schema.get("additionalProperties"), dict):
                validate_schema_document(
                    value, schema["additionalProperties"], root=root, path=f"{path}.{name}"
                )


def schema_document(name):
    path = schema_root() / f"{name}.schema.json"
    return json.loads(path.read_text(encoding="utf-8"))


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
                          "name": "sentaurus_internal_semantics_residual",
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


def transition(solver, topology="sketch", bias=-1.0, branch_classification="multiplication_like"):
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
        "branch_classification": branch_classification,
        "branch_threshold_version": THRESHOLD,
        "convergence_metadata": {"iterations": 4, "residual_norm": 1.0e-10},
        "stdout": "", "stderr": "",
    }


def comparison_report():
    vela = transition("vela", branch_classification="unidentified")
    sentaurus = transition("sentaurus", branch_classification="unidentified")
    checkpoint = {
        "topology": "sketch", "bias_V": -1.0, "classification": "common_exact",
        "vela": copy.deepcopy(vela), "sentaurus": copy.deepcopy(sentaurus),
        "branch_classification": "multiplication_like",
        "branch_threshold_version": THRESHOLD,
        "branch_ratio_evidence": {
            "vela_anode_current_A_per_um": -1.0e-8,
            "sentaurus_anode_current_A_per_um": -1.0e-8,
            "absolute_vela_over_sentaurus": 1.0,
            "geometric_zero": False,
            "threshold_version": THRESHOLD,
        },
        "terminal_current_ratio": {"classification": "available", "value": 1.0},
        "maximum_field_ratio": {"classification": "available", "value": 1.0},
        "native_source_ratio": {"classification": "available", "value": 1.0},
        "reconstructed_source_ratio": {"classification": "available", "value": 1.0},
        "vela_one_volt_current_growth": {
            "classification": "unavailable", "value": None,
            "reason": "next exact one-volt checkpoint unavailable",
        },
        "sentaurus_one_volt_current_growth": {
            "classification": "unavailable", "value": None,
            "reason": "next exact one-volt checkpoint unavailable",
        },
        "gap_closure": {
            "status": "closed",
            "tolerance_dex": 1.0e-10,
            "gaps": [
                closed_gap(quantity)
                for quantity in (
                    "terminal_current", "maximum_field",
                    "native_source", "reconstructed_source",
                )
            ],
        },
        "contact_current_conservation": {
            "unit": "A/um",
            "vela": contact_conservation(vela),
            "sentaurus": contact_conservation(sentaurus),
        },
    }
    metadata = {"DiagnosticDisclaimer": DISCLAIMER, "SolverTermination": "Every recorded solver failure transition is explicitly marked.", "BVExtrapolation": "No physical breakdown voltage (BV) is extrapolated."}
    names = ("terminal_current.png", "one_volt_growth.png", "maximum_field.png", "source_integrals.png", "topology.png")
    series = {
        "terminal_current.png": [{"solver": solver, "topology": "sketch", "quantity": "terminal_current"} for solver in ("vela", "sentaurus")],
        "maximum_field.png": [{"solver": solver, "topology": "sketch", "quantity": "maximum_field"} for solver in ("vela", "sentaurus")],
        "source_integrals.png": [{"solver": solver, "topology": "sketch", "quantity": quantity} for solver in ("vela", "sentaurus") for quantity in ("native_source", "reconstructed_source")],
        "one_volt_growth.png": [], "topology.png": [],
    }
    figures = {name: {"sha256": SHA_A, "width_px": 900, "height_px": 504, "metadata": metadata, "series_identities": series[name], "failure_transition_markers": []} for name in names}
    return {
        "schema": "vela.pn2d_minimal6_bv_comparison.v1",
        "diagnostic_disclaimer": DISCLAIMER, "interpolation": "forbidden",
        "branch_threshold_version": THRESHOLD,
        "solver_configurations": {
            "vela": solver_configuration("vela"),
            "sentaurus": solver_configuration("sentaurus"),
        },
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
        "artifact_hashes": {name: SHA_A for name in names}, "input_artifacts": {},
        "figure_contract": {"schema": "vela.pn2d_minimal6_figure_contract.v1", "figures": figures}, "closure": {
            "status": "closed",
            "eligible_gaps": 4,
            "rule": "each eligible gap records named contributions and a typed residual; non-common points are side-only",
        },
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
        root = schema_root()
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

    def test_formula_actual_bias_tolerance_matches_runtime_and_schema(self):
        document = schema_document("vela.pn2d_minimal6_formula_difference.v1")
        for index, state in enumerate(state_matrix()):
            expected = state["requested_bias_V"]
            for sign in (-1.0, 1.0):
                for deviation in (1.0e-13, 1.0e-12):
                    report = formula_report()
                    report["state_matrix"][index]["actual_bias_V"] = expected + sign * deviation
                    with self.subTest(
                        index=index, sign=sign, deviation=deviation, validator="runtime"
                    ):
                        self.assertIsNone(schemas.validate_formula_difference_v1(report))
                    with self.subTest(
                        index=index, sign=sign, deviation=deviation, validator="schema"
                    ):
                        validate_schema_document(report, document)

            for sign in (-1.0, 1.0):
                outside = formula_report()
                outside["state_matrix"][index]["actual_bias_V"] = expected + sign * 1.01e-12
                with self.subTest(index=index, sign=sign, validator="runtime-outside"):
                    with self.assertRaises(ValueError):
                        schemas.validate_formula_difference_v1(outside)
                with self.subTest(index=index, sign=sign, validator="schema-outside"):
                    with self.assertRaises(ValueError):
                        validate_schema_document(outside, document)

    def test_transition_bias_tolerance_uses_inclusive_nonzero_bounds(self):
        boundaries = {
            -12.0: (-12.000000000001, -11.999999999999),
            -19.0: (-19.000000000001, -18.999999999999),
        }

        def validate(kind, target, actual):
            report = comparison_report() if kind == "comparison" else sweep_manifest()
            if kind == "comparison":
                row = report["accepted_transitions"]["vela"][0]
                validator = schemas.validate_bv_comparison_v1
            else:
                row = report["accepted_checkpoints"][0]
                validator = schemas.validate_sweep_manifest_v1
            row["target_bias_V"] = target
            row["actual_bias_V"] = actual
            return validator(report)

        for target, bounds in boundaries.items():
            accepted = (*bounds, target - 1.0e-13, target + 1.0e-13)
            for actual in accepted:
                for kind in ("comparison", "sweep"):
                    with self.subTest(kind=kind, target=target, actual=actual):
                        self.assertIsNone(validate(kind, target, actual))
            for sign in (-1.0, 1.0):
                actual = target + sign * 1.01e-12
                for kind in ("comparison", "sweep"):
                    with self.subTest(kind=kind, target=target, actual=actual):
                        with self.assertRaisesRegex(ValueError, "not an exact checkpoint"):
                            validate(kind, target, actual)

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

    def test_formula_runtime_requires_paths_interactions_and_named_residual_identities(self):
        interaction = {
            "first_factor": "mobility", "second_factor": "alpha_law",
            "path_identity": "forward_adjacent", "baseline": 1.0,
            "a_only": 2.0, "b_only": 3.0, "both": 6.0, "interaction_dex": 0.0,
        }
        golden = formula_report()
        golden["waterfall_paths"][0]["interactions"] = [copy.deepcopy(interaction)]
        golden["interactions"] = [{
            "topology": "sketch", "bias_V": 0.0, **copy.deepcopy(interaction)
        }]
        self.assertIsNone(schemas.validate_formula_difference_v1(golden))

        mutations = []
        empty_forward = copy.deepcopy(golden)
        empty_forward["waterfall_paths"][0]["forward"] = {}
        mutations.append(empty_forward)
        malformed_interaction = copy.deepcopy(golden)
        del malformed_interaction["interactions"][0]["first_factor"]
        mutations.append(malformed_interaction)
        wrong_interaction_identity = copy.deepcopy(golden)
        wrong_interaction_identity["interactions"][0]["bias_V"] = -12.0
        wrong_interaction_identity["interactions"].append(
            copy.deepcopy(wrong_interaction_identity["interactions"][0])
        )
        mutations.append(wrong_interaction_identity)
        unnamed_residual = copy.deepcopy(golden)
        del unnamed_residual["sentaurus_internal_semantics_residual"][0]["name"]
        mutations.append(unnamed_residual)
        wrong_residual_name = copy.deepcopy(golden)
        wrong_residual_name["sentaurus_internal_semantics_residual"][0]["name"] = "residual"
        mutations.append(wrong_residual_name)
        duplicate_residual = copy.deepcopy(golden)
        duplicate_residual["sentaurus_internal_semantics_residual"][-1] = copy.deepcopy(
            duplicate_residual["sentaurus_internal_semantics_residual"][0]
        )
        mutations.append(duplicate_residual)
        for index, mutation in enumerate(mutations):
            with self.subTest(index=index), self.assertRaises(ValueError):
                schemas.validate_formula_difference_v1(mutation)

    def test_comparison_runtime_requires_transition_branch_threshold_and_convergence(self):
        failure = {
            "solver": "vela", "topology": "mirror", "start_bias_V": 0.0,
            "target_bias_V": -1.0, "status": "rejected", "observables": None,
            "branch_classification": "unidentified", "branch_threshold_version": THRESHOLD,
            "convergence_metadata": {"exit_code": 1, "reason": "Newton failure"},
        }
        golden = comparison_report()
        golden["failed_transitions"] = [failure]
        golden["failure_transitions"] = [copy.deepcopy(failure)]
        marker = {key: failure[key] for key in ("solver", "topology", "start_bias_V", "target_bias_V")}
        for figure in golden["figure_contract"]["figures"].values():
            figure["failure_transition_markers"] = [copy.deepcopy(marker)]
        self.assertIsNone(schemas.validate_bv_comparison_v1(golden))
        for collection, index, field in (
            ("accepted_transitions", 0, "branch_classification"),
            ("accepted_transitions", 0, "branch_threshold_version"),
            ("accepted_transitions", 0, "convergence_metadata"),
            ("failed_transitions", 0, "branch_classification"),
            ("failed_transitions", 0, "branch_threshold_version"),
            ("failed_transitions", 0, "convergence_metadata"),
        ):
            invalid = copy.deepcopy(golden)
            rows = invalid[collection]["vela"] if collection == "accepted_transitions" else invalid[collection]
            del rows[index][field]
            with self.subTest(collection=collection, field=field), self.assertRaises(ValueError):
                schemas.validate_bv_comparison_v1(invalid)

    def test_comparison_runtime_ties_every_transition_version_to_report(self):
        failure = {
            "solver": "vela", "topology": "mirror", "start_bias_V": 0.0,
            "target_bias_V": -1.0, "status": "rejected", "observables": None,
            "branch_classification": "unidentified", "branch_threshold_version": THRESHOLD,
            "convergence_metadata": {"exit_code": 1, "reason": "Newton failure"},
        }
        golden = comparison_report()
        golden["failed_transitions"] = [failure]
        golden["failure_transitions"] = [copy.deepcopy(failure)]
        mutations = []
        for solver in ("vela", "sentaurus"):
            invalid = copy.deepcopy(golden)
            invalid["accepted_transitions"][solver][0]["branch_threshold_version"] = "contradictory"
            mutations.append((f"accepted-{solver}", invalid))
        invalid = copy.deepcopy(golden)
        invalid["failed_transitions"][0]["branch_threshold_version"] = "contradictory"
        invalid["failure_transitions"][0]["branch_threshold_version"] = "contradictory"
        mutations.append(("failed", invalid))
        invalid = copy.deepcopy(golden)
        invalid["checkpoints"][0]["branch_threshold_version"] = "contradictory"
        mutations.append(("checkpoint", invalid))
        for solver in ("vela", "sentaurus"):
            invalid = copy.deepcopy(golden)
            invalid["checkpoints"][0][solver]["branch_threshold_version"] = "contradictory"
            mutations.append((f"checkpoint-{solver}", invalid))
        for location, invalid in mutations:
            with self.subTest(location=location), self.assertRaisesRegex(
                ValueError, "branch threshold version mismatch"
            ):
                schemas.validate_bv_comparison_v1(invalid)

    def test_sweep_runtime_ties_every_transition_version_to_package(self):
        failure = {
            "solver": "vela", "topology": "mirror", "start_bias_V": 0.0,
            "target_bias_V": -1.0, "status": "rejected", "observables": None,
            "branch_classification": "unidentified", "branch_threshold_version": THRESHOLD,
            "convergence_metadata": {"exit_code": 1, "reason": "Newton failure"},
        }
        golden = sweep_manifest()
        golden["failed_transition"] = copy.deepcopy(failure)
        golden["failed_transitions"] = [copy.deepcopy(failure)]
        mutations = []
        invalid = copy.deepcopy(golden)
        invalid["accepted_checkpoints"][0]["branch_threshold_version"] = "contradictory"
        mutations.append(("accepted", invalid))
        invalid = copy.deepcopy(golden)
        invalid["failed_transitions"][0]["branch_threshold_version"] = "contradictory"
        mutations.append(("failed-list", invalid))
        invalid = copy.deepcopy(golden)
        invalid["failed_transition"]["branch_threshold_version"] = "contradictory"
        mutations.append(("failed-singular", invalid))
        for location, invalid in mutations:
            with self.subTest(location=location), self.assertRaisesRegex(
                ValueError, "branch threshold version mismatch"
            ):
                schemas.validate_sweep_manifest_v1(invalid)

    def test_figure_schema_requires_semantic_metadata_series_and_markers(self):
        comparison = comparison_report()
        schema = schema_document("vela.pn2d_minimal6_bv_comparison.v1")
        mutations = (
            lambda value: value["figure_contract"]["figures"]["terminal_current.png"]["metadata"].update(
                DiagnosticDisclaimer="not the diagnostic disclaimer"
            ),
            lambda value: value["figure_contract"]["figures"]["terminal_current.png"].update(
                series_identities=[{"solver": "vela"}]
            ),
            lambda value: value["figure_contract"]["figures"]["terminal_current.png"].update(
                failure_transition_markers=[{"solver": "vela"}]
            ),
        )
        for mutate in mutations:
            with self.subTest(mutate=mutate):
                invalid = copy.deepcopy(comparison)
                mutate(invalid)
                with self.assertRaises(ValueError):
                    validate_schema_document(invalid, schema)
    def test_schema_documents_execute_against_golden_and_invalid_fixtures(self):
        formula = formula_report()
        formula_schema = schema_document("vela.pn2d_minimal6_formula_difference.v1")
        comparison = comparison_report()
        comparison_schema = schema_document("vela.pn2d_minimal6_bv_comparison.v1")
        sweep = sweep_manifest()
        sweep_schema = schema_document("vela.pn2d_minimal6_sweep_manifest.v1")
        for instance, document in (
            (formula, formula_schema), (comparison, comparison_schema), (sweep, sweep_schema)
        ):
            validate_schema_document(instance, document)

        formula_mutations = []
        wrong_order = copy.deepcopy(formula)
        wrong_order["state_matrix"][0], wrong_order["state_matrix"][1] = (
            wrong_order["state_matrix"][1], wrong_order["state_matrix"][0]
        )
        formula_mutations.append(wrong_order)
        extra_state_property = copy.deepcopy(formula)
        extra_state_property["state_matrix"][0]["label"] = "free-form"
        formula_mutations.append(extra_state_property)
        wrong_actual_bias = copy.deepcopy(formula)
        wrong_actual_bias["state_matrix"][0]["actual_bias_V"] = 1.000001e-12
        formula_mutations.append(wrong_actual_bias)
        wrong_residual_name = copy.deepcopy(formula)
        wrong_residual_name["sentaurus_internal_semantics_residual"][0]["name"] = "residual"
        formula_mutations.append(wrong_residual_name)
        malformed_interaction = copy.deepcopy(formula)
        malformed_interaction["waterfall_paths"][0]["interactions"] = [{"first_factor": "mobility"}]
        formula_mutations.append(malformed_interaction)
        for index, invalid in enumerate(formula_mutations):
            with self.subTest(schema="formula", index=index), self.assertRaises(ValueError):
                validate_schema_document(invalid, formula_schema)

        for location in ("accepted", "checkpoint", "solver_configuration"):
            invalid = copy.deepcopy(comparison)
            if location == "accepted":
                invalid["accepted_transitions"]["vela"][0]["raw_state"] = {"psi": [0.0]}
            elif location == "checkpoint":
                invalid["checkpoints"][0]["raw_states"] = []
            else:
                invalid["solver_configurations"]["vela"]["state_payload"] = {}
            with self.subTest(schema="comparison", location=location), self.assertRaises(ValueError):
                validate_schema_document(invalid, comparison_schema)

        for field in ("branch_classification", "branch_threshold_version", "convergence_metadata"):
            invalid = copy.deepcopy(comparison)
            del invalid["accepted_transitions"]["vela"][0][field]
            with self.subTest(schema="comparison", field=field), self.assertRaises(ValueError):
                validate_schema_document(invalid, comparison_schema)
        failed = {
            "solver": "vela", "topology": "mirror", "start_bias_V": 0.0,
            "target_bias_V": -1.0, "status": "rejected", "observables": None,
            "state_path": None, "state_sha256": None,
            "branch_classification": "unidentified", "branch_threshold_version": THRESHOLD,
            "convergence_metadata": {"exit_code": 1},
        }
        comparison_with_failure = copy.deepcopy(comparison)
        comparison_with_failure["failed_transitions"] = [failed]
        comparison_with_failure["failure_transitions"] = [copy.deepcopy(failed)]
        validate_schema_document(comparison_with_failure, comparison_schema)
        for field in ("branch_classification", "branch_threshold_version", "convergence_metadata"):
            invalid = copy.deepcopy(comparison_with_failure)
            del invalid["failed_transitions"][0][field]
            with self.subTest(schema="comparison-failure", field=field), self.assertRaises(ValueError):
                validate_schema_document(invalid, comparison_schema)

        for location in ("checkpoint", "segment"):
            invalid = copy.deepcopy(sweep)
            if location == "checkpoint":
                invalid["accepted_checkpoints"][0]["raw_state"] = {"psi": [0.0]}
            else:
                invalid["segments"] = [{"solver": "vela", "raw_states": []}]
            with self.subTest(schema="sweep", location=location), self.assertRaises(ValueError):
                validate_schema_document(invalid, sweep_schema)

        sweep_with_failure = copy.deepcopy(sweep)
        sweep_with_failure["failed_transition"] = failed
        sweep_with_failure["failed_transitions"] = [copy.deepcopy(failed)]
        validate_schema_document(sweep_with_failure, sweep_schema)
        for field in ("branch_classification", "branch_threshold_version", "convergence_metadata"):
            invalid = copy.deepcopy(sweep_with_failure)
            del invalid["failed_transitions"][0][field]
            with self.subTest(schema="sweep-failure", field=field), self.assertRaises(ValueError):
                validate_schema_document(invalid, sweep_schema)

        nonfinite = copy.deepcopy(comparison)
        nonfinite["accepted_transitions"]["vela"][0]["observables"]["max_field_V_per_m"] = math.nan
        with self.assertRaises(ValueError):
            validate_schema_document(nonfinite, comparison_schema)


if __name__ == "__main__":
    unittest.main()
