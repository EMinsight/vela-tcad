"""Strict runtime validation for Minimal6 diagnostic report contracts."""

from collections.abc import Mapping
import hashlib
import json
import math
import re


DISCLAIMER = "minimal6 diagnostic sweep; not a physical BV curve"
_BRANCHES = {"leakage_like", "multiplication_like", "unidentified"}
_OBSERVABLES = {
    "anode_current_A_per_um",
    "cathode_current_A_per_um",
    "max_field_V_per_m",
    "native_source_integral_s_inv_per_cm",
    "reconstructed_source_integral_s_inv_per_cm",
}
_EXPECTED_STATES = {
    (topology, bias)
    for topology in ("sketch", "mirror")
    for bias in (0.0, -12.0, -19.0)
}
_SHA256 = re.compile(r"[0-9a-fA-F]{64}\Z")
_RAW_STATE_KEYS = {"raw_state", "raw_states", "state_payload", "state_vectors", "state_values"}
_BIAS_TOLERANCE_V = 1.0e-12
_CONTACT_TOLERANCE_RELATIVE = 1.0e-9
_CONTACT_TOLERANCE_FLOOR_A_PER_UM = 1.0e-18
_BRANCH_THRESHOLD_VERSION = "v1: multiplication=[0.1,10], leakage<=1e-3"
_GEOMETRIC_ZERO_SOURCE_INTEGRAL = 1.0e-285
_GAP_TOLERANCE_DEX = 1.0e-10
_GAP_QUANTITIES = {
    "terminal_current": ("terminal_current_ratio", "anode_current_A_per_um", True),
    "maximum_field": ("maximum_field_ratio", "max_field_V_per_m", False),
    "native_source": ("native_source_ratio", "native_source_integral_s_inv_per_cm", False),
    "reconstructed_source": ("reconstructed_source_ratio", "reconstructed_source_integral_s_inv_per_cm", False),
}


def _finite_tree(value, path="report") -> None:
    if isinstance(value, float) and not math.isfinite(value):
        raise ValueError(f"non-finite values are forbidden at {path}")
    if isinstance(value, Mapping):
        for key, item in value.items():
            _finite_tree(item, f"{path}.{key}")
    elif isinstance(value, list):
        for index, item in enumerate(value):
            _finite_tree(item, f"{path}[{index}]")


def _require_mapping(value, name: str, *, nonempty: bool = False) -> Mapping:
    if not isinstance(value, Mapping) or (nonempty and not value):
        qualifier = "non-empty " if nonempty else ""
        raise ValueError(f"{name} must be a {qualifier}object")
    return value


def _require_list(value, name: str) -> list:
    if not isinstance(value, list):
        raise ValueError(f"{name} must be an array")
    return value


def _require_fields(report: Mapping, schema: str, fields: set[str]) -> None:
    if report.get("schema") != schema:
        raise ValueError("invalid report schema")
    if report.get("diagnostic_disclaimer") != DISCLAIMER:
        raise ValueError("missing diagnostic disclaimer")
    missing = sorted(fields - set(report))
    if missing:
        raise ValueError(f"missing required report fields: {missing}")
    _finite_tree(report)


def _number(value, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{name} must be numeric")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{name} must be finite")
    return result

def _within_bias_tolerance(expected: float, actual: float) -> bool:
    return expected - _BIAS_TOLERANCE_V <= actual <= expected + _BIAS_TOLERANCE_V


def _sha(value, name: str) -> None:
    if not isinstance(value, str) or not _SHA256.fullmatch(value):
        raise ValueError(f"{name} must be a SHA-256 hex digest")


def _state_identity(row: Mapping, name: str) -> tuple[str, float]:
    topology = row.get("topology_id", row.get("topology"))
    bias = row.get("requested_bias_V", row.get("bias_V"))
    if topology not in {"sketch", "mirror"}:
        raise ValueError(f"{name} has invalid topology")
    bias_value = _number(bias, f"{name}.bias_V")
    if "actual_bias_V" in row:
        actual = _number(row["actual_bias_V"], f"{name}.actual_bias_V")
        if not _within_bias_tolerance(bias_value, actual):
            raise ValueError(f"{name} is not an exact state")
    return str(topology), bias_value


def _reject_raw_states(value, path="report") -> None:
    if isinstance(value, Mapping):
        forbidden = _RAW_STATE_KEYS & set(value)
        if forbidden:
            raise ValueError(f"embedded raw states are forbidden at {path}: {sorted(forbidden)}")
        for key, item in value.items():
            _reject_raw_states(item, f"{path}.{key}")
    elif isinstance(value, list):
        for index, item in enumerate(value):
            _reject_raw_states(item, f"{path}[{index}]")


def _observables(value, name: str) -> Mapping:
    values = _require_mapping(value, name)
    if set(values) != _OBSERVABLES:
        missing = sorted(_OBSERVABLES - set(values))
        extra = sorted(set(values) - _OBSERVABLES)
        raise ValueError(f"{name} has invalid named observables; missing={missing}, extra={extra}")
    for key, item in values.items():
        _number(item, f"{name}.{key}")
    return values


def _branch_evidence(row: Mapping, name: str) -> str:
    classification = row.get("branch_classification")
    if classification not in _BRANCHES:
        raise ValueError(f"{name} lacks typed branch classification")
    version = row.get("branch_threshold_version")
    if not isinstance(version, str) or not version:
        raise ValueError(f"{name} lacks branch-threshold provenance")
    _require_mapping(row.get("convergence_metadata"), f"{name}.convergence_metadata", nonempty=True)
    return version


def _canonical_sha(value) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _solver_configuration(value, name: str) -> None:
    config = _require_mapping(value, name)
    if set(config) != {"template", "topology_input_sha256", "deck_sha256", "configuration_sha256"}:
        raise ValueError(f"{name} has invalid configuration provenance fields")
    template = _require_mapping(config["template"], f"{name}.template", nonempty=True)
    _sha(template.get("sha256"), f"{name}.template.sha256")
    topology_hashes = _require_mapping(config["topology_input_sha256"], f"{name}.topology_input_sha256")
    if set(topology_hashes) != {"sketch", "mirror"}:
        raise ValueError(f"{name} requires sketch/mirror topology hashes")
    for topology, entries in topology_hashes.items():
        entries = _require_mapping(entries, f"{name}.{topology}", nonempty=True)
        for key, digest in entries.items():
            if not isinstance(key, str) or not key:
                raise ValueError(f"{name}.{topology} has an invalid input name")
            _sha(digest, f"{name}.{topology}.{key}")
    for index, digest in enumerate(_require_list(config["deck_sha256"], f"{name}.deck_sha256")):
        _sha(digest, f"{name}.deck_sha256[{index}]")
    supplied = config["configuration_sha256"]
    _sha(supplied, f"{name}.configuration_sha256")
    canonical = {key: config[key] for key in ("template", "topology_input_sha256", "deck_sha256")}
    if supplied != _canonical_sha(canonical):
        raise ValueError(f"{name}.configuration_sha256 does not match configuration content")


def _contact_conservation(value, name: str) -> None:
    evidence = _require_mapping(value, name)
    expected_fields = {
        "anode_current_A_per_um", "cathode_current_A_per_um",
        "signed_residual_A_per_um", "tolerance_A_per_um", "tolerance_formula", "classification",
    }
    if set(evidence) != expected_fields:
        raise ValueError(f"{name} has invalid conservation fields")
    anode = _number(evidence["anode_current_A_per_um"], f"{name}.anode_current_A_per_um")
    cathode = _number(evidence["cathode_current_A_per_um"], f"{name}.cathode_current_A_per_um")
    residual = _number(evidence["signed_residual_A_per_um"], f"{name}.signed_residual_A_per_um")
    tolerance = _number(evidence["tolerance_A_per_um"], f"{name}.tolerance_A_per_um")
    expected_residual = anode + cathode
    expected_tolerance = max(
        _CONTACT_TOLERANCE_FLOOR_A_PER_UM,
        max(abs(anode), abs(cathode)) * _CONTACT_TOLERANCE_RELATIVE,
    )
    if residual != expected_residual or tolerance != expected_tolerance:
        raise ValueError(f"{name} conservation residual/tolerance was not independently recomputed")
    if evidence["tolerance_formula"] != "max(1e-18 A/um, 1e-9 * max(|Anode|, |Cathode|))":
        raise ValueError(f"{name} has a noncanonical tolerance formula")
    expected_classification = "conserved" if abs(residual) <= tolerance else "not_conserved"
    if evidence["classification"] != expected_classification:
        raise ValueError(f"{name} conservation classification is inconsistent")


def _geometric_zero_pair(*rows: Mapping) -> bool:
    return any(
        abs(_number(row["observables"]["native_source_integral_s_inv_per_cm"], "native_source"))
        <= _GEOMETRIC_ZERO_SOURCE_INTEGRAL
        and abs(_number(row["observables"]["reconstructed_source_integral_s_inv_per_cm"], "reconstructed_source"))
        <= _GEOMETRIC_ZERO_SOURCE_INTEGRAL
        for row in rows
    )


def _expected_ratio(numerator: float, denominator: float, *, geometric_zero: bool, absolute: bool) -> dict:
    if geometric_zero:
        return {"classification": "geometric_zero", "value": None}
    left, right = (abs(numerator), abs(denominator)) if absolute else (numerator, denominator)
    if right == 0.0:
        return {"classification": "zero_denominator", "value": None}
    if left == 0.0:
        return {"classification": "zero_numerator", "value": 0.0}
    return {"classification": "available", "value": left / right}


def _validate_ratio_record(value, expected: dict, name: str) -> None:
    ratio = _require_mapping(value, name)
    if set(ratio) != {"classification", "value"}:
        raise ValueError(f"{name} has invalid ratio fields")
    if ratio.get("classification") != expected["classification"]:
        raise ValueError(f"{name} classification does not match nested observables")
    actual_value = ratio.get("value")
    expected_value = expected["value"]
    if expected_value is None:
        if actual_value is not None:
            raise ValueError(f"{name} unavailable ratio must have null value")
    elif _number(actual_value, f"{name}.value") != expected_value:
        raise ValueError(f"{name} value does not match nested observables")


def _validate_checkpoint_branch_and_gaps(row: Mapping, name: str, version: str) -> int:
    vela = _require_mapping(row["vela"], f"{name}.vela")
    sentaurus = _require_mapping(row["sentaurus"], f"{name}.sentaurus")
    vela_observables = _require_mapping(vela["observables"], f"{name}.vela.observables")
    sentaurus_observables = _require_mapping(sentaurus["observables"], f"{name}.sentaurus.observables")
    vela_current = _number(vela_observables["anode_current_A_per_um"], f"{name}.vela.anode_current")
    sentaurus_current = _number(sentaurus_observables["anode_current_A_per_um"], f"{name}.sentaurus.anode_current")
    geometric_zero = _geometric_zero_pair(vela, sentaurus)
    absolute_ratio = None if sentaurus_current == 0.0 else abs(vela_current / sentaurus_current)
    expected_evidence = {
        "vela_anode_current_A_per_um": vela_current,
        "sentaurus_anode_current_A_per_um": sentaurus_current,
        "absolute_vela_over_sentaurus": absolute_ratio,
        "geometric_zero": geometric_zero,
        "threshold_version": version,
    }
    evidence = _require_mapping(row.get("branch_ratio_evidence"), f"{name}.branch_ratio_evidence")
    if set(evidence) != set(expected_evidence) or evidence != expected_evidence:
        raise ValueError(f"{name} branch ratio evidence does not match nested observables")
    if geometric_zero or sentaurus_current == 0.0 or vela_current == 0.0:
        expected_branch = "unidentified"
    elif 0.1 <= absolute_ratio <= 10.0:
        expected_branch = "multiplication_like"
    elif absolute_ratio <= 1.0e-3:
        expected_branch = "leakage_like"
    else:
        expected_branch = "unidentified"
    if row.get("branch_classification") != expected_branch:
        raise ValueError(f"{name} branch classification does not match recomputed evidence")

    expected_ratios = {}
    for quantity, (field, observable, absolute) in _GAP_QUANTITIES.items():
        expected = _expected_ratio(
            _number(vela_observables[observable], f"{name}.vela.{observable}"),
            _number(sentaurus_observables[observable], f"{name}.sentaurus.{observable}"),
            geometric_zero=geometric_zero,
            absolute=absolute,
        )
        _validate_ratio_record(row.get(field), expected, f"{name}.{field}")
        expected_ratios[quantity] = expected

    closure = _require_mapping(row.get("gap_closure"), f"{name}.gap_closure")
    if set(closure) != {"status", "tolerance_dex", "gaps"} or closure.get("status") != "closed":
        raise ValueError(f"{name}.gap_closure has invalid closure fields")
    if _number(closure.get("tolerance_dex"), f"{name}.gap_closure.tolerance_dex") != _GAP_TOLERANCE_DEX:
        raise ValueError(f"{name}.gap_closure has noncanonical tolerance")
    gaps = _require_list(closure.get("gaps"), f"{name}.gap_closure.gaps")
    by_quantity = {}
    for index, gap_value in enumerate(gaps):
        gap = _require_mapping(gap_value, f"{name}.gap_closure.gaps[{index}]")
        quantity = gap.get("quantity")
        if quantity not in _GAP_QUANTITIES or quantity in by_quantity:
            raise ValueError(f"{name}.gap_closure has missing, duplicate, or unknown quantities")
        by_quantity[quantity] = gap
    if set(by_quantity) != set(_GAP_QUANTITIES) or len(gaps) != len(_GAP_QUANTITIES):
        raise ValueError(f"{name}.gap_closure must contain exactly four named quantities")

    eligible = 0
    for quantity, expected_ratio in expected_ratios.items():
        gap = by_quantity[quantity]
        expected_fields = {
            "quantity", "classification", "log_gap_dex", "named_contributions",
            "residual", "closure_error_dex",
        }
        if set(gap) != expected_fields or gap.get("classification") != expected_ratio["classification"]:
            raise ValueError(f"{name}.{quantity} gap classification/fields are inconsistent")
        ratio_value = expected_ratio["value"]
        is_eligible = expected_ratio["classification"] == "available" and ratio_value > 0.0
        if not is_eligible:
            residual = _require_mapping(gap.get("residual"), f"{name}.{quantity}.residual")
            expected_residual = {
                "name": "cross_solver_semantics_residual",
                "classification": "unidentifiable",
                "value_dex": None,
            }
            if (gap.get("log_gap_dex") is not None or gap.get("named_contributions") != []
                    or residual != expected_residual or gap.get("closure_error_dex") is not None):
                raise ValueError(f"{name}.{quantity} ineligible gap carries fabricated closure")
            continue
        eligible += 1
        expected_log_gap = math.log10(abs(ratio_value))
        log_gap = _number(gap.get("log_gap_dex"), f"{name}.{quantity}.log_gap_dex")
        if log_gap != expected_log_gap:
            raise ValueError(f"{name}.{quantity} log gap does not match the serialized ratio")
        contributions = _require_list(gap.get("named_contributions"), f"{name}.{quantity}.named_contributions")
        names, contribution_sum = set(), 0.0
        if not contributions:
            raise ValueError(f"{name}.{quantity} eligible gap lacks named contributions")
        for contribution_index, contribution_value in enumerate(contributions):
            contribution = _require_mapping(contribution_value, f"{name}.{quantity}.named_contributions[{contribution_index}]")
            if set(contribution) != {"name", "contribution_dex"}:
                raise ValueError(f"{name}.{quantity} contribution has invalid fields")
            contribution_name = contribution.get("name")
            if not isinstance(contribution_name, str) or not contribution_name or contribution_name in names:
                raise ValueError(f"{name}.{quantity} contribution names must be unique and non-empty")
            names.add(contribution_name)
            contribution_sum += _number(contribution.get("contribution_dex"), f"{name}.{quantity}.contribution_dex")
        residual = _require_mapping(gap.get("residual"), f"{name}.{quantity}.residual")
        if set(residual) != {"name", "classification", "value_dex"}:
            raise ValueError(f"{name}.{quantity} residual has invalid fields")
        if residual.get("name") != "cross_solver_semantics_residual" or residual.get("classification") != "available":
            raise ValueError(f"{name}.{quantity} residual must be named and available")
        residual_value = _number(residual.get("value_dex"), f"{name}.{quantity}.residual.value_dex")
        recomputed_error = log_gap - (contribution_sum + residual_value)
        closure_error = _number(gap.get("closure_error_dex"), f"{name}.{quantity}.closure_error_dex")
        if closure_error != recomputed_error or abs(recomputed_error) > _GAP_TOLERANCE_DEX:
            raise ValueError(f"{name}.{quantity} gap closure exceeds tolerance or was tampered")
    return eligible


def _path_definition(value, name: str) -> tuple[str, ...]:
    definition = _require_mapping(value, name, nonempty=True)
    missing = {"order", "contributions"} - set(definition)
    if missing:
        raise ValueError(f"{name} lacks explicit path fields: {sorted(missing)}")
    order = _require_list(definition["order"], f"{name}.order")
    if not order or any(not isinstance(factor, str) or not factor for factor in order):
        raise ValueError(f"{name}.order must contain named factors")
    if len(order) != len(set(order)):
        raise ValueError(f"{name}.order contains duplicate factors")
    _require_list(definition["contributions"], f"{name}.contributions")
    return tuple(order)


def _interaction_record(value, name: str, *, expected_identity=None, require_identity=False):
    row = _require_mapping(value, name)
    required = {
        "first_factor", "second_factor", "path_identity", "baseline", "a_only",
        "b_only", "both", "interaction_dex",
    }
    missing = required - set(row)
    if missing:
        raise ValueError(f"{name} lacks interaction fields: {sorted(missing)}")
    for field in ("first_factor", "second_factor", "path_identity"):
        if not isinstance(row[field], str) or not row[field]:
            raise ValueError(f"{name}.{field} must be a non-empty string")
    for field in ("baseline", "a_only", "b_only", "both", "interaction_dex"):
        _number(row[field], f"{name}.{field}")
    identity = expected_identity
    if require_identity:
        identity = _state_identity(row, name)
    elif "topology" in row or "bias_V" in row:
        supplied = _state_identity(row, name)
        if expected_identity is not None and supplied != expected_identity:
            raise ValueError(f"{name} state identity does not match its waterfall path")
        identity = supplied
    return identity, row["first_factor"], row["second_factor"], row["path_identity"]


def _quantity_ledger_result(value, state_sha256: str, name: str) -> None:
    result = _require_mapping(value, name)
    expected_fields = {"state_sha256", "status", "dominant_factor", "ranking", "closure"}
    if set(result) != expected_fields:
        raise ValueError(f"{name} has invalid ledger-result fields")
    _sha(result.get("state_sha256"), f"{name}.state_sha256")
    if result["state_sha256"] != state_sha256:
        raise ValueError(f"{name} state hash does not match its checkpoint")
    if result.get("status") != "available":
        raise ValueError(f"{name} status is inconsistent")
    ranking = _require_list(result.get("ranking"), f"{name}.ranking")
    if (not ranking or any(not isinstance(factor, str) or not factor for factor in ranking)
            or len(ranking) != len(set(ranking))):
        raise ValueError(f"{name}.ranking must be non-empty, named, and unique")
    if result.get("dominant_factor") != ranking[0]:
        raise ValueError(f"{name}.dominant_factor must equal ranking[0]")
    closure = _require_mapping(result.get("closure"), f"{name}.closure")
    if set(closure) != {"status", "tolerance_dex", "closure_error_dex"} or closure.get("status") != "closed":
        raise ValueError(f"{name}.closure has invalid fields or status")
    tolerance = _number(closure.get("tolerance_dex"), f"{name}.closure.tolerance_dex")
    error = _number(closure.get("closure_error_dex"), f"{name}.closure.closure_error_dex")
    if tolerance != _GAP_TOLERANCE_DEX or abs(error) > tolerance:
        raise ValueError(f"{name}.closure error exceeds its canonical tolerance")


def _accepted_transition(row, name: str, *, strict_package: bool) -> tuple[str, str, float]:
    row = _require_mapping(row, name)
    solver = row.get("solver")
    topology = row.get("topology")
    if solver not in {"vela", "sentaurus"} or topology not in {"sketch", "mirror"}:
        raise ValueError(f"{name} lacks typed solver/topology identity")
    if row.get("status") != "accepted":
        raise ValueError(f"{name} must have accepted status")
    target = _number(row.get("target_bias_V"), f"{name}.target_bias_V")
    actual = _number(row.get("actual_bias_V"), f"{name}.actual_bias_V")
    if not _within_bias_tolerance(target, actual):
        raise ValueError(f"{name} is not an exact checkpoint")
    _observables(row.get("observables"), f"{name}.observables")
    _branch_evidence(row, name)
    if strict_package:
        if not isinstance(row.get("state_path"), str) or not row["state_path"]:
            raise ValueError(f"{name} lacks a state path")
        _sha(row.get("state_sha256"), f"{name}.state_sha256")
    if "quantity_ledger_result" in row:
        _quantity_ledger_result(
            row["quantity_ledger_result"], row.get("state_sha256"),
            f"{name}.quantity_ledger_result",
        )
    return str(solver), str(topology), target


def _failed_transition(row, name: str, *, strict_package: bool) -> tuple[str, str, float]:
    row = _require_mapping(row, name)
    solver, topology = row.get("solver"), row.get("topology")
    if solver not in {"vela", "sentaurus"} or topology not in {"sketch", "mirror"}:
        raise ValueError(f"{name} lacks typed solver/topology identity")
    if row.get("status") != "rejected" or row.get("observables") is not None:
        raise ValueError("failed transition must preserve no fabricated observables")
    target = _number(row.get("target_bias_V"), f"{name}.target_bias_V")
    _branch_evidence(row, name)
    if row["branch_classification"] != "unidentified":
        raise ValueError("failed transition branch must be unidentified")
    return str(solver), str(topology), target


_FIGURE_NAMES = (
    "terminal_current.png", "one_volt_growth.png", "maximum_field.png",
    "source_integrals.png", "topology.png",
)
_FIGURE_METADATA = {
    "DiagnosticDisclaimer": DISCLAIMER,
    "SolverTermination": "Every recorded solver failure transition is explicitly marked.",
    "BVExtrapolation": "No physical breakdown voltage (BV) is extrapolated.",
}


def _validate_figure_contract(report: Mapping) -> None:
    contract = _require_mapping(report.get("figure_contract"), "figure_contract")
    if set(contract) != {"schema", "figures"} or contract.get("schema") != "vela.pn2d_minimal6_figure_contract.v1":
        raise ValueError("figure_contract has an invalid schema")
    entries = _require_mapping(contract.get("figures"), "figure_contract.figures")
    if set(entries) != set(_FIGURE_NAMES):
        raise ValueError("figure_contract requires the exact five figures")
    markers = [{key: row[key] for key in ("solver", "topology", "start_bias_V", "target_bias_V")} for row in report["failure_transitions"]]
    accepted = report["accepted_transitions"]
    def identities(quantity):
        return [{"solver": solver, "topology": topology, "quantity": quantity}
                for solver in ("vela", "sentaurus") for topology in ("sketch", "mirror")
                if any(row["topology"] == topology for row in accepted[solver])]
    expected = {
        "terminal_current.png": identities("terminal_current"),
        "maximum_field.png": identities("maximum_field"),
        "source_integrals.png": [
            {"solver": solver, "topology": topology, "quantity": quantity}
            for solver in ("vela", "sentaurus") for topology in ("sketch", "mirror")
            for quantity in ("native_source", "reconstructed_source")
            if any(row["topology"] == topology for row in accepted[solver])
        ],
        "one_volt_growth.png": [
            {"solver": solver, "quantity": "one_volt_growth"}
            for solver, field in (("vela", "vela_one_volt_current_growth"), ("sentaurus", "sentaurus_one_volt_current_growth"))
            if any(row[field]["classification"] == "available" for row in report["checkpoints"])
        ],
        "topology.png": [
            {"solver": solver, "quantity": "terminal_current_sketch_over_mirror"}
            for solver in ("vela", "sentaurus")
            if any(row["solver"] == solver and row["terminal_current_sketch_over_mirror"]["classification"] == "available" for row in report["topology_sensitivity"])
        ],
    }
    hashes = _require_mapping(report["artifact_hashes"], "artifact_hashes")
    entry_fields = {"sha256", "width_px", "height_px", "metadata", "series_identities", "failure_transition_markers"}
    for name in _FIGURE_NAMES:
        entry = _require_mapping(entries[name], f"figure_contract.figures.{name}")
        if set(entry) != entry_fields:
            raise ValueError(f"figure contract entry has invalid fields: {name}")
        _sha(entry.get("sha256"), f"figure_contract.figures.{name}.sha256")
        if entry["sha256"] != hashes.get(name):
            raise ValueError(f"figure contract hash mismatch: {name}")
        if entry["width_px"] != 900 or entry["height_px"] != 504 or entry["width_px"] < 640 or entry["height_px"] < 360:
            raise ValueError(f"figure contract dimensions mismatch: {name}")
        if entry["metadata"] != _FIGURE_METADATA:
            raise ValueError(f"figure contract metadata mismatch: {name}")
        if entry["series_identities"] != expected[name]:
            raise ValueError(f"figure contract series identities mismatch: {name}")
        if entry["failure_transition_markers"] != markers:
            raise ValueError(f"figure contract failure markers mismatch: {name}")


def validate_formula_difference_v1(report: dict) -> None:
    fields = {
        "input_provenance", "audit_provenance", "state_matrix", "row_counts",
        "waterfall_paths", "interactions", "dominance_rules",
        "sentaurus_internal_semantics_residual", "vela_parameter_agreement",
        "artifact_hashes", "records",
    }
    _require_fields(report, "vela.pn2d_minimal6_formula_difference.v1", fields)
    _require_mapping(report["input_provenance"], "input_provenance", nonempty=True)
    _require_mapping(report["audit_provenance"], "audit_provenance", nonempty=True)
    identities = {
        _state_identity(_require_mapping(row, f"state_matrix[{index}]"), f"state_matrix[{index}]")
        for index, row in enumerate(_require_list(report["state_matrix"], "state_matrix"))
    }
    if identities != _EXPECTED_STATES or len(report["state_matrix"]) != 6:
        raise ValueError("formula report requires the exact six-state matrix")
    if report["row_counts"] != {"node": 36, "edge": 54, "triangle": 24}:
        raise ValueError("formula report requires exact 36/54/24 row counts")
    paths = _require_list(report["waterfall_paths"], "waterfall_paths")
    path_ids = set()
    for index, path in enumerate(paths):
        path = _require_mapping(path, f"waterfall_paths[{index}]")
        missing = {"topology", "bias_V", "forward", "reverse", "interactions", "residual_dex", "status"} - set(path)
        if missing:
            raise ValueError(f"waterfall path lacks definitions: {sorted(missing)}")
        identity = _state_identity(path, f"waterfall_paths[{index}]")
        path_ids.add(identity)
        forward = _path_definition(path["forward"], f"waterfall_paths[{index}].forward")
        reverse = _path_definition(path["reverse"], f"waterfall_paths[{index}].reverse")
        if set(forward) != set(reverse):
            raise ValueError("forward and reverse paths must cover the same named factors")
        interactions = _require_list(path["interactions"], f"waterfall_paths[{index}].interactions")
        seen_path_interactions = set()
        for interaction_index, interaction in enumerate(interactions):
            key = _interaction_record(
                interaction,
                f"waterfall_paths[{index}].interactions[{interaction_index}]",
                expected_identity=identity,
            )
            if key in seen_path_interactions:
                raise ValueError("duplicate waterfall interaction record")
            seen_path_interactions.add(key)
    if path_ids != _EXPECTED_STATES or len(paths) != 6:
        raise ValueError("waterfall paths must cover the exact six-state matrix")
    interactions = _require_list(report["interactions"], "interactions")
    seen_interactions = set()
    for index, interaction in enumerate(interactions):
        key = _interaction_record(interaction, f"interactions[{index}]", require_identity=True)
        if key in seen_interactions:
            raise ValueError("duplicate state interaction record")
        seen_interactions.add(key)
    _require_mapping(report["dominance_rules"], "dominance_rules", nonempty=True)
    residuals = _require_list(report["sentaurus_internal_semantics_residual"], "sentaurus_internal_semantics_residual")
    residual_ids = set()
    for index, residual in enumerate(residuals):
        residual = _require_mapping(residual, f"sentaurus_internal_semantics_residual[{index}]")
        if residual.get("name") != "sentaurus_internal_semantics_residual":
            raise ValueError("residual record has the wrong name")
        residual_ids.add(
            _state_identity(residual, f"sentaurus_internal_semantics_residual[{index}]")
        )
    if len(residuals) != 6 or residual_ids != _EXPECTED_STATES:
        raise ValueError("named residual must cover all six states")
    hashes = _require_mapping(report["artifact_hashes"], "artifact_hashes", nonempty=True)
    _sha(hashes.get("state_manifest_sha256"), "artifact_hashes.state_manifest_sha256")
    _require_list(report["records"], "records")


def validate_bv_comparison_v1(report: dict) -> None:
    fields = {
        "solver_configurations", "accepted_transitions", "failed_transitions",
        "checkpoints", "terminal_currents", "maximum_fields", "source_integrals",
        "convergence_metadata", "curve_artifact_hashes", "records",
        "branch_threshold_version", "deepest_common_bias_V", "missing_tails",
        "topology_sensitivity", "fixed_state_recheck", "failure_transitions",
        "artifact_hashes", "input_artifacts", "closure", "interpolation",
    }
    _require_fields(report, "vela.pn2d_minimal6_bv_comparison.v1", fields)
    _reject_raw_states(report)
    if report["interpolation"] != "forbidden":
        raise ValueError("comparison must forbid interpolation")
    version = report["branch_threshold_version"]
    if version != _BRANCH_THRESHOLD_VERSION:
        raise ValueError("missing or noncanonical branch threshold version")
    configurations = _require_mapping(report["solver_configurations"], "solver_configurations")
    if set(configurations) != {"vela", "sentaurus"}:
        raise ValueError("solver_configurations must cover both solvers")
    for solver in ("vela", "sentaurus"):
        _solver_configuration(configurations[solver], f"solver_configurations.{solver}")
    accepted = _require_mapping(report["accepted_transitions"], "accepted_transitions")
    if set(accepted) != {"vela", "sentaurus"}:
        raise ValueError("accepted_transitions must be keyed by both solvers")
    for solver in ("vela", "sentaurus"):
        seen = set()
        for index, row in enumerate(_require_list(accepted[solver], f"accepted_transitions.{solver}")):
            identity = _accepted_transition(row, f"accepted_transitions.{solver}[{index}]", strict_package=True)
            if identity[0] != solver or identity[1:] in seen:
                raise ValueError("accepted transition solver/identity mismatch or duplicate")
            if row["branch_threshold_version"] != version:
                raise ValueError("accepted transition branch threshold version mismatch")
            seen.add(identity[1:])
    for index, row in enumerate(_require_list(report["failed_transitions"], "failed_transitions")):
        _failed_transition(row, f"failed_transitions[{index}]", strict_package=False)
        if row["branch_threshold_version"] != version:
            raise ValueError("failed transition branch threshold version mismatch")
    if report["failure_transitions"] != report["failed_transitions"]:
        raise ValueError("failure transition aliases must preserve identical evidence")
    checkpoints = _require_list(report["checkpoints"], "checkpoints")
    seen = set()
    eligible_gaps = 0
    for index, row in enumerate(checkpoints):
        row = _require_mapping(row, f"checkpoints[{index}]")
        identity = _state_identity(row, f"checkpoints[{index}]")
        if identity in seen:
            raise ValueError("duplicate exact comparison checkpoint")
        seen.add(identity)
        if row.get("branch_classification") not in _BRANCHES:
            raise ValueError("checkpoint lacks typed branch classification")
        if row.get("branch_threshold_version") != version:
            raise ValueError("checkpoint branch threshold version mismatch")
        conservation = _require_mapping(row.get("contact_current_conservation"), f"checkpoints[{index}].contact_current_conservation")
        if set(conservation) != {"unit", "vela", "sentaurus"} or conservation.get("unit") != "A/um":
            raise ValueError("checkpoint contact conservation must use A/um for both solvers")
        for solver in ("vela", "sentaurus"):
            _contact_conservation(conservation[solver], f"checkpoints[{index}].contact_current_conservation.{solver}")
        for solver in ("vela", "sentaurus"):
            nested = _accepted_transition(row.get(solver), f"checkpoints[{index}].{solver}", strict_package=True)
            if row[solver]["branch_classification"] != "unidentified":
                raise ValueError("nested accepted transition branch evidence is noncanonical")
            observables = row[solver]["observables"]
            evidence = conservation[solver]
            if (evidence["anode_current_A_per_um"] != observables["anode_current_A_per_um"] or
                    evidence["cathode_current_A_per_um"] != observables["cathode_current_A_per_um"]):
                raise ValueError("contact conservation currents do not match checkpoint observables")
            if row[solver]["branch_threshold_version"] != version:
                raise ValueError("checkpoint nested branch threshold version mismatch")
            if nested != (solver, identity[0], identity[1]):
                raise ValueError("checkpoint nested solver identity mismatch")
        eligible_gaps += _validate_checkpoint_branch_and_gaps(
            row, f"checkpoints[{index}]", version
        )
    if report["artifact_hashes"] or "figure_contract" in report:
        _validate_figure_contract(report)
    for name in ("terminal_currents", "maximum_fields", "source_integrals", "records"):
        _require_list(report[name], name)
    _require_mapping(report["convergence_metadata"], "convergence_metadata", nonempty=True)
    hashes = _require_mapping(report["curve_artifact_hashes"], "curve_artifact_hashes")
    for solver in ("vela_manifest", "sentaurus_manifest"):
        _sha(hashes.get(solver), f"curve_artifact_hashes.{solver}")
    top_closure = _require_mapping(report["closure"], "closure")
    if set(top_closure) != {"status", "eligible_gaps", "rule"} or top_closure.get("status") != "closed":
        raise ValueError("comparison closure has invalid fields or status")
    if (isinstance(top_closure.get("eligible_gaps"), bool)
            or not isinstance(top_closure.get("eligible_gaps"), int)
            or top_closure["eligible_gaps"] != eligible_gaps):
        raise ValueError("comparison eligible gap count does not match recomputed gaps")
    if not isinstance(top_closure.get("rule"), str) or not top_closure["rule"]:
        raise ValueError("comparison closure lacks its named eligibility rule")


def validate_sweep_manifest_v1(report: dict) -> None:
    fields = {
        "targets_V", "template", "topology_input_sha256", "segments",
        "sentaurus_segments", "accepted_checkpoints", "failed_transition",
        "failed_transitions", "interpolation", "branch_threshold_version",
    }
    _require_fields(report, "vela.pn2d_minimal6_sweep_manifest.v1", fields)
    _reject_raw_states(report)
    if report["interpolation"] != "forbidden":
        raise ValueError("sweep manifest must forbid interpolation")
    targets = [_number(value, "targets_V") for value in _require_list(report["targets_V"], "targets_V")]
    expected = [float(-index) for index in range(len(targets))]
    if not targets or targets != expected or len(targets) > 21:
        raise ValueError("sweep manifest targets must be the exact integer prefix 0..-20 V")
    template = _require_mapping(report["template"], "template", nonempty=True)
    _sha(template.get("sha256"), "template.sha256")
    topology_hashes = _require_mapping(report["topology_input_sha256"], "topology_input_sha256", nonempty=True)
    for topology, hashes in topology_hashes.items():
        if topology not in {"sketch", "mirror"}:
            raise ValueError("topology_input_sha256 has an invalid topology")
        for name, value in _require_mapping(hashes, f"topology_input_sha256.{topology}", nonempty=True).items():
            _sha(value, f"topology_input_sha256.{topology}.{name}")
    _require_list(report["segments"], "segments")
    _require_list(report["sentaurus_segments"], "sentaurus_segments")
    version = report["branch_threshold_version"]
    if version != _BRANCH_THRESHOLD_VERSION:
        raise ValueError("missing or noncanonical branch threshold version")
    accepted_rows = _require_list(report["accepted_checkpoints"], "accepted_checkpoints")
    seen = set()
    for index, row in enumerate(accepted_rows):
        identity = _accepted_transition(row, f"accepted_checkpoints[{index}]", strict_package=True)
        if identity in seen:
            raise ValueError("duplicate accepted exact checkpoint")
        seen.add(identity)
        if row["branch_threshold_version"] != version:
            raise ValueError("accepted checkpoint branch threshold version mismatch")
    failures = _require_list(report["failed_transitions"], "failed_transitions")
    for index, row in enumerate(failures):
        _failed_transition(row, f"failed_transitions[{index}]", strict_package=True)
        if row["branch_threshold_version"] != version:
            raise ValueError("failed transition branch threshold version mismatch")
    first = report["failed_transition"]
    if first is not None:
        _failed_transition(first, "failed_transition", strict_package=True)
        if first["branch_threshold_version"] != version:
            raise ValueError("failed transition branch threshold version mismatch")
        if first not in failures:
            raise ValueError("failed_transition must be preserved in failed_transitions")
