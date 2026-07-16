"""Strict runtime validation for Minimal6 diagnostic report contracts."""

from collections.abc import Mapping
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
        if not bias_value - 1.0e-12 <= actual <= bias_value + 1.0e-12:
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
    if abs(target - actual) > 1.0e-12:
        raise ValueError(f"{name} is not an exact checkpoint")
    _observables(row.get("observables"), f"{name}.observables")
    _branch_evidence(row, name)
    if strict_package:
        if not isinstance(row.get("state_path"), str) or not row["state_path"]:
            raise ValueError(f"{name} lacks a state path")
        _sha(row.get("state_sha256"), f"{name}.state_sha256")
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
    if not isinstance(version, str) or not version:
        raise ValueError("missing branch threshold version")
    accepted = _require_mapping(report["accepted_transitions"], "accepted_transitions")
    if set(accepted) != {"vela", "sentaurus"}:
        raise ValueError("accepted_transitions must be keyed by both solvers")
    for solver in ("vela", "sentaurus"):
        seen = set()
        for index, row in enumerate(_require_list(accepted[solver], f"accepted_transitions.{solver}")):
            identity = _accepted_transition(row, f"accepted_transitions.{solver}[{index}]", strict_package=False)
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
        for solver in ("vela", "sentaurus"):
            nested = _accepted_transition(row.get(solver), f"checkpoints[{index}].{solver}", strict_package=False)
            if row[solver]["branch_threshold_version"] != version:
                raise ValueError("checkpoint nested branch threshold version mismatch")
            if nested != (solver, identity[0], identity[1]):
                raise ValueError("checkpoint nested solver identity mismatch")
    for name in ("terminal_currents", "maximum_fields", "source_integrals", "records"):
        _require_list(report[name], name)
    _require_mapping(report["convergence_metadata"], "convergence_metadata", nonempty=True)
    hashes = _require_mapping(report["curve_artifact_hashes"], "curve_artifact_hashes")
    for solver in ("vela_manifest", "sentaurus_manifest"):
        _sha(hashes.get(solver), f"curve_artifact_hashes.{solver}")
    if _require_mapping(report["closure"], "closure").get("status") != "closed":
        raise ValueError("comparison closure is not closed")


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
    if not isinstance(version, str) or not version:
        raise ValueError("missing branch threshold version")
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
