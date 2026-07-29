#!/usr/bin/env python3
"""Fail-closed contract for paired PN2D BV process observations."""

from __future__ import annotations

import hashlib
import json
import math
import re
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any


SCHEMA_ID = "vela.pn2d_bv_process_run.v1"
EXACT_BIAS_TOLERANCE_V = 1.0e-10
BRANCHES = frozenset(
    {
        "avalanche_off",
        "iic_postprocess",
        "avalanche_on",
        "avalanche_on_aval_derivatives",
    }
)
SUPPORT_KINDS = frozenset(
    {
        "physical_node",
        "contact_support_vertex",
        "cell",
        "element_local_edge",
        "element_local_vertex",
        "contact",
    }
)
PROVENANCE = frozenset(
    {
        "native",
        "operator_replay",
        "reconstructed",
        "solver_used",
        "postprocessed",
    }
)
CARRIERS = frozenset({"electron", "hole", "total", "none"})
CANONICAL_UNITS = frozenset(
    {
        "um",
        "V",
        "cm^-3",
        "V/cm",
        "cm^2/(V s)",
        "cm/s",
        "A/cm^2",
        "cm^-1",
        "cm^-3 s^-1",
        "A/um",
        "1",
    }
)
SUPPORT_CENTERING = {
    "physical_node": "vertex",
    "contact_support_vertex": "vertex",
    "cell": "cell",
    "element_local_edge": "element_edge",
    "element_local_vertex": "element_vertex",
    "contact": "contact",
}
QUANTITY_UNITS = {
    "coordinate": "um",
    "potential": "V",
    "quasi_fermi": "V",
    "density": "cm^-3",
    "charge_density": "cm^-3",
    "doping": "cm^-3",
    "electric_field": "V/cm",
    "quasi_fermi_gradient": "V/cm",
    "mobility": "cm^2/(V s)",
    "velocity": "cm/s",
    "current_density": "A/cm^2",
    "avalanche_alpha": "cm^-1",
    "avalanche_generation": "cm^-3 s^-1",
    "srh_recombination": "cm^-3 s^-1",
    "terminal_current": "A/um",
    "integrated_source": "A/um",
    "ionization_integral": "1",
}
HASH_RE = re.compile(r"^[0-9a-f]{64}$")
REQUIRED_RUN_FIELDS = frozenset(
    {
        "schema",
        "status",
        "outcome",
        "run_id",
        "simulator",
        "release",
        "missing_value_policy",
        "input_hashes",
        "normalized_output_hashes",
        "branch_records",
        "field_records",
        "aggregate_records",
        "newton_attempt_records",
    }
)
BRANCH_RECORD_FIELDS = frozenset(
    {"branch", "requested_biases_V", "bias_records"}
)
BIAS_RECORD_FIELDS = frozenset(
    {
        "requested_bias_V",
        "actual_bias_V",
        "snapshot_tdr",
        "currentplot",
        "process_record",
    }
)
FIELD_RECORD_FIELDS = frozenset(
    {
        "branch",
        "requested_bias_V",
        "actual_bias_V",
        "support_kind",
        "support_key",
        "centering",
        "provenance",
        "carrier",
        "quantity",
        "components",
        "unit",
        "values",
        "coordinates_um",
        "connectivity",
        "source",
    }
)
AGGREGATE_RECORD_FIELDS = frozenset(
    {
        "branch",
        "requested_bias_V",
        "actual_bias_V",
        "carrier",
        "quantity",
        "unit",
        "value",
        "provenance",
        "source",
    }
)
NEWTON_ATTEMPT_RECORD_FIELDS = frozenset(
    {
        "branch",
        "attempt_id",
        "requested_bias_V",
        "actual_bias_V",
        "status",
        "reason",
        "source",
    }
)


class ProcessContractError(ValueError):
    """A typed process-data contract violation."""

    def __init__(self, reason: str, detail: str):
        self.reason = reason
        super().__init__(f"{reason}: {detail}")


def _fail(reason: str, detail: str) -> None:
    raise ProcessContractError(reason, detail)


def _require(condition: bool, reason: str, detail: str) -> None:
    if not condition:
        _fail(reason, detail)


def _validate_field_set(
    record: Mapping[str, Any],
    allowed: frozenset[str],
    label: str,
) -> None:
    unexpected = set(record) - allowed
    _require(
        not unexpected,
        "unexpected_record_field",
        f"{label}:{sorted(unexpected)}",
    )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def normalized_records_sha256(records: Sequence[Mapping[str, Any]]) -> str:
    """Return an order-independent hash without changing row semantics."""

    serialized = [
        json.dumps(record, sort_keys=True, separators=(",", ":"), allow_nan=False)
        for record in records
    ]
    payload = ("\n".join(sorted(serialized)) + "\n").encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _finite(value: Any, label: str) -> float:
    _require(
        isinstance(value, (int, float)) and not isinstance(value, bool),
        "non_numeric_value",
        label,
    )
    result = float(value)
    _require(math.isfinite(result), "nonfinite_value", label)
    return result


def _validate_hash_map(
    hashes: Any,
    label: str,
    *,
    base_dir: Path | None,
) -> None:
    _require(
        isinstance(hashes, Mapping) and bool(hashes),
        "missing_hashes",
        label,
    )
    for name, digest in hashes.items():
        _require(
            isinstance(name, str) and bool(name),
            "invalid_hash_name",
            label,
        )
        _require(
            isinstance(digest, str) and HASH_RE.fullmatch(digest) is not None,
            "invalid_sha256",
            f"{label}:{name}",
        )
        if base_dir is not None:
            path = base_dir / name
            _require(path.is_file(), "hashed_file_missing", str(path))
            _require(
                sha256(path) == digest,
                "hash_drift",
                str(path),
            )


def _validate_artifact(
    artifact: Any,
    label: str,
    *,
    base_dir: Path | None,
) -> None:
    _require(isinstance(artifact, Mapping), "missing_artifact", label)
    _require(
        set(artifact) == {"path", "sha256"},
        "artifact_field_set",
        label,
    )
    path_text = artifact.get("path")
    digest = artifact.get("sha256")
    _require(
        isinstance(path_text, str) and bool(path_text),
        "missing_artifact_path",
        label,
    )
    _require(
        isinstance(digest, str) and HASH_RE.fullmatch(digest) is not None,
        "invalid_sha256",
        label,
    )
    if base_dir is not None:
        path = base_dir / path_text
        _require(path.is_file(), "artifact_missing", str(path))
        _require(sha256(path) == digest, "hash_drift", str(path))


def _validate_exact_bias(
    record: Mapping[str, Any],
    label: str,
) -> tuple[float, float]:
    requested = _finite(record.get("requested_bias_V"), f"{label}.requested")
    actual = _finite(record.get("actual_bias_V"), f"{label}.actual")
    _require(
        abs(requested - actual) <= EXACT_BIAS_TOLERANCE_V,
        "nearest_bias_substitution",
        f"{label}: requested={requested:.17g}, actual={actual:.17g}",
    )
    return requested, actual


def _validate_source(source: Any, label: str) -> None:
    _require(isinstance(source, Mapping), "missing_source", label)
    _require(
        set(source) == {"file", "dataset", "index"},
        "source_field_set",
        label,
    )
    _require(
        isinstance(source.get("file"), str) and bool(source["file"]),
        "missing_source_file",
        label,
    )
    _require(
        isinstance(source.get("dataset"), str) and bool(source["dataset"]),
        "missing_source_dataset",
        label,
    )
    _require(
        isinstance(source.get("index"), int)
        and not isinstance(source["index"], bool)
        and source["index"] >= 0,
        "invalid_source_index",
        label,
    )


def _validate_carrier_quantity(record: Mapping[str, Any], label: str) -> None:
    carrier = record.get("carrier")
    _require(carrier in CARRIERS, "missing_or_unknown_carrier", label)
    quantity = record.get("quantity")
    _require(
        isinstance(quantity, str) and bool(quantity),
        "missing_quantity",
        label,
    )
    unit = record.get("unit")
    _require(unit in CANONICAL_UNITS, "unknown_unit", f"{label}:{unit!r}")
    expected_unit = QUANTITY_UNITS.get(quantity)
    _require(
        expected_unit is not None,
        "unknown_quantity",
        f"{label}:{quantity}",
    )
    _require(
        unit == expected_unit,
        "wrong_unit",
        f"{label}:{quantity} requires {expected_unit}, got {unit}",
    )
    carrier_specific = {
        "density",
        "quasi_fermi",
        "quasi_fermi_gradient",
        "mobility",
        "velocity",
        "avalanche_alpha",
    }
    if quantity in carrier_specific:
        _require(
            carrier in {"electron", "hole"},
            "missing_or_unknown_carrier",
            f"{label}:{quantity}",
        )


def _validate_field_records(records: Any) -> None:
    _require(isinstance(records, Sequence), "missing_field_records", "field_records")
    keys: set[tuple[Any, ...]] = set()
    for index, record in enumerate(records):
        label = f"field_records[{index}]"
        _require(isinstance(record, Mapping), "invalid_field_record", label)
        _validate_field_set(record, FIELD_RECORD_FIELDS, label)
        branch = record.get("branch")
        _require(branch in BRANCHES, "unknown_branch", label)
        _, actual = _validate_exact_bias(record, label)
        support = record.get("support_kind")
        _require(support in SUPPORT_KINDS, "unknown_support", label)
        support_key = record.get("support_key")
        _require(
            isinstance(support_key, str) and bool(support_key),
            "missing_support_key",
            label,
        )
        centering = record.get("centering")
        _require(
            centering == SUPPORT_CENTERING[support],
            "wrong_centering",
            f"{label}:{support} requires {SUPPORT_CENTERING[support]}",
        )
        provenance = record.get("provenance")
        _require(provenance in PROVENANCE, "unknown_provenance", label)
        _validate_carrier_quantity(record, label)
        components = record.get("components")
        values = record.get("values")
        _require(
            isinstance(components, Sequence)
            and not isinstance(components, (str, bytes))
            and bool(components),
            "missing_components",
            label,
        )
        _require(
            all(
                isinstance(component, str) and bool(component)
                for component in components
            ),
            "invalid_component",
            label,
        )
        _require(
            len(set(components)) == len(components),
            "duplicate_components",
            label,
        )
        _require(
            isinstance(values, Sequence)
            and not isinstance(values, (str, bytes))
            and len(values) == len(components),
            "implicit_zero_fill",
            label,
        )
        for value_index, value in enumerate(values):
            _finite(value, f"{label}.values[{value_index}]")
        if (
            support == "element_local_edge"
            and provenance == "native"
            and record.get("quantity") == "current_density"
        ):
            _fail(
                "unsupported_native_edge_claim",
                label,
            )
        coordinates = record.get("coordinates_um")
        connectivity = record.get("connectivity")
        if support in {"physical_node", "contact_support_vertex"}:
            _require(
                isinstance(coordinates, Sequence)
                and not isinstance(coordinates, (str, bytes))
                and len(coordinates) >= 2,
                "missing_support_coordinates",
                label,
            )
        if support in {"cell", "element_local_edge", "element_local_vertex"}:
            _require(
                isinstance(connectivity, Sequence)
                and not isinstance(connectivity, (str, bytes))
                and bool(connectivity),
                "missing_support_connectivity",
                label,
            )
        if coordinates is not None:
            for coordinate_index, coordinate in enumerate(coordinates):
                _finite(coordinate, f"{label}.coordinates_um[{coordinate_index}]")
        if connectivity is not None:
            _require(
                all(
                    isinstance(node, int)
                    and not isinstance(node, bool)
                    and node >= 0
                    for node in connectivity
                ),
                "invalid_support_connectivity",
                label,
            )
        key = (
            branch,
            actual,
            support,
            support_key,
            provenance,
            record.get("carrier"),
            record.get("quantity"),
            tuple(components),
        )
        _require(key not in keys, "duplicate_support_key", label)
        keys.add(key)
        _validate_source(record.get("source"), label)


def _validate_aggregate_records(records: Any) -> None:
    _require(
        isinstance(records, Sequence),
        "missing_aggregate_records",
        "aggregate_records",
    )
    keys: set[tuple[Any, ...]] = set()
    for index, record in enumerate(records):
        label = f"aggregate_records[{index}]"
        _require(isinstance(record, Mapping), "invalid_aggregate_record", label)
        _validate_field_set(record, AGGREGATE_RECORD_FIELDS, label)
        _require(record.get("branch") in BRANCHES, "unknown_branch", label)
        _, actual = _validate_exact_bias(record, label)
        _require(
            record.get("provenance") in PROVENANCE,
            "unknown_provenance",
            label,
        )
        _validate_carrier_quantity(record, label)
        _finite(record.get("value"), f"{label}.value")
        key = (
            record.get("branch"),
            actual,
            record.get("carrier"),
            record.get("quantity"),
            record.get("provenance"),
        )
        _require(key not in keys, "duplicate_aggregate_key", label)
        keys.add(key)
        _validate_source(record.get("source"), label)


def _validate_newton_attempt_records(records: Any) -> None:
    _require(
        isinstance(records, Sequence),
        "missing_newton_attempt_records",
        "newton_attempt_records",
    )
    keys: set[tuple[str, str]] = set()
    for index, record in enumerate(records):
        label = f"newton_attempt_records[{index}]"
        _require(isinstance(record, Mapping), "invalid_newton_attempt", label)
        _validate_field_set(record, NEWTON_ATTEMPT_RECORD_FIELDS, label)
        branch = record.get("branch")
        _require(branch in BRANCHES, "unknown_branch", label)
        _validate_exact_bias(record, label)
        attempt_id = record.get("attempt_id")
        _require(
            isinstance(attempt_id, str) and bool(attempt_id),
            "missing_attempt_id",
            label,
        )
        key = (branch, attempt_id)
        _require(key not in keys, "duplicate_attempt_id", label)
        keys.add(key)
        _require(
            record.get("status") in {"accepted", "rejected", "failed"},
            "unknown_attempt_status",
            label,
        )
        _require(
            isinstance(record.get("reason"), str) and bool(record["reason"]),
            "missing_attempt_reason",
            label,
        )
        _validate_source(record.get("source"), label)


def validate_process_run(
    manifest: Mapping[str, Any],
    *,
    base_dir: Path | None = None,
) -> None:
    """Validate manifest semantics and optionally close every file hash."""

    _require(
        manifest.get("schema") == SCHEMA_ID,
        "schema_mismatch",
        str(manifest.get("schema")),
    )
    _require(
        set(manifest) == REQUIRED_RUN_FIELDS,
        "run_field_set_mismatch",
        str(sorted(set(manifest) ^ REQUIRED_RUN_FIELDS)),
    )
    _require(manifest.get("status") in {"passed", "failed"}, "invalid_status", "")
    for label in ("outcome", "run_id", "release"):
        _require(
            isinstance(manifest.get(label), str) and bool(manifest[label]),
            f"missing_{label}",
            label,
        )
    _require(
        manifest.get("simulator") in {"vela", "sentaurus"},
        "unknown_simulator",
        str(manifest.get("simulator")),
    )
    _require(
        manifest.get("missing_value_policy") == "reject",
        "implicit_zero_fill",
        "missing_value_policy must be reject",
    )
    _validate_hash_map(
        manifest.get("input_hashes"),
        "input_hashes",
        base_dir=base_dir,
    )
    _validate_hash_map(
        manifest.get("normalized_output_hashes"),
        "normalized_output_hashes",
        base_dir=base_dir,
    )

    branches = manifest.get("branch_records")
    _require(
        isinstance(branches, Sequence)
        and not isinstance(branches, (str, bytes))
        and bool(branches),
        "missing_branch_records",
        "",
    )
    seen_branches: set[str] = set()
    declared_biases: set[tuple[str, float]] = set()
    for branch_index, branch_record in enumerate(branches):
        label = f"branch_records[{branch_index}]"
        _require(isinstance(branch_record, Mapping), "invalid_branch_record", label)
        _validate_field_set(branch_record, BRANCH_RECORD_FIELDS, label)
        branch = branch_record.get("branch")
        _require(branch in BRANCHES, "unknown_branch", label)
        _require(branch not in seen_branches, "duplicate_branch", str(branch))
        seen_branches.add(branch)
        requested_values = branch_record.get("requested_biases_V")
        _require(
            isinstance(requested_values, Sequence)
            and not isinstance(requested_values, (str, bytes))
            and bool(requested_values),
            "missing_requested_biases",
            label,
        )
        requested = tuple(_finite(value, label) for value in requested_values)
        _require(
            len(set(requested)) == len(requested),
            "duplicate_requested_bias",
            label,
        )
        bias_records = branch_record.get("bias_records")
        _require(
            isinstance(bias_records, Sequence)
            and not isinstance(bias_records, (str, bytes)),
            "missing_bias_records",
            label,
        )
        observed: list[float] = []
        for bias_index, bias_record in enumerate(bias_records):
            bias_label = f"{label}.bias_records[{bias_index}]"
            _require(
                isinstance(bias_record, Mapping),
                "invalid_bias_record",
                bias_label,
            )
            _validate_field_set(bias_record, BIAS_RECORD_FIELDS, bias_label)
            requested_bias, actual_bias = _validate_exact_bias(
                bias_record,
                bias_label,
            )
            observed.append(requested_bias)
            declared_biases.add((branch, actual_bias))
            for artifact_name in ("snapshot_tdr", "currentplot", "process_record"):
                _validate_artifact(
                    bias_record.get(artifact_name),
                    f"{bias_label}.{artifact_name}",
                    base_dir=base_dir,
                )
        _require(
            len(set(observed)) == len(observed),
            "duplicate_bias_row",
            label,
        )
        _require(
            set(observed) == set(requested),
            "missing_bias_row",
            f"{label}: requested={requested}, observed={tuple(observed)}",
        )

    for record_group in (
        manifest.get("field_records", ()),
        manifest.get("aggregate_records", ()),
        manifest.get("newton_attempt_records", ()),
    ):
        for record in record_group:
            if isinstance(record, Mapping):
                pair = (
                    record.get("branch"),
                    float(record.get("actual_bias_V", math.nan)),
                )
                _require(
                    pair in declared_biases,
                    "undeclared_bias_record",
                    str(pair),
                )

    _validate_field_records(manifest.get("field_records"))
    _validate_aggregate_records(manifest.get("aggregate_records"))
    _validate_newton_attempt_records(manifest.get("newton_attempt_records"))
