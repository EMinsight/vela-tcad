#!/usr/bin/env python3
"""Export exact-bias Sentaurus states on both PN2D minimal6 topologies."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path
from typing import Callable, Sequence


REPO = Path(__file__).resolve().parents[1]
if __package__ in (None, ""):
    sys.path.insert(0, str(REPO))

from scripts.pn2d_minimal6_topology import load_topology  # noqa: E402
from scripts.run_pn2d_minimal6_sentaurus_gate import (  # noqa: E402
    DEFAULT_REMOTE_ROOT,
    MODELS_SOURCE,
    TOPOLOGY_FIXTURE,
    build_gate_bundle,
)
from scripts.run_sentaurus_vm_reference import (  # noqa: E402
    default_windows_openssh,
    run_checked,
    write_manifest,
)
from scripts.sentaurus_import import parse_quoted_list, parse_values_block  # noqa: E402


SCHEMA = "vela.pn2d_minimal6_states.v1"
V2_SCHEMA = "vela.pn2d_minimal6_states.v2"
V2_SCHEMA_PATH = REPO / "schemas" / "vela.pn2d_minimal6_states.v2.schema.json"
V2_FIELD_CONTRACT = (
    ("sentaurus_electrostatic_potential", "ElectrostaticPotential", 1, "V",
     "electrostatic_potential"),
    ("sentaurus_electron_quasi_fermi_potential", "eQuasiFermiPotential", 1, "V",
     "carrier_quasi_fermi_potential"),
    ("sentaurus_hole_quasi_fermi_potential", "hQuasiFermiPotential", 1, "V",
     "carrier_quasi_fermi_potential"),
    ("sentaurus_electron_density", "eDensity", 1, "cm^-3", "carrier_density"),
    ("sentaurus_hole_density", "hDensity", 1, "cm^-3", "carrier_density"),
    ("sentaurus_electric_field", "ElectricField", 2, "V*cm^-1", "electric_field"),
    ("sentaurus_electron_current_density", "eCurrentDensity", 2, "A*cm^-2",
     "carrier_current_density"),
    ("sentaurus_hole_current_density", "hCurrentDensity", 2, "A*cm^-2",
     "carrier_current_density"),
    ("sentaurus_electron_mobility", "eMobility", 1, "cm^2*V^-1*s^-1",
     "carrier_mobility"),
    ("sentaurus_hole_mobility", "hMobility", 1, "cm^2*V^-1*s^-1",
     "carrier_mobility"),
    ("sentaurus_electron_avalanche_coefficient", "eAlphaAvalanche", 1, "cm^-1",
     "impact_ionization_coefficient"),
    ("sentaurus_hole_avalanche_coefficient", "hAlphaAvalanche", 1, "cm^-1",
     "impact_ionization_coefficient"),
    ("sentaurus_native_avalanche_generation", "ImpactIonization", 1,
     "cm^-3*s^-1", "native_avalanche_generation"),
    ("sentaurus_electron_speed", "eVelocity", 1, "cm*s^-1", "carrier_speed"),
    ("sentaurus_hole_speed", "hVelocity", 1, "cm*s^-1", "carrier_speed"),
    ("sentaurus_electron_ionization_integral", "eIonIntegral", 1, "1",
     "path_ionization_integral"),
    ("sentaurus_hole_ionization_integral", "hIonIntegral", 1, "1",
     "path_ionization_integral"),
    ("sentaurus_mean_ionization_integral", "MeanIonIntegral", 1, "1",
     "path_ionization_integral"),
)
REQUIRED_TOPOLOGIES = ("sketch", "mirror")
REQUIRED_BIASES = (0.0, -12.0, -19.0)
BIAS_TOLERANCE_V = 1.0e-12
LOCAL_RECOVERY_RUN_ID = "minimal6_states_live_20260713_v2"
LOCAL_RECOVERY_MANIFEST_SHA256 = (
    "b44ad95d5df6d57383ba3d5b292818568e358d67f0fc0424ee72f95b673e8aaa"
)
AUDIT_PRODUCER = "build-release/pn2d_minimal6_operator_audit.exe"
AUDIT_PATH_OPTIONS = (
    ("--mesh", "mesh.json"),
    ("--doping", "doping.csv"),
    ("--state", "state.csv"),
    ("--config", "audit.json"),
    ("--node-out", "vela_node_state.csv"),
    ("--edge-out", "vela_edge_audit.csv"),
    ("--triangle-out", "vela_triangle_audit.csv"),
)
AUDIT_INPUT_COUNT = 4
SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")

COORDINATE_TOLERANCE_UM = 1.0e-12
SOURCE_DECK = (
    REPO / "reference_tcad" / "pn2d_sentaurus2018_minimal6" / "source"
    / "pn2d_minimal6_state_sdevice.cmd"
)
DEFAULT_OUTPUT_DIR = (
    REPO / "build-release" / "reference_tcad" / "pn2d_sentaurus2018_minimal6"
    / "state_exports"
)
RECOVERY_VALIDATION_SCHEMA = "vela.pn2d_minimal6_recovery_validation.v1"
SEALED_SOURCE_VALIDATION_PATH = "source_recovery_validation.json"
LOCAL_RECOVERY_VALIDATION_PATH = (
    REPO / "build-release" / "reference_tcad" / "pn2d_sentaurus2018_minimal6"
    / "recovery_validation" / LOCAL_RECOVERY_RUN_ID / "recovery_validation.json"
)
LOCAL_RECOVERY_VALIDATION_SHA256 = (
    "9466ee2db317fda4707254403fa33c2e1ebee666f8fdc72b776fa4a8ab689ec3"
)
REQUIRED_RECOVERY_FLAGS = {
    "outputs_complete": True,
    "exact_state_matrix": True,
    "field_contract": True,
    "member_hashes_verified": True,
}

DEFAULT_IMPORTER = REPO / "build-release" / (
    "sentaurus_import.exe" if os.name == "nt" else "sentaurus_import"
)
_SAFE_RUN_ID = re.compile(r"^[A-Za-z0-9_.-]+$")
_SAFE_REMOTE_COMPONENT = re.compile(r"^[A-Za-z0-9_./~:-]+$")

_FIELD_CONTRACT = {
    "ElectrostaticPotential": (1, "V"),
    "eQuasiFermiPotential": (1, "V"),
    "hQuasiFermiPotential": (1, "V"),
    "eDensity": (1, "cm^-3"),
    "hDensity": (1, "cm^-3"),
    "ElectricField": (2, "V*cm^-1"),
    "eCurrentDensity": (2, "A*cm^-2"),
    "hCurrentDensity": (2, "A*cm^-2"),
    "eMobility": (1, "cm^2*V^-1*s^-1"),
    "hMobility": (1, "cm^2*V^-1*s^-1"),
    "eAlphaAvalanche": (1, "cm^-1"),
    "hAlphaAvalanche": (1, "cm^-1"),
    "ImpactIonization": (1, "cm^-3*s^-1"),
    "eVelocity": (1, "cm*s^-1"),
    "hVelocity": (1, "cm*s^-1"),
    "eIonIntegral": (1, "1"),
    "hIonIntegral": (1, "1"),
    "MeanIonIntegral": (1, "1"),
}

def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()

def collect_member_hashes(root: Path) -> dict[str, str]:
    root = Path(root)
    return {path.relative_to(root).as_posix(): _sha256(path)
            for path in sorted(root.rglob("*")) if path.is_file()}

def validate_member_hashes(root: Path, expected: dict[str, str]) -> None:
    actual = collect_member_hashes(root)
    if set(actual) != set(expected):
        raise ValueError("archive member set mismatch")
    for name, digest in expected.items():
        if actual[name] != digest:
            raise ValueError(f"archive member hash mismatch: {name}")


def _source_member_hashes(root: Path, validation_path: Path) -> dict[str, str]:
    root = root.resolve()
    validation_path = validation_path.resolve()
    members: dict[str, str] = {}
    for path in sorted(root.rglob("*")):
        if path.is_file() and path.resolve() != validation_path:
            members[path.relative_to(root).as_posix()] = _sha256(path)
    return members


def validate_source_recovery(
    source_root: Path,
    source_manifest: dict[str, object],
    expected_manifest_sha256: str,
    *,
    source_kind: str,
    expected_source_validation_sha256: str | None = None,
) -> tuple[dict[str, object], Path, str]:
    source_root = source_root.resolve()
    if source_kind == "local_recovery":
        required_root = (DEFAULT_OUTPUT_DIR / LOCAL_RECOVERY_RUN_ID).resolve()
        if source_root != required_root:
            raise ValueError("local recovery archive root identity mismatch")
        validation_path = LOCAL_RECOVERY_VALIDATION_PATH.resolve()
        expected_validation_hash = LOCAL_RECOVERY_VALIDATION_SHA256
    elif source_kind == "regenerated":
        validation_path = (source_root / "recovery_validation.json").resolve()
        if expected_source_validation_sha256 is None:
            if not validation_path.is_file():
                raise ValueError("source member ledger recovery validation is missing")
            expected_validation_hash = _sha256(validation_path)
        else:
            expected_validation_hash = expected_source_validation_sha256.lower()
    else:
        raise ValueError(f"unsupported source kind: {source_kind}")

    if not SHA256_PATTERN.fullmatch(expected_validation_hash):
        raise ValueError("source member ledger recovery validation hash is invalid")
    if (
        not validation_path.is_file()
        or _sha256(validation_path) != expected_validation_hash
    ):
        raise ValueError("source member ledger recovery validation hash mismatch")
    try:
        validation = json.loads(validation_path.read_text(encoding="utf-8-sig"))
    except json.JSONDecodeError as error:
        raise ValueError("source member ledger recovery validation is invalid JSON") from error
    if not isinstance(validation, dict):
        raise ValueError("source member ledger recovery validation must be an object")
    if validation.get("schema") != RECOVERY_VALIDATION_SCHEMA:
        raise ValueError("source member ledger recovery validation schema mismatch")
    if validation.get("run_id") != source_manifest.get("run_id"):
        raise ValueError("source member ledger recovery validation run mismatch")
    archive_root = validation.get("archive_root")
    if not isinstance(archive_root, str) or Path(archive_root).resolve() != source_root:
        raise ValueError("source member ledger recovery validation root mismatch")
    if (
        validation.get("archive_manifest_sha256")
        != expected_manifest_sha256.lower()
    ):
        raise ValueError("source member ledger recovery validation manifest binding mismatch")
    if validation.get("validation") != REQUIRED_RECOVERY_FLAGS:
        raise ValueError("source member ledger recovery validation flags are not all true")

    ledger = validation.get("member_sha256")
    if not isinstance(ledger, dict) or not ledger:
        raise ValueError("source member ledger is missing")
    if not all(
        isinstance(path, str)
        and isinstance(digest, str)
        and SHA256_PATTERN.fullmatch(digest)
        for path, digest in ledger.items()
    ):
        raise ValueError("source member ledger contains an invalid entry")
    if validation.get("member_count") != len(ledger):
        raise ValueError("source member ledger count mismatch")
    actual = _source_member_hashes(source_root, validation_path)
    if set(actual) != set(ledger):
        missing = sorted(set(ledger) - set(actual))
        extra = sorted(set(actual) - set(ledger))
        raise ValueError(
            f"source member ledger set mismatch; missing={missing}, extra={extra}"
        )
    for relative, digest in ledger.items():
        if actual[relative] != digest:
            raise ValueError(f"source member ledger hash mismatch: {relative}")
    return validation, validation_path, expected_validation_hash


def _legacy_expected_matrix() -> tuple[tuple[str, float], ...]:
    return tuple(
        (topology, bias)
        for topology in REQUIRED_TOPOLOGIES
        for bias in REQUIRED_BIASES
    )


def _normalize_expected_matrix(
    expected_matrix: Sequence[Sequence[object]], *, label: str = "expected matrix"
) -> tuple[tuple[str, float], ...]:
    matrix: list[tuple[str, float]] = []
    for entry in expected_matrix:
        if isinstance(entry, (str, bytes)) or not isinstance(entry, Sequence) or len(entry) != 2:
            raise ValueError(f"{label} entry must be a topology/bias pair")
        topology = entry[0]
        if not isinstance(topology, str) or topology not in REQUIRED_TOPOLOGIES:
            raise ValueError(f"{label} contains invalid topology: {topology!r}")
        try:
            bias = float(entry[1])
        except (TypeError, ValueError) as error:
            raise ValueError(f"{label} contains an invalid bias") from error
        if not math.isfinite(bias):
            raise ValueError(f"{label} contains a non-finite bias")
        key = (topology, bias)
        if key in matrix:
            raise ValueError(f"{label} contains duplicate state {topology} at {bias:g} V")
        matrix.append(key)
    if not matrix:
        raise ValueError(f"{label} must not be empty")
    return tuple(matrix)


def _expected_matrix_document(
    expected_matrix: Sequence[Sequence[object]],
) -> list[list[object]]:
    return [[topology, bias] for topology, bias in expected_matrix]


def _manifest_expected_matrix(manifest: dict[str, object]) -> tuple[tuple[str, float], ...]:
    raw = manifest.get("expected_matrix")
    if not isinstance(raw, list):
        raise ValueError("prepared manifest is missing expected_matrix")
    expected = _normalize_expected_matrix(raw, label="manifest expected_matrix")
    canonical = tuple(sorted(expected))
    if raw != _expected_matrix_document(canonical):
        raise ValueError("manifest expected_matrix is not canonical")
    return canonical


def _prepared_manifest_expected_matrix(
    manifest: dict[str, object],
) -> tuple[tuple[str, float], ...]:
    expected = _manifest_expected_matrix(manifest)
    prepared: list[tuple[str, float]] = []
    states = manifest.get("states")
    if not isinstance(states, list):
        raise ValueError("prepared manifest states must be a list")
    for state in states:
        if not isinstance(state, dict):
            raise ValueError("prepared manifest state must be an object")
        topology = state.get("topology_id")
        if not isinstance(topology, str) or topology not in REQUIRED_TOPOLOGIES:
            raise ValueError(f"prepared state has invalid topology: {topology!r}")
        try:
            requested = float(state.get("requested_bias_V"))
        except (TypeError, ValueError) as error:
            raise ValueError("prepared state has invalid requested bias") from error
        if not math.isfinite(requested):
            raise ValueError("prepared state has non-finite requested bias")
        key = (topology, requested)
        if key in prepared:
            raise ValueError(f"prepared manifest has duplicate state {topology} at {requested:g} V")
        prepared.append(key)
    if set(prepared) != set(expected) or len(prepared) != len(expected):
        raise ValueError("prepared manifest states do not match expected_matrix")
    return expected


def _validate_sentaurus_version_provenance(manifest: dict[str, object]) -> None:
    version = manifest.get("sentaurus_version")
    if not isinstance(version, str) or not version:
        raise ValueError("completed manifest is missing sentaurus_version")
    for state in manifest.get("states", []):
        if not isinstance(state, dict) or state.get("sentaurus_version") != version:
            raise ValueError("completed manifest has mixed Sentaurus versions")


def validate_recovered_archive(
    root: Path, expected_manifest_sha256: str,
    expected_matrix: Sequence[Sequence[object]] | None = None,
) -> dict[str, object]:
    root = Path(root)
    manifest_path = root / "manifest.json"
    if not manifest_path.is_file():
        raise ValueError("recovered archive is missing manifest.json")
    if _sha256(manifest_path).lower() != expected_manifest_sha256.lower():
        raise ValueError("recovered archive manifest hash mismatch")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("outputs_complete") is not True:
        raise ValueError("recovered archive outputs are incomplete")
    if "expected_matrix" in manifest:
        recorded = _manifest_expected_matrix(manifest)
        if expected_matrix is not None and recorded != _normalize_expected_matrix(expected_matrix):
            raise ValueError("recovered archive expected_matrix does not match caller declaration")
        expected_matrix = recorded
    validate_state_matrix(manifest.get("states", []), expected_matrix)
    if "sentaurus_version" in manifest:
        _validate_sentaurus_version_provenance(manifest)
    return manifest


def _bias_tag(bias_V: float) -> str:
    sign = "m" if bias_V < 0.0 else "p"
    magnitude = format(abs(bias_V), ".17g").replace(".", "p")
    return f"{sign}{magnitude}V"


def validate_final_bias(requested_bias_V: float, actual_bias_V: float) -> float:
    requested = float(requested_bias_V)
    actual = float(actual_bias_V)
    if not math.isfinite(requested):
        raise ValueError(f"requested Anode contact voltage is not finite: {requested_bias_V}")
    if not math.isfinite(actual):
        raise ValueError(f"final Anode contact voltage is not finite: {actual_bias_V}")
    if abs(actual - requested) > BIAS_TOLERANCE_V:
        raise ValueError(
            f"final Anode contact voltage {actual:.17g} V does not match requested "
            f"{requested:.17g} V within 1e-12 V"
        )
    return actual


def validate_state_matrix(
    states: Sequence[dict[str, object]],
    expected_matrix: Sequence[Sequence[object]] | None = None,
) -> list[tuple[str, float]]:
    required = (
        _normalize_expected_matrix(expected_matrix)
        if expected_matrix is not None else _legacy_expected_matrix()
    )
    matrix: list[tuple[str, float]] = []
    for state in states:
        topology = str(state.get("topology_id", ""))
        requested = float(state.get("requested_bias_V", math.nan))
        if state.get("status") != "passed":
            raise ValueError(f"state {topology} at {requested:g} V is not passed")
        try:
            actual = validate_final_bias(requested, float(state.get("actual_bias_V", math.nan)))
        except ValueError as error:
            raise ValueError(f"{topology} at {requested:g} V: {error}") from error
        key = (topology, requested)
        if key in matrix:
            raise ValueError(f"duplicate state {topology} at {requested:g} V")
        matrix.append(key)
        if actual != actual:  # Defensive; validate_final_bias already rejects NaN.
            raise ValueError("unreachable non-finite bias")
    if set(matrix) != set(required) or len(matrix) != len(required):
        missing = sorted(set(required) - set(matrix))
        extra = sorted(set(matrix) - set(required))
        if expected_matrix is None:
            raise ValueError(f"exact six-state matrix mismatch; missing={missing}, extra={extra}")
        raise ValueError(
            f"exact declared state matrix mismatch; missing={missing}, extra={extra}"
        )
    return matrix


def validate_field_manifest(manifest: dict[str, object]) -> dict[str, dict[str, object]]:
    fields = manifest.get("fields")
    if not isinstance(fields, list):
        raise ValueError("field_manifest.json must contain a fields list")
    selected: dict[str, dict[str, object]] = {}
    for name, (components, unit) in _FIELD_CONTRACT.items():
        candidates = [field for field in fields if isinstance(field, dict) and field.get("name") == name]
        expected = {
            "region": 0,
            "components": components,
            "unit": unit,
            "mapping_status": "complete",
            "global_node_mapping": "global_vertex_order",
        }
        match = next(
            (field for field in candidates if all(field.get(key) == value for key, value in expected.items())),
            None,
        )
        if match is None:
            contract = ", ".join(f"{key}={value}" for key, value in expected.items())
            raise ValueError(f"{name} must satisfy {contract}; got {candidates}")
        selected[name] = match
    return selected


def validate_v2_field_manifest(
    manifest: dict[str, object],
) -> dict[str, dict[str, object]]:
    fields = manifest.get("fields")
    if not isinstance(fields, list):
        raise ValueError("field_manifest.json must contain a fields list")
    selected: dict[str, dict[str, object]] = {}
    raw_aliases: set[str] = set()
    for normalized, raw, components, unit, semantic_role in V2_FIELD_CONTRACT:
        if raw in raw_aliases:
            raise ValueError(f"duplicate v2 raw alias in contract: {raw}")
        raw_aliases.add(raw)
        candidates = [
            field for field in fields
            if isinstance(field, dict) and field.get("name") == raw
        ]
        expected = {
            "region": 0,
            "components": components,
            "unit": unit,
            "mapping_status": "complete",
            "global_node_mapping": "global_vertex_order",
        }
        matches = [
            candidate for candidate in candidates
            if all(candidate.get(key) == value for key, value in expected.items())
        ]
        if len(matches) > 1:
            raise ValueError(f"duplicate alias for required field {raw}")
        if len(matches) != 1:
            contract = ", ".join(f"{key}={value}" for key, value in expected.items())
            raise ValueError(f"{raw} must satisfy {contract}; got {candidates}")
        selected[normalized] = {
            "normalized_name": normalized,
            "raw_name": raw,
            "components": components,
            "unit": unit,
            "semantic_role": semantic_role,
        }
    return selected


def _schema_error(path: str, message: str) -> None:
    raise ValueError(f"v2 schema validation failed at {path}: {message}")


def _validate_json_schema(instance: object, schema: dict[str, object], path: str = "$") -> None:
    expected_type = schema.get("type")
    type_matches = {
        "object": isinstance(instance, dict),
        "array": isinstance(instance, list),
        "string": isinstance(instance, str),
        "number": isinstance(instance, (int, float)) and not isinstance(instance, bool),
        "integer": isinstance(instance, int) and not isinstance(instance, bool),
        "boolean": isinstance(instance, bool),
    }
    if expected_type in type_matches and not type_matches[expected_type]:
        _schema_error(path, f"expected {expected_type}")
    if expected_type == "number" and not math.isfinite(float(instance)):
        _schema_error(path, "number must be finite")
    if "const" in schema and instance != schema["const"]:
        _schema_error(path, f"expected constant {schema['const']!r}")
    if "enum" in schema and instance not in schema["enum"]:
        _schema_error(path, f"value {instance!r} is not in enum")
    if isinstance(instance, str):
        if len(instance) < int(schema.get("minLength", 0)):
            _schema_error(path, "string is too short")
        pattern = schema.get("pattern")
        if isinstance(pattern, str) and re.search(pattern, instance) is None:
            _schema_error(path, f"string does not match {pattern}")
    if isinstance(instance, list):
        if len(instance) < int(schema.get("minItems", 0)):
            _schema_error(path, "array has too few items")
        if "maxItems" in schema and len(instance) > int(schema["maxItems"]):
            _schema_error(path, "array has too many items")
        if schema.get("uniqueItems") is True:
            encoded = [json.dumps(item, sort_keys=True) for item in instance]
            if len(encoded) != len(set(encoded)):
                _schema_error(path, "array items must be unique")
        item_schema = schema.get("items")
        if isinstance(item_schema, dict):
            for index, item in enumerate(instance):
                _validate_json_schema(item, item_schema, f"{path}[{index}]")
    if isinstance(instance, dict):
        required = schema.get("required", [])
        if isinstance(required, list):
            missing = [name for name in required if name not in instance]
            if missing:
                _schema_error(path, f"missing required properties {missing}")
        properties = schema.get("properties", {})
        if not isinstance(properties, dict):
            _schema_error(path, "schema properties must be an object")
        for name, value in instance.items():
            child = properties.get(name)
            if isinstance(child, dict):
                _validate_json_schema(value, child, f"{path}.{name}")
            elif schema.get("additionalProperties") is False:
                _schema_error(path, f"unexpected property {name}")


def _load_v2_schema() -> dict[str, object]:
    try:
        schema = json.loads(V2_SCHEMA_PATH.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"cannot load v2 state schema: {error}") from error
    if (
        not isinstance(schema, dict)
        or schema.get("$schema") != "https://json-schema.org/draft/2020-12/schema"
        or schema.get("title") != V2_SCHEMA
    ):
        raise ValueError("unexpected v2 state schema document")
    return schema


def _expected_topology_contract(topology_id: str) -> dict[str, object]:
    topology = load_topology(TOPOLOGY_FIXTURE, topology_id)
    edges = {
        tuple(sorted(edge))
        for triangle in topology.triangles
        for edge in (
            (triangle[0], triangle[1]),
            (triangle[1], triangle[2]),
            (triangle[2], triangle[0]),
        )
    }
    return {
        "nodes": len(topology.nodes),
        "triangles": len(topology.triangles),
        "edges": len(edges),
        "contact_edges": {
            name: list(edge) for name, edge in topology.contacts.items()
        },
        "triangle_connectivity": [list(triangle) for triangle in topology.triangles],
    }


def _resolve_inside(root: Path, value: str, label: str) -> Path:
    path = Path(value)
    resolved = path.resolve() if path.is_absolute() else (root / path).resolve()
    if not resolved.is_relative_to(root.resolve()):
        raise ValueError(f"{label} escapes recovered archive root")
    return resolved


def _validate_recovery_states(
    root: Path, manifest: dict[str, object]
) -> dict[tuple[str, float], dict[str, object]]:
    validated: dict[tuple[str, float], dict[str, object]] = {}
    for state in manifest["states"]:
        topology = str(state["topology_id"])
        bias = float(state["requested_bias_V"])
        tag = _bias_tag(bias)
        if state.get("bias_tag") != tag:
            raise ValueError(f"{topology} at {bias:g} V has wrong bias tag")
        expected_topology = _expected_topology_contract(topology)
        if state.get("topology_contract") != expected_topology:
            raise ValueError(f"{topology} at {bias:g} V has wrong topology contract")
        export_dir = _resolve_inside(root, str(state.get("export_dir", "")), "export_dir")
        field_manifest = export_dir / "field_manifest.json"
        if not field_manifest.is_file():
            raise ValueError(f"{topology} at {bias:g} V is missing field_manifest.json")
        try:
            raw_fields = json.loads(field_manifest.read_text(encoding="utf-8"))
        except json.JSONDecodeError as error:
            raise ValueError(f"{topology} at {bias:g} V has invalid field manifest") from error
        fields = validate_v2_field_manifest(raw_fields)
        state_root = (root / "states" / topology / tag).resolve()
        if not state_root.is_dir() or not export_dir.is_relative_to(state_root):
            raise ValueError(f"{topology} at {bias:g} V has invalid state artifact root")
        members = collect_member_hashes(state_root)
        if not members:
            raise ValueError(f"{topology} at {bias:g} V has no raw artifacts")
        validated[(topology, bias)] = {
            "state_root": state_root,
            "export_dir": export_dir,
            "fields": fields,
            "member_sha256": members,
        }
    return validated


def _validate_v2_manifest_shape(manifest: dict[str, object]) -> None:
    _validate_json_schema(manifest, _load_v2_schema())
    if manifest["source_kind"] == "local_recovery" and (
        manifest["source_run_id"] != LOCAL_RECOVERY_RUN_ID
        or manifest["source_manifest_sha256"]
        != LOCAL_RECOVERY_MANIFEST_SHA256
    ):
        raise ValueError(
            "local recovery must use the named local recovery run "
            f"{LOCAL_RECOVERY_RUN_ID} with its prescribed manifest hash"
        )
    validate_state_matrix(manifest["states"])
    for state in manifest["states"]:
        topology = str(state["topology_id"])
        bias = float(state["requested_bias_V"])
        if state["bias_tag"] != _bias_tag(bias):
            raise ValueError(f"{topology} at {bias:g} V has wrong bias tag")
        if state["topology_contract"] != _expected_topology_contract(topology):
            raise ValueError(f"{topology} at {bias:g} V has wrong topology contract")
        state_root = Path("states") / topology / state["bias_tag"]
        canonical_paths = {
            "export_dir": (state_root / "export").as_posix(),
            "field_manifest": (
                state_root / "export" / "field_manifest.json"
            ).as_posix(),
            "state_csv": (state_root / "export" / "state.csv").as_posix(),
        }
        for name, expected in canonical_paths.items():
            if state.get(name) != expected:
                raise ValueError(
                    f"{topology} at {bias:g} V has non-canonical state path: {name}"
                )
        expected_fields = [
            {
                "normalized_name": normalized,
                "raw_name": raw,
                "components": components,
                "unit": unit,
                "semantic_role": role,
            }
            for normalized, raw, components, unit, role in V2_FIELD_CONTRACT
        ]
        if state["fields"] != expected_fields:
            raise ValueError(f"{topology} at {bias:g} V has wrong v2 field identities")


def _validate_audit_provenance(
    manifest: dict[str, object],
    artifacts_by_state: dict[tuple[str, float], dict[str, str]],
) -> None:
    provenance = manifest.get("audit_provenance")
    if not isinstance(provenance, dict):
        raise ValueError("audit provenance is required")
    producer = provenance.get("producer")
    producer_hash = provenance.get("producer_sha256")
    if producer != AUDIT_PRODUCER or not isinstance(producer_hash, str):
        raise ValueError("audit provenance has an invalid executable reference")
    producer_path = REPO / producer
    if (
        not SHA256_PATTERN.fullmatch(producer_hash)
        or not producer_path.is_file()
        or _sha256(producer_path) != producer_hash
    ):
        raise ValueError("audit provenance executable hash mismatch")
    source_commit = provenance.get("task4_source_commit")
    if not isinstance(source_commit, str) or re.fullmatch(
        r"[0-9a-f]{40}", source_commit
    ) is None:
        raise ValueError("audit provenance has an invalid source commit")
    if not isinstance(provenance.get("replay_environment"), str) or not str(
        provenance["replay_environment"]
    ):
        raise ValueError("audit provenance has no replay environment")
    replays = provenance.get("replays")
    if not isinstance(replays, list) or len(replays) != 6:
        raise ValueError("audit provenance must contain exactly six replays")

    seen: set[tuple[str, float]] = set()
    for replay in replays:
        if not isinstance(replay, dict):
            raise ValueError("audit provenance replay must be an object")
        try:
            key = (str(replay["topology_id"]), float(replay["bias_V"]))
        except (KeyError, TypeError, ValueError) as error:
            raise ValueError("audit provenance replay identity is invalid") from error
        if key not in artifacts_by_state or key in seen:
            raise ValueError("audit provenance replay matrix is not exact")
        seen.add(key)
        if replay.get("producer") != producer or replay.get("exit_code") != 0:
            raise ValueError("audit provenance replay execution is invalid")
        arguments = replay.get("arguments")
        if not isinstance(arguments, list) or len(arguments) != 2 * len(
            AUDIT_PATH_OPTIONS
        ):
            raise ValueError("audit provenance arguments must be a 14-item array")
        state_root = (
            Path("states") / key[0] / _bias_tag(key[1]) / "export"
        )
        argument_paths: list[str] = []
        canonical_paths: list[str] = []
        for index, (option, filename) in enumerate(AUDIT_PATH_OPTIONS):
            value = arguments[2 * index + 1]
            if arguments[2 * index] != option or not isinstance(value, str):
                raise ValueError("audit provenance argument order is invalid")
            canonical = (state_root / filename).as_posix()
            normalized = value.replace("\\", "/")
            if normalized != canonical and not normalized.endswith("/" + canonical):
                raise ValueError("audit provenance path is not canonical for its state")
            if canonical not in artifacts_by_state[key]:
                raise ValueError("audit provenance path is absent from artifact ledger")
            argument_paths.append(value)
            canonical_paths.append(canonical)
        if replay.get("command") != " ".join([producer, *arguments]):
            raise ValueError("audit provenance command does not match arguments")
        for label, start, stop in (
            ("input", 0, AUDIT_INPUT_COUNT),
            ("output", AUDIT_INPUT_COUNT, len(AUDIT_PATH_OPTIONS)),
        ):
            hashes = replay.get(f"{label}_sha256")
            expected_keys = set(argument_paths[start:stop])
            if not isinstance(hashes, dict) or set(hashes) != expected_keys:
                raise ValueError(
                    f"audit provenance {label} hash references are not exact"
                )
            for index in range(start, stop):
                source_path = argument_paths[index]
                canonical = canonical_paths[index]
                if hashes[source_path] != artifacts_by_state[key][canonical]:
                    raise ValueError(
                        f"audit provenance {label} hash mismatch: {canonical}"
                    )
    if seen != set(artifacts_by_state):
        raise ValueError("audit provenance replay matrix is not exact")


def seal_recovered_archive(
    source_root: Path,
    sealed_root: Path,
    expected_manifest_sha256: str,
    *,
    run_id: str,
    source_kind: str = "local_recovery",
    expected_source_validation_sha256: str | None = None,
) -> dict[str, object]:
    source_root = Path(source_root).resolve()
    sealed_root = Path(sealed_root).resolve()
    if source_root.parent != sealed_root.parent:
        raise ValueError("sealed v2 archive must be a sibling of the source archive")
    if not _SAFE_RUN_ID.fullmatch(run_id) or sealed_root.name != run_id:
        raise ValueError("sealed v2 run ID is invalid or does not match its directory")
    if sealed_root.exists():
        raise FileExistsError(f"sealed v2 archive already exists: {sealed_root}")
    source_manifest = validate_recovered_archive(source_root, expected_manifest_sha256)
    if source_kind == "local_recovery":
        if (
            source_manifest.get("run_id") != LOCAL_RECOVERY_RUN_ID
            or expected_manifest_sha256.lower()
            != LOCAL_RECOVERY_MANIFEST_SHA256
        ):
            raise ValueError(
                "local recovery must use the named local recovery run "
                f"{LOCAL_RECOVERY_RUN_ID} with its prescribed manifest hash"
            )
    elif source_kind != "regenerated":
        raise ValueError(f"unsupported source kind: {source_kind}")
    if not isinstance(source_manifest.get("task4_provenance"), dict):
        raise ValueError("audit provenance is required in the source manifest")

    if source_manifest.get("schema") != SCHEMA:
        raise ValueError("recovered archive must preserve the v1 schema")
    validated_states = _validate_recovery_states(source_root, source_manifest)
    source_validation, source_validation_path, source_validation_hash = (
        validate_source_recovery(
            source_root, source_manifest, expected_manifest_sha256,
            source_kind=source_kind,
            expected_source_validation_sha256=expected_source_validation_sha256,
        )
    )

    sealed_root.mkdir(parents=True)
    shutil.copyfile(source_root / "manifest.json", sealed_root / "source_manifest.json")
    shutil.copytree(source_root / "states", sealed_root / "states")
    shutil.copyfile(
        source_validation_path,
        sealed_root / SEALED_SOURCE_VALIDATION_PATH,
    )
    states: list[dict[str, object]] = []
    for source_state in source_manifest["states"]:
        topology = str(source_state["topology_id"])
        bias = float(source_state["requested_bias_V"])
        tag = _bias_tag(bias)
        state_root = sealed_root / "states" / topology / tag
        export_dir = state_root / "export"
        fields = list(validated_states[(topology, bias)]["fields"].values())
        raw_artifacts = [
            {"path": path, "sha256": digest}
            for path, digest in collect_member_hashes(state_root).items()
        ]
        for artifact in raw_artifacts:
            artifact["path"] = (
                Path("states") / topology / tag / artifact["path"]
            ).as_posix()
        state = {
            "topology_id": topology,
            "requested_bias_V": bias,
            "actual_bias_V": float(source_state["actual_bias_V"]),
            "bias_tag": tag,
            "status": "passed",
            "topology_contract": source_state["topology_contract"],
            "export_dir": export_dir.relative_to(sealed_root).as_posix(),
            "field_manifest": (
                export_dir / "field_manifest.json"
            ).relative_to(sealed_root).as_posix(),
            "fields": fields,
            "raw_artifacts": raw_artifacts,
        }
        state_csv = export_dir / "state.csv"
        if state_csv.is_file():
            state["state_csv"] = state_csv.relative_to(sealed_root).as_posix()
        states.append(state)
    manifest: dict[str, object] = {
        "schema": V2_SCHEMA,
        "run_id": run_id,
        "source_kind": source_kind,
        "source_schema": SCHEMA,
        "source_run_id": str(source_manifest.get("run_id", "")),
        "source_manifest_path": "source_manifest.json",
        "source_manifest_sha256": expected_manifest_sha256.lower(),
        "bias_tolerance_V": BIAS_TOLERANCE_V,
        "outputs_complete": True,
        "states": states,
        "source_validation": {
            "kind": source_kind,
            "schema": RECOVERY_VALIDATION_SCHEMA,
            "path": SEALED_SOURCE_VALIDATION_PATH,
            "sha256": source_validation_hash,
            "archive_root": str(source_validation["archive_root"]),
            "archive_manifest_sha256": str(
                source_validation["archive_manifest_sha256"]
            ),
            "member_count": int(source_validation["member_count"]),
            "validation": source_validation["validation"],
        },
    }
    manifest["audit_provenance"] = source_manifest["task4_provenance"]
    _validate_v2_manifest_shape(manifest)
    write_manifest(sealed_root / "manifest.json", manifest)
    return validate_sealed_archive(sealed_root)


def _validate_sealed_source_recovery(
    root: Path, manifest: dict[str, object]
) -> None:
    identity = manifest["source_validation"]
    if identity["kind"] != manifest["source_kind"]:
        raise ValueError("sealed source member ledger kind mismatch")
    if (
        identity["archive_manifest_sha256"]
        != manifest["source_manifest_sha256"]
    ):
        raise ValueError("sealed source member ledger manifest binding mismatch")
    if manifest["source_kind"] == "local_recovery":
        required_root = str((DEFAULT_OUTPUT_DIR / LOCAL_RECOVERY_RUN_ID).resolve())
        if (
            identity["sha256"] != LOCAL_RECOVERY_VALIDATION_SHA256
            or Path(str(identity["archive_root"])).resolve()
            != Path(required_root).resolve()
        ):
            raise ValueError("sealed local recovery validation identity mismatch")
    validation_path = _resolve_inside(
        root, str(identity["path"]), "source recovery validation path"
    )
    if (
        not validation_path.is_file()
        or _sha256(validation_path) != identity["sha256"]
    ):
        raise ValueError("sealed source member ledger validation hash mismatch")
    try:
        validation = json.loads(validation_path.read_text(encoding="utf-8-sig"))
    except json.JSONDecodeError as error:
        raise ValueError("sealed source member ledger validation is invalid") from error
    summary = {
        "schema": validation.get("schema"),
        "archive_root": validation.get("archive_root"),
        "archive_manifest_sha256": validation.get("archive_manifest_sha256"),
        "member_count": validation.get("member_count"),
        "validation": validation.get("validation"),
    }
    expected_summary = {
        "schema": identity["schema"],
        "archive_root": identity["archive_root"],
        "archive_manifest_sha256": identity["archive_manifest_sha256"],
        "member_count": identity["member_count"],
        "validation": identity["validation"],
    }
    if summary != expected_summary or validation.get("run_id") != manifest["source_run_id"]:
        raise ValueError("sealed source member ledger identity mismatch")
    if identity["validation"] != REQUIRED_RECOVERY_FLAGS:
        raise ValueError("sealed source member ledger validation flags are not all true")
    ledger = validation.get("member_sha256")
    if not isinstance(ledger, dict) or identity["member_count"] != len(ledger):
        raise ValueError("sealed source member ledger count mismatch")
    source_manifest_path = _resolve_inside(
        root, str(manifest["source_manifest_path"]), "source_manifest_path"
    )
    actual = {"manifest.json": _sha256(source_manifest_path)}
    states_root = root / "states"
    for path in sorted(states_root.rglob("*")):
        if path.is_file():
            actual[path.relative_to(root).as_posix()] = _sha256(path)
    if set(actual) != set(ledger):
        missing = sorted(set(ledger) - set(actual))
        extra = sorted(set(actual) - set(ledger))
        raise ValueError(
            f"sealed source member ledger set mismatch; missing={missing}, extra={extra}"
        )
    for relative, digest in ledger.items():
        if not isinstance(digest, str) or actual[relative] != digest:
            raise ValueError(f"sealed source member ledger hash mismatch: {relative}")


def validate_sealed_archive(root: Path) -> dict[str, object]:
    root = Path(root).resolve()
    manifest_path = root / "manifest.json"
    if not manifest_path.is_file():
        raise ValueError("sealed archive is missing manifest.json")
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as error:
        raise ValueError("sealed archive manifest is invalid JSON") from error
    _validate_v2_manifest_shape(manifest)

    source_path = _resolve_inside(
        root, str(manifest["source_manifest_path"]), "source_manifest_path"
    )
    if not source_path.is_file() or _sha256(source_path) != manifest["source_manifest_sha256"]:
        raise ValueError("source manifest hash mismatch")
    source_manifest = json.loads(source_path.read_text(encoding="utf-8"))
    if (
        source_manifest.get("schema") != SCHEMA
        or source_manifest.get("run_id") != manifest["source_run_id"]
        or source_manifest.get("outputs_complete") is not True
    ):
        raise ValueError("source manifest identity mismatch")
    if source_manifest.get("task4_provenance") != manifest["audit_provenance"]:
        raise ValueError(
            "audit provenance does not match the source manifest"
        )

    _validate_sealed_source_recovery(root, manifest)
    artifacts_by_state: dict[tuple[str, float], dict[str, str]] = {}
    for state in manifest["states"]:
        topology = str(state["topology_id"])
        bias = float(state["requested_bias_V"])
        tag = str(state["bias_tag"])
        state_root = (root / "states" / topology / tag).resolve()
        artifact_paths = [
            str(artifact["path"]) for artifact in state["raw_artifacts"]
        ]
        if len(artifact_paths) != len(set(artifact_paths)):
            raise ValueError(
                f"duplicate artifact path: {topology} at {bias:g} V"
            )
        recorded = {
            str(artifact["path"]): str(artifact["sha256"])
            for artifact in state["raw_artifacts"]
        }
        actual = {
            (state_root / relative).relative_to(root).as_posix(): digest
            for relative, digest in collect_member_hashes(state_root).items()
        }
        if set(recorded) != set(actual):
            raise ValueError(
                f"raw artifact member set mismatch: {topology} at {bias:g} V"
            )
        required_members = {
            str(state["field_manifest"]),
            str(state["state_csv"]),
        }
        if not required_members.issubset(recorded):
            raise ValueError(
                f"canonical state path is absent from ledger: {topology} at {bias:g} V"
            )
        export_prefix = str(state["export_dir"]) + "/"
        if not any(path.startswith(export_prefix) for path in recorded):
            raise ValueError(
                "canonical export_dir does not contain the artifact ledger"
            )
        artifacts_by_state[(topology, bias)] = recorded
        for relative, digest in recorded.items():
            path = _resolve_inside(root, relative, "raw artifact path")
            if not path.is_relative_to(state_root) or actual[relative] != digest:
                raise ValueError(f"raw artifact hash mismatch: {relative}")
        field_manifest = _resolve_inside(
            root, str(state["field_manifest"]), "field_manifest"
        )
        identities = list(validate_v2_field_manifest(
            json.loads(field_manifest.read_text(encoding="utf-8"))
        ).values())
        if identities != state["fields"]:
            raise ValueError(
                f"v2 field identity mismatch: {topology} at {bias:g} V"
            )
    _validate_audit_provenance(manifest, artifacts_by_state)
    return manifest



def canonical_minimal6_coordinates() -> dict[int, tuple[float, float]]:
    topology = load_topology(TOPOLOGY_FIXTURE, "sketch")
    return {
        label - 1: tuple(topology.nodes[label])
        for label in sorted(topology.nodes)
    }


def _read_csv(path: Path, required: set[str]) -> list[dict[str, str]]:
    if not path.is_file():
        raise ValueError(f"missing neutral export file: {path.name}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        columns = set(reader.fieldnames or [])
        if not required.issubset(columns):
            raise ValueError(f"{path.name} missing columns: {sorted(required - columns)}")
        return list(reader)


def _canonical_source_ids(
    export_dir: Path, coordinates: dict[int, tuple[float, float]]
) -> dict[int, int]:
    rows = _read_csv(export_dir / "nodes.csv", {"id", "x_um", "y_um"})
    if len(rows) != len(coordinates):
        raise ValueError(f"expected {len(coordinates)} canonical nodes, found {len(rows)}")
    source_for_canonical: dict[int, int] = {}
    used_sources: set[int] = set()
    for canonical, expected in coordinates.items():
        matches = [
            int(row["id"])
            for row in rows
            if abs(float(row["x_um"]) - expected[0]) < COORDINATE_TOLERANCE_UM
            and abs(float(row["y_um"]) - expected[1]) < COORDINATE_TOLERANCE_UM
        ]
        if len(matches) != 1:
            raise ValueError(
                f"canonical node {canonical} requires one exact coordinate mapping; got {matches}"
            )
        if matches[0] in used_sources:
            raise ValueError("canonical node mapping reuses a source node")
        source_for_canonical[canonical] = matches[0]
        used_sources.add(matches[0])
    return source_for_canonical


def _read_scalar_field(export_dir: Path, name: str) -> dict[int, float]:
    path = export_dir / "fields" / f"{name}_region0.csv"
    if not path.is_file():
        role = {"eQuasiFermiPotential": "phin", "hQuasiFermiPotential": "phip"}.get(name)
        suffix = f" ({role})" if role else ""
        raise ValueError(f"missing required field {name}{suffix}; density-derived QF is forbidden")
    rows = _read_csv(path, {"node_id", "component0"})
    values: dict[int, float] = {}
    for row in rows:
        node_id = int(row["node_id"])
        value = float(row["component0"])
        if node_id in values or not math.isfinite(value):
            raise ValueError(f"{name} has duplicate or non-finite node {node_id}")
        values[node_id] = value
    return values


def write_state_csv(
    export_dir: Path, coordinates: dict[int, tuple[float, float]] | None = None
) -> Path:
    export_dir = Path(export_dir)
    manifest_path = export_dir / "field_manifest.json"
    if not manifest_path.is_file():
        raise ValueError("missing neutral export file: field_manifest.json")
    validate_field_manifest(json.loads(manifest_path.read_text(encoding="utf-8")))
    raw = {
        "psi_V": _read_scalar_field(export_dir, "ElectrostaticPotential"),
        "phin_V": _read_scalar_field(export_dir, "eQuasiFermiPotential"),
        "phip_V": _read_scalar_field(export_dir, "hQuasiFermiPotential"),
        "n_m3": _read_scalar_field(export_dir, "eDensity"),
        "p_m3": _read_scalar_field(export_dir, "hDensity"),
    }
    source_ids = _canonical_source_ids(export_dir, coordinates or canonical_minimal6_coordinates())
    output = export_dir / "state.csv"
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["node_id", "psi_V", "phin_V", "phip_V", "n_m3", "p_m3"])
        for canonical in sorted(source_ids):
            source = source_ids[canonical]
            try:
                row = [raw[name][source] for name in ("psi_V", "phin_V", "phip_V", "n_m3", "p_m3")]
            except KeyError as error:
                raise ValueError(f"state field mapping is missing source node {source}") from error
            row[3] *= 1.0e6
            row[4] *= 1.0e6
            writer.writerow([canonical, *(format(value, ".17g") for value in row)])
    return output


def _render_deck(target_bias_V: float, tag: str, destination: Path) -> None:
    source = SOURCE_DECK.read_text(encoding="utf-8")
    rendered = source.replace("__TARGET_BIAS_V__", format(target_bias_V, ".17g"))
    rendered = rendered.replace("__BIAS_TAG__", tag)
    if "__TARGET_BIAS_V__" in rendered or "__BIAS_TAG__" in rendered:
        raise ValueError("state deck placeholder replacement is incomplete")
    destination.write_text(rendered, encoding="utf-8", newline="\n")


def _validated_requested_matrix(
    topology_ids: Sequence[str], biases: Sequence[float]
) -> tuple[tuple[str, ...], tuple[float, ...]]:
    topologies = tuple(topology_ids)
    if not topologies or any(
        not isinstance(topology, str) or topology not in REQUIRED_TOPOLOGIES
        for topology in topologies
    ):
        raise ValueError(f"topologies must be nonempty members of {REQUIRED_TOPOLOGIES}")
    if len(topologies) != len(set(topologies)):
        raise ValueError("topologies must not contain duplicates")
    try:
        bias_values = tuple(float(value) for value in biases)
    except (TypeError, ValueError) as error:
        raise ValueError("biases must be numeric") from error
    if not bias_values or not all(math.isfinite(bias) for bias in bias_values):
        raise ValueError("biases must be nonempty and finite")
    if len(bias_values) != len(set(bias_values)):
        raise ValueError("biases must not contain duplicates")
    return topologies, bias_values


def prepare_exports(
    *, topology_ids: Sequence[str], biases: Sequence[float], run_id: str,
    output_dir: Path, ssh_target: str, remote_root: str = DEFAULT_REMOTE_ROOT,
    importer: Path = DEFAULT_IMPORTER,
) -> dict[str, object]:
    topologies, bias_values = _validated_requested_matrix(topology_ids, biases)
    if not _SAFE_RUN_ID.fullmatch(run_id):
        raise ValueError("run ID contains unsupported characters")
    if not _SAFE_REMOTE_COMPONENT.fullmatch(remote_root):
        raise ValueError("remote root contains unsupported characters")
    if not SOURCE_DECK.is_file() or not MODELS_SOURCE.is_file():
        raise FileNotFoundError("missing minimal6 state deck or model parameters")
    run_root = (Path(output_dir) / run_id).resolve()
    states: list[dict[str, object]] = []
    for topology_id in topologies:
        for bias in bias_values:
            tag = _bias_tag(bias)
            state_root = run_root / "states" / topology_id / tag
            bundle = state_root / "source"
            artifacts = state_root / "artifacts"
            neutral = state_root / "export"
            gate_bundle = build_gate_bundle(topology_id, bundle)
            gate_deck = bundle / "pn2d_minimal6_gate_sdevice.cmd"
            gate_deck.unlink()
            deck_name = f"pn2d_minimal6_state_{tag}_sdevice.cmd"
            _render_deck(bias, tag, bundle / deck_name)
            staged = ["pn2d_minimal6.grd", "pn2d_minimal6.dat", "models.par", deck_name]
            remote_dir = f"{remote_root.rstrip('/')}/{run_id}/{topology_id}/{tag}"
            plot_stem = f"pn2d_minimal6_state_{tag}"
            final_tdr_name = f"{plot_stem}.tdr"
            current_plt_name = f"{plot_stem}.plt"
            log_name = f"{plot_stem}_des.log"
            stdout_name = f"run_{plot_stem}.out"
            returned_files = [
                final_tdr_name,
                current_plt_name,
                log_name,
                "pn2d_minimal6.tdr",
                "pn2d_minimal6.grd",
                "pn2d_minimal6.dat",
                "run_tdx_dfise_to_tdr.out",
                stdout_name,
            ]
            states.append({
                "topology_id": topology_id,
                "requested_bias_V": bias,
                "bias_tag": tag,
                "bundle_dir": str(bundle),
                "artifacts_dir": str(artifacts),
                "export_dir": str(neutral),
                "remote_dir": remote_dir,
                "deck_name": deck_name,
                "staged_files": staged,
                "topology_contract": gate_bundle["topology_contract"],
                "remote_commands": [
                    f"cd {remote_dir} && tdx -d pn2d_minimal6.grd pn2d_minimal6.dat "
                    "pn2d_minimal6.tdr > run_tdx_dfise_to_tdr.out 2>&1",
                    f"cd {remote_dir} && sdevice {deck_name} > {stdout_name} 2>&1",
                ],
                "final_tdr_name": final_tdr_name,
                "current_plt_name": current_plt_name,
                "log_name": log_name,
                "stdout_name": stdout_name,
                "returned_files": returned_files,
                "status": "prepared",
            })
    expected_matrix = tuple(sorted(
        (topology, bias) for topology in topologies for bias in bias_values
    ))
    manifest: dict[str, object] = {
        "schema": SCHEMA,
        "run_id": run_id,
        "ssh_target": ssh_target,
        "remote_root": remote_root,
        "importer": str(Path(importer).resolve()),
        "bias_tolerance_V": BIAS_TOLERANCE_V,
        "expected_matrix": _expected_matrix_document(expected_matrix),
        "sentaurus_version": None,
        "outputs_complete": False,
        "states": states,
        "manifest_path": str(run_root / "manifest.json"),
    }
    write_manifest(Path(str(manifest["manifest_path"])), manifest)
    return manifest


def _parse_final_anode_bias(plt_path: Path) -> float:
    text = plt_path.read_text(errors="replace")
    datasets = parse_quoted_list(text, "datasets")
    if "Anode OuterVoltage" not in datasets:
        raise ValueError(f"{plt_path.name} lacks Anode OuterVoltage")
    rows = parse_values_block(text, len(datasets))
    if not rows:
        raise ValueError(f"{plt_path.name} has no bias rows")
    return float(rows[-1][datasets.index("Anode OuterVoltage")])


def _live_executor(
    state: dict[str, object], *, ssh_bin: str, scp_bin: str,
    ssh_target: str, importer: Path,
) -> dict[str, object]:
    if not importer.is_file():
        raise FileNotFoundError(f"Sentaurus importer is not built: {importer}")
    bundle = Path(str(state["bundle_dir"]))
    artifacts = Path(str(state["artifacts_dir"]))
    neutral = Path(str(state["export_dir"]))
    artifacts.mkdir(parents=True, exist_ok=True)
    remote_dir = str(state["remote_dir"])
    returned = [str(name) for name in state["returned_files"]]

    def return_argv(name: str) -> list[str]:
        return [
            scp_bin,
            f"{ssh_target}:{remote_dir}/{name}",
            str(artifacts) + os.sep,
        ]

    try:
        run_checked([ssh_bin, ssh_target, f"mkdir -p {remote_dir}"])
        for name in state["staged_files"]:
            run_checked([scp_bin, str(bundle / str(name)), f"{ssh_target}:{remote_dir}/"])
        for command in state["remote_commands"]:
            run_checked([ssh_bin, ssh_target, str(command)])
    except Exception:
        recovery_errors: list[str] = []
        for name in returned:
            try:
                run_checked(return_argv(name))
            except Exception as recovery_error:
                recovery_errors.append(f"{name}: {recovery_error}")
        state["artifact_recovery_errors"] = recovery_errors
        raise
    for name in returned:
        run_checked(return_argv(name))
    final_tdr = artifacts / str(state["final_tdr_name"])
    # The current C++ CLI performs the export-neutral operation when --export-dir is present.
    run_checked([
        str(importer), "--tdr", str(final_tdr), "--export-dir", str(neutral),
        "--compensated-doping-policy", "reported",
    ])
    version = subprocess.run(
        [ssh_bin, ssh_target, "sdevice -version"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    if not version:
        raise ValueError("remote Sentaurus version command produced no output")
    return {
        "actual_bias_V": _parse_final_anode_bias(artifacts / str(state["current_plt_name"])),
        "export_dir": str(neutral),
        "sentaurus_version": version,
    }


def run_exports(
    manifest: dict[str, object],
    *, executor: Callable[[dict[str, object]], dict[str, object] | None],
) -> None:
    expected_matrix = _prepared_manifest_expected_matrix(manifest)
    expected_document = _expected_matrix_document(expected_matrix)
    manifest_path = Path(str(manifest["manifest_path"]))
    manifest["outputs_complete"] = False
    write_manifest(manifest_path, manifest)
    for state in manifest["states"]:
        try:
            result = executor(state) or {}
            if "actual_bias_V" not in result:
                raise ValueError("executor result is missing actual_bias_V")
            try:
                actual_bias = float(result["actual_bias_V"])
            except (TypeError, ValueError) as error:
                raise ValueError(
                    f"executor result has invalid actual_bias_V: {result['actual_bias_V']!r}"
                ) from error
            state["actual_bias_V"] = validate_final_bias(
                float(state["requested_bias_V"]), actual_bias
            )
            export_dir = Path(str(result.get("export_dir", state["export_dir"])))
            if not export_dir.is_dir():
                raise ValueError(f"missing neutral export directory: {export_dir}")
            version = result.get("sentaurus_version")
            if not isinstance(version, str) or not version:
                raise ValueError("executor result is missing sentaurus_version")
            recorded_version = manifest.get("sentaurus_version")
            if recorded_version is not None and recorded_version != version:
                raise ValueError(
                    f"mixed Sentaurus versions: {recorded_version!r} and {version!r}"
                )
            manifest["sentaurus_version"] = version
            state["sentaurus_version"] = version
            state["state_csv"] = str(write_state_csv(export_dir))
            state["field_manifest"] = str(export_dir / "field_manifest.json")
            state["member_sha256"] = collect_member_hashes(export_dir)
            state["status"] = "passed"
            if manifest.get("expected_matrix") != expected_document:
                raise ValueError("manifest expected_matrix changed after preparation")
            write_manifest(manifest_path, manifest)
        except Exception as error:
            manifest["expected_matrix"] = expected_document
            state["status"] = "failed"
            state["error"] = str(error)
            manifest["error"] = str(error)
            write_manifest(manifest_path, manifest)
            raise
    validate_state_matrix(manifest["states"], expected_matrix)
    if manifest.get("expected_matrix") != expected_document:
        raise ValueError("manifest expected_matrix changed after preparation")
    _validate_sentaurus_version_provenance(manifest)
    manifest["outputs_complete"] = True
    manifest.pop("error", None)
    write_manifest(manifest_path, manifest)


def _parse_csv_values(raw: str) -> tuple[float, ...]:
    values = tuple(float(value.strip()) for value in raw.split(",") if value.strip())
    if not values:
        raise ValueError("bias list is empty")
    return values


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topologies", default="sketch,mirror")
    parser.add_argument("--biases", default="0,-12,-19")
    parser.add_argument("--run-id", default=None)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--ssh-target", default="sentaurus")
    parser.add_argument("--ssh-bin", default=None)
    parser.add_argument("--scp-bin", default=None)
    parser.add_argument("--remote-root", default=DEFAULT_REMOTE_ROOT)
    parser.add_argument("--importer", type=Path, default=DEFAULT_IMPORTER)
    parser.add_argument("--dry-run", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    manifest: dict[str, object] | None = None
    try:
        run_id = args.run_id or datetime.now().strftime("minimal6_states_%Y%m%d_%H%M%S")
        manifest = prepare_exports(
            topology_ids=tuple(value.strip() for value in args.topologies.split(",") if value.strip()),
            biases=_parse_csv_values(args.biases),
            run_id=run_id,
            output_dir=args.output_dir,
            ssh_target=args.ssh_target,
            remote_root=args.remote_root,
            importer=args.importer,
        )
        if not args.dry_run:
            ssh_bin = args.ssh_bin or default_windows_openssh("ssh")
            scp_bin = args.scp_bin or default_windows_openssh("scp")
            run_exports(
                manifest,
                executor=lambda state: _live_executor(
                    state, ssh_bin=ssh_bin, scp_bin=scp_bin,
                    ssh_target=args.ssh_target, importer=args.importer.resolve(),
                ),
            )
        print(json.dumps(manifest, indent=2))
        return 0
    except Exception as error:  # noqa: BLE001 - partial manifest is the contract.
        if manifest is not None:
            manifest["outputs_complete"] = False
            manifest["error"] = str(error)
            write_manifest(Path(str(manifest["manifest_path"])), manifest)
        print(str(error), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
