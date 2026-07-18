"""Hash-bound, diagnostic-only canonical inputs for the Minimal6 inverse audit."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path, PureWindowsPath
from typing import Any

try:
    from .inverse_contracts import Observation, SampleStatus, SupportKind, classify_numeric_sample
except ImportError:
    from scripts.pn2d_minimal6_diagnostics.inverse_contracts import (  # type: ignore
        Observation, SampleStatus, SupportKind, classify_numeric_sample,
    )


COMMON_BIASES = tuple(float(-value) for value in range(1, 21))
TOPOLOGIES = ("sketch", "mirror")
DISCOVERY_KEYS = tuple(("sketch", value) for value in (-1.0, -4.0, -8.0,
                                                       -12.0, -16.0, -19.0, -20.0))
REQUIRED_SENTARUS_FIELDS = {
    "ElectrostaticPotential": (1, "V"), "eQuasiFermiPotential": (1, "V"),
    "hQuasiFermiPotential": (1, "V"), "eDensity": (1, "cm^-3"),
    "hDensity": (1, "cm^-3"), "ElectricField": (2, "V*cm^-1"),
    "eCurrentDensity": (2, "A*cm^-2"), "hCurrentDensity": (2, "A*cm^-2"),
    "eAlphaAvalanche": (1, "cm^-1"), "hAlphaAvalanche": (1, "cm^-1"),
    "ImpactIonization": (1, "cm^-3*s^-1"),
}
SUPPLEMENTAL_FIELDS = {
    "eMobility": (1, "cm^2*V^-1*s^-1"), "hMobility": (1, "cm^2*V^-1*s^-1"),
    "eVelocity": (1, "cm*s^-1"), "hVelocity": (1, "cm*s^-1"),
}

_SHA256 = re.compile(r"[0-9a-f]{64}")
_EXPECTED_KEYS = {(topology, bias) for topology in TOPOLOGIES for bias in COMMON_BIASES}
_CONVERSIONS = {
    "V": (1.0, "V", "V_to_V"), "cm^-3": (1.0e6, "m^-3", "cm^-3_to_m^-3"),
    "V*cm^-1": (100.0, "V/m", "V*cm^-1_to_V/m"),
    "A*cm^-2": (1.0e4, "A/m^2", "A*cm^-2_to_A/m^2"),
    "cm^2*V^-1*s^-1": (1.0e-4, "m^2*V^-1*s^-1", "cm^2*V^-1*s^-1_to_m^2*V^-1*s^-1"),
    "cm^-1": (100.0, "m^-1", "cm^-1_to_m^-1"),
    "cm*s^-1": (1.0e-2, "m*s^-1", "cm*s^-1_to_m*s^-1"),
    "cm^-3*s^-1": (1.0e6, "m^-3*s^-1", "cm^-3*s^-1_to_m^-3*s^-1"),
}


@dataclass(frozen=True)
class InputBundle:
    common_keys: tuple[tuple[str, float], ...]
    discovery_keys: tuple[tuple[str, float], ...]
    holdout_keys: tuple[tuple[str, float], ...]
    observations: tuple[Observation, ...]
    fields: tuple[tuple[str, int, str], ...]
    executable_hashes: tuple[str, ...]
    tracked_source_hashes: tuple[tuple[str, str], ...]
    input_hashes: tuple[tuple[str, str], ...]


@dataclass(frozen=True)
class _State:
    topology: str
    bias: float
    path: Path
    digest: str
    coordinate_frame: str
    orientation: str
    rows: tuple[dict[str, str], ...]


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _require_sha256(value: Any, label: str) -> str:
    if not isinstance(value, str) or not _SHA256.fullmatch(value.lower()):
        raise ValueError(f"{label} must be a SHA-256")
    return value.lower()


def _checked_path(root: Path, value: Any, label: str) -> Path:
    if not isinstance(value, str) or not value or Path(value).is_absolute() or PureWindowsPath(value).is_absolute():
        raise ValueError(f"{label} escapes root")
    candidate = (root / value).resolve()
    try:
        candidate.relative_to(root)
    except ValueError as error:
        raise ValueError(f"{label} escapes root") from error
    if not candidate.is_file():
        raise ValueError(f"{label} is not a file")
    return candidate


def _load_manifest(root_value: str | Path, solver: str) -> tuple[Path, dict[str, Any], dict[str, str]]:
    root = Path(root_value).resolve()
    if not root.is_dir():
        raise ValueError(f"{solver} root is not a directory")
    manifest_path = _checked_path(root, "manifest.json", f"{solver} manifest")
    seal_path = _checked_path(root, "seal.json", f"{solver} seal")
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        seal = json.loads(seal_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"{solver} manifest schema mismatch") from error
    if not isinstance(manifest, dict) or not isinstance(seal, dict):
        raise ValueError(f"{solver} manifest schema mismatch")
    expected = _require_sha256(seal.get("manifest_sha256"), f"{solver} seal manifest_sha256")
    if _sha256(manifest_path) != expected:
        raise ValueError(f"{solver} manifest hash mismatch")
    if manifest.get("schema") != "vela.pn2d_minimal6_inverse_input.v1" or manifest.get("solver") != solver:
        raise ValueError(f"{solver} manifest schema mismatch")
    ledger = manifest.get("member_sha256")
    if not isinstance(ledger, dict) or not ledger:
        raise ValueError(f"{solver} member hash ledger mismatch")
    checked: dict[str, str] = {}
    for relative, digest in ledger.items():
        path = _checked_path(root, relative, f"{solver} member")
        expected_digest = _require_sha256(digest, f"{solver} member hash")
        if _sha256(path) != expected_digest:
            raise ValueError(f"{solver} state hash mismatch")
        checked[relative] = expected_digest
    return root, manifest, checked


def _field_contract(manifest: dict[str, Any], expected: dict[str, tuple[int, str]], label: str) -> tuple[tuple[str, int, str], ...]:
    items = manifest.get("fields")
    if not isinstance(items, list):
        raise ValueError(f"{label} field schema mismatch")
    actual: dict[str, tuple[int, str]] = {}
    for item in items:
        if not isinstance(item, dict) or not isinstance(item.get("name"), str):
            raise ValueError(f"{label} field schema mismatch")
        name = item["name"]
        if name in actual or item.get("support_kind") != SupportKind.NODE.value:
            raise ValueError(f"{label} field support mismatch")
        components, unit = item.get("components"), item.get("unit")
        if not isinstance(components, int) or not isinstance(unit, str):
            raise ValueError(f"{label} field schema mismatch")
        actual[name] = (components, unit)
    if set(actual) != set(expected):
        raise ValueError(f"{label} field mismatch")
    for name, contract in expected.items():
        components, unit = actual[name]
        if components != contract[0]:
            raise ValueError(f"{label} field component mismatch")
        if unit != contract[1]:
            raise ValueError(f"{label} field unit mismatch")
    return tuple((name, *expected[name]) for name in sorted(expected))


def _states(root: Path, manifest: dict[str, Any], checked: dict[str, str], label: str) -> dict[tuple[str, float], _State]:
    tolerance = manifest.get("bias_tolerance_V")
    if not isinstance(tolerance, (int, float)) or not math.isfinite(tolerance) or tolerance > 1.0e-12 or tolerance < 0.0:
        raise ValueError(f"{label} bias tolerance mismatch")
    rows = manifest.get("states")
    if not isinstance(rows, list):
        raise ValueError(f"{label} state schema mismatch")
    result: dict[tuple[str, float], _State] = {}
    for item in rows:
        if not isinstance(item, dict):
            raise ValueError(f"{label} state schema mismatch")
        topology, requested, actual = item.get("topology"), item.get("requested_bias_V"), item.get("actual_bias_V")
        if topology not in TOPOLOGIES:
            raise ValueError(f"{label} checkpoint matrix mismatch")
        if not isinstance(requested, (int, float)) or not isinstance(actual, (int, float)):
            raise ValueError(f"{label} state schema mismatch")
        requested, actual = float(requested), float(actual)
        if not math.isfinite(requested) or not math.isfinite(actual) or abs(actual - requested) > 1.0e-12:
            raise ValueError(f"{label} exact bias mismatch")
        key = (topology, requested)
        if key not in _EXPECTED_KEYS:
            raise ValueError(f"{label} checkpoint mismatch")
        if key in result:
            raise ValueError(f"{label} duplicate checkpoint")
        if item.get("support_kind") != SupportKind.NODE.value:
            raise ValueError(f"{label} support mismatch")
        coordinate_frame, orientation = item.get("coordinate_frame"), item.get("orientation")
        if not isinstance(coordinate_frame, str) or not coordinate_frame or not isinstance(orientation, str) or not orientation:
            raise ValueError(f"{label} provenance mismatch")
        relative = item.get("state_path")
        path = _checked_path(root, relative, f"{label} state")
        digest = _require_sha256(item.get("state_sha256"), f"{label} state hash")
        if relative not in checked or checked[relative] != digest or _sha256(path) != digest:
            raise ValueError(f"{label} state hash mismatch")
        try:
            with path.open(newline="", encoding="utf-8") as handle:
                data = tuple(csv.DictReader(handle))
        except (OSError, csv.Error) as error:
            raise ValueError(f"{label} state CSV mismatch") from error
        if not data or not all(row.get("canonical_node_id") and row.get("x_um") and row.get("y_um") for row in data):
            raise ValueError(f"{label} canonical node mismatch")
        node_ids = [row["canonical_node_id"] for row in data]
        if len(node_ids) != len(set(node_ids)):
            raise ValueError(f"{label} duplicate canonical node")
        result[key] = _State(topology, requested, path, digest, coordinate_frame, orientation, data)
    if set(result) != _EXPECTED_KEYS:
        raise ValueError(f"{label} checkpoint matrix mismatch")
    return result


def _provenance(manifest: dict[str, Any], label: str) -> tuple[tuple[str, ...], tuple[tuple[str, str], ...]]:
    provenance = manifest.get("provenance")
    if not isinstance(provenance, dict):
        raise ValueError(f"{label} provenance mismatch")
    executable = _require_sha256(provenance.get("executable_sha256"), f"{label} executable hash")
    source = provenance.get("tracked_source_sha256")
    if not isinstance(source, dict) or not source:
        raise ValueError(f"{label} tracked-source provenance mismatch")
    sources = tuple(sorted((str(path), _require_sha256(value, f"{label} tracked source hash")) for path, value in source.items()))
    return (executable,), sources


def _value(row: dict[str, str], header: str) -> float | None:
    value = row.get(header)
    if value is None or value == "":
        return None
    try:
        return float(value)
    except ValueError as error:
        raise ValueError(f"non-numeric field value for {header}") from error


def _validate_headers(states: dict[tuple[str, float], _State], fields: dict[str, tuple[int, str]], label: str) -> None:
    required = {"canonical_node_id", "x_um", "y_um"}
    required.update(f"{name}_component{component}" for name, (components, _) in fields.items() for component in range(components))
    for state in states.values():
        missing = required - set(state.rows[0])
        if missing:
            raise ValueError(f"{label} field component header mismatch")


def _observations(states: dict[tuple[str, float], _State], fields: dict[str, tuple[int, str]], *, solver: str) -> tuple[Observation, ...]:
    observations: list[Observation] = []
    for key in sorted(states):
        state = states[key]
        for row in state.rows:
            node_id = row["canonical_node_id"]
            if not node_id.isdecimal() or str(int(node_id)) != node_id:
                raise ValueError("canonical node identifier mismatch")
            for component, raw_value in (("x", _value(row, "x_um")), ("y", _value(row, "y_um"))):
                if raw_value is None or not math.isfinite(raw_value):
                    raise ValueError("coordinate mismatch")
                observations.append(Observation(
                    solver=solver, topology=state.topology, bias_V=state.bias,
                    support_kind=SupportKind.NODE, support_id=node_id,
                    quantity="coordinate", component=component, raw_value=raw_value,
                    raw_unit="um", value_si=raw_value * 1.0e-6, unit_si="m",
                    coordinate_frame=state.coordinate_frame, orientation=state.orientation,
                    conversion="um_to_m", status=SampleStatus.VALID,
                    source_path=state.path.as_posix(), source_sha256=state.digest,
                ))
            for quantity, (components, raw_unit) in fields.items():
                scale, unit_si, conversion = _CONVERSIONS[raw_unit]
                for component_index in range(components):
                    raw_value = _value(row, f"{quantity}_component{component_index}")
                    observations.append(Observation(
                        solver=solver, topology=state.topology, bias_V=state.bias,
                        support_kind=SupportKind.NODE, support_id=row["canonical_node_id"],
                        quantity=quantity, component=f"component{component_index}", raw_value=raw_value,
                        raw_unit=raw_unit, value_si=None if raw_value is None else raw_value * scale,
                        unit_si=unit_si, coordinate_frame=state.coordinate_frame,
                        orientation=state.orientation, conversion=conversion,
                        status=classify_numeric_sample(raw_value, floor=0.0),
                        source_path=state.path.as_posix(), source_sha256=state.digest,
                    ))
    return tuple(observations)


def load_input_bundle(vela_root: str | Path, sentaurus_root: str | Path, supplemental_root: str | Path) -> InputBundle:
    vela_base, vela_manifest, vela_hashes = _load_manifest(vela_root, "vela")
    sentaurus_base, sentaurus_manifest, sentaurus_hashes = _load_manifest(sentaurus_root, "sentaurus")
    supplemental_base, supplemental_manifest, supplemental_hashes = _load_manifest(supplemental_root, "supplemental")
    vela_fields = _field_contract(vela_manifest, REQUIRED_SENTARUS_FIELDS, "vela")
    sentaurus_fields = _field_contract(sentaurus_manifest, REQUIRED_SENTARUS_FIELDS, "sentaurus")
    supplemental_fields = _field_contract(supplemental_manifest, SUPPLEMENTAL_FIELDS, "supplemental")
    vela_states = _states(vela_base, vela_manifest, vela_hashes, "vela")
    sentaurus_states = _states(sentaurus_base, sentaurus_manifest, sentaurus_hashes, "sentaurus")
    supplemental_states = _states(supplemental_base, supplemental_manifest, supplemental_hashes, "supplemental")
    _validate_headers(vela_states, REQUIRED_SENTARUS_FIELDS, "vela")
    _validate_headers(sentaurus_states, REQUIRED_SENTARUS_FIELDS, "sentaurus")
    _validate_headers(supplemental_states, SUPPLEMENTAL_FIELDS, "supplemental")
    if set(vela_states) != set(sentaurus_states):
        raise ValueError("Vela/Sentaurus checkpoint matrix mismatch")
    if set(supplemental_states) != set(sentaurus_states):
        raise ValueError("supplemental checkpoint matrix mismatch")
    for key, sentaurus_state in sentaurus_states.items():
        if {row["canonical_node_id"] for row in supplemental_states[key].rows} != {row["canonical_node_id"] for row in sentaurus_state.rows}:
            raise ValueError("supplemental canonical-node mismatch")
    executable: set[str] = set()
    sources: set[tuple[str, str]] = set()
    for manifest, label in ((vela_manifest, "vela"), (sentaurus_manifest, "sentaurus"), (supplemental_manifest, "supplemental")):
        current_executable, current_sources = _provenance(manifest, label)
        executable.update(current_executable)
        sources.update((f"{label}:{path}", digest) for path, digest in current_sources)
    common = tuple(sorted(_EXPECTED_KEYS))
    observations = (_observations(vela_states, REQUIRED_SENTARUS_FIELDS, solver="vela")
                    + _observations(sentaurus_states, REQUIRED_SENTARUS_FIELDS, solver="sentaurus")
                    + _observations(supplemental_states, SUPPLEMENTAL_FIELDS, solver="sentaurus"))
    return InputBundle(common, DISCOVERY_KEYS, tuple(key for key in common if key not in DISCOVERY_KEYS), observations,
                       tuple(sorted(set(vela_fields + sentaurus_fields + supplemental_fields))), tuple(sorted(executable)),
                       tuple(sorted(sources)), tuple(sorted(
                           [(f"vela:{path}", digest) for path, digest in vela_hashes.items()] +
                           [(f"sentaurus:{path}", digest) for path, digest in sentaurus_hashes.items()] +
                           [(f"supplemental:{path}", digest) for path, digest in supplemental_hashes.items()])))


def canonical_observations(bundle: InputBundle) -> tuple[Observation, ...]:
    return bundle.observations


def field_inventory(bundle: InputBundle) -> dict[str, dict[str, object]]:
    return {name: {"components": components, "unit": unit, "support_kind": SupportKind.NODE.value}
            for name, components, unit in bundle.fields}


def write_input_manifest(bundle: InputBundle, path: str | Path) -> None:
    target = Path(path)
    target.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "schema": "vela.pn2d_minimal6_inverse_inputs.v1",
        "common_keys": [[topology, bias] for topology, bias in bundle.common_keys],
        "discovery_keys": [[topology, bias] for topology, bias in bundle.discovery_keys],
        "holdout_keys": [[topology, bias] for topology, bias in bundle.holdout_keys],
        "field_inventory": field_inventory(bundle), "executable_sha256": list(bundle.executable_hashes),
        "tracked_source_sha256": {item[0]: item[1] for item in bundle.tracked_source_hashes},
        "input_sha256": {item[0]: item[1] for item in bundle.input_hashes},
    }
    target.write_text(json.dumps(payload, sort_keys=True, indent=2) + "\n", encoding="utf-8")
