"""Fail-closed canonical sealer for live Minimal6 inverse-audit inputs."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import os
from pathlib import Path, PureWindowsPath
import re
import shutil
import subprocess
import tempfile
from typing import Callable, Iterable

try:
    from .inverse_inputs import (
        COMMON_BIASES, REQUIRED_SENTARUS_FIELDS, SUPPLEMENTAL_FIELDS,
        TOPOLOGIES, load_input_bundle,
    )
except ImportError:
    from scripts.pn2d_minimal6_diagnostics.inverse_inputs import (  # type: ignore
        COMMON_BIASES, REQUIRED_SENTARUS_FIELDS, SUPPLEMENTAL_FIELDS,
        TOPOLOGIES, load_input_bundle,
    )


BIAS_TOLERANCE_V = 1.0e-12
COORDINATE_TOLERANCE_UM = 1.0e-12
ORIENTATION = "+x,+y;canonical_mesh_node_ids"
FRAME = "minimal6_cartesian"
SUPPLEMENTAL_VERSION = "O-2018.06-SP2"
_PHASE_BASE = re.compile(r"[0-9a-fA-F]{7,40}")
_ImporterRunner = Callable[[Path, Path], None]
_REPO = Path(__file__).resolve().parents[2]
_TRACKED_SOURCE_PATHS = (
    "include/vela/discretization/ScharfetterGummel.h",
    "include/vela/physics/ImpactIonizationModel.h",
    "include/vela/equation/AssemblerUtils.h",
    "src/equation/CoupledDDAssembler.cpp",
    "src/simulation/DCSweep.cpp",
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _read_json(path: Path, label: str) -> dict:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"{label} is not valid JSON") from error
    if not isinstance(value, dict):
        raise ValueError(f"{label} must be an object")
    return value


def _write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, sort_keys=True, indent=2) + "\n", encoding="utf-8")


def _root(path: str | Path, label: str) -> Path:
    result = Path(path).resolve()
    if not result.is_dir():
        raise ValueError(f"{label} is not a directory")
    return result


def _checked(root: Path, value: object, label: str) -> Path:
    if (not isinstance(value, str) or not value or Path(value).is_absolute()
            or PureWindowsPath(value).is_absolute()):
        raise ValueError(f"{label} escapes root")
    path = (root / value).resolve()
    try:
        path.relative_to(root)
    except ValueError as error:
        raise ValueError(f"{label} escapes root") from error
    if not path.is_file():
        raise ValueError(f"{label} is not a file")
    return path


def _require_hash(path: Path, expected: object, label: str) -> None:
    if not isinstance(expected, str) or not re.fullmatch(r"[0-9a-fA-F]{64}", expected):
        raise ValueError(f"{label} hash is invalid")
    if _sha256(path) != expected.lower():
        raise ValueError(f"{label} hash mismatch")


def _finite(value: object, label: str) -> float:
    if not isinstance(value, (int, float)) or not math.isfinite(float(value)):
        raise ValueError(f"{label} must be finite")
    return float(value)


def _bias_tag(bias: float) -> str:
    magnitude = int(round(abs(bias)))
    if bias != -float(magnitude) or magnitude not in range(1, 21):
        raise ValueError("bias tag requires an exact -1..-20 V bias")
    return f"m{magnitude}V"


def _matrix(include_zero: bool) -> set[tuple[str, float]]:
    biases = (0.0,) + COMMON_BIASES if include_zero else COMMON_BIASES
    return {(topology, bias) for topology in TOPOLOGIES for bias in biases}


def _validate_sweep(root: Path, solver: str, *, include_zero: bool) -> tuple[dict, dict, dict]:
    manifest_path = root / "sweep_manifest.json"
    manifest = _read_json(manifest_path, f"{solver} sweep manifest")
    if manifest.get("schema") != "vela.pn2d_minimal6_sweep_manifest.v1":
        raise ValueError(f"{solver} sweep schema mismatch")
    if manifest.get("interpolation") != "forbidden":
        raise ValueError(f"{solver} interpolation must be forbidden")
    if manifest.get("failed_transition") is not None or manifest.get("failed_transitions") not in (None, []):
        raise ValueError(f"{solver} sweep contains failed transitions")
    expected = _matrix(include_zero)
    accepted = manifest.get("accepted_checkpoints")
    if not isinstance(accepted, list):
        raise ValueError(f"{solver} accepted checkpoint schema mismatch")
    states = {}
    for row in accepted:
        if not isinstance(row, dict) or row.get("solver") != solver or row.get("status") != "accepted":
            raise ValueError(f"{solver} checkpoint status mismatch")
        topology = row.get("topology")
        requested = _finite(row.get("target_bias_V"), f"{solver} target bias")
        actual = _finite(row.get("actual_bias_V"), f"{solver} actual bias")
        if topology not in TOPOLOGIES or abs(requested - actual) > BIAS_TOLERANCE_V:
            raise ValueError(f"{solver} exact bias mismatch")
        key = (topology, requested)
        if key not in expected or key in states:
            raise ValueError(f"{solver} checkpoint matrix mismatch")
        state = _checked(root, row.get("state_path"), f"{solver} state")
        _require_hash(state, row.get("state_sha256"), f"{solver} state")
        states[key] = (row, state)
    if set(states) != expected:
        raise ValueError(f"{solver} checkpoint matrix mismatch")
    segment_name = "segments" if solver == "vela" else "sentaurus_segments"
    segments = manifest.get(segment_name)
    if not isinstance(segments, list):
        raise ValueError(f"{solver} segment schema mismatch")
    decks = {}
    for row in segments:
        if not isinstance(row, dict) or row.get("solver") != solver or row.get("status") != "accepted":
            raise ValueError(f"{solver} segment status mismatch")
        topology = row.get("topology")
        bias = _finite(row.get("target_bias_V"), f"{solver} segment bias")
        key = (topology, bias)
        if key not in expected or key in decks:
            raise ValueError(f"{solver} segment matrix mismatch")
        deck = _checked(root, row.get("deck"), f"{solver} deck")
        _require_hash(deck, row.get("deck_sha256"), f"{solver} deck")
        decks[key] = deck
    if set(decks) != expected:
        raise ValueError(f"{solver} segment matrix mismatch")
    return manifest, states, decks


def _mesh_contract(path: Path) -> dict:
    value = _read_json(path, "Vela topology mesh")
    if value.get("coordinate_unit") != "um":
        raise ValueError("Vela topology mesh coordinate unit mismatch")
    nodes = value.get("nodes")
    if not isinstance(nodes, list) or len(nodes) != 6:
        raise ValueError("Vela topology mesh node mismatch")
    coordinates = {}
    for row in nodes:
        if not isinstance(row, dict) or not isinstance(row.get("id"), int):
            raise ValueError("Vela topology mesh node mismatch")
        node = row["id"]
        if node in coordinates:
            raise ValueError("Vela topology mesh duplicate node")
        coordinates[node] = (_finite(row.get("x"), "mesh x"), _finite(row.get("y"), "mesh y"))
    if set(coordinates) != set(range(6)):
        raise ValueError("Vela topology mesh canonical node mismatch")
    triangles = value.get("triangles")
    if not isinstance(triangles, list) or len(triangles) != 4:
        raise ValueError("Vela topology mesh triangle mismatch")
    triangle_ids = set()
    topology = []
    for row in triangles:
        ids = row.get("node_ids") if isinstance(row, dict) else None
        element_id = row.get("id") if isinstance(row, dict) else None
        if not isinstance(element_id, int) or element_id in triangle_ids:
            raise ValueError("Vela topology mesh duplicate triangle id")
        triangle_ids.add(element_id)
        if (not isinstance(ids, list) or len(ids) != 3 or len(set(ids)) != 3
                or any(node not in coordinates for node in ids)):
            raise ValueError("Vela topology mesh triangle mismatch")
        triangle = tuple(sorted(ids))
        if triangle in topology:
            raise ValueError("Vela topology mesh duplicate triangle")
        topology.append(triangle)
    if len(topology) != 4:
        raise ValueError("Vela topology mesh duplicate triangle")
    contacts = value.get("contacts")
    if not isinstance(contacts, list) or len(contacts) != 2:
        raise ValueError("Vela topology mesh contact mismatch")
    contact_ids = set()
    contact_names = set()
    contact_map = {}
    for row in contacts:
        if not isinstance(row, dict) or not isinstance(row.get("id"), int):
            raise ValueError("Vela topology mesh contact mismatch")
        contact_id = row["id"]
        name = row.get("name")
        node_ids = row.get("node_ids")
        if contact_id in contact_ids:
            raise ValueError("Vela topology mesh duplicate contact id")
        if not isinstance(name, str) or name in contact_names:
            raise ValueError("Vela topology mesh duplicate contact name")
        if (not isinstance(node_ids, list) or len(node_ids) != 2
                or len(set(node_ids)) != 2
                or any(node not in coordinates for node in node_ids)):
            raise ValueError("Vela topology mesh duplicate contact node")
        contact_ids.add(contact_id)
        contact_names.add(name)
        contact_map[name] = set(node_ids)
    if contact_map != {"Anode": {0, 4}, "Cathode": {2, 3}}:
        raise ValueError("Vela topology mesh contact mismatch")
    return {"coordinates": coordinates, "triangles": set(topology), "contacts": contact_map}


def _read_vela_state(path: Path) -> dict[int, dict[str, float]]:
    expected = ["node_id", "psi", "phin", "phip", "electrons_m3", "holes_m3"]
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != expected:
            raise ValueError("Vela state CSV header mismatch")
        rows = {}
        for row in reader:
            try:
                node = int(row["node_id"])
            except ValueError as error:
                raise ValueError("Vela state node mismatch") from error
            if str(node) != row["node_id"] or node in rows:
                raise ValueError("Vela state duplicate node")
            values = {name: _finite(float(row[name]), f"Vela state {name}") for name in expected[1:]}
            rows[node] = values
    if set(rows) != set(range(6)):
        raise ValueError("Vela state canonical node mismatch")
    return rows


def _validate_supplemental(root: Path) -> tuple[dict, dict]:
    manifest = _read_json(root / "manifest.json", "supplemental manifest")
    if manifest.get("schema") != "vela.pn2d_minimal6_states.v1":
        raise ValueError("supplemental schema mismatch")
    if manifest.get("outputs_complete") is not True:
        raise ValueError("supplemental outputs are incomplete")
    if manifest.get("sentaurus_version") != SUPPLEMENTAL_VERSION:
        raise ValueError("supplemental Sentaurus version mismatch")
    raw_matrix = manifest.get("expected_matrix")
    if type(raw_matrix) is not list or len(raw_matrix) != 40:
        raise ValueError("supplemental expected matrix mismatch")
    declared = set()
    for row in raw_matrix:
        if type(row) is not list or len(row) != 2:
            raise ValueError("supplemental expected matrix mismatch")
        topology, bias = row
        if (type(topology) is not str
                or type(bias) is not float
                or not math.isfinite(bias)):
            raise ValueError("supplemental expected matrix mismatch")
        key = (topology, bias)
        if key in declared:
            raise ValueError("supplemental expected matrix mismatch")
        declared.add(key)
    if declared != _matrix(False):
        raise ValueError("supplemental expected matrix mismatch")
    rows = manifest.get("states")
    if not isinstance(rows, list):
        raise ValueError("supplemental state schema mismatch")
    result = {}
    for row in rows:
        if not isinstance(row, dict) or row.get("status") != "passed":
            raise ValueError("supplemental state status mismatch")
        topology = row.get("topology_id")
        requested = _finite(row.get("requested_bias_V"), "supplemental requested bias")
        actual = _finite(row.get("actual_bias_V"), "supplemental actual bias")
        key = (topology, requested)
        if key not in _matrix(False) or key in result or abs(requested - actual) > BIAS_TOLERANCE_V:
            raise ValueError("supplemental exact matrix mismatch")
        tag = _bias_tag(requested)
        if row.get("bias_tag") != tag or row.get("sentaurus_version") != SUPPLEMENTAL_VERSION:
            raise ValueError("supplemental state provenance mismatch")
        state_root = (root / "states" / topology / tag).resolve()
        try:
            state_root.relative_to(root)
        except ValueError as error:
            raise ValueError("supplemental state root escapes run root") from error
        export = state_root / "export"
        artifacts = state_root / "artifacts"
        if not export.is_dir() or not artifacts.is_dir():
            raise ValueError("supplemental state artifacts are missing")
        ledger = row.get("member_sha256")
        if not isinstance(ledger, dict) or not ledger:
            raise ValueError("supplemental member hash ledger mismatch")
        actual_members = {path.relative_to(export).as_posix() for path in export.rglob("*") if path.is_file()}
        if set(ledger) != actual_members:
            raise ValueError("supplemental member hash ledger is incomplete")
        for relative, digest in ledger.items():
            member = _checked(export, relative, "supplemental export member")
            _require_hash(member, digest, "supplemental export member")
        tdr = _checked(artifacts, row.get("final_tdr_name"), "supplemental final TDR")
        result[key] = (row, tdr, dict(ledger))
    if set(result) != _matrix(False):
        raise ValueError("supplemental exact matrix mismatch")
    return manifest, result


def _read_rows(path: Path, required: Iterable[str], label: str) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fields = set(reader.fieldnames or [])
        missing = set(required) - fields
        if missing:
            raise ValueError(f"{label} missing columns {sorted(missing)}")
        return list(reader)


def _validate_import(export: Path, contract: dict, mesh_contract: dict, bias: float) -> dict:
    manifest = _read_json(export / "field_manifest.json", "imported field manifest")
    entries = manifest.get("fields")
    if not isinstance(entries, list):
        raise ValueError("imported field schema mismatch")
    by_name = {}
    contact_entries = []
    for entry in entries:
        if not isinstance(entry, dict):
            continue
        name = entry.get("name")
        if name == "ContactExternalVoltage":
            contact_entries.append(entry)
            continue
        if entry.get("region") != 0:
            continue
        name = entry.get("name")
        if name in contract:
            if name in by_name:
                raise ValueError("imported duplicate field")
            components, unit = contract[name]
            if entry.get("components") != components:
                raise ValueError(f"imported {name} component mismatch")
            if entry.get("unit") != unit:
                raise ValueError(f"imported {name} unit mismatch")
            if entry.get("mapping_status") != "complete" or entry.get("global_vertex_count") != 6:
                raise ValueError(f"imported {name} mapping mismatch")
            by_name[name] = entry
    if set(by_name) != set(contract):
        raise ValueError("imported required field mismatch")
    node_rows = _read_rows(export / "nodes.csv", {"id", "x_um", "y_um"}, "imported nodes")
    if len(node_rows) != 6:
        raise ValueError("imported node count mismatch")
    source_coordinates = {}
    for row in node_rows:
        node = int(row["id"])
        if node in source_coordinates:
            raise ValueError("imported duplicate node")
        source_coordinates[node] = (_finite(float(row["x_um"]), "imported x"),
                                    _finite(float(row["y_um"]), "imported y"))
    source_for_canonical = {}
    used = set()
    for canonical, expected in mesh_contract["coordinates"].items():
        matches = [node for node, xy in source_coordinates.items()
                   if abs(xy[0] - expected[0]) <= COORDINATE_TOLERANCE_UM
                   and abs(xy[1] - expected[1]) <= COORDINATE_TOLERANCE_UM]
        if len(matches) != 1 or matches[0] in used:
            raise ValueError("imported canonical coordinate mismatch")
        source_for_canonical[canonical] = matches[0]
        used.add(matches[0])
    canonical_for_source = {source: canonical for canonical, source in source_for_canonical.items()}
    element_rows = _read_rows(
        export / "elements.csv", {"id", "node0", "node1", "node2"}, "imported elements")
    if len(element_rows) != 4:
        raise ValueError("imported topology element count mismatch")
    element_ids = [int(row["id"]) for row in element_rows]
    if len(set(element_ids)) != len(element_ids):
        raise ValueError("imported duplicate element id")
    source_triangles = [tuple(int(row[f"node{i}"]) for i in range(3)) for row in element_rows]
    if any(len(set(nodes)) != 3 for nodes in source_triangles):
        raise ValueError("imported topology duplicate element node")
    try:
        imported_triangles = {tuple(sorted(canonical_for_source[int(row[f"node{i}"])] for i in range(3)))
                              for row in element_rows}
    except KeyError as error:
        raise ValueError("imported topology references an unknown node") from error
    if imported_triangles != mesh_contract["triangles"]:
        raise ValueError("imported topology mismatch")
    contact_rows = _read_rows(
        export / "contacts.csv", {"name", "node_ids", "region"}, "imported contacts")
    if len(contact_rows) != 2:
        raise ValueError("imported contact topology count mismatch")
    contacts = {}
    for row in contact_rows:
        name = row["name"]
        if name in contacts:
            raise ValueError("imported duplicate contact name")
        if row["region"] != "R.Si":
            raise ValueError("imported contact region mismatch")
        source_node_ids = [value for value in row["node_ids"].split(";") if value]
        if len(source_node_ids) != 2 or len(set(source_node_ids)) != 2:
            raise ValueError("imported duplicate contact node")
        try:
            contacts[name] = {canonical_for_source[int(value)] for value in source_node_ids}
        except (KeyError, ValueError) as error:
            raise ValueError("imported contact references an unknown node") from error
    if contacts != mesh_contract["contacts"]:
        raise ValueError("imported contact topology mismatch")
    expected_contact_regions = {1: "Cathode", 2: "Anode"}
    if len(contact_entries) != 2:
        raise ValueError("imported contact voltage manifest count mismatch")
    entries_by_region = {}
    for entry in contact_entries:
        region = entry.get("region")
        if region not in expected_contact_regions or region in entries_by_region:
            raise ValueError("imported contact voltage manifest region mismatch")
        if (entry.get("components") != 1
                or entry.get("unit") != "V"
                or entry.get("mapping_status") != "scalar"
                or entry.get("region_name") != expected_contact_regions[region]
                or entry.get("global_node_mapping") != "contact_scalar"
                or entry.get("values") != 1
                or entry.get("global_vertex_count") != 6):
            raise ValueError("imported contact voltage manifest contract mismatch")
        entries_by_region[region] = entry
    expected_paths = {export / "fields" / f"ContactExternalVoltage_region{region}.csv"
                      for region in expected_contact_regions}
    actual_paths = set((export / "fields").glob("ContactExternalVoltage_region*.csv"))
    if actual_paths != expected_paths:
        raise ValueError("imported contact voltage file set mismatch")
    for region, role in expected_contact_regions.items():
        path = export / "fields" / f"ContactExternalVoltage_region{region}.csv"
        rows = _read_rows(path, {"node_id", "component0"}, "imported contact voltage")
        if len(rows) != 1:
            raise ValueError("imported contact voltage row mismatch")
        try:
            source = int(rows[0]["node_id"])
            canonical = canonical_for_source[source]
        except (KeyError, ValueError) as error:
            raise ValueError("imported contact voltage node mismatch") from error
        if canonical not in contacts[role]:
            raise ValueError("imported contact voltage role mismatch")
        value = _finite(float(rows[0]["component0"]), "imported contact voltage")
        expected = bias if role == "Anode" else 0.0
        if abs(value - expected) > BIAS_TOLERANCE_V:
            raise ValueError("imported contact bias role mismatch")
    fields = {}
    for name, (components, _) in contract.items():
        rows = _read_rows(export / "fields" / f"{name}_region0.csv",
                          {"node_id", *(f"component{i}" for i in range(components))},
                          f"imported {name}")
        values = {}
        for row in rows:
            source = int(row["node_id"])
            if source in values:
                raise ValueError(f"imported {name} duplicate node")
            values[source] = tuple(_finite(float(row[f"component{i}"]), f"imported {name}")
                                   for i in range(components))
        if set(values) != set(source_coordinates):
            raise ValueError(f"imported {name} node mismatch")
        fields[name] = {canonical: values[source] for canonical, source in source_for_canonical.items()}
    return {"fields": fields, "coordinates": mesh_contract["coordinates"]}


def _copy(source: Path, root: Path, relative: str) -> None:
    target = root / relative
    target.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, target)


def _fmt(value: float | None) -> str:
    return "" if value is None else format(value, ".17g")


def _write_state(path: Path, coordinates: dict, fields: dict, contract: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    header = ["canonical_node_id", "x_um", "y_um"]
    for name, (components, _) in contract.items():
        header.extend(f"{name}_component{component}" for component in range(components))
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(header)
        for node in sorted(coordinates):
            row = [str(node), _fmt(coordinates[node][0]), _fmt(coordinates[node][1])]
            for name, (components, _) in contract.items():
                values = fields[name][node]
                if len(values) != components:
                    raise ValueError("canonical field component mismatch")
                row.extend(_fmt(value) for value in values)
            writer.writerow(row)


def _seal_root(root: Path, solver: str, contract: dict, records: list[dict],
               executable: Path, source_manifest: Path, phase_base: str,
               extra_provenance: dict,
               expected_executable_sha256: str | None = None) -> None:
    executable_sha256 = _sha256(executable)
    if expected_executable_sha256 is not None and executable_sha256 != expected_executable_sha256:
        raise ValueError("sealed executable SHA mismatch")
    _copy(source_manifest, root, "source/source_manifest.json")
    _copy(executable, root, f"source/executables/{executable.name}")
    _copy(Path(__file__), root, "source/tooling/inverse_input_sealer.py")
    for relative in _TRACKED_SOURCE_PATHS:
        source = _REPO / relative
        if not source.is_file():
            raise ValueError(f"tracked source is missing: {relative}")
        _copy(source, root, f"source/tracked/{relative}")
    states = []
    for record in records:
        relative = f"states/{record['topology']}/{_bias_tag(record['bias'])}.csv"
        path = root / relative
        _write_state(path, record["coordinates"], record["fields"], contract)
        for source, evidence_relative in record["evidence"]:
            _copy(source, root, evidence_relative)
        states.append({
            "topology": record["topology"], "requested_bias_V": record["bias"],
            "actual_bias_V": record["bias"], "state_path": relative,
            "state_sha256": _sha256(path), "support_kind": "node",
            "coordinate_frame": FRAME, "orientation": ORIENTATION,
        })
    ledger = {path.relative_to(root).as_posix(): _sha256(path)
              for path in sorted(root.rglob("*")) if path.is_file()}
    executable_member = f"source/executables/{executable.name}"
    if ledger.get(executable_member) != executable_sha256:
        raise ValueError("sealed executable SHA mismatch")
    tracked = {f"source/tracked/{relative}": ledger[f"source/tracked/{relative}"]
               for relative in _TRACKED_SOURCE_PATHS}
    tracked.update({
        "source/source_manifest.json": ledger["source/source_manifest.json"],
        "source/tooling/inverse_input_sealer.py": ledger["source/tooling/inverse_input_sealer.py"],
    })
    provenance = {
        "executable_sha256": executable_sha256,
        "tracked_source_sha256": tracked,
        "phase_base": phase_base,
        "tracked_source_binding": "current_workspace_bytes; phase_base verified separately",
        **extra_provenance,
    }
    manifest = {
        "schema": "vela.pn2d_minimal6_inverse_input.v1", "solver": solver,
        "bias_tolerance_V": BIAS_TOLERANCE_V,
        "fields": [{"name": name, "components": components, "unit": unit,
                    "support_kind": "node"}
                   for name, (components, unit) in contract.items()],
        "states": states, "member_sha256": ledger, "provenance": provenance,
    }
    manifest_path = root / "manifest.json"
    _write_json(manifest_path, manifest)
    _write_json(root / "seal.json", {"manifest_sha256": _sha256(manifest_path)})


def _run_import(importer: Path, tdr: Path, output: Path,
                runner: _ImporterRunner | None) -> None:
    try:
        if runner is not None:
            runner(tdr, output)
        else:
            result = subprocess.run(
                [str(importer), "--tdr", str(tdr), "--export-dir", str(output),
                 "--compensated-doping-policy", "reported"],
                text=True, capture_output=True, check=False,
            )
            if result.returncode != 0:
                raise RuntimeError(result.stderr or result.stdout or f"exit {result.returncode}")
    except Exception as error:
        raise ValueError(f"importer failed for {tdr.name}: {error}") from error
    if not output.is_dir():
        raise ValueError("importer did not create its export directory")

def _require_reimport_matches_ledger(export: Path, ledger: dict) -> None:
    actual = {path.relative_to(export).as_posix(): _sha256(path)
              for path in export.rglob("*") if path.is_file()}
    trusted_members = set(ledger)
    actual_members = set(actual)
    if "state.csv" in trusted_members:
        if (actual_members - {"state.csv"}
                != trusted_members - {"state.csv"}):
            raise ValueError("supplemental reimport mismatch: member set")
        try:
            from scripts.export_pn2d_minimal6_states import write_state_csv
            write_state_csv(export)
        except Exception as error:
            raise ValueError("supplemental reimport mismatch: state.csv regeneration") from error
        actual = {path.relative_to(export).as_posix(): _sha256(path)
                  for path in export.rglob("*") if path.is_file()}
        actual_members = set(actual)
    if actual_members != trusted_members:
        raise ValueError("supplemental reimport mismatch: member set")
    for relative, digest in ledger.items():
        if (not isinstance(digest, str)
                or actual[relative] != digest.lower()):
            raise ValueError(f"supplemental reimport mismatch: {relative}")

def _require_importer_sha(importer: Path, expected: str) -> None:
    try:
        actual = _sha256(importer)
    except OSError as error:
        raise ValueError("importer SHA changed during sealing") from error
    if actual != expected:
        raise ValueError("importer SHA changed during sealing")


def seal_inverse_input_roots(
    vela_sweep_root: str | Path,
    sentaurus_sweep_root: str | Path,
    supplemental_root: str | Path,
    output_root: str | Path,
    *,
    importer: str | Path,
    vela_executable: str | Path,
    phase_base: str,
    importer_runner: _ImporterRunner | None = None,
) -> dict[str, Path]:
    """Validate live roots and atomically create three canonical sealed roots."""
    destination = Path(output_root).resolve()
    if destination.exists():
        raise ValueError(f"output root already exists: {destination}")
    if not isinstance(phase_base, str) or not _PHASE_BASE.fullmatch(phase_base):
        raise ValueError("phase base must be a hexadecimal commit identifier")
    vela_root = _root(vela_sweep_root, "Vela sweep root")
    sentaurus_root = _root(sentaurus_sweep_root, "Sentaurus sweep root")
    supplemental_base = _root(supplemental_root, "supplemental root")
    importer_path = Path(importer).resolve()
    vela_executable_path = Path(vela_executable).resolve()
    if not importer_path.is_file() or not vela_executable_path.is_file():
        raise ValueError("declared executable is not a file")
    importer_sha256 = _sha256(importer_path)

    vela_manifest, vela_states, vela_decks = _validate_sweep(vela_root, "vela", include_zero=False)
    _, sent_states, sent_decks = _validate_sweep(
        sentaurus_root, "sentaurus", include_zero=True)
    supplemental_manifest, supplemental_states = _validate_supplemental(supplemental_base)
    declared_importer = supplemental_manifest.get("importer")
    if (not isinstance(declared_importer, str)
            or not Path(declared_importer).is_absolute()
            or os.path.normcase(str(Path(declared_importer).resolve()))
            != os.path.normcase(str(importer_path))):
        raise ValueError("supplemental importer identity mismatch")

    meshes = {}
    mesh_paths = {}
    doping_paths = {}
    topology_hashes = vela_manifest.get("topology_input_sha256")
    if not isinstance(topology_hashes, dict):
        raise ValueError("Vela topology input hash schema mismatch")
    for topology in TOPOLOGIES:
        mesh_path = _checked(vela_root, f"inputs/{topology}/mesh.json",
                             "Vela topology mesh")
        doping_path = _checked(vela_root, f"inputs/{topology}/doping.csv",
                               "Vela topology doping")
        declared = topology_hashes.get(topology)
        if not isinstance(declared, dict):
            raise ValueError("Vela topology input hash schema mismatch")
        _require_hash(mesh_path, declared.get("mesh.json"), "Vela topology mesh")
        _require_hash(doping_path, declared.get("doping.csv"), "Vela topology doping")
        meshes[topology] = _mesh_contract(mesh_path)
        mesh_paths[topology], doping_paths[topology] = mesh_path, doping_path

    vela_records = []
    for topology in TOPOLOGIES:
        for bias in COMMON_BIASES:
            _, state_path = vela_states[(topology, bias)]
            raw = _read_vela_state(state_path)
            fields = {}
            for name, (components, _) in REQUIRED_SENTARUS_FIELDS.items():
                values = {}
                for node in range(6):
                    source = raw[node]
                    mapped = {
                        "ElectrostaticPotential": (source["psi"],),
                        "eQuasiFermiPotential": (source["phin"],),
                        "hQuasiFermiPotential": (source["phip"],),
                        "eDensity": (source["electrons_m3"] / 1.0e6,),
                        "hDensity": (source["holes_m3"] / 1.0e6,),
                    }.get(name, tuple(None for _ in range(components)))
                    values[node] = mapped
                fields[name] = values
            tag = _bias_tag(bias)
            vela_records.append({
                "topology": topology, "bias": bias, "coordinates": meshes[topology]["coordinates"],
                "fields": fields,
                "evidence": [
                    (state_path, f"source/states/{topology}/{tag}.csv"),
                    (vela_decks[(topology, bias)], f"source/decks/{topology}/{tag}.json"),
                    (mesh_paths[topology], f"source/topologies/{topology}/mesh.json"),
                    (doping_paths[topology], f"source/topologies/{topology}/doping.csv"),
                ],
            })

    destination.parent.mkdir(parents=True, exist_ok=True)
    stage = Path(tempfile.mkdtemp(prefix=f".{destination.name}.staging-", dir=destination.parent))
    work = stage / "_work"
    try:
        fixed_importer = stage / "_tools" / importer_path.name
        fixed_importer.parent.mkdir(parents=True)
        _require_importer_sha(importer_path, importer_sha256)
        shutil.copyfile(importer_path, fixed_importer)
        _require_importer_sha(fixed_importer, importer_sha256)
        _require_importer_sha(importer_path, importer_sha256)
        sent_records = []
        supplemental_records = []
        for solver_label, state_map, deck_map, contract, target_records in (
            ("sentaurus", sent_states, sent_decks, REQUIRED_SENTARUS_FIELDS, sent_records),
            ("supplemental", supplemental_states, None, SUPPLEMENTAL_FIELDS, supplemental_records),
        ):
            for topology in TOPOLOGIES:
                for bias in COMMON_BIASES:
                    tag = _bias_tag(bias)
                    if solver_label == "sentaurus":
                        _, tdr = state_map[(topology, bias)]
                        trusted_ledger = None
                    else:
                        _, tdr, trusted_ledger = state_map[(topology, bias)]
                    export = work / solver_label / topology / tag
                    export.parent.mkdir(parents=True, exist_ok=True)
                    _require_importer_sha(importer_path, importer_sha256)
                    _require_importer_sha(fixed_importer, importer_sha256)
                    _run_import(fixed_importer, tdr, export, importer_runner)
                    _require_importer_sha(fixed_importer, importer_sha256)
                    _require_importer_sha(importer_path, importer_sha256)
                    if trusted_ledger is not None:
                        _require_reimport_matches_ledger(export, trusted_ledger)
                    imported = _validate_import(export, contract, meshes[topology], bias)
                    evidence = [(tdr, f"source/tdr/{topology}/{tag}.tdr")]
                    if deck_map is not None:
                        evidence.append((deck_map[(topology, bias)],
                                         f"source/decks/{topology}/{tag}.cmd"))
                    for name in ("field_manifest.json", "nodes.csv", "elements.csv", "contacts.csv"):
                        evidence.append((export / name,
                                         f"source/reimport/{topology}/{tag}/{name}"))
                    for name in contract:
                        evidence.append((export / "fields" / f"{name}_region0.csv",
                                         f"source/reimport/{topology}/{tag}/fields/{name}_region0.csv"))
                    target_records.append({
                        "topology": topology, "bias": bias,
                        "coordinates": imported["coordinates"], "fields": imported["fields"],
                        "evidence": evidence,
                    })
        _require_importer_sha(importer_path, importer_sha256)
        _require_importer_sha(fixed_importer, importer_sha256)
        _seal_root(stage / "vela", "vela", REQUIRED_SENTARUS_FIELDS, vela_records,
                   vela_executable_path, vela_root / "sweep_manifest.json", phase_base,
                   {"execution_binding_status": "post_hoc_observed",
                    "scientific_limitation": "unavailable native nodal E/J/alpha/source remain blank"})
        _seal_root(stage / "sentaurus", "sentaurus", REQUIRED_SENTARUS_FIELDS, sent_records,
                   fixed_importer, sentaurus_root / "sweep_manifest.json", phase_base,
                   {"remote_solver_binding_status": "not_declared_by_source_manifest"},
                   expected_executable_sha256=importer_sha256)
        _seal_root(stage / "supplemental", "supplemental", SUPPLEMENTAL_FIELDS,
                   supplemental_records, fixed_importer, supplemental_base / "manifest.json", phase_base,
                   {"remote_solver_binding_status": "release_declared_by_supplemental_manifest",
                    "sentaurus_version": supplemental_manifest["sentaurus_version"],
                    "declared_importer": str(importer_path),
                    "importer_sha256": importer_sha256},
                   expected_executable_sha256=importer_sha256)
        shutil.rmtree(work)
        load_input_bundle(stage / "vela", stage / "sentaurus", stage / "supplemental")
        shutil.rmtree(fixed_importer.parent)
        os.replace(stage, destination)
    except Exception:
        shutil.rmtree(stage, ignore_errors=True)
        raise
    return {name: destination / name for name in ("vela", "sentaurus", "supplemental")}
