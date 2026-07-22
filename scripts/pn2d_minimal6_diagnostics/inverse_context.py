"""Ledger-bound topology and Vela production physics context."""

from __future__ import annotations

import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping


TOPOLOGIES = ("sketch", "mirror")
HEADER_RELATIVE = "source/tracked/include/vela/physics/ImpactIonizationModel.h"
MESH_RELATIVE = "source/topologies/{topology}/mesh.json"
_NUMBER = r"([0-9.eE+-]+)"
_BOLTZMANN_EV_PER_K = 8.617333262145e-5


@dataclass(frozen=True)
class DiagnosticContext:
    mesh_by_topology: dict[str, dict[str, dict[str, tuple[str, ...]]]]
    van_overstraeten_parameters: dict[str, object]
    thermal_voltage_V: float
    provenance: dict[str, object]


def _finite(value: object, label: str, *, positive: bool = False) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{label} must be finite") from error
    if not math.isfinite(result) or (positive and result <= 0.0):
        raise ValueError(f"{label} must be finite and positive")
    return result


def _sealed_path(root: Path, hashes: Mapping[str, str], relative: str) -> Path:
    if relative not in hashes:
        raise ValueError(f"Vela sealed context lacks {relative}")
    path = (root / relative).resolve()
    try:
        path.relative_to(root.resolve())
    except ValueError as error:
        raise ValueError("Vela sealed context path escapes root") from error
    if not path.is_file():
        raise ValueError(f"Vela sealed context lacks {relative}")
    return path


def _node_key(value: str) -> tuple[int, object]:
    return (0, int(value)) if value.isdecimal() else (1, value)


def _load_mesh(path: Path, topology: str) -> dict[str, dict[str, tuple[str, ...]]]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ValueError(f"sealed {topology} mesh is not valid JSON") from error
    if not isinstance(payload, dict) or payload.get("coordinate_unit") != "um":
        raise ValueError(f"sealed {topology} mesh coordinate contract mismatch")
    raw_nodes, raw_triangles = payload.get("nodes"), payload.get("triangles")
    if not isinstance(raw_nodes, list) or not isinstance(raw_triangles, list):
        raise ValueError(f"sealed {topology} mesh topology contract mismatch")
    coordinates: dict[str, tuple[float, float]] = {}
    for item in raw_nodes:
        if not isinstance(item, dict):
            raise ValueError(f"sealed {topology} mesh node contract mismatch")
        node = str(item.get("id"))
        if node in coordinates:
            raise ValueError(f"sealed {topology} mesh has duplicate node")
        coordinates[node] = (
            _finite(item.get("x"), f"{topology} mesh x"),
            _finite(item.get("y"), f"{topology} mesh y"),
        )
    triangles: dict[str, tuple[str, ...]] = {}
    for item in raw_triangles:
        if not isinstance(item, dict) or not isinstance(item.get("node_ids"), list):
            raise ValueError(f"sealed {topology} mesh triangle contract mismatch")
        cell = str(item.get("id"))
        nodes = tuple(str(node) for node in item["node_ids"])
        if cell in triangles or len(nodes) != 3 or len(set(nodes)) != 3:
            raise ValueError(f"sealed {topology} mesh triangle contract mismatch")
        if any(node not in coordinates for node in nodes):
            raise ValueError(f"sealed {topology} mesh triangle references unknown node")
        points = tuple(coordinates[node] for node in nodes)
        twice_area = (
            (points[1][0] - points[0][0]) * (points[2][1] - points[0][1])
            - (points[2][0] - points[0][0]) * (points[1][1] - points[0][1])
        )
        if not math.isfinite(twice_area) or abs(twice_area) <= 1.0e-300:
            raise ValueError(f"sealed {topology} mesh has degenerate triangle")
        triangles[cell] = nodes
    if not triangles:
        raise ValueError(f"sealed {topology} mesh has no triangles")
    edge_pairs = {
        tuple(sorted((nodes[index], nodes[(index + 1) % 3]), key=_node_key))
        for nodes in triangles.values()
        for index in range(3)
    }
    edges = {
        str(index): pair
        for index, pair in enumerate(sorted(edge_pairs, key=lambda pair: (
            _node_key(pair[0]), _node_key(pair[1])
        )))
    }
    return {
        "triangles": dict(sorted(triangles.items(), key=lambda item: _node_key(item[0]))),
        "edges": edges,
    }


def _header_value(text: str, name: str) -> float:
    matches = re.findall(rf"\bReal\s+{re.escape(name)}\s*=\s*{_NUMBER}\s*;", text)
    if len(matches) != 1:
        raise ValueError(f"Vela production header lacks unique numeric {name}")
    return _finite(matches[0], f"Vela production {name}", positive=True)


def _load_parameters(path: Path) -> tuple[dict[str, object], float]:
    text = path.read_text(encoding="utf-8")
    values = {
        name: _header_value(text, name)
        for name in (
            "electronALow", "electronAHigh", "electronBLow", "electronBHigh",
            "holeALow", "holeAHigh", "holeBLow", "holeBHigh", "switchField",
            "phononEnergy", "referenceTemperature_K", "temperature_K",
        )
    }
    reference = values["referenceTemperature_K"]
    temperature = values["temperature_K"]
    phonon = values["phononEnergy"]
    gamma = (
        math.tanh(phonon / (2.0 * _BOLTZMANN_EV_PER_K * reference))
        / math.tanh(phonon / (2.0 * _BOLTZMANN_EV_PER_K * temperature))
    )
    parameters: dict[str, object] = {
        "gamma": gamma,
        "switch_field_V_m": values["switchField"],
        "electron": {
            "low": (values["electronALow"], values["electronBLow"]),
            "high": (values["electronAHigh"], values["electronBHigh"]),
        },
        "hole": {
            "low": (values["holeALow"], values["holeBLow"]),
            "high": (values["holeAHigh"], values["holeBHigh"]),
        },
        "phonon_energy_eV": phonon,
        "reference_temperature_K": reference,
        "temperature_K": temperature,
    }
    return parameters, _BOLTZMANN_EV_PER_K * temperature


def load_sealed_vela_context(
    root: str | Path, member_hashes: Mapping[str, str],
) -> DiagnosticContext:
    base = Path(root).resolve()
    mesh = {}
    mesh_provenance = {}
    for topology in TOPOLOGIES:
        relative = MESH_RELATIVE.format(topology=topology)
        path = _sealed_path(base, member_hashes, relative)
        mesh[topology] = _load_mesh(path, topology)
        mesh_provenance[topology] = {
            "logical_id": f"vela:{relative}", "sha256": member_hashes[relative],
        }
    header_path = _sealed_path(base, member_hashes, HEADER_RELATIVE)
    parameters, thermal_voltage = _load_parameters(header_path)
    provenance = {
        "parameter_identity": "sealed_vela_production_defaults",
        "parameter_source": f"vela:{HEADER_RELATIVE}",
        "parameter_sha256": member_hashes[HEADER_RELATIVE],
        "mesh_sources": mesh_provenance,
    }
    return DiagnosticContext(
        mesh, parameters, thermal_voltage, provenance,
    )
