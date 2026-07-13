#!/usr/bin/env python3
"""Canonical fixed-state PN2D minimal6 topologies and DF-ISE text writers."""

from __future__ import annotations

import argparse
import json
import re
import sys
from collections import Counter
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.compare_sentaurus_tdr_tdx import parse_dat, parse_grd


_NODES = {
    1: (0.0, 0.5),
    2: (1.0, 0.5),
    3: (2.0, 0.5),
    4: (2.0, 0.0),
    5: (0.0, 0.0),
    6: (1.0, 0.0),
}
_CONTACTS = {"Anode": (1, 5), "Cathode": (3, 4)}
_ACCEPTORS = {1: 1.0e17, 2: 1.0e17, 3: 0.0, 4: 0.0, 5: 1.0e17, 6: 1.0e17}
_DONORS = {1: 0.0, 2: 1.0e17, 3: 1.0e17, 4: 1.0e17, 5: 0.0, 6: 1.0e17}
_TRIANGLES = {
    "sketch": ((1, 5, 2), (5, 6, 2), (2, 6, 4), (2, 4, 3)),
    "mirror": ((1, 5, 6), (1, 6, 2), (2, 6, 3), (6, 4, 3)),
}
_MIRROR_LABELS = {1: 5, 2: 6, 3: 4, 4: 3, 5: 1, 6: 2}
_REGION_OWNERSHIP = {"R.Si": [0, 1, 2, 3], "Cathode": [4], "Anode": [5]}
_REGION_ORDER = ["R.Si", "Cathode", "Anode"]
_MATERIALS = ["Silicon", "Contact", "Contact"]


@dataclass(frozen=True)
class Topology:
    topology_id: str
    nodes: Mapping[int, tuple[float, float]]
    triangles: tuple[tuple[int, int, int], ...]
    contacts: Mapping[str, tuple[int, int]]
    acceptors_cm3: Mapping[int, float]
    donors_cm3: Mapping[int, float]

    def __post_init__(self) -> None:
        object.__setattr__(self, "nodes", MappingProxyType(dict(self.nodes)))
        object.__setattr__(
            self,
            "triangles",
            tuple(tuple(triangle) for triangle in self.triangles),
        )
        object.__setattr__(self, "contacts", MappingProxyType(dict(self.contacts)))
        object.__setattr__(
            self, "acceptors_cm3", MappingProxyType(dict(self.acceptors_cm3))
        )
        object.__setattr__(self, "donors_cm3", MappingProxyType(dict(self.donors_cm3)))


@dataclass(frozen=True)
class TopologySummary:
    nodes: int
    triangles: int
    edges: int
    contact_edges: dict[str, tuple[int, int]]


def canonical_edges(
    triangles: Iterable[tuple[int, int, int]],
) -> list[tuple[int, int]]:
    return sorted(
        {
            tuple(sorted(edge))
            for triangle in triangles
            for edge in (
                (triangle[0], triangle[1]),
                (triangle[1], triangle[2]),
                (triangle[2], triangle[0]),
            )
        }
    )


def canonical_triangle(triangle: tuple[int, int, int]) -> tuple[int, int, int]:
    """Normalize a triangle cyclically without changing its orientation."""
    return min(triangle[index:] + triangle[:index] for index in range(3))


_CANONICAL_NODE_KEYS = {str(node_id) for node_id in _NODES}


def _reject_duplicate_json_keys(
    pairs: list[tuple[str, object]],
) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key: {key}")
        result[key] = value
    return result


def _require_canonical_node_keys(data: dict[str, object], name: str) -> None:
    values = data.get(name)
    if not isinstance(values, dict) or set(values) != _CANONICAL_NODE_KEYS:
        raise ValueError(
            f"{name} keys must be exactly {sorted(_CANONICAL_NODE_KEYS)}"
        )


def _as_int_keyed(values: dict[str, object]) -> dict[int, object]:
    return {int(key): value for key, value in values.items()}


def load_topology(path: Path, topology_id: str) -> Topology:
    data = json.loads(
        path.read_text(encoding="utf-8"),
        object_pairs_hook=_reject_duplicate_json_keys,
    )
    if data.get("schema") != "vela.pn2d_minimal6_topologies.v1":
        raise ValueError("unexpected minimal6 topology schema")
    if data.get("length_unit") != "um" or data.get("doping_unit") != "cm^-3":
        raise ValueError("unexpected minimal6 topology units")
    for name in ("nodes", "acceptors_cm3", "donors_cm3"):
        _require_canonical_node_keys(data, name)
    try:
        triangles = data["topologies"][topology_id]
    except KeyError as error:
        raise ValueError(f"unknown topology ID: {topology_id}") from error
    topology = Topology(
        topology_id=topology_id,
        nodes={
            node_id: tuple(float(value) for value in point)
            for node_id, point in _as_int_keyed(data["nodes"]).items()
        },
        triangles=[tuple(int(node_id) for node_id in triangle) for triangle in triangles],
        contacts={
            name: tuple(int(node_id) for node_id in edge)
            for name, edge in data["contacts"].items()
        },
        acceptors_cm3={
            node_id: float(value)
            for node_id, value in _as_int_keyed(data["acceptors_cm3"]).items()
        },
        donors_cm3={
            node_id: float(value)
            for node_id, value in _as_int_keyed(data["donors_cm3"]).items()
        },
    )
    validate_topology(topology)
    return topology


def _signed_area2(topology: Topology, triangle: tuple[int, int, int]) -> float:
    (ax, ay), (bx, by), (cx, cy) = (
        topology.nodes[node_id] for node_id in triangle
    )
    return (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)


def validate_topology(topology: Topology) -> TopologySummary:
    if topology.topology_id not in _TRIANGLES:
        raise ValueError(f"unknown topology ID: {topology.topology_id}")
    if topology.nodes != _NODES:
        raise ValueError("minimal6 nodes must use the six approved coordinates")
    if topology.contacts != _CONTACTS:
        raise ValueError("minimal6 contacts must remain unchanged")
    if topology.acceptors_cm3 != _ACCEPTORS or topology.donors_cm3 != _DONORS:
        raise ValueError("minimal6 doping must remain canonical")
    if (
        topology.acceptors_cm3[2] != topology.donors_cm3[2]
        or topology.acceptors_cm3[6] != topology.donors_cm3[6]
    ):
        raise ValueError("nodes 2 and 6 must remain compensated")
    if topology.triangles != _TRIANGLES[topology.topology_id]:
        raise ValueError("minimal6 triangle connectivity must remain canonical")
    if len(topology.nodes) != 6 or len(topology.triangles) != 4:
        raise ValueError("minimal6 node or triangle count mismatch")
    if len(canonical_edges(topology.triangles)) != 9:
        raise ValueError("minimal6 edge count mismatch")
    if any(_signed_area2(topology, triangle) <= 0.0 for triangle in topology.triangles):
        raise ValueError("minimal6 triangles must be CCW")
    if topology.topology_id == "mirror":
        reflected_nodes = {
            _MIRROR_LABELS[node_id]: (x, 0.5 - y)
            for node_id, (x, y) in _NODES.items()
        }
        reflected_triangles = {
            canonical_triangle(
                tuple(
                    reversed(
                        tuple(_MIRROR_LABELS[node_id] for node_id in triangle)
                    )
                )
            )
            for triangle in _TRIANGLES["sketch"]
        }
        if topology.nodes != reflected_nodes or {
            canonical_triangle(triangle) for triangle in topology.triangles
        } != reflected_triangles:
            raise ValueError("mirror topology must be the labelled vertical reflection")
    return TopologySummary(
        nodes=len(topology.nodes),
        triangles=len(topology.triangles),
        edges=len(canonical_edges(topology.triangles)),
        contact_edges=topology.contacts,
    )


def _edge_reference(
    edge: tuple[int, int], edge_ids: dict[tuple[int, int], int]
) -> int:
    canonical = tuple(sorted(edge))
    edge_id = edge_ids[canonical]
    return edge_id if edge == canonical else -(edge_id + 1)


def _edge_locations(
    triangles: Iterable[tuple[int, int, int]], edges: Sequence[tuple[int, int]]
) -> list[str]:
    incidence = Counter(
        tuple(sorted(edge))
        for triangle in triangles
        for edge in (
            (triangle[0], triangle[1]),
            (triangle[1], triangle[2]),
            (triangle[2], triangle[0]),
        )
    )
    if any(incidence[edge] not in (1, 2) for edge in edges):
        raise ValueError("minimal6 edges must have one or two triangle owners")
    return ["e" if incidence[edge] == 1 else "i" for edge in edges]


def write_dfise_grid(topology: Topology, path: Path) -> None:
    validate_topology(topology)
    edges = canonical_edges(topology.triangles)
    edge_ids = {edge: index for index, edge in enumerate(edges)}
    locations = _edge_locations(topology.triangles, edges)
    triangle_records = [
        "2 "
        + " ".join(
            str(_edge_reference(edge, edge_ids))
            for edge in (
                (triangle[0], triangle[1]),
                (triangle[1], triangle[2]),
                (triangle[2], triangle[0]),
            )
        )
        for triangle in topology.triangles
    ]
    vertices = "\n".join(
        f" {x:.15e} {y:.15e}" for _, (x, y) in sorted(topology.nodes.items())
    )
    edge_rows = "\n".join(f" {start - 1} {end - 1}" for start, end in edges)
    elements = "\n".join(
        [f" {record}" for record in triangle_records]
        + [
            f" 1 {topology.contacts['Cathode'][0] - 1} {topology.contacts['Cathode'][1] - 1}",
            f" 1 {topology.contacts['Anode'][0] - 1} {topology.contacts['Anode'][1] - 1}",
        ]
    )
    text = f"""DF-ISE text

Info {{
  version = 1.1
  type    = grid
  dimension   = 2
  nb_vertices = 6
  nb_edges = 9
  nb_faces = 0
  nb_elements = 6
  nb_regions = 3
  regions = [ "R.Si" "Cathode" "Anode" ]
  materials = [ Silicon Contact Contact ]
}}

Data {{

  CoordSystem {{
    translate = [  0.000000000000000e+00 0.000000000000000e+00 0.000000000000000e+00 ]
    transform = [  1.000000000000000e+00 0.000000000000000e+00 0.000000000000000e+00
 0.000000000000000e+00 1.000000000000000e+00 0.000000000000000e+00
 0.000000000000000e+00 0.000000000000000e+00 1.000000000000000e+00
 ]
  }}

  Vertices (6) {{
{vertices}
  }}

  Edges (9) {{
{edge_rows}
  }}

  Locations (9) {{
{''.join(locations)}
  }}

  Elements (6) {{
{elements}
  }}

  Region ("R.Si") {{
    material = Silicon
    Elements (4) {{
 0 1 2 3
    }}
  }}

  Region ("Cathode") {{
    material = Contact
    Elements (1) {{
 4
    }}
  }}

  Region ("Anode") {{
    material = Contact
    Elements (1) {{
 5
    }}
  }}

}}
"""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def _doping_datasets(topology: Topology) -> dict[str, list[float]]:
    return {
        "DopingConcentration": [
            topology.donors_cm3[node_id] - topology.acceptors_cm3[node_id]
            for node_id in sorted(topology.nodes)
        ],
        "PhosphorusActiveConcentration": [
            topology.donors_cm3[node_id] for node_id in sorted(topology.nodes)
        ],
        "BoronActiveConcentration": [
            topology.acceptors_cm3[node_id] for node_id in sorted(topology.nodes)
        ],
    }


def write_dfise_doping(topology: Topology, path: Path) -> None:
    validate_topology(topology)
    datasets = _doping_datasets(topology)
    blocks = []
    for name, values in datasets.items():
        value_rows = "\n".join(f" {value:.15e}" for value in values)
        blocks.append(
            f"""  Dataset ("{name}") {{
    function  = {name}
    type      = scalar
    dimension = 1
    location  = vertex
    validity  = [ "R.Si" ]
    Values (6) {{
{value_rows}
    }}
  }}"""
        )
    names = " ".join(f'"{name}"' for name in datasets)
    functions = " ".join(datasets)
    dataset_blocks = "\n\n".join(blocks)
    text = f"""DF-ISE text

Info {{
  version = 1.0
  type    = dataset
  dimension   = 2
  nb_vertices = 6
  nb_edges    = 9
  nb_faces    = 0
  nb_elements = 6
  nb_regions  = 3
  datasets    = [ {names} ]
  functions   = [ {functions} ]
}}

Data {{

{dataset_blocks}

}}
"""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


_NUMBER_RE = re.compile(r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?")


def _metadata_assignment(body: str, name: str) -> str | None:
    match = re.search(rf"\b{re.escape(name)}\s*=\s*([^\s\]\}}]+)", body)
    return match.group(1) if match else None


def _metadata_list(body: str, name: str) -> list[str]:
    match = re.search(rf"\b{re.escape(name)}\s*=\s*\[(.*?)\]", body, re.S)
    if not match:
        return []
    raw = match.group(1)
    token_pattern = re.compile(r'\s*(?:"([^"]+)"|([A-Za-z0-9_.:+-]+))')
    tokens: list[str] = []
    position = 0
    while position < len(raw):
        if not raw[position:].strip():
            break
        token = token_pattern.match(raw, position)
        if not token:
            return []
        tokens.append(token.group(1) or token.group(2))
        position = token.end()
    return tokens


def _metadata_numbers(body: str, name: str) -> tuple[float, ...]:
    match = re.search(rf"\b{re.escape(name)}\s*=\s*\[(.*?)\]", body, re.S)
    if not match:
        return ()
    return tuple(float(value) for value in _NUMBER_RE.findall(match.group(1)))


_INTEGER_RE = re.compile(r"[-+]?\d+")


def _metadata_assignment_names(body: str) -> list[str]:
    return re.findall(r"\b([A-Za-z_][A-Za-z0-9_]*)\s*=", body)


def _strict_numeric_tokens(body: str, *, integer: bool) -> list[str] | None:
    tokens = body.split()
    pattern = _INTEGER_RE if integer else _NUMBER_RE
    if not tokens or any(not pattern.fullmatch(token) for token in tokens):
        return None
    return tokens


def _strict_numeric_rows(
    text: str,
    block_name: str,
    expected_count: int,
    allowed_columns: set[int],
    *,
    integer: bool,
) -> bool:
    match = re.search(
        rf"\b{re.escape(block_name)}\s*\((\d+)\)\s*"
        rf"\{{(?P<body>.*?)\n\s*\}}",
        text,
        re.S,
    )
    if not match or int(match.group(1)) != expected_count:
        return False
    rows = [
        line.split()
        for line in match.group("body").splitlines()
        if line.strip()
    ]
    pattern = _INTEGER_RE if integer else _NUMBER_RE
    return (
        len(rows) == expected_count
        and all(len(row) in allowed_columns for row in rows)
        and all(pattern.fullmatch(token) for row in rows for token in row)
    )


def _strict_grid_metadata(text: str) -> bool:
    info_matches = list(re.finditer(r"\bInfo\s*\{(?P<body>.*?)\n\}", text, re.S))
    coord_matches = list(
        re.finditer(r"\bCoordSystem\s*\{(?P<body>.*?)\n\s*\}", text, re.S)
    )
    region_names = re.findall(r'Region\s*\("([^"]+)"\)\s*\{', text)
    if (
        len(info_matches) != 1
        or len(coord_matches) != 1
        or len(re.findall(r"\bData\s*\{", text)) != 1
        or region_names != _REGION_ORDER
        or not re.search(r"\}\s*\Z", text)
    ):
        return False
    info_match = info_matches[0]
    coord_match = coord_matches[0]
    info = info_match.group("body")
    expected_scalars = {
        "version": "1.1",
        "type": "grid",
        "dimension": "2",
        "nb_vertices": "6",
        "nb_edges": "9",
        "nb_faces": "0",
        "nb_elements": "6",
        "nb_regions": "3",
    }
    if any(
        _metadata_assignment(info, name) != value
        for name, value in expected_scalars.items()
    ):
        return False
    if _metadata_assignment_names(info) != [
        *expected_scalars,
        "regions",
        "materials",
    ]:
        return False
    if _metadata_list(info, "regions") != _REGION_ORDER:
        return False
    if _metadata_list(info, "materials") != _MATERIALS:
        return False

    coord = coord_match.group("body")
    if _metadata_assignment_names(coord) != ["translate", "transform"]:
        return False
    if _metadata_numbers(coord, "translate") != (0.0, 0.0, 0.0):
        return False
    if _metadata_numbers(coord, "transform") != (
        1.0,
        0.0,
        0.0,
        0.0,
        1.0,
        0.0,
        0.0,
        0.0,
        1.0,
    ):
        return False

    if not _strict_numeric_rows(text, "Vertices", 6, {2}, integer=False):
        return False
    if not _strict_numeric_rows(text, "Edges", 9, {2}, integer=True):
        return False
    if not _strict_numeric_rows(text, "Elements", 6, {3, 4}, integer=True):
        return False
    locations_match = re.search(r"\bLocations\s*\((\d+)\)", text)
    if not locations_match or int(locations_match.group(1)) != 9:
        return False
    expected_materials = dict(zip(_REGION_ORDER, _MATERIALS, strict=True))
    for region_name, element_ids in _REGION_OWNERSHIP.items():
        match = re.search(
            rf'Region\s*\("{re.escape(region_name)}"\)\s*\{{\s*'
            rf'material\s*=\s*{re.escape(expected_materials[region_name])}\s+'
            rf'Elements\s*\((\d+)\)\s*\{{(?P<values>.*?)\n\s*\}}',
            text,
            re.S,
        )
        if not match or int(match.group(1)) != len(element_ids):
            return False
        tokens = _strict_numeric_tokens(match.group("values"), integer=True)
        if tokens is None or len(tokens) != len(element_ids):
            return False
    return True


def _strict_dataset_metadata(text: str, expected_names: list[str]) -> bool:
    info_matches = list(re.finditer(r"\bInfo\s*\{(?P<body>.*?)\n\}", text, re.S))
    dataset_names = re.findall(r'Dataset\s*\("([^"]+)"\)\s*\{', text)
    if (
        len(info_matches) != 1
        or len(re.findall(r"\bData\s*\{", text)) != 1
        or dataset_names != expected_names
        or not re.search(r"\}\s*\Z", text)
    ):
        return False
    info = info_matches[0].group("body")
    expected_scalars = {
        "version": "1.0",
        "type": "dataset",
        "dimension": "2",
        "nb_vertices": "6",
        "nb_edges": "9",
        "nb_faces": "0",
        "nb_elements": "6",
        "nb_regions": "3",
    }
    if any(
        _metadata_assignment(info, name) != value
        for name, value in expected_scalars.items()
    ):
        return False
    if _metadata_assignment_names(info) != [
        *expected_scalars,
        "datasets",
        "functions",
    ]:
        return False
    if _metadata_list(info, "datasets") != expected_names:
        return False
    if _metadata_list(info, "functions") != expected_names:
        return False

    for name in expected_names:
        match = re.search(
            rf'Dataset\s*\("{re.escape(name)}"\)\s*\{{'
            rf'(?P<header>.*?)\bValues\s*\((\d+)\)\s*'
            rf'\{{(?P<values>.*?)\n\s*\}}',
            text,
            re.S,
        )
        if not match or int(match.group(2)) != 6:
            return False
        header = match.group("header")
        if _metadata_assignment_names(header) != [
            "function",
            "type",
            "dimension",
            "location",
            "validity",
        ]:
            return False
        values = _strict_numeric_tokens(match.group("values"), integer=False)
        if values is None or len(values) != 6:
            return False
        if _metadata_assignment(header, "function") != name:
            return False
        if _metadata_assignment(header, "type") != "scalar":
            return False
        if _metadata_assignment(header, "dimension") != "1":
            return False
        if _metadata_assignment(header, "location") != "vertex":
            return False
        if _metadata_list(header, "validity") != ["R.Si"]:
            return False
    return True


def _decode_triangle(
    record: list[int], edges: list[list[int]]
) -> tuple[int, int, int]:
    if record[0] != 2 or len(record) != 4:
        raise ValueError("expected a triangular DF-ISE element")
    oriented = []
    for reference in record[1:]:
        edge_index = reference if reference >= 0 else -reference - 1
        if edge_index < 0 or edge_index >= len(edges):
            raise ValueError(f"DF-ISE edge reference {reference} is out of range")
        if len(edges[edge_index]) != 2:
            raise ValueError(f"DF-ISE edge reference {reference} is not a segment")
        edge = (
            edges[edge_index]
            if reference >= 0
            else list(reversed(edges[edge_index]))
        )
        oriented.append(edge)
    if (
        oriented[0][1] != oriented[1][0]
        or oriented[1][1] != oriented[2][0]
        or oriented[2][1] != oriented[0][0]
    ):
        raise ValueError("DF-ISE triangle edge loop is disconnected")
    return tuple(
        node_id + 1
        for node_id in (oriented[0][0], oriented[0][1], oriented[1][1])
    )


def _decode_contact(record: list[int]) -> tuple[int, int] | None:
    if len(record) != 3 or record[0] != 1:
        return None
    return (record[1] + 1, record[2] + 1)


def _parse_locations(path: Path) -> list[str]:
    text = path.read_text(encoding="utf-8", errors="replace")
    match = re.search(
        r"\bLocations\s*\((\d+)\)\s*\{(?P<body>.*?)\n\s*\}", text, re.S
    )
    if not match:
        return []
    raw_locations = re.sub(r"\s+", "", match.group("body"))
    if not re.fullmatch(r"[ei]+", raw_locations):
        return []
    locations = list(raw_locations)
    if len(locations) != int(match.group(1)):
        return []
    return locations


def validate_dfise_roundtrip(
    topology: Topology, grd: Path, dat: Path
) -> dict[str, object]:
    validate_topology(topology)
    grid_text = grd.read_text(encoding="utf-8", errors="replace")
    dataset_text = dat.read_text(encoding="utf-8", errors="replace")
    grid = parse_grd(grd)
    datasets = parse_dat(dat)
    grid_metadata_matches = _strict_grid_metadata(grid_text)
    edges = grid["edges"]
    elements = grid["elements"]
    silicon_ids = grid["region_elements"].get("R.Si", [])
    reconstructed = tuple(
        _decode_triangle(elements[element_id], edges) for element_id in silicon_ids
    )
    region_ownership = {
        name: grid["region_elements"].get(name, []) for name in _REGION_ORDER
    }
    contacts = {
        name: len(region_ownership[name]) for name in ("Anode", "Cathode")
    }
    contact_edges = {}
    for name in ("Anode", "Cathode"):
        owned = region_ownership[name]
        contact_edges[name] = (
            _decode_contact(elements[owned[0]])
            if len(owned) == 1 and owned[0] < len(elements)
            else None
        )

    vertices = {
        node_id: tuple(vertex)
        for node_id, vertex in enumerate(grid["vertices"], start=1)
    }
    parsed_edges = [
        (int(edge[0]) + 1, int(edge[1]) + 1)
        for edge in edges
        if len(edge) == 2
    ]
    locations = _parse_locations(grd)
    location_counts = {"e": locations.count("e"), "i": locations.count("i")}
    expected_edges = canonical_edges(topology.triangles)
    expected_locations = _edge_locations(topology.triangles, expected_edges)
    expected_datasets = _doping_datasets(topology)
    dataset_metadata_matches = _strict_dataset_metadata(
        dataset_text, list(expected_datasets)
    )
    doping_matches = (
        set(datasets) == set(expected_datasets)
        and all(
            datasets[name].get("count") == 6
            and datasets[name].get("dimension") == 1
            and datasets[name].get("values") == values
            for name, values in expected_datasets.items()
        )
    )
    passed = (
        grid_metadata_matches
        and dataset_metadata_matches
        and grid["vertex_count"] == 6
        and vertices == topology.nodes
        and len(edges) == 9
        and parsed_edges == expected_edges
        and locations == expected_locations
        and grid["element_count"] == 6
        and len(elements) == 6
        and grid["regions"] == _REGION_ORDER
        and grid["materials"] == _MATERIALS
        and region_ownership == _REGION_OWNERSHIP
        and reconstructed == topology.triangles
        and contacts == {"Anode": 1, "Cathode": 1}
        and contact_edges == topology.contacts
        and doping_matches
    )
    return {
        "passed": passed,
        "grid_metadata_matches": grid_metadata_matches,
        "dataset_metadata_matches": dataset_metadata_matches,
        "vertices": vertices,
        "edge_count": len(edges),
        "edges": parsed_edges,
        "location_counts": location_counts,
        "locations": locations,
        "triangles": reconstructed,
        "silicon_triangle_count": len(reconstructed),
        "contact_element_counts": contacts,
        "contact_edges": contact_edges,
        "region_ownership": region_ownership,
    }


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Validate or emit a canonical PN2D minimal6 topology."
    )
    parser.add_argument(
        "--fixture",
        type=Path,
        default=Path(__file__).resolve().parents[1]
        / "reference_tcad"
        / "pn2d_sentaurus2018_minimal6"
        / "source"
        / "minimal6_topologies.json",
    )
    parser.add_argument("--topology", choices=sorted(_TRIANGLES), default="sketch")
    parser.add_argument("--grd", type=Path)
    parser.add_argument("--dat", type=Path)
    return parser


def main() -> int:
    parser = _build_parser()
    args = parser.parse_args()
    topology = load_topology(args.fixture, args.topology)
    if (args.grd is None) != (args.dat is None):
        parser.error("--grd and --dat must be provided together")
    if args.grd is not None:
        write_dfise_grid(topology, args.grd)
        write_dfise_doping(topology, args.dat)
        report = validate_dfise_roundtrip(topology, args.grd, args.dat)
        print(json.dumps(report, indent=2, sort_keys=True))
        return 0 if report["passed"] else 1
    summary = validate_topology(topology)
    print(
        json.dumps(
            {
                "nodes": summary.nodes,
                "triangles": summary.triangles,
                "edges": summary.edges,
                "contact_edges": dict(summary.contact_edges),
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
