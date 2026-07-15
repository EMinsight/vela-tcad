from __future__ import annotations
import math

def project_vector_to_edge(vector, start, end):
    dx,dy = float(end[0])-float(start[0]), float(end[1])-float(start[1])
    length = math.hypot(dx,dy)
    if length == 0.0: raise ValueError("edge has zero length")
    return (float(vector[0])*dx+float(vector[1])*dy)/length
def integrate_nodal_field(values, weights):
    if len(values) != len(weights): raise ValueError("values and weights must align")
    return sum(float(value) * float(weight) for value, weight in zip(values, weights))

def map_local_sources_to_nodes(triangles, sources):
    if len(triangles) != len(sources): raise ValueError("triangles and sources must align")
    result = {}
    for nodes, source in zip(triangles, sources):
        if len(nodes) != 3 or len(set(nodes)) != 3: raise ValueError("local source needs three unique nodes")
        share = float(source) / 3.0
        for node in nodes: result[node] = result.get(node, 0.0) + share
    if abs(sum(result.values()) - sum(float(source) for source in sources)) > 1e-12: raise AssertionError("node mapping is not conservative")
    return result
def integrate_cell_field(points, values, *, partial_volume_fraction=1.0):
    if len(points) != 3 or len(values) != 3: raise ValueError("cell field requires one triangle")
    (x0,y0),(x1,y1),(x2,y2) = ((float(x),float(y)) for x,y in points)
    twice_area = (x1-x0)*(y2-y0)-(x2-x0)*(y1-y0)
    if twice_area <= 0.0: raise ValueError("triangle must be strictly CCW")
    fraction = float(partial_volume_fraction)
    if fraction < 0.0 or fraction > 1.0: raise ValueError("partial-volume fraction must be in [0, 1]")
    return fraction * 0.5 * twice_area * sum(float(value) for value in values) / 3.0
def node_scalar_to_cells(node_values, cells, *, quantity=None):
    """P1 node-to-cell averaging with explicit normalized node weights."""
    if quantity in {"ImpactIonization", "AvalancheGeneration"}:
        raise ValueError("native avalanche generation must be integrated, not averaged")
    values, weights = [], []
    for nodes in cells:
        ids = tuple(nodes)
        if len(ids) != 3 or len(set(ids)) != 3 or any(node not in node_values for node in ids):
            raise ValueError("cell needs three distinct nodal values")
        weight = {node: 1.0 / 3.0 for node in ids}
        weights.append(weight)
        values.append(sum(float(node_values[node]) * weight[node] for node in ids))
    return {"values": values, "weights": weights}


def node_vector_to_edges(node_vectors, edges, coordinates):
    """Project endpoint-averaged nodal vectors onto directed edges."""
    values, weights = [], []
    for start, end in edges:
        if (start == end or start not in node_vectors or end not in node_vectors
                or start not in coordinates or end not in coordinates):
            raise ValueError("edge needs two distinct nodal vectors and coordinates")
        weight = {start: 0.5, end: 0.5}
        vector = (
            0.5 * (float(node_vectors[start][0]) + float(node_vectors[end][0])),
            0.5 * (float(node_vectors[start][1]) + float(node_vectors[end][1])),
        )
        values.append(project_vector_to_edge(vector, coordinates[start], coordinates[end]))
        weights.append(weight)
    return {"values": values, "weights": weights}


def local_edge_sources_to_nodes(edges, sources):
    """Conservatively distribute each extensive local-edge source to endpoints."""
    if len(edges) != len(sources):
        raise ValueError("local edges and sources must align")
    values, weights = {}, []
    for edge, source in zip(edges, sources):
        if len(edge) != 2 or edge[0] == edge[1]:
            raise ValueError("local source needs two distinct edge endpoints")
        start, end = edge
        weight = {start: 0.5, end: 0.5}
        for node in (start, end):
            values[node] = values.get(node, 0.0) + float(source) * weight[node]
        weights.append(weight)
    if abs(sum(values.values()) - sum(float(source) for source in sources)) > 1.0e-12:
        raise AssertionError("local-edge mapping is not conservative")
    return {"values": values, "weights": weights}
def edge_scalar_to_cells(edge_values, cell_edges):
    """Average three edge-centered intensive values onto each triangle with explicit weights."""
    values, weights = [], []
    for edges in cell_edges:
        ids = tuple(edges)
        if len(ids) != 3 or len(set(ids)) != 3 or any(edge not in edge_values for edge in ids):
            raise ValueError("cell needs three distinct edge values")
        weight = {edge: 1.0 / 3.0 for edge in ids}
        weights.append(weight)
        values.append(sum(float(edge_values[edge]) * weight[edge] for edge in ids))
    return {"values": values, "weights": weights}