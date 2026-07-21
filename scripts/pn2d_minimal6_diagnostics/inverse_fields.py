"""Diagnostic-only recovery of electrostatic gradients and field candidates."""

from __future__ import annotations

from dataclasses import dataclass
import math
import statistics
from typing import Iterable, Mapping, Sequence

try:
    from .inverse_contracts import (
        AcceptanceThresholds,
        Identifiability,
        Observation,
        SampleStatus,
        SupportKind,
    )
except ImportError:
    from scripts.pn2d_minimal6_diagnostics.inverse_contracts import (  # type: ignore
        AcceptanceThresholds,
        Identifiability,
        Observation,
        SampleStatus,
        SupportKind,
    )


Vector = tuple[float, float]
Point = tuple[float, float]


@dataclass(frozen=True)
class VectorErrorResult:
    magnitude_status: SampleStatus
    direction_status: SampleStatus
    relative_magnitude_error: float | None
    angle_deg: float | None


@dataclass(frozen=True)
class FieldCandidateSample:
    candidate: str
    solver: str
    topology: str
    bias_V: float
    support_kind: SupportKind
    support_id: int | str
    candidate_value: Vector | None
    reference_value: Vector | None
    unit_si: str
    coordinate_frame: str
    orientation: str
    candidate_support_transform: str
    reference_support_transform: str
    error: VectorErrorResult


@dataclass(frozen=True)
class FieldCandidateResult:
    candidate: str
    solver: str
    topology: str
    bias_V: float
    support_kind: SupportKind
    samples: tuple[FieldCandidateSample, ...]
    magnitude_valid_count: int
    direction_valid_count: int
    median_relative_magnitude_error: float | None
    median_angle_deg: float | None
    classification: Identifiability


def _stable_key(value: object) -> tuple[str, str]:
    return type(value).__name__, str(value)


def _items(values):
    if isinstance(values, Mapping):
        return tuple(sorted(values.items(), key=lambda item: _stable_key(item[0])))
    return tuple(enumerate(values))


def _finite_point(value, label: str) -> Point:
    try:
        if len(value) != 2:
            raise ValueError
        result = float(value[0]), float(value[1])
    except (TypeError, ValueError, IndexError) as error:
        raise ValueError(f"{label} must have two numeric components") from error
    if not all(math.isfinite(component) for component in result):
        raise ValueError(f"{label} must be finite")
    return result


def _finite_vector(value, label: str) -> Vector:
    return _finite_point(value, label)


def _edge_geometry(start: Point, end: Point) -> tuple[float, float, float]:
    dx, dy = end[0] - start[0], end[1] - start[1]
    if not math.isfinite(dx) or not math.isfinite(dy):
        raise ValueError("edge has nonfinite geometry")
    length = math.hypot(dx, dy)
    if not math.isfinite(length):
        raise ValueError("edge has nonfinite geometry")
    if length == 0.0:
        raise ValueError("edge has zero length")
    return dx, dy, length


def _triangle_area(points: Sequence[Point]) -> float:
    if len(points) != 3:
        raise ValueError("triangle needs exactly three points")
    (x0, y0), (x1, y1), (x2, y2) = points
    twice_area = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)
    if not math.isfinite(twice_area) or abs(twice_area) <= 1.0e-300:
        raise ValueError("degenerate triangle")
    return 0.5 * abs(twice_area)


def triangle_gradient(points, values) -> Vector:
    """Return the exact P1 scalar gradient on one non-degenerate triangle."""
    if len(points) != 3 or len(values) != 3:
        raise ValueError("triangle gradient needs three points and values")
    (x0, y0), (x1, y1), (x2, y2) = (
        _finite_point(point, "triangle point") for point in points
    )
    try:
        f0, f1, f2 = (float(value) for value in values)
    except (TypeError, ValueError) as error:
        raise ValueError("triangle values must be numeric") from error
    if not all(math.isfinite(value) for value in (f0, f1, f2)):
        raise ValueError("triangle values must be finite")
    det = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)
    if not math.isfinite(det):
        raise ValueError("nonfinite triangle geometry")
    if abs(det) <= 1e-300:
        raise ValueError("degenerate triangle")
    gx = ((f1 - f0) * (y2 - y0) - (f2 - f0) * (y1 - y0)) / det
    gy = ((x1 - x0) * (f2 - f0) - (x2 - x0) * (f1 - f0)) / det
    if not math.isfinite(gx) or not math.isfinite(gy):
        raise ValueError("nonfinite triangle gradient")
    return gx, gy


def _cell_geometry(cells, coordinates):
    cell_items = _items(cells)
    coordinate_map = dict(_items(coordinates))
    normalized = {}
    areas = {}
    for cell_id, nodes_value in cell_items:
        nodes = tuple(nodes_value)
        if len(nodes) != 3 or len(set(nodes)) != 3:
            raise ValueError("cell needs three distinct nodes")
        if any(node not in coordinate_map for node in nodes):
            raise ValueError("cell coordinate is missing")
        points = tuple(_finite_point(coordinate_map[node], "node coordinate") for node in nodes)
        normalized[cell_id] = nodes
        areas[cell_id] = _triangle_area(points)
    return normalized, areas, coordinate_map


def _weighted_vector(cell_ids, cell_vectors, areas):
    total = sum(areas[cell_id] for cell_id in cell_ids)
    if not math.isfinite(total) or total <= 0.0:
        raise ValueError("support has no positive adjacent-cell area")
    x = sum(areas[cell_id] * cell_vectors[cell_id][0] for cell_id in cell_ids) / total
    y = sum(areas[cell_id] * cell_vectors[cell_id][1] for cell_id in cell_ids) / total
    return (x, y), tuple((cell_id, areas[cell_id] / total) for cell_id in cell_ids)


def cell_to_node_vectors(cell_vectors, cells, coordinates):
    """Area-weight cell vectors onto incident nodes, retaining explicit weights."""
    normalized_cells, areas, _ = _cell_geometry(cells, coordinates)
    vectors = dict(_items(cell_vectors))
    if set(vectors) != set(normalized_cells):
        raise ValueError("cell vectors and topology must have identical support")
    vectors = {cell_id: _finite_vector(vector, "cell vector")
               for cell_id, vector in vectors.items()}
    node_ids = sorted({node for nodes in normalized_cells.values() for node in nodes}, key=_stable_key)
    values, weights = {}, {}
    for node_id in node_ids:
        adjacent = tuple(cell_id for cell_id in sorted(normalized_cells, key=_stable_key)
                         if node_id in normalized_cells[cell_id])
        values[node_id], weights[node_id] = _weighted_vector(adjacent, vectors, areas)
    return {"values": values, "weights": weights}


def cell_to_edge_vectors(cell_vectors, cells, edges, coordinates):
    """Area-weight cell vectors onto explicitly directed adjacent edges."""
    normalized_cells, areas, coordinate_map = _cell_geometry(cells, coordinates)
    vectors = dict(_items(cell_vectors))
    if set(vectors) != set(normalized_cells):
        raise ValueError("cell vectors and topology must have identical support")
    vectors = {cell_id: _finite_vector(vector, "cell vector")
               for cell_id, vector in vectors.items()}
    values, weights = {}, {}
    for edge_id, edge_value in _items(edges):
        edge = tuple(edge_value)
        if len(edge) != 2 or edge[0] == edge[1]:
            raise ValueError("edge needs two distinct directed endpoints")
        if any(node not in coordinate_map for node in edge):
            raise ValueError("edge coordinate is missing")
        start, end = (_finite_point(coordinate_map[node], "edge coordinate") for node in edge)
        _edge_geometry(start, end)
        adjacent = tuple(cell_id for cell_id in sorted(normalized_cells, key=_stable_key)
                         if edge[0] in normalized_cells[cell_id] and edge[1] in normalized_cells[cell_id])
        if not adjacent:
            raise ValueError("edge has no adjacent cell")
        values[edge_id], weights[edge_id] = _weighted_vector(adjacent, vectors, areas)
    return {"values": values, "weights": weights}


def edge_scalar_difference(start_value, end_value, start, end) -> float:
    """Return ``(end_value - start_value) / length`` along the declared edge."""
    try:
        first, second = float(start_value), float(end_value)
    except (TypeError, ValueError) as error:
        raise ValueError("edge scalar values must be numeric") from error
    if not math.isfinite(first) or not math.isfinite(second):
        raise ValueError("edge scalar values must be finite")
    x0, y0 = _finite_point(start, "edge start")
    x1, y1 = _finite_point(end, "edge end")
    _, _, length = _edge_geometry((x0, y0), (x1, y1))
    result = (second - first) / length
    if not math.isfinite(result):
        raise ValueError("nonfinite edge scalar difference")
    return result


def mirror_vector(vector) -> Vector:
    """Transform a global Cartesian vector through the x-coordinate mirror."""
    x, y = _finite_vector(vector, "vector")
    return -x, y


def vector_error(candidate, reference, *, reference_floor: float) -> VectorErrorResult:
    """Compare vector magnitude and direction without component-wise ratios."""
    try:
        floor = float(reference_floor)
    except (TypeError, ValueError) as error:
        raise ValueError("reference floor must be finite and non-negative") from error
    if not math.isfinite(floor) or floor < 0.0:
        raise ValueError("reference floor must be finite and non-negative")
    if candidate is None or reference is None:
        return VectorErrorResult(SampleStatus.MISSING_FIELD, SampleStatus.MISSING_FIELD,
                                 None, None)
    try:
        if len(candidate) != 2 or len(reference) != 2:
            raise ValueError
        candidate_vector = float(candidate[0]), float(candidate[1])
        reference_vector = float(reference[0]), float(reference[1])
    except (TypeError, ValueError, IndexError):
        return VectorErrorResult(SampleStatus.NONFINITE, SampleStatus.NONFINITE,
                                 None, None)
    if not all(math.isfinite(value) for value in candidate_vector + reference_vector):
        return VectorErrorResult(SampleStatus.NONFINITE, SampleStatus.NONFINITE,
                                 None, None)
    candidate_magnitude = math.hypot(*candidate_vector)
    reference_magnitude = math.hypot(*reference_vector)
    if not math.isfinite(candidate_magnitude) or not math.isfinite(reference_magnitude):
        return VectorErrorResult(SampleStatus.NONFINITE, SampleStatus.NONFINITE,
                                 None, None)
    if reference_magnitude == 0.0:
        return VectorErrorResult(SampleStatus.GEOMETRIC_ZERO,
                                 SampleStatus.DIRECTION_UNDEFINED, None, None)
    if reference_magnitude <= floor:
        return VectorErrorResult(SampleStatus.BELOW_FLOOR,
                                 SampleStatus.DIRECTION_UNDEFINED, None, None)
    relative_error = abs(candidate_magnitude - reference_magnitude) / reference_magnitude
    if not math.isfinite(relative_error):
        return VectorErrorResult(SampleStatus.NONFINITE, SampleStatus.NONFINITE,
                                 None, None)
    if candidate_magnitude == 0.0:
        return VectorErrorResult(SampleStatus.VALID, SampleStatus.DIRECTION_UNDEFINED,
                                 relative_error, None)
    candidate_unit = (candidate_vector[0] / candidate_magnitude,
                      candidate_vector[1] / candidate_magnitude)
    reference_unit = (reference_vector[0] / reference_magnitude,
                      reference_vector[1] / reference_magnitude)
    if not all(math.isfinite(value) for value in candidate_unit + reference_unit):
        return VectorErrorResult(SampleStatus.VALID, SampleStatus.NONFINITE,
                                 relative_error, None)
    cosine = (candidate_unit[0] * reference_unit[0]
              + candidate_unit[1] * reference_unit[1])
    if not math.isfinite(cosine):
        return VectorErrorResult(SampleStatus.VALID, SampleStatus.NONFINITE,
                                 relative_error, None)
    cosine = max(-1.0, min(1.0, cosine))
    angle = math.degrees(math.acos(cosine))
    if not math.isfinite(angle):
        return VectorErrorResult(SampleStatus.VALID, SampleStatus.NONFINITE,
                                 relative_error, None)
    return VectorErrorResult(SampleStatus.VALID, SampleStatus.VALID,
                             relative_error, angle)


def _topology(mesh, topology):
    if not isinstance(mesh, Mapping):
        raise ValueError("topology mesh must be a mapping")
    selected = mesh
    if "triangles" not in selected or "edges" not in selected:
        selected = mesh.get(topology)
    if not isinstance(selected, Mapping) or "triangles" not in selected or "edges" not in selected:
        raise ValueError(f"topology mesh is missing {topology}")
    cells = {cell_id: tuple(str(node) for node in nodes)
             for cell_id, nodes in _items(selected["triangles"])}
    edges = {edge_id: tuple(str(node) for node in nodes)
             for edge_id, nodes in _items(selected["edges"])}
    return cells, edges


def _usable(observation: Observation | None, expected_unit: str) -> tuple[float | None, SampleStatus]:
    if observation is None:
        return None, SampleStatus.MISSING_FIELD
    if observation.status is not SampleStatus.VALID:
        return None, observation.status
    if observation.value_si is None:
        return None, SampleStatus.MISSING_FIELD
    if observation.unit_si != expected_unit:
        return None, SampleStatus.INVALID_UNIT
    try:
        value = float(observation.value_si)
    except (TypeError, ValueError):
        return None, SampleStatus.NONFINITE
    if not math.isfinite(value) or observation.status is SampleStatus.NONFINITE:
        return None, SampleStatus.NONFINITE
    return value, SampleStatus.VALID


def _invalid_error(status: SampleStatus) -> VectorErrorResult:
    return VectorErrorResult(status, status, None, None)


def _sample(candidate, state, support_kind, support_id, candidate_value, reference_value,
            frame, orientation, candidate_transform, reference_transform, reference_floor,
            invalid_status=None):
    error = (_invalid_error(invalid_status) if invalid_status is not None
             else vector_error(candidate_value, reference_value, reference_floor=reference_floor))
    return FieldCandidateSample(
        candidate=candidate, solver=state[0], topology=state[1], bias_V=state[2],
        support_kind=support_kind, support_id=support_id,
        candidate_value=candidate_value, reference_value=reference_value,
        unit_si="V/m", coordinate_frame=frame, orientation=orientation,
        candidate_support_transform=candidate_transform,
        reference_support_transform=reference_transform, error=error,
    )


def _summarize(candidate, state, support_kind, samples, thresholds):
    magnitude = [sample.error.relative_magnitude_error for sample in samples
                 if sample.error.magnitude_status is SampleStatus.VALID
                 and sample.error.relative_magnitude_error is not None]
    angles = [sample.error.angle_deg for sample in samples
              if sample.error.direction_status is SampleStatus.VALID
              and sample.error.angle_deg is not None]
    median_magnitude = statistics.median(magnitude) if magnitude else None
    median_angle = statistics.median(angles) if angles else None
    if median_magnitude is None or median_angle is None:
        classification = Identifiability.INSUFFICIENT_DATA
    elif (median_magnitude <= thresholds.field_median_relative
          and median_angle <= thresholds.field_median_angle_deg):
        classification = Identifiability.IDENTIFIED
    else:
        classification = Identifiability.REJECTED
    return FieldCandidateResult(
        candidate=candidate, solver=state[0], topology=state[1], bias_V=state[2],
        support_kind=support_kind, samples=tuple(samples),
        magnitude_valid_count=len(magnitude), direction_valid_count=len(angles),
        median_relative_magnitude_error=median_magnitude,
        median_angle_deg=median_angle, classification=classification,
    )


def evaluate_field_candidates(
    observations: Iterable[Observation],
    mesh: Mapping,
    *,
    reference_floor: float,
    thresholds: AcceptanceThresholds | None = None,
) -> tuple[FieldCandidateResult, ...]:
    """Evaluate four fixed electrostatic-field candidates on compatible support.

    Coordinates and physical fields must already be canonical SI node
    observations. Triangle and edge identities come from ``mesh``. Direct
    nodal electric fields are explicitly averaged to cell/edge support before
    comparison; no cross-support comparison is implicit.
    """
    canonical_limits = AcceptanceThresholds()
    if thresholds is None:
        limits = canonical_limits
    elif not isinstance(thresholds, AcceptanceThresholds):
        raise ValueError("thresholds must use the fixed inverse-audit contract")
    elif (thresholds.field_median_relative != canonical_limits.field_median_relative
          or thresholds.field_median_angle_deg != canonical_limits.field_median_angle_deg):
        raise ValueError("field thresholds are immutable at 2 percent and 1 degree")
    else:
        limits = thresholds
    rows = tuple(observations)
    relevant = [row for row in rows if row.support_kind is SupportKind.NODE
                and row.quantity in {"coordinate", "ElectrostaticPotential", "ElectricField"}]
    seen = set()
    groups = {}
    for row in relevant:
        if row.key in seen:
            raise ValueError("duplicate field observation key")
        seen.add(row.key)
        state = row.solver, row.topology, float(row.bias_V)
        groups.setdefault(state, []).append(row)

    results = []
    for state in sorted(groups, key=lambda item: (item[0], item[1], item[2])):
        state_rows = groups[state]
        frames = {row.coordinate_frame for row in state_rows}
        orientations = {row.orientation for row in state_rows}
        if len(frames) != 1 or len(orientations) != 1:
            raise ValueError("field observations have incompatible coordinate frames or orientations")
        frame, orientation = next(iter(frames)), next(iter(orientations))
        index = {(str(row.support_id), row.quantity, row.component): row for row in state_rows}
        node_ids = sorted({key[0] for key in index if key[1] == "coordinate"}, key=_stable_key)
        coordinates = {}
        for node_id in node_ids:
            x, x_status = _usable(index.get((node_id, "coordinate", "x")), "m")
            y, y_status = _usable(index.get((node_id, "coordinate", "y")), "m")
            if x_status is not SampleStatus.VALID or y_status is not SampleStatus.VALID:
                raise ValueError("canonical SI node coordinates are required")
            coordinates[node_id] = (x, y)

        cells, edges = _topology(mesh, state[1])
        normalized_cells, areas, _ = _cell_geometry(cells, coordinates)
        for edge_id, edge in edges.items():
            if len(edge) != 2 or edge[0] == edge[1] or any(node not in coordinates for node in edge):
                raise ValueError("invalid directed topology edge")
            if not any(edge[0] in nodes and edge[1] in nodes for nodes in normalized_cells.values()):
                raise ValueError("edge has no adjacent cell")
            start, end = coordinates[edge[0]], coordinates[edge[1]]
            _edge_geometry(start, end)

        potentials, potential_status = {}, {}
        fields, field_status = {}, {}
        for node_id in node_ids:
            potentials[node_id], potential_status[node_id] = _usable(
                index.get((node_id, "ElectrostaticPotential", "component0")), "V")
            ex, ex_status = _usable(index.get((node_id, "ElectricField", "component0")), "V/m")
            ey, ey_status = _usable(index.get((node_id, "ElectricField", "component1")), "V/m")
            status = ex_status if ex_status is not SampleStatus.VALID else ey_status
            fields[node_id] = None if status is not SampleStatus.VALID else (ex, ey)
            field_status[node_id] = status

        gradients, gradient_status = {}, {}
        for cell_id, nodes in normalized_cells.items():
            statuses = [potential_status.get(node, SampleStatus.MISSING_FIELD) for node in nodes]
            invalid = next((status for status in statuses if status is not SampleStatus.VALID), None)
            if invalid is not None:
                gradients[cell_id], gradient_status[cell_id] = None, invalid
                continue
            gradient = triangle_gradient(tuple(coordinates[node] for node in nodes),
                                         tuple(potentials[node] for node in nodes))
            gradients[cell_id], gradient_status[cell_id] = (-gradient[0], -gradient[1]), SampleStatus.VALID

        triangle_samples = []
        for cell_id in sorted(normalized_cells, key=_stable_key):
            nodes = normalized_cells[cell_id]
            reference_invalid = next((field_status[node] for node in nodes
                                      if field_status[node] is not SampleStatus.VALID), None)
            reference = None if reference_invalid else (
                sum(fields[node][0] for node in nodes) / 3.0,
                sum(fields[node][1] for node in nodes) / 3.0,
            )
            invalid = gradient_status[cell_id] if gradient_status[cell_id] is not SampleStatus.VALID else reference_invalid
            triangle_samples.append(_sample(
                "triangle_minus_grad_psi", state, SupportKind.CELL, cell_id,
                gradients[cell_id], reference, frame, orientation,
                "P1_triangle_minus_gradient", "node_to_cell_equal_P1_weights",
                reference_floor, invalid,
            ))

        node_samples = []
        for node_id in node_ids:
            adjacent = tuple(cell_id for cell_id in sorted(normalized_cells, key=_stable_key)
                             if node_id in normalized_cells[cell_id])
            invalid = next((gradient_status[cell_id] for cell_id in adjacent
                            if gradient_status[cell_id] is not SampleStatus.VALID), None)
            candidate_value = None
            if invalid is None:
                candidate_value, _ = _weighted_vector(adjacent, gradients, areas)
            reference_invalid = field_status[node_id]
            if invalid is None and reference_invalid is not SampleStatus.VALID:
                invalid = reference_invalid
            node_samples.append(_sample(
                "node_area_weighted_minus_grad_psi", state, SupportKind.NODE, node_id,
                candidate_value, fields[node_id], frame, orientation,
                "cell_to_node_area_weighted", "native_node_vector",
                reference_floor, invalid,
            ))

        edge_vector_samples, edge_scalar_samples = [], []
        for edge_id in sorted(edges, key=_stable_key):
            start_id, end_id = edges[edge_id]
            adjacent = tuple(cell_id for cell_id in sorted(normalized_cells, key=_stable_key)
                             if start_id in normalized_cells[cell_id] and end_id in normalized_cells[cell_id])
            invalid = next((gradient_status[cell_id] for cell_id in adjacent
                            if gradient_status[cell_id] is not SampleStatus.VALID), None)
            candidate_value = None
            if invalid is None:
                candidate_value, _ = _weighted_vector(adjacent, gradients, areas)
            reference_invalid = next((field_status[node] for node in (start_id, end_id)
                                      if field_status[node] is not SampleStatus.VALID), None)
            reference = None if reference_invalid else (
                0.5 * (fields[start_id][0] + fields[end_id][0]),
                0.5 * (fields[start_id][1] + fields[end_id][1]),
            )
            vector_invalid = invalid if invalid is not None else reference_invalid
            edge_vector_samples.append(_sample(
                "edge_area_weighted_minus_grad_psi", state, SupportKind.EDGE, edge_id,
                candidate_value, reference, frame, orientation,
                "cell_to_edge_adjacent_area_weighted", "node_to_edge_endpoint_average",
                reference_floor, vector_invalid,
            ))

            scalar_invalid = next((potential_status[node] for node in (start_id, end_id)
                                   if potential_status[node] is not SampleStatus.VALID), None)
            if scalar_invalid is None and reference_invalid is not None:
                scalar_invalid = reference_invalid
            signed_vector = projected_reference = None
            if scalar_invalid is None:
                start, end = coordinates[start_id], coordinates[end_id]
                dx, dy, length = _edge_geometry(start, end)
                tangent = dx / length, dy / length
                signed = -edge_scalar_difference(potentials[start_id], potentials[end_id], start, end)
                reference_projection = reference[0] * tangent[0] + reference[1] * tangent[1]
                signed_vector = signed * tangent[0], signed * tangent[1]
                projected_reference = reference_projection * tangent[0], reference_projection * tangent[1]
            edge_scalar_samples.append(_sample(
                "signed_edge_minus_delta_psi_over_h", state, SupportKind.EDGE, edge_id,
                signed_vector, projected_reference, frame, orientation,
                "directed_edge_minus_delta_psi_over_length",
                "node_to_edge_average_then_directed_tangent_projection",
                reference_floor, scalar_invalid,
            ))

        result_groups = (
            ("triangle_minus_grad_psi", SupportKind.CELL, triangle_samples),
            ("node_area_weighted_minus_grad_psi", SupportKind.NODE, node_samples),
            ("edge_area_weighted_minus_grad_psi", SupportKind.EDGE, edge_vector_samples),
            ("signed_edge_minus_delta_psi_over_h", SupportKind.EDGE, edge_scalar_samples),
        )
        results.extend(_summarize(name, state, support, samples, limits)
                       for name, support, samples in result_groups)

    return tuple(sorted(results, key=lambda item: (
        item.solver, item.topology, item.bias_V, item.candidate,
    )))
