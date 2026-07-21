"""Diagnostic-only recovery of quasi-Fermi and current-density semantics.

All scalar state is consumed in canonical SI units.  Directed edge quantities
retain their declared start-to-end orientation; nodal, edge, and cell support
are transformed explicitly before comparison.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
import statistics
from typing import Iterable, Mapping

try:
    from .inverse_contracts import (
        AcceptanceThresholds,
        Identifiability,
        Observation,
        SampleStatus,
        SupportKind,
    )
    from .inverse_fields import (
        _usable,
        cell_to_edge_vectors,
        cell_to_node_vectors,
        edge_scalar_difference,
        triangle_gradient,
    )
    from .inverse_inputs import DISCOVERY_KEYS
except ImportError:
    from scripts.pn2d_minimal6_diagnostics.inverse_contracts import (  # type: ignore
        AcceptanceThresholds,
        Identifiability,
        Observation,
        SampleStatus,
        SupportKind,
    )
    from scripts.pn2d_minimal6_diagnostics.inverse_fields import (  # type: ignore
        _usable,
        cell_to_edge_vectors,
        cell_to_node_vectors,
        edge_scalar_difference,
        triangle_gradient,
    )
    from scripts.pn2d_minimal6_diagnostics.inverse_inputs import DISCOVERY_KEYS  # type: ignore


Vector = tuple[float, float]
Point = tuple[float, float]
ELEMENTARY_CHARGE_C = 1.602176634e-19


@dataclass(frozen=True)
class TransportVectorError:
    magnitude_status: SampleStatus
    direction_status: SampleStatus
    abs_log10_error: float | None
    angle_deg: float | None


@dataclass(frozen=True)
class TransportMetricSummary:
    errors: tuple[TransportVectorError, ...]
    valid_count: int
    direction_valid_count: int
    median_abs_log10_error: float | None
    p95_abs_log10_error: float | None
    median_angle_deg: float | None
    classification: Identifiability


@dataclass(frozen=True)
class TransportConfoundingRecord:
    candidate: str
    carrier: str
    topology: str
    bias_V: float
    support_kind: SupportKind
    support_id: int | str
    status: SampleStatus
    classification: Identifiability
    missing_inputs: tuple[str, ...]
    observable: str | None
    observable_value: Vector | None
    observable_unit_si: str | None


@dataclass(frozen=True)
class TransportCandidateSample:
    candidate: str
    solver: str
    topology: str
    bias_V: float
    carrier: str
    split: str
    support_kind: SupportKind
    support_id: int | str
    candidate_value: Vector | None
    reference_value: Vector | None
    unit_si: str
    coordinate_frame: str
    orientation: str
    candidate_support_transform: str
    reference_support_transform: str
    error: TransportVectorError


@dataclass(frozen=True)
class TransportCandidateResult:
    candidate: str
    solver: str
    topology: str
    bias_V: float
    carrier: str
    split: str
    support_kind: SupportKind
    unit_si: str
    samples: tuple[TransportCandidateSample, ...]
    valid_count: int
    direction_valid_count: int
    median_abs_log10_error: float | None
    p95_abs_log10_error: float | None
    median_angle_deg: float | None
    classification: Identifiability
    confoundings: tuple[TransportConfoundingRecord, ...]


def _carrier_sign(carrier: str) -> float:
    if carrier == "electron":
        return -1.0
    if carrier == "hole":
        return 1.0
    raise ValueError("carrier must be 'electron' or 'hole'")


def _finite_scalar(value, label: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{label} must be finite") from error
    if not math.isfinite(result):
        raise ValueError(f"{label} must be finite")
    return result


def _finite_vector(value, label: str) -> Vector:
    try:
        if len(value) != 2:
            raise ValueError
        result = float(value[0]), float(value[1])
    except (TypeError, ValueError, IndexError) as error:
        raise ValueError(f"{label} must have two finite components") from error
    if not all(math.isfinite(component) for component in result):
        raise ValueError(f"{label} must be finite")
    return result


def qf_current_density(
    carrier,
    density,
    mobility,
    gradient,
    *,
    q=ELEMENTARY_CHARGE_C,
) -> Vector:
    """Return conventional current density from a quasi-Fermi gradient.

    Electron and hole quasi-Fermi potentials have different conventional-
    current signs: ``Jn = -q mu_n n grad(phi_n)`` and
    ``Jp = +q mu_p p grad(phi_p)``.
    """
    sign = _carrier_sign(carrier)
    density_value = _finite_scalar(density, "density")
    mobility_value = _finite_scalar(mobility, "mobility")
    charge = _finite_scalar(q, "q")
    gradient_value = _finite_vector(gradient, "gradient")
    if density_value < 0.0 or mobility_value < 0.0:
        raise ValueError("density and mobility must be non-negative")
    if charge <= 0.0:
        raise ValueError("q must be positive")
    factor = sign * charge * mobility_value * density_value
    result = factor * gradient_value[0], factor * gradient_value[1]
    if not all(math.isfinite(component) for component in result):
        raise ValueError("current density must be finite")
    return result


def current_inverted_qf_gradient(
    carrier,
    density,
    mobility,
    current,
    *,
    q=ELEMENTARY_CHARGE_C,
    density_floor=0.0,
    mobility_floor=0.0,
) -> Vector:
    """Invert a current vector into a diagnostic quasi-Fermi gradient.

    The inversion is diagnostic evidence only.  It is not a production
    transport formula and refuses every zero/floor denominator.
    """
    sign = _carrier_sign(carrier)
    density_value = _finite_scalar(density, "density")
    mobility_value = _finite_scalar(mobility, "mobility")
    charge = _finite_scalar(q, "q")
    density_limit = _finite_scalar(density_floor, "density floor")
    mobility_limit = _finite_scalar(mobility_floor, "mobility floor")
    current_value = _finite_vector(current, "current")
    if density_limit < 0.0 or mobility_limit < 0.0:
        raise ValueError("floors must be non-negative")
    if density_value <= density_limit:
        raise ValueError("density is at or below floor")
    if mobility_value <= mobility_limit:
        raise ValueError("mobility is at or below floor")
    if charge <= 0.0:
        raise ValueError("q must be positive")
    denominator = sign * charge * mobility_value * density_value
    if not math.isfinite(denominator) or denominator == 0.0:
        raise ValueError("current inversion denominator is invalid")
    result = current_value[0] / denominator, current_value[1] / denominator
    if not all(math.isfinite(component) for component in result):
        raise ValueError("inverted gradient must be finite")
    return result


def _edge_tangent(start, end) -> Vector:
    start_value = _finite_vector(start, "edge start")
    end_value = _finite_vector(end, "edge end")
    dx = end_value[0] - start_value[0]
    dy = end_value[1] - start_value[1]
    if not math.isfinite(dx) or not math.isfinite(dy):
        raise ValueError("edge geometry must be finite")
    length = math.hypot(dx, dy)
    if not math.isfinite(length):
        raise ValueError("edge geometry must be finite")
    if length == 0.0:
        raise ValueError("edge has zero length")
    return dx / length, dy / length


def project_vector_to_edge(vector, start, end) -> float:
    """Project a Cartesian vector onto the directed start-to-end tangent."""
    vector_value = _finite_vector(vector, "vector")
    tangent = _edge_tangent(start, end)
    result = vector_value[0] * tangent[0] + vector_value[1] * tangent[1]
    if not math.isfinite(result):
        raise ValueError("edge projection must be finite")
    return result


def reconstruct_edge_vector(signed_component, start, end) -> Vector:
    """Reconstruct a tangent vector from a signed directed-edge component."""
    component = _finite_scalar(signed_component, "signed edge component")
    tangent = _edge_tangent(start, end)
    result = component * tangent[0], component * tangent[1]
    if not all(math.isfinite(value) for value in result):
        raise ValueError("reconstructed edge vector must be finite")
    return result


def _transport_vector_error(candidate, reference, *, reference_floor: float) -> TransportVectorError:
    floor = _finite_scalar(reference_floor, "reference floor")
    if floor < 0.0:
        raise ValueError("reference floor must be non-negative")
    if candidate is None or reference is None:
        return TransportVectorError(
            SampleStatus.MISSING_FIELD, SampleStatus.MISSING_FIELD, None, None
        )
    try:
        candidate_value = _finite_vector(candidate, "candidate")
        reference_value = _finite_vector(reference, "reference")
    except ValueError:
        return TransportVectorError(
            SampleStatus.NONFINITE, SampleStatus.NONFINITE, None, None
        )
    candidate_magnitude = math.hypot(*candidate_value)
    reference_magnitude = math.hypot(*reference_value)
    if not math.isfinite(candidate_magnitude) or not math.isfinite(reference_magnitude):
        return TransportVectorError(
            SampleStatus.NONFINITE, SampleStatus.NONFINITE, None, None
        )
    if reference_magnitude == 0.0:
        return TransportVectorError(
            SampleStatus.GEOMETRIC_ZERO, SampleStatus.DIRECTION_UNDEFINED, None, None
        )
    if reference_magnitude <= floor or candidate_magnitude <= floor:
        return TransportVectorError(
            SampleStatus.BELOW_FLOOR, SampleStatus.DIRECTION_UNDEFINED, None, None
        )
    log_error = abs(math.log10(candidate_magnitude / reference_magnitude))
    cosine = (
        candidate_value[0] * reference_value[0]
        + candidate_value[1] * reference_value[1]
    ) / (candidate_magnitude * reference_magnitude)
    if not math.isfinite(log_error) or not math.isfinite(cosine):
        return TransportVectorError(
            SampleStatus.NONFINITE, SampleStatus.NONFINITE, None, None
        )
    angle = math.degrees(math.acos(max(-1.0, min(1.0, cosine))))
    return TransportVectorError(SampleStatus.VALID, SampleStatus.VALID, log_error, angle)


def _fixed_thresholds(thresholds: AcceptanceThresholds | None) -> AcceptanceThresholds:
    canonical = AcceptanceThresholds()
    if thresholds is None:
        return canonical
    if not isinstance(thresholds, AcceptanceThresholds):
        raise ValueError("thresholds must use the fixed inverse-audit contract")
    if (
        thresholds.gradient_median_abs_dex != canonical.gradient_median_abs_dex
        or thresholds.gradient_p95_abs_dex != canonical.gradient_p95_abs_dex
        or thresholds.gradient_median_angle_deg != canonical.gradient_median_angle_deg
    ):
        raise ValueError("transport thresholds are immutable at 0.1/0.3 dex and 5 degrees")
    return thresholds


def _nearest_rank_p95(values: list[float]) -> float:
    ordered = sorted(values)
    return ordered[max(0, math.ceil(0.95 * len(ordered)) - 1)]


def summarize_transport_errors(
    errors: Iterable[TransportVectorError],
    *,
    thresholds: AcceptanceThresholds | None = None,
) -> TransportMetricSummary:
    """Summarize fixed transport magnitude/direction gates without fitting."""
    limits = _fixed_thresholds(thresholds)
    rows = tuple(errors)
    magnitudes = [
        row.abs_log10_error for row in rows
        if row.magnitude_status is SampleStatus.VALID and row.abs_log10_error is not None
    ]
    angles = [
        row.angle_deg for row in rows
        if row.direction_status is SampleStatus.VALID and row.angle_deg is not None
    ]
    median = statistics.median(magnitudes) if magnitudes else None
    p95 = _nearest_rank_p95(magnitudes) if magnitudes else None
    median_angle = statistics.median(angles) if angles else None
    if median is None or p95 is None or median_angle is None:
        classification = Identifiability.INSUFFICIENT_DATA
    elif (
        median <= limits.gradient_median_abs_dex
        and p95 <= limits.gradient_p95_abs_dex
        and median_angle <= limits.gradient_median_angle_deg
    ):
        classification = Identifiability.IDENTIFIED
    else:
        classification = Identifiability.REJECTED
    return TransportMetricSummary(
        rows, len(magnitudes), len(angles), median, p95, median_angle, classification
    )


def _stable_key(value: object) -> tuple[str, str]:
    return type(value).__name__, str(value)


def _items(values):
    if isinstance(values, Mapping):
        return tuple(sorted(values.items(), key=lambda item: _stable_key(item[0])))
    return tuple(enumerate(values))


def _topology(mesh, topology):
    if not isinstance(mesh, Mapping):
        raise ValueError("topology mesh must be a mapping")
    selected = mesh
    if "triangles" not in selected or "edges" not in selected:
        selected = mesh.get(topology)
    if not isinstance(selected, Mapping) or "triangles" not in selected or "edges" not in selected:
        raise ValueError(f"topology mesh is missing {topology}")
    cells = {
        cell_id: tuple(str(node) for node in nodes)
        for cell_id, nodes in _items(selected["triangles"])
    }
    edges = {
        edge_id: tuple(str(node) for node in nodes)
        for edge_id, nodes in _items(selected["edges"])
    }
    return cells, edges


def _average(values: Iterable[float]) -> float:
    rows = tuple(values)
    return sum(rows) / len(rows)


def _average_vectors(values: Iterable[Vector]) -> Vector:
    rows = tuple(values)
    return _average(value[0] for value in rows), _average(value[1] for value in rows)


def _status_value(index, node, quantity, component, unit):
    return _usable(index.get((node, quantity, component)), unit)


def _first_invalid(items):
    return next((status for status in items if status is not SampleStatus.VALID), None)


def _quantity_contract(carrier):
    prefix = "e" if carrier == "electron" else "h"
    return (
        f"{prefix}QuasiFermiPotential",
        f"{prefix}Density",
        f"{prefix}Mobility",
        f"{prefix}CurrentDensity",
    )


def _bernoulli(value: float) -> float:
    if abs(value) < 1.0e-8:
        return 1.0 - 0.5 * value + value * value / 12.0
    if value > 50.0:
        return value * math.exp(-value)
    if value < -50.0:
        return -value
    return value / math.expm1(value)


def _sg_current(carrier, density0, density1, dpsi, thermal_voltage, mobility, length, q):
    u = dpsi / thermal_voltage
    if carrier == "electron":
        difference = _bernoulli(u) * density1 - _bernoulli(-u) * density0
    else:
        difference = _bernoulli(-u) * density1 - _bernoulli(u) * density0
    result = q * mobility * thermal_voltage * difference / length
    if not math.isfinite(result):
        raise ValueError("SG current must be finite")
    return result


def _confounding(
    candidate,
    carrier,
    state,
    support_kind,
    support_id,
    status,
    missing_inputs,
    *,
    observable=None,
    observable_value=None,
    observable_unit=None,
):
    classification = (
        Identifiability.CONFOUNDED
        if observable_value is not None
        else Identifiability.INSUFFICIENT_DATA
    )
    return TransportConfoundingRecord(
        candidate, carrier, state[1], state[2], support_kind, support_id,
        status, classification, tuple(missing_inputs), observable,
        observable_value, observable_unit,
    )


def _sample(
    candidate,
    carrier,
    state,
    split,
    support_kind,
    support_id,
    candidate_value,
    reference_value,
    unit,
    frame,
    orientation,
    candidate_transform,
    reference_transform,
    reference_floor,
    invalid_status=None,
):
    error = (
        TransportVectorError(invalid_status, invalid_status, None, None)
        if invalid_status is not None
        else _transport_vector_error(
            candidate_value, reference_value, reference_floor=reference_floor
        )
    )
    return TransportCandidateSample(
        candidate, state[0], state[1], state[2], carrier, split, support_kind,
        support_id, candidate_value, reference_value, unit, frame, orientation,
        candidate_transform, reference_transform, error,
    )


def _result(candidate, carrier, state, split, support_kind, unit, samples, confoundings, limits):
    summary = summarize_transport_errors(
        (sample.error for sample in samples), thresholds=limits
    )
    records = tuple(confoundings)
    classification = summary.classification
    if records:
        if any(record.classification is Identifiability.CONFOUNDED for record in records):
            classification = Identifiability.CONFOUNDED
        else:
            classification = Identifiability.INSUFFICIENT_DATA
    return TransportCandidateResult(
        candidate, state[0], state[1], state[2], carrier, split, support_kind,
        unit, tuple(samples), summary.valid_count, summary.direction_valid_count,
        summary.median_abs_log10_error, summary.p95_abs_log10_error,
        summary.median_angle_deg, classification, records,
    )


def evaluate_transport_candidates(
    observations: Iterable[Observation],
    mesh: Mapping,
    *,
    density_floor: float,
    current_floor: float,
    thermal_voltage_V: float | None = None,
    thresholds: AcceptanceThresholds | None = None,
    q: float = ELEMENTARY_CHARGE_C,
) -> tuple[TransportCandidateResult, ...]:
    """Evaluate carrier-resolved QF, inverse, SG, and drift/diffusion candidates.

    ``current_inverted_qf_gradient`` rows are diagnostic candidates only.  The
    immutable discovery set is imported from the sealed input contract; every
    other checkpoint is a holdout and no caller-supplied split is accepted.
    """
    limits = _fixed_thresholds(thresholds)
    density_limit = _finite_scalar(density_floor, "density floor")
    current_limit = _finite_scalar(current_floor, "current floor")
    charge = _finite_scalar(q, "q")
    if density_limit < 0.0 or current_limit < 0.0:
        raise ValueError("floors must be non-negative")
    if charge <= 0.0:
        raise ValueError("q must be positive")
    thermal_voltage = None
    if thermal_voltage_V is not None:
        thermal_voltage = _finite_scalar(thermal_voltage_V, "thermal voltage")
        if thermal_voltage <= 0.0:
            raise ValueError("thermal voltage must be positive")

    relevant_quantities = {
        "coordinate", "ElectrostaticPotential", "eQuasiFermiPotential",
        "hQuasiFermiPotential", "eDensity", "hDensity", "eMobility",
        "hMobility", "eCurrentDensity", "hCurrentDensity",
    }
    groups = {}
    seen = set()
    for row in observations:
        if row.support_kind is not SupportKind.NODE or row.quantity not in relevant_quantities:
            continue
        if row.key in seen:
            raise ValueError("duplicate transport observation key")
        seen.add(row.key)
        state = row.solver, row.topology, float(row.bias_V)
        groups.setdefault(state, []).append(row)

    results = []
    discovery = set(DISCOVERY_KEYS)
    for state in sorted(groups, key=lambda value: (value[0], value[1], value[2])):
        rows = groups[state]
        frames = {row.coordinate_frame for row in rows}
        orientations = {row.orientation for row in rows}
        if len(frames) != 1 or len(orientations) != 1:
            raise ValueError("transport observations have incompatible coordinate frames or orientations")
        frame = next(iter(frames))
        native_orientation = next(iter(orientations))
        split = "discovery" if (state[1], state[2]) in discovery else "holdout"
        index = {(str(row.support_id), row.quantity, row.component): row for row in rows}
        node_ids = sorted(
            {key[0] for key in index if key[1] == "coordinate"}, key=_stable_key
        )
        coordinates = {}
        for node in node_ids:
            x, x_status = _status_value(index, node, "coordinate", "x", "m")
            y, y_status = _status_value(index, node, "coordinate", "y", "m")
            if x_status is not SampleStatus.VALID or y_status is not SampleStatus.VALID:
                raise ValueError("canonical SI node coordinates are required")
            coordinates[node] = (x, y)

        cells, edges = _topology(mesh, state[1])
        for nodes in cells.values():
            if len(nodes) != 3 or any(node not in coordinates for node in nodes):
                raise ValueError("invalid triangle topology")
        for endpoints in edges.values():
            if len(endpoints) != 2 or any(node not in coordinates for node in endpoints):
                raise ValueError("invalid directed edge topology")
            _edge_tangent(coordinates[endpoints[0]], coordinates[endpoints[1]])

        psi, psi_status = {}, {}
        for node in node_ids:
            psi[node], psi_status[node] = _status_value(
                index, node, "ElectrostaticPotential", "component0", "V"
            )

        for carrier in ("electron", "hole"):
            qf_quantity, density_quantity, mobility_quantity, current_quantity = _quantity_contract(carrier)
            qf, qf_status = {}, {}
            density, density_status = {}, {}
            mobility, mobility_status = {}, {}
            current, current_status = {}, {}
            for node in node_ids:
                qf[node], qf_status[node] = _status_value(
                    index, node, qf_quantity, "component0", "V"
                )
                density[node], density_status[node] = _status_value(
                    index, node, density_quantity, "component0", "m^-3"
                )
                mobility[node], mobility_status[node] = _status_value(
                    index, node, mobility_quantity, "component0", "m^2*V^-1*s^-1"
                )
                jx, jx_status = _status_value(
                    index, node, current_quantity, "component0", "A/m^2"
                )
                jy, jy_status = _status_value(
                    index, node, current_quantity, "component1", "A/m^2"
                )
                invalid = _first_invalid((jx_status, jy_status))
                current[node] = None if invalid else (jx, jy)
                current_status[node] = invalid or SampleStatus.VALID
                if density_status[node] is SampleStatus.VALID and density[node] <= density_limit:
                    density_status[node] = SampleStatus.BELOW_FLOOR
                if mobility_status[node] is SampleStatus.VALID and mobility[node] <= 0.0:
                    mobility_status[node] = SampleStatus.BELOW_FLOOR
                if current_status[node] is SampleStatus.VALID and math.hypot(*current[node]) <= current_limit:
                    current_status[node] = SampleStatus.BELOW_FLOOR

            cell_gradients, cell_gradient_status = {}, {}
            for cell_id in sorted(cells, key=_stable_key):
                nodes = cells[cell_id]
                invalid = _first_invalid(qf_status[node] for node in nodes)
                if invalid:
                    cell_gradients[cell_id], cell_gradient_status[cell_id] = None, invalid
                else:
                    cell_gradients[cell_id] = triangle_gradient(
                        tuple(coordinates[node] for node in nodes),
                        tuple(qf[node] for node in nodes),
                    )
                    cell_gradient_status[cell_id] = SampleStatus.VALID

            usable_cell_gradients = {
                cell_id: gradient for cell_id, gradient in cell_gradients.items()
                if gradient is not None
            }
            node_gradient_values = {}
            edge_gradient_values = {}
            if len(usable_cell_gradients) == len(cells):
                node_gradient_values = cell_to_node_vectors(
                    usable_cell_gradients, cells, coordinates
                )["values"]
                edge_gradient_values = cell_to_edge_vectors(
                    usable_cell_gradients, cells, edges, coordinates
                )["values"]

            specifications = (
                ("triangle_qf_gradient_current", SupportKind.CELL),
                ("node_area_weighted_qf_gradient_current", SupportKind.NODE),
                ("edge_area_weighted_qf_gradient_current", SupportKind.EDGE),
                ("signed_edge_qf_difference_current", SupportKind.EDGE),
            )
            for candidate, support in specifications:
                samples, confoundings = [], []
                support_items = (
                    tuple(sorted(cells, key=_stable_key)) if support is SupportKind.CELL
                    else tuple(node_ids) if support is SupportKind.NODE
                    else tuple(sorted(edges, key=_stable_key))
                )
                for support_id in support_items:
                    if support is SupportKind.CELL:
                        nodes = cells[support_id]
                        gradient = cell_gradients[support_id]
                        transform = "P1_triangle_qf_gradient"
                    elif support is SupportKind.NODE:
                        nodes = (support_id,)
                        gradient = node_gradient_values.get(support_id)
                        transform = "cell_qf_gradient_to_node_area_weighted"
                    else:
                        nodes = edges[support_id]
                        if candidate == "signed_edge_qf_difference_current":
                            invalid_qf = _first_invalid(qf_status[node] for node in nodes)
                            gradient = None
                            if invalid_qf is None:
                                start, end = (coordinates[node] for node in nodes)
                                signed_gradient = edge_scalar_difference(
                                    qf[nodes[0]], qf[nodes[1]], start, end
                                )
                                gradient = reconstruct_edge_vector(signed_gradient, start, end)
                            transform = "directed_edge_delta_qf_over_length"
                        else:
                            gradient = edge_gradient_values.get(support_id)
                            transform = "cell_qf_gradient_to_edge_area_weighted"
                    statuses = (
                        tuple(qf_status[node] for node in nodes)
                        + tuple(density_status[node] for node in nodes)
                        + tuple(mobility_status[node] for node in nodes)
                        + tuple(current_status[node] for node in nodes)
                    )
                    invalid = _first_invalid(statuses)
                    density_value = None if _first_invalid(density_status[node] for node in nodes) else _average(density[node] for node in nodes)
                    mobility_value = None if _first_invalid(mobility_status[node] for node in nodes) else _average(mobility[node] for node in nodes)
                    reference = None if _first_invalid(current_status[node] for node in nodes) else _average_vectors(current[node] for node in nodes)
                    if (
                        candidate == "signed_edge_qf_difference_current"
                        and reference is not None
                    ):
                        start, end = (coordinates[node] for node in nodes)
                        reference = reconstruct_edge_vector(project_vector_to_edge(reference, start, end), start, end)
                    candidate_value = None
                    if invalid is None and gradient is not None:
                        candidate_value = qf_current_density(
                            carrier, density_value, mobility_value, gradient, q=charge
                        )
                    else:
                        missing = []
                        if _first_invalid(density_status[node] for node in nodes):
                            missing.append("density")
                        if _first_invalid(mobility_status[node] for node in nodes):
                            missing.append("mobility")
                        if _first_invalid(current_status[node] for node in nodes):
                            missing.append("current")
                        if gradient is None:
                            missing.append("qf_gradient")
                        observable = observable_value = observable_unit = None
                        if (
                            missing == ["mobility"] and reference is not None
                            and density_value is not None
                        ):
                            denominator = _carrier_sign(carrier) * charge * density_value
                            observable_value = (
                                reference[0] / denominator,
                                reference[1] / denominator,
                            )
                            observable = "mu_times_grad_qf"
                            observable_unit = "m/s"
                        confoundings.append(_confounding(
                            candidate, carrier, state, support, support_id,
                            invalid or SampleStatus.MISSING_FIELD, missing,
                            observable=observable, observable_value=observable_value,
                            observable_unit=observable_unit,
                        ))
                    samples.append(_sample(
                        candidate, carrier, state, split, support, support_id,
                        candidate_value, reference, "A/m^2", frame,
                        native_orientation if support is not SupportKind.EDGE
                        else f"{native_orientation};edge=start_to_end",
                        transform,
                        "native_node_current_vector" if support is SupportKind.NODE
                        else "node_current_to_cell_equal_P1_weights" if support is SupportKind.CELL
                        else "node_current_to_edge_average_then_signed_projection"
                        if candidate == "signed_edge_qf_difference_current"
                        else "node_current_to_edge_average_preserving_vector",
                        current_limit, invalid,
                    ))
                results.append(_result(
                    candidate, carrier, state, split, support, "A/m^2",
                    samples, confoundings, limits,
                ))

            inverse_candidate = "current_inverted_qf_gradient"
            inverse_samples, inverse_confoundings = [], []
            for node in node_ids:
                reference = node_gradient_values.get(node)
                statuses = (density_status[node], mobility_status[node], current_status[node])
                invalid = _first_invalid(statuses)
                candidate_value = None
                if invalid is None and reference is not None:
                    candidate_value = current_inverted_qf_gradient(
                        carrier, density[node], mobility[node], current[node], q=charge,
                        density_floor=density_limit, mobility_floor=0.0,
                    )
                else:
                    missing = tuple(
                        name for name, status in zip(
                            ("density", "mobility", "current"), statuses
                        ) if status is not SampleStatus.VALID
                    )
                    inverse_confoundings.append(_confounding(
                        inverse_candidate, carrier, state, SupportKind.NODE, node,
                        invalid or SampleStatus.MISSING_FIELD, missing,
                    ))
                inverse_samples.append(_sample(
                    inverse_candidate, carrier, state, split, SupportKind.NODE,
                    node, candidate_value, reference, "V/m", frame,
                    native_orientation, "native_current_divided_by_q_mu_density",
                    "cell_qf_gradient_to_node_area_weighted", 0.0, invalid,
                ))
            results.append(_result(
                inverse_candidate, carrier, state, split, SupportKind.NODE, "V/m",
                inverse_samples, inverse_confoundings, limits,
            ))

            if thermal_voltage is not None:
                for candidate in (
                    "signed_edge_sg_density_current",
                    "signed_edge_drift_diffusion_current",
                ):
                    samples, confoundings = [], []
                    for edge_id in sorted(edges, key=_stable_key):
                        start_node, end_node = edges[edge_id]
                        nodes = (start_node, end_node)
                        statuses = (
                            tuple(psi_status[node] for node in nodes)
                            + tuple(density_status[node] for node in nodes)
                            + tuple(mobility_status[node] for node in nodes)
                            + tuple(current_status[node] for node in nodes)
                        )
                        invalid = _first_invalid(statuses)
                        start, end = coordinates[start_node], coordinates[end_node]
                        reference = None
                        candidate_value = None
                        if _first_invalid(current_status[node] for node in nodes) is None:
                            native = _average_vectors(current[node] for node in nodes)
                            reference = reconstruct_edge_vector(
                                project_vector_to_edge(native, start, end), start, end
                            )
                        if invalid is None:
                            length = math.dist(start, end)
                            mobility_value = _average(mobility[node] for node in nodes)
                            if candidate == "signed_edge_sg_density_current":
                                scalar = _sg_current(
                                    carrier, density[start_node], density[end_node],
                                    psi[end_node] - psi[start_node], thermal_voltage,
                                    mobility_value, length, charge,
                                )
                                transform = "density_scharfetter_gummel_then_directed_tangent"
                            else:
                                grad_psi = edge_scalar_difference(
                                    psi[start_node], psi[end_node], start, end
                                )
                                grad_density = edge_scalar_difference(
                                    density[start_node], density[end_node], start, end
                                )
                                density_value = _average(density[node] for node in nodes)
                                if carrier == "electron":
                                    scalar = charge * mobility_value * (
                                        -density_value * grad_psi
                                        + thermal_voltage * grad_density
                                    )
                                else:
                                    scalar = charge * mobility_value * (
                                        density_value * grad_psi
                                        - thermal_voltage * grad_density
                                    )
                                transform = "edge_drift_plus_diffusion_then_directed_tangent"
                            candidate_value = reconstruct_edge_vector(scalar, start, end)
                        else:
                            missing = []
                            for name, status_group in (
                                ("potential", tuple(psi_status[node] for node in nodes)),
                                ("density", tuple(density_status[node] for node in nodes)),
                                ("mobility", tuple(mobility_status[node] for node in nodes)),
                                ("current", tuple(current_status[node] for node in nodes)),
                            ):
                                if _first_invalid(status_group):
                                    missing.append(name)
                            confoundings.append(_confounding(
                                candidate, carrier, state, SupportKind.EDGE, edge_id,
                                invalid, missing,
                            ))
                            transform = (
                                "density_scharfetter_gummel_then_directed_tangent"
                                if candidate == "signed_edge_sg_density_current"
                                else "edge_drift_plus_diffusion_then_directed_tangent"
                            )
                        samples.append(_sample(
                            candidate, carrier, state, split, SupportKind.EDGE,
                            edge_id, candidate_value, reference, "A/m^2", frame,
                            f"{native_orientation};edge=start_to_end", transform,
                            "node_current_to_edge_average_then_signed_projection",
                            current_limit, invalid,
                        ))
                    results.append(_result(
                        candidate, carrier, state, split, SupportKind.EDGE, "A/m^2",
                        samples, confoundings, limits,
                    ))

    return tuple(sorted(results, key=lambda item: (
        item.solver, item.topology, item.bias_V, item.carrier,
        item.candidate, item.support_kind.value,
    )))
