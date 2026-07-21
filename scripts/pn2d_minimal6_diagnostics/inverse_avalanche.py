"""Diagnostic-only avalanche-driver inversion and source reconstruction.

All public numerical inputs use SI units unless a name explicitly states a
different convention.  Intensive generation (``m^-3*s^-1``), two-dimensional
integrals per unit out-of-plane depth (``m^-1*s^-1``), and finite-depth rates
(``s^-1``) are deliberately represented by different fields.
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
    from .inverse_fields import _topology, _usable, cell_to_node_vectors, triangle_gradient
    from .inverse_inputs import DISCOVERY_KEYS
    from .support import local_edge_sources_to_nodes, map_local_sources_to_nodes
except ImportError:
    from scripts.pn2d_minimal6_diagnostics.inverse_contracts import (  # type: ignore
        AcceptanceThresholds,
        Identifiability,
        Observation,
        SampleStatus,
        SupportKind,
    )
    from scripts.pn2d_minimal6_diagnostics.inverse_fields import (  # type: ignore
        _topology, _usable, cell_to_node_vectors, triangle_gradient,
    )
    from scripts.pn2d_minimal6_diagnostics.inverse_inputs import DISCOVERY_KEYS  # type: ignore
    from scripts.pn2d_minimal6_diagnostics.support import (  # type: ignore
        local_edge_sources_to_nodes,
        map_local_sources_to_nodes,
    )


Vector = tuple[float, float]
ELEMENTARY_CHARGE_C = 1.602176634e-19


@dataclass(frozen=True)
class GenerationError:
    status: SampleStatus
    abs_log10_error: float | None


@dataclass(frozen=True)
class GenerationMetricSummary:
    local_errors: tuple[GenerationError, ...]
    integrated_error: GenerationError
    local_valid_count: int
    local_median_abs_log10_error: float | None
    local_max_abs_log10_error: float | None
    integrated_abs_log10_error: float | None
    classification: Identifiability


@dataclass(frozen=True)
class GenerationSupportReconstruction:
    """Separate intensive, per-depth, finite-depth, and mapped source layers."""

    native_nodal_generation_m3_s: Mapping[object, float]
    candidate_cell_generation_m3_s: Mapping[object, float]
    native_cell_integrals_per_m_s: Mapping[object, float]
    candidate_cell_integrals_per_m_s: Mapping[object, float]
    native_integrated_per_m_s: float
    candidate_integrated_per_m_s: float
    native_one_cm_depth_s_inv: float
    candidate_one_cm_depth_s_inv: float
    candidate_node_mapped_per_m_s: Mapping[object, float]
    vela_edge_partial_sources_per_m_s: Mapping[object, float]
    vela_node_mapped_per_m_s: Mapping[object, float]
    depth_m: float
    native_nodal_unit: str = "m^-3*s^-1"
    candidate_cell_unit: str = "m^-3*s^-1"
    cell_integral_unit: str = "m^-1*s^-1"
    depth_integral_unit: str = "s^-1"


@dataclass(frozen=True)
class AvalancheExclusion:
    candidate: str
    solver: str
    topology: str
    bias_V: float
    carrier: str
    support_kind: SupportKind
    support_id: int | str
    status: SampleStatus
    missing_inputs: tuple[str, ...]


@dataclass(frozen=True)
class AvalancheCandidateSample:
    candidate: str
    solver: str
    topology: str
    bias_V: float
    split: str
    support_kind: SupportKind
    support_id: int | str
    electron_driver_V_m: float | None
    hole_driver_V_m: float | None
    electron_reference_driver_V_m: float | None
    hole_reference_driver_V_m: float | None
    electron_alpha_m_inv: float | None
    hole_alpha_m_inv: float | None
    candidate_generation_m3_s: float | None
    reference_generation_m3_s: float | None
    error: GenerationError


@dataclass(frozen=True)
class AvalancheCandidateResult:
    candidate: str
    solver: str
    topology: str
    bias_V: float
    split: str
    support_kind: SupportKind
    samples: tuple[AvalancheCandidateSample, ...]
    supports: GenerationSupportReconstruction | None
    summary: GenerationMetricSummary
    exclusions: tuple[AvalancheExclusion, ...]


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


def invert_van_overstraeten_alpha(
    alpha,
    *,
    prefactor,
    critical_field,
    gamma,
    branch: str | None = None,
    switch_field: float | None = None,
    numerical_floor: float = 0.0,
) -> tuple[float | None, SampleStatus]:
    """Invert ``alpha=gamma*a*exp(-gamma*b/F)`` on one declared branch.

    ``prefactor`` and ``alpha`` must use the same inverse-length unit, while
    ``critical_field``, the returned field, and ``switch_field`` use the same
    electric-field unit.  Invalid observations return a typed status; invalid
    model configuration raises ``ValueError``.
    """

    prefactor_value = _finite_scalar(prefactor, "prefactor")
    critical_value = _finite_scalar(critical_field, "critical field")
    gamma_value = _finite_scalar(gamma, "gamma")
    floor = _finite_scalar(numerical_floor, "numerical floor")
    if prefactor_value <= 0.0 or critical_value <= 0.0 or gamma_value <= 0.0:
        raise ValueError("prefactor, critical field, and gamma must be positive")
    if floor < 0.0:
        raise ValueError("numerical floor must be non-negative")
    if (branch is None) != (switch_field is None):
        raise ValueError("branch and switch field must be declared together")
    if branch is not None and branch not in {"low", "high"}:
        raise ValueError("branch must be 'low' or 'high'")
    switch = None
    if switch_field is not None:
        switch = _finite_scalar(switch_field, "switch field")
        if switch <= 0.0:
            raise ValueError("switch field must be positive")

    if alpha is None:
        return None, SampleStatus.MISSING_FIELD
    try:
        alpha_value = float(alpha)
    except (TypeError, ValueError):
        return None, SampleStatus.NONFINITE
    if not math.isfinite(alpha_value):
        return None, SampleStatus.NONFINITE
    if alpha_value <= floor:
        return None, SampleStatus.BELOW_FLOOR

    ceiling = gamma_value * prefactor_value
    if not math.isfinite(ceiling):
        raise ValueError("gamma times prefactor must be finite")
    if alpha_value >= ceiling:
        return None, SampleStatus.BRANCH_AMBIGUOUS
    ratio = alpha_value / ceiling
    if ratio == 0.0:
        return None, SampleStatus.EXPONENTIAL_UNDERFLOW
    logarithm = math.log(ratio)
    if not math.isfinite(logarithm):
        return None, SampleStatus.EXPONENTIAL_UNDERFLOW
    field = -gamma_value * critical_value / logarithm
    if not math.isfinite(field) or field <= 0.0:
        return None, SampleStatus.EXPONENTIAL_UNDERFLOW
    if branch == "low" and field >= switch:
        return None, SampleStatus.BRANCH_AMBIGUOUS
    if branch == "high" and field < switch:
        return None, SampleStatus.BRANCH_AMBIGUOUS
    return field, SampleStatus.VALID


def current_aligned_magnitude(driver, current, *, current_floor: float = 0.0) -> float:
    """Return ``|driver dot J_hat|`` without depending on current sign."""

    vector = _finite_vector(driver, "driver")
    current_vector = _finite_vector(current, "current")
    floor = _finite_scalar(current_floor, "current floor")
    if floor < 0.0:
        raise ValueError("current floor must be non-negative")
    magnitude = math.hypot(*current_vector)
    if magnitude <= floor:
        raise ValueError("current direction is undefined at or below floor")
    return abs((vector[0] * current_vector[0] + vector[1] * current_vector[1]) / magnitude)


def impact_generation(
    electron_alpha,
    electron_current,
    hole_alpha,
    hole_current,
    *,
    q: float = ELEMENTARY_CHARGE_C,
) -> float:
    """Return ``(alpha_n*|Jn| + alpha_p*|Jp|)/q`` in ``m^-3*s^-1``.

    Alpha is in ``m^-1``, conventional current density is in ``A/m^2``, and
    ``q`` is in coulombs.  Only current magnitudes enter the source.
    """

    alpha_n = _finite_scalar(electron_alpha, "electron alpha")
    alpha_p = _finite_scalar(hole_alpha, "hole alpha")
    current_n = _finite_vector(electron_current, "electron current")
    current_p = _finite_vector(hole_current, "hole current")
    charge = _finite_scalar(q, "elementary charge")
    if alpha_n < 0.0 or alpha_p < 0.0:
        raise ValueError("ionization coefficients must be non-negative")
    if charge <= 0.0:
        raise ValueError("elementary charge must be positive")
    result = (
        alpha_n * math.hypot(*current_n)
        + alpha_p * math.hypot(*current_p)
    ) / charge
    if not math.isfinite(result):
        raise ValueError("impact generation is nonfinite")
    return result


def _triangle_area(points) -> float:
    if len(points) != 3:
        raise ValueError("triangle needs exactly three coordinates")
    (x0, y0), (x1, y1), (x2, y2) = (
        _finite_vector(point, "coordinate") for point in points
    )
    twice_area = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)
    if not math.isfinite(twice_area) or twice_area <= 0.0:
        raise ValueError("generation integration requires strictly CCW triangles")
    return 0.5 * twice_area


def reconstruct_generation_supports(
    *,
    coordinates_m: Mapping,
    triangles: Mapping,
    native_nodal_generation_m3_s: Mapping,
    candidate_cell_generation_m3_s: Mapping,
    edges: Mapping | None = None,
    vela_edge_partial_sources_per_m_s: Mapping | None = None,
    depth_m: float = 0.01,
) -> GenerationSupportReconstruction:
    """Integrate and map native, cell, and Vela partial-volume source layers."""

    depth = _finite_scalar(depth_m, "out-of-plane depth")
    if depth <= 0.0:
        raise ValueError("out-of-plane depth must be positive")
    cell_ids = tuple(triangles)
    if set(cell_ids) != set(candidate_cell_generation_m3_s):
        raise ValueError("candidate cell generation and triangles must have identical support")

    native = {
        node: _finite_scalar(value, "native nodal generation")
        for node, value in native_nodal_generation_m3_s.items()
    }
    candidate = {
        cell: _finite_scalar(value, "candidate cell generation")
        for cell, value in candidate_cell_generation_m3_s.items()
    }
    if any(value < 0.0 for value in (*native.values(), *candidate.values())):
        raise ValueError("generation must be non-negative")

    normalized_triangles = {}
    native_cell_integrals = {}
    candidate_cell_integrals = {}
    for cell in cell_ids:
        nodes = tuple(triangles[cell])
        if len(nodes) != 3 or len(set(nodes)) != 3:
            raise ValueError("generation cell needs three distinct nodes")
        if any(node not in coordinates_m or node not in native for node in nodes):
            raise ValueError("generation cell lacks a coordinate or native nodal value")
        normalized_triangles[cell] = nodes
        area_m2 = _triangle_area(tuple(coordinates_m[node] for node in nodes))
        native_cell_integrals[cell] = (
            area_m2 * sum(native[node] for node in nodes) / 3.0
        )
        candidate_cell_integrals[cell] = area_m2 * candidate[cell]

    candidate_node_mapping = map_local_sources_to_nodes(
        tuple(normalized_triangles[cell] for cell in cell_ids),
        tuple(candidate_cell_integrals[cell] for cell in cell_ids),
    )

    edge_map = dict(edges or {})
    vela_edge_sources = {
        edge: _finite_scalar(value, "Vela edge partial-volume source")
        for edge, value in dict(vela_edge_partial_sources_per_m_s or {}).items()
    }
    if set(edge_map) != set(vela_edge_sources):
        raise ValueError("Vela edge sources and edges must have identical support")
    if any(value < 0.0 for value in vela_edge_sources.values()):
        raise ValueError("Vela edge partial-volume sources must be non-negative")
    vela_node_mapping = (
        local_edge_sources_to_nodes(
            tuple(edge_map[edge] for edge in edge_map),
            tuple(vela_edge_sources[edge] for edge in edge_map),
        )["values"]
        if edge_map else {}
    )

    native_total = sum(native_cell_integrals.values())
    candidate_total = sum(candidate_cell_integrals.values())
    return GenerationSupportReconstruction(
        native_nodal_generation_m3_s=native,
        candidate_cell_generation_m3_s=candidate,
        native_cell_integrals_per_m_s=native_cell_integrals,
        candidate_cell_integrals_per_m_s=candidate_cell_integrals,
        native_integrated_per_m_s=native_total,
        candidate_integrated_per_m_s=candidate_total,
        native_one_cm_depth_s_inv=native_total * depth,
        candidate_one_cm_depth_s_inv=candidate_total * depth,
        candidate_node_mapped_per_m_s=candidate_node_mapping,
        vela_edge_partial_sources_per_m_s=vela_edge_sources,
        vela_node_mapped_per_m_s=vela_node_mapping,
        depth_m=depth,
    )


def generation_error(candidate, reference, *, floor: float) -> GenerationError:
    """Return an absolute log10 source error with typed floor handling."""

    floor_value = _finite_scalar(floor, "generation floor")
    if floor_value < 0.0:
        raise ValueError("generation floor must be non-negative")
    if candidate is None or reference is None:
        return GenerationError(SampleStatus.MISSING_FIELD, None)
    try:
        candidate_value, reference_value = float(candidate), float(reference)
    except (TypeError, ValueError):
        return GenerationError(SampleStatus.NONFINITE, None)
    if not math.isfinite(candidate_value) or not math.isfinite(reference_value):
        return GenerationError(SampleStatus.NONFINITE, None)
    if candidate_value < 0.0 or reference_value < 0.0:
        return GenerationError(SampleStatus.NONFINITE, None)
    if reference_value <= floor_value or candidate_value <= floor_value:
        return GenerationError(SampleStatus.BELOW_FLOOR, None)
    error = abs(math.log10(candidate_value / reference_value))
    if not math.isfinite(error):
        return GenerationError(SampleStatus.NONFINITE, None)
    return GenerationError(SampleStatus.VALID, error)


def summarize_generation_errors(
    local_errors,
    integrated_error: GenerationError,
    *,
    thresholds: AcceptanceThresholds | None = None,
) -> GenerationMetricSummary:
    """Apply the fixed 0.3-dex local and 0.1-dex integrated gates."""

    canonical = AcceptanceThresholds()
    if thresholds is None:
        limits = canonical
    elif not isinstance(thresholds, AcceptanceThresholds):
        raise ValueError("thresholds must use the fixed inverse-audit contract")
    elif (
        thresholds.local_generation_abs_dex != canonical.local_generation_abs_dex
        or thresholds.integrated_generation_abs_dex
        != canonical.integrated_generation_abs_dex
    ):
        raise ValueError("generation thresholds are immutable at 0.3 and 0.1 dex")
    else:
        limits = thresholds

    errors = tuple(local_errors)
    values = [
        float(error.abs_log10_error)
        for error in errors
        if error.status is SampleStatus.VALID
        and error.abs_log10_error is not None
        and math.isfinite(float(error.abs_log10_error))
    ]
    integrated_value = (
        float(integrated_error.abs_log10_error)
        if integrated_error.status is SampleStatus.VALID
        and integrated_error.abs_log10_error is not None
        and math.isfinite(float(integrated_error.abs_log10_error))
        else None
    )
    median = statistics.median(values) if values else None
    maximum = max(values) if values else None
    if not values or integrated_value is None:
        classification = Identifiability.INSUFFICIENT_DATA
    elif (
        maximum <= limits.local_generation_abs_dex
        and integrated_value <= limits.integrated_generation_abs_dex
    ):
        classification = Identifiability.IDENTIFIED
    else:
        classification = Identifiability.REJECTED
    return GenerationMetricSummary(
        local_errors=errors,
        integrated_error=integrated_error,
        local_valid_count=len(values),
        local_median_abs_log10_error=median,
        local_max_abs_log10_error=maximum,
        integrated_abs_log10_error=integrated_value,
        classification=classification,
    )

def _stable_key(value: object) -> tuple[str, str]:
    return type(value).__name__, str(value)


def _first_invalid(statuses: Iterable[SampleStatus]) -> SampleStatus | None:
    for status in statuses:
        if status is not SampleStatus.VALID:
            return status
    return None


def _average(values: Iterable[float]) -> float:
    rows = tuple(values)
    if not rows:
        raise ValueError("cannot average empty support")
    return sum(rows) / len(rows)


def _average_vectors(values: Iterable[Vector]) -> Vector:
    rows = tuple(values)
    return _average(value[0] for value in rows), _average(value[1] for value in rows)


def _avalanche_parameters(parameters, carrier: str):
    if not isinstance(parameters, Mapping):
        raise ValueError("avalanche parameters must be a mapping")
    gamma = _finite_scalar(parameters.get("gamma"), "gamma")
    switch = _finite_scalar(parameters.get("switch_field_V_m"), "switch field")
    carrier_parameters = parameters.get(carrier)
    if gamma <= 0.0 or switch <= 0.0 or not isinstance(carrier_parameters, Mapping):
        raise ValueError("invalid avalanche parameter contract")
    branches = {}
    for branch in ("low", "high"):
        values = carrier_parameters.get(branch)
        try:
            if len(values) != 2:
                raise ValueError
            prefactor = _finite_scalar(values[0], f"{carrier} {branch} prefactor")
            critical = _finite_scalar(values[1], f"{carrier} {branch} critical field")
        except (TypeError, ValueError, IndexError) as error:
            raise ValueError(f"invalid {carrier} {branch} avalanche parameters") from error
        if prefactor <= 0.0 or critical <= 0.0:
            raise ValueError("avalanche branch parameters must be positive")
        branches[branch] = prefactor, critical
    return gamma, switch, branches


def _forward_van_overstraeten(driver, parameters, carrier: str):
    gamma, switch, branches = _avalanche_parameters(parameters, carrier)
    if driver is None:
        return None, SampleStatus.MISSING_FIELD
    try:
        field = float(driver)
    except (TypeError, ValueError):
        return None, SampleStatus.NONFINITE
    if not math.isfinite(field):
        return None, SampleStatus.NONFINITE
    if field <= 0.0:
        return None, SampleStatus.BELOW_FLOOR
    branch = "low" if field < switch else "high"
    prefactor, critical = branches[branch]
    alpha = gamma * prefactor * math.exp(-gamma * critical / field)
    if alpha == 0.0:
        return None, SampleStatus.EXPONENTIAL_UNDERFLOW
    if not math.isfinite(alpha):
        return None, SampleStatus.NONFINITE
    return alpha, SampleStatus.VALID


def _invert_piecewise_alpha(alpha, parameters, carrier: str, numerical_floor: float):
    gamma, switch, branches = _avalanche_parameters(parameters, carrier)
    recovered = []
    statuses = []
    for branch in ("low", "high"):
        value, status = invert_van_overstraeten_alpha(
            alpha,
            prefactor=branches[branch][0],
            critical_field=branches[branch][1],
            gamma=gamma,
            branch=branch,
            switch_field=switch,
            numerical_floor=numerical_floor,
        )
        statuses.append(status)
        if status is SampleStatus.VALID:
            recovered.append(value)
    if len(recovered) == 1:
        return recovered[0], SampleStatus.VALID
    if len(recovered) > 1:
        return None, SampleStatus.BRANCH_AMBIGUOUS
    for status in statuses:
        if status in {
            SampleStatus.MISSING_FIELD,
            SampleStatus.NONFINITE,
            SampleStatus.BELOW_FLOOR,
            SampleStatus.EXPONENTIAL_UNDERFLOW,
        }:
            return None, status
    return None, SampleStatus.BRANCH_AMBIGUOUS


def _node_scalar(index, node, quantity, unit):
    return _usable(index.get((node, quantity, "component0")), unit)


def _node_vector(index, node, quantity, unit):
    x, x_status = _usable(index.get((node, quantity, "component0")), unit)
    y, y_status = _usable(index.get((node, quantity, "component1")), unit)
    invalid = _first_invalid((x_status, y_status))
    return (None if invalid else (x, y)), (invalid or SampleStatus.VALID)


def _candidate_driver(candidate, electric_field, qf_gradient, current, density, reference_density):
    if candidate == "electric_field_magnitude":
        return (math.hypot(*electric_field), SampleStatus.VALID) if electric_field is not None else (None, SampleStatus.MISSING_FIELD)
    if candidate == "qf_gradient_magnitude":
        return (math.hypot(*qf_gradient), SampleStatus.VALID) if qf_gradient is not None else (None, SampleStatus.MISSING_FIELD)
    if candidate == "electric_field_current_aligned":
        if electric_field is None:
            return None, SampleStatus.MISSING_FIELD
        if current is None:
            return None, SampleStatus.MISSING_FIELD
        return current_aligned_magnitude(electric_field, current), SampleStatus.VALID
    if candidate == "qf_gradient_current_aligned":
        if qf_gradient is None:
            return None, SampleStatus.MISSING_FIELD
        if current is None:
            return None, SampleStatus.MISSING_FIELD
        return current_aligned_magnitude(qf_gradient, current), SampleStatus.VALID
    if candidate == "density_interpolated_qf_electric":
        if electric_field is None or qf_gradient is None or density is None:
            return None, SampleStatus.MISSING_FIELD
        weight = density / (density + reference_density)
        return (
            weight * math.hypot(*qf_gradient)
            + (1.0 - weight) * math.hypot(*electric_field),
            SampleStatus.VALID,
        )
    raise ValueError(f"unknown avalanche candidate {candidate}")


def evaluate_avalanche_candidates(
    observations,
    mesh,
    *,
    parameters,
    generation_floor,
    current_floor=0.0,
    reference_densities_m3=None,
    q=ELEMENTARY_CHARGE_C,
    vela_edge_partial_sources_per_state=None,
    thresholds=None,
    depth_m=0.01,
) -> tuple[AvalancheCandidateResult, ...]:
    """Evaluate fixed avalanche drivers on canonical nodal observations.

    Native impact generation is used only as the comparison reference.  Every
    candidate source is reconstructed from its driver, the declared
    Van Overstraeten parameters, and native current magnitude.
    """

    generation_limit = _finite_scalar(generation_floor, "generation floor")
    current_limit = _finite_scalar(current_floor, "current floor")
    charge = _finite_scalar(q, "q")
    depth = _finite_scalar(depth_m, "out-of-plane depth")
    if generation_limit < 0.0 or current_limit < 0.0:
        raise ValueError("floors must be non-negative")
    if charge <= 0.0 or depth <= 0.0:
        raise ValueError("q and out-of-plane depth must be positive")
    # Validate both carriers before touching observations.
    for carrier in ("electron", "hole"):
        _avalanche_parameters(parameters, carrier)

    reference_densities = None
    if reference_densities_m3 is not None:
        if not isinstance(reference_densities_m3, Mapping):
            raise ValueError("reference densities must be a mapping")
        reference_densities = {}
        for carrier in ("electron", "hole"):
            value = _finite_scalar(
                reference_densities_m3.get(carrier), f"{carrier} reference density"
            )
            if value <= 0.0:
                raise ValueError("reference densities must be positive")
            reference_densities[carrier] = value

    relevant = {
        "coordinate", "ElectricField", "eQuasiFermiPotential",
        "hQuasiFermiPotential", "eCurrentDensity", "hCurrentDensity",
        "eDensity", "hDensity", "eAlphaAvalanche", "hAlphaAvalanche",
        "ImpactIonization",
    }
    groups = {}
    seen = set()
    for row in observations:
        if row.support_kind is not SupportKind.NODE or row.quantity not in relevant:
            continue
        if row.key in seen:
            raise ValueError("duplicate avalanche observation key")
        seen.add(row.key)
        state = row.solver, row.topology, float(row.bias_V)
        groups.setdefault(state, []).append(row)

    candidate_names = [
        "electric_field_magnitude",
        "qf_gradient_magnitude",
        "electric_field_current_aligned",
        "qf_gradient_current_aligned",
    ]
    if reference_densities is not None:
        candidate_names.append("density_interpolated_qf_electric")

    source_by_state = dict(vela_edge_partial_sources_per_state or {})
    discovery = set(DISCOVERY_KEYS)
    results = []
    for state in sorted(groups, key=lambda value: (value[0], value[1], value[2])):
        rows = groups[state]
        frames = {row.coordinate_frame for row in rows}
        orientations = {row.orientation for row in rows}
        if len(frames) != 1 or len(orientations) != 1:
            raise ValueError("avalanche observations have incompatible frames or orientations")
        index = {
            (str(row.support_id), row.quantity, row.component): row for row in rows
        }
        node_ids = sorted(
            {key[0] for key in index if key[1] == "coordinate"}, key=_stable_key
        )
        coordinates = {}
        for node in node_ids:
            x, x_status = _usable(index.get((node, "coordinate", "x")), "m")
            y_row = index.get((node, "coordinate", "y"))
            y, y_status = _usable(y_row, "m")
            if _first_invalid((x_status, y_status)):
                raise ValueError("canonical SI node coordinates are required")
            coordinates[node] = (x, y)

        cells, edges = _topology(mesh, state[1])
        if any(node not in coordinates for nodes in cells.values() for node in nodes):
            raise ValueError("avalanche topology references an unknown node")

        node_fields = {}
        node_field_status = {}
        native_generation = {}
        native_generation_status = {}
        carrier_data = {}
        for node in node_ids:
            node_fields[node], node_field_status[node] = _node_vector(
                index, node, "ElectricField", "V/m"
            )
            native_generation[node], native_generation_status[node] = _node_scalar(
                index, node, "ImpactIonization", "m^-3*s^-1"
            )

        for carrier in ("electron", "hole"):
            prefix = "e" if carrier == "electron" else "h"
            qf = {}
            qf_status = {}
            currents = {}
            current_status = {}
            densities = {}
            density_status = {}
            native_alpha = {}
            alpha_status = {}
            reference_driver = {}
            reference_driver_status = {}
            for node in node_ids:
                qf[node], qf_status[node] = _node_scalar(
                    index, node, f"{prefix}QuasiFermiPotential", "V"
                )
                currents[node], current_status[node] = _node_vector(
                    index, node, f"{prefix}CurrentDensity", "A/m^2"
                )
                if (
                    current_status[node] is SampleStatus.VALID
                    and math.hypot(*currents[node]) <= current_limit
                ):
                    currents[node] = None
                    current_status[node] = SampleStatus.BELOW_FLOOR
                densities[node], density_status[node] = _node_scalar(
                    index, node, f"{prefix}Density", "m^-3"
                )
                if density_status[node] is SampleStatus.VALID and densities[node] < 0.0:
                    densities[node] = None
                    density_status[node] = SampleStatus.NONFINITE
                native_alpha[node], alpha_status[node] = _node_scalar(
                    index, node, f"{prefix}AlphaAvalanche", "m^-1"
                )
                if alpha_status[node] is SampleStatus.VALID:
                    (
                        reference_driver[node],
                        reference_driver_status[node],
                    ) = _invert_piecewise_alpha(
                        native_alpha[node], parameters, carrier, 0.0
                    )
                else:
                    reference_driver[node] = None
                    reference_driver_status[node] = alpha_status[node]

            cell_gradients = {}
            complete_qf = all(
                qf_status[node] is SampleStatus.VALID for node in node_ids
            )
            if complete_qf:
                for cell_id in sorted(cells, key=_stable_key):
                    nodes = cells[cell_id]
                    cell_gradients[cell_id] = triangle_gradient(
                        tuple(coordinates[node] for node in nodes),
                        tuple(qf[node] for node in nodes),
                    )
                node_gradients = cell_to_node_vectors(
                    cell_gradients, cells, coordinates
                )["values"]
                node_gradient_status = {
                    node: SampleStatus.VALID for node in node_ids
                }
            else:
                missing_status = _first_invalid(qf_status[node] for node in node_ids)
                node_gradients = {node: None for node in node_ids}
                node_gradient_status = {
                    node: missing_status or SampleStatus.MISSING_FIELD
                    for node in node_ids
                }
            carrier_data[carrier] = {
                "current": currents,
                "current_status": current_status,
                "density": densities,
                "density_status": density_status,
                "gradient": node_gradients,
                "gradient_status": node_gradient_status,
                "reference_driver": reference_driver,
                "reference_driver_status": reference_driver_status,
            }

        split = "discovery" if (state[1], state[2]) in discovery else "holdout"
        for candidate in candidate_names:
            samples = []
            exclusions = []
            candidate_node_generation = {}
            for node in node_ids:
                driver_values = {}
                alpha_values = {}
                invalid_statuses = []
                for carrier in ("electron", "hole"):
                    data = carrier_data[carrier]
                    field = (
                        node_fields[node]
                        if node_field_status[node] is SampleStatus.VALID else None
                    )
                    gradient = (
                        data["gradient"][node]
                        if data["gradient_status"][node] is SampleStatus.VALID else None
                    )
                    current = (
                        data["current"][node]
                        if data["current_status"][node] is SampleStatus.VALID else None
                    )
                    density = (
                        data["density"][node]
                        if data["density_status"][node] is SampleStatus.VALID else None
                    )
                    reference_density = (
                        reference_densities[carrier]
                        if reference_densities is not None else None
                    )
                    try:
                        driver, driver_status = _candidate_driver(
                            candidate, field, gradient, current, density,
                            reference_density,
                        )
                    except ValueError as error:
                        if "direction" not in str(error):
                            raise
                        driver, driver_status = None, SampleStatus.DIRECTION_UNDEFINED

                    dependency_statuses = []
                    if candidate.startswith("electric_field") or candidate.startswith("density_interpolated"):
                        dependency_statuses.append(node_field_status[node])
                    if candidate.startswith("qf_gradient") or candidate.startswith("density_interpolated"):
                        dependency_statuses.append(data["gradient_status"][node])
                    if candidate.endswith("current_aligned"):
                        dependency_statuses.append(data["current_status"][node])
                    if candidate.startswith("density_interpolated"):
                        dependency_statuses.append(data["density_status"][node])
                    dependency_invalid = _first_invalid(dependency_statuses)
                    if dependency_invalid is not None:
                        driver, driver_status = None, dependency_invalid

                    alpha, forward_status = _forward_van_overstraeten(
                        driver, parameters, carrier
                    )
                    status = _first_invalid((
                        driver_status,
                        forward_status,
                        data["current_status"][node],
                        data["reference_driver_status"][node],
                    ))
                    if status is not None:
                        missing = []
                        if driver_status is not SampleStatus.VALID:
                            missing.append("driving_force")
                        if data["current_status"][node] is not SampleStatus.VALID:
                            missing.append("current")
                        if data["reference_driver_status"][node] is not SampleStatus.VALID:
                            missing.append("native_alpha_branch")
                        exclusions.append(AvalancheExclusion(
                            candidate, state[0], state[1], state[2], carrier,
                            SupportKind.NODE, node, status, tuple(missing),
                        ))
                        invalid_statuses.append(status)
                        driver_values[carrier] = None
                        alpha_values[carrier] = None
                    else:
                        driver_values[carrier] = driver
                        alpha_values[carrier] = alpha

                reference_value = (
                    native_generation[node]
                    if native_generation_status[node] is SampleStatus.VALID else None
                )
                invalid = _first_invalid(invalid_statuses)
                candidate_value = None
                if invalid is None:
                    candidate_value = impact_generation(
                        alpha_values["electron"],
                        carrier_data["electron"]["current"][node],
                        alpha_values["hole"],
                        carrier_data["hole"]["current"][node],
                        q=charge,
                    )
                    candidate_node_generation[node] = candidate_value
                error = (
                    GenerationError(invalid, None)
                    if invalid is not None
                    else generation_error(
                        candidate_value, reference_value, floor=generation_limit
                    )
                )
                samples.append(AvalancheCandidateSample(
                    candidate=candidate, solver=state[0], topology=state[1],
                    bias_V=state[2], split=split, support_kind=SupportKind.NODE,
                    support_id=node,
                    electron_driver_V_m=driver_values.get("electron"),
                    hole_driver_V_m=driver_values.get("hole"),
                    electron_reference_driver_V_m=carrier_data["electron"]["reference_driver"][node],
                    hole_reference_driver_V_m=carrier_data["hole"]["reference_driver"][node],
                    electron_alpha_m_inv=alpha_values.get("electron"),
                    hole_alpha_m_inv=alpha_values.get("hole"),
                    candidate_generation_m3_s=candidate_value,
                    reference_generation_m3_s=reference_value,
                    error=error,
                ))

            supports = None
            integrated_error = GenerationError(SampleStatus.MISSING_FIELD, None)
            if (
                len(candidate_node_generation) == len(node_ids)
                and all(status is SampleStatus.VALID for status in native_generation_status.values())
            ):
                candidate_cells = {
                    cell_id: _average(
                        candidate_node_generation[node] for node in cells[cell_id]
                    )
                    for cell_id in sorted(cells, key=_stable_key)
                }
                state_edge_sources = source_by_state.get(state)
                supports = reconstruct_generation_supports(
                    coordinates_m=coordinates,
                    triangles=cells,
                    native_nodal_generation_m3_s=native_generation,
                    candidate_cell_generation_m3_s=candidate_cells,
                    edges=edges if state_edge_sources is not None else None,
                    vela_edge_partial_sources_per_m_s=state_edge_sources,
                    depth_m=depth,
                )
                integrated_error = generation_error(
                    supports.candidate_integrated_per_m_s,
                    supports.native_integrated_per_m_s,
                    floor=generation_limit,
                )
            summary = summarize_generation_errors(
                tuple(sample.error for sample in samples),
                integrated_error,
                thresholds=thresholds,
            )
            results.append(AvalancheCandidateResult(
                candidate=candidate, solver=state[0], topology=state[1],
                bias_V=state[2], split=split, support_kind=SupportKind.NODE,
                samples=tuple(samples), supports=supports, summary=summary,
                exclusions=tuple(exclusions),
            ))

    return tuple(sorted(results, key=lambda item: (
        item.solver, item.topology, item.bias_V, item.candidate,
        item.support_kind.value,
    )))
