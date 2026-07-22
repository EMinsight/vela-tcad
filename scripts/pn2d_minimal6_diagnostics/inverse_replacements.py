"""Causal replacements and formula identifiability for the Minimal6 audit.

This module is diagnostic-only.  It deliberately keeps whole-state replay out
of the formula dependency graph and accepts no locally fitted candidate scale.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Callable, Iterable, Mapping
import math
import statistics
from typing import Any

try:
    from .counterfactual import DependencyCounterfactualEngine, interaction_dex
    from .inverse_contracts import (
        AcceptanceThresholds,
        Identifiability,
        SampleStatus,
        SupportKind,
    )
    from .inverse_inputs import DISCOVERY_KEYS
except ImportError:
    from scripts.pn2d_minimal6_diagnostics.counterfactual import (  # type: ignore
        DependencyCounterfactualEngine,
        interaction_dex,
    )
    from scripts.pn2d_minimal6_diagnostics.inverse_contracts import (  # type: ignore
        AcceptanceThresholds,
        Identifiability,
        SampleStatus,
        SupportKind,
    )
    from scripts.pn2d_minimal6_diagnostics.inverse_inputs import DISCOVERY_KEYS  # type: ignore


INVERSE_DEPENDENCIES: dict[str, tuple[str, ...]] = {
    "gradient_recovery": (),
    "mobility": ("gradient_recovery",),
    "current_semantics": ("mobility",),
    "impact_driving_field": ("current_semantics",),
    "alpha_law": ("impact_driving_field",),
    "geometric_integration": ("alpha_law",),
    "source_to_node_mapping": ("geometric_integration",),
}

CLOSURE_TOLERANCE_DEX = 1.0e-10
NONUNIQUE_TOLERANCE_DEX = 1.0e-10
_STATE_LAYERS = ("potential", "carrier_state", "quasi_fermi_state")
_COMPATIBILITY_FIELDS = (
    "support_kind", "support_id", "unit_si", "carrier", "topology", "bias_V"
)
_CONTEXT_FIELDS = ("topology", "bias_V", "carrier", "support_kind", "support_id")
_FORBIDDEN_FIT_DIMENSIONS = {
    "bias", "bias_V", "node", "edge", "cell", "support", "support_id",
    "carrier", "topology",
}


def _field(record: object, name: str, default: Any = None) -> Any:
    if isinstance(record, Mapping):
        return record.get(name, default)
    return getattr(record, name, default)


def _status(value: object) -> SampleStatus:
    if isinstance(value, SampleStatus):
        return value
    try:
        return SampleStatus(str(value))
    except ValueError as error:
        raise ValueError(f"unknown sample status {value!r}") from error


def _support(value: object) -> str:
    if isinstance(value, SupportKind):
        return value.value
    try:
        return SupportKind(str(value)).value
    except ValueError as error:
        raise ValueError(f"unknown support_kind {value!r}") from error


def _finite_number(value: object, label: str, *, positive: bool) -> float:
    if isinstance(value, bool):
        raise TypeError(f"{label} must be numeric, not boolean")
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{label} must be finite") from error
    if not math.isfinite(result):
        raise ValueError(f"{label} must be finite")
    if positive and result <= 0.0:
        raise ValueError(f"{label} must be finite and nonzero positive")
    return result


def _normalize_operand(factor: str, record: object, *, positive: bool = True) -> dict:
    declared_factor = _field(record, "factor")
    if declared_factor != factor:
        raise ValueError(f"operand factor mismatch for {factor}")
    missing = []
    for name in _COMPATIBILITY_FIELDS:
        if isinstance(record, Mapping):
            present = name in record
        else:
            present = hasattr(record, name)
        if not present:
            missing.append(name)
    if missing:
        raise ValueError(f"{factor} lacks compatibility metadata: {', '.join(missing)}")
    status = _status(_field(record, "status", SampleStatus.VALID))
    normalized = {
        "factor": factor,
        "status": status,
        "support_kind": _support(_field(record, "support_kind")),
        "support_id": _field(record, "support_id"),
        "unit_si": str(_field(record, "unit_si")),
        "carrier": _field(record, "carrier"),
        "topology": str(_field(record, "topology")),
        "bias_V": _finite_number(_field(record, "bias_V"), f"{factor} bias_V", positive=False),
        "value": _field(record, "value"),
    }
    if status is SampleStatus.VALID:
        normalized["value"] = _finite_number(
            normalized["value"], f"{factor} value", positive=positive
        )
    return normalized


def _validate_operand_sets(baseline_values, replacement_values, *, positive=True):
    expected = set(INVERSE_DEPENDENCIES)
    for label, values in (("baseline", baseline_values),
                          ("replacement", replacement_values)):
        if not isinstance(values, Mapping) or set(values) != expected:
            raise ValueError(f"{label} values must cover exactly the inverse dependency graph")
    baseline, replacement = {}, {}
    for factor in INVERSE_DEPENDENCIES:
        baseline[factor] = _normalize_operand(
            factor, baseline_values[factor], positive=positive
        )
        replacement[factor] = _normalize_operand(
            factor, replacement_values[factor], positive=positive
        )
        for field in _COMPATIBILITY_FIELDS:
            if baseline[factor][field] != replacement[factor][field]:
                raise ValueError(
                    f"{factor} {field} mismatch: "
                    f"{baseline[factor][field]!r} != {replacement[factor][field]!r}"
                )
    return baseline, replacement


def _validate_global_context(records: Mapping[str, Mapping[str, object]], label: str) -> None:
    for field in _CONTEXT_FIELDS:
        values = [(name, row[field]) for name, row in records.items()
                  if row[field] is not None]
        if not values:
            continue
        expected = values[0][1]
        mismatched = [(name, value) for name, value in values[1:]
                      if value != expected]
        if mismatched:
            details = ", ".join(f"{name}={value!r}" for name, value in mismatched)
            raise ValueError(
                f"global {label} {field} mismatch: expected {expected!r}; {details}"
            )


def _chain_operators() -> dict[str, Callable]:
    operators: dict[str, Callable] = {}
    for factor, parents in INVERSE_DEPENDENCIES.items():
        if not parents:
            operators[factor] = lambda _inputs, raw: raw
        else:
            parent = parents[0]
            operators[factor] = (
                lambda inputs, raw, parent=parent: inputs[parent] * raw
            )
    return operators


def _typed_unavailable(factor: str, status: SampleStatus) -> dict:
    return {
        "status": status.value,
        "unavailable_factor": factor,
        "dependency_order": list(INVERSE_DEPENDENCIES),
        "baseline": None,
        "one_factor": [],
        "forward": [],
        "reverse": [],
        "full_replacement": None,
        "adjacent_interactions": [],
        "closure": None,
    }


def run_replacement_matrix(
    baseline_values: Mapping[str, object],
    replacement_values: Mapping[str, object],
    *,
    direct_target: float | None = None,
) -> dict:
    """Run isolated, forward-staged, and reverse-restoration replacements.

    Each factor value is a positive multiplicative operator carrying explicit
    support, unit, carrier, and state metadata.  The arithmetic is delegated to
    :class:`DependencyCounterfactualEngine`; this wrapper owns the physical DAG
    and the immutable closure gate.
    """
    baseline, replacement = _validate_operand_sets(
        baseline_values, replacement_values
    )
    _validate_global_context(baseline, "replacement matrix")
    _validate_global_context(replacement, "replacement matrix")
    if direct_target is None:
        raise ValueError("direct_target is required for independent replacement closure")
    target = _finite_number(direct_target, "direct_target", positive=True)
    for factor in INVERSE_DEPENDENCIES:
        for operand_set in (baseline, replacement):
            status = operand_set[factor]["status"]
            if status is not SampleStatus.VALID:
                return _typed_unavailable(factor, status)

    engine = DependencyCounterfactualEngine(
        dependencies=INVERSE_DEPENDENCIES,
        baseline_values={name: row["value"] for name, row in baseline.items()},
        replacement_values={name: row["value"] for name, row in replacement.items()},
        operators=_chain_operators(),
        output_factor="source_to_node_mapping",
    )
    order = tuple(INVERSE_DEPENDENCIES)
    baseline_stages = engine._full(set())
    baseline_output = baseline_stages["source_to_node_mapping"]

    one_factor = []
    for factor in order:
        stages = engine._full({factor})
        value = stages["source_to_node_mapping"]
        one_factor.append({
            "factor": factor,
            "value": value,
            "delta_dex": math.log10(value / baseline_output),
            "stage_values": stages,
            "changed_stages": list(engine.downstream[factor]),
        })

    replaced: set[str] = set()
    current = baseline_output
    forward = []
    for factor in order:
        replaced.add(factor)
        stages = engine._full(replaced)
        value = stages["source_to_node_mapping"]
        forward.append({
            "factor": factor,
            "replaced_factors": [name for name in order if name in replaced],
            "value": value,
            "incremental_dex": math.log10(value / current),
            "stage_values": stages,
            "recomputed_stages": list(engine.downstream[factor]),
        })
        current = value
    full_replacement = current

    restored = set(order)
    current = full_replacement
    reverse = []
    for factor in reversed(order):
        restored.remove(factor)
        stages = engine._full(restored)
        value = stages["source_to_node_mapping"]
        reverse.append({
            "factor": factor,
            "remaining_replacements": [name for name in order if name in restored],
            "value": value,
            "incremental_dex": math.log10(value / current),
            "stage_values": stages,
            "recomputed_stages": list(engine.downstream[factor]),
        })
        current = value

    forward_gap = math.log10(full_replacement / baseline_output)
    reverse_gap = math.log10(baseline_output / full_replacement)
    closure = {
        "tolerance_dex": CLOSURE_TOLERANCE_DEX,
        "forward_abs_dex": abs(
            forward_gap - sum(row["incremental_dex"] for row in forward)
        ),
        "reverse_abs_dex": abs(
            reverse_gap - sum(row["incremental_dex"] for row in reverse)
        ),
        "direct_abs_dex": abs(math.log10(full_replacement / target)),
    }
    failed = {name: value for name, value in closure.items()
              if name.endswith("_abs_dex") and value > CLOSURE_TOLERANCE_DEX}
    if failed:
        details = ", ".join(f"{name}={value:.17g}" for name, value in failed.items())
        raise ValueError(
            f"replacement closure exceeds immutable {CLOSURE_TOLERANCE_DEX:.1e} dex: {details}"
        )

    adjacent_interactions = []
    for first, second in zip(order, order[1:]):
        adjacent_interactions.append({
            "first_factor": first,
            "second_factor": second,
            "interaction_dex": interaction_dex(
                baseline=baseline_output,
                a_only=engine.evaluate_replacements({first}),
                b_only=engine.evaluate_replacements({second}),
                both=engine.evaluate_replacements({first, second}),
            ),
        })
    return {
        "status": SampleStatus.VALID.value,
        "dependency_order": list(order),
        "baseline": baseline_output,
        "baseline_stage_values": baseline_stages,
        "one_factor": one_factor,
        "forward": forward,
        "reverse": reverse,
        "full_replacement": full_replacement,
        "direct_target": target,
        "adjacent_interactions": adjacent_interactions,
        "closure": closure,
    }


def _candidate_split(record: object) -> str:
    explicit = _field(record, "split")
    if explicit in ("discovery", "holdout"):
        return str(explicit)
    key = (str(_field(record, "topology")), float(_field(record, "bias_V")))
    return "discovery" if key in DISCOVERY_KEYS else "holdout"


def _signed_vector_prediction(candidate, reference, status: SampleStatus) -> float | None:
    if status is not SampleStatus.VALID or candidate is None or reference is None:
        return None
    try:
        candidate_magnitude = math.hypot(float(candidate[0]), float(candidate[1]))
        reference_magnitude = math.hypot(float(reference[0]), float(reference[1]))
    except (TypeError, ValueError, IndexError):
        return None
    if (not math.isfinite(candidate_magnitude)
            or not math.isfinite(reference_magnitude)
            or candidate_magnitude <= 0.0 or reference_magnitude <= 0.0):
        return None
    return math.log10(candidate_magnitude / reference_magnitude)


def _signed_scalar_prediction(candidate, reference, status: SampleStatus) -> float | None:
    if status is not SampleStatus.VALID or candidate is None or reference is None:
        return None
    try:
        candidate_value, reference_value = float(candidate), float(reference)
    except (TypeError, ValueError):
        return None
    if (not math.isfinite(candidate_value) or not math.isfinite(reference_value)
            or candidate_value <= 0.0 or reference_value <= 0.0):
        return None
    return math.log10(candidate_value / reference_value)


def _formula_evidence(sample: object, *, factor: str, metric: str,
                      status: SampleStatus, error: float | None,
                      prediction_dex: float | None,
                      carrier: str | None = None,
                      support_kind: object | None = None,
                      support_id: object | None = None,
                      candidate: str | None = None,
                      split: str | None = None) -> dict:
    return {
        "record_kind": "formula_candidate",
        "candidate": candidate or str(_field(sample, "candidate")),
        "solver": str(_field(sample, "solver")),
        "factor": factor,
        "split": split or _candidate_split(sample),
        "topology": str(_field(sample, "topology")),
        "bias_V": float(_field(sample, "bias_V")),
        "support_kind": support_kind if support_kind is not None else _field(sample, "support_kind"),
        "support_id": support_id if support_id is not None else _field(sample, "support_id"),
        "carrier": carrier if carrier is not None else _field(sample, "carrier"),
        "metric": metric,
        "error": error,
        "prediction_dex": prediction_dex,
        "status": status,
    }


def _flatten_field_sample(sample: object) -> list[dict]:
    error = _field(sample, "error")
    magnitude_status = _status(_field(error, "magnitude_status"))
    direction_status = _status(_field(error, "direction_status"))
    prediction = _signed_vector_prediction(
        _field(sample, "candidate_value"), _field(sample, "reference_value"),
        magnitude_status,
    )
    return [
        _formula_evidence(
            sample, factor="gradient_recovery", metric="field_magnitude_relative",
            status=magnitude_status, error=_field(error, "relative_magnitude_error"),
            prediction_dex=prediction,
        ),
        _formula_evidence(
            sample, factor="gradient_recovery", metric="field_direction_deg",
            status=direction_status, error=_field(error, "angle_deg"),
            prediction_dex=prediction,
        ),
    ]


def _flatten_transport_sample(sample: object) -> list[dict]:
    error = _field(sample, "error")
    magnitude_status = _status(_field(error, "magnitude_status"))
    direction_status = _status(_field(error, "direction_status"))
    prediction = _signed_vector_prediction(
        _field(sample, "candidate_value"), _field(sample, "reference_value"),
        magnitude_status,
    )
    carrier = str(_field(sample, "carrier"))
    return [
        _formula_evidence(
            sample, factor="current_semantics", metric="transport_abs_dex",
            status=magnitude_status, error=_field(error, "abs_log10_error"),
            prediction_dex=prediction, carrier=carrier,
        ),
        _formula_evidence(
            sample, factor="current_semantics", metric="transport_direction_deg",
            status=direction_status, error=_field(error, "angle_deg"),
            prediction_dex=prediction, carrier=carrier,
        ),
    ]


def _flatten_transport_confounding(record: object, result: object) -> dict:
    missing_inputs = tuple(str(value) for value in (_field(record, "missing_inputs", ()) or ()))
    missing = ("mobility",) if "mobility" in missing_inputs else ()
    row = _formula_evidence(
        record, factor="current_semantics", metric="transport_abs_dex",
        status=_status(_field(record, "status")), error=None,
        prediction_dex=None, carrier=str(_field(record, "carrier")),
        candidate=str(_field(record, "candidate")), split=_candidate_split(result),
    )
    row["solver"] = str(_field(result, "solver"))
    row["missing_independent_factors"] = missing
    return row


def _flatten_avalanche_sample(sample: object) -> list[dict]:
    error = _field(sample, "error")
    status = _status(_field(error, "status"))
    prediction = _signed_scalar_prediction(
        _field(sample, "candidate_generation_m3_s"),
        _field(sample, "reference_generation_m3_s"), status,
    )
    return [_formula_evidence(
        sample, factor="impact_driving_field", metric="local_generation_abs_dex",
        status=status, error=_field(error, "abs_log10_error"),
        prediction_dex=prediction,
    )]


def _flatten_avalanche_integrated(result: object) -> dict:
    summary = _field(result, "summary")
    error = _field(summary, "integrated_error")
    status = _status(_field(error, "status"))
    supports = _field(result, "supports")
    prediction = None
    if supports is not None:
        prediction = _signed_scalar_prediction(
            _field(supports, "candidate_integrated_per_m_s"),
            _field(supports, "native_integrated_per_m_s"), status,
        )
    return _formula_evidence(
        result, factor="impact_driving_field",
        metric="integrated_generation_abs_dex", status=status,
        error=_field(error, "abs_log10_error"), prediction_dex=prediction,
        support_kind=SupportKind.INTEGRATED, support_id="integrated",
        split=_candidate_split(result),
    )


def _flatten_error_record(record: object) -> list[dict] | None:
    name = type(record).__name__
    if name == "VectorErrorResult":
        return [
            {"status": _field(record, "magnitude_status"),
             "error": _field(record, "relative_magnitude_error")},
            {"status": _field(record, "direction_status"),
             "error": _field(record, "angle_deg")},
        ]
    if name == "TransportVectorError":
        return [
            {"status": _field(record, "magnitude_status"),
             "error": _field(record, "abs_log10_error")},
            {"status": _field(record, "direction_status"),
             "error": _field(record, "angle_deg")},
        ]
    if name == "GenerationError":
        return [{"status": _field(record, "status"),
                 "error": _field(record, "abs_log10_error")}]
    return None


def _flatten_evidence(records: Iterable[object]) -> list[object]:
    flattened: list[object] = []
    for record in records:
        if isinstance(record, Mapping) or isinstance(record, (int, float)):
            flattened.append(record)
            continue
        name = type(record).__name__
        errors = _flatten_error_record(record)
        if errors is not None:
            flattened.extend(errors)
        elif name == "FieldCandidateSample":
            flattened.extend(_flatten_field_sample(record))
        elif name == "TransportCandidateSample":
            flattened.extend(_flatten_transport_sample(record))
        elif name == "AvalancheCandidateSample":
            flattened.extend(_flatten_avalanche_sample(record))
        elif name == "FieldCandidateResult":
            for sample in _field(record, "samples", ()):
                flattened.extend(_flatten_field_sample(sample))
        elif name == "TransportCandidateResult":
            for sample in _field(record, "samples", ()):
                flattened.extend(_flatten_transport_sample(sample))
            for confounding in _field(record, "confoundings", ()):
                flattened.append(_flatten_transport_confounding(confounding, record))
        elif name == "AvalancheCandidateResult":
            for sample in _field(record, "samples", ()):
                flattened.extend(_flatten_avalanche_sample(sample))
            flattened.append(_flatten_avalanche_integrated(record))
        else:
            flattened.append(record)
    return flattened
def replacement_from_evidence(records: Iterable[object]) -> dict:
    """Build a fail-closed replacement result from typed candidate evidence.

    Only factors carrying finite signed prediction evidence are observed.
    Absent independent factors remain typed unavailable; they are never filled
    with identity values and no direct closure target is synthesized.
    """
    flattened = [
        row for row in _flatten_evidence(records)
        if isinstance(row, Mapping)
        and _field(row, "record_kind", "formula_candidate") == "formula_candidate"
    ]
    candidates = sorted({
        str(_field(row, "candidate")) for row in flattened
        if _field(row, "candidate") is not None
    })
    observed: dict[str, float] = {}
    for factor in INVERSE_DEPENDENCIES:
        predictions = []
        for row in flattened:
            if _field(row, "factor") != factor:
                continue
            if _status(_field(row, "status")) is not SampleStatus.VALID:
                continue
            value = _field(row, "prediction_dex")
            if isinstance(value, (int, float)) and not isinstance(value, bool):
                value = float(value)
                if math.isfinite(value):
                    predictions.append(value)
        if predictions:
            observed[factor] = statistics.median(predictions)
    unavailable = next(
        (factor for factor in INVERSE_DEPENDENCIES if factor not in observed),
        None,
    )
    if unavailable is None:
        raise ValueError(
            "typed evidence lacks an independently observed direct replacement target"
        )
    result = _typed_unavailable(unavailable, SampleStatus.MISSING_FIELD)
    result.update({
        "evidence_source": "typed_candidate_evidence",
        "evidence_candidates": candidates,
        "observed_prediction_dex": observed,
    })
    return result

def _percentile(sorted_values: list[float], percentile: float) -> float:
    if len(sorted_values) == 1:
        return sorted_values[0]
    position = (len(sorted_values) - 1) * percentile
    lower = math.floor(position)
    upper = math.ceil(position)
    fraction = position - lower
    return sorted_values[lower] + fraction * (
        sorted_values[upper] - sorted_values[lower]
    )


def metric_summary(records: Iterable[object], *, value_field: str = "error") -> dict:
    """Summarize finite, explicitly valid errors without coercing typed masks."""
    values: list[float] = []
    statuses: Counter[str] = Counter()
    total = 0
    for record in _flatten_evidence(records):
        total += 1
        if isinstance(record, (int, float)) and not isinstance(record, bool):
            status = SampleStatus.VALID
            raw = record
        else:
            if isinstance(record, Mapping) and "status" not in record:
                raise ValueError("metric record lacks typed status")
            if not isinstance(record, Mapping) and not hasattr(record, "status"):
                raise ValueError("metric record lacks typed status")
            status = _status(_field(record, "status"))
            raw = _field(record, value_field)
        statuses[status.value] += 1
        if status is not SampleStatus.VALID:
            continue
        value = _finite_number(raw, "valid metric error", positive=False)
        values.append(abs(value))
    values.sort()
    if values:
        middle = len(values) // 2
        median = (values[middle] if len(values) % 2 else
                  0.5 * (values[middle - 1] + values[middle]))
        p95 = _percentile(values, 0.95)
    else:
        median = p95 = None
    return {
        "total_count": total,
        "valid_count": len(values),
        "median_abs_error": median,
        "p95_abs_error": p95,
        "status_counts": dict(sorted(statuses.items())),
    }


def _candidate_rows(candidate: str, records: Iterable[object]) -> tuple[list[object], list[object]]:
    materialized = list(records)
    controls = [row for row in materialized
                if _field(row, "record_kind", "formula_candidate") == "localization_control"]
    rows = [row for row in materialized
            if _field(row, "record_kind", "formula_candidate") == "formula_candidate"
            and str(_field(row, "candidate")) == candidate]
    if not rows and any(str(_field(row, "candidate")) == candidate for row in controls):
        raise ValueError("whole-state localization control is excluded from formula classification")
    if not rows:
        raise ValueError(f"candidate {candidate!r} has no formula evidence")
    return rows, materialized


def _metric_limit(metric: str, thresholds: AcceptanceThresholds) -> tuple[float, float]:
    if metric in {"field_magnitude_relative", "field_relative"}:
        return thresholds.field_median_relative, math.inf
    if metric in {"field_direction_deg", "field_angle_deg"}:
        return thresholds.field_median_angle_deg, math.inf
    if metric in {"gradient_abs_dex", "current_abs_dex", "transport_abs_dex"}:
        return thresholds.gradient_median_abs_dex, thresholds.gradient_p95_abs_dex
    if metric in {"gradient_direction_deg", "current_direction_deg",
                  "transport_direction_deg"}:
        return thresholds.gradient_median_angle_deg, math.inf
    if metric == "integrated_generation_abs_dex":
        return thresholds.integrated_generation_abs_dex, thresholds.integrated_generation_abs_dex
    if metric == "local_generation_abs_dex":
        return thresholds.local_generation_abs_dex, thresholds.local_generation_abs_dex
    if metric == "replacement_closure_abs_dex":
        return CLOSURE_TOLERANCE_DEX, CLOSURE_TOLERANCE_DEX
    raise ValueError(f"unknown candidate metric {metric!r}")


def _fit_leakage(rows: Iterable[object]) -> tuple[str, ...]:
    leaking = set()
    for row in rows:
        scope = _field(row, "fit_scope")
        if scope not in (None, "global", "discovery_global"):
            leaking.add(str(scope))
        dimensions = _field(row, "fit_dimensions", ()) or ()
        leaking.update(str(value) for value in dimensions
                       if str(value) in _FORBIDDEN_FIT_DIMENSIONS)
    return tuple(sorted(leaking))


def _valid_rows(rows: Iterable[object]) -> list[object]:
    return [row for row in rows if _status(_field(row, "status")) is SampleStatus.VALID]


def _base_classification(candidate: str, records: Iterable[object], *,
                         thresholds: AcceptanceThresholds,
                         minimum_valid_samples: int) -> Identifiability:
    rows, _ = _candidate_rows(candidate, records)
    if minimum_valid_samples < 1:
        raise ValueError("minimum_valid_samples must be positive")
    missing_factors = {
        str(factor) for row in rows
        for factor in (_field(row, "missing_independent_factors", ()) or ())
    }
    if missing_factors:
        return Identifiability.CONFOUNDED

    valid = _valid_rows(rows)
    discovery = [row for row in valid if _field(row, "split") == "discovery"]
    holdout = [row for row in valid if _field(row, "split") == "holdout"]
    if (len(discovery) < minimum_valid_samples
            or len(holdout) < minimum_valid_samples):
        return Identifiability.INSUFFICIENT_DATA
    if _fit_leakage(rows):
        return Identifiability.REJECTED

    metrics = sorted({str(_field(row, "metric")) for row in valid})
    if not metrics:
        return Identifiability.INSUFFICIENT_DATA
    for metric in metrics:
        combined_rows = [row for row in rows if str(_field(row, "metric")) == metric]
        holdout_rows = [row for row in combined_rows if _field(row, "split") == "holdout"]
        combined_summary = metric_summary(combined_rows)
        holdout_summary = metric_summary(holdout_rows)
        if (combined_summary["valid_count"] < minimum_valid_samples
                or holdout_summary["valid_count"] < minimum_valid_samples):
            return Identifiability.INSUFFICIENT_DATA
        median_limit, p95_limit = _metric_limit(metric, thresholds)
        for summary in (combined_summary, holdout_summary):
            if (summary["median_abs_error"] > median_limit
                    or summary["p95_abs_error"] > p95_limit):
                return Identifiability.REJECTED
    return Identifiability.IDENTIFIED


def _sample_identity(row: object) -> tuple:
    support = _field(row, "support_kind")
    if isinstance(support, SupportKind):
        support = support.value
    return (
        str(_field(row, "factor")), str(_field(row, "metric")),
        str(_field(row, "solver")), str(_field(row, "split")),
        str(_field(row, "topology")), float(_field(row, "bias_V")),
        str(support), str(_field(row, "support_id")),
        str(_field(row, "carrier")), str(_field(row, "component", "")),
    )


def _indistinguishable(first: str, second: str, records: Iterable[object]) -> bool:
    predictions = {}
    for candidate in (first, second):
        rows = [row for row in records
                if str(_field(row, "candidate")) == candidate
                and _field(row, "record_kind", "formula_candidate") == "formula_candidate"
                and _status(_field(row, "status")) is SampleStatus.VALID]
        by_key = {}
        for row in rows:
            raw = _field(row, "prediction_dex")
            if raw is None:
                return False
            value = _finite_number(raw, "valid candidate prediction_dex", positive=False)
            key = _sample_identity(row)
            if key in by_key:
                raise ValueError(f"duplicate candidate sample identity for {candidate}: {key}")
            by_key[key] = value
        predictions[candidate] = by_key
    first_values, second_values = predictions[first], predictions[second]
    if not first_values or set(first_values) != set(second_values):
        return False
    return all(abs(first_values[key] - second_values[key]) <= NONUNIQUE_TOLERANCE_DEX
               for key in first_values)


def classify_candidate(
    candidate: str,
    records: Iterable[object],
    *,
    thresholds: AcceptanceThresholds | None = None,
    minimum_valid_samples: int = 1,
) -> Identifiability:
    """Classify one formula candidate using combined and holdout evidence."""
    materialized = _flatten_evidence(records)
    thresholds = thresholds or AcceptanceThresholds()
    base = _base_classification(
        candidate, materialized, thresholds=thresholds,
        minimum_valid_samples=minimum_valid_samples,
    )
    if base is not Identifiability.IDENTIFIED:
        return base
    peers = sorted({str(_field(row, "candidate")) for row in materialized
                    if _field(row, "record_kind", "formula_candidate") == "formula_candidate"
                    and str(_field(row, "candidate")) != candidate})
    for peer in peers:
        peer_base = _base_classification(
            peer, materialized, thresholds=thresholds,
            minimum_valid_samples=minimum_valid_samples,
        )
        if (peer_base is Identifiability.IDENTIFIED
                and _indistinguishable(candidate, peer, materialized)):
            return Identifiability.CONSISTENT_NONUNIQUE
    return base


def _summaries_by_metric(rows: Iterable[object]) -> dict[str, dict]:
    materialized = list(rows)
    metrics = sorted({str(_field(row, "metric")) for row in materialized
                      if _field(row, "metric") is not None})
    return {metric: metric_summary(
        [row for row in materialized if str(_field(row, "metric")) == metric]
    ) for metric in metrics}


def _classification_reason(classification: Identifiability, rows: list[object]) -> str:
    if classification is Identifiability.IDENTIFIED:
        return "all combined and holdout gates pass"
    if classification is Identifiability.CONSISTENT_NONUNIQUE:
        return "passing prediction is indistinguishable within 1e-10 dex"
    if classification is Identifiability.CONFOUNDED:
        missing = sorted({str(factor) for row in rows
                          for factor in (_field(row, "missing_independent_factors", ()) or ())})
        return f"independent factors are absent: {', '.join(missing)}"
    if classification is Identifiability.INSUFFICIENT_DATA:
        return "discovery or holdout lacks valid typed support"
    leakage = _fit_leakage(rows)
    if leakage:
        return f"local fit leakage is forbidden: {', '.join(leakage)}"
    return "combined or holdout acceptance gate failed"


def rank_candidates(
    records: Iterable[object],
    *,
    thresholds: AcceptanceThresholds | None = None,
    minimum_valid_samples: int = 1,
) -> list[dict]:
    """Rank formula candidates by discovery evidence only.

    Holdout evidence affects the classification gate, never the ranking score.
    Localization controls are omitted before candidate names are collected.
    """
    materialized = [row for row in _flatten_evidence(records)
                    if _field(row, "record_kind", "formula_candidate") == "formula_candidate"]
    if not materialized:
        return []
    thresholds = thresholds or AcceptanceThresholds()
    candidates = sorted({str(_field(row, "candidate")) for row in materialized})
    ranked = []
    for candidate in candidates:
        rows = [row for row in materialized if str(_field(row, "candidate")) == candidate]
        discovery = [row for row in rows if _field(row, "split") == "discovery"]
        holdout = [row for row in rows if _field(row, "split") == "holdout"]
        classification = classify_candidate(
            candidate, materialized, thresholds=thresholds,
            minimum_valid_samples=minimum_valid_samples,
        )
        discovery_valid = _valid_rows(discovery)
        score = (sum(abs(_finite_number(_field(row, "error"),
                                        "valid discovery error", positive=False))
                     for row in discovery_valid) / len(discovery_valid)
                 if discovery_valid else math.inf)
        ranked.append({
            "candidate": candidate,
            "classification": classification,
            "reason": _classification_reason(classification, rows),
            "discovery_score": score,
            "discovery_metrics": _summaries_by_metric(discovery),
            "combined_metrics": _summaries_by_metric(rows),
            "holdout_metrics": _summaries_by_metric(holdout),
        })
    return sorted(
        ranked, key=lambda row: (row["discovery_score"], row["candidate"])
    )


def run_state_localization_control(
    baseline_state: Mapping[str, object],
    replacement_state: Mapping[str, object],
    *,
    evaluate: Callable[[Mapping[str, float]], float],
) -> dict:
    """Replay all observed state layers together for localization only."""
    if set(baseline_state) != set(_STATE_LAYERS) or set(replacement_state) != set(_STATE_LAYERS):
        raise ValueError("whole-state control requires potential, carrier, and quasi-Fermi state")
    baseline, replacement = {}, {}
    unavailable: tuple[str, SampleStatus] | None = None
    for layer in _STATE_LAYERS:
        baseline[layer] = _normalize_operand(layer, baseline_state[layer], positive=False)
        replacement[layer] = _normalize_operand(layer, replacement_state[layer], positive=False)
        for field in _COMPATIBILITY_FIELDS:
            if baseline[layer][field] != replacement[layer][field]:
                raise ValueError(f"{layer} {field} mismatch in localization control")
        for values in (baseline, replacement):
            status = values[layer]["status"]
            if status is not SampleStatus.VALID and unavailable is None:
                unavailable = layer, status
    _validate_global_context(baseline, "localization control")
    _validate_global_context(replacement, "localization control")
    if unavailable is not None:
        layer, status = unavailable
        return {
            "record_kind": "localization_control",
            "candidate": "whole_state_sentaurus_replay",
            "classification": "localization_control",
            "status": status.value,
            "unavailable_layer": layer,
            "eligible_for_candidate_ranking": False,
            "eligible_for_formula_classification": False,
            "layers_replaced": list(_STATE_LAYERS),
        }
    baseline_output = _finite_number(
        evaluate({name: baseline[name]["value"] for name in _STATE_LAYERS}),
        "localization baseline output", positive=False,
    )
    replay_output = _finite_number(
        evaluate({name: replacement[name]["value"] for name in _STATE_LAYERS}),
        "localization replay output", positive=False,
    )
    localization_dex = (
        math.log10(replay_output / baseline_output)
        if baseline_output > 0.0 and replay_output > 0.0 else None
    )
    return {
        "record_kind": "localization_control",
        "candidate": "whole_state_sentaurus_replay",
        "classification": "localization_control",
        "status": SampleStatus.VALID.value,
        "layers_replaced": list(_STATE_LAYERS),
        "baseline_output": baseline_output,
        "replay_output": replay_output,
        "downstream_localization_dex": localization_dex,
        "eligible_for_candidate_ranking": False,
        "eligible_for_formula_classification": False,
    }
