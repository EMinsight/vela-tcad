from __future__ import annotations

import math
from collections.abc import Callable, Iterable, Mapping

from scripts.export_pn2d_minimal6_states import validate_state_matrix

FACTOR_DEPENDENCIES: dict[str, tuple[str, ...]] = {
    "ni_eff/BGN": (),
    "gradient_recovery": ("ni_eff/BGN",),
    "mobility": ("gradient_recovery",),
    "current_semantics": ("mobility",),
    "impact_driving_field": ("current_semantics",),
    "alpha_law": ("impact_driving_field",),
    "partial_volume": ("alpha_law",),
    "source_to_node_mapping": ("partial_volume",),
}


def validate_formula_input(manifest: dict) -> dict:
    if manifest.get("outputs_complete") is not True:
        raise ValueError("state package is incomplete")
    validate_state_matrix(manifest.get("states", []))
    return {"row_counts": {"node": 36, "edge": 54, "triangle": 24}, "sentaurus_internal_semantics_residual": None}


def validate_dependency_dag(dependencies: Mapping[str, Iterable[str]]) -> list[str]:
    names = set(dependencies)
    normalized = {factor: tuple(parents) for factor, parents in dependencies.items()}
    for factor, parents in normalized.items():
        unknown = set(parents) - names
        if unknown:
            raise ValueError(f"{factor} has undeclared dependencies: {sorted(unknown)}")
    order: list[str] = []
    visiting, visited = set(), set()
    def visit(factor: str) -> None:
        if factor in visited:
            return
        if factor in visiting:
            raise ValueError(f"counterfactual dependency cycle at {factor}")
        visiting.add(factor)
        for parent in normalized[factor]:
            visit(parent)
        visiting.remove(factor)
        visited.add(factor)
        order.append(factor)
    for factor in normalized:
        visit(factor)
    return order


def _validate_operator_value(value, path: str) -> None:
    if value is None:
        return
    if isinstance(value, bool):
        raise TypeError(f"{path} must not contain booleans")
    if isinstance(value, (int, float)):
        if not math.isfinite(float(value)):
            raise ValueError(f"{path} must be finite")
        return
    if isinstance(value, Mapping):
        for name, child in value.items():
            _validate_operator_value(child, f"{path}.{name}")
        return
    if isinstance(value, (tuple, list)):
        for index, child in enumerate(value):
            _validate_operator_value(child, f"{path}[{index}]")
        return
    raise TypeError(f"{path} has unsupported type {type(value).__name__}")

class _DeclaredInputs(Mapping[str, float]):
    def __init__(self, factor: str, allowed: tuple[str, ...], values: Mapping[str, float]):
        self._factor = factor
        self._allowed = allowed
        self._values = values

    def __getitem__(self, name: str) -> float:
        if name not in self._allowed:
            raise ValueError(f"{self._factor} accessed undeclared dependency {name}")
        return self._values[name]

    def __iter__(self):
        return iter(self._allowed)

    def __len__(self) -> int:
        return len(self._allowed)


class DependencyCounterfactualEngine:
    """Replay one named operator replacement through its declared DAG dependents."""

    def __init__(self, *, dependencies, baseline_values, replacement_values,
                 operators, output_factor):
        self.dependencies = {name: tuple(parents) for name, parents in dependencies.items()}
        self.order = validate_dependency_dag(self.dependencies)
        names = set(self.dependencies)
        for label, values in (("baseline", baseline_values),
                              ("replacement", replacement_values),
                              ("operator", operators)):
            if set(values) != names:
                raise ValueError(f"{label} values must cover exactly the declared operators")
        if output_factor not in names:
            raise ValueError("output factor is not declared")
        self.baseline_values = dict(baseline_values)
        self.replacement_values = dict(replacement_values)
        for name, value in (*self.baseline_values.items(), *self.replacement_values.items()):
            _validate_operator_value(value, f"operator input {name}")
        self.operators = dict(operators)
        self.output_factor = output_factor
        self.downstream = {}
        for changed in self.order:
            affected = {changed}
            for name in self.order:
                if any(parent in affected for parent in self.dependencies[name]):
                    affected.add(name)
            self.downstream[changed] = tuple(name for name in self.order if name in affected)

    def _compute(self, name: str, replaced: set[str], cache: dict):
        raw = self.replacement_values[name] if name in replaced else self.baseline_values[name]
        inputs = _DeclaredInputs(name, self.dependencies[name], cache)
        value = self.operators[name](inputs, raw)
        _validate_operator_value(value, f"operator result {name}")
        return value

    def _full(self, replaced: set[str]) -> dict[str, float]:
        unknown = set(replaced) - set(self.dependencies)
        if unknown:
            raise ValueError(f"undeclared replacement operators: {sorted(unknown)}")
        cache = {}
        for name in self.order:
            cache[name] = self._compute(name, replaced, cache)
        return cache

    def _advance(self, replaced: set[str], cache: dict[str, float], changed: str):
        recomputed = self.downstream[changed]
        updated = dict(cache)
        for name in recomputed:
            updated[name] = self._compute(name, replaced, updated)
        return updated, list(recomputed)

    def evaluate_replacements(self, replaced: set[str]) -> float:
        value = self._full(set(replaced))[self.output_factor]
        if value <= 0.0:
            raise ValueError("counterfactual output must be positive")
        return value

    def _path(self, sequence: Iterable[str]) -> dict:
        sequence = list(sequence)
        replaced: set[str] = set()
        cache = self._full(replaced)
        current = cache[self.output_factor]
        if current <= 0.0:
            raise ValueError("counterfactual baseline must be positive")
        contributions = []
        for name in sequence:
            if name in replaced:
                raise ValueError(f"operator {name} was replaced more than once")
            replaced.add(name)
            cache, recomputed = self._advance(replaced, cache, name)
            next_value = cache[self.output_factor]
            if next_value <= 0.0:
                raise ValueError("counterfactual output must be positive")
            contributions.append({
                "factor": name,
                "contribution_dex": math.log10(next_value / current),
                "recomputed": recomputed,
            })
            current = next_value
        return {"order": sequence, "contributions": contributions, "result": current}

    def evaluate_paths(self, *, native: float) -> dict:
        native = float(native)
        baseline = self.evaluate_replacements(set())
        if not math.isfinite(native) or native <= 0.0:
            raise ValueError("native source must be finite and positive")
        forward = self._path(self.order)
        reverse = self._path(reversed(self.order))
        native_gap = math.log10(native / baseline)
        residual = native_gap - sum(row["contribution_dex"] for row in forward["contributions"])
        assert_counterfactual_closure(
            native_gap_dex=native_gap,
            contributions_dex=(row["contribution_dex"] for row in forward["contributions"]),
            residual_dex=residual,
        )
        assert_counterfactual_closure(
            native_gap_dex=native_gap,
            contributions_dex=(row["contribution_dex"] for row in reverse["contributions"]),
            residual_dex=residual,
        )
        return {
            "dependency_order": list(self.order),
            "forward": forward,
            "reverse": reverse,
            "native_gap_dex": native_gap,
            "residual_dex": residual,
        }
def _numeric_tuple(raw, factor: str) -> tuple[float, ...]:
    if not isinstance(raw, (tuple, list)) or not raw:
        raise ValueError(f"{factor} requires a nonempty numeric sequence")
    values = tuple(float(value) for value in raw)
    if any(not math.isfinite(value) for value in values):
        raise ValueError(f"{factor} contains a non-finite value")
    return values


def _stage_with(parent: dict, key: str, raw, factor: str) -> dict:
    result = dict(parent)
    values = _numeric_tuple(raw, factor)
    if "density" in result and len(values) != len(result["density"]):
        raise ValueError(f"{factor} length does not match ni_eff/BGN")
    result[key] = values
    return result


def _formula_ni(_inputs, raw):
    return {"density": _numeric_tuple(raw, "ni_eff/BGN")}


def _formula_gradient(inputs, raw):
    return _stage_with(inputs["ni_eff/BGN"], "gradient", raw, "gradient_recovery")


def _formula_mobility(inputs, raw):
    return _stage_with(inputs["gradient_recovery"], "mobility", raw, "mobility")


def _formula_current(inputs, raw):
    stage = dict(inputs["mobility"])
    if raw is None:
        stage["flux"] = tuple(
            mobility * density * gradient
            for mobility, density, gradient in zip(
                stage["mobility"], stage["density"], stage["gradient"]
            )
        )
    else:
        stage["flux"] = _numeric_tuple(raw, "current_semantics")
        if len(stage["flux"]) != len(stage["density"]):
            raise ValueError("current_semantics length does not match ni_eff/BGN")
    return stage


def _formula_impact_field(inputs, raw):
    return _stage_with(inputs["current_semantics"], "impact_field", raw, "impact_driving_field")


def _alpha_from_parameters(field: float, parameters) -> float:
    if not isinstance(parameters, (tuple, list)) or len(parameters) not in (2, 5):
        raise ValueError("alpha_law parameters require (a,b) or (a_low,b_low,a_high,b_high,switch)")
    values = tuple(float(value) for value in parameters)
    if len(values) == 2:
        a_value, b_value = values
    else:
        a_low, b_low, a_high, b_high, switch = values
        a_value, b_value = (a_low, b_low) if abs(field) < switch else (a_high, b_high)
    magnitude = abs(float(field))
    if magnitude <= 0.0:
        return 0.0
    return a_value * math.exp(-b_value / magnitude)


def _formula_alpha(inputs, raw):
    stage = dict(inputs["impact_driving_field"])
    if not isinstance(raw, (tuple, list)) or len(raw) != len(stage["impact_field"]):
        raise ValueError("alpha_law requires one parameter tuple per contribution")
    alpha_values = []
    for field, parameters in zip(stage["impact_field"], raw):
        if isinstance(parameters, (int, float)) and not isinstance(parameters, bool):
            alpha = float(parameters)
            if not math.isfinite(alpha) or alpha < 0.0:
                raise ValueError("direct alpha_law values must be finite and nonnegative")
        else:
            alpha = _alpha_from_parameters(field, parameters)
        alpha_values.append(alpha)
    stage["alpha"] = tuple(alpha_values)
    return stage


def _formula_partial_volume(inputs, raw):
    stage = dict(inputs["alpha_law"])
    volumes = _numeric_tuple(raw, "partial_volume")
    if len(volumes) != len(stage["flux"]):
        raise ValueError("partial_volume length does not match source contributions")
    stage["local_source"] = tuple(
        alpha * abs(flux) * volume
        for alpha, flux, volume in zip(stage["alpha"], stage["flux"], volumes)
    )
    return stage


def _formula_source_mapping(inputs, raw):
    stage = inputs["partial_volume"]
    weights = _numeric_tuple(raw, "source_to_node_mapping")
    if len(weights) != len(stage["local_source"]):
        raise ValueError("source_to_node_mapping length does not match source contributions")
    return sum(source * weight for source, weight in zip(stage["local_source"], weights))


FORMULA_OPERATORS = {
    "ni_eff/BGN": _formula_ni,
    "gradient_recovery": _formula_gradient,
    "mobility": _formula_mobility,
    "current_semantics": _formula_current,
    "impact_driving_field": _formula_impact_field,
    "alpha_law": _formula_alpha,
    "partial_volume": _formula_partial_volume,
    "source_to_node_mapping": _formula_source_mapping,
}


def make_formula_operator_engine(*, baseline_values, replacement_values):
    return DependencyCounterfactualEngine(
        dependencies=FACTOR_DEPENDENCIES,
        baseline_values=baseline_values,
        replacement_values=replacement_values,
        operators=FORMULA_OPERATORS,
        output_factor="source_to_node_mapping",
    )


def evaluate_formula_counterfactual(*, native: float, baseline_values,
                                    replacement_values,
                                    unavailable_reasons=None) -> dict:
    unavailable_reasons = dict(unavailable_reasons or {})
    supplied = dict(replacement_values)
    completed = {
        factor: supplied.get(factor, baseline_values[factor])
        for factor in FACTOR_DEPENDENCIES
    }
    engine = make_formula_operator_engine(
        baseline_values=baseline_values,
        replacement_values=completed,
    )
    result = engine.evaluate_paths(native=native)
    result["factor_availability"] = [
        ({"factor": factor, "status": "available"}
         if factor in supplied else
         {"factor": factor, "status": "unavailable",
          "reason": unavailable_reasons.get(factor, "required operator input is absent")})
        for factor in FACTOR_DEPENDENCIES
    ]
    result["engine"] = engine
    return result


def evaluate_counterfactual_paths(*, native: float, baseline: float, factors: Mapping[str, float], dependencies: Mapping[str, Iterable[str]]) -> dict:
    if set(factors) != set(dependencies):
        raise ValueError("counterfactual factors and dependencies must name the same operators")
    if native <= 0.0 or baseline <= 0.0 or any(value <= 0.0 for value in factors.values()):
        raise ValueError("counterfactual sources and factors must be positive")
    order = validate_dependency_dag(dependencies)
    def path(sequence: Iterable[str]) -> dict:
        sequence = list(sequence)
        current = float(baseline)
        contributions = []
        for name in sequence:
            next_value = current * float(factors[name])
            contributions.append({"factor": name, "contribution_dex": math.log10(next_value / current)})
            current = next_value
        return {"order": list(sequence), "contributions": contributions, "result": current}
    forward = path(order)
    reverse = path(reversed(order))
    native_gap = math.log10(native / baseline)
    residual = native_gap - sum(row["contribution_dex"] for row in forward["contributions"])
    return {"dependency_order": order, "forward": forward, "reverse": reverse, "native_gap_dex": native_gap, "residual_dex": residual}


def interaction_dex(*, baseline: float, a_only: float, b_only: float, both: float) -> float:
    if any(value <= 0.0 for value in (baseline, a_only, b_only, both)):
        raise ValueError("interaction sources must be positive")
    return math.log10(both) - math.log10(a_only) - math.log10(b_only) + math.log10(baseline)


def build_adjacent_interactions(forward: list[dict], reverse: list[dict], evaluate_replacements: Callable[[set[str]], float], *, threshold_dex: float = 0.3) -> list[dict]:
    forward_by_name = {row["factor"]: float(row["contribution_dex"]) for row in forward}
    reverse_by_name = {row["factor"]: float(row["contribution_dex"]) for row in reverse}
    result = []
    for first, second in zip(forward, forward[1:]):
        a, b = first["factor"], second["factor"]
        if max(abs(forward_by_name[a] - reverse_by_name[a]), abs(forward_by_name[b] - reverse_by_name[b])) <= threshold_dex:
            continue
        baseline = evaluate_replacements(set())
        a_only = evaluate_replacements({a})
        b_only = evaluate_replacements({b})
        both = evaluate_replacements({a, b})
        result.append({
            "first_factor": a, "second_factor": b, "path_identity": "forward_adjacent",
            "baseline": baseline, "a_only": a_only, "b_only": b_only, "both": both,
            "interaction_dex": interaction_dex(baseline=baseline, a_only=a_only, b_only=b_only, both=both),
        })
    for first, second in zip(reverse, reverse[1:]):
        a, b = first["factor"], second["factor"]
        if max(abs(forward_by_name[a] - reverse_by_name[a]), abs(forward_by_name[b] - reverse_by_name[b])) <= threshold_dex:
            continue
        baseline = evaluate_replacements(set())
        a_only = evaluate_replacements({a})
        b_only = evaluate_replacements({b})
        both = evaluate_replacements({a, b})
        result.append({
            "first_factor": a, "second_factor": b, "path_identity": "reverse_adjacent",
            "baseline": baseline, "a_only": a_only, "b_only": b_only, "both": both,
            "interaction_dex": interaction_dex(baseline=baseline, a_only=a_only, b_only=b_only, both=both),
        })
    return result


def symmetric_contributions(forward: list[dict], reverse: list[dict]) -> dict[str, float]:
    forward_by_name = {row["factor"]: float(row["contribution_dex"]) for row in forward}
    reverse_by_name = {row["factor"]: float(row["contribution_dex"]) for row in reverse}
    if set(forward_by_name) != set(reverse_by_name):
        raise ValueError("forward and reverse paths do not cover the same factors")
    return {name: 0.5 * (forward_by_name[name] + reverse_by_name[name]) for name in forward_by_name}


def assert_counterfactual_closure(*, native_gap_dex: float, contributions_dex: Iterable[float], residual_dex: float, tolerance_dex: float = 1.0e-10) -> None:
    closure = sum(float(value) for value in contributions_dex) + float(residual_dex)
    if abs(float(native_gap_dex) - closure) > float(tolerance_dex):
        raise ValueError(f"counterfactual closure failed: gap={native_gap_dex}, closure={closure}")


def score_dominance(states: list[dict]) -> dict:
    eligible = [state for state in states if float(state["bias_V"]) in (-12.0, -19.0) and abs(float(state["native_gap_dex"])) > 0.0]
    required = {(topology, bias) for topology in ("sketch", "mirror") for bias in (-12.0, -19.0)}
    present = {(state["topology"], float(state["bias_V"])) for state in eligible}
    if present != required:
        return {"status": "insufficient_data", "reason": "dominance requires both topologies at -12 V and -19 V"}
    for state in eligible:
        unavailable = [
            row for row in state.get("factor_availability", [])
            if row.get("status") != "available"
        ]
        if unavailable:
            names = ", ".join(sorted(str(row.get("factor")) for row in unavailable))
            return {"status": "insufficient_data", "reason": f"unavailable factors: {names}"}
        if abs(float(state["residual_dex"])) > 0.25 * abs(float(state["native_gap_dex"])):
            return {"status": "insufficient_data", "reason": "residual exceeds 25 percent of source gap"}
    common = set.intersection(*(set(state["symmetric_contributions"]) for state in eligible))
    if not common:
        return {"status": "insufficient_data", "reason": "no common factor contributions"}
    scores = {factor: sum(abs(float(state["symmetric_contributions"][factor])) for state in eligible) / len(eligible) for factor in common}
    dominant = max(scores, key=lambda factor: (scores[factor], factor))
    return {"status": "available", "dominant_factor": dominant, "scores": scores}


def validate_field_units(fields: Iterable[Mapping[str, str]], expected: Mapping[str, str]) -> None:
    by_name: dict[str, set[str]] = {}
    for field in fields:
        by_name.setdefault(str(field["name"]), set()).add(str(field["unit"]))
    for name, unit in expected.items():
        if unit not in by_name.get(name, set()):
            raise ValueError(f"field {name} lacks required unit {unit}")


def validate_source_anchor_kind(source_family: str, source_kind: str, *, native: bool) -> None:
    native_family = "sentaurus_native_avalanche_generation"
    expected_kind = "sentaurus" if native else "derived"
    if native and source_family != native_family:
        raise ValueError("only raw ImpactIonization may be labelled native")
    if not native and source_family == native_family:
        raise ValueError("native source family cannot be labelled reconstructed")
    if source_kind != expected_kind:
        raise ValueError(
            f"source family {source_family} has inconsistent SourceKind {source_kind}"
        )
def native_source_anchor(values, *, volume_m3):
    if volume_m3 is None or volume_m3 <= 0.0:
        return {"status": "insufficient_data", "reason": "native generation lacks a physical volume", "value": None}
    return {"status": "available", "reason": None, "value": sum(float(value) for value in values) * float(volume_m3)}


def integrate_native_nodal_per_unit_depth(mesh: dict, values_by_node: Mapping[int, float]) -> dict:
    coordinate_unit = str(mesh.get("coordinate_unit", "m")).strip().lower()
    coordinate_scale_to_m = {
        "m": 1.0,
        "meter": 1.0,
        "metre": 1.0,
        "um": 1.0e-6,
        "micrometer": 1.0e-6,
        "micrometre": 1.0e-6,
    }.get(coordinate_unit)
    if coordinate_scale_to_m is None:
        raise ValueError(
            f"unsupported native source mesh coordinate_unit {coordinate_unit}"
        )
    coordinates = {
        int(node["id"]): (
            float(node["x"]) * coordinate_scale_to_m,
            float(node["y"]) * coordinate_scale_to_m,
        )
        for node in mesh["nodes"]
    }
    if len(coordinates) != len(mesh["nodes"]):
        raise ValueError("native source mesh contains duplicate node IDs")
    total = 0.0
    for triangle in mesh["triangles"]:
        ids = [int(value) for value in triangle["node_ids"]]
        if len(ids) != 3 or any(node not in coordinates or node not in values_by_node for node in ids):
            raise ValueError("triangle lacks node coordinate or native source")
        (x0, y0), (x1, y1), (x2, y2) = (coordinates[node] for node in ids)
        twice_area_m2 = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)
        if twice_area_m2 <= 0.0:
            raise ValueError("native source mesh requires CCW triangles")
        total += 0.5 * twice_area_m2 * 1.0e4 * sum(float(values_by_node[node]) for node in ids) / 3.0
    return {"status": "available", "value_s_inv_per_unit_depth": total, "depth_convention": "unit_out_of_plane_length_cm"}


def integrate_vela_reconstructed_per_unit_depth(rows: Iterable[Mapping[str, str]]) -> float:
    total_per_m_s = 0.0
    for row in rows:
        for key, value in row.items():
            if key.endswith("_electron_source_integral_per_m_s") or key.endswith("_hole_source_integral_per_m_s"):
                total_per_m_s += float(value)
    return total_per_m_s * 1.0e-2


def sentaurus_alpha_current_nodal(electron_alpha_cm_inv, electron_current_A_per_cm2, hole_alpha_cm_inv, hole_current_A_per_cm2, *, elementary_charge_C: float = 1.602176634e-19):
    if elementary_charge_C <= 0.0:
        raise ValueError("elementary charge must be positive")
    keys = set(electron_alpha_cm_inv) | set(electron_current_A_per_cm2) | set(hole_alpha_cm_inv) | set(hole_current_A_per_cm2)
    result = {}
    for node in keys:
        ea, ha = float(electron_alpha_cm_inv.get(node, 0.0)), float(hole_alpha_cm_inv.get(node, 0.0))
        ex, ey = electron_current_A_per_cm2.get(node, (0.0, 0.0))
        hx, hy = hole_current_A_per_cm2.get(node, (0.0, 0.0))
        result[node] = (ea * math.hypot(float(ex), float(ey)) + ha * math.hypot(float(hx), float(hy))) / elementary_charge_C
    return result


def source_log_gap(native: float, reconstructed: float) -> dict:
    native, reconstructed = abs(float(native)), abs(float(reconstructed))
    if native <= 1.0e-285 and reconstructed <= 1.0e-285:
        return {"classification": "geometric_zero", "dex": None}
    if native <= 0.0 or reconstructed <= 0.0:
        return {"classification": "unavailable", "dex": None}
    return {"classification": "available", "dex": math.log10(native / reconstructed)}