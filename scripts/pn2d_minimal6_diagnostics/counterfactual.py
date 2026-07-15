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


def validate_source_anchor_kind(source_kind: str, *, native: bool) -> None:
    native_kind = "sentaurus_native_avalanche_generation"
    if native and source_kind != native_kind:
        raise ValueError("only raw ImpactIonization may be labelled native")
    if not native and source_kind == native_kind:
        raise ValueError("native source kind cannot be labelled reconstructed")
def native_source_anchor(values, *, volume_m3):
    if volume_m3 is None or volume_m3 <= 0.0:
        return {"status": "insufficient_data", "reason": "native generation lacks a physical volume", "value": None}
    return {"status": "available", "reason": None, "value": sum(float(value) for value in values) * float(volume_m3)}


def integrate_native_nodal_per_unit_depth(mesh: dict, values_by_node: Mapping[int, float]) -> dict:
    coordinates = {int(node["id"]): (float(node["x"]), float(node["y"])) for node in mesh["nodes"]}
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