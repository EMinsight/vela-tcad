"""Deterministic helpers for the Minimal6 internal-QFP SG replacement audit."""

from __future__ import annotations

import math
from collections.abc import Mapping


ELEMENTARY_CHARGE_C = 1.602176634e-19
INTERNAL_NODE_IDS = (1, 5)


def _carrier(value: str) -> str:
    carrier = value.strip().lower()
    if carrier not in {"electron", "hole"}:
        raise ValueError(f"unsupported carrier {value!r}")
    return carrier


def production_bernoulli(value: float) -> float:
    value = float(value)
    if not math.isfinite(value):
        raise ValueError("Bernoulli input must be finite")
    if abs(value) < 1.0e-10:
        return 1.0 - value * 0.5 + value * value / 12.0
    if value > 500.0:
        return value * math.exp(-value)
    if value < -500.0:
        return -value
    return value / math.expm1(value)


def _limited_exp(value: float) -> float:
    return math.exp(max(-500.0, min(500.0, float(value))))


def density_sg_flux(
    carrier: str,
    density0_m3: float,
    density1_m3: float,
    psi0_V: float,
    psi1_V: float,
    thermal_voltage_V: float,
    coefficient_m_per_s: float,
) -> float:
    carrier = _carrier(carrier)
    u = (float(psi1_V) - float(psi0_V)) / float(thermal_voltage_V)
    if carrier == "electron":
        signed = (
            production_bernoulli(-u) * float(density0_m3)
            - production_bernoulli(u) * float(density1_m3)
        )
    else:
        signed = (
            production_bernoulli(u) * float(density0_m3)
            - production_bernoulli(-u) * float(density1_m3)
        )
    return float(coefficient_m_per_s) * signed


def qf_sg_flux(
    carrier: str,
    ni0_m3: float,
    ni1_m3: float,
    psi0_V: float,
    psi1_V: float,
    qf0_V: float,
    qf1_V: float,
    thermal_voltage_V: float,
    coefficient_m_per_s: float,
) -> float:
    carrier = _carrier(carrier)
    ni0, ni1 = float(ni0_m3), float(ni1_m3)
    psi0, psi1 = float(psi0_V), float(psi1_V)
    qf0, qf1 = float(qf0_V), float(qf1_V)
    vt, coef = float(thermal_voltage_V), float(coefficient_m_per_s)
    if qf0 == qf1:
        return 0.0
    if ni0 <= 0.0 or ni1 <= 0.0:
        if carrier == "electron":
            density0 = ni0 * _limited_exp((psi0 - qf0) / vt)
            density1 = ni1 * _limited_exp((psi1 - qf1) / vt)
        else:
            density0 = ni0 * _limited_exp((qf0 - psi0) / vt)
            density1 = ni1 * _limited_exp((qf1 - psi1) / vt)
        return density_sg_flux(
            carrier, density0, density1, psi0, psi1, vt, coef
        )
    if carrier == "electron":
        eta = (psi1 - psi0) / vt + math.log(ni1 / ni0)
        density0 = ni0 * _limited_exp((psi0 - qf0) / vt)
        density1 = ni1 * _limited_exp((psi1 - qf1) / vt)
        signed = (
            production_bernoulli(-eta) * density0
            - production_bernoulli(eta) * density1
        )
    else:
        eta = (psi1 - psi0) / vt + math.log(ni0 / ni1)
        density0 = ni0 * _limited_exp((qf0 - psi0) / vt)
        density1 = ni1 * _limited_exp((qf1 - psi1) / vt)
        signed = (
            production_bernoulli(eta) * density0
            - production_bernoulli(-eta) * density1
        )
    return coef * signed


def replace_internal_qfp(
    vela_state: Mapping[int, Mapping[str, float]],
    sentaurus_state: Mapping[int, Mapping[str, float]],
    *,
    replace_electron: bool,
    replace_hole: bool,
) -> dict[int, dict[str, float]]:
    result = {int(node): dict(values) for node, values in vela_state.items()}
    for node in INTERNAL_NODE_IDS:
        if node not in result or node not in sentaurus_state:
            raise ValueError(f"replacement input lacks internal node {node}")
        if replace_electron:
            result[node]["phin_V"] = float(sentaurus_state[node]["phin_V"])
        if replace_hole:
            result[node]["phip_V"] = float(sentaurus_state[node]["phip_V"])
    return result


def continuity_flux_from_current_proxy(
    carrier: str, current_tangent_A_per_m2: float
) -> float:
    sign = -1.0 if _carrier(carrier) == "electron" else 1.0
    return sign * float(current_tangent_A_per_m2) / ELEMENTARY_CHARGE_C


def symmetric_relative_residual(candidate: float, reference: float) -> float:
    denominator = abs(candidate) + abs(reference)
    return 0.0 if denominator == 0.0 else abs(candidate - reference) / denominator


def absolute_log10_error(candidate: float, reference: float) -> float | None:
    if candidate == 0.0 or reference == 0.0:
        return None
    return abs(math.log10(abs(candidate) / abs(reference)))
