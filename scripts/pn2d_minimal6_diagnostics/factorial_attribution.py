"""Exact replacement-order and Shapley attribution helpers."""

from __future__ import annotations

import itertools
import math
from collections.abc import Mapping, Sequence


FACTORS = ("low_field", "drive", "support")


def replacement_order_contributions(
    values: Mapping[frozenset[str], float],
    factors: Sequence[str] = FACTORS,
) -> list[dict[str, object]]:
    expected = {
        frozenset(combo)
        for size in range(len(factors) + 1)
        for combo in itertools.combinations(factors, size)
    }
    if set(values) != expected:
        raise ValueError("factorial lattice is incomplete")
    if not all(math.isfinite(float(value)) for value in values.values()):
        raise ValueError("factorial values must be finite")
    output: list[dict[str, object]] = []
    for order in itertools.permutations(factors):
        active: frozenset[str] = frozenset()
        for step, factor in enumerate(order, start=1):
            updated = active | {factor}
            output.append(
                {
                    "order": ">".join(order),
                    "step": step,
                    "factor": factor,
                    "before": values[active],
                    "after": values[updated],
                    "increment": values[updated] - values[active],
                }
            )
            active = updated
    return output


def shapley_contributions(
    values: Mapping[frozenset[str], float],
    factors: Sequence[str] = FACTORS,
) -> dict[str, float]:
    increments = replacement_order_contributions(values, factors)
    order_count = math.factorial(len(factors))
    return {
        factor: sum(
            float(row["increment"])
            for row in increments
            if row["factor"] == factor
        )
        / order_count
        for factor in factors
    }


def interaction_remainder(
    values: Mapping[frozenset[str], float],
    factors: Sequence[str] = FACTORS,
) -> float:
    baseline = values[frozenset()]
    target = values[frozenset(factors)]
    main_effect = sum(
        values[frozenset((factor,))] - baseline for factor in factors
    )
    return target - baseline - main_effect


__all__ = [
    "FACTORS",
    "interaction_remainder",
    "replacement_order_contributions",
    "shapley_contributions",
]
