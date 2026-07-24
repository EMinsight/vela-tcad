"""Typed contracts and arithmetic for Minimal6 high-field box replay."""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence

from .mobility_diagnosis import FIELD, field_limited_mobility


BRANCH_ID = "sentaurus_lowfield_element_electric_field"
RECONSTRUCTION_LABEL = "box_operator_reconstruction"
CELL_MAPPING = {
    "mirror": (0, 1, 2, 3),
    "sketch": (0, 3, 2, 1),
}
REQUIRED_SAMPLE_FIELDS = (
    "topology",
    "bias_V",
    "carrier",
    "vela_triangle_id",
    "sentaurus_region_cell_id",
    "low_field_source_sha256",
    "electric_field_source_sha256",
    "saturation_velocity_m_per_s",
    "field_beta",
    "kappa",
    "status",
    "reconstruction_label",
)


def inverted_effective_field(
    carrier: str,
    low_field_mobility: float,
    final_mobility: float,
) -> float:
    if carrier not in FIELD:
        raise ValueError(f"unsupported carrier: {carrier}")
    low = float(low_field_mobility)
    final = float(final_mobility)
    if not math.isfinite(low) or low <= 0.0:
        raise ValueError("low-field mobility must be finite and positive")
    if not math.isfinite(final) or final <= 0.0:
        raise ValueError("final mobility must be finite and positive")
    if final > low:
        raise ValueError("final mobility exceeds low-field mobility")
    beta = FIELD[carrier]["beta"]
    power = max(0.0, (low / final) ** beta - 1.0)
    return (
        FIELD[carrier]["saturation_velocity"]
        / low
        * power ** (1.0 / beta)
    )


def coefficient_weighted_mobility(
    values: Sequence[tuple[float, float]],
) -> dict[str, object]:
    if not values:
        raise ValueError("coefficient-weighted mobility requires samples")
    normalized: list[tuple[float, float]] = []
    for kappa_value, mobility_value in values:
        kappa = float(kappa_value)
        mobility = float(mobility_value)
        if not math.isfinite(kappa):
            raise ValueError("kappa must be finite")
        if not math.isfinite(mobility) or mobility <= 0.0:
            raise ValueError("mobility must be finite and positive")
        normalized.append((kappa, mobility))
    kappa_sum = sum(kappa for kappa, _ in normalized)
    if kappa_sum == 0.0:
        return {
            "status": "geometric_zero",
            "mobility_m2_per_Vs": None,
            "kappa_sum": 0.0,
        }
    weighted = sum(
        kappa * mobility for kappa, mobility in normalized
    ) / kappa_sum
    if not math.isfinite(weighted) or weighted <= 0.0:
        raise ValueError("weighted mobility must be finite and positive")
    return {
        "status": "valid",
        "mobility_m2_per_Vs": weighted,
        "kappa_sum": kappa_sum,
    }


def validate_sample_record(record: Mapping[str, object]) -> None:
    missing = [
        field for field in REQUIRED_SAMPLE_FIELDS if field not in record
    ]
    if missing:
        raise ValueError(
            "sample record lacks required fields: " + ", ".join(missing)
        )
    topology = str(record["topology"])
    if topology not in CELL_MAPPING:
        raise ValueError(f"unsupported topology: {topology}")
    if record["reconstruction_label"] != RECONSTRUCTION_LABEL:
        raise ValueError("sample reconstruction label mismatch")
    if record["status"] not in (
        "valid",
        "geometric_zero",
        "reference_missing",
    ):
        raise ValueError("sample status is not typed")


__all__ = [
    "BRANCH_ID",
    "CELL_MAPPING",
    "FIELD",
    "RECONSTRUCTION_LABEL",
    "REQUIRED_SAMPLE_FIELDS",
    "coefficient_weighted_mobility",
    "field_limited_mobility",
    "inverted_effective_field",
    "validate_sample_record",
]
