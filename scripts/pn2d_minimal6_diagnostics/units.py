"""Dimension-checked unit conversion for Minimal6 diagnostic quantities."""

from dataclasses import dataclass
import math


@dataclass(frozen=True)
class _Unit:
    dimension: str
    to_si: float


_UNITS = {
    "V/m": _Unit("electric_field", 1.0),
    "V/cm": _Unit("electric_field", 1.0e2),
    "A/m^2": _Unit("current_density", 1.0),
    "A/cm^2": _Unit("current_density", 1.0e4),
    "A*cm^-2": _Unit("current_density", 1.0e4),
    "m^-3": _Unit("number_density", 1.0),
    "cm^-3": _Unit("number_density", 1.0e6),
    "m^2/(V s)": _Unit("mobility", 1.0),
    "cm^2/(V s)": _Unit("mobility", 1.0e-4),
    "m^2*V^-1*s^-1": _Unit("mobility", 1.0),
    "cm^2*V^-1*s^-1": _Unit("mobility", 1.0e-4),
    "m^-1": _Unit("inverse_length", 1.0),
    "cm^-1": _Unit("inverse_length", 1.0e2),
    "m/s": _Unit("speed", 1.0),
    "cm/s": _Unit("speed", 1.0e-2),
    "m*s^-1": _Unit("speed", 1.0),
    "cm*s^-1": _Unit("speed", 1.0e-2),
    "m^-3*s^-1": _Unit("volume_generation_rate", 1.0),
    "cm^-3*s^-1": _Unit("volume_generation_rate", 1.0e6),
    "1": _Unit("dimensionless", 1.0),
}


def convert_value(value: float, source_unit: str, target_unit: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise TypeError("converted value must be numeric")
    if not math.isfinite(float(value)):
        raise ValueError("converted value must be finite")
    try:
        source = _UNITS[source_unit]
        target = _UNITS[target_unit]
    except KeyError as error:
        raise ValueError(f"unsupported unit: {error.args[0]}") from error
    if source.dimension != target.dimension:
        raise ValueError(f"dimension mismatch: {source_unit} -> {target_unit}")
    return float(value) * source.to_si / target.to_si
