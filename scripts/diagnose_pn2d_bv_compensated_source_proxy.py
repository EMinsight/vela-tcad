#!/usr/bin/env python3
"""Compare PN2D BV baseline vs compensated-junction source proxy factors.

The diagnostic focuses on the three horizontal junction cuts used by the
coarse7x3 BV debug artifacts. It intentionally does not alter solver state.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import statistics
import sys
from decimal import Decimal, localcontext
from pathlib import Path
from typing import Any


ELEMENTARY_CHARGE_C = 1.602176634e-19
BIASES = [-12.0, -19.0, -20.0]
Y_CUTS = [0.0, 0.25, 0.5]
REPLAY_VARIANT_MATRIX = {
    "legacy_density_gradient": {
        "doping_strategy": "legacy",
        "compensated_doping_policy": "dominant_signed_region",
        "current_variant": "density_gradient",
        "current_approximation": "density_gradient",
    },
    "legacy_gss_midpoint": {
        "doping_strategy": "legacy",
        "compensated_doping_policy": "dominant_signed_region",
        "current_variant": "gss_midpoint",
        "current_approximation": "cell_reconstructed",
    },
    "reported_density_gradient": {
        "doping_strategy": "reported",
        "compensated_doping_policy": "reported",
        "current_variant": "density_gradient",
        "current_approximation": "density_gradient",
    },
    "reported_gss_midpoint": {
        "doping_strategy": "reported",
        "compensated_doping_policy": "reported",
        "current_variant": "gss_midpoint",
        "current_approximation": "cell_reconstructed",
    },
}
REPLAY_VARIANTS = tuple(REPLAY_VARIANT_MATRIX)
BOLTZMANN_OVER_CHARGE_V_K = 8.617333262145e-5
ELECTRON_SG_FIELDS = (
    "electron_sg_ni0",
    "electron_sg_ni1",
    "electron_sg_n0",
    "electron_sg_n1",
    "electron_sg_psi0",
    "electron_sg_psi1",
    "electron_sg_phin0",
    "electron_sg_phin1",
    "electron_sg_eta",
    "electron_sg_b_minus_eta",
    "electron_sg_b_eta",
    "electron_sg_coef",
    "electron_sg_left_term",
    "electron_sg_right_term",
    "electron_sg_signed_difference",
    "electron_sg_reconstructed_flux_native",
    "electron_sg_stable_factorized_flux_native",
    "electron_sg_production_signed_flux_native",
    "electron_sg_cancellation_condition",
    "electron_sg_node0_exponent_clamped_low",
    "electron_sg_node0_exponent_clamped_high",
    "electron_sg_node1_exponent_clamped_low",
    "electron_sg_node1_exponent_clamped_high",
    "electron_sg_include_ni_gradient_drift",
    "electron_sg_flat_qf_short_circuit",
    "electron_sg_reconstruction_relative_error",
    "electron_sg_high_precision_reference_flux_native",
    "electron_sg_production_vs_high_precision_reference_relative_error",
    "electron_sg_stable_vs_high_precision_reference_relative_error",
    "electron_sg_production_signed_continuity_particle_flux_m2_s",
    "electron_sg_production_abs_continuity_particle_flux_m2_s",
    "electron_sg_production_signed_conventional_current_density_A_per_m2",
    "electron_sg_production_signed_conventional_current_density_A_per_cm2",
)
REQUIRED_ENRICHED_FIELDS = ELECTRON_SG_FIELDS + (
    "sentaurus_e_psi0_V",
    "sentaurus_e_psi1_V",
    "sentaurus_e_phin0_V",
    "sentaurus_e_phin1_V",
    "sentaurus_e_density0_m3",
    "sentaurus_e_density1_m3",
    "sentaurus_e_mobility0_m2_V_s",
    "sentaurus_e_mobility1_m2_V_s",
    "sentaurus_e_ni_inferred0_m3",
    "sentaurus_e_ni_inferred1_m3",
    "sentaurus_e_alpha0_m_inv",
    "sentaurus_e_alpha1_m_inv",
    "sentaurus_e_alpha_edge_average_m_inv",
    "vela_e_over_sentaurus_alpha_abs_ratio",
    "sentaurus_e_current_edge_signed_A_cm2",
    "sentaurus_edge_length_m",
    "sentaurus_e_continuity_edge_signed_flux_m2_s",
    "sentaurus_e_sg_vela_mobility_signed_flux_m2_s",
    "sentaurus_e_sg_vela_mobility_abs_flux_m2_s",
    "sentaurus_e_sg_vela_mobility_conventional_current_A_cm2",
    "sentaurus_e_sg_sentaurus_mobility_signed_flux_m2_s",
    "sentaurus_e_sg_sentaurus_mobility_abs_flux_m2_s",
    "sentaurus_e_sg_sentaurus_mobility_conventional_current_A_cm2",
    "vela_e_sg_production_canonical_signed_flux_m2_s",
    "vela_e_over_sentaurus_vector_abs_ratio",
    "sentaurus_e_sg_vela_mobility_over_vector_abs_ratio",
    "sentaurus_e_sg_sentaurus_mobility_over_vector_abs_ratio",
    "vela_e_source_integral_physical_m_inv_s",
    "sentaurus_e_source_on_vela_area_physical_m_inv_s",
    "vela_e_over_sentaurus_source_abs_ratio",
    "vela_e_source_closure_ratio",
)
SENTAURUS_SCALAR_FIELD_SPECS = {
    "ElectrostaticPotential": "V",
    "eQuasiFermiPotential": "V",
    "eDensity": "cm^-3",
    "eMobility": "cm^2*V^-1*s^-1",
    "eAlphaAvalanche": "cm^-1",
}
EDGE_BY_SIDE = {
    0.0: {"left": 9, "right": 13},
    0.25: {"left": 12, "right": 16},
    0.5: {"left": 34, "right": 37},
}
VELA_X_COLUMNS = [2.0 / 3.0, 1.0, 4.0 / 3.0]
SENTAURUS_X_COLUMNS = [0.75, 1.0, 1.25]
NODE_FIELDS = ["Potential", "ElectronQuasiFermi", "HoleQuasiFermi", "Electrons", "Holes"]
FIELD_TO_OUTPUT = {
    "Potential": "psi",
    "ElectronQuasiFermi": "phin",
    "HoleQuasiFermi": "phip",
    "Electrons": "electrons",
    "Holes": "holes",
}
RATIO_FIELDS = [
    "edge_source_integral",
    "electron_source_integral",
    "hole_source_integral",
    "electron_alpha_m_inv",
    "hole_alpha_m_inv",
    "electron_flux_proxy",
    "hole_flux_proxy",
    "electron_raw_flux_proxy",
    "hole_raw_flux_proxy",
    "electron_mobility_m2_V_s",
    "hole_mobility_m2_V_s",
    "electron_density_mid_m3",
    "hole_density_mid_m3",
    "edge_area_proxy_m2",
]


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--variants-root", type=Path)
    parser.add_argument("--baseline-report-root", type=Path)
    parser.add_argument("--probe-root", type=Path)
    parser.add_argument("--sentaurus-root", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    has_legacy_pair = args.baseline_report_root is not None and args.probe_root is not None
    if args.variants_root is None and not has_legacy_pair:
        parser.error(
            "provide --variants-root or both --baseline-report-root and --probe-root"
        )
    if (args.baseline_report_root is None) != (args.probe_root is None):
        parser.error("--baseline-report-root and --probe-root must be provided together")
    return args


def long_path(path: Path) -> str:
    resolved = str(path.resolve())
    if os.name == "nt" and not resolved.startswith("\\\\?\\"):
        return "\\\\?\\" + resolved
    return resolved


def open_path(path: Path, mode: str, **kwargs: Any) -> Any:
    return open(long_path(path), mode, **kwargs)


def read_csv(path: Path) -> list[dict[str, str]]:
    with open_path(path, "r", newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise SystemExit(f"no rows to write: {path}")
    with open_path(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


def finite_float(raw: Any, default: float = math.nan) -> float:
    if raw in (None, ""):
        return default
    try:
        value = float(raw)
    except (TypeError, ValueError):
        return default
    return value if math.isfinite(value) else default


def safe_ratio(numerator: float, denominator: float) -> float | None:
    if not math.isfinite(numerator) or not math.isfinite(denominator) or denominator == 0.0:
        return None
    return numerator / denominator


def abs_ratio(numerator: float, denominator: float) -> float | None:
    if not math.isfinite(numerator) or not math.isfinite(denominator) or denominator == 0.0:
        return None
    return abs(numerator) / abs(denominator)


def log10_abs(value: float | None) -> float | None:
    if value is None or not math.isfinite(value) or value == 0.0:
        return None
    return math.log10(abs(value))


def clean_json(value: Any) -> Any:
    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        return None
    if isinstance(value, dict):
        return {key: clean_json(item) for key, item in value.items()}
    if isinstance(value, list):
        return [clean_json(item) for item in value]
    return value


def production_bernoulli(value: float) -> float:
    """Match Vela's production Bernoulli branch thresholds."""
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


def _finite_saturated_float(value: Decimal) -> float:
    maximum = Decimal(str(sys.float_info.max))
    if value.is_nan():
        raise ValueError("Decimal SG reference produced NaN")
    if value >= maximum:
        return sys.float_info.max
    if value <= -maximum:
        return -sys.float_info.max
    result = float(value)
    if not math.isfinite(result):
        return math.copysign(sys.float_info.max, result)
    return result


def _decimal_bernoulli(value: Decimal) -> Decimal:
    if abs(value) < Decimal("1e-30"):
        return Decimal(1) - value / 2 + value * value / 12
    if value > 500:
        return value * (-value).exp()
    if value < -500:
        return -value
    return value / (value.exp() - 1)


def replay_electron_variable_ni_sg(
    *,
    ni0: float,
    ni1: float,
    psi0: float,
    psi1: float,
    phin0: float,
    phin1: float,
    vt: float,
    mobility_m2_V_s: float,
    length_m: float,
    include_ni_gradient_drift: bool = True,
) -> dict[str, Any]:
    """Replay Vela's variable-ni electron SG flux with a Decimal reference."""
    raw = {
        "ni0": ni0,
        "ni1": ni1,
        "psi0": psi0,
        "psi1": psi1,
        "phin0": phin0,
        "phin1": phin1,
        "vt": vt,
        "mobility_m2_V_s": mobility_m2_V_s,
        "length_m": length_m,
    }
    invalid = [name for name, value in raw.items() if not math.isfinite(float(value))]
    if invalid:
        raise ValueError(f"SG replay inputs must be finite: {invalid}")
    if vt <= 0.0 or length_m <= 0.0:
        raise ValueError("SG replay vt and length_m must be positive")
    if ni0 < 0.0 or ni1 < 0.0 or mobility_m2_V_s < 0.0:
        raise ValueError("SG replay ni and mobility must be non-negative")

    coefficient = mobility_m2_V_s * vt / length_m
    endpoint_exponent0 = (psi0 - phin0) / vt
    endpoint_exponent1 = (psi1 - phin1) / vt
    clamped_exponent0 = max(-500.0, min(500.0, endpoint_exponent0))
    clamped_exponent1 = max(-500.0, min(500.0, endpoint_exponent1))
    density0 = ni0 * math.exp(clamped_exponent0)
    density1 = ni1 * math.exp(clamped_exponent1)
    eta = (psi1 - psi0) / vt
    if ni0 > 0.0 and ni1 > 0.0 and include_ni_gradient_drift:
        eta += math.log(ni1 / ni0)
    b_minus_eta = production_bernoulli(-eta)
    b_eta = production_bernoulli(eta)
    left_term = b_minus_eta * density0
    right_term = b_eta * density1
    flat_qf = phin0 == phin1
    signed_difference = 0.0 if flat_qf else left_term - right_term
    double_flux = coefficient * signed_difference
    term_scale = abs(left_term) + abs(right_term)
    if term_scale == 0.0:
        cancellation_condition = 0.0
    elif signed_difference == 0.0:
        cancellation_condition = sys.float_info.max
    else:
        cancellation_condition = min(
            sys.float_info.max,
            term_scale / abs(signed_difference),
        )

    with localcontext() as context:
        context.prec = 100
        dec = {name: Decimal(str(value)) for name, value in raw.items()}
        dec_endpoint0 = (dec["psi0"] - dec["phin0"]) / dec["vt"]
        dec_endpoint1 = (dec["psi1"] - dec["phin1"]) / dec["vt"]
        dec_clamped0 = max(Decimal(-500), min(Decimal(500), dec_endpoint0))
        dec_clamped1 = max(Decimal(-500), min(Decimal(500), dec_endpoint1))
        dec_density0 = dec["ni0"] * dec_clamped0.exp()
        dec_density1 = dec["ni1"] * dec_clamped1.exp()
        dec_eta = (dec["psi1"] - dec["psi0"]) / dec["vt"]
        if ni0 > 0.0 and ni1 > 0.0 and include_ni_gradient_drift:
            dec_eta += (dec["ni1"] / dec["ni0"]).ln()
        dec_left = _decimal_bernoulli(-dec_eta) * dec_density0
        dec_right = _decimal_bernoulli(dec_eta) * dec_density1
        dec_coefficient = (
            dec["mobility_m2_V_s"] * dec["vt"] / dec["length_m"]
        )
        if flat_qf:
            decimal_flux_raw = Decimal(0)
        elif ni0 > 0.0 and ni1 > 0.0:
            log_left_over_right = (
                (dec["phin1"] - dec["phin0"]) / dec["vt"]
                + (dec_clamped0 - dec_endpoint0)
                - (dec_clamped1 - dec_endpoint1)
            )
            if not include_ni_gradient_drift:
                log_left_over_right += (dec["ni0"] / dec["ni1"]).ln()
            decimal_flux_raw = (
                dec_coefficient
                * dec_right
                * (log_left_over_right.exp() - Decimal(1))
            )
        else:
            decimal_flux_raw = dec_coefficient * (dec_left - dec_right)
        decimal_term_scale_raw = abs(dec_coefficient) * (
            abs(dec_left) + abs(dec_right)
        )
        decimal_flux = _finite_saturated_float(decimal_flux_raw)
        decimal_term_scale = _finite_saturated_float(decimal_term_scale_raw)

    reference_scale = max(abs(decimal_flux), 1.0e-300)
    double_highprec_relative_error = abs(double_flux - decimal_flux) / reference_scale
    if not math.isfinite(double_highprec_relative_error):
        double_highprec_relative_error = sys.float_info.max
    return {
        "ni0": ni0,
        "ni1": ni1,
        "n0": density0,
        "n1": density1,
        "psi0": psi0,
        "psi1": psi1,
        "phin0": phin0,
        "phin1": phin1,
        "eta": eta,
        "bernoulli_minus_eta": b_minus_eta,
        "bernoulli_eta": b_eta,
        "coef": coefficient,
        "left_term": left_term,
        "right_term": right_term,
        "signed_difference": signed_difference,
        "double_flux_m2_s": double_flux,
        "decimal_flux_m2_s": decimal_flux,
        "decimal_reference_term_scale": decimal_term_scale,
        "double_highprec_relative_error": double_highprec_relative_error,
        "cancellation_condition": cancellation_condition,
        "node0_exponent_clamped_low": endpoint_exponent0 < -500.0,
        "node0_exponent_clamped_high": endpoint_exponent0 > 500.0,
        "node1_exponent_clamped_low": endpoint_exponent1 < -500.0,
        "node1_exponent_clamped_high": endpoint_exponent1 > 500.0,
        "include_ni_gradient_drift": include_ni_gradient_drift,
        "flat_qf_short_circuit": flat_qf,
    }


def require_finite_fields(
    row: dict[str, Any],
    fields: tuple[str, ...] | list[str],
    *,
    context: str,
) -> None:
    """Require numeric, finite values instead of silently emitting JSON nulls."""
    missing = [field for field in fields if field not in row or row[field] in (None, "")]
    if missing:
        raise ValueError(f"{context}: missing required fields: {missing}")
    nonfinite: list[str] = []
    for field in fields:
        try:
            value = float(row[field])
        except (TypeError, ValueError):
            nonfinite.append(field)
            continue
        if not math.isfinite(value):
            nonfinite.append(field)
    if nonfinite:
        raise ValueError(f"{context}: non-finite required fields: {nonfinite}")


def validate_72_row_keys(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Validate the exact 3 bias x 3 cut x 2 side x 4 variant matrix."""
    expected = {
        (variant, float(bias), float(y_um), side)
        for variant in REPLAY_VARIANTS
        for bias in BIASES
        for y_um in Y_CUTS
        for side in ("left", "right")
    }
    keys: list[tuple[str, float, float, str]] = []
    for index, row in enumerate(rows):
        try:
            key = (
                str(row["variant"]),
                float(row["bias_V"]),
                float(row["y_um"]),
                str(row["side"]),
            )
        except (KeyError, TypeError, ValueError) as exc:
            raise ValueError(f"row {index}: invalid replay key") from exc
        keys.append(key)
    if len(keys) != len(set(keys)):
        raise ValueError("duplicate row key in compensated SG replay matrix")
    if len(keys) != 72:
        raise ValueError(f"expected 72 replay rows, got {len(keys)}")
    actual = set(keys)
    if actual != expected:
        missing = sorted(expected - actual)
        extra = sorted(actual - expected)
        raise ValueError(f"replay row key mismatch: missing={missing}, extra={extra}")
    return rows


def standard_variant_inputs(variants_root: Path) -> dict[str, dict[str, Any]]:
    """Resolve standardized current-HEAD two-by-two replay inputs."""
    result: dict[str, dict[str, Any]] = {}
    for name, metadata in REPLAY_VARIANT_MATRIX.items():
        variant_root = variants_root / name
        run_root = variant_root / "run"
        current_variant = metadata["current_variant"]
        result[name] = {
            "name": name,
            "variant_root": variant_root,
            "run_root": run_root,
            "doping_csv": variant_root / "imported" / "vela" / "doping.csv",
            "deck_path": run_root / f"simulation_pn2d_bv_{current_variant}.json",
            "sg_csv": run_root / f"sg_avalanche_edges_{current_variant}.csv",
            "vtk_root": run_root / "vtk" / current_variant,
            "implementation": "current_head",
            **metadata,
        }
    return result


def variant_run_status(spec: dict[str, Any]) -> dict[str, Any]:
    result = {
        "variant": str(spec["name"]),
        "doping_strategy": str(spec["doping_strategy"]),
        "current_variant": str(spec["current_variant"]),
        "run_status": "missing_deck",
        "last_converged_bias_V": None,
        "handoff_stage": "",
    }
    deck_path = Path(spec["deck_path"])
    if not deck_path.is_file():
        return result
    deck = json.loads(deck_path.read_text(encoding="utf-8-sig"))
    output_csv = Path(str(deck.get("output_csv", "")))
    if not output_csv.is_absolute():
        output_csv = deck_path.parent / output_csv
    if not output_csv.is_file():
        result["run_status"] = "prepared"
        return result
    rows = read_csv(output_csv)
    if not rows:
        result["run_status"] = "empty"
        return result
    converged_rows = [
        row for row in rows
        if str(row.get("converged", "")).strip().lower() in {"1", "true", "yes"}
    ]
    if converged_rows:
        result["last_converged_bias_V"] = float(converged_rows[-1]["bias_V"])
    result["handoff_stage"] = str(rows[-1].get("handoff_stage", ""))
    last_converged = result["last_converged_bias_V"]
    result["run_status"] = (
        "complete"
        if last_converged is not None and last_converged <= -20.0 + 1.0e-12
        else "partial"
    )
    return result


def validate_manifest_vector_components(
    manifest: dict[str, Any],
    field_name: str,
    *,
    expected_components: int,
) -> dict[str, Any]:
    """Return a manifest field only when it has the required vector arity."""
    fields = manifest.get("fields")
    if not isinstance(fields, list):
        raise ValueError("field manifest is missing a fields list")
    name_matches = [
        field for field in fields
        if isinstance(field, dict) and str(field.get("name", "")) == field_name
    ]
    if not name_matches:
        raise ValueError(f"field manifest missing {field_name}")
    matches: list[dict[str, Any]] = []
    for field in name_matches:
        try:
            components = int(field.get("components"))
        except (TypeError, ValueError):
            continue
        if components == expected_components:
            matches.append(field)
    if not matches:
        actual = [field.get("components") for field in name_matches]
        raise ValueError(
            f"{field_name} must have components={expected_components}, got {actual}"
        )
    if len(matches) != 1:
        raise ValueError(
            f"field manifest has duplicate {field_name} components={expected_components} entries"
        )
    field = matches[0]
    return field


def _require_manifest_field(
    manifest: dict[str, Any],
    field_name: str,
    *,
    components: int,
    unit: str,
) -> dict[str, Any]:
    fields = manifest.get("fields")
    if not isinstance(fields, list):
        raise ValueError("field manifest is missing a fields list")
    matches = [
        field for field in fields
        if (
            isinstance(field, dict)
            and str(field.get("name", "")) == field_name
            and int(field.get("components", -1)) == components
            and int(field.get("region", -1)) == 0
        )
    ]
    if len(matches) != 1:
        raise ValueError(
            f"field manifest requires exactly one region0 {field_name} "
            f"with components={components}"
        )
    field = matches[0]
    if str(field.get("unit", "")) != unit:
        raise ValueError(f"{field_name} unit must be {unit}")
    if str(field.get("mapping_status", "")) != "complete":
        raise ValueError(f"{field_name} mapping_status must be complete")
    if str(field.get("global_node_mapping", "")) != "global_vertex_order":
        raise ValueError(f"{field_name} must use global_vertex_order mapping")
    return field


def project_endpoint_current_to_canonical_edge(
    *,
    point0: tuple[float, float],
    point1: tuple[float, float],
    current0_A_cm2: tuple[float, float],
    current1_A_cm2: tuple[float, float],
) -> dict[str, Any]:
    """Average a two-component conventional current and project low-x/low-y to high."""
    raw_values = (*point0, *point1, *current0_A_cm2, *current1_A_cm2)
    if not all(math.isfinite(float(value)) for value in raw_values):
        raise ValueError("edge projection inputs must be finite")
    p0 = (float(point0[0]), float(point0[1]))
    p1 = (float(point1[0]), float(point1[1]))
    j0 = (float(current0_A_cm2[0]), float(current0_A_cm2[1]))
    j1 = (float(current1_A_cm2[0]), float(current1_A_cm2[1]))
    reversed_input = p1 < p0
    if reversed_input:
        p0, p1 = p1, p0
        j0, j1 = j1, j0
    dx = p1[0] - p0[0]
    dy = p1[1] - p0[1]
    length = math.hypot(dx, dy)
    if length <= 0.0:
        raise ValueError("canonical edge projection requires distinct endpoints")
    tangent = (dx / length, dy / length)
    average_current = (0.5 * (j0[0] + j1[0]), 0.5 * (j0[1] + j1[1]))
    conventional_current = (
        average_current[0] * tangent[0] + average_current[1] * tangent[1]
    )
    electron_continuity_flux = (
        -conventional_current * 1.0e4 / ELEMENTARY_CHARGE_C
    )
    return {
        "canonical_point0": p0,
        "canonical_point1": p1,
        "canonical_tangent": tangent,
        "canonical_length_coordinate_units": length,
        "input_orientation_reversed": reversed_input,
        "conventional_current_A_cm2": conventional_current,
        "conventional_current_abs_A_cm2": abs(conventional_current),
        "electron_continuity_flux_m2_s": electron_continuity_flux,
        "electron_continuity_flux_abs_m2_s": abs(electron_continuity_flux),
    }


def classify_root_cause(evidence: dict[str, Any]) -> dict[str, Any]:
    """Apply the ordered compensated-SG root-cause rules with explicit evidence."""
    numeric_fields = (
        "double_highprec_relative_error",
        "cancellation_condition",
        "sent_state_gap_recovery",
        "sent_state_replay_residual_dex",
        "sent_state_vector_residual_dex",
        "raw_edge_residual_dex",
        "source_residual_dex",
        "alpha_residual_dex",
        "source_closure_residual_dex",
        "terminal_residual_dex",
    )
    metrics: dict[str, float] = {}
    for name in numeric_fields:
        raw = evidence.get(name)
        if raw is None:
            continue
        try:
            value = float(raw)
        except (TypeError, ValueError) as exc:
            raise ValueError(f"classifier evidence {name} must be numeric") from exc
        if not math.isfinite(value):
            raise ValueError(f"classifier evidence {name} must be finite")
        metrics[name] = value

    clamp = bool(evidence.get("any_exponent_clamped")) or any(
        bool(evidence.get(name))
        for name in (
            "node0_exponent_clamped_low",
            "node0_exponent_clamped_high",
            "node1_exponent_clamped_low",
            "node1_exponent_clamped_high",
        )
    )
    error = metrics.get("double_highprec_relative_error")
    cancellation = metrics.get("cancellation_condition")
    if (
        error is not None
        and error > 1.0e-6
        and (clamp or (cancellation is not None and cancellation >= 1.0e12))
    ):
        classification = "variable_ni_sg_numerical_stability"
        rule = (
            "double/high-precision relative error > 1e-6 with exponent clamp "
            "or cancellation condition >= 1e12"
        )
    elif (
        metrics.get("sent_state_gap_recovery", -math.inf) >= 0.8
        and metrics.get("sent_state_replay_residual_dex", math.inf) <= 0.1
    ):
        classification = "vela_internal_state_branch"
        rule = "Sentaurus-state replay recovers >= 0.8 of gap with <= 0.1 dex residual"
    elif metrics.get("sent_state_vector_residual_dex", -math.inf) > 0.2:
        classification = "sg_discretization_ni_or_current_semantics"
        rule = "Sentaurus-state SG replay differs from vector projection by > 0.2 dex"
    elif (
        metrics.get("raw_edge_residual_dex", math.inf) <= 0.1
        and metrics.get("alpha_residual_dex", -math.inf) > 0.1
    ):
        classification = "impact_coefficient_or_source_semantics"
        rule = (
            "raw edge agrees within 0.1 dex but impact coefficient differs "
            "by > 0.1 dex"
        )
    elif (
        metrics.get("raw_edge_residual_dex", math.inf) <= 0.1
        and metrics.get("alpha_residual_dex", math.inf) <= 0.1
        and (
            metrics.get("source_residual_dex", -math.inf) > 0.2
            or metrics.get("terminal_residual_dex", -math.inf) > 0.2
        )
    ):
        classification = "ownership_support_mapping"
        rule = (
            "raw edge and alpha agree but source or terminal differs by > 0.2 dex"
        )
    elif (
        evidence.get("coarse_only") is True
        and evidence.get("main_comparison_supports_same_failure") is False
    ):
        classification = "coarse_artifact"
        rule = "later main comparison explicitly does not reproduce the coarse-only failure"
    else:
        classification = "inconclusive"
        rule = "no ordered root-cause threshold was met"

    return {
        "classification": classification,
        "rule": rule,
        "thresholds": {
            "double_highprec_relative_error": 1.0e-6,
            "cancellation_condition": 1.0e12,
            "sent_state_gap_recovery": 0.8,
            "sent_state_replay_residual_dex": 0.1,
            "sent_state_vector_residual_dex": 0.2,
            "raw_edge_residual_dex": 0.1,
            "alpha_residual_dex": 0.1,
            "source_or_terminal_residual_dex": 0.2,
        },
        "evidence": dict(evidence),
    }


def _read_sentaurus_component_field(
    path: Path,
    *,
    components: int,
) -> dict[int, float | tuple[float, float]]:
    rows = read_csv(path)
    if not rows:
        raise ValueError(f"empty Sentaurus field CSV: {path}")
    result: dict[int, float | tuple[float, float]] = {}
    required_columns = tuple(f"component{index}" for index in range(components))
    for row_index, row in enumerate(rows):
        raw_node = row.get("node_id", row.get("id"))
        if raw_node in (None, ""):
            raise ValueError(f"{path}: row {row_index} missing node id")
        node_id = int(raw_node)
        if node_id in result:
            raise ValueError(f"{path}: duplicate node id {node_id}")
        values: list[float] = []
        for column in required_columns:
            if row.get(column) in (None, ""):
                raise ValueError(f"{path}: node {node_id} missing {column}")
            value = float(row[column])
            if not math.isfinite(value):
                raise ValueError(f"{path}: node {node_id} has non-finite {column}")
            values.append(value)
        result[node_id] = values[0] if components == 1 else (values[0], values[1])
    return result


def load_sentaurus_electron_state(export_dir: Path) -> dict[str, Any]:
    """Load strict scalar endpoints plus the two-component eCurrentDensity field."""
    manifest_path = export_dir / "field_manifest.json"
    if not manifest_path.exists():
        raise FileNotFoundError(manifest_path)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8-sig"))
    vector_field = validate_manifest_vector_components(
        manifest,
        "eCurrentDensity",
        expected_components=2,
    )
    if int(vector_field.get("region", -1)) != 0:
        raise ValueError("eCurrentDensity vector must be region0")
    if str(vector_field.get("unit", "")) != "A*cm^-2":
        raise ValueError("eCurrentDensity unit must be A*cm^-2")
    if str(vector_field.get("mapping_status", "")) != "complete":
        raise ValueError("eCurrentDensity mapping_status must be complete")
    if str(vector_field.get("global_node_mapping", "")) != "global_vertex_order":
        raise ValueError("eCurrentDensity must use global_vertex_order mapping")
    for field_name, unit in SENTAURUS_SCALAR_FIELD_SPECS.items():
        _require_manifest_field(manifest, field_name, components=1, unit=unit)
    fields_dir = export_dir / "fields"
    scalar_specs = {
        "psi_V": ("ElectrostaticPotential", 1.0),
        "phin_V": ("eQuasiFermiPotential", 1.0),
        "density_m3": ("eDensity", 1.0e6),
        "mobility_m2_V_s": ("eMobility", 1.0e-4),
        "alpha_m_inv": ("eAlphaAvalanche", 1.0e2),
    }
    state: dict[str, Any] = {}
    expected_nodes: set[int] | None = None
    for output_name, (field_name, scale) in scalar_specs.items():
        raw = _read_sentaurus_component_field(
            fields_dir / f"{field_name}_region0.csv",
            components=1,
        )
        values = {node: float(value) * scale for node, value in raw.items()}
        nodes = set(values)
        if expected_nodes is None:
            expected_nodes = nodes
        elif nodes != expected_nodes:
            raise ValueError(f"{field_name} node ids do not match Sentaurus scalar state")
        state[output_name] = values
    current = _read_sentaurus_component_field(
        fields_dir / "eCurrentDensity_region0.csv",
        components=2,
    )
    if set(current) != (expected_nodes or set()):
        raise ValueError("eCurrentDensity node ids do not match Sentaurus scalar state")
    state["current_A_cm2"] = current
    return state


def _inferred_electron_ni_m3(
    density_m3: float,
    psi_V: float,
    phin_V: float,
    vt: float,
) -> float:
    if density_m3 <= 0.0:
        raise ValueError("Sentaurus electron density must be positive")
    exponent = max(-500.0, min(500.0, (psi_V - phin_V) / vt))
    value = density_m3 / math.exp(exponent)
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError("inferred Sentaurus electron ni must be positive and finite")
    return value


def _finite_abs_ratio(numerator: float, denominator: float) -> float:
    if not math.isfinite(numerator) or not math.isfinite(denominator):
        raise ValueError("ratio operands must be finite")
    if denominator == 0.0:
        return 1.0 if numerator == 0.0 else sys.float_info.max
    value = abs(numerator) / abs(denominator)
    if not math.isfinite(value):
        return sys.float_info.max
    return max(value, 1.0e-300)


def _append_replay_columns(
    output: dict[str, Any],
    prefix: str,
    replay: dict[str, Any],
) -> None:
    for name, value in replay.items():
        if isinstance(value, bool):
            output[f"{prefix}_{name}"] = int(value)
        elif isinstance(value, (int, float)):
            numeric = float(value)
            if not math.isfinite(numeric):
                raise ValueError(f"{prefix}_{name} is non-finite")
            output[f"{prefix}_{name}"] = numeric
    signed_flux = float(replay["double_flux_m2_s"])
    output[f"{prefix}_signed_flux_m2_s"] = signed_flux
    output[f"{prefix}_abs_flux_m2_s"] = abs(signed_flux)
    output[f"{prefix}_conventional_current_A_cm2"] = (
        -ELEMENTARY_CHARGE_C * signed_flux / 1.0e4
    )


def enrich_edge_with_sentaurus_replay(
    *,
    edge_row: dict[str, Any],
    sentaurus_state: dict[str, Any],
    sentaurus_node0: dict[str, Any],
    sentaurus_node1: dict[str, Any],
    temperature_K: float,
    unit_system: str = "tcad_internal",
) -> dict[str, Any]:
    """Append strict Vela SG decomposition and canonical Sentaurus replays."""
    if not math.isfinite(temperature_K) or temperature_K <= 0.0:
        raise ValueError("temperature_K must be positive and finite")
    if unit_system != "tcad_internal":
        raise ValueError("compensated SG replay requires tcad_internal unit_system")
    require_finite_fields(
        edge_row,
        (
            "edge_length_m",
            "edge_area_proxy_m2",
            "electron_mobility_m2_V_s",
            "electron_alpha_m_inv",
            "electron_source_integral",
            "x0_um", "y0_um", "x1_um", "y1_um",
        ),
        context="SG edge row",
    )
    output: dict[str, Any] = {}
    for field in ELECTRON_SG_FIELDS:
        if edge_row.get(field) in (None, ""):
            raise ValueError(f"SG edge row missing {field}")
        value = float(edge_row[field])
        if not math.isfinite(value):
            raise ValueError(f"SG edge row has non-finite {field}")
        output[field] = value

    node0 = dict(sentaurus_node0)
    node1 = dict(sentaurus_node1)
    point0 = (float(node0["x_um"]), float(node0["y_um"]))
    point1 = (float(node1["x_um"]), float(node1["y_um"]))
    sentaurus_reversed = point1 < point0
    if sentaurus_reversed:
        node0, node1 = node1, node0
        point0, point1 = point1, point0
    id0 = int(node0["id"])
    id1 = int(node1["id"])
    try:
        psi0 = float(sentaurus_state["psi_V"][id0])
        psi1 = float(sentaurus_state["psi_V"][id1])
        phin0 = float(sentaurus_state["phin_V"][id0])
        phin1 = float(sentaurus_state["phin_V"][id1])
        density0 = float(sentaurus_state["density_m3"][id0])
        density1 = float(sentaurus_state["density_m3"][id1])
        mobility0 = float(sentaurus_state["mobility_m2_V_s"][id0])
        mobility1 = float(sentaurus_state["mobility_m2_V_s"][id1])
        alpha0 = float(sentaurus_state["alpha_m_inv"][id0])
        alpha1 = float(sentaurus_state["alpha_m_inv"][id1])
        current0 = sentaurus_state["current_A_cm2"][id0]
        current1 = sentaurus_state["current_A_cm2"][id1]
    except KeyError as exc:
        raise ValueError(f"Sentaurus endpoint state missing {exc}") from exc
    endpoint_values = (
        psi0, psi1, phin0, phin1, density0, density1, mobility0, mobility1,
        alpha0, alpha1,
        *current0, *current1,
    )
    if not all(math.isfinite(float(value)) for value in endpoint_values):
        raise ValueError("Sentaurus endpoint state contains non-finite values")

    vt = BOLTZMANN_OVER_CHARGE_V_K * temperature_K
    ni0 = _inferred_electron_ni_m3(density0, psi0, phin0, vt)
    ni1 = _inferred_electron_ni_m3(density1, psi1, phin1, vt)
    projection = project_endpoint_current_to_canonical_edge(
        point0=point0,
        point1=point1,
        current0_A_cm2=current0,
        current1_A_cm2=current1,
    )
    edge_length_m = (
        float(projection["canonical_length_coordinate_units"]) * 1.0e-6
    )
    vela_mobility = float(edge_row["electron_mobility_m2_V_s"])
    sentaurus_mobility = 0.5 * (mobility0 + mobility1)
    vela_replay = replay_electron_variable_ni_sg(
        ni0=ni0,
        ni1=ni1,
        psi0=psi0,
        psi1=psi1,
        phin0=phin0,
        phin1=phin1,
        vt=vt,
        mobility_m2_V_s=vela_mobility,
        length_m=edge_length_m,
    )
    sentaurus_replay = replay_electron_variable_ni_sg(
        ni0=ni0,
        ni1=ni1,
        psi0=psi0,
        psi1=psi1,
        phin0=phin0,
        phin1=phin1,
        vt=vt,
        mobility_m2_V_s=sentaurus_mobility,
        length_m=edge_length_m,
    )
    output.update({
        "sentaurus_input_orientation_reversed": int(sentaurus_reversed),
        "sentaurus_e_psi0_V": psi0,
        "sentaurus_e_psi1_V": psi1,
        "sentaurus_e_phin0_V": phin0,
        "sentaurus_e_phin1_V": phin1,
        "sentaurus_e_density0_m3": density0,
        "sentaurus_e_density1_m3": density1,
        "sentaurus_e_mobility0_m2_V_s": mobility0,
        "sentaurus_e_mobility1_m2_V_s": mobility1,
        "sentaurus_e_ni_inferred0_m3": ni0,
        "sentaurus_e_ni_inferred1_m3": ni1,
        "sentaurus_e_alpha0_m_inv": alpha0,
        "sentaurus_e_alpha1_m_inv": alpha1,
        "sentaurus_edge_length_m": edge_length_m,
        "sentaurus_e_current_edge_signed_A_cm2": projection["conventional_current_A_cm2"],
        "sentaurus_e_current_edge_abs_A_cm2": projection["conventional_current_abs_A_cm2"],
        "sentaurus_e_continuity_edge_signed_flux_m2_s": projection["electron_continuity_flux_m2_s"],
        "sentaurus_e_continuity_edge_abs_flux_m2_s": projection["electron_continuity_flux_abs_m2_s"],
    })
    _append_replay_columns(
        output,
        "sentaurus_e_sg_vela_mobility",
        vela_replay,
    )
    _append_replay_columns(
        output,
        "sentaurus_e_sg_sentaurus_mobility",
        sentaurus_replay,
    )

    vela_point0 = (float(edge_row["x0_um"]), float(edge_row["y0_um"]))
    vela_point1 = (float(edge_row["x1_um"]), float(edge_row["y1_um"]))
    vela_orientation_sign = -1.0 if vela_point1 < vela_point0 else 1.0
    vela_production = (
        float(edge_row["electron_sg_production_signed_continuity_particle_flux_m2_s"])
        * vela_orientation_sign
    )
    sent_vector_flux = float(projection["electron_continuity_flux_m2_s"])
    vela_replay_flux = float(vela_replay["double_flux_m2_s"])
    sentaurus_replay_flux = float(sentaurus_replay["double_flux_m2_s"])
    output.update({
        "vela_edge_orientation_sign_to_canonical": vela_orientation_sign,
        "vela_e_sg_production_canonical_signed_flux_m2_s": vela_production,
        "vela_e_over_sentaurus_vector_abs_ratio": _finite_abs_ratio(
            vela_production, sent_vector_flux
        ),
        "sentaurus_e_sg_vela_mobility_over_vector_abs_ratio": _finite_abs_ratio(
            vela_replay_flux, sent_vector_flux
        ),
        "sentaurus_e_sg_sentaurus_mobility_over_vector_abs_ratio": _finite_abs_ratio(
            sentaurus_replay_flux, sent_vector_flux
        ),
    })
    source = float(edge_row["electron_source_integral"])
    alpha = float(edge_row["electron_alpha_m_inv"])
    area = float(edge_row["edge_area_proxy_m2"])
    sentaurus_alpha_average = 0.5 * (alpha0 + alpha1)
    output["sentaurus_e_alpha_edge_average_m_inv"] = sentaurus_alpha_average
    output["vela_e_over_sentaurus_alpha_abs_ratio"] = _finite_abs_ratio(
        alpha, sentaurus_alpha_average
    )
    # TCAD-native source uses cm^-1, cm^-2 s^-1 and um^2. Converting the
    # assembled 2-D source per device depth gives 100 * 1e4 * 1e-12 = 1e-6.
    source_physical = source * 1.0e-6
    output["vela_e_source_integral_physical_m_inv_s"] = source_physical
    sentaurus_source_physical = (
        sentaurus_alpha_average
        * abs(float(projection["electron_continuity_flux_m2_s"]))
        * area
    )
    output["sentaurus_e_source_on_vela_area_physical_m_inv_s"] = sentaurus_source_physical
    output["vela_e_over_sentaurus_source_abs_ratio"] = _finite_abs_ratio(
        source_physical, sentaurus_source_physical
    )
    output["vela_e_source_closure_ratio"] = _finite_abs_ratio(
        source_physical,
        alpha * abs(vela_production) * area,
    )
    require_finite_fields(
        output,
        REQUIRED_ENRICHED_FIELDS,
        context="enriched SG edge",
    )
    return output


def validate_enriched_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    validate_72_row_keys(rows)
    for index, row in enumerate(rows):
        require_finite_fields(
            row,
            REQUIRED_ENRICHED_FIELDS + CURRENT_DISCRETIZATION_RATIO_FIELDS,
            context=f"row {index}",
        )
    return rows


def classify_doping(donors: float, acceptors: float) -> tuple[str, float]:
    net = donors - acceptors
    threshold = 1.0e-6 * max(abs(donors), abs(acceptors), 1.0)
    if abs(net) <= threshold:
        return "compensated", net
    if net > 0.0:
        return "n", net
    return "p", net


def load_doping(path: Path) -> dict[int, dict[str, Any]]:
    result: dict[int, dict[str, Any]] = {}
    for row in read_csv(path):
        node_id = int(finite_float(row.get("node_id", row.get("id"))))
        donors = finite_float(row.get("donors_cm3"), 0.0)
        acceptors = finite_float(row.get("acceptors_cm3"), 0.0)
        kind, net = classify_doping(donors, acceptors)
        result[node_id] = {
            "donors_cm3": donors,
            "acceptors_cm3": acceptors,
            "net_doping_cm3": net,
            "type": kind,
        }
    return result


def parse_vtk(
    path: Path,
    *, coordinate_scale_to_um: float = 1.0e6,
) -> dict[str, Any]:
    lines = path.read_text(encoding="utf-8").splitlines()
    points: list[tuple[float, float, float]] = []
    scalars: dict[str, list[float]] = {}
    section: str | None = None
    section_count = 0
    index = 0
    while index < len(lines):
        parts = lines[index].split()
        if not parts:
            index += 1
            continue
        if parts[0] == "POINTS":
            count = int(parts[1])
            index += 1
            values: list[float] = []
            while len(values) < 3 * count:
                values.extend(float(item) for item in lines[index].split())
                index += 1
            points = [(values[i] * coordinate_scale_to_um, values[i + 1] * coordinate_scale_to_um, values[i + 2] * coordinate_scale_to_um) for i in range(0, len(values), 3)]
            continue
        if parts[0] == "POINT_DATA":
            section = "point"
            section_count = int(parts[1])
            index += 1
            continue
        if parts[0] == "CELL_DATA":
            section = "cell"
            section_count = int(parts[1])
            index += 1
            continue
        if parts[0] == "SCALARS" and section == "point":
            name = parts[1]
            index += 1
            if lines[index].strip().startswith("LOOKUP_TABLE"):
                index += 1
            values = []
            while len(values) < section_count:
                values.extend(float(item) for item in lines[index].split())
                index += 1
            scalars[name] = values[:section_count]
            continue
        if parts[0] == "VECTORS" and section in {"point", "cell"}:
            index += 1 + section_count
            continue
        index += 1
    missing = [name for name in NODE_FIELDS if name not in scalars]
    if missing:
        raise SystemExit(f"missing VTK scalar(s) {missing} in {path}")
    return {"points": points, "scalars": scalars}


def nearest_node(points: list[tuple[float, float, float]], x_um: float, y_um: float) -> int:
    best: tuple[float, int] | None = None
    for node_id, (x, y, _z) in enumerate(points):
        distance = (x - x_um) ** 2 + (y - y_um) ** 2
        if best is None or distance < best[0]:
            best = (distance, node_id)
    if best is None:
        raise SystemExit("no VTK points loaded")
    return best[1]


def vtk_for_bias(root: Path, prefix: str, bias: float) -> Path:
    index = int(round(abs(bias) / 0.05))
    exact = root / f"{prefix}_{index:04d}_{bias:g}V.vtk"
    if exact.exists():
        return exact
    matches = sorted(root.glob(f"{prefix}_{index:04d}_*.vtk"))
    if not matches:
        raise SystemExit(f"no VTK file found for bias {bias} in {root}")
    return matches[0]


def load_sg_edges(path: Path) -> dict[tuple[float, int], dict[str, str]]:
    result: dict[tuple[float, int], dict[str, str]] = {}
    for row in read_csv(path):
        bias = round(finite_float(row.get("bias_V")), 10)
        edge_id = int(finite_float(row.get("edge_id")))
        key = (bias, edge_id)
        if key in result:
            raise ValueError(
                f"duplicate SG edge row for bias={bias:g}, edge_id={edge_id}"
            )
        result[key] = row
    return result

def unique_sg_edge_for_nodes(
    edges: Mapping[tuple[float, int], dict[str, str]],
    bias: float,
    node0: int,
    node1: int,
) -> tuple[int, dict[str, str]]:
    bias_key = round(float(bias), 10)
    target_nodes = {int(node0), int(node1)}
    matches: list[tuple[int, dict[str, str]]] = []
    for (row_bias, edge_id), row in edges.items():
        if row_bias != bias_key:
            continue
        row_nodes = {
            int(finite_float(row.get("node0"))),
            int(finite_float(row.get("node1"))),
        }
        if row_nodes == target_nodes:
            matches.append((edge_id, row))
    if len(matches) != 1:
        raise ValueError(
            f"expected one SG edge for bias={bias:g}, nodes={sorted(target_nodes)}; "
            f"found {len(matches)}"
        )
    return matches[0]


def load_sentaurus_nodes(sentaurus_root: Path, bias: float) -> tuple[list[dict[str, Any]], dict[int, dict[str, Any]]]:
    bias_name = f"sentaurus_{bias:g}v"
    directory = sentaurus_root / bias_name
    if not directory.exists():
        raise SystemExit(f"missing Sentaurus export directory: {directory}")
    nodes = []
    for row in read_csv(directory / "nodes.csv"):
        nodes.append({"id": int(row["id"]), "x_um": finite_float(row["x_um"]), "y_um": finite_float(row["y_um"])})
    doping = load_doping(directory / "doping.csv")
    return nodes, doping


def nearest_sentaurus_node(nodes: list[dict[str, Any]], x_um: float, y_um: float) -> dict[str, Any]:
    best: tuple[float, dict[str, Any]] | None = None
    for node in nodes:
        distance = (node["x_um"] - x_um) ** 2 + (node["y_um"] - y_um) ** 2
        if best is None or distance < best[0]:
            best = (distance, node)
    if best is None:
        raise SystemExit("no Sentaurus nodes loaded")
    return best[1]


def row_side_nodes(state: dict[str, Any], side: str, y_um: float) -> tuple[int, int]:
    x0, x1 = (VELA_X_COLUMNS[0], VELA_X_COLUMNS[1]) if side == "left" else (VELA_X_COLUMNS[1], VELA_X_COLUMNS[2])
    points = state["points"]
    return nearest_node(points, x0, y_um), nearest_node(points, x1, y_um)


def scalar_drop(state: dict[str, Any], field: str, node0: int, node1: int) -> float:
    values = state["scalars"][field]
    return values[node1] - values[node0]


def scalar_mid(state: dict[str, Any], field: str, node0: int, node1: int) -> float:
    values = state["scalars"][field]
    return 0.5 * (values[node0] + values[node1])


def endpoint_values(state: dict[str, Any], field: str, node0: int, node1: int) -> tuple[float, float]:
    values = state["scalars"][field]
    return values[node0], values[node1]


def build_legacy_detail_rows(args: argparse.Namespace) -> list[dict[str, Any]]:
    baseline_mesh_doping = load_doping(args.baseline_report_root.parent.parent / "imported_reference" / "vela" / "doping.csv")
    probe_doping = load_doping(args.probe_root / "doping_compensated_x1_column.csv")
    variants = {
        "baseline": {
            "sg": load_sg_edges(args.baseline_report_root / "sg_avalanche_edges_density_gradient_0p05.csv"),
            "vtk_root": args.baseline_report_root / "vtk_density_gradient_0p05",
            "vtk_prefix": "dc_sweep",
            "doping": baseline_mesh_doping,
        },
        "compensated_probe": {
            "sg": load_sg_edges(args.probe_root / "sg_avalanche_edges_compensated_junction_0p05.csv"),
            "vtk_root": args.probe_root / "vtk_compensated_junction_0p05",
            "vtk_prefix": "dc_sweep",
            "doping": probe_doping,
        },
    }

    rows: list[dict[str, Any]] = []
    for bias in BIASES:
        sentaurus_nodes, sentaurus_doping = load_sentaurus_nodes(args.sentaurus_root, bias)
        states = {
            name: parse_vtk(vtk_for_bias(data["vtk_root"], data["vtk_prefix"], bias))
            for name, data in variants.items()
        }
        for y_um in Y_CUTS:
            sent_left = nearest_sentaurus_node(sentaurus_nodes, SENTAURUS_X_COLUMNS[0], y_um)
            sent_mid = nearest_sentaurus_node(sentaurus_nodes, SENTAURUS_X_COLUMNS[1], y_um)
            sent_right = nearest_sentaurus_node(sentaurus_nodes, SENTAURUS_X_COLUMNS[2], y_um)
            sent_edge_types = {
                "left": f"{sentaurus_doping[sent_left['id']]['type']}-{sentaurus_doping[sent_mid['id']]['type']}",
                "right": f"{sentaurus_doping[sent_mid['id']]['type']}-{sentaurus_doping[sent_right['id']]['type']}",
            }
            for variant_name, variant in variants.items():
                state = states[variant_name]
                for side in ["left", "right"]:
                    edge_id = EDGE_BY_SIDE[y_um][side]
                    edge_row = variant["sg"].get((round(bias, 10), edge_id))
                    if edge_row is None:
                        raise SystemExit(f"missing SG edge row for {variant_name} bias={bias} edge={edge_id}")
                    node0, node1 = row_side_nodes(state, side, y_um)
                    doping0 = variant["doping"][node0]
                    doping1 = variant["doping"][node1]
                    item: dict[str, Any] = {
                        "variant": variant_name,
                        "bias_V": bias,
                        "y_um": y_um,
                        "side": side,
                        "edge_id": edge_id,
                        "node0": node0,
                        "node1": node1,
                        "node0_type": doping0["type"],
                        "node1_type": doping1["type"],
                        "edge_type": f"{doping0['type']}-{doping1['type']}",
                        "node0_net_doping_cm3": doping0["net_doping_cm3"],
                        "node1_net_doping_cm3": doping1["net_doping_cm3"],
                        "sentaurus_edge_type": sent_edge_types[side],
                        "sentaurus_nearest_left_node": sent_left["id"],
                        "sentaurus_nearest_mid_node": sent_mid["id"],
                        "sentaurus_nearest_right_node": sent_right["id"],
                    }
                    for field, output_name in FIELD_TO_OUTPUT.items():
                        item[f"{output_name}_drop_V"] = scalar_drop(state, field, node0, node1)
                        value0, value1 = endpoint_values(state, field, node0, node1)
                        item[f"{output_name}0"] = value0
                        item[f"{output_name}1"] = value1
                    item["electron_density_mid_m3"] = scalar_mid(state, "Electrons", node0, node1)
                    item["hole_density_mid_m3"] = scalar_mid(state, "Holes", node0, node1)
                    item["electron_density_endpoint_abs_ratio"] = abs_ratio(item["electrons1"], item["electrons0"])
                    item["hole_density_endpoint_abs_ratio"] = abs_ratio(item["holes1"], item["holes0"])
                    for field in [
                        "edge_length_m",
                        "edge_couple_m",
                        "edge_area_proxy_m2",
                        "electric_field_V_per_m",
                        "electron_impact_field_V_per_m",
                        "hole_impact_field_V_per_m",
                        "electron_alpha_m_inv",
                        "hole_alpha_m_inv",
                        "electron_mobility_m2_V_s",
                        "hole_mobility_m2_V_s",
                        "electron_flux_proxy",
                        "hole_flux_proxy",
                        "electron_raw_flux_proxy",
                        "hole_raw_flux_proxy",
                        "electron_reconstructed_flux_proxy",
                        "hole_reconstructed_flux_proxy",
                        "electron_final_over_raw_flux_proxy",
                        "hole_final_over_raw_flux_proxy",
                        "electron_source_integral",
                        "hole_source_integral",
                        "edge_source_integral",
                    ]:
                        item[field] = finite_float(edge_row.get(field))
                    rows.append(item)
    add_pair_ratios(rows)
    add_probe_over_baseline(rows)
    return rows


CURRENT_DISCRETIZATION_RATIO_FIELDS = (
    "edge_source_integral",
    "electron_source_integral",
    "hole_source_integral",
    "electron_alpha_m_inv",
    "hole_alpha_m_inv",
    "electron_flux_proxy",
    "hole_flux_proxy",
    "electron_mobility_m2_V_s",
    "hole_mobility_m2_V_s",
)
STANDARD_REPLAY_RATIO_FIELDS = (
    "vela_e_sg_production_canonical_signed_flux_m2_s",
    "sentaurus_e_continuity_edge_signed_flux_m2_s",
    "sentaurus_e_sg_vela_mobility_signed_flux_m2_s",
    "sentaurus_e_sg_sentaurus_mobility_signed_flux_m2_s",
    "sentaurus_e_alpha_edge_average_m_inv",
    "vela_e_over_sentaurus_alpha_abs_ratio",
    "vela_e_source_integral_physical_m_inv_s",
    "sentaurus_e_source_on_vela_area_physical_m_inv_s",
    "vela_e_over_sentaurus_source_abs_ratio",
    "vela_e_source_closure_ratio",
) + CURRENT_DISCRETIZATION_RATIO_FIELDS


def _log_gap_from_abs_ratio(ratio: float) -> float:
    value = float(ratio)
    if not math.isfinite(value) or value <= 0.0:
        raise ValueError("SG replay ratio must be positive and finite")
    return abs(math.log10(value))


def _standard_row_classification(enriched: dict[str, Any]) -> dict[str, Any]:
    raw_gap = _log_gap_from_abs_ratio(
        enriched["vela_e_over_sentaurus_vector_abs_ratio"]
    )
    replay_gap = _log_gap_from_abs_ratio(
        enriched["sentaurus_e_sg_vela_mobility_over_vector_abs_ratio"]
    )
    recovery = (
        1.0 - replay_gap / raw_gap
        if raw_gap > 1.0e-15
        else (1.0 if replay_gap <= 1.0e-15 else 0.0)
    )
    evidence = {
        "double_highprec_relative_error": float(
            enriched[
                "electron_sg_production_vs_high_precision_reference_relative_error"
            ]
        ),
        "cancellation_condition": float(
            enriched["electron_sg_cancellation_condition"]
        ),
        "any_exponent_clamped": any(
            bool(float(enriched[name]))
            for name in (
                "electron_sg_node0_exponent_clamped_low",
                "electron_sg_node0_exponent_clamped_high",
                "electron_sg_node1_exponent_clamped_low",
                "electron_sg_node1_exponent_clamped_high",
            )
        ),
        "raw_edge_residual_dex": raw_gap,
        "sent_state_vector_residual_dex": replay_gap,
        "sent_state_replay_residual_dex": replay_gap,
        "sent_state_gap_recovery": recovery,
        "source_residual_dex": _log_gap_from_abs_ratio(
            enriched["vela_e_over_sentaurus_source_abs_ratio"]
        ),
        "sentaurus_source_residual_dex": _log_gap_from_abs_ratio(
            enriched["vela_e_over_sentaurus_source_abs_ratio"]
        ),
        "alpha_residual_dex": _log_gap_from_abs_ratio(
            enriched["vela_e_over_sentaurus_alpha_abs_ratio"]
        ),
        "source_closure_residual_dex": _log_gap_from_abs_ratio(
            enriched["vela_e_source_closure_ratio"]
        ),
    }
    result = classify_root_cause(evidence)
    return {
        "root_cause_classification": result["classification"],
        "root_cause_rule": result["rule"],
        "classifier_gap_recovery": recovery,
        "classifier_raw_gap_dex": raw_gap,
        "classifier_sent_state_residual_dex": replay_gap,
        "classifier_source_residual_dex": evidence["source_residual_dex"],
        "classifier_alpha_residual_dex": evidence["alpha_residual_dex"],
        "classifier_source_closure_residual_dex": evidence["source_closure_residual_dex"],
        "root_cause_evidence_json": json.dumps(
            clean_json(evidence), sort_keys=True, separators=(",", ":")
        ),
        **evidence,
    }


def _append_standard_ratios(rows: list[dict[str, Any]]) -> None:
    by_pair: dict[tuple[str, float, float], dict[str, dict[str, Any]]] = {}
    for row in rows:
        metadata = REPLAY_VARIANT_MATRIX[str(row["variant"])]
        row.setdefault("doping_strategy", metadata["doping_strategy"])
        row.setdefault("current_variant", metadata["current_variant"])
        by_pair.setdefault(
            (row["variant"], row["bias_V"], row["y_um"]), {}
        )[row["side"]] = row
    for pair in by_pair.values():
        left = pair["left"]
        right = pair["right"]
        for field in STANDARD_REPLAY_RATIO_FIELDS:
            ratio = abs_ratio(float(right[field]), float(left[field]))
            left[f"right_over_left_{field}"] = ratio
            right[f"right_over_left_{field}"] = ratio

    by_current_pair: dict[
        tuple[str, float, float, str], dict[str, dict[str, Any]]
    ] = {}
    for row in rows:
        by_current_pair.setdefault(
            (
                row["doping_strategy"],
                row["bias_V"],
                row["y_um"],
                row["side"],
            ),
            {},
        )[row["current_variant"]] = row
    for pair in by_current_pair.values():
        density = pair["density_gradient"]
        gss = pair["gss_midpoint"]
        for field in CURRENT_DISCRETIZATION_RATIO_FIELDS:
            ratio = abs_ratio(float(gss[field]), float(density[field]))
            if ratio is None or not math.isfinite(ratio):
                raise ValueError(f"non-finite GSS/density ratio for {field}")
            key = f"gss_midpoint_over_density_gradient_{field}"
            density[key] = ratio
            gss[key] = ratio

    by_doping_pair: dict[
        tuple[str, float, float, str], dict[str, dict[str, Any]]
    ] = {}
    for row in rows:
        by_doping_pair.setdefault(
            (
                row["current_variant"],
                row["bias_V"],
                row["y_um"],
                row["side"],
            ),
            {},
        )[row["doping_strategy"]] = row
    for pair in by_doping_pair.values():
        if "legacy" not in pair or "reported" not in pair:
            continue
        legacy = pair["legacy"]
        reported = pair["reported"]
        for field in STANDARD_REPLAY_RATIO_FIELDS:
            ratio = abs_ratio(float(reported[field]), float(legacy[field]))
            legacy[f"reported_over_legacy_{field}"] = ratio
            reported[f"reported_over_legacy_{field}"] = ratio


def current_discretization_pair_rows(
    rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    pairs: dict[tuple[str, float, float, str], dict[str, Any]] = {}
    for row in rows:
        key = (
            str(row["doping_strategy"]),
            float(row["bias_V"]),
            float(row["y_um"]),
            str(row["side"]),
        )
        if str(row["current_variant"]) != "gss_midpoint":
            continue
        item = {
            "doping_strategy": key[0],
            "bias_V": key[1],
            "y_um": key[2],
            "side": key[3],
        }
        for field in CURRENT_DISCRETIZATION_RATIO_FIELDS:
            ratio_key = f"gss_midpoint_over_density_gradient_{field}"
            item[ratio_key] = float(row[ratio_key])
        pairs[key] = item
    return [pairs[key] for key in sorted(pairs)]


def build_standard_detail_rows(args: argparse.Namespace) -> list[dict[str, Any]]:
    specs = standard_variant_inputs(args.variants_root)
    variants: dict[str, dict[str, Any]] = {}
    for name, spec in specs.items():
        variants[name] = {
            "sg": load_sg_edges(spec["sg_csv"]),
            "vtk_root": spec["vtk_root"],
            "vtk_prefix": "dc_sweep",
            "doping": load_doping(spec["doping_csv"]),
            "metadata": spec,
        }

    rows: list[dict[str, Any]] = []
    for bias in BIASES:
        export_dir = args.sentaurus_root / f"sentaurus_{bias:g}v"
        sentaurus_nodes, sentaurus_doping = load_sentaurus_nodes(
            args.sentaurus_root, bias
        )
        sentaurus_state = load_sentaurus_electron_state(export_dir)
        states = {
            name: parse_vtk(
                vtk_for_bias(data["vtk_root"], data["vtk_prefix"], bias),
                coordinate_scale_to_um=1.0,
            )
            for name, data in variants.items()
        }
        for y_um in Y_CUTS:
            sent_left = nearest_sentaurus_node(
                sentaurus_nodes, SENTAURUS_X_COLUMNS[0], y_um
            )
            sent_mid = nearest_sentaurus_node(
                sentaurus_nodes, SENTAURUS_X_COLUMNS[1], y_um
            )
            sent_right = nearest_sentaurus_node(
                sentaurus_nodes, SENTAURUS_X_COLUMNS[2], y_um
            )
            sent_pairs = {
                "left": (sent_left, sent_mid),
                "right": (sent_mid, sent_right),
            }
            for variant_name, variant in variants.items():
                state = states[variant_name]
                for side in ("left", "right"):
                    node0, node1 = row_side_nodes(state, side, y_um)
                    try:
                        edge_id, edge_row = unique_sg_edge_for_nodes(
                            variant["sg"], bias, node0, node1
                        )
                    except ValueError as exc:
                        raise ValueError(f"{variant_name}: {exc}") from exc
                    doping0 = variant["doping"][node0]
                    doping1 = variant["doping"][node1]
                    sent0, sent1 = sent_pairs[side]
                    item: dict[str, Any] = {
                        "variant": variant_name,
                        "doping_strategy": variant["metadata"]["doping_strategy"],
                        "compensated_doping_policy": variant["metadata"]["compensated_doping_policy"],
                        "current_variant": variant["metadata"]["current_variant"],
                        "current_approximation": variant["metadata"]["current_approximation"],
                        "bias_V": bias,
                        "y_um": y_um,
                        "side": side,
                        "edge_id": edge_id,
                        "node0": node0,
                        "node1": node1,
                        "node0_type": doping0["type"],
                        "node1_type": doping1["type"],
                        "edge_type": f"{doping0['type']}-{doping1['type']}",
                        "node0_net_doping_cm3": doping0["net_doping_cm3"],
                        "node1_net_doping_cm3": doping1["net_doping_cm3"],
                        "sentaurus_edge_type": (
                            f"{sentaurus_doping[sent0['id']]['type']}-"
                            f"{sentaurus_doping[sent1['id']]['type']}"
                        ),
                        "sentaurus_nearest_node0": sent0["id"],
                        "sentaurus_nearest_node1": sent1["id"],
                    }
                    for field, output_name in FIELD_TO_OUTPUT.items():
                        item[f"{output_name}_drop_V"] = scalar_drop(
                            state, field, node0, node1
                        )
                        value0, value1 = endpoint_values(
                            state, field, node0, node1
                        )
                        item[f"{output_name}0"] = value0
                        item[f"{output_name}1"] = value1
                    for key, raw in edge_row.items():
                        if key in item:
                            continue
                        if key == "edge_class":
                            item[key] = str(raw)
                            continue
                        value = finite_float(raw)
                        if not math.isfinite(value):
                            raise ValueError(
                                f"non-finite SG edge field {key} for edge {edge_id}"
                            )
                        item[key] = value
                    enriched = enrich_edge_with_sentaurus_replay(
                        edge_row=edge_row,
                        sentaurus_state=sentaurus_state,
                        sentaurus_node0=sent0,
                        sentaurus_node1=sent1,
                        temperature_K=300.0,
                        unit_system="tcad_internal",
                    )
                    item.update(enriched)
                    item.update(_standard_row_classification(enriched))
                    rows.append(item)
    validate_enriched_rows(rows)
    _append_standard_ratios(rows)
    return rows


def build_detail_rows(args: argparse.Namespace) -> list[dict[str, Any]]:
    if args.variants_root is not None:
        return build_standard_detail_rows(args)
    return build_legacy_detail_rows(args)


def add_pair_ratios(rows: list[dict[str, Any]]) -> None:
    by_key: dict[tuple[str, float, float], dict[str, dict[str, Any]]] = {}
    for row in rows:
        by_key.setdefault((row["variant"], row["bias_V"], row["y_um"]), {})[row["side"]] = row
    for (_variant, _bias, _y), pair in by_key.items():
        left = pair.get("left")
        right = pair.get("right")
        if left is None or right is None:
            continue
        for field in RATIO_FIELDS + ["psi_drop_V", "phin_drop_V", "phip_drop_V"]:
            ratio = abs_ratio(right.get(field, math.nan), left.get(field, math.nan))
            left[f"right_over_left_{field}"] = ratio
            right[f"right_over_left_{field}"] = ratio


def add_probe_over_baseline(rows: list[dict[str, Any]]) -> None:
    by_key: dict[tuple[float, float, str], dict[str, dict[str, Any]]] = {}
    for row in rows:
        by_key.setdefault((row["bias_V"], row["y_um"], row["side"]), {})[row["variant"]] = row
    for (_bias, _y, _side), pair in by_key.items():
        baseline = pair.get("baseline")
        probe = pair.get("compensated_probe")
        if baseline is None or probe is None:
            continue
        for field in RATIO_FIELDS + ["psi_drop_V", "phin_drop_V", "phip_drop_V"]:
            ratio = abs_ratio(probe.get(field, math.nan), baseline.get(field, math.nan))
            baseline[f"probe_over_baseline_{field}"] = ratio
            probe[f"probe_over_baseline_{field}"] = ratio


def median(values: list[float]) -> float | None:
    finite = [value for value in values if value is not None and math.isfinite(value)]
    return statistics.median(finite) if finite else None


def summarize(rows: list[dict[str, Any]]) -> dict[str, Any]:
    aggregate: list[dict[str, Any]] = []
    dominant_by_bias: list[dict[str, Any]] = []
    for variant in ["baseline", "compensated_probe"]:
        for bias in BIASES:
            subset = [row for row in rows if row["variant"] == variant and row["bias_V"] == bias and row["side"] == "right"]
            item: dict[str, Any] = {"variant": variant, "bias_V": bias}
            for field in [
                "edge_source_integral",
                "electron_source_integral",
                "hole_source_integral",
                "phin_drop_V",
                "electron_alpha_m_inv",
                "hole_alpha_m_inv",
                "electron_flux_proxy",
                "hole_flux_proxy",
                "electron_raw_flux_proxy",
                "hole_raw_flux_proxy",
                "electron_mobility_m2_V_s",
                "hole_mobility_m2_V_s",
                "electron_density_mid_m3",
                "hole_density_mid_m3",
            ]:
                item[f"median_right_over_left_{field}"] = median([
                    row.get(f"right_over_left_{field}") for row in subset
                ])
            aggregate.append(item)

    for bias in BIASES:
        subset = [row for row in rows if row["variant"] == "compensated_probe" and row["bias_V"] == bias and row["side"] == "right"]
        source_ratio = median([row.get("right_over_left_edge_source_integral") for row in subset])
        source_log = log10_abs(source_ratio)
        channel = {
            "electron_source_right_left_ratio": median([row.get("right_over_left_electron_source_integral") for row in subset]),
            "hole_source_right_left_ratio": median([row.get("right_over_left_hole_source_integral") for row in subset]),
            "electron_alpha_right_left_ratio": median([row.get("right_over_left_electron_alpha_m_inv") for row in subset]),
            "electron_flux_right_left_ratio": median([row.get("right_over_left_electron_flux_proxy") for row in subset]),
            "electron_raw_flux_right_left_ratio": median([row.get("right_over_left_electron_raw_flux_proxy") for row in subset]),
            "electron_mobility_right_left_ratio": median([row.get("right_over_left_electron_mobility_m2_V_s") for row in subset]),
            "electron_density_mid_right_left_ratio": median([row.get("right_over_left_electron_density_mid_m3") for row in subset]),
            "hole_alpha_right_left_ratio": median([row.get("right_over_left_hole_alpha_m_inv") for row in subset]),
            "hole_flux_right_left_ratio": median([row.get("right_over_left_hole_flux_proxy") for row in subset]),
            "edge_area_right_left_ratio": median([row.get("right_over_left_edge_area_proxy_m2") for row in subset]),
        }
        channel["electron_alpha_x_flux_right_left_ratio"] = (
            channel["electron_alpha_right_left_ratio"] * channel["electron_flux_right_left_ratio"]
            if channel["electron_alpha_right_left_ratio"] is not None and channel["electron_flux_right_left_ratio"] is not None
            else None
        )
        channel["hole_alpha_x_flux_right_left_ratio"] = (
            channel["hole_alpha_right_left_ratio"] * channel["hole_flux_right_left_ratio"]
            if channel["hole_alpha_right_left_ratio"] is not None and channel["hole_flux_right_left_ratio"] is not None
            else None
        )
        if (channel["electron_source_right_left_ratio"] or 0.0) > 1.0 and (channel["hole_source_right_left_ratio"] or math.inf) < 1.0:
            channel["dominant_physical_reading"] = "electron source is right-heavy while hole source is left-heavy; residual right bias follows electron SG flux proxy moderated by alpha"
        elif (channel["hole_source_right_left_ratio"] or 0.0) > 1.0 and (channel["electron_source_right_left_ratio"] or math.inf) < 1.0:
            channel["dominant_physical_reading"] = "hole source is right-heavy while electron source is left-heavy"
        else:
            channel["dominant_physical_reading"] = "both carrier source channels have the same right/left direction or one channel is unavailable"

        candidates = []
        for label, fields in [
            ("electron_flux_proxy", ["right_over_left_electron_flux_proxy"]),
            ("electron_raw_flux_proxy", ["right_over_left_electron_raw_flux_proxy"]),
            ("electron_alpha", ["right_over_left_electron_alpha_m_inv"]),
            ("electron_density_mid", ["right_over_left_electron_density_mid_m3"]),
            ("electron_mobility", ["right_over_left_electron_mobility_m2_V_s"]),
            ("hole_flux_proxy", ["right_over_left_hole_flux_proxy"]),
            ("hole_alpha", ["right_over_left_hole_alpha_m_inv"]),
            ("edge_area", ["right_over_left_edge_area_proxy_m2"]),
            ("electron_alpha_x_flux", ["right_over_left_electron_alpha_m_inv", "right_over_left_electron_flux_proxy"]),
            ("hole_alpha_x_flux", ["right_over_left_hole_alpha_m_inv", "right_over_left_hole_flux_proxy"]),
        ]:
            product = 1.0
            ok = True
            for field in fields:
                value = median([row.get(field) for row in subset])
                if value is None or not math.isfinite(value):
                    ok = False
                    break
                product *= value
            if not ok:
                continue
            candidate_log = log10_abs(product)
            distance = abs(candidate_log - source_log) if candidate_log is not None and source_log is not None else None
            candidates.append({
                "factor": label,
                "median_ratio_or_product": product,
                "log10_ratio_or_product": candidate_log,
                "distance_to_source_log10": distance,
            })
        candidates.sort(key=lambda item: math.inf if item["distance_to_source_log10"] is None else item["distance_to_source_log10"])
        dominant_by_bias.append({
            "bias_V": bias,
            "source_right_left_ratio": source_ratio,
            "source_log10_ratio": source_log,
            "channel_decomposition": channel,
            "closest_factor": candidates[0] if candidates else None,
            "ranked_factors": candidates,
        })
    return {"aggregate": aggregate, "dominant_by_bias": dominant_by_bias}


def summarize_standard(rows: list[dict[str, Any]]) -> dict[str, Any]:
    aggregate: list[dict[str, Any]] = []
    classifications: list[dict[str, Any]] = []
    classification_counts: dict[str, int] = {}
    evidence_fields = (
        "double_highprec_relative_error",
        "cancellation_condition",
        "sent_state_gap_recovery",
        "sent_state_replay_residual_dex",
        "sent_state_vector_residual_dex",
        "raw_edge_residual_dex",
        "source_residual_dex",
        "sentaurus_source_residual_dex",
        "alpha_residual_dex",
        "source_closure_residual_dex",
    )
    for row in rows:
        name = str(row["root_cause_classification"])
        classification_counts[name] = classification_counts.get(name, 0) + 1
    for variant in REPLAY_VARIANTS:
        for bias in BIASES:
            subset = [
                row for row in rows
                if row["variant"] == variant and row["bias_V"] == bias
            ]
            metadata = REPLAY_VARIANT_MATRIX[variant]
            item: dict[str, Any] = {
                "variant": variant,
                "doping_strategy": metadata["doping_strategy"],
                "current_variant": metadata["current_variant"],
                "current_approximation": metadata["current_approximation"],
                "bias_V": bias,
                "edge_count": len(subset),
            }
            aggregate_evidence: dict[str, Any] = {}
            for field in evidence_fields:
                value = median([
                    float(row.get(field, 0.0)) for row in subset
                ])
                item[f"median_{field}"] = value
                if value is not None:
                    aggregate_evidence[field] = value
            aggregate_evidence["any_exponent_clamped"] = any(
                bool(float(row.get(name, 0.0)))
                for row in subset
                for name in (
                    "electron_sg_node0_exponent_clamped_low",
                    "electron_sg_node0_exponent_clamped_high",
                    "electron_sg_node1_exponent_clamped_low",
                    "electron_sg_node1_exponent_clamped_high",
                )
            )
            classification = classify_root_cause(aggregate_evidence)
            item["classification"] = classification["classification"]
            item["classification_rule"] = classification["rule"]
            for field in CURRENT_DISCRETIZATION_RATIO_FIELDS:
                item[f"median_right_over_left_{field}"] = median([
                    float(row[f"right_over_left_{field}"]) for row in subset
                ])
            aggregate.append(item)
            classifications.append({
                "variant": variant,
                "bias_V": bias,
                "classification": classification["classification"],
                "rule": classification["rule"],
                "evidence": aggregate_evidence,
            })
    non_inconclusive = {
        key: value for key, value in classification_counts.items()
        if key != "inconclusive"
    }
    dominant = (
        max(non_inconclusive, key=lambda key: (non_inconclusive[key], key))
        if non_inconclusive else "inconclusive"
    )
    return {
        "schema": "vela.pn2d_bv_compensated_sg_replay.summary.v2",
        "row_count": len(rows),
        "variants": {
            name: {
                "implementation": "current_head",
                **metadata,
            }
            for name, metadata in REPLAY_VARIANT_MATRIX.items()
        },
        "classification_counts": classification_counts,
        "dominant_classification": dominant,
        "aggregate": aggregate,
        "classifications": classifications,
    }


def write_standard_report(path: Path, summary: dict[str, Any]) -> None:
    lines = [
        "# PN2D BV Compensated SG Same-Edge Replay",
        "",
        "The matrix compares legacy/reported compensated-doping policies with "
        "density-gradient SG current and the GSS Bernoulli-midpoint current proxy.",
        "",
        f"- Rows: `{summary['row_count']}`",
        f"- Dominant coarse classification: `{summary['dominant_classification']}`",
        "",
        "## Structured Root-Cause Classifications",
        "",
        "| variant | bias (V) | edges | raw gap (dex) | Sent-state residual (dex) | gap recovery | classification |",
        "|---|---:|---:|---:|---:|---:|---|",
    ]
    for item in summary["aggregate"]:
        lines.append(
            "| {variant} | {bias:g} | {count} | {raw:.6g} | {residual:.6g} | "
            "{recovery:.6g} | {classification} |".format(
                variant=item["variant"],
                bias=item["bias_V"],
                count=item["edge_count"],
                raw=item.get("median_raw_edge_residual_dex") or 0.0,
                residual=item.get("median_sent_state_replay_residual_dex") or 0.0,
                recovery=item.get("median_sent_state_gap_recovery") or 0.0,
                classification=item["classification"],
            )
        )
    lines.append("")
    run_statuses = summary.get("run_statuses", [])
    if run_statuses:
        lines.extend([
            "## Run Status",
            "",
            "| variant | status | last converged bias (V) | handoff stage |",
            "|---|---|---:|---|",
        ])
        for status in run_statuses:
            lines.append(
                "| {variant} | {run_status} | {last_bias} | {handoff} |".format(
                    variant=status["variant"],
                    run_status=status["run_status"],
                    last_bias=status["last_converged_bias_V"],
                    handoff=status["handoff_stage"],
                )
            )
        lines.append("")
    current_pairs = summary.get("current_discretization_pairs", [])
    if current_pairs:
        lines.extend([
            "## GSS Midpoint Over Density Gradient",
            "",
            "| doping | bias (V) | y (um) | side | total source | electron source | hole source | e-alpha | e-flux | e-mobility |",
            "|---|---:|---:|---|---:|---:|---:|---:|---:|---:|",
        ])
        for pair in current_pairs:
            lines.append(
                "| {doping_strategy} | {bias_V:g} | {y_um:g} | {side} | "
                "{total:.6g} | {electron:.6g} | {hole:.6g} | {alpha:.6g} | "
                "{flux:.6g} | {mobility:.6g} |".format(
                    **pair,
                    total=pair["gss_midpoint_over_density_gradient_edge_source_integral"],
                    electron=pair["gss_midpoint_over_density_gradient_electron_source_integral"],
                    hole=pair["gss_midpoint_over_density_gradient_hole_source_integral"],
                    alpha=pair["gss_midpoint_over_density_gradient_electron_alpha_m_inv"],
                    flux=pair["gss_midpoint_over_density_gradient_electron_flux_proxy"],
                    mobility=pair["gss_midpoint_over_density_gradient_electron_mobility_m2_V_s"],
                )
            )
        lines.append("")
    lines.extend([
        "",
        "`sentaurus_e_source_on_vela_area_physical_m_inv_s` is a same-area "
        "proxy (endpoint alpha arithmetic mean times projected vector current "
        "on the Vela edge area), not a native Sentaurus source discretization.",
        "",
        "This report is diagnostic only. A coarse classification is not promoted "
        "to a solver fix until the main-mesh five-anchor gate is evaluated.",
    ])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_report(path: Path, summary: dict[str, Any], rows: list[dict[str, Any]]) -> None:
    lines: list[str] = []
    lines.append("# PN2D BV Compensated Junction Source Proxy Compare")
    lines.append("")
    lines.append("## Summary")
    lines.append("")
    lines.append("This diagnostic compares the baseline density-gradient BV run with the compensated-junction probe at -12 V, -19 V, and -20 V. It classifies nodes from donors-acceptors and decomposes the remaining right-heavy source into QF drops, carrier densities, mobilities, SG flux proxies, alpha, and source integrals.")
    lines.append("")
    lines.append("Direct `phin/phip` clamp or zeroing is intentionally not part of this diagnostic. Doping classification is used only as metadata for artifact alignment and source-proxy interpretation.")
    lines.append("")
    lines.append("## Median Right/Left Ratios")
    lines.append("")
    lines.append("| variant | bias | source | phin drop | e-alpha | e-flux proxy | e-raw flux | e-density mid | e-mobility |")
    lines.append("|---|---:|---:|---:|---:|---:|---:|---:|---:|")
    for item in summary["aggregate"]:
        lines.append("| {variant} | {bias:g} | {source:.6g} | {phin:.6g} | {alpha:.6g} | {flux:.6g} | {raw:.6g} | {density:.6g} | {mob:.6g} |".format(
            variant=item["variant"],
            bias=item["bias_V"],
            source=item.get("median_right_over_left_edge_source_integral") or math.nan,
            phin=item.get("median_right_over_left_phin_drop_V") or math.nan,
            alpha=item.get("median_right_over_left_electron_alpha_m_inv") or math.nan,
            flux=item.get("median_right_over_left_electron_flux_proxy") or math.nan,
            raw=item.get("median_right_over_left_electron_raw_flux_proxy") or math.nan,
            density=item.get("median_right_over_left_electron_density_mid_m3") or math.nan,
            mob=item.get("median_right_over_left_electron_mobility_m2_V_s") or math.nan,
        ))
    lines.append("")
    lines.append("## Channel Source Decomposition For Compensated Probe")
    lines.append("")
    lines.append("| bias | total source R/L | electron source R/L | hole source R/L | e-alpha | e-flux proxy | e-alpha x flux | e-mobility | reading |")
    lines.append("|---:|---:|---:|---:|---:|---:|---:|---:|---|")
    for item in summary["dominant_by_bias"]:
        channel = item.get("channel_decomposition") or {}
        lines.append("| {bias:g} | {source:.6g} | {esrc:.6g} | {hsrc:.6g} | {ealpha:.6g} | {eflux:.6g} | {eaf:.6g} | {emob:.6g} | {reading} |".format(
            bias=item["bias_V"],
            source=item.get("source_right_left_ratio") or math.nan,
            esrc=channel.get("electron_source_right_left_ratio") or math.nan,
            hsrc=channel.get("hole_source_right_left_ratio") or math.nan,
            ealpha=channel.get("electron_alpha_right_left_ratio") or math.nan,
            eflux=channel.get("electron_flux_right_left_ratio") or math.nan,
            eaf=channel.get("electron_alpha_x_flux_right_left_ratio") or math.nan,
            emob=channel.get("electron_mobility_right_left_ratio") or math.nan,
            reading=channel.get("dominant_physical_reading", ""),
        ))
    lines.append("")
    lines.append("## Scalar Closest-Factor Ranking For Compensated Probe")
    lines.append("")
    lines.append("This table is a scalar log-distance screen only. Use it with the channel decomposition above so that a numerically close factor is not mistaken for the contributing carrier channel.")
    lines.append("")
    lines.append("| bias | source right/left | closest factor | factor ratio/product | log10 distance |")
    lines.append("|---:|---:|---|---:|---:|")
    for item in summary["dominant_by_bias"]:
        closest = item.get("closest_factor") or {}
        lines.append("| {bias:g} | {source:.6g} | {factor} | {ratio:.6g} | {distance:.6g} |".format(
            bias=item["bias_V"],
            source=item.get("source_right_left_ratio") or math.nan,
            factor=closest.get("factor", ""),
            ratio=closest.get("median_ratio_or_product") or math.nan,
            distance=closest.get("distance_to_source_log10") or math.nan,
        ))
    lines.append("")
    lines.append("## Interpretation")
    lines.append("")
    probe_agg = [item for item in summary["aggregate"] if item["variant"] == "compensated_probe"]
    max_phin_ratio = max((item.get("median_right_over_left_phin_drop_V") or 0.0) for item in probe_agg)
    source_ratios = [item.get("median_right_over_left_edge_source_integral") for item in probe_agg]
    lines.append(f"- The compensated probe keeps left/right `phin` drops balanced: max median right/left ratio is `{max_phin_ratio:.6g}`.")
    lines.append("- Remaining source right/left ratios are `{}` for -12/-19/-20 V.".format(
        ", ".join(f"{value:.6g}" for value in source_ratios if value is not None)
    ))
    lines.append("- Channel decomposition shows the residual right-heavy source is carried by the electron source channel; the hole source channel is left-heavy at all three inspected biases.")
    lines.append("- The electron right/left source ratio is driven mainly by the electron SG flux proxy / raw flux proxy, while electron alpha is below 1 and mobility is close to 1 after compensation.")
    lines.append("- This points the next debug target at density-gradient SG current/source construction and carrier-density/flux proxy selection, not at a QF hard limiter.")
    lines.append("")
    lines.append("## Artifacts")
    lines.append("")
    lines.append("- `compensated_source_proxy_compare.csv`")
    lines.append("- `compensated_source_proxy_compare_summary.json`")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    rows = build_detail_rows(args)
    if args.variants_root is not None:
        validate_enriched_rows(rows)
        _append_standard_ratios(rows)
        current_pairs = current_discretization_pair_rows(rows)
        if len(current_pairs) != 36:
            raise ValueError(
                f"expected 36 GSS/density paired rows, got {len(current_pairs)}"
            )
        run_statuses = [
            variant_run_status(spec)
            for spec in standard_variant_inputs(args.variants_root).values()
        ]
        summary = summarize_standard(rows)
        summary["current_discretization_pairs"] = current_pairs
        summary["run_statuses"] = run_statuses
        csv_path = args.out_dir / "compensated_sg_replay.csv"
        json_path = args.out_dir / "compensated_sg_replay.json"
        report_path = args.out_dir / "compensated_sg_replay_report.md"
        write_csv(csv_path, rows)
        json_path.write_text(json.dumps(clean_json({
            "schema": "vela.pn2d_bv_compensated_sg_replay.v2",
            "row_count": len(rows),
            "summary": summary,
            "classifications": summary["classifications"],
            "current_discretization_pairs": current_pairs,
            "run_statuses": run_statuses,
            "rows": rows,
        }), indent=2) + "\n", encoding="utf-8")
        write_standard_report(report_path, summary)
        print(json.dumps({"csv": str(csv_path), "json": str(json_path), "report": str(report_path), "rows": len(rows)}, indent=2))
        return
    summary = summarize(rows)
    csv_path = args.out_dir / "compensated_source_proxy_compare.csv"
    json_path = args.out_dir / "compensated_source_proxy_compare_summary.json"
    report_path = args.out_dir / "compensated_source_proxy_compare_report_20260709.md"
    write_csv(csv_path, rows)
    json_path.write_text(json.dumps(clean_json(summary), indent=2), encoding="utf-8")
    write_report(report_path, summary, rows)
    print(json.dumps({
        "csv": str(csv_path),
        "json": str(json_path),
        "report": str(report_path),
        "rows": len(rows),
    }, indent=2))


if __name__ == "__main__":
    main()
