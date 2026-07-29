#!/usr/bin/env python3
"""Fail-closed PN2D avalanche-on curve and knee parity verifier."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


EXACT_BIAS_TOLERANCE_V = 1.0e-10
PLOTTING_FLOOR_A_PER_UM = 1.0e-30
GLOBAL_BIASES_V = tuple(float(-value) for value in range(21))
KNEE_BIASES_V = (
    -18.0,
    -18.5,
    -19.0,
    -19.25,
    -19.5,
    -19.7,
    -19.8,
    -19.85,
    -19.9,
    -19.95,
    -20.0,
)
CURVE_NAMES = (
    "vela_on",
    "vela_off",
    "sentaurus_on",
    "sentaurus_off",
)
CURRENT_COLUMNS = (
    "current_A_per_um",
    "current_total_A_per_um",
    "terminal_current_A_per_um",
    "vela_current_total_A_per_um",
    "sentaurus_anode_total_current_A",
    "Anode TotalCurrent",
    "current_total",
)


class CurveContractError(ValueError):
    """A typed curve-input or metric-contract failure."""

    def __init__(self, reason: str, detail: str):
        self.reason = reason
        self.detail = detail
        super().__init__(f"{reason}: {detail}")


@dataclass(frozen=True)
class CurvePoint:
    bias_V: float
    current_A_per_um: float
    converged: bool = True
    failure_reason: str = ""
    source_index: int = 0
    electron_closure_relative: float | None = None
    hole_closure_relative: float | None = None
    terminal_pair_closure_A_per_um: float | None = None
    internal_kcl_relative: float | None = None


def _finite(value: Any, label: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise CurveContractError("non_numeric_value", label) from error
    if not math.isfinite(result):
        raise CurveContractError("nonfinite_value", label)
    return result


def _optional_finite(row: Mapping[str, str], name: str) -> float | None:
    value = row.get(name, "")
    return None if value == "" else _finite(value, name)


def _as_bool(value: str | None) -> bool:
    if value is None or value == "":
        return True
    normalized = value.strip().lower()
    if normalized in {"1", "true", "yes", "converged", "accepted", "passed"}:
        return True
    if normalized in {"0", "false", "no", "failed", "rejected"}:
        return False
    raise CurveContractError("invalid_convergence_flag", value)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_curve(path: Path) -> list[CurvePoint]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise CurveContractError("empty_curve", str(path))
    fieldnames = set(rows[0])
    if {"carrier", "quantity", "value"}.issubset(fieldnames):
        rows = [
            row
            for row in rows
            if row["carrier"] == "total"
            and row["quantity"] == "terminal_current"
            and row.get("provenance", "native") == "native"
        ]
        if not rows:
            raise CurveContractError("missing_terminal_current_rows", str(path))
        fieldnames = set(rows[0])
    bias_column = next(
        (name for name in ("bias_V", "actual_bias_V") if name in fieldnames),
        None,
    )
    current_column = (
        "value"
        if {"carrier", "quantity", "value"}.issubset(fieldnames)
        else next(
            (
                name
                for name in CURRENT_COLUMNS
                if name in fieldnames
                and any(row.get(name, "") != "" for row in rows)
            ),
            None,
        )
    )
    if bias_column is None:
        raise CurveContractError("missing_bias_column", str(path))
    if current_column is None:
        raise CurveContractError("missing_current_column", str(path))
    points: list[CurvePoint] = []
    for index, row in enumerate(rows):
        points.append(
            CurvePoint(
                bias_V=_finite(row[bias_column], f"{path}:{index}:bias"),
                current_A_per_um=_finite(
                    row[current_column],
                    f"{path}:{index}:current",
                ),
                converged=_as_bool(row.get("converged")),
                failure_reason=(
                    row.get("failure_reason")
                    or row.get("newton_failure_class")
                    or row.get("outcome")
                    or ""
                ),
                source_index=index,
                electron_closure_relative=_optional_finite(
                    row,
                    "electron_closure_relative",
                ),
                hole_closure_relative=_optional_finite(
                    row,
                    "hole_closure_relative",
                ),
                terminal_pair_closure_A_per_um=_optional_finite(
                    row,
                    "terminal_pair_closure_A_per_um",
                ),
                internal_kcl_relative=_optional_finite(
                    row,
                    "internal_kcl_relative",
                ),
            )
        )
    return points


def first_solver_failure(
    points: Sequence[CurvePoint],
) -> dict[str, Any] | None:
    for point in points:
        if not point.converged:
            return {
                "bias_V": point.bias_V,
                "failure_reason": point.failure_reason or "unclassified_solver_failure",
                "source_index": point.source_index,
            }
    return None


def validate_curve_points(
    points: Sequence[CurvePoint],
    curve_name: str,
) -> None:
    if not points:
        raise CurveContractError("empty_curve", curve_name)
    for index, point in enumerate(points):
        _finite(point.bias_V, f"{curve_name}:{index}:bias")
        _finite(point.current_A_per_um, f"{curve_name}:{index}:current")
        for field in (
            "electron_closure_relative",
            "hole_closure_relative",
            "terminal_pair_closure_A_per_um",
            "internal_kcl_relative",
        ):
            value = getattr(point, field)
            if value is not None:
                _finite(value, f"{curve_name}:{index}:{field}")


def exact_curve_index(
    points: Sequence[CurvePoint],
    targets: Sequence[float],
    *,
    curve_name: str,
) -> tuple[dict[float, CurvePoint], list[float]]:
    result: dict[float, CurvePoint] = {}
    missing: list[float] = []
    for target in targets:
        matches = [
            point
            for point in points
            if point.converged
            and abs(point.bias_V - target) <= EXACT_BIAS_TOLERANCE_V
        ]
        if len(matches) > 1:
            raise CurveContractError(
                "duplicate_exact_bias",
                f"{curve_name}:{target:g}",
            )
        if not matches:
            missing.append(target)
            continue
        result[target] = matches[0]
    return result, missing


def percentile(values: Sequence[float], fraction: float) -> float:
    if not values:
        raise CurveContractError("empty_metric", "percentile")
    ordered = sorted(values)
    rank = fraction * (len(ordered) - 1)
    lower = math.floor(rank)
    upper = math.ceil(rank)
    if lower == upper:
        return ordered[lower]
    weight = rank - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def log_current(point: CurvePoint, label: str) -> float:
    magnitude = abs(point.current_A_per_um)
    if magnitude <= PLOTTING_FLOOR_A_PER_UM:
        raise CurveContractError(
            "numerical_floor_row",
            f"{label}:{point.bias_V:g}",
        )
    return math.log10(magnitude)


def adjacent_slopes(
    curve: Mapping[float, CurvePoint],
    biases: Sequence[float],
    label: str,
) -> list[dict[str, float]]:
    rows: list[dict[str, float]] = []
    for left_bias, right_bias in zip(biases, biases[1:]):
        left_y = log_current(curve[left_bias], label)
        right_y = log_current(curve[right_bias], label)
        delta = abs(right_bias) - abs(left_bias)
        if delta <= 0.0:
            raise CurveContractError("invalid_bias_order", label)
        rows.append(
            {
                "left_bias_V": left_bias,
                "right_bias_V": right_bias,
                "midpoint_bias_V": -(abs(left_bias) + abs(right_bias)) / 2.0,
                "slope_dex_per_V": (right_y - left_y) / delta,
            }
        )
    return rows


def slope_knee(slopes: Sequence[Mapping[str, float]]) -> float | None:
    for index, row in enumerate(slopes):
        if row["slope_dex_per_V"] < 1.0:
            continue
        if (
            index + 1 < len(slopes)
            and slopes[index + 1]["slope_dex_per_V"] < 1.0
        ):
            continue
        if index == 0:
            return float(row["midpoint_bias_V"])
        previous = float(slopes[index - 1]["slope_dex_per_V"])
        current = float(row["slope_dex_per_V"])
        left = abs(float(slopes[index - 1]["midpoint_bias_V"]))
        right = abs(float(row["midpoint_bias_V"]))
        if current == previous:
            return -right
        fraction = (1.0 - previous) / (current - previous)
        return -(left + min(max(fraction, 0.0), 1.0) * (right - left))
    return None


def solve_3x3(matrix: list[list[float]], vector: list[float]) -> list[float]:
    augmented = [row[:] + [value] for row, value in zip(matrix, vector, strict=True)]
    for column in range(3):
        pivot = max(range(column, 3), key=lambda row: abs(augmented[row][column]))
        if abs(augmented[pivot][column]) <= 1.0e-18:
            raise CurveContractError("singular_knee_fit", str(column))
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        scale = augmented[column][column]
        augmented[column] = [value / scale for value in augmented[column]]
        for row in range(3):
            if row == column:
                continue
            factor = augmented[row][column]
            augmented[row] = [
                value - factor * pivot_value
                for value, pivot_value in zip(
                    augmented[row],
                    augmented[column],
                    strict=True,
                )
            ]
    return [augmented[index][3] for index in range(3)]


def continuous_breakpoint(
    curve: Mapping[float, CurvePoint],
    biases: Sequence[float],
    label: str,
) -> float:
    x_values = [abs(bias) for bias in biases]
    y_values = [log_current(curve[bias], label) for bias in biases]
    lower = x_values[2]
    upper = x_values[-3]
    steps = max(1, math.ceil((upper - lower) / 0.001))
    best_knot = lower
    best_sse = math.inf
    for step in range(steps + 1):
        knot = lower + (upper - lower) * step / steps
        design = [[1.0, x, max(0.0, x - knot)] for x in x_values]
        normal = [
            [
                sum(row[left] * row[right] for row in design)
                for right in range(3)
            ]
            for left in range(3)
        ]
        rhs = [
            sum(row[column] * y for row, y in zip(design, y_values, strict=True))
            for column in range(3)
        ]
        coefficients = solve_3x3(normal, rhs)
        sse = sum(
            (
                y
                - sum(
                    coefficient * value
                    for coefficient, value in zip(
                        coefficients,
                        row,
                        strict=True,
                    )
                )
            )
            ** 2
            for row, y in zip(design, y_values, strict=True)
        )
        if sse < best_sse:
            best_sse = sse
            best_knot = knot
    return -best_knot


def curvature_knee(slopes: Sequence[Mapping[str, float]]) -> float:
    if len(slopes) < 2:
        raise CurveContractError("insufficient_knee_slopes", "")
    changes = [
        float(right["slope_dex_per_V"]) - float(left["slope_dex_per_V"])
        for left, right in zip(slopes, slopes[1:])
    ]
    index = max(range(len(changes)), key=changes.__getitem__) + 1
    return float(slopes[index]["midpoint_bias_V"])


def monotonic(points: Mapping[float, CurvePoint], biases: Sequence[float]) -> bool:
    magnitudes = [abs(points[bias].current_A_per_um) for bias in biases]
    return all(
        right >= left
        for left, right in zip(magnitudes, magnitudes[1:])
    )


def closure_guard(points: Mapping[float, CurvePoint]) -> dict[str, Any]:
    fields = (
        ("electron_closure_relative", 1.0e-5),
        ("hole_closure_relative", 1.0e-5),
        ("terminal_pair_closure_A_per_um", 1.0e-20),
        ("internal_kcl_relative", 1.0e-8),
    )
    missing: list[str] = []
    maxima: dict[str, float] = {}
    passed = True
    for name, threshold in fields:
        values = [getattr(point, name) for point in points.values()]
        if any(value is None for value in values):
            missing.append(name)
            passed = False
            continue
        maximum = max(abs(float(value)) for value in values if value is not None)
        maxima[name] = maximum
        passed = passed and maximum <= threshold
    return {"passed": passed, "missing_fields": missing, "maxima": maxima}


def curve_error_metrics(
    left: Mapping[float, CurvePoint],
    right: Mapping[float, CurvePoint],
    biases: Sequence[float],
) -> tuple[list[dict[str, float]], dict[str, float]]:
    rows: list[dict[str, float]] = []
    errors: list[float] = []
    for bias in biases:
        error = abs(
            log_current(left[bias], "left") - log_current(right[bias], "right")
        )
        errors.append(error)
        rows.append({"bias_V": bias, "absolute_log_error_dex": error})
    return rows, {
        "median_absolute_log_error_dex": statistics.median(errors),
        "p95_absolute_log_error_dex": percentile(errors, 0.95),
        "maximum_absolute_log_error_dex": max(errors),
    }


def gain_error_metrics(
    curves: Mapping[str, Mapping[float, CurvePoint]],
    biases: Sequence[float],
) -> dict[str, float]:
    errors = []
    for bias in biases:
        vela_gain = abs(curves["vela_on"][bias].current_A_per_um) / abs(
            curves["vela_off"][bias].current_A_per_um
        )
        sentaurus_gain = abs(
            curves["sentaurus_on"][bias].current_A_per_um
        ) / abs(curves["sentaurus_off"][bias].current_A_per_um)
        if vela_gain <= 0.0 or sentaurus_gain <= 0.0:
            raise CurveContractError("invalid_avalanche_gain", str(bias))
        errors.append(abs(math.log10(vela_gain / sentaurus_gain)))
    return {
        "median_log_error_dex": statistics.median(errors),
        "maximum_log_error_dex": max(errors),
    }


def comparison_rows(
    curves: Mapping[str, Mapping[float, CurvePoint]],
    biases: Sequence[float],
) -> list[dict[str, float]]:
    rows = []
    for bias in biases:
        vela_on = curves["vela_on"][bias].current_A_per_um
        vela_off = curves["vela_off"][bias].current_A_per_um
        sentaurus_on = curves["sentaurus_on"][bias].current_A_per_um
        sentaurus_off = curves["sentaurus_off"][bias].current_A_per_um
        rows.append(
            {
                "bias_V": bias,
                "vela_on_A_per_um": vela_on,
                "vela_off_A_per_um": vela_off,
                "sentaurus_on_A_per_um": sentaurus_on,
                "sentaurus_off_A_per_um": sentaurus_off,
                "absolute_log_current_error_dex": abs(
                    math.log10(abs(vela_on / sentaurus_on))
                ),
                "vela_gain": abs(vela_on / vela_off),
                "sentaurus_gain": abs(sentaurus_on / sentaurus_off),
            }
        )
    return rows


def analyze_curves(
    raw_curves: Mapping[str, Sequence[CurvePoint]],
) -> dict[str, Any]:
    if set(raw_curves) != set(CURVE_NAMES):
        raise CurveContractError("curve_set_mismatch", str(sorted(raw_curves)))
    failures = {
        name: failure
        for name, points in raw_curves.items()
        if (failure := first_solver_failure(points)) is not None
    }
    indexed: dict[str, dict[float, CurvePoint]] = {}
    missing: dict[str, dict[str, list[float]]] = {}
    for name, points in raw_curves.items():
        validate_curve_points(points, name)
        global_rows, missing_global = exact_curve_index(
            points,
            GLOBAL_BIASES_V,
            curve_name=name,
        )
        knee_rows, missing_knee = exact_curve_index(
            points,
            KNEE_BIASES_V,
            curve_name=name,
        )
        indexed[name] = {**global_rows, **knee_rows}
        if missing_global or missing_knee:
            missing[name] = {
                "global_biases_V": missing_global,
                "knee_biases_V": missing_knee,
            }
    if missing:
        outcome = "solver_first_failure" if failures else "incomplete_exact_lattice"
        return {
            "schema": "vela.pn2d_avalanche_on_bv_parity.v1",
            "outcome": outcome,
            "missing_exact_rows": missing,
            "solver_failures": failures,
        }

    global_rows, global_metrics = curve_error_metrics(
        indexed["vela_on"],
        indexed["sentaurus_on"],
        GLOBAL_BIASES_V,
    )
    knee_rows, knee_metrics = curve_error_metrics(
        indexed["vela_on"],
        indexed["sentaurus_on"],
        KNEE_BIASES_V,
    )
    vela_slopes = adjacent_slopes(indexed["vela_on"], KNEE_BIASES_V, "vela_on")
    sentaurus_slopes = adjacent_slopes(
        indexed["sentaurus_on"],
        KNEE_BIASES_V,
        "sentaurus_on",
    )
    slope_errors = [
        float(vela["slope_dex_per_V"]) - float(sentaurus["slope_dex_per_V"])
        for vela, sentaurus in zip(vela_slopes, sentaurus_slopes, strict=True)
    ]
    slope_rmse = math.sqrt(
        sum(error * error for error in slope_errors) / len(slope_errors)
    )
    knee_estimators: dict[str, dict[str, float | None]] = {}
    for name, slopes in (
        ("vela", vela_slopes),
        ("sentaurus", sentaurus_slopes),
    ):
        curve = indexed[f"{name}_on"]
        knee_estimators[name] = {
            "V_slope": slope_knee(slopes),
            "V_break": continuous_breakpoint(curve, KNEE_BIASES_V, f"{name}_on"),
            "V_curvature": curvature_knee(slopes),
        }
    ill_conditioned = any(
        metrics["V_slope"] is None
        or abs(float(metrics["V_slope"]) - float(metrics["V_break"])) > 0.20
        for metrics in knee_estimators.values()
    )
    gain = gain_error_metrics(indexed, KNEE_BIASES_V)
    closure = closure_guard(indexed["vela_on"])
    monotonicity = all(
        monotonic(indexed[name], GLOBAL_BIASES_V)
        and monotonic(indexed[name], KNEE_BIASES_V)
        for name in ("vela_on", "sentaurus_on")
    )
    gates = {
        "global_median": global_metrics["median_absolute_log_error_dex"] <= 0.05,
        "global_p95": global_metrics["p95_absolute_log_error_dex"] <= 0.10,
        "global_maximum": global_metrics["maximum_absolute_log_error_dex"] <= 0.15,
        "knee_median": knee_metrics["median_absolute_log_error_dex"] <= 0.05,
        "knee_maximum": knee_metrics["maximum_absolute_log_error_dex"] <= 0.10,
        "slope_knee": (
            knee_estimators["vela"]["V_slope"] is not None
            and knee_estimators["sentaurus"]["V_slope"] is not None
            and abs(
                float(knee_estimators["vela"]["V_slope"])
                - float(knee_estimators["sentaurus"]["V_slope"])
            )
            <= 0.10
        ),
        "break_knee": abs(
            float(knee_estimators["vela"]["V_break"])
            - float(knee_estimators["sentaurus"]["V_break"])
        )
        <= 0.10,
        "slope_rmse": slope_rmse <= 0.20,
        "monotonicity": monotonicity,
        "gain_median": gain["median_log_error_dex"] <= 0.05,
        "gain_maximum": gain["maximum_log_error_dex"] <= 0.10,
        "closure": closure["passed"],
    }
    if ill_conditioned:
        outcome = "ill_conditioned_knee_metric"
    elif all(gates.values()):
        outcome = "curve_knee_parity_passed"
    else:
        outcome = "curve_knee_gate_failed"
    return {
        "schema": "vela.pn2d_avalanche_on_bv_parity.v1",
        "outcome": outcome,
        "global_metrics": global_metrics,
        "knee_metrics": knee_metrics,
        "gain_metrics": gain,
        "knee_estimators": knee_estimators,
        "adjacent_slope_rmse_dex_per_V": slope_rmse,
        "closure": closure,
        "gates": gates,
        "curve_rows": comparison_rows(
            indexed,
            tuple(sorted(set(GLOBAL_BIASES_V) | set(KNEE_BIASES_V), reverse=True)),
        ),
        "global_error_rows": global_rows,
        "knee_error_rows": knee_rows,
        "vela_slopes": vela_slopes,
        "sentaurus_slopes": sentaurus_slopes,
    }


def write_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="ascii", newline="\n")
        return
    with path.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_svg(
    path: Path,
    series: Mapping[str, Sequence[tuple[float, float]]],
    title: str,
) -> None:
    width, height = 900, 520
    margin = 55
    all_points = [point for points in series.values() for point in points]
    if not all_points:
        return
    xs = [point[0] for point in all_points]
    ys = [point[1] for point in all_points]
    x_min, x_max = min(xs), max(xs)
    y_min, y_max = min(ys), max(ys)
    if y_min == y_max:
        y_min -= 1.0
        y_max += 1.0

    def transform(point: tuple[float, float]) -> tuple[float, float]:
        x, y = point
        px = margin + (x - x_min) / max(x_max - x_min, 1.0) * (
            width - 2 * margin
        )
        py = height - margin - (y - y_min) / (y_max - y_min) * (
            height - 2 * margin
        )
        return px, py

    colors = ("#2563eb", "#dc2626", "#059669", "#7c3aed")
    polylines = []
    legend = []
    for index, (label, points) in enumerate(series.items()):
        color = colors[index % len(colors)]
        coordinates = " ".join(
            f"{x:.2f},{y:.2f}" for x, y in map(transform, points)
        )
        polylines.append(
            f'<polyline fill="none" stroke="{color}" stroke-width="2" '
            f'points="{coordinates}"/>'
        )
        legend.append(
            f'<text x="{margin + index * 180}" y="30" fill="{color}">{label}</text>'
        )
    svg = f"""\
<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">
<rect width="100%" height="100%" fill="white"/>
<text x="{width / 2}" y="20" text-anchor="middle">{title}</text>
<line x1="{margin}" y1="{height-margin}" x2="{width-margin}" y2="{height-margin}" stroke="black"/>
<line x1="{margin}" y1="{margin}" x2="{margin}" y2="{height-margin}" stroke="black"/>
{''.join(legend)}
{''.join(polylines)}
</svg>
"""
    path.write_text(svg, encoding="ascii", newline="\n")


def manifest_failure(
    manifest: Mapping[str, Any],
    curve_name: str,
) -> dict[str, Any] | None:
    for key in ("first_failure",):
        value = manifest.get(key)
        if isinstance(value, Mapping):
            return dict(value)
    if not curve_name.endswith("_on"):
        return None
    for key in ("avalanche_on_run_a", "avalanche_on_run_b"):
        value = manifest.get(key)
        if isinstance(value, Mapping) and isinstance(value.get("first_failure"), Mapping):
            return dict(value["first_failure"])
    return None


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in CURVE_NAMES:
        option = name.replace("_", "-")
        parser.add_argument(f"--{option}-csv", type=Path, required=True)
        parser.add_argument(f"--{option}-manifest", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args(argv)
    curves: dict[str, list[CurvePoint]] = {}
    inputs: dict[str, Any] = {}
    manifest_failures: dict[str, Any] = {}
    for name in CURVE_NAMES:
        csv_path = getattr(args, f"{name}_csv").resolve()
        manifest_path = getattr(args, f"{name}_manifest").resolve()
        curves[name] = load_curve(csv_path)
        manifest = json.loads(manifest_path.read_text(encoding="ascii"))
        failure = manifest_failure(manifest, name)
        if failure is not None:
            manifest_failures[name] = failure
        inputs[name] = {
            "csv": str(csv_path),
            "csv_sha256": sha256(csv_path),
            "manifest": str(manifest_path),
            "manifest_sha256": sha256(manifest_path),
        }
    result = analyze_curves(curves)
    if result["outcome"] == "incomplete_exact_lattice" and manifest_failures:
        result["outcome"] = "solver_first_failure"
    result["manifest_solver_failures"] = manifest_failures
    result["inputs"] = inputs

    output = args.output_root.resolve()
    output.mkdir(parents=True, exist_ok=True)
    write_csv(
        output / "curve_points.csv",
        result.get("curve_rows", []),
    )
    slope_rows = [
        {"simulator": "vela", **row} for row in result.get("vela_slopes", [])
    ] + [
        {"simulator": "sentaurus", **row}
        for row in result.get("sentaurus_slopes", [])
    ]
    write_csv(output / "slope_points.csv", slope_rows)
    (output / "knee_metrics.json").write_text(
        json.dumps(result.get("knee_estimators", {}), indent=2, sort_keys=True)
        + "\n",
        encoding="ascii",
        newline="\n",
    )
    (output / "acceptance.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n",
        encoding="ascii",
        newline="\n",
    )
    if result.get("vela_slopes") and result.get("sentaurus_slopes"):
        write_svg(
            output / "local_slope_comparison.svg",
            {
                "Vela": [
                    (abs(row["midpoint_bias_V"]), row["slope_dex_per_V"])
                    for row in result["vela_slopes"]
                ],
                "Sentaurus": [
                    (abs(row["midpoint_bias_V"]), row["slope_dex_per_V"])
                    for row in result["sentaurus_slopes"]
                ],
            },
            "PN2D BV local log-current slope",
        )
    if result.get("curve_rows"):
        rows = result["curve_rows"]
        write_svg(
            output / "linear_current_comparison.svg",
            {
                "Vela on": [
                    (abs(row["bias_V"]), abs(row["vela_on_A_per_um"]))
                    for row in rows
                ],
                "Sentaurus on": [
                    (abs(row["bias_V"]), abs(row["sentaurus_on_A_per_um"]))
                    for row in rows
                ],
            },
            "PN2D BV terminal current",
        )
        write_svg(
            output / "log_current_comparison.svg",
            {
                "Vela on": [
                    (
                        abs(row["bias_V"]),
                        math.log10(abs(row["vela_on_A_per_um"])),
                    )
                    for row in rows
                ],
                "Sentaurus on": [
                    (
                        abs(row["bias_V"]),
                        math.log10(abs(row["sentaurus_on_A_per_um"])),
                    )
                    for row in rows
                ],
            },
            "PN2D BV log terminal current",
        )
        write_svg(
            output / "gain_comparison.svg",
            {
                "Vela gain": [
                    (abs(row["bias_V"]), row["vela_gain"]) for row in rows
                ],
                "Sentaurus gain": [
                    (abs(row["bias_V"]), row["sentaurus_gain"]) for row in rows
                ],
            },
            "PN2D avalanche gain",
        )
    print(json.dumps(result, indent=2, sort_keys=True))
    successful = {
        "curve_knee_parity_passed",
        "solver_first_failure",
        "incomplete_exact_lattice",
        "curve_knee_gate_failed",
        "ill_conditioned_knee_metric",
    }
    return 0 if result["outcome"] in successful else 1


if __name__ == "__main__":
    raise SystemExit(main())
