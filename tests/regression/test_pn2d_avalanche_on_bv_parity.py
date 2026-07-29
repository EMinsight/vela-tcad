#!/usr/bin/env python3

from __future__ import annotations

import math
import unittest

from scripts.analyze_pn2d_avalanche_on_bv_parity import (
    GLOBAL_BIASES_V,
    KNEE_BIASES_V,
    CurveContractError,
    CurvePoint,
    analyze_curves,
)


ALL_BIASES = tuple(
    sorted(set(GLOBAL_BIASES_V) | set(KNEE_BIASES_V), reverse=True)
)


def current(bias: float, shift: float = 0.0) -> float:
    magnitude = abs(bias)
    log_value = -17.0 + 0.02 * magnitude
    log_value += 1.7 * max(0.0, magnitude - (19.35 + shift)) ** 2
    return 10.0**log_value


def curve(
    *,
    shift: float = 0.0,
    scale: float = 1.0,
    nonmonotonic: bool = False,
) -> list[CurvePoint]:
    rows = []
    for bias in ALL_BIASES:
        value = current(bias, shift) * scale
        if nonmonotonic and bias == -19.8:
            value = current(-19.7, shift) * scale * 0.5
        rows.append(
            CurvePoint(
                bias_V=bias,
                current_A_per_um=value,
                electron_closure_relative=1.0e-8,
                hole_closure_relative=1.0e-8,
                terminal_pair_closure_A_per_um=1.0e-24,
                internal_kcl_relative=1.0e-10,
            )
        )
    return rows


def paired(**vela_on_options: object) -> dict[str, list[CurvePoint]]:
    baseline = curve()
    return {
        "vela_on": curve(**vela_on_options),
        "vela_off": baseline,
        "sentaurus_on": baseline,
        "sentaurus_off": baseline,
    }


class PN2DAvalancheOnBVParityTest(unittest.TestCase):
    def test_identical_curves_pass(self) -> None:
        result = analyze_curves(paired())
        self.assertEqual(result["outcome"], "curve_knee_parity_passed")
        self.assertTrue(all(result["gates"].values()))

    def test_known_vertical_mismatch_fails(self) -> None:
        result = analyze_curves(paired(scale=10.0**0.2))
        self.assertEqual(result["outcome"], "curve_knee_gate_failed")
        self.assertFalse(result["gates"]["global_maximum"])

    def test_known_knee_shift_is_detected(self) -> None:
        result = analyze_curves(paired(shift=0.10))
        detected = abs(
            float(result["knee_estimators"]["vela"]["V_slope"])
            - float(result["knee_estimators"]["sentaurus"]["V_slope"])
        )
        self.assertAlmostEqual(detected, 0.10, delta=0.01)

    def test_correct_knee_with_wrong_post_knee_slope_fails(self) -> None:
        rows = curve()
        modified = []
        for point in rows:
            value = point.current_A_per_um
            if point.bias_V <= -19.8:
                value *= 10.0 ** (2.0 * (abs(point.bias_V) - 19.8))
            modified.append(
                CurvePoint(
                    **{
                        **point.__dict__,
                        "current_A_per_um": value,
                    }
                )
            )
        curves = paired()
        curves["vela_on"] = modified
        result = analyze_curves(curves)
        self.assertEqual(result["outcome"], "curve_knee_gate_failed")
        self.assertFalse(result["gates"]["slope_rmse"])

    def test_nonmonotonic_curve_fails(self) -> None:
        result = analyze_curves(paired(nonmonotonic=True))
        self.assertEqual(result["outcome"], "curve_knee_gate_failed")
        self.assertFalse(result["gates"]["monotonicity"])

    def test_ill_conditioned_knee_metric_is_typed(self) -> None:
        rows = curve()
        shaped = []
        for point in rows:
            magnitude = abs(point.bias_V)
            if magnitude < 18.0:
                log_value = -17.0 + 0.02 * magnitude
            elif magnitude <= 19.7:
                log_value = -17.0 + 1.2 * (magnitude - 18.0)
            else:
                log_value = -17.0 + 1.2 * 1.7 + 4.0 * (magnitude - 19.7)
            shaped.append(
                CurvePoint(
                    **{
                        **point.__dict__,
                        "current_A_per_um": 10.0**log_value,
                    }
                )
            )
        curves = paired()
        curves["vela_on"] = shaped
        curves["sentaurus_on"] = shaped
        result = analyze_curves(curves)
        self.assertEqual(result["outcome"], "ill_conditioned_knee_metric")

    def test_missing_closure_fields_fail_gate(self) -> None:
        curves = paired()
        curves["vela_on"] = [
            CurvePoint(
                bias_V=point.bias_V,
                current_A_per_um=point.current_A_per_um,
            )
            for point in curves["vela_on"]
        ]
        result = analyze_curves(curves)
        self.assertEqual(result["outcome"], "curve_knee_gate_failed")
        self.assertFalse(result["gates"]["closure"])
        self.assertEqual(len(result["closure"]["missing_fields"]), 4)

    def test_row_order_does_not_change_results(self) -> None:
        forward = paired()
        reverse = {
            name: list(reversed(points)) for name, points in forward.items()
        }
        self.assertEqual(analyze_curves(forward), analyze_curves(reverse))

    def test_duplicate_and_missing_rows_fail_closed(self) -> None:
        duplicate = paired()
        duplicate["vela_on"].append(duplicate["vela_on"][0])
        with self.assertRaises(CurveContractError) as captured:
            analyze_curves(duplicate)
        self.assertEqual(captured.exception.reason, "duplicate_exact_bias")

        missing = paired()
        missing["vela_on"] = [
            point for point in missing["vela_on"] if point.bias_V != -19.85
        ]
        result = analyze_curves(missing)
        self.assertEqual(result["outcome"], "incomplete_exact_lattice")
        self.assertEqual(
            result["missing_exact_rows"]["vela_on"]["knee_biases_V"],
            [-19.85],
        )

    def test_solver_failure_takes_precedence_over_missing_rows(self) -> None:
        curves = paired()
        curves["vela_on"] = [
            point for point in curves["vela_on"] if point.bias_V >= -19.5
        ]
        curves["vela_on"].append(
            CurvePoint(
                bias_V=-19.69375,
                current_A_per_um=0.0,
                converged=False,
                failure_reason="max_iterations",
            )
        )
        result = analyze_curves(curves)
        self.assertEqual(result["outcome"], "solver_first_failure")
        self.assertEqual(
            result["solver_failures"]["vela_on"]["failure_reason"],
            "max_iterations",
        )

    def test_numerical_floor_row_fails_typed(self) -> None:
        curves = paired()
        first = curves["vela_on"][0]
        curves["vela_on"][0] = CurvePoint(
            **{**first.__dict__, "current_A_per_um": 1.0e-31}
        )
        with self.assertRaises(CurveContractError) as captured:
            analyze_curves(curves)
        self.assertEqual(captured.exception.reason, "numerical_floor_row")

    def test_nonfinite_row_fails_typed(self) -> None:
        curves = paired()
        first = curves["vela_on"][0]
        curves["vela_on"][0] = CurvePoint(
            **{**first.__dict__, "current_A_per_um": math.nan}
        )
        with self.assertRaises(CurveContractError):
            analyze_curves(curves)


if __name__ == "__main__":
    unittest.main()
