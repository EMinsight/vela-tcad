import math
import unittest
from dataclasses import replace

from scripts.pn2d_minimal6_diagnostics.inverse_contracts import (
    AcceptanceThresholds,
    Identifiability,
    Observation,
    SampleStatus,
    SupportKind,
)
from scripts.pn2d_minimal6_diagnostics.inverse_transport import (
    TransportVectorError,
    current_inverted_qf_gradient,
    evaluate_transport_candidates,
    project_vector_to_edge,
    qf_current_density,
    reconstruct_edge_vector,
    summarize_transport_errors,
)


class InverseTransportTest(unittest.TestCase):
    def test_qf_current_density_uses_explicit_carrier_signs_and_si_dimensions(self):
        self.assertEqual(
            qf_current_density("electron", 2.0, 3.0, (4.0, -5.0), q=1.0),
            (-24.0, 30.0),
        )
        self.assertEqual(
            qf_current_density("hole", 2.0, 3.0, (4.0, -5.0), q=1.0),
            (24.0, -30.0),
        )
        self.assertEqual(
            qf_current_density("electron", 2.0, 3.0, (0.0, 0.0), q=1.0),
            (0.0, 0.0),
        )

        # C * (m^2/V/s) * (1/m^3) * (V/m) = A/m^2.
        value = qf_current_density(
            "hole", 2.0e21, 0.15, (40.0, -50.0), q=1.602176634e-19
        )
        self.assertAlmostEqual(value[0], 1922.6119608)
        self.assertAlmostEqual(value[1], -2403.264951)

    def test_qf_current_density_rejects_unknown_negative_and_nonfinite_inputs(self):
        with self.assertRaisesRegex(ValueError, "carrier"):
            qf_current_density("ambipolar", 2.0, 3.0, (4.0, -5.0))
        for density, mobility in ((-1.0, 1.0), (1.0, -1.0)):
            with self.subTest(density=density, mobility=mobility):
                with self.assertRaisesRegex(ValueError, "non-negative"):
                    qf_current_density("electron", density, mobility, (1.0, 0.0))
        for bad in (math.nan, math.inf):
            with self.subTest(bad=bad):
                with self.assertRaisesRegex(ValueError, "finite"):
                    qf_current_density("electron", 1.0, 1.0, (bad, 0.0))

    def test_current_inversion_is_carrier_specific_and_fails_closed_at_floors(self):
        self.assertEqual(
            current_inverted_qf_gradient(
                "electron", 2.0, 3.0, (-24.0, 30.0), q=1.0
            ),
            (4.0, -5.0),
        )
        self.assertEqual(
            current_inverted_qf_gradient(
                "hole", 2.0, 3.0, (24.0, -30.0), q=1.0
            ),
            (4.0, -5.0),
        )
        with self.assertRaisesRegex(ValueError, "density is at or below floor"):
            current_inverted_qf_gradient(
                "electron", 1.0e-12, 3.0, (1.0, 0.0), q=1.0,
                density_floor=1.0e-12,
            )
        with self.assertRaisesRegex(ValueError, "mobility is at or below floor"):
            current_inverted_qf_gradient(
                "electron", 2.0, 0.0, (1.0, 0.0), q=1.0,
                mobility_floor=0.0,
            )
        with self.assertRaisesRegex(ValueError, "current must be finite"):
            current_inverted_qf_gradient(
                "electron", 2.0, 3.0, (math.inf, 0.0), q=1.0
            )

    def test_projection_and_reconstruction_preserve_directed_edge_orientation(self):
        self.assertEqual(project_vector_to_edge((3.0, 4.0), (0.0, 0.0), (2.0, 0.0)), 3.0)
        self.assertEqual(project_vector_to_edge((3.0, 4.0), (2.0, 0.0), (0.0, 0.0)), -3.0)
        self.assertEqual(reconstruct_edge_vector(3.0, (0.0, 0.0), (2.0, 0.0)), (3.0, 0.0))
        self.assertEqual(reconstruct_edge_vector(3.0, (2.0, 0.0), (0.0, 0.0)), (-3.0, 0.0))
        diagonal = reconstruct_edge_vector(5.0, (0.0, 0.0), (3.0, 4.0))
        self.assertEqual(diagonal, (3.0, 4.0))
        self.assertAlmostEqual(project_vector_to_edge(diagonal, (0.0, 0.0), (3.0, 4.0)), 5.0)
        with self.assertRaisesRegex(ValueError, "zero length"):
            project_vector_to_edge((1.0, 0.0), (0.0, 0.0), (0.0, 0.0))
        with self.assertRaisesRegex(ValueError, "finite"):
            reconstruct_edge_vector(math.inf, (0.0, 0.0), (1.0, 0.0))

    def test_fixed_transport_gates_include_exact_boundaries_and_reject_overrides(self):
        boundary = summarize_transport_errors(
            (
                TransportVectorError(SampleStatus.VALID, SampleStatus.VALID, 0.0, 5.0),
                TransportVectorError(SampleStatus.VALID, SampleStatus.VALID, 0.1, 5.0),
                TransportVectorError(SampleStatus.VALID, SampleStatus.VALID, 0.3, 5.0),
            )
        )
        self.assertEqual(boundary.median_abs_log10_error, 0.1)
        self.assertEqual(boundary.p95_abs_log10_error, 0.3)
        self.assertEqual(boundary.median_angle_deg, 5.0)
        self.assertIs(boundary.classification, Identifiability.IDENTIFIED)

        above_median = summarize_transport_errors(
            (
                TransportVectorError(SampleStatus.VALID, SampleStatus.VALID, 0.1000001, 0.0),
            )
        )
        self.assertIs(above_median.classification, Identifiability.REJECTED)
        above_p95 = summarize_transport_errors(
            (
                TransportVectorError(SampleStatus.VALID, SampleStatus.VALID, 0.0, 0.0),
                TransportVectorError(SampleStatus.VALID, SampleStatus.VALID, 0.1, 0.0),
                TransportVectorError(SampleStatus.VALID, SampleStatus.VALID, 0.3000001, 0.0),
            )
        )
        self.assertIs(above_p95.classification, Identifiability.REJECTED)
        above_angle = summarize_transport_errors(
            (TransportVectorError(SampleStatus.VALID, SampleStatus.VALID, 0.0, 5.000001),)
        )
        self.assertIs(above_angle.classification, Identifiability.REJECTED)

        with self.assertRaisesRegex(ValueError, "immutable"):
            summarize_transport_errors(
                boundary.errors,
                thresholds=AcceptanceThresholds(gradient_median_abs_dex=0.2),
            )

    def test_evaluator_recovers_all_carrier_resolved_current_semantics(self):
        observations = self.make_affine_observations(solver="sentaurus", mobility=True)
        results = evaluate_transport_candidates(
            observations, self.mesh(), density_floor=1.0e-30, current_floor=1.0e-30,
            thermal_voltage_V=0.025, q=1.0,
        )

        expected = {
            "triangle_qf_gradient_current",
            "node_area_weighted_qf_gradient_current",
            "edge_area_weighted_qf_gradient_current",
            "signed_edge_qf_difference_current",
            "current_inverted_qf_gradient",
            "signed_edge_sg_density_current",
            "signed_edge_drift_diffusion_current",
        }
        self.assertEqual({item.candidate for item in results}, expected)
        self.assertEqual({item.carrier for item in results}, {"electron", "hole"})
        self.assertTrue(all(item.split == "discovery" for item in results))
        self.assertTrue(all(item.classification is Identifiability.IDENTIFIED for item in results))
        self.assertTrue(all(item.unit_si in {"A/m^2", "V/m"} for item in results))
        self.assertTrue(all(item.median_abs_log10_error == 0.0 for item in results))
        self.assertTrue(all(item.median_angle_deg == 0.0 for item in results))

        electron_node = next(
            item for item in results
            if item.carrier == "electron"
            and item.candidate == "node_area_weighted_qf_gradient_current"
        )
        self.assertIs(electron_node.support_kind, SupportKind.NODE)
        self.assertEqual(electron_node.samples[0].candidate_value, (-16.0, 12.0))
        self.assertEqual(electron_node.samples[0].reference_value, (-16.0, 12.0))
        self.assertEqual(electron_node.samples[0].candidate_support_transform,
                         "cell_qf_gradient_to_node_area_weighted")
        self.assertEqual(electron_node.samples[0].reference_support_transform,
                         "native_node_current_vector")

        signed_edge = next(
            item for item in results
            if item.carrier == "hole"
            and item.candidate == "signed_edge_qf_difference_current"
        )
        self.assertIs(signed_edge.support_kind, SupportKind.EDGE)
        self.assertTrue(all(sample.orientation == "global_xy;edge=start_to_end"
                            for sample in signed_edge.samples))

        reversed_observations = tuple(reversed(observations))
        self.assertEqual(
            results,
            evaluate_transport_candidates(
                reversed_observations, self.mesh(), density_floor=1.0e-30,
                current_floor=1.0e-30, thermal_voltage_V=0.025, q=1.0,
            ),
        )

    def test_zero_field_sg_diffusion_signs_follow_both_carrier_qf_definitions(self):
        rows = []
        thermal_voltage = 0.5
        for row in self.make_affine_observations(solver="sentaurus", mobility=True):
            node_density = 3.0 if str(row.support_id) == "1" else 1.0
            if row.quantity == "ElectrostaticPotential":
                row = replace(row, raw_value=0.0, value_si=0.0)
            elif row.quantity in {"eDensity", "hDensity"}:
                row = replace(row, raw_value=node_density, value_si=node_density)
            elif row.quantity in {"eQuasiFermiPotential", "hQuasiFermiPotential"}:
                qf = -thermal_voltage * math.log(node_density)
                row = replace(row, raw_value=qf, value_si=qf)
            elif row.quantity in {"eCurrentDensity", "hCurrentDensity"}:
                current_x = 1.0 if row.quantity == "eCurrentDensity" else -1.0
                current_value = (
                    current_x if row.component == "component0" else 0.0
                )
                row = replace(row, raw_value=current_value, value_si=current_value)
            rows.append(row)

        results = evaluate_transport_candidates(
            rows, self.mesh(), density_floor=1.0e-30, current_floor=1.0e-30,
            thermal_voltage_V=thermal_voltage, q=1.0,
        )
        sg = {
            item.carrier: next(
                sample for sample in item.samples if sample.support_id == "edge-01"
            )
            for item in results if item.candidate == "signed_edge_sg_density_current"
        }
        self.assertEqual(sg["electron"].candidate_value, (1.0, 0.0))
        self.assertEqual(sg["hole"].candidate_value, (-1.0, 0.0))

    def test_candidates_remain_computable_when_current_reference_is_below_floor(self):
        rows = [
            replace(row, raw_value=0.0, value_si=0.0)
            if row.quantity in {"eCurrentDensity", "hCurrentDensity"} else row
            for row in self.make_affine_observations(solver="sentaurus", mobility=True)
        ]
        results = evaluate_transport_candidates(
            rows, self.mesh(), density_floor=1.0e-30, current_floor=1.0e-12,
            thermal_voltage_V=0.025, q=1.0,
        )
        inverse = [
            item for item in results
            if item.candidate == "current_inverted_qf_gradient"
        ]
        self.assertTrue(inverse)
        self.assertTrue(all(
            sample.candidate_value is None
            for item in inverse for sample in item.samples
        ))
        self.assertTrue(all(
            sample.error.magnitude_status is SampleStatus.BELOW_FLOOR
            for item in inverse for sample in item.samples
        ))
        independently_computable = [
            item for item in results if item.candidate != "current_inverted_qf_gradient"
        ]
        self.assertTrue(independently_computable)
        self.assertTrue(all(
            sample.candidate_value is not None
            for item in independently_computable for sample in item.samples
        ))
        self.assertTrue(any(
            sample.candidate_value != (0.0, 0.0)
            for item in independently_computable for sample in item.samples
        ))
        self.assertTrue(all(
            sample.error.magnitude_status is SampleStatus.BELOW_FLOOR
            and sample.error.direction_status is SampleStatus.DIRECTION_UNDEFINED
            for item in independently_computable for sample in item.samples
        ))
        self.assertTrue(all(item.valid_count == 0 for item in independently_computable))
        self.assertTrue(all(item.classification is Identifiability.INSUFFICIENT_DATA
                            for item in independently_computable))

    def test_flat_qf_produces_explicit_zero_reconstructed_qf_current(self):
        rows = []
        for row in self.make_affine_observations(solver="sentaurus", mobility=True):
            if row.quantity in {"eQuasiFermiPotential", "hQuasiFermiPotential"}:
                row = replace(row, raw_value=2.0, value_si=2.0)
            if row.quantity in {"eCurrentDensity", "hCurrentDensity"}:
                row = replace(row, raw_value=0.0, value_si=0.0)
            rows.append(row)
        results = evaluate_transport_candidates(
            rows, self.mesh(), density_floor=1.0e-30, current_floor=1.0e-30,
            q=1.0,
        )
        qf_candidates = {
            "triangle_qf_gradient_current",
            "node_area_weighted_qf_gradient_current",
            "edge_area_weighted_qf_gradient_current",
            "signed_edge_qf_difference_current",
        }
        samples = [
            sample for item in results if item.candidate in qf_candidates
            for sample in item.samples
        ]
        self.assertTrue(samples)
        self.assertTrue(all(sample.candidate_value == (0.0, 0.0)
                            for sample in samples))
        self.assertTrue(all(
            sample.error.magnitude_status is SampleStatus.BELOW_FLOOR
            and sample.error.direction_status is SampleStatus.DIRECTION_UNDEFINED
            for sample in samples
        ))

    def test_missing_mobility_is_typed_confounded_with_mu_times_gradient_observable(self):
        results = evaluate_transport_candidates(
            self.make_affine_observations(solver="vela", mobility=False), self.mesh(),
            density_floor=1.0e-30, current_floor=1.0e-30, q=1.0,
        )
        direct = next(
            item for item in results
            if item.carrier == "electron"
            and item.candidate == "node_area_weighted_qf_gradient_current"
        )
        self.assertIs(direct.classification, Identifiability.CONFOUNDED)
        self.assertTrue(direct.confoundings)
        record = direct.confoundings[0]
        self.assertIs(record.status, SampleStatus.MISSING_FIELD)
        self.assertIs(record.classification, Identifiability.CONFOUNDED)
        self.assertEqual(record.observable, "mu_times_grad_qf")
        self.assertEqual(record.observable_value, (8.0, -6.0))
        self.assertEqual(record.observable_unit_si, "m/s")
        self.assertEqual(record.missing_inputs, ("mobility",))

    def test_density_floor_and_nonfinite_current_never_divide_or_fabricate_zero(self):
        rows = list(self.make_affine_observations(solver="sentaurus", mobility=True))
        rows = [
            replace(row, raw_value=0.0, value_si=0.0)
            if row.quantity == "eDensity" else row
            for row in rows
        ]
        results = evaluate_transport_candidates(
            rows, self.mesh(), density_floor=1.0e-12, current_floor=1.0e-30,
            q=1.0,
        )
        inverse = next(
            item for item in results
            if item.carrier == "electron"
            and item.candidate == "current_inverted_qf_gradient"
        )
        self.assertIs(inverse.classification, Identifiability.INSUFFICIENT_DATA)
        self.assertTrue(inverse.confoundings)
        self.assertTrue(all(record.status is SampleStatus.BELOW_FLOOR
                            for record in inverse.confoundings))
        self.assertTrue(all(record.observable_value is None
                            for record in inverse.confoundings))
        self.assertTrue(all(sample.candidate_value is None for sample in inverse.samples))

        nonfinite_rows = [
            replace(row, raw_value=math.inf, value_si=math.inf,
                    status=SampleStatus.NONFINITE)
            if row.quantity == "hCurrentDensity" and row.component == "component0" else row
            for row in self.make_affine_observations(solver="sentaurus", mobility=True)
        ]
        nonfinite = evaluate_transport_candidates(
            nonfinite_rows, self.mesh(), density_floor=1.0e-30, current_floor=1.0e-30,
            q=1.0,
        )
        hole_inverse = next(
            item for item in nonfinite
            if item.carrier == "hole" and item.candidate == "current_inverted_qf_gradient"
        )
        self.assertTrue(any(record.status is SampleStatus.NONFINITE
                            for record in hole_inverse.confoundings))

    @staticmethod
    def mesh():
        return {
            "triangles": {"cell-0": ("0", "1", "2")},
            "edges": {
                "edge-01": ("0", "1"),
                "edge-02": ("0", "2"),
                "edge-12": ("1", "2"),
            },
        }

    @staticmethod
    def make_affine_observations(*, solver, mobility):
        points = {"0": (0.0, 0.0), "1": (2.0, 0.0), "2": (0.0, 1.0)}
        gradient = (4.0, -3.0)
        rows = []

        def row(node, quantity, component, value, unit):
            return Observation(
                solver=solver, topology="sketch", bias_V=-12.0,
                support_kind=SupportKind.NODE, support_id=node, quantity=quantity,
                component=component, raw_value=value, raw_unit=unit, value_si=value,
                unit_si=unit, coordinate_frame="canonical_xy",
                orientation="global_xy", conversion="identity",
                status=SampleStatus.VALID, source_path="state.csv",
                source_sha256="0" * 64,
            )

        for node, (x, y) in points.items():
            qf = gradient[0] * x + gradient[1] * y + 2.0
            rows.extend((
                row(node, "coordinate", "x", x, "m"),
                row(node, "coordinate", "y", y, "m"),
                row(node, "ElectrostaticPotential", "component0", qf, "V"),
                row(node, "eQuasiFermiPotential", "component0", qf, "V"),
                row(node, "hQuasiFermiPotential", "component0", qf, "V"),
                row(node, "eDensity", "component0", 2.0, "m^-3"),
                row(node, "hDensity", "component0", 2.0, "m^-3"),
                row(node, "eCurrentDensity", "component0", -16.0, "A/m^2"),
                row(node, "eCurrentDensity", "component1", 12.0, "A/m^2"),
                row(node, "hCurrentDensity", "component0", 16.0, "A/m^2"),
                row(node, "hCurrentDensity", "component1", -12.0, "A/m^2"),
            ))
            if mobility:
                rows.extend((
                    row(node, "eMobility", "component0", 2.0, "m^2*V^-1*s^-1"),
                    row(node, "hMobility", "component0", 2.0, "m^2*V^-1*s^-1"),
                ))
        return tuple(rows)


if __name__ == "__main__":
    unittest.main()
