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
from scripts.pn2d_minimal6_diagnostics.inverse_fields import (
    _usable,
    cell_to_edge_vectors,
    cell_to_node_vectors,
    edge_scalar_difference,
    evaluate_field_candidates,
    mirror_vector,
    triangle_gradient,
    vector_error,
)


class InverseFieldsTest(unittest.TestCase):
    def test_usable_preserves_every_preexisting_nonvalid_status(self):
        valid = self.make_observations(reference_scale=1.0)[0]

        for status in SampleStatus:
            if status is SampleStatus.VALID:
                continue
            with self.subTest(status=status):
                value, result_status = _usable(
                    replace(valid, status=status, value_si=1.0, unit_si=valid.unit_si), valid.unit_si
                )
                self.assertIsNone(value)
                self.assertIs(result_status, status)

    def test_triangle_gradient_is_exact_and_rejects_singular_geometry(self):
        points = ((0.0, 0.0), (2.0, 0.0), (0.0, 1.0))
        values = tuple(3.0 * x - 4.0 * y + 2.0 for x, y in points)
        self.assertEqual(triangle_gradient(points, values), (3.0, -4.0))

        with self.assertRaisesRegex(ValueError, "degenerate triangle"):
            triangle_gradient(((0.0, 0.0), (1.0, 0.0), (2.0, 0.0)), (0.0, 1.0, 2.0))


        with self.assertRaisesRegex(ValueError, "nonfinite triangle geometry"):
            triangle_gradient(((0.0, 0.0), (1.0e308, 0.0), (0.0, 1.0e308)),
                              (0.0, 1.0, 0.0))
    def test_area_weighted_recovery_preserves_support_identity(self):
        coordinates = {
            "0": (0.0, 0.0),
            "1": (2.0, 0.0),
            "2": (0.0, 1.0),
            "3": (0.0, -2.0),
        }
        cells = {
            "small": ("0", "1", "2"),
            "large": ("0", "3", "1"),
        }
        vectors = {"small": (1.0, 2.0), "large": (4.0, -1.0)}

        nodes = cell_to_node_vectors(vectors, cells, coordinates)
        self.assertEqual(list(nodes["values"]), ["0", "1", "2", "3"])
        self.assertEqual(nodes["values"]["0"], (3.0, 0.0))
        self.assertEqual(nodes["weights"]["0"], (("large", 2.0 / 3.0), ("small", 1.0 / 3.0)))
        self.assertEqual(nodes["values"]["2"], (1.0, 2.0))

        edges = {"shared": ("0", "1"), "outer": ("1", "2")}
        recovered_edges = cell_to_edge_vectors(vectors, cells, edges, coordinates)
        self.assertEqual(list(recovered_edges["values"]), ["outer", "shared"])
        self.assertEqual(recovered_edges["values"]["shared"], (3.0, 0.0))
        self.assertEqual(recovered_edges["weights"]["shared"],
                         (("large", 2.0 / 3.0), ("small", 1.0 / 3.0)))

    def test_mirror_and_directed_edge_difference_have_deterministic_signs(self):
        self.assertEqual(mirror_vector((3.0, -4.0)), (-3.0, -4.0))
        self.assertEqual(edge_scalar_difference(2.0, 8.0, (0.0, 0.0), (2.0, 0.0)), 3.0)
        self.assertEqual(edge_scalar_difference(8.0, 2.0, (2.0, 0.0), (0.0, 0.0)), -3.0)
        with self.assertRaisesRegex(ValueError, "zero length"):
            edge_scalar_difference(2.0, 8.0, (0.0, 0.0), (0.0, 0.0))

        overflowing_edges = (
            ((-1.0e308, 0.0), (1.0e308, 0.0)),
            ((0.0, 0.0), (1.3e308, 1.3e308)),
        )
        for start, end in overflowing_edges:
            with self.subTest(start=start, end=end):
                with self.assertRaises(ValueError):
                    edge_scalar_difference(2.0, 8.0, start, end)

    def test_vector_error_types_zero_floor_nonfinite_and_undefined_direction(self):
        exact = vector_error((-3.0, 4.0), (-3.0, 4.0), reference_floor=1.0e-12)
        self.assertEqual(exact.magnitude_status, SampleStatus.VALID)
        self.assertEqual(exact.direction_status, SampleStatus.VALID)
        self.assertEqual(exact.relative_magnitude_error, 0.0)
        self.assertEqual(exact.angle_deg, 0.0)

        zero_reference = vector_error((1.0, 0.0), (0.0, 0.0), reference_floor=1.0e-12)
        self.assertEqual(zero_reference.magnitude_status, SampleStatus.GEOMETRIC_ZERO)
        self.assertEqual(zero_reference.direction_status, SampleStatus.DIRECTION_UNDEFINED)
        self.assertIsNone(zero_reference.relative_magnitude_error)
        self.assertIsNone(zero_reference.angle_deg)

        below_floor = vector_error((1.0, 0.0), (1.0e-15, 0.0), reference_floor=1.0e-12)
        self.assertEqual(below_floor.magnitude_status, SampleStatus.BELOW_FLOOR)
        self.assertEqual(below_floor.direction_status, SampleStatus.DIRECTION_UNDEFINED)

        at_floor = vector_error((1.0, 0.0), (1.0e-12, 0.0), reference_floor=1.0e-12)
        self.assertEqual(at_floor.magnitude_status, SampleStatus.BELOW_FLOOR)
        self.assertEqual(at_floor.direction_status, SampleStatus.DIRECTION_UNDEFINED)

        zero_candidate = vector_error((0.0, 0.0), (1.0, 0.0), reference_floor=1.0e-12)
        self.assertEqual(zero_candidate.magnitude_status, SampleStatus.VALID)
        self.assertEqual(zero_candidate.relative_magnitude_error, 1.0)
        self.assertEqual(zero_candidate.direction_status, SampleStatus.DIRECTION_UNDEFINED)

        nonfinite = vector_error((math.inf, 0.0), (1.0, 0.0), reference_floor=1.0e-12)
        self.assertEqual(nonfinite.magnitude_status, SampleStatus.NONFINITE)
        self.assertEqual(nonfinite.direction_status, SampleStatus.NONFINITE)

    def test_vector_error_handles_extreme_finite_vectors_safely(self):
        huge = vector_error(
            (1.0e308, 1.0e308), (1.0e308, -1.0e308), reference_floor=0.0
        )
        self.assertIs(huge.magnitude_status, SampleStatus.VALID)
        self.assertIs(huge.direction_status, SampleStatus.VALID)
        self.assertEqual(huge.relative_magnitude_error, 0.0)
        self.assertEqual(huge.angle_deg, 90.0)

        subnormal = math.nextafter(0.0, 1.0)
        tiny = vector_error(
            (subnormal, 0.0), (subnormal, 0.0), reference_floor=0.0
        )
        self.assertIs(tiny.magnitude_status, SampleStatus.VALID)
        self.assertIs(tiny.direction_status, SampleStatus.VALID)
        self.assertEqual(tiny.relative_magnitude_error, 0.0)
        self.assertEqual(tiny.angle_deg, 0.0)

    def test_evaluator_emits_all_candidates_and_applies_fixed_field_thresholds(self):
        observations = self.make_observations(reference_scale=1.0)
        mesh = {
            "triangles": {"cell-0": ("0", "1", "2")},
            "edges": {
                "edge-01": ("0", "1"),
                "edge-02": ("0", "2"),
                "edge-12": ("1", "2"),
            },
        }
        results = evaluate_field_candidates(observations, mesh, reference_floor=1.0e-12)
        self.assertEqual([item.candidate for item in results], [
            "edge_area_weighted_minus_grad_psi",
            "node_area_weighted_minus_grad_psi",
            "signed_edge_minus_delta_psi_over_h",
            "triangle_minus_grad_psi",
        ])
        self.assertTrue(all(item.classification is Identifiability.IDENTIFIED for item in results))
        self.assertTrue(all(item.median_relative_magnitude_error == 0.0 for item in results))
        self.assertTrue(all(item.median_angle_deg == 0.0 for item in results))
        triangle = next(item for item in results if item.candidate == "triangle_minus_grad_psi")
        self.assertEqual(triangle.support_kind, SupportKind.CELL)
        self.assertEqual(triangle.samples[0].candidate_value, (-3.0, 4.0))
        self.assertEqual(triangle.samples[0].reference_value, (-3.0, 4.0))

        outside_gate = evaluate_field_candidates(
            self.make_observations(reference_scale=0.97), mesh, reference_floor=1.0e-12,
            thresholds=AcceptanceThresholds(),
        )
        self.assertTrue(all(item.classification is Identifiability.REJECTED
                            for item in outside_gate))
        self.assertTrue(all(item.median_relative_magnitude_error > 0.02
                            for item in outside_gate))
        self.assertTrue(all(item.median_angle_deg <= 1.0 for item in outside_gate))
        magnitude_boundary = evaluate_field_candidates(
            self.make_observations(reference_scale=math.nextafter(1.0 / 1.02, 1.0)), mesh,
            reference_floor=1.0e-12,
        )
        magnitude_triangle = next(
            item for item in magnitude_boundary if item.support_kind is SupportKind.CELL
        )
        self.assertAlmostEqual(magnitude_triangle.median_relative_magnitude_error, 0.02, places=14)
        self.assertIs(magnitude_triangle.classification, Identifiability.IDENTIFIED)

        magnitude_above = evaluate_field_candidates(
            self.make_observations(reference_scale=1.0 / 1.0200000001), mesh,
            reference_floor=1.0e-12,
        )
        magnitude_above_triangle = next(
            item for item in magnitude_above if item.support_kind is SupportKind.CELL
        )
        self.assertGreater(magnitude_above_triangle.median_relative_magnitude_error, 0.02)
        self.assertIs(magnitude_above_triangle.classification, Identifiability.REJECTED)

        direction_boundary = evaluate_field_candidates(
            self.make_observations(reference_scale=1.0, reference_angle_deg=1.0), mesh,
            reference_floor=1.0e-12,
        )
        direction_triangle = next(
            item for item in direction_boundary if item.support_kind is SupportKind.CELL
        )
        self.assertAlmostEqual(direction_triangle.median_angle_deg, 1.0, places=12)
        self.assertIs(direction_triangle.classification, Identifiability.IDENTIFIED)

        direction_above = evaluate_field_candidates(
            self.make_observations(reference_scale=1.0, reference_angle_deg=1.0000001), mesh,
            reference_floor=1.0e-12,
        )
        direction_above_triangle = next(
            item for item in direction_above if item.support_kind is SupportKind.CELL
        )
        self.assertGreater(direction_above_triangle.median_angle_deg, 1.0)
        self.assertIs(direction_above_triangle.classification, Identifiability.REJECTED)

        with self.assertRaises(ValueError):
            evaluate_field_candidates(
                observations, mesh, reference_floor=1.0e-12,
                thresholds=AcceptanceThresholds(field_median_relative=0.03),
            )

    @staticmethod
    def make_observations(*, reference_scale, reference_angle_deg=0.0):
        points = {"0": (0.0, 0.0), "1": (2.0, 0.0), "2": (0.0, 1.0)}
        theta = math.radians(reference_angle_deg)
        reference_x = reference_scale * (-3.0 * math.cos(theta) - 4.0 * math.sin(theta))
        reference_y = reference_scale * (-3.0 * math.sin(theta) + 4.0 * math.cos(theta))
        rows = []

        def row(node, quantity, component, value, unit):
            return Observation(
                solver="sentaurus", topology="sketch", bias_V=-12.0,
                support_kind=SupportKind.NODE, support_id=node, quantity=quantity,
                component=component, raw_value=value, raw_unit=unit, value_si=value,
                unit_si=unit, coordinate_frame="canonical_xy", orientation="global_xy",
                conversion="identity", status=SampleStatus.VALID, source_path="state.csv",
                source_sha256="0" * 64,
            )

        for node, (x, y) in points.items():
            rows.extend((
                row(node, "coordinate", "x", x, "m"),
                row(node, "coordinate", "y", y, "m"),
                row(node, "ElectrostaticPotential", "component0", 3.0 * x - 4.0 * y + 2.0, "V"),
                row(node, "ElectricField", "component0", reference_x, "V/m"),
                row(node, "ElectricField", "component1", reference_y, "V/m"),
            ))
        return tuple(rows)


if __name__ == "__main__":
    unittest.main()
