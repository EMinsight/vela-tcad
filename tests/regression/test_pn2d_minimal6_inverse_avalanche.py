import math
import unittest
from dataclasses import replace

from scripts.pn2d_minimal6_diagnostics.inverse_avalanche import (
    GenerationError,
    current_aligned_magnitude,
    evaluate_avalanche_candidates,
    impact_generation,
    invert_van_overstraeten_alpha,
    reconstruct_generation_supports,
    summarize_generation_errors,
)
from scripts.pn2d_minimal6_diagnostics.inverse_contracts import (
    AcceptanceThresholds,
    Identifiability,
    Observation,
    SampleStatus,
    SupportKind,
)
from scripts.pn2d_minimal6_diagnostics.physics import van_overstraeten_alpha


class InverseAvalancheTest(unittest.TestCase):
    def test_forward_inverse_identity_on_every_carrier_branch(self):
        switch = 4.0e5
        parameters = {
            "electron": {
                "low": (7.03e5, 1.231e6, 2.0e5),
                "high": (7.03e5, 1.231e6, 8.0e5),
            },
            "hole": {
                "low": (1.582e6, 2.036e6, 2.0e5),
                "high": (6.71e5, 1.693e6, 8.0e5),
            },
        }
        for carrier, branches in parameters.items():
            low_a, low_b, _ = branches["low"]
            high_a, high_b, _ = branches["high"]
            for branch, (prefactor, critical_field, field) in branches.items():
                with self.subTest(carrier=carrier, branch=branch):
                    alpha = van_overstraeten_alpha(
                        field, low_a, low_b, high_a, high_b, switch
                    )
                    recovered, status = invert_van_overstraeten_alpha(
                        alpha,
                        prefactor=prefactor,
                        critical_field=critical_field,
                        gamma=1.0,
                        branch=branch,
                        switch_field=switch,
                    )
                    self.assertIs(status, SampleStatus.VALID)
                    self.assertAlmostEqual(recovered / field, 1.0, delta=1.0e-12)

    def test_temperature_factor_inverse_is_exact(self):
        field = 6.75e7
        prefactor = 7.03e7
        critical_field = 1.231e8
        gamma = 1.17
        alpha = gamma * prefactor * math.exp(-gamma * critical_field / field)
        recovered, status = invert_van_overstraeten_alpha(
            alpha,
            prefactor=prefactor,
            critical_field=critical_field,
            gamma=gamma,
        )
        self.assertIs(status, SampleStatus.VALID)
        self.assertAlmostEqual(recovered / field, 1.0, delta=1.0e-12)

    def test_inverse_fails_closed_with_typed_floor_branch_underflow_and_nonfinite(self):
        arguments = dict(prefactor=10.0, critical_field=20.0, gamma=2.0)
        self.assertEqual(
            invert_van_overstraeten_alpha(0.0, **arguments),
            (None, SampleStatus.BELOW_FLOOR),
        )
        self.assertEqual(
            invert_van_overstraeten_alpha(None, **arguments),
            (None, SampleStatus.MISSING_FIELD),
        )
        self.assertEqual(
            invert_van_overstraeten_alpha(math.inf, **arguments),
            (None, SampleStatus.NONFINITE),
        )
        self.assertEqual(
            invert_van_overstraeten_alpha(20.0, **arguments),
            (None, SampleStatus.BRANCH_AMBIGUOUS),
        )
        self.assertEqual(
            invert_van_overstraeten_alpha(
                5.0e-324,
                prefactor=1.0e308,
                critical_field=1.0,
                gamma=1.0,
            ),
            (None, SampleStatus.EXPONENTIAL_UNDERFLOW),
        )
        self.assertEqual(
            invert_van_overstraeten_alpha(
                1.0e-15, numerical_floor=1.0e-12, **arguments
            ),
            (None, SampleStatus.BELOW_FLOOR),
        )

    def test_inverse_rejects_fields_inconsistent_with_declared_branch(self):
        low_inconsistent_alpha = 10.0 * math.exp(-20.0 / 8.0)
        self.assertEqual(
            invert_van_overstraeten_alpha(
                low_inconsistent_alpha,
                prefactor=10.0,
                critical_field=20.0,
                gamma=1.0,
                branch="low",
                switch_field=5.0,
            ),
            (None, SampleStatus.BRANCH_AMBIGUOUS),
        )
        high_inconsistent_alpha = 10.0 * math.exp(-20.0 / 2.0)
        self.assertEqual(
            invert_van_overstraeten_alpha(
                high_inconsistent_alpha,
                prefactor=10.0,
                critical_field=20.0,
                gamma=1.0,
                branch="high",
                switch_field=5.0,
            ),
            (None, SampleStatus.BRANCH_AMBIGUOUS),
        )

    def test_current_projection_uses_magnitude_and_never_invents_zero_direction(self):
        self.assertEqual(current_aligned_magnitude((3.0, 4.0), (0.0, 24.0)), 4.0)
        self.assertEqual(current_aligned_magnitude((-3.0, -4.0), (0.0, -24.0)), 4.0)
        with self.assertRaisesRegex(ValueError, "direction"):
            current_aligned_magnitude((3.0, 4.0), (0.0, 0.0))
        with self.assertRaisesRegex(ValueError, "finite"):
            current_aligned_magnitude((math.inf, 0.0), (1.0, 0.0))

    def test_impact_generation_is_sign_independent_and_has_si_source_dimensions(self):
        self.assertEqual(
            impact_generation(2.0, (3.0, 4.0), 7.0, (0.0, 24.0), q=2.0),
            89.0,
        )
        self.assertEqual(
            impact_generation(2.0, (-3.0, -4.0), 7.0, (0.0, -24.0), q=2.0),
            89.0,
        )
        # (m^-1 * A*m^-2) / C = m^-3*s^-1.
        value = impact_generation(
            2.0e5,
            (3.0e4, 4.0e4),
            7.0e5,
            (0.0, 2.4e5),
            q=1.602176634e-19,
        )
        # (2e5*5e4 + 7e5*2.4e5) / q = 1.78e11 / q.
        self.assertAlmostEqual(value, 1.78e11 / 1.602176634e-19)
        with self.assertRaisesRegex(ValueError, "non-negative"):
            impact_generation(-1.0, (1.0, 0.0), 1.0, (1.0, 0.0))

    def test_two_triangle_reconstruction_keeps_local_integrated_and_depth_units_distinct(self):
        # A 1 cm by 1 cm square expressed in canonical metres.
        coordinates_m = {
            "0": (0.0, 0.0),
            "1": (0.01, 0.0),
            "2": (0.01, 0.01),
            "3": (0.0, 0.01),
        }
        triangles = {
            "a": ("0", "1", "2"),
            "b": ("0", "2", "3"),
        }
        reconstruction = reconstruct_generation_supports(
            coordinates_m=coordinates_m,
            triangles=triangles,
            native_nodal_generation_m3_s={node: 6.0e6 for node in coordinates_m},
            candidate_cell_generation_m3_s={"a": 6.0e6, "b": 6.0e6},
            edges={"e01": ("0", "1"), "e23": ("2", "3")},
            vela_edge_partial_sources_per_m_s={"e01": 2.0, "e23": 4.0},
            depth_m=0.01,
        )

        self.assertEqual(reconstruction.native_nodal_unit, "m^-3*s^-1")
        self.assertEqual(reconstruction.cell_integral_unit, "m^-1*s^-1")
        self.assertEqual(reconstruction.depth_integral_unit, "s^-1")
        self.assertEqual(reconstruction.native_cell_integrals_per_m_s, {"a": 300.0, "b": 300.0})
        self.assertEqual(reconstruction.candidate_cell_integrals_per_m_s, {"a": 300.0, "b": 300.0})
        self.assertEqual(reconstruction.native_integrated_per_m_s, 600.0)
        self.assertEqual(reconstruction.candidate_integrated_per_m_s, 600.0)
        self.assertEqual(reconstruction.native_one_cm_depth_s_inv, 6.0)
        self.assertEqual(reconstruction.candidate_one_cm_depth_s_inv, 6.0)
        self.assertEqual(
            reconstruction.candidate_node_mapped_per_m_s,
            {"0": 200.0, "1": 100.0, "2": 200.0, "3": 100.0},
        )
        self.assertEqual(
            reconstruction.vela_node_mapped_per_m_s,
            {"0": 1.0, "1": 1.0, "2": 2.0, "3": 2.0},
        )
        self.assertEqual(sum(reconstruction.candidate_node_mapped_per_m_s.values()), 600.0)
        self.assertEqual(sum(reconstruction.vela_node_mapped_per_m_s.values()), 6.0)

    def test_generation_gates_include_boundaries_and_reject_overrides(self):
        boundary = summarize_generation_errors(
            (
                GenerationError(SampleStatus.VALID, 0.0),
                GenerationError(SampleStatus.VALID, 0.3),
            ),
            GenerationError(SampleStatus.VALID, 0.1),
        )
        self.assertEqual(boundary.local_median_abs_log10_error, 0.15)
        self.assertEqual(boundary.local_max_abs_log10_error, 0.3)
        self.assertEqual(boundary.integrated_abs_log10_error, 0.1)
        self.assertIs(boundary.classification, Identifiability.IDENTIFIED)

        local_failure = summarize_generation_errors(
            (GenerationError(SampleStatus.VALID, 0.3000001),),
            GenerationError(SampleStatus.VALID, 0.0),
        )
        self.assertIs(local_failure.classification, Identifiability.REJECTED)
        integrated_failure = summarize_generation_errors(
            (GenerationError(SampleStatus.VALID, 0.0),),
            GenerationError(SampleStatus.VALID, 0.1000001),
        )
        self.assertIs(integrated_failure.classification, Identifiability.REJECTED)
        unavailable = summarize_generation_errors(
            (GenerationError(SampleStatus.BELOW_FLOOR, None),),
            GenerationError(SampleStatus.BELOW_FLOOR, None),
        )
        self.assertIs(unavailable.classification, Identifiability.INSUFFICIENT_DATA)
        with self.assertRaisesRegex(ValueError, "immutable"):
            summarize_generation_errors(
                (GenerationError(SampleStatus.VALID, 0.0),),
                GenerationError(SampleStatus.VALID, 0.0),
                thresholds=AcceptanceThresholds(local_generation_abs_dex=0.4),
            )

    def test_evaluator_emits_fixed_candidates_splits_and_deterministic_generation_records(self):
        observations = self.make_avalanche_observations()
        kwargs = dict(
            mesh=self.mesh(),
            parameters=self.parameters(),
            generation_floor=1.0e-30,
            current_floor=1.0e-30,
            reference_densities_m3={"electron": 2.0, "hole": 2.0},
            q=1.0,
            vela_edge_partial_sources_per_state={
                ("sentaurus", "sketch", -12.0): {"edge-01": 2.0},
            },
        )
        results = evaluate_avalanche_candidates(observations, **kwargs)
        expected_candidates = {
            "electric_field_magnitude",
            "qf_gradient_magnitude",
            "electric_field_current_aligned",
            "qf_gradient_current_aligned",
            "density_interpolated_qf_electric",
        }
        self.assertEqual({item.candidate for item in results}, expected_candidates)
        self.assertEqual({item.split for item in results}, {"discovery", "holdout"})
        self.assertTrue(all(item.summary.classification is Identifiability.IDENTIFIED
                            for item in results))
        self.assertTrue(all(item.summary.local_max_abs_log10_error == 0.0
                            for item in results))
        self.assertTrue(all(item.summary.integrated_abs_log10_error == 0.0
                            for item in results))

        discovery = next(
            item for item in results
            if item.candidate == "electric_field_magnitude"
            and item.split == "discovery"
        )
        alpha = 10.0 * math.exp(-4.0)
        generation = alpha * 15.0
        self.assertAlmostEqual(discovery.samples[0].candidate_generation_m3_s, generation)
        self.assertAlmostEqual(discovery.supports.native_integrated_per_m_s,
                               0.5 * generation)
        self.assertAlmostEqual(discovery.supports.candidate_integrated_per_m_s,
                               0.5 * generation)
        self.assertEqual(discovery.supports.depth_m, 0.01)
        self.assertEqual(discovery.supports.vela_edge_partial_sources_per_m_s,
                         {"edge-01": 2.0})
        self.assertEqual(discovery.supports.vela_node_mapped_per_m_s,
                         {"0": 1.0, "1": 1.0})
        self.assertEqual(
            results,
            evaluate_avalanche_candidates(tuple(reversed(observations)), **kwargs),
        )

    def test_evaluator_never_uses_native_generation_to_construct_candidate(self):
        kwargs = dict(
            mesh=self.mesh(), parameters=self.parameters(), generation_floor=1.0e-30,
            current_floor=1.0e-30, q=1.0,
        )
        baseline = evaluate_avalanche_candidates(
            self.make_avalanche_observations(states=(("sketch", -12.0),)), **kwargs
        )
        changed_rows = tuple(
            replace(row, raw_value=row.raw_value * 10.0, value_si=row.value_si * 10.0)
            if row.quantity == "ImpactIonization" else row
            for row in self.make_avalanche_observations(states=(("sketch", -12.0),))
        )
        changed = evaluate_avalanche_candidates(changed_rows, **kwargs)
        self.assertEqual(
            [sample.candidate_generation_m3_s for item in baseline for sample in item.samples],
            [sample.candidate_generation_m3_s for item in changed for sample in item.samples],
        )
        self.assertTrue(all(item.summary.classification is Identifiability.REJECTED
                            for item in changed))
        self.assertTrue(all(item.summary.local_max_abs_log10_error == 1.0
                            for item in changed))

    def test_evaluator_types_missing_qf_and_omits_undeclared_density_interpolation(self):
        observations = tuple(
            replace(row, raw_value=None, value_si=None, status=SampleStatus.MISSING_FIELD)
            if row.quantity in {"eQuasiFermiPotential", "hQuasiFermiPotential"}
            else row
            for row in self.make_avalanche_observations(states=(("sketch", -12.0),))
        )
        results = evaluate_avalanche_candidates(
            observations, self.mesh(), parameters=self.parameters(),
            generation_floor=1.0e-30, current_floor=1.0e-30, q=1.0,
        )
        self.assertNotIn("density_interpolated_qf_electric",
                         {item.candidate for item in results})
        electric = next(item for item in results
                        if item.candidate == "electric_field_magnitude")
        self.assertIs(electric.summary.classification, Identifiability.IDENTIFIED)
        qf = next(item for item in results
                  if item.candidate == "qf_gradient_magnitude")
        self.assertIs(qf.summary.classification, Identifiability.INSUFFICIENT_DATA)
        self.assertTrue(qf.exclusions)
        self.assertTrue(all(record.status is SampleStatus.MISSING_FIELD
                            for record in qf.exclusions))

    @staticmethod
    def mesh():
        return {
            "triangles": {"cell-0": ("0", "1", "2")},
            "edges": {"edge-01": ("0", "1")},
        }

    @staticmethod
    def parameters():
        return {
            "gamma": 1.0,
            "switch_field_V_m": 10.0,
            "electron": {"low": (10.0, 20.0), "high": (8.0, 16.0)},
            "hole": {"low": (10.0, 20.0), "high": (8.0, 16.0)},
        }

    @staticmethod
    def make_avalanche_observations(
        *, states=(("sketch", -12.0), ("mirror", -12.0))
    ):
        points = {"0": (0.0, 0.0), "1": (1.0, 0.0), "2": (0.0, 1.0)}
        alpha = 10.0 * math.exp(-4.0)
        generation = alpha * 15.0
        rows = []

        def row(topology, bias, node, quantity, component, value, unit):
            return Observation(
                solver="sentaurus", topology=topology, bias_V=bias,
                support_kind=SupportKind.NODE, support_id=node, quantity=quantity,
                component=component, raw_value=value, raw_unit=unit, value_si=value,
                unit_si=unit, coordinate_frame="canonical_xy", orientation="global_xy",
                conversion="identity", status=SampleStatus.VALID,
                source_path="state.csv", source_sha256="0" * 64,
            )

        for topology, bias in states:
            for node, (x, y) in points.items():
                qf = 3.0 * x + 4.0 * y
                rows.extend((
                    row(topology, bias, node, "coordinate", "x", x, "m"),
                    row(topology, bias, node, "coordinate", "y", y, "m"),
                    row(topology, bias, node, "ElectricField", "component0", 3.0, "V/m"),
                    row(topology, bias, node, "ElectricField", "component1", 4.0, "V/m"),
                    row(topology, bias, node, "eQuasiFermiPotential", "component0", qf, "V"),
                    row(topology, bias, node, "hQuasiFermiPotential", "component0", qf, "V"),
                    row(topology, bias, node, "eCurrentDensity", "component0", 3.0, "A/m^2"),
                    row(topology, bias, node, "eCurrentDensity", "component1", 4.0, "A/m^2"),
                    row(topology, bias, node, "hCurrentDensity", "component0", 6.0, "A/m^2"),
                    row(topology, bias, node, "hCurrentDensity", "component1", 8.0, "A/m^2"),
                    row(topology, bias, node, "eDensity", "component0", 2.0, "m^-3"),
                    row(topology, bias, node, "hDensity", "component0", 2.0, "m^-3"),
                    row(topology, bias, node, "eAlphaAvalanche", "component0", alpha, "m^-1"),
                    row(topology, bias, node, "hAlphaAvalanche", "component0", alpha, "m^-1"),
                    row(topology, bias, node, "ImpactIonization", "component0", generation,
                        "m^-3*s^-1"),
                ))
        return tuple(rows)


if __name__ == "__main__":
    unittest.main()
