import tempfile
import unittest
from pathlib import Path
from scripts.pn2d_minimal6_diagnostics.geometry import p1_gradient
from scripts.pn2d_minimal6_diagnostics.support import project_vector_to_edge, integrate_nodal_field, integrate_cell_field, map_local_sources_to_nodes, node_scalar_to_cells, node_vector_to_edges, edge_scalar_to_cells, local_edge_sources_to_nodes
from scripts.pn2d_minimal6_diagnostics.physics import van_overstraeten_alpha, invert_alpha, infer_ni_eff, parse_van_overstraeten_de_man, invert_piecewise_alpha, compare_van_overstraeten_parameters

class DiagnosticPhysicsTest(unittest.TestCase):
    def test_p1_gradient_recovers_affine_field_in_ccw_triangle(self):
        # f(x,y)=2x-3y+5 on a CCW triangle.
        gradient = p1_gradient(((0.,0.), (2.,0.), (0.,1.)), (5.,9.,2.))
        self.assertEqual(gradient, (2., -3.))
    def test_p1_gradient_is_affine_on_both_canonical_topologies(self):
        coordinates = {0:(0.,0.), 1:(1.,0.), 2:(2.,0.), 3:(0.,1.), 4:(1.,1.), 5:(2.,1.)}
        for triangles in (((0, 1, 4), (0, 4, 3), (1, 2, 5), (1, 5, 4)),
                          ((0, 1, 3), (1, 4, 3), (1, 2, 4), (2, 5, 4))):
            for triangle in triangles:
                values = [2. * coordinates[node][0] - 3. * coordinates[node][1] + 5. for node in triangle]
                self.assertEqual(p1_gradient([coordinates[node] for node in triangle], values), (2., -3.))
    def test_p1_gradient_rejects_clockwise_triangle(self):
        with self.assertRaises(ValueError):
            p1_gradient(((0.,0.), (0.,1.), (2.,0.)), (5.,2.,9.))
    def test_van_overstraeten_switch_and_inverse_domain(self):
        low = van_overstraeten_alpha(1.0e5, 1.0e6, 2.0e5, 2.0e6, 4.0e5, 3.0e5)
        high = van_overstraeten_alpha(5.0e5, 1.0e6, 2.0e5, 2.0e6, 4.0e5, 3.0e5)
        self.assertGreater(high, low)
        self.assertAlmostEqual(invert_alpha(low, 1.0e6, 2.0e5), 1.0e5)
        with self.assertRaises(ValueError): invert_alpha(0.0, 1.0e6, 2.0e5)
    def test_piecewise_alpha_inversion_returns_all_valid_branches(self):
        candidates = invert_piecewise_alpha(3.5e5, low_a_cm_inv=1.e6, low_b_v_per_cm=2.e5, high_a_cm_inv=4.e5, high_b_v_per_cm=1.e5, switch_v_per_cm=3.e5)
        self.assertEqual([item["branch"] for item in candidates], ["low", "high"])
        self.assertLess(candidates[0]["field_v_per_cm"], 3.e5)
        self.assertGreaterEqual(candidates[1]["field_v_per_cm"], 3.e5)
    def test_parses_versioned_van_overstraeten_parameters(self):
        text = """vanOverstraetendeMan * Impact Ionization:
{
    a(low) = 7.03e5 , 1.582e6
    a(high) = 7.03e5 , 6.71e5
    b(low) = 1.231e6 , 2.036e6
    b(high) = 1.231e6 , 1.693e6
    E0 = 4.0e5 , 4.0e5
    hbarOmega = 0.063 , 0.063
}
"""
        with tempfile.TemporaryDirectory() as temp:
            path = Path(temp) / "models.par"
            path.write_text(text, encoding="utf-8")
            parsed = parse_van_overstraeten_de_man(path)
        self.assertEqual(parsed["electron"]["a_low_cm_inv"], 7.03e5)
        self.assertEqual(parsed["hole"]["b_high_v_per_cm"], 1.693e6)
        self.assertEqual(parsed["switch_field_v_per_cm"], 4.0e5)
        self.assertEqual(len(parsed["sha256"]), 64)
    def test_parameter_agreement_requires_numeric_vela_coefficients(self):
        parsed = {
            "electron": {"a_low_cm_inv": 1., "a_high_cm_inv": 2., "b_low_v_per_cm": 3., "b_high_v_per_cm": 4.},
            "hole": {"a_low_cm_inv": 5., "a_high_cm_inv": 6., "b_low_v_per_cm": 7., "b_high_v_per_cm": 8.},
            "switch_field_v_per_cm": 9.,
        }
        self.assertEqual(compare_van_overstraeten_parameters(parsed, None)["status"], "unavailable")
        production = {carrier: dict(values) for carrier, values in parsed.items() if carrier in ("electron", "hole")}
        production["switch_field_v_per_cm"] = 9.
        self.assertEqual(compare_van_overstraeten_parameters(parsed, production)["status"], "available")
        production["electron"]["a_low_cm_inv"] = 1.1
        self.assertEqual(compare_van_overstraeten_parameters(parsed, production)["status"], "mismatch")
    def test_integral_and_node_mapping_are_conservative(self):
        self.assertEqual(integrate_nodal_field((2., 4., 6.), (1., 2., 3.)), 28.)
        mapped = map_local_sources_to_nodes(((0, 1, 2),), (12.,))
        self.assertEqual(mapped, {0: 4., 1: 4., 2: 4.})
    def test_named_support_conversions_return_normalized_weights(self):
        cells = node_scalar_to_cells({0:3., 1:6., 2:9.}, ((0,1,2),))
        self.assertEqual(cells["values"], [6.])
        self.assertEqual(cells["weights"], [{0:1./3., 1:1./3., 2:1./3.}])
        with self.assertRaises(ValueError):
            node_scalar_to_cells({0:3., 1:6., 2:9.}, ((0,1,2),), quantity="ImpactIonization")
        edges = node_vector_to_edges({0:(2.,0.), 1:(4.,0.)}, ((0,1),), {0:(0.,0.), 1:(1.,0.)})
        self.assertEqual(edges["values"], [3.])
        self.assertEqual(edges["weights"], [{0:0.5, 1:0.5}])
        nodes = local_edge_sources_to_nodes(((0,1), (1,2), (2,0)), (3., 3., 3.))
        self.assertEqual(nodes["values"], {0:3., 1:3., 2:3.})
        self.assertEqual(nodes["weights"], [{0:0.5, 1:0.5}, {1:0.5, 2:0.5}, {2:0.5, 0:0.5}])
    def test_edge_to_cell_conversion_is_weighted_and_conservative(self):
        converted = edge_scalar_to_cells({10:3., 11:6., 12:9.}, ((10,11,12),))
        self.assertEqual(converted["values"], [6.])
        self.assertEqual(converted["weights"], [{10:1./3., 11:1./3., 12:1./3.}])
    def test_infer_ni_eff_reports_electron_hole_consistency(self):
        result = infer_ni_eff(psi_V=0.2, phin_V=0.3, phip_V=0.1, n_cm3=1.0e12, p_cm3=1.0e12, thermal_voltage_V=0.1)
        self.assertAlmostEqual(result["electron_cm3"] / 1.0e12, 1.0 / 2.718281828459045)
        self.assertAlmostEqual(result["hole_cm3"] / 1.0e12, 1.0 / 2.718281828459045)
        self.assertAlmostEqual(result["relative_residual"], 0.0)
    def test_cell_field_integrates_triangle_average_times_area(self):
        self.assertAlmostEqual(integrate_cell_field(((0.,0.),(2.,0.),(0.,1.)), (3.,6.,9.)), 6.0)
        self.assertAlmostEqual(integrate_cell_field(((0.,0.),(2.,0.),(0.,1.)), (3.,6.,9.), partial_volume_fraction=0.25), 1.5)
    def test_signed_vector_projection_obeys_edge_direction(self):
        self.assertEqual(project_vector_to_edge((3.,4.), (0.,0.), (1.,0.)), 3.)
        self.assertEqual(project_vector_to_edge((3.,4.), (1.,0.), (0.,0.)), -3.)

if __name__ == '__main__': unittest.main()
