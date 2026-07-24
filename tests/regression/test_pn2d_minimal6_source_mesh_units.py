import unittest

from scripts.pn2d_minimal6_diagnostics.counterfactual import (
    integrate_native_nodal_per_unit_depth,
)


class NativeSourceMeshUnitTests(unittest.TestCase):
    def test_meter_and_micrometer_meshes_integrate_identically(self) -> None:
        meter_mesh = {
            "coordinate_unit": "m",
            "nodes": [
                {"id": 0, "x": 0.0, "y": 0.0},
                {"id": 1, "x": 1.0e-6, "y": 0.0},
                {"id": 2, "x": 0.0, "y": 1.0e-6},
            ],
            "triangles": [{"node_ids": [0, 1, 2]}],
        }
        micrometer_mesh = {
            "coordinate_unit": "um",
            "nodes": [
                {"id": 0, "x": 0.0, "y": 0.0},
                {"id": 1, "x": 1.0, "y": 0.0},
                {"id": 2, "x": 0.0, "y": 1.0},
            ],
            "triangles": [{"node_ids": [0, 1, 2]}],
        }
        values = {0: 3.0, 1: 6.0, 2: 9.0}

        meter = integrate_native_nodal_per_unit_depth(meter_mesh, values)
        micrometer = integrate_native_nodal_per_unit_depth(
            micrometer_mesh, values
        )

        self.assertAlmostEqual(
            meter["value_s_inv_per_unit_depth"], 3.0e-8
        )
        self.assertEqual(
            meter["value_s_inv_per_unit_depth"],
            micrometer["value_s_inv_per_unit_depth"],
        )


if __name__ == "__main__":
    unittest.main()
