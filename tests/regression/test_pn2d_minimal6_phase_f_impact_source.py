import unittest

from scripts.pn2d_minimal6_diagnostics.phase_f_self_consistent import (
    _triangle_source_per_cm_s,
)


class PhaseFImpactSourceTests(unittest.TestCase):
    def test_triangle_source_uses_physical_per_meter_to_per_cm_conversion(
        self,
    ) -> None:
        row = {}
        for local, (volume, electron, hole) in enumerate(
            (
                (1.0e-13, 2.0, 3.0),
                (1.0e-13, 5.0, 7.0),
                (0.0, 0.0, 0.0),
            )
        ):
            prefix = f"local_edge{local}_"
            row[prefix + "truncated_partial_volume_m2"] = str(volume)
            row[prefix + "electron_source_integral_per_m_s"] = str(electron)
            row[prefix + "hole_source_integral_per_m_s"] = str(hole)

        self.assertEqual(_triangle_source_per_cm_s([row]), 0.17)

    def test_triangle_source_rejects_nonzero_geometric_zero(self) -> None:
        row = {}
        for local in range(3):
            prefix = f"local_edge{local}_"
            row[prefix + "truncated_partial_volume_m2"] = "0"
            row[prefix + "electron_source_integral_per_m_s"] = (
                "1" if local == 2 else "0"
            )
            row[prefix + "hole_source_integral_per_m_s"] = "0"

        with self.assertRaisesRegex(ValueError, "geometric-zero"):
            _triangle_source_per_cm_s([row])


if __name__ == "__main__":
    unittest.main()
