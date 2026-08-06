"""Regression coverage for BVmethods native-to-SI postprocessing units."""

from __future__ import annotations

import unittest

from scripts.analyze_bvmethods_nmos_iic_postprocess import (
    Q_C,
    avalanche_current_from_native_source,
    current_density_from_native_particle_flux,
)
from scripts.compare_bvmethods_nmos_iic_multibias_fields import vela_edge_value


class BVMethodsNmosIicUnitTests(unittest.TestCase):
    def test_native_source_integral_to_current_per_um_has_two_length_factors(self) -> None:
        # alpha: cm^-1 -> m^-1, flux: cm^-2 -> m^-2, area: um^2 -> m^2
        # combine to 1e-6; per-m device depth -> per-um contributes another 1e-6.
        self.assertAlmostEqual(
            avalanche_current_from_native_source(1.0),
            Q_C * 1.0e-12,
        )

    def test_native_particle_flux_to_current_density_uses_cm2_to_m2(self) -> None:
        self.assertAlmostEqual(
            current_density_from_native_particle_flux(1.0),
            Q_C * 1.0e4,
        )

    def test_edge_generation_converts_native_line_source_before_dividing_by_si_area(self) -> None:
        row = {"edge_source_integral": "8.0", "edge_area_proxy_m2": "2e-12"}
        self.assertAlmostEqual(vela_edge_value(row, "avalanche_generation"), 4.0e6)


if __name__ == "__main__":
    unittest.main()
