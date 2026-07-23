import unittest

from scripts.verify_pn2d_minimal6_transport_element_closure import close


class Minimal6TransportElementVerifyTest(unittest.TestCase):
    def test_relative_comparison_scales_with_magnitude(self):
        self.assertTrue(close(1.0e20, 1.0e20 * (1.0 + 1.0e-13)))
        self.assertFalse(close(1.0e20, 1.0e20 * (1.0 + 1.0e-8)))

    def test_absolute_comparison_handles_zero(self):
        self.assertTrue(close(0.0, 1.0e-13))
        self.assertFalse(close(0.0, 1.0e-4))


if __name__ == "__main__":
    unittest.main()
