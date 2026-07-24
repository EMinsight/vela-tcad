import unittest

from scripts.pn2d_minimal6_diagnostics.factorial_attribution import (
    interaction_remainder,
    replacement_order_contributions,
    shapley_contributions,
)


class Minimal6FactorialAttributionTest(unittest.TestCase):
    def test_three_factor_shapley_closes_for_interacting_lattice(self) -> None:
        factors = ("low_field", "drive", "support")
        values = {}
        for mask in range(8):
            active = frozenset(
                factor
                for bit, factor in enumerate(factors)
                if mask & (1 << bit)
            )
            values[active] = (
                0.2
                + 0.1 * ("low_field" in active)
                + 0.3 * ("drive" in active)
                - 0.2 * ("support" in active)
                + 0.4
                * (
                    "drive" in active
                    and "support" in active
                )
            )
        orders = replacement_order_contributions(values)
        shapley = shapley_contributions(values)
        closure = (
            values[frozenset(factors)]
            - values[frozenset()]
            - sum(shapley.values())
        )
        self.assertEqual(len(orders), 18)
        self.assertAlmostEqual(closure, 0.0)
        self.assertAlmostEqual(interaction_remainder(values), 0.4)

    def test_incomplete_lattice_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "incomplete"):
            shapley_contributions({frozenset(): 0.0})


if __name__ == "__main__":
    unittest.main()
