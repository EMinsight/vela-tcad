import unittest

from scripts.pn2d_minimal6_diagnostics.phase_e_continuity_residual import (
    _node_divergence,
    _norm,
)


class PhaseEContinuityResidualTest(unittest.TestCase):
    def test_directed_edge_flux_is_added_to_node0_and_subtracted_from_node1(self):
        edges = [
            {
                "edge_id": "0",
                "node0": "0",
                "node1": "1",
                "electron_flux": "2.0",
            },
            {
                "edge_id": "1",
                "node0": "1",
                "node1": "2",
                "electron_flux": "-3.0",
            },
        ]
        divergence, absolute = _node_divergence(edges, "electron")
        self.assertEqual(divergence, {0: 2.0, 1: -5.0, 2: 3.0})
        self.assertEqual(absolute, {0: 2.0, 1: 5.0, 2: 3.0})

    def test_edge_mobility_replay_is_linear_in_mobility(self):
        edges = [
            {
                "edge_id": "7",
                "node0": "1",
                "node1": "5",
                "electron_flux": "4.0",
            }
        ]
        divergence, _ = _node_divergence(edges, "electron", {7: 0.25})
        self.assertEqual(divergence, {1: 1.0, 5: -1.0})

    def test_term_balance_normalization_uses_absolute_component_sum(self):
        self.assertAlmostEqual(_norm(2.0, [4.0, -1.0, -1.0]), 1.0 / 3.0)


if __name__ == "__main__":
    unittest.main()
