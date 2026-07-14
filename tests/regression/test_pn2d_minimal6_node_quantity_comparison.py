import importlib.util
from pathlib import Path
import unittest


REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "plot_pn2d_minimal6_node_quantity_comparison.py"


def load_module():
    spec = importlib.util.spec_from_file_location("minimal6_node_compare", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


comparison = load_module()


class Minimal6NodeQuantityComparisonTest(unittest.TestCase):
    def test_edge_quantities_are_aggregated_to_six_nodes(self):
        rows = [
            {
                "node0": "1", "node1": "2",
                "vela_electron_impact_field_V_per_m": "2",
                "vela_hole_impact_field_V_per_m": "4",
                "vela_electron_alpha_per_m": "6",
                "vela_hole_alpha_per_m": "8",
                "vela_electron_flux_per_m2_s": "10",
                "vela_hole_flux_per_m2_s": "12",
                "vela_electron_edge_source_per_s": "14",
                "vela_hole_edge_source_per_s": "16",
            }
        ]
        result = comparison.aggregate_edges(rows, set(range(1, 7)))
        self.assertEqual(result[1]["incident_edge_count"], 1)
        self.assertEqual(result[1]["vela_electron_alpha_mean_per_m"], 6)
        self.assertEqual(result[1]["vela_electron_source_node_per_s"], 7)
        self.assertEqual(result[2]["vela_total_source_node_per_s"], 15)
        self.assertEqual(result[6]["incident_edge_count"], 0)

    def test_join_requires_exact_six_nodes_per_topology_bias(self):
        rows = [
            {"topology_id": "sketch", "bias_V": "0", "node_id": str(node)}
            for node in range(1, 6)
        ]
        with self.assertRaisesRegex(ValueError, "exactly nodes 1..6"):
            comparison.validate_node_groups(rows)


if __name__ == "__main__":
    unittest.main()
