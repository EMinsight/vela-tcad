from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from scripts.diagnose_pn2d_imported_state_qfp_update import deck, impact_config, internal_nodes
from scripts.verify_pn2d_imported_state_qfp_update import derive_causality


class ImportedStateQfpUpdateContractTest(unittest.TestCase):
    def setUp(self) -> None:
        self.base = {
            "contacts": [
                {"name": "Anode", "bias": 0.0},
                {"name": "Cathode", "bias": 0.0},
            ],
            "solver": {"mobility": "constant", "recombination": ["srh"]},
        }

    def test_production_triangle_is_frozen_canonical_gss_configuration(self) -> None:
        impact = impact_config("production_triangle")
        self.assertEqual(impact["model"], "van_overstraeten")
        self.assertEqual(impact["driving_force"], "quasi_fermi_gradient")
        self.assertEqual(impact["current_approximation"], "cell_reconstructed")
        self.assertEqual(impact["source_mapping_mode"], "triangle_gss_gradqf_truncated")
        self.assertEqual(impact["source_geometry_scale"], 1.0)

    def test_element_edge_candidate_stays_opt_in_and_unfitted(self) -> None:
        impact = impact_config("element_edge_opt_in")
        self.assertEqual(impact["model"], "van_overstraeten")
        self.assertEqual(impact["driving_force"], "quasi_fermi_gradient")
        self.assertEqual(impact["current_approximation"], "element_edge_sg_gss_laux")
        self.assertEqual(impact["source_mapping_mode"], "element_vertex_box_measure")
        self.assertNotIn("A_scale", impact)
        self.assertNotIn("B_scale", impact)
        self.assertNotIn("source_geometry_scale", impact)

    def test_controls_change_only_the_requested_source_family(self) -> None:
        full = deck(self.base, 20, Path("fields"), Path("full.csv"), "newton_carrier_term_probe", "element_edge_opt_in")
        avalanche_off = deck(self.base, 20, Path("fields"), Path("off.csv"), "newton_carrier_term_probe", "element_edge_opt_in", "avalanche_off")
        srh_off = deck(self.base, 20, Path("fields"), Path("srh.csv"), "newton_carrier_term_probe", "element_edge_opt_in", "srh_off")
        self.assertNotIn("impact_ionization", avalanche_off["solver"])
        self.assertEqual(avalanche_off["solver"]["recombination"], ["srh"])
        self.assertEqual(srh_off["solver"]["recombination"], [])
        self.assertEqual(srh_off["solver"]["impact_ionization"], full["solver"]["impact_ionization"])
        self.assertEqual(full["contacts"][0]["bias"], -20.0)
        self.assertEqual(full["contacts"][1]["bias"], 0.0)
        self.assertEqual(
            full["solver"]["mobility"],
            {"model": "masetti_field", "high_field_driving_force": "quasi_fermi_gradient"},
        )

    def test_internal_nodes_exclude_all_contact_nodes(self) -> None:
        mesh = {
            "nodes": [{"id": value} for value in range(5)],
            "contacts": [{"node_ids": [0, 1]}, {"node_ids": [4]}],
        }
        with tempfile.TemporaryDirectory() as temporary:
            path = Path(temporary) / "mesh.json"
            path.write_text(json.dumps(mesh), encoding="utf-8")
            self.assertEqual(internal_nodes({"mesh_file": str(path)}), {2, 3})


    def test_independent_causality_ignores_equal_low_signal_then_rejects_opposite_directions(self) -> None:
        residuals = []
        updates = []
        topologies = ("minimal6_mirror", "minimal6_sketch", "coarse7x3")
        for topology in topologies:
            for bias in (1, 10, 20):
                for carrier in ("electron", "hole"):
                    baseline = 1.0
                    if bias == 1:
                        candidate = baseline
                    elif topology == "coarse7x3":
                        candidate = 0.5
                    else:
                        candidate = 1.5
                    for variant, value in (
                        ("production_triangle", baseline),
                        ("element_edge_opt_in", candidate),
                    ):
                        residuals.append({
                            "topology": topology, "bias_V": str(-bias),
                            "variant": variant, "carrier": carrier,
                            "is_boundary": "0",
                            "final_residual_normalized": str(value),
                        })
                        for mode in ("carrier_only", "coupled"):
                            updates.append({
                                "topology": topology, "bias_V": str(-bias),
                                "variant": variant, "mode": mode,
                                "carrier": carrier, "delta_qfp_V": str(value),
                            })
        groups, authorized = derive_causality(residuals, updates)
        self.assertFalse(authorized)
        first_bias_directions = {
            (row["topology"], row["direction"])
            for row in groups if row["bias_V"] == -10
        }
        self.assertIn(("coarse7x3", "improved"), first_bias_directions)
        self.assertIn(("minimal6_mirror", "worsened"), first_bias_directions)
        self.assertTrue(all(row["direction"] == "equal" for row in groups if row["bias_V"] == -1))


if __name__ == "__main__":
    unittest.main()
