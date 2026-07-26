from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

from scripts.diagnose_pn2d_imported_state_qfp_update import deck, impact_config, internal_nodes


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


if __name__ == "__main__":
    unittest.main()
