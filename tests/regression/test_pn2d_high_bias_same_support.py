from __future__ import annotations

import unittest

from scripts.analyze_pn2d_high_bias_same_support import (
    hotspot_evidence,
    support_class_summary,
)


class HighBiasSameSupportTests(unittest.TestCase):
    @staticmethod
    def row(bias: float, cell: int, generation: float, density: float) -> dict[str, object]:
        return {
            "bias_V": bias,
            "cell_id": cell,
            "node_ids": "0;1;2",
            "contact_class": "interior" if cell == 1 else "contact_adjacent",
            "active_region": 1 if cell == 1 else 0,
            "mean_n_cm3": density,
            "mean_p_cm3": 10.0,
            "efield_V_cm": 100.0,
            "grad_qf_n_V_cm": 100.0,
            "grad_qf_p_V_cm": 100.0,
            "mu_n_cm2_Vs": 1000.0,
            "mu_p_cm2_Vs": 500.0,
            "current_n_A_cm2": density,
            "current_p_A_cm2": 1.0,
            "max_alpha_n_cm_inv": 2.0,
            "max_alpha_p_cm_inv": 1.0,
            "max_generation_total_cm3_s": generation,
        }

    def test_fixed_minus20_hotspot_preserves_same_cell_process_order(self) -> None:
        rows = [
            self.row(-19.0, 1, 2.0, 2.0),
            self.row(-20.0, 1, 8.0, 4.0),
            self.row(-19.0, 2, 1.0, 1.0),
            self.row(-20.0, 2, 2.0, 1.0),
        ]
        chain, summary = hotspot_evidence(rows)
        self.assertEqual(summary["hotspot_cell_id"], 1)
        self.assertEqual(summary["first_material_stage"], "density")
        self.assertEqual([row["cell_id"] for row in chain], [1, 1])
        self.assertEqual(summary["metric_ratios_m20_over_m19"]["mean_n_cm3"], 2.0)

    def test_support_summary_keeps_contact_and_active_classes_separate(self) -> None:
        rows = [
            self.row(-20.0, 1, 8.0, 4.0),
            self.row(-20.0, 2, 2.0, 1.0),
        ]
        summary = support_class_summary(rows)
        keys = {
            (row["contact_class"], row["active_region"], row["cell_count"])
            for row in summary
        }
        self.assertEqual(
            keys,
            {("interior", 1, 1), ("contact_adjacent", 0, 1)},
        )


if __name__ == "__main__":
    unittest.main()
