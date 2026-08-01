from __future__ import annotations

import unittest

from scripts.run_pn2d_forward_node_volume_policy_acceptance import anchor_metrics


class ForwardNodeVolumePolicyAcceptanceTests(unittest.TestCase):
    def test_anchor_metrics_report_improvement_and_median(self) -> None:
        sentaurus = {bias: float(bias) for bias in (1, 2, 5, 10, 15, 20)}
        baseline = {float(bias): 1.01 * bias for bias in sentaurus}
        candidate = {float(bias): 1.001 * bias for bias in sentaurus}
        metrics = anchor_metrics(candidate, baseline, sentaurus)
        self.assertAlmostEqual(metrics["candidate_median_relative_error"], 0.001)
        self.assertLess(metrics["maximum_error_degradation_over_barycentric"], 0.0)


if __name__ == "__main__":
    unittest.main()
