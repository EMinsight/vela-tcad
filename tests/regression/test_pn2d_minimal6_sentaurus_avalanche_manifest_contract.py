#!/usr/bin/env python3
"""Tests for cross-root Sentaurus avalanche evidence contracts."""

from __future__ import annotations

import copy
import unittest

from scripts.compare_pn2d_minimal6_sentaurus_avalanche_corrected_controls import (
    EXPECTED_SOURCE_VARIANTS,
    validate_source_manifests,
)


def source_manifest(label: str) -> dict:
    variants = EXPECTED_SOURCE_VARIANTS[label]
    topologies = {}
    for topology, tdr_hash in (
        ("mirror", "mirror_tdr"),
        ("sketch", "sketch_tdr"),
    ):
        topologies[topology] = {
            variant: {
                "status": "passed",
                "bundle_sha256": {
                    "pn2d_minimal6.tdr": tdr_hash,
                    "models.par": "models_hash",
                },
            }
            for variant in variants
        }
    return {
        "status": "passed",
        "sentaurus_release": "O-2018.06-SP2",
        "variants": list(variants),
        "topologies": topologies,
    }


class SentaurusAvalancheManifestContractTest(unittest.TestCase):
    def setUp(self) -> None:
        self.manifests = {
            label: source_manifest(label)
            for label in EXPECTED_SOURCE_VARIANTS
        }

    def test_matching_release_and_static_bundles_pass(self) -> None:
        release, hashes = validate_source_manifests(self.manifests)
        self.assertEqual(release, "O-2018.06-SP2")
        self.assertEqual(
            hashes["mirror"]["pn2d_minimal6.tdr"],
            "mirror_tdr",
        )

    def test_release_mismatch_is_rejected(self) -> None:
        changed = copy.deepcopy(self.manifests)
        changed["contact"]["sentaurus_release"] = "P-2019.03"
        with self.assertRaisesRegex(ValueError, "releases differ"):
            validate_source_manifests(changed)

    def test_static_bundle_mismatch_is_rejected(self) -> None:
        changed = copy.deepcopy(self.manifests)
        changed["contact"]["topologies"]["mirror"][
            "grad_qf_use_qf_contacts"
        ]["bundle_sha256"]["models.par"] = "different_models"
        with self.assertRaisesRegex(ValueError, "static bundle hash mismatch"):
            validate_source_manifests(changed)

    def test_variant_matrix_mismatch_is_rejected(self) -> None:
        changed = copy.deepcopy(self.manifests)
        changed["contact"]["variants"].reverse()
        with self.assertRaisesRegex(ValueError, "variant matrix mismatch"):
            validate_source_manifests(changed)


if __name__ == "__main__":
    unittest.main()
