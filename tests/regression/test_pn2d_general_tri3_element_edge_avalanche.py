#!/usr/bin/env python3
"""RED contracts for the general-Tri3 avalanche follow-up."""

from __future__ import annotations

import copy
import unittest
from pathlib import Path

from scripts.pn2d_general_tri3_contract import (
    EXACT_BIASES_V,
    SCHEMA_ID,
    SENTAURUS_RELEASE,
    validate_cell_mapping,
    validate_driver_contract,
    validate_record,
    validate_source_manifests,
)


REPO = Path(__file__).resolve().parents[2]
ANALYZER = REPO / "scripts" / "diagnose_pn2d_general_tri3_element_edge_avalanche.py"


def valid_record() -> dict:
    return {
        "case": "coarse7x3",
        "topology": "native",
        "bias_V": -10.0,
        "carrier": "electron",
        "vela_cell_id": 4,
        "sentaurus_region_cell_id": 7,
        "node_ids": [1, 2, 8],
        "local_edge": 0,
        "edge_start": 1,
        "edge_end": 2,
        "edge_orientation": 1,
        "contact_adjacent": False,
        "interior_element": True,
        "triangle_angles_deg": [50.0, 60.0, 70.0],
        "signed_area_um2": 0.1,
        "cell_orientation": "ccw",
        "read_coefficient": 0.4,
        "read_measure_um2": 0.03,
        "electric_field_source_sha256": "electric",
        "qfp_gradient_source_sha256": "qfp",
        "low_field_mobility_source_sha256": "mu0",
        "final_mobility_source_sha256": "mu",
        "coefficient_driving_force": "quasi_fermi_gradient",
        "current_density_approximation": "element_edge_sg",
        "observation_label": "box_operator_reconstruction",
        "support_status": "valid",
        "unit": "A/cm2",
        "absolute_error_dex": 0.01,
        "global_default_driver": "quasi_fermi_gradient",
        "effective_driver": "quasi_fermi_gradient",
        "use_quasi_fermi_at_contacts": False,
    }


def source_manifest() -> dict:
    return {
        "schema": SCHEMA_ID,
        "status": "passed",
        "sentaurus_release": SENTAURUS_RELEASE,
        "exact_biases_V": list(EXACT_BIASES_V),
        "case_hashes": {
            "coarse7x3": {
                "tdr": "mesh",
                "models.par": "models",
            }
        },
    }


class GeneralTri3ContractTest(unittest.TestCase):
    def test_valid_record_and_manifests_pass(self) -> None:
        record = valid_record()
        validate_record(record)
        validate_driver_contract(record)
        hashes = validate_source_manifests(
            {"implicit": source_manifest(), "explicit": source_manifest()}
        )
        self.assertEqual(hashes["coarse7x3"]["tdr"], "mesh")

    def test_wrong_carrier_or_edge_orientation_is_rejected(self) -> None:
        record = valid_record()
        record["edge_orientation"] = 0
        with self.assertRaisesRegex(ValueError, "carrier or edge orientation"):
            validate_record(record)

    def test_wrong_cell_permutation_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "wrong cell permutation"):
            validate_cell_mapping(
                [valid_record()],
                {("native", 4): 6},
            )

    def test_contact_fallback_is_not_global_default(self) -> None:
        record = valid_record()
        record.update(
            {
                "contact_adjacent": True,
                "interior_element": False,
                "effective_driver": "quasi_fermi_gradient",
            }
        )
        with self.assertRaisesRegex(ValueError, "contact fallback"):
            validate_driver_contract(record)

    def test_geometric_zero_cannot_have_finite_dex(self) -> None:
        record = valid_record()
        record.update(
            {
                "support_status": "geometric_zero",
                "read_coefficient": 0.0,
                "absolute_error_dex": 0.0,
            }
        )
        with self.assertRaisesRegex(ValueError, "finite dex"):
            validate_record(record)

    def test_reconstructed_edge_cannot_be_labeled_native(self) -> None:
        record = valid_record()
        record["observation_label"] = "native_element"
        with self.assertRaisesRegex(ValueError, "mislabeled as native"):
            validate_record(record)

    def test_static_bundle_hash_mismatch_is_rejected(self) -> None:
        changed = copy.deepcopy(source_manifest())
        changed["case_hashes"]["coarse7x3"]["models.par"] = "different"
        with self.assertRaisesRegex(ValueError, "static bundle hash mismatch"):
            validate_source_manifests(
                {"implicit": source_manifest(), "explicit": changed}
            )

    def test_sentaurus_release_mismatch_is_rejected(self) -> None:
        changed = copy.deepcopy(source_manifest())
        changed["sentaurus_release"] = "P-2019.03"
        with self.assertRaisesRegex(ValueError, "releases differ"):
            validate_source_manifests(
                {"implicit": source_manifest(), "explicit": changed}
            )

    def test_analyzer_entrypoint_is_present(self) -> None:
        self.assertTrue(
            ANALYZER.is_file(),
            "general-Tri3 avalanche analyzer is not implemented",
        )


if __name__ == "__main__":
    unittest.main()
