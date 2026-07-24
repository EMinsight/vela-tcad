import unittest

from scripts.pn2d_minimal6_diagnostics.highfield_box_replay import (
    BRANCH_ID,
    CELL_MAPPING,
    RECONSTRUCTION_LABEL,
    REQUIRED_SAMPLE_FIELDS,
    coefficient_weighted_mobility,
    validate_sample_record,
)


class Minimal6HighfieldBoxCurrentContractTest(unittest.TestCase):
    def test_typed_branch_and_reconstruction_label_are_frozen(self) -> None:
        self.assertEqual(
            BRANCH_ID,
            "sentaurus_lowfield_element_electric_field",
        )
        self.assertEqual(
            RECONSTRUCTION_LABEL,
            "box_operator_reconstruction",
        )

    def test_region_cell_mapping_is_frozen(self) -> None:
        self.assertEqual(CELL_MAPPING["mirror"], (0, 1, 2, 3))
        self.assertEqual(CELL_MAPPING["sketch"], (0, 3, 2, 1))

    def test_coefficient_weighted_mobility(self) -> None:
        result = coefficient_weighted_mobility(
            [(0.25, 0.10), (0.75, 0.20)]
        )
        self.assertEqual(result["status"], "valid")
        self.assertAlmostEqual(result["mobility_m2_per_Vs"], 0.175)
        self.assertAlmostEqual(result["kappa_sum"], 1.0)

    def test_geometric_zero_is_typed(self) -> None:
        result = coefficient_weighted_mobility(
            [(0.0, 0.10), (0.0, 0.20)]
        )
        self.assertEqual(result["status"], "geometric_zero")
        self.assertIsNone(result["mobility_m2_per_Vs"])
        self.assertEqual(result["kappa_sum"], 0.0)

    def test_sample_schema_rejects_missing_support_provenance(self) -> None:
        record = {field: "sealed" for field in REQUIRED_SAMPLE_FIELDS}
        record["topology"] = "mirror"
        record["status"] = "valid"
        record["reconstruction_label"] = RECONSTRUCTION_LABEL
        validate_sample_record(record)
        del record["electric_field_source_sha256"]
        with self.assertRaisesRegex(
            ValueError,
            "electric_field_source_sha256",
        ):
            validate_sample_record(record)


if __name__ == "__main__":
    unittest.main()
