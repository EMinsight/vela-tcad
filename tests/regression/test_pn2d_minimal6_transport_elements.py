import csv
import tempfile
import unittest
from pathlib import Path

from scripts.export_pn2d_minimal6_transport_elements import (
    FIELD_SPECS,
    cell_field,
    read_cell_values,
)


def valid_field(name):
    components, unit = FIELD_SPECS[name]
    return {
        "index": 11,
        "name": name,
        "support_kind": "cell",
        "location_type": 3,
        "components": components,
        "values": 4,
        "raw_value_count": 4 * components,
        "global_node_mapping": "region_cell_order",
        "mapping_status": "complete",
        "unit": unit,
        "csv_file": f"{name}_region0_cells.csv",
    }


class Minimal6TransportElementExportTest(unittest.TestCase):
    def test_contract_accepts_all_seven_native_element_fields(self):
        for name in FIELD_SPECS:
            selected = cell_field({"fields": [valid_field(name)]}, name)
            self.assertEqual(selected["name"], name)

    def test_contract_rejects_node_or_wrong_shape(self):
        for key, value in (
            ("support_kind", "node"),
            ("location_type", 0),
            ("values", 6),
            ("mapping_status", "partial"),
        ):
            field = valid_field("ElectricField")
            field[key] = value
            with self.assertRaises(ValueError):
                cell_field({"fields": [field]}, "ElectricField")

    def test_reads_scalar_and_vector_cells(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for components in (1, 2):
                path = root / f"field{components}.csv"
                with path.open("w", newline="", encoding="utf-8") as handle:
                    writer = csv.writer(handle)
                    writer.writerow(
                        ["cell_id"]
                        + [
                            f"component{index}"
                            for index in range(components)
                        ]
                    )
                    for cell in range(4):
                        writer.writerow(
                            [cell]
                            + [
                                cell + index / 10.0
                                for index in range(components)
                            ]
                        )
                rows = read_cell_values(path, components)
                self.assertEqual([row[0] for row in rows], [0, 1, 2, 3])
                self.assertEqual(len(rows[0][1]), components)


if __name__ == "__main__":
    unittest.main()
