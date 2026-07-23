import csv
import tempfile
import unittest
from pathlib import Path

from scripts.export_pn2d_minimal6_element_currents import (
    cell_field,
    read_cell_vectors,
)


def valid_manifest():
    return {
        "fields": [
            {
                "index": 11,
                "name": "eCurrentDensity",
                "support_kind": "cell",
                "location_type": 3,
                "components": 2,
                "values": 4,
                "raw_value_count": 8,
                "global_node_mapping": "region_cell_order",
                "mapping_status": "complete",
                "unit": "A*cm^-2",
                "csv_file": "eCurrentDensity_region0_cells.csv",
            },
            {
                "index": 13,
                "name": "eCurrentDensity",
                "support_kind": "node",
                "location_type": 0,
                "components": 2,
                "values": 6,
                "raw_value_count": 12,
                "global_node_mapping": "global_vertex_order",
                "mapping_status": "complete",
                "unit": "A*cm^-2",
                "csv_file": "eCurrentDensity_region0.csv",
            },
        ]
    }


class Minimal6ElementCurrentExportTest(unittest.TestCase):
    def test_selects_only_complete_native_cell_vector(self):
        field = cell_field(valid_manifest(), "eCurrentDensity")
        self.assertEqual(field["index"], 11)
        self.assertEqual(field["csv_file"], "eCurrentDensity_region0_cells.csv")

    def test_rejects_node_only_current(self):
        manifest = valid_manifest()
        manifest["fields"] = manifest["fields"][1:]
        with self.assertRaisesRegex(ValueError, "expected one cell field"):
            cell_field(manifest, "eCurrentDensity")

    def test_rejects_wrong_location_or_shape(self):
        for key, value in (
            ("location_type", 0),
            ("components", 1),
            ("values", 6),
            ("raw_value_count", 12),
        ):
            manifest = valid_manifest()
            manifest["fields"][0][key] = value
            with self.assertRaisesRegex(ValueError, key):
                cell_field(manifest, "eCurrentDensity")

    def test_reads_canonical_four_cell_vectors(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "cells.csv"
            with path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle)
                writer.writerow(["cell_id", "component0", "component1"])
                for cell in range(4):
                    writer.writerow([cell, cell + 0.25, -cell - 0.5])
            vectors = read_cell_vectors(path)
        self.assertEqual(len(vectors), 4)
        self.assertEqual(vectors[0], (0, 0.25, -0.5))
        self.assertEqual(vectors[3], (3, 3.25, -3.5))


if __name__ == "__main__":
    unittest.main()
