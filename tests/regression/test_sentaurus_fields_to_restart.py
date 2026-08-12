import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]


class SentaurusFieldsToRestartTest(unittest.TestCase):
    def test_merges_regions_and_converts_density_to_m3(self) -> None:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fields = root / "fields"
            fields.mkdir()
            (root / "nodes.csv").write_text("node_id,x,y\n0,0,0\n1,1,0\n")
            entries = []
            for name in ("ElectrostaticPotential", "eQuasiFermiPotential",
                         "hQuasiFermiPotential", "eQuantumPotential"):
                for region, node in enumerate((0, 1)):
                    filename = f"{name}_region{region}.csv"
                    (fields / filename).write_text(f"node_id,component0\n{node},{node + 0.1}\n")
                    entries.append({"name": name, "csv_file": filename,
                                    "unit": "V", "mapping_status": "complete"})
            for name in ("eDensity", "hDensity"):
                filename = f"{name}_region0.csv"
                (fields / filename).write_text("node_id,component0\n0,2e10\n")
                entries.append({"name": name, "csv_file": filename,
                                "unit": "cm^-3", "mapping_status": "complete"})
            (root / "field_manifest.json").write_text(json.dumps({"fields": entries}))
            output = root / "state.csv"
            subprocess.run([sys.executable, str(REPO / "scripts" / "sentaurus_fields_to_restart.py"),
                            "--export-dir", str(root), "--output", str(output)], check=True)
            with output.open(newline="") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual(len(rows), 2)
            self.assertEqual(float(rows[0]["electrons_m3"]), 2e16)
            self.assertEqual(float(rows[1]["electrons_m3"]), 0.0)
            self.assertEqual(float(rows[1]["phin"]), 0.0)
            self.assertEqual(float(rows[1]["phip"]), 0.0)


if __name__ == "__main__":
    unittest.main()
