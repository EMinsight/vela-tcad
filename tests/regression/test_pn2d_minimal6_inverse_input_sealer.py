import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]


class PN2DMinimal6InverseInputSealerTest(unittest.TestCase):
    def test_report_and_verifier_help_work_outside_repository(self) -> None:
        scripts = (
            REPO / "scripts" / "diagnose_pn2d_minimal6_physics_inverse_audit.py",
            REPO / "scripts" / "verify_pn2d_minimal6_physics_inverse_audit.py",
            REPO / "scripts" / "seal_pn2d_minimal6_inverse_inputs.py",
        )
        with tempfile.TemporaryDirectory() as tmp:
            for script in scripts:
                result = subprocess.run(
                    [sys.executable, str(script), "--help"], cwd=tmp,
                    text=True, capture_output=True, check=False,
                )
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn("usage:", result.stdout.lower())


if __name__ == "__main__":
    unittest.main()
