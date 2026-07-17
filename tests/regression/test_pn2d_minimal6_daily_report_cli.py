import subprocess
import sys
import unittest
from pathlib import Path


REPO = Path(__file__).resolve().parents[2]
SCRIPT = REPO / "scripts" / "render_pn2d_minimal6_daily_report_figures.py"


class DailyReportFigureCliTest(unittest.TestCase):
    def test_script_entrypoint_loads_repository_modules(self):
        result = subprocess.run(
            [sys.executable, str(SCRIPT), "--help"],
            cwd=REPO,
            text=True,
            capture_output=True,
            check=False,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("--comparison-report", result.stdout)
        self.assertIn("--out-dir", result.stdout)


if __name__ == "__main__":
    unittest.main()
