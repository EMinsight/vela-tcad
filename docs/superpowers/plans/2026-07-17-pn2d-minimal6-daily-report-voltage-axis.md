# PN2D Minimal6 Daily-Report Voltage Axis Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Generate three reproducible daily-report PNGs whose applied-voltage axis runs from `-1 V` on the left to `-20 V` on the right without changing any sealed Task 8 artifact.

**Architecture:** Add one standalone presentation-layer renderer that first verifies the existing Task 8 comparison package, selects only exact common checkpoints, and renders the three approved figures into a new directory. A separate manifest binds every PNG to the source report hash and records descending numeric axis limits; its verifier rerenders into a temporary directory and compares pixels.

**Tech Stack:** Python 3, Matplotlib, Pillow, `unittest`, existing `scripts.compare_pn2d_minimal6_diagnostic_sweeps` verification and plotting helpers.

## Global Constraints

- Generate exactly `terminal_current.png`, `maximum_field.png`, and `source_integrals.png`.
- Display signed applied voltage in the order `-1, -2, ..., -20 V` from left to right.
- Do not modify or overwrite the Task 8 comparison directory, report, figures, hashes, data, units, series identities, disclaimer, or scientific conclusions.
- Reject invalid source packages, missing exact common checkpoints, non-finite values, or non-descending display limits.
- Preserve the established 900 x 504 pixel figure size.

---

### Task 1: Daily-report renderer, contract, and final figures

**Files:**
- Create: `scripts/render_pn2d_minimal6_daily_report_figures.py`
- Create: `tests/regression/test_pn2d_minimal6_daily_report_figures.py`
- Generate: `build-release/pn2d-minimal6-daily-report-20260717/terminal_current.png`
- Generate: `build-release/pn2d-minimal6-daily-report-20260717/maximum_field.png`
- Generate: `build-release/pn2d-minimal6-daily-report-20260717/source_integrals.png`
- Generate: `build-release/pn2d-minimal6-daily-report-20260717/daily_report_figure_manifest.json`

**Interfaces:**
- Consumes: `Path` to a verified `sweep_comparison.json` with exact common `checkpoints` and the original Task 8 figure contract.
- Produces: `render_daily_report_figures(report_path: Path, out_dir: Path) -> dict[str, Any]`.
- Produces: `verify_daily_report_figures(manifest_path: Path) -> bool`.
- CLI: `python scripts/render_pn2d_minimal6_daily_report_figures.py --comparison-report <sweep_comparison.json> --out-dir <directory>`.

- [ ] **Step 1: Write the failing contract tests**

Create `tests/regression/test_pn2d_minimal6_daily_report_figures.py`. Reuse the existing sweep-comparison test helpers to create a fully verified synthetic comparison package, then require the new API and contract:

```python
import hashlib
import json
import tempfile
import unittest
from pathlib import Path

from PIL import Image
from tests.regression.test_pn2d_minimal6_sweep_comparison import (
    checkpoint, manifest, write_comparison_package,
)
from scripts.render_pn2d_minimal6_daily_report_figures import (
    DAILY_FIGURE_NAMES,
    render_daily_report_figures,
    verify_daily_report_figures,
)


class DailyReportFigureTest(unittest.TestCase):
    def comparison_package(self, root: Path) -> Path:
        biases = (-1.0, -2.0, -20.0)
        rows = {
            solver: [
                checkpoint(solver, topology, bias, bias, 1.0, 2.0, 3.0, 4.0)
                for topology in ("sketch", "mirror") for bias in biases
            ]
            for solver in ("vela", "sentaurus")
        }
        write_comparison_package(
            root,
            manifest("vela", rows["vela"]),
            manifest("sentaurus", rows["sentaurus"]),
            fixed_state_report={},
        )
        return root / "sweep_comparison.json"

    def test_contract_has_only_daily_figures_and_descending_voltage_axis(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            source = self.comparison_package(root / "source")
            output = root / "daily"
            contract = render_daily_report_figures(source, output)
            self.assertEqual(set(contract["figures"]), set(DAILY_FIGURE_NAMES))
            self.assertEqual(set(DAILY_FIGURE_NAMES), {
                "terminal_current.png", "maximum_field.png", "source_integrals.png"})
            for name, entry in contract["figures"].items():
                self.assertEqual(entry["x_quantity"], "applied_bias_V")
                self.assertEqual(entry["x_axis_order"], "decreasing_left_to_right")
                self.assertGreater(entry["x_limits_V"][0], entry["x_limits_V"][1])
                self.assertEqual(entry["x_limits_V"], [-1.0, -20.0])
                with Image.open(output / name) as image:
                    self.assertEqual(image.size, (900, 504))
                    self.assertEqual(image.info["VoltageAxisOrder"],
                                     "decreasing_left_to_right")
            self.assertTrue(verify_daily_report_figures(
                output / "daily_report_figure_manifest.json"))

    def test_rerender_is_byte_deterministic_and_source_is_unchanged(self):
        with tempfile.TemporaryDirectory() as temp:
            root = Path(temp)
            source = self.comparison_package(root / "source")
            before = hashlib.sha256(source.read_bytes()).hexdigest()
            render_daily_report_figures(source, root / "first")
            render_daily_report_figures(source, root / "second")
            self.assertEqual(before, hashlib.sha256(source.read_bytes()).hexdigest())
            for name in (*DAILY_FIGURE_NAMES, "daily_report_figure_manifest.json"):
                self.assertEqual((root / "first" / name).read_bytes(),
                                 (root / "second" / name).read_bytes())
```

- [ ] **Step 2: Run the tests and verify RED**

Run:

```powershell
D:\msys64\ucrt64\bin\python.exe -m unittest tests.regression.test_pn2d_minimal6_daily_report_figures -v
```

Expected: FAIL with `ModuleNotFoundError: No module named 'scripts.render_pn2d_minimal6_daily_report_figures'`.

- [ ] **Step 3: Implement the minimal standalone renderer**

Create `scripts/render_pn2d_minimal6_daily_report_figures.py` with these exact public constants and functions:

```python
DAILY_FIGURE_NAMES = (
    "terminal_current.png",
    "maximum_field.png",
    "source_integrals.png",
)
SCHEMA = "vela.pn2d_minimal6_daily_report_figures.v1"
AXIS_ORDER = "decreasing_left_to_right"

def render_daily_report_figures(report_path: Path, out_dir: Path) -> dict[str, Any]:
    comparison.verify_comparison_artifacts(report_path)
    report = json.loads(report_path.read_text(encoding="utf-8"),
                        parse_constant=_reject_nonfinite)
    common_biases = sorted({float(row["bias_V"]) for row in report["checkpoints"]},
                           reverse=True)
    if not common_biases:
        raise ValueError("daily-report figures require exact common checkpoints")
    x_limits = [common_biases[0], common_biases[-1]]
    if not x_limits[0] > x_limits[1]:
        raise ValueError("daily-report voltage axis must decrease left to right")
    if x_limits != [-1.0, -20.0]:
        raise ValueError("daily-report voltage axis must span -1 V to -20 V")
    daily_report = exact_common_checkpoint_view(report, common_biases)
    figures = render_three_figures(daily_report, out_dir, x_limits)
    contract = build_manifest(
        report_path=report_path.resolve(),
        source_sha256=sha256(report_path.read_bytes()),
        figures=figures,
    )
    write_manifest(out_dir / "daily_report_figure_manifest.json", contract)
    return contract

def verify_daily_report_figures(manifest_path: Path) -> bool:
    contract = load_and_validate_manifest(manifest_path)
    report_path = Path(contract["source_comparison"]["path"])
    comparison.verify_comparison_artifacts(report_path)
    verify_bound_source_hash(report_path, contract["source_comparison"]["sha256"])
    verify_png_contracts(manifest_path.parent, contract)
    rerender_and_compare_pixels(report_path, manifest_path.parent, contract)
    return True
```

The written manifest must have this exact shape:

```json
{
  "schema": "vela.pn2d_minimal6_daily_report_figures.v1",
  "source_comparison": {
    "path": "absolute path to sweep_comparison.json",
    "sha256": "64 lowercase hex characters"
  },
  "figures": {
    "terminal_current.png": {
      "source_comparison_path": "same absolute path",
      "source_comparison_sha256": "same SHA-256",
      "sha256": "PNG SHA-256",
      "width_px": 900,
      "height_px": 504,
      "x_quantity": "applied_bias_V",
      "x_axis_order": "decreasing_left_to_right",
      "x_limits_V": [-1.0, -20.0],
      "series_identities": []
    }
  }
}
```

Populate `series_identities` from the actual plotted Vela/Sentaurus and
sketch/mirror series. Use the established comparison helpers for titles,
labels, grid, disclaimer, failure markers, dimensions, and palette. Reject any
selected observable for which `math.isfinite(float(value))` is false.

The CLI must resolve both paths, call `render_daily_report_figures`, immediately
call `verify_daily_report_figures`, print the manifest path, and return zero.

- [ ] **Step 4: Run the focused tests and verify GREEN**

Run:

```powershell
D:\msys64\ucrt64\bin\python.exe -m unittest tests.regression.test_pn2d_minimal6_daily_report_figures -v
```

Expected: both tests PASS.

- [ ] **Step 5: Run the existing sweep-comparison regression suite**

Run:

```powershell
D:\msys64\ucrt64\bin\python.exe -m unittest tests.regression.test_pn2d_minimal6_sweep_comparison -v
```

Expected: all existing tests PASS, demonstrating that the main Task 8 renderer
and verification contract remain unchanged.

- [ ] **Step 6: Generate and verify the final daily-report figures**

Run:

```powershell
D:\msys64\ucrt64\bin\python.exe scripts\render_pn2d_minimal6_daily_report_figures.py `
  --comparison-report build-release\pn2d-minimal6-comparison-task8-fresh-20260717-a\sweep_comparison.json `
  --out-dir build-release\pn2d-minimal6-daily-report-20260717
```

Expected: exit code 0 and output ending in
`daily_report_figure_manifest.json`.

Run the CLI a second time with `--out-dir
build-release\pn2d-minimal6-daily-report-20260717-repeat`; compare all four
files byte-for-byte. Expected: identical files.

- [ ] **Step 7: Inspect the three final images and verify source integrity**

Open all three PNGs at their final size. Confirm `-1 V` is the leftmost labeled
voltage, `-20 V` is the rightmost, legends and disclaimer are readable, and no
text is mirrored or clipped. Re-run:

```powershell
D:\msys64\ucrt64\bin\python.exe -c "from pathlib import Path; from scripts.render_pn2d_minimal6_daily_report_figures import verify_daily_report_figures; assert verify_daily_report_figures(Path(r'build-release\pn2d-minimal6-daily-report-20260717\daily_report_figure_manifest.json'))"
```

Expected: exit code 0. Then run `git diff --check`; expected: no output.

- [ ] **Step 8: Commit the implementation**

```powershell
D:\msys64\usr\bin\git.exe add scripts/render_pn2d_minimal6_daily_report_figures.py tests/regression/test_pn2d_minimal6_daily_report_figures.py
D:\msys64\usr\bin\git.exe commit -m "Add decreasing-axis daily-report figures"
```

Do not add `build-release` outputs or the pre-existing untracked
`docs/validation/figures/` directory.
