#!/usr/bin/env python3
"""Freeze and verify the TransportModels DD/DG parity-improvement baseline."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
BASELINE = REPO_ROOT / "build-release/reference_tcad/transportmodels_sentaurus2022/vela_baseline"
DEFAULT_JSON = REPO_ROOT / "docs/validation/transportmodels_dg_parity_baseline_2026-08-20.json"
DEFAULT_MD = REPO_ROOT / "docs/validation/transportmodels_dg_parity_baseline_2026-08-20.md"

CURVES = {
    "dd_idvg": (
        BASELINE / "generated/reference_curves/transportmodels_sentaurus2022_dd_idvg_reference.csv",
        BASELINE / "workflow_dd_vector_run01/dd_idvg_curve_comparison_candidate.csv",
    ),
    "dd_idvd": (
        BASELINE / "generated/reference_curves/transportmodels_sentaurus2022_dd_idvd_reference.csv",
        BASELINE / "workflow_dd_vector_run01/dd_idvd_curve_comparison_candidate.csv",
    ),
    "dg_idvg": (
        BASELINE / "generated/reference_curves/transportmodels_sentaurus2022_dg_idvg_reference.csv",
        BASELINE / "workflow_dg_outer80_resume_m036_run01/dg_idvg_curve_comparison_candidate.csv",
    ),
    "dg_idvd": (
        BASELINE / "generated/reference_curves/transportmodels_sentaurus2022_dg_idvd_reference.csv",
        BASELINE / "workflow_dg_idvd_resume_1p9_direct2p0_fixed1p5_outer200_run01/dg_idvd_curve_comparison_candidate.csv",
    ),
}

ARTIFACTS = {
    "dd_workflow_manifest": BASELINE / "workflow_dd_vector_run01/workflow_manifest.json",
    "dg_idvg_candidate": CURVES["dg_idvg"][1],
    "dg_idvd_candidate": CURVES["dg_idvd"][1],
    "dg_idvd_2V_state": BASELINE / "workflow_dg_idvd_resume_1p9_direct2p0_fixed1p5_outer200_run01/dg_idvd_curve_state_bias_2p000000.csv",
    "dg_final_workflow_manifest": BASELINE / "workflow_dg_idvd_resume_1p9_direct2p0_fixed1p5_outer200_run01/workflow_manifest.json",
    "dg_bias_regime_analysis": BASELINE / "workflow_dg_idvd_resume_1p9_direct2p0_fixed1p5_outer200_run01/reports/dg_bias_regime_analysis.json",
    "dg_spatial_summary": REPO_ROOT / "docs/progress_report_2026Q3/2026-08-19_transportmodels_dg_daily_report/figures/spatial_fields/transportmodels_dg_spatial_comparison_summary.json",
}

TARGETS = {
    "dd_control_idvd_max_relative_error": 0.02,
    "dd_control_idvg_on_max_relative_error": 0.20,
    "dg_idvd_max_relative_error": 0.05,
    "dg_idvd_endpoint_relative_error": 0.03,
    "dg_idvg_on_max_relative_error": 0.10,
    "dg_idvg_transition_max_absolute_log_error_dex": 0.15,
    "dg_quantum_potential_surface_p95_absolute_error_mV": 20.0,
    "dg_electron_density_surface_p95_absolute_log_error_dex": 0.20,
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest().upper()


def load_curve(path: Path, current_column: str) -> list[tuple[float, float]]:
    rows: list[tuple[float, float]] = []
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        for row in csv.DictReader(stream):
            rows.append((float(row["bias_V"]), abs(float(row[current_column]))))
    if len(rows) != 21 or len({bias for bias, _ in rows}) != 21:
        raise RuntimeError(f"Expected 21 unique biases in {path}, got {len(rows)}")
    return rows


def aligned_curve(name: str) -> list[tuple[float, float, float]]:
    reference_path, candidate_path = CURVES[name]
    reference = load_curve(reference_path, "current_total")
    candidate = load_curve(candidate_path, "current_total_A_per_um")
    result: list[tuple[float, float, float]] = []
    for (rbias, rvalue), (cbias, cvalue) in zip(reference, candidate):
        if not math.isclose(rbias, cbias, rel_tol=0.0, abs_tol=1.0e-12):
            raise RuntimeError(f"Bias mismatch for {name}: {rbias} vs {cbias}")
        if not all(math.isfinite(value) for value in (rvalue, cvalue)):
            raise RuntimeError(f"Non-finite current for {name} at {rbias} V")
        result.append((rbias, rvalue, cvalue))
    return result


def error_summary(rows: list[tuple[float, float, float]]) -> dict[str, float | int]:
    nonzero = [(bias, ref, vela) for bias, ref, vela in rows if ref > 0.0]
    rel = [abs(vela - ref) / ref for _, ref, vela in nonzero]
    dex = [abs(math.log10(vela / ref)) for _, ref, vela in nonzero if vela > 0.0]
    return {
        "points": len(rows),
        "median_relative_error": statistics.median(rel),
        "max_relative_error": max(rel),
        "median_absolute_log_error_dex": statistics.median(dex),
        "max_absolute_log_error_dex": max(dex),
        "endpoint_relative_error": abs(rows[-1][2] - rows[-1][1]) / rows[-1][1],
    }


def region_summary(
    rows: list[tuple[float, float, float]], predicate
) -> dict[str, float | int]:
    selected = [(bias, ref, vela) for bias, ref, vela in rows if predicate(bias)]
    rel = [abs(vela - ref) / ref for _, ref, vela in selected]
    dex = [abs(math.log10(vela / ref)) for _, ref, vela in selected]
    return {
        "points": len(selected),
        "median_absolute_log_error_dex": statistics.median(dex),
        "max_absolute_log_error_dex": max(dex),
        "max_relative_error": max(rel),
    }


def percentile_from_spatial_summary(
    summary_path: Path, field_key: str, source_field: str
) -> float:
    # The plotting summary freezes median/RMSE/max and the 99th-percentile color
    # limit. Phase 1 will add dedicated percentile metrics; use the current 99th
    # percentile as a conservative pre-calibration proxy in the frozen baseline.
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    return float(summary["fields"][field_key][source_field])


def build_snapshot() -> dict[str, object]:
    for path in [item for pair in CURVES.values() for item in pair] + list(ARTIFACTS.values()):
        if not path.is_file():
            raise FileNotFoundError(path)

    curves = {name: aligned_curve(name) for name in CURVES}
    metrics = {
        "dd_idvd": error_summary([row for row in curves["dd_idvd"] if row[0] > 0.0]),
        "dd_idvg_on": region_summary(curves["dd_idvg"], lambda bias: bias >= 0.28),
        "dg_idvd": error_summary([row for row in curves["dg_idvd"] if row[0] > 0.0]),
        "dg_idvg_off": region_summary(curves["dg_idvg"], lambda bias: bias <= -0.68),
        "dg_idvg_transition": region_summary(curves["dg_idvg"], lambda bias: -0.52 <= bias <= 0.12),
        "dg_idvg_on": region_summary(curves["dg_idvg"], lambda bias: bias >= 0.28),
        "dg_quantum_potential_surface_p99_absolute_error_mV": percentile_from_spatial_summary(
            ARTIFACTS["dg_spatial_summary"],
            "electron_quantum_potential",
            "difference_color_limit_99pct",
        ),
        "dg_electron_density_surface_p99_absolute_log_error_dex": percentile_from_spatial_summary(
            ARTIFACTS["dg_spatial_summary"],
            "electron_density",
            "difference_color_limit_99pct",
        ),
    }
    return {
        "schema": "vela.transportmodels.dg_parity_baseline.v1",
        "scope": {
            "sentaurus_version": "T-2022.03-SP2",
            "device": "TransportModels MOS",
            "curves": ["DD Id-Vg", "DD Id-Vd", "DG Id-Vg", "DG Id-Vd"],
            "curve_points_each": 21,
            "spatial_oracle": {"gate_bias_V": 1.0, "drain_bias_V": 2.0},
        },
        "status": "frozen",
        "artifact_hashes": {
            name: {"path": str(path.resolve()), "sha256": sha256(path)}
            for name, path in ARTIFACTS.items()
        },
        "curve_hashes": {
            name: {
                "reference_sha256": sha256(reference),
                "candidate_sha256": sha256(candidate),
            }
            for name, (reference, candidate) in CURVES.items()
        },
        "current_metrics": metrics,
        "target_acceptance": TARGETS,
        "qualification_policy": {
            "deep_off_idvg": "excluded from ordinary relative-error gates",
            "curve_bias_alignment_tolerance_V": 1.0e-12,
            "dd_is_non_regression_control": True,
            "targets_are_not_claimed_as_current_passes": True,
        },
    }


def render_markdown(snapshot: dict[str, object]) -> str:
    metrics = snapshot["current_metrics"]
    assert isinstance(metrics, dict)
    return f"""# TransportModels DG parity-improvement frozen baseline

Date: 2026-08-20

Status: **frozen and reproducible**. This baseline is the phase-0 control for
improving Vela/Sentaurus DG parity. It does not claim that the improvement
targets already pass.

## Frozen curve metrics

| Metric | Current value | Improvement target |
|---|---:|---:|
| DD Id-Vd maximum relative error | {metrics['dd_idvd']['max_relative_error']:.6%} | <= {TARGETS['dd_control_idvd_max_relative_error']:.2%} |
| DD Id-Vg on-state maximum relative error | {metrics['dd_idvg_on']['max_relative_error']:.6%} | <= {TARGETS['dd_control_idvg_on_max_relative_error']:.2%} |
| DG Id-Vd maximum relative error | {metrics['dg_idvd']['max_relative_error']:.6%} | <= {TARGETS['dg_idvd_max_relative_error']:.2%} |
| DG Id-Vd endpoint relative error | {metrics['dg_idvd']['endpoint_relative_error']:.6%} | <= {TARGETS['dg_idvd_endpoint_relative_error']:.2%} |
| DG Id-Vg on-state maximum relative error | {metrics['dg_idvg_on']['max_relative_error']:.6%} | <= {TARGETS['dg_idvg_on_max_relative_error']:.2%} |
| DG Id-Vg transition maximum log error | {metrics['dg_idvg_transition']['max_absolute_log_error_dex']:.6f} dex | <= {TARGETS['dg_idvg_transition_max_absolute_log_error_dex']:.3f} dex |
| DG Qn surface 99th-percentile absolute error | {metrics['dg_quantum_potential_surface_p99_absolute_error_mV']:.6f} mV | phase-1 target uses p95 <= {TARGETS['dg_quantum_potential_surface_p95_absolute_error_mV']:.1f} mV |
| DG electron-density surface 99th-percentile log error | {metrics['dg_electron_density_surface_p99_absolute_log_error_dex']:.6f} dex | phase-1 target uses p95 <= {TARGETS['dg_electron_density_surface_p95_absolute_log_error_dex']:.2f} dex |

## Reproduction

```powershell
D:\\msys64\\ucrt64\\bin\\python.exe scripts\\freeze_transportmodels_dg_parity_baseline.py --check
```

The JSON companion records absolute paths and SHA-256 values for every frozen
candidate, final state, workflow manifest, regional analysis, and spatial
summary. Curve files must contain 21 unique finite points aligned to the
Sentaurus bias lattice within 1e-12 V.

## Policy

- DD remains a non-regression control while DG is modified.
- Deep-off Id-Vg is excluded from ordinary relative-error gates.
- No later phase may relax nonlinear convergence tolerances to pass an
  amplitude target.
- A changed artifact hash requires an explicitly regenerated phase-0 baseline,
  not silent acceptance by `--check`.
"""


def compare_frozen(expected: dict[str, object], actual: dict[str, object]) -> list[str]:
    failures: list[str] = []
    for section in ("artifact_hashes", "curve_hashes", "current_metrics", "target_acceptance"):
        if expected.get(section) != actual.get(section):
            failures.append(section)
    return failures


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-json", type=Path, default=DEFAULT_JSON)
    parser.add_argument("--output-md", type=Path, default=DEFAULT_MD)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    output_json = args.output_json.resolve()
    output_md = args.output_md.resolve()
    snapshot = build_snapshot()
    if args.check:
        if not output_json.is_file():
            raise FileNotFoundError(output_json)
        expected = json.loads(output_json.read_text(encoding="utf-8"))
        failures = compare_frozen(expected, snapshot)
        if failures:
            print("Frozen baseline mismatch: " + ", ".join(failures))
            return 1
        print("TransportModels DG parity baseline check: PASS")
        return 0

    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(snapshot, indent=2) + "\n", encoding="utf-8")
    output_md.write_text(render_markdown(snapshot), encoding="utf-8")
    print(output_json)
    print(output_md)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
