#!/usr/bin/env python3
"""Build the portable technical-report artifact for the mobility unit audit."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


TITLE = "PN2D Minimal6 mobility unit root-cause audit"


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _number(value: str) -> int | float:
    parsed = float(value)
    return int(parsed) if parsed.is_integer() else parsed


def build_artifact(root: Path) -> dict[str, object]:
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    summary = [
        {
            "support": row["support"],
            "carrier": row["carrier"],
            "branch": row["branch"],
            "n": int(row["sample_count"]),
            "median_abs_dex": _number(row["median_abs_log10_error_dex"]),
            "p95_abs_dex": _number(row["p95_abs_log10_error_dex"]),
            "maximum_abs_dex": _number(row["maximum_abs_log10_error_dex"]),
        }
        for row in _rows(root / "summary.csv")
    ]
    code_audit = [
        {
            "component": "Mobility default",
            "declared_physical_unit": "m2/(V s)",
            "unit_scaled_internal_unit": "cm2/(V s)",
            "conversion_status": "converted",
        },
        {
            "component": "QFP-gradient field",
            "declared_physical_unit": "V/m",
            "unit_scaled_internal_unit": "V/cm",
            "conversion_status": "converted",
        },
        {
            "component": "Saturation velocity",
            "declared_physical_unit": "m/s",
            "unit_scaled_internal_unit": "cm/s",
            "conversion_status": "not converted",
        },
        {
            "component": "High-field ratio mu E / vsat",
            "declared_physical_unit": "dimensionless",
            "unit_scaled_internal_unit": "dimensionless",
            "conversion_status": "100x too large",
        },
    ]
    sources = [
        {
            "id": "root_cause_output",
            "label": "Deterministic 40-state mobility unit audit",
            "path": (
                "build-release/pn2d-minimal6-mobility-unit-root-cause-"
                "20260723-a/manifest.json"
            ),
        },
        {
            "id": "summary_query",
            "label": "Mobility unit root-cause summary",
            "path": (
                "build-release/pn2d-minimal6-mobility-unit-root-cause-"
                "20260723-a/summary.csv"
            ),
            "query": {
                "engine": "duckdb",
                "language": "sql",
                "description": (
                    "Read all support-aligned mobility error branches in "
                    "deterministic output order."
                ),
                "sql": (
                    "SELECT support, carrier, branch, sample_count AS n, "
                    "median_abs_log10_error_dex AS median_abs_dex, "
                    "p95_abs_log10_error_dex AS p95_abs_dex, "
                    "maximum_abs_log10_error_dex AS maximum_abs_dex "
                    "FROM read_csv_auto('summary.csv', header = true)"
                ),
            },
        },
        {
            "id": "production_units",
            "label": "Vela unit and mobility production implementation",
            "path": "src/physics/MobilityModel.cpp",
        },
        {
            "id": "native_elements",
            "label": "Native Sentaurus element transport export",
            "path": (
                "build-release/pn2d-minimal6-transport-elements-20260723-b/"
                "export/transport_element_values.csv"
            ),
        },
    ]
    return {
        "surface": "report",
        "manifest": {
            "version": 1,
            "surface": "report",
            "title": TITLE,
            "generatedAt": "2026-07-23T00:00:00+08:00",
            "cards": [],
            "charts": [],
            "tables": [
                {
                    "id": "summary_table",
                    "title": "Support-aligned mobility error decomposition",
                    "description": (
                        "40 states; direct local-edge and native element "
                        "comparisons; absolute log10 mobility-ratio error"
                    ),
                    "dataset": "summary",
                    "sourceId": "summary_query",
                    "columns": [
                        {"field": "support", "label": "Support"},
                        {"field": "carrier", "label": "Carrier"},
                        {"field": "branch", "label": "Branch"},
                        {"field": "n", "label": "N", "format": "number"},
                        {
                            "field": "median_abs_dex",
                            "label": "Median abs error (dex)",
                            "format": "number",
                        },
                        {
                            "field": "p95_abs_dex",
                            "label": "P95 abs error (dex)",
                            "format": "number",
                        },
                        {
                            "field": "maximum_abs_dex",
                            "label": "Maximum abs error (dex)",
                            "format": "number",
                        },
                    ],
                    "defaultSort": {"field": "support", "direction": "asc"},
                },
                {
                    "id": "code_audit_table",
                    "title": "Unit-scaling path",
                    "description": (
                        "Physical declarations versus internal arithmetic "
                        "used by the high-field limiter"
                    ),
                    "dataset": "code_audit",
                    "sourceId": "production_units",
                    "columns": [
                        {"field": "component", "label": "Component"},
                        {
                            "field": "declared_physical_unit",
                            "label": "Declared physical unit",
                        },
                        {
                            "field": "unit_scaled_internal_unit",
                            "label": "Internal unit",
                        },
                        {
                            "field": "conversion_status",
                            "label": "Conversion status",
                        },
                    ],
                    "defaultSort": {"field": "component", "direction": "asc"},
                },
            ],
            "blocks": [
                {"id": "title", "type": "markdown", "body": f"# {TITLE}"},
                {
                    "id": "technical_summary",
                    "type": "markdown",
                    "sourceId": "root_cause_output",
                    "body": (
                        "## Technical summary\n\n"
                        "The dominant Vela mobility discrepancy is a "
                        "saturation-velocity conversion defect in "
                        "`unit_scaling`, not a failure of the Masetti "
                        "low-field formula. Mobility and QFP-gradient field "
                        "are converted to cm2/(V s) and V/cm, but a value "
                        "declared as m/s remains numerically unchanged and is "
                        "therefore consumed as cm/s. The effective physical "
                        "saturation velocity is 100 times too small."
                    ),
                },
                {
                    "id": "direct_finding",
                    "type": "markdown",
                    "sourceId": "root_cause_output",
                    "body": (
                        "## Direct C++ values identify the exact legacy "
                        "interpretation\n\n"
                        "Across 480 local-edge samples per carrier, direct "
                        "C++ mobility matches the unconverted-velocity branch "
                        "at a median error of 9.64e-17 dex for both carriers. "
                        "The correctly converted branch differs from the "
                        "current C++ output by 1.875 dex for electrons and "
                        "1.790 dex for holes. This machine-precision closure "
                        "turns the unit explanation from a scale hypothesis "
                        "into a verified arithmetic reproduction."
                    ),
                },
                {
                    "id": "summary_block",
                    "type": "table",
                    "tableId": "summary_table",
                },
                {
                    "id": "element_finding",
                    "type": "markdown",
                    "sourceId": "root_cause_output",
                    "body": (
                        "## Native element support confirms the physical "
                        "correction\n\n"
                        "Using the native Sentaurus element QFP-gradient "
                        "support and cell-average doping, the legacy branch "
                        "has median errors of 1.877 dex (electron) and 1.839 "
                        "dex (hole). Correct velocity conversion reduces "
                        "those medians to 0.0527 and 0.0478 dex. Support "
                        "alignment is essential: averaging Vela local-edge "
                        "mobilities is a different operator because a "
                        "zero-QFP-difference edge retains low-field mobility."
                    ),
                },
                {
                    "id": "code_audit_finding",
                    "type": "markdown",
                    "body": (
                        "## The conversion gap enters only through velocity\n\n"
                        "The high-field ratio is `mu * E / vsat`. In the "
                        "unit-scaled path, the numerator is expressed in "
                        "cm/s. Both default values and explicit JSON values "
                        "named `*_saturation_velocity_m_s` bypass a velocity "
                        "conversion, so the denominator remains numerically "
                        "an m/s value while being interpreted as cm/s."
                    ),
                },
                {
                    "id": "code_audit_block",
                    "type": "table",
                    "tableId": "code_audit_table",
                },
                {
                    "id": "scope",
                    "type": "markdown",
                    "sourceId": "root_cause_output",
                    "body": (
                        "## Scope, data, and definitions\n\n"
                        "The audit covers 40 exact states: mirror and sketch "
                        "topologies at -1 through -20 V. Direct evidence has "
                        "960 carrier-local-edge rows; native element evidence "
                        "has 320 carrier-element rows. The metric is the "
                        "absolute base-10 logarithm of the positive mobility "
                        "ratio. The element reconstruction uses affine cell "
                        "QFP-gradient magnitude and arithmetic cell-average "
                        "net doping."
                    ),
                },
                {
                    "id": "method",
                    "type": "markdown",
                    "body": (
                        "## Methodology and robustness\n\n"
                        "The replay keeps Masetti parameters, field exponent, "
                        "doping, QFP state, and support fixed. The two branches "
                        "differ only by saturation-velocity scale: 1.0 for "
                        "physical m/s and 0.01 for the current unit-scaled "
                        "interpretation. No production formula or solver state "
                        "was modified. Exact row counts and input/output hashes "
                        "are sealed in the diagnostic manifest."
                    ),
                },
                {
                    "id": "limitations",
                    "type": "markdown",
                    "body": (
                        "## Limitations and uncertainty\n\n"
                        "The fixed-state audit already exposes direct triangle "
                        "local-edge mobility but not the separate global SG "
                        "edge mobility. Sentaurus native element mobility is "
                        "not a directed-edge flux. Residual element errors up "
                        "to 0.313 dex for electrons and 0.184 dex for holes "
                        "can include doping interpolation, temperature "
                        "dependence, and solver-specific element evaluation."
                    ),
                },
                {
                    "id": "next_steps",
                    "type": "markdown",
                    "body": (
                        "## Recommended next steps\n\n"
                        "1. Add an explicit velocity conversion and apply it "
                        "to default and JSON-provided saturation velocities.\n"
                        "2. Add legacy-SI versus unit-scaled mobility parity "
                        "tests at fixed physical doping and QFP field.\n"
                        "3. Re-run the 40-state fixed-state and "
                        "self-consistent replacement audits.\n"
                        "4. Resume directed-edge SG current inversion only "
                        "after mobility parity is restored."
                    ),
                },
                {
                    "id": "further_questions",
                    "type": "markdown",
                    "body": (
                        "## Further question\n\n"
                        "After unit parity is restored, how much of the "
                        "remaining approximately 0.05 dex median element gap "
                        "comes from doping interpolation versus "
                        "Sentaurus-specific temperature dependence?"
                    ),
                },
            ],
            "sources": sources,
        },
        "snapshot": {
            "version": 1,
            "generatedAt": "2026-07-23T00:00:00+08:00",
            "status": "ready",
            "datasets": {
                "summary": summary,
                "code_audit": code_audit,
            },
        },
        "sources": sources,
        "package_info": {
            "report_shape": "technical",
            "chart_omission_reason": (
                "Exact tables preserve mixed-support branch semantics; a "
                "single chart would visually combine incomparable supports."
            ),
            "source_manifest": manifest["experiment"],
        },
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    artifact = build_artifact(args.root.resolve())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(artifact, indent=2, ensure_ascii=True) + "\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
