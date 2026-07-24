#!/usr/bin/env python3
"""Build the validator-complete portable mobility root-cause artifact."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from build_pn2d_minimal6_mobility_unit_root_cause_artifact import (
    build_artifact,
)


def build_portable_artifact(root: Path) -> dict[str, object]:
    artifact = build_artifact(root)
    manifest = artifact["manifest"]
    snapshot = artifact["snapshot"]
    summary = snapshot["datasets"]["summary"]
    chart_rows = [
        {
            "carrier": row["carrier"],
            "velocity_interpretation": (
                "Legacy, no velocity conversion"
                if row["branch"] == "legacy_cell_average_doping"
                else "Correct m/s conversion"
            ),
            "median_abs_dex": row["median_abs_dex"],
            "p95_abs_dex": row["p95_abs_dex"],
            "n": row["n"],
        }
        for row in summary
        if row["support"] == "sentaurus_native_element"
        and row["branch"]
        in ("legacy_cell_average_doping", "correct_cell_average_doping")
    ]
    snapshot["datasets"]["native_element_chart"] = chart_rows
    manifest["charts"] = [
        {
            "id": "native_element_error_chart",
            "title": "Native element mobility error",
            "description": (
                "Median absolute log10 error versus Sentaurus; "
                "160 elements per carrier and branch"
            ),
            "type": "bar",
            "dataset": "native_element_chart",
            "sourceId": "summary_query",
            "encodings": {
                "x": {"field": "carrier", "type": "nominal"},
                "y": {"field": "median_abs_dex", "type": "quantitative"},
                "color": {
                    "field": "velocity_interpretation",
                    "type": "nominal",
                },
            },
        }
    ]
    blocks = manifest["blocks"]
    finding_index = next(
        index
        for index, block in enumerate(blocks)
        if block["id"] == "element_finding"
    )
    blocks.insert(
        finding_index + 1,
        {
            "id": "native_element_chart_block",
            "type": "chart",
            "chartId": "native_element_error_chart",
        },
    )
    for source in manifest["sources"]:
        if source["id"] != "production_units":
            continue
        source["query"] = {
            "engine": "duckdb",
            "language": "sql",
            "description": (
                "Materialize the four reviewed unit-path audit rows used by "
                "the report table."
            ),
            "sql": (
                "SELECT * FROM (VALUES "
                "('Mobility default','m2/(V s)','cm2/(V s)','converted'),"
                "('QFP-gradient field','V/m','V/cm','converted'),"
                "('Saturation velocity','m/s','cm/s','not converted'),"
                "('High-field ratio mu E / vsat','dimensionless',"
                "'dimensionless','100x too large')) AS t("
                "component,declared_physical_unit,"
                "unit_scaled_internal_unit,conversion_status)"
            ),
        }
    artifact["sources"] = manifest["sources"]
    artifact["package_info"] = {
        "report_shape": "technical",
        "chart_map": (
            "Native-element finding | compare legacy versus correct velocity "
            "interpretation | grouped bar | x=carrier, y=median_abs_dex, "
            "color=velocity_interpretation | hard two-root palette | report.html"
        ),
    }
    return artifact


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    artifact = build_portable_artifact(args.root.resolve())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(artifact, indent=2, ensure_ascii=True) + "\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
