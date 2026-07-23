#!/usr/bin/env python3
"""Build the canonical portable-report artifact for the mobility diagnosis."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


TITLE = "Minimal6 mobility and quasi-Fermi-gradient diagnosis"


def _round(value, digits=12):
    return None if value is None else round(float(value), digits)


def _combined(rows):
    return [row for row in rows if row["topology"] == "combined"]


def build_artifact(report_path: Path) -> dict:
    report = json.loads(report_path.read_text(encoding="utf-8"))
    conclusions = report["conclusions"]
    mobility_rows = []
    for row in _combined(report["mobility_summary"]):
        mobility_rows.append({
            "carrier": row["carrier"],
            "branch": row["branch"],
            "n": row["sample_count"],
            "median_mobility": _round(row["median_mobility_m2_per_Vs"]),
            "median_abs_dex": _round(row["median_abs_log10_error_vs_sentaurus"]),
            "p95_abs_dex": _round(row["p95_abs_log10_error_vs_sentaurus"]),
            "median_relative_error": _round(row["median_relative_error_vs_sentaurus"]),
        })
    qf_rows = []
    for row in _combined(report["qf_gradient_summary"]):
        qf_rows.append({
            "support": row["support"],
            "carrier": row["carrier"],
            "branch": row["branch"],
            "n_valid": row["magnitude_valid_count"],
            "median_abs_dex": _round(row["median_abs_log10_error"]),
            "p95_abs_dex": _round(row["p95_abs_log10_error"]),
            "median_angle_deg": _round(row["median_angle_deg"]),
            "sign_agreement": _round(row["sign_agreement_fraction"]),
        })
    orientation_rows = [{
        "carrier": row["carrier"],
        "transform": row["transform"],
        "n": row["valid_count"],
        "median_angle_deg": _round(row["median_angle_deg"]),
        "p95_angle_deg": _round(row["p95_angle_deg"]),
    } for row in report["orientation_control_summary"]]

    sources = [
        {
            "id": "diagnosis_output",
            "label": "Deterministic mobility diagnosis",
            "path": "build-release/pn2d-minimal6-mobility-diagnosis-20260723-b/mobility_diagnosis.json",
        },
        {
            "id": "sealed_inputs",
            "label": "Sealed Minimal6 inverse inputs",
            "path": "build-release/pn2d-minimal6-inverse-inputs-20260722-a",
        },
        {
            "id": "production_mobility",
            "label": "Vela production mobility implementation",
            "path": "src/physics/MobilityModel.cpp",
        },
        {
            "id": "mobility_summary_query",
            "label": "Combined edge-mobility summary",
            "path": "build-release/pn2d-minimal6-mobility-diagnosis-20260723-b/mobility_summary.csv",
            "query": {
                "engine": "duckdb",
                "language": "sql",
                "description": "Select combined-topology mobility branches for the chart and exact table.",
                "sql": "SELECT carrier, branch, sample_count AS n, median_mobility_m2_per_Vs AS median_mobility, median_abs_log10_error_vs_sentaurus AS median_abs_dex, p95_abs_log10_error_vs_sentaurus AS p95_abs_dex, median_relative_error_vs_sentaurus AS median_relative_error FROM read_csv_auto('mobility_summary.csv', header = true) WHERE topology = 'combined'",
            },
        },
        {
            "id": "qf_summary_query",
            "label": "Combined quasi-Fermi-gradient summary",
            "path": "build-release/pn2d-minimal6-mobility-diagnosis-20260723-b/qf_gradient_summary.csv",
            "query": {
                "engine": "duckdb",
                "language": "sql",
                "description": "Select combined-topology quasi-Fermi-gradient branches for the exact table.",
                "sql": "SELECT support, carrier, branch, magnitude_valid_count AS n_valid, median_abs_log10_error AS median_abs_dex, p95_abs_log10_error AS p95_abs_dex, median_angle_deg, sign_agreement_fraction AS sign_agreement FROM read_csv_auto('qf_gradient_summary.csv', header = true) WHERE topology = 'combined'",
            },
        },
        {
            "id": "orientation_summary_query",
            "label": "Signed-axis orientation controls",
            "path": "build-release/pn2d-minimal6-mobility-diagnosis-20260723-b/orientation_control_summary.csv",
            "query": {
                "engine": "duckdb",
                "language": "sql",
                "description": "Select all fixed signed-axis orientation controls for the exact table.",
                "sql": "SELECT carrier, transform, valid_count AS n, median_angle_deg, p95_angle_deg FROM read_csv_auto('orientation_control_summary.csv', header = true)",
            },
        },
    ]
    summary = (
        "## Technical summary\n\n"
        f"The earlier near-90-degree aggregate is a sign-pooling artifact: its reproduced median is "
        f"{conclusions['legacy_mixed_sign_pooled_median_angle_deg']:.6f} degrees. Using the data-supported "
        f"negative sign for both exported carrier quasi-Fermi potentials reduces the pooled median to "
        f"{conclusions['corrected_same_sign_pooled_median_angle_deg']:.6f} degrees; identity is the best "
        "coordinate transform for both carriers. Mobility changes magnitude, not this direction result."
    )
    return {
        "surface": "report",
        "manifest": {
            "version": 1,
            "surface": "report",
            "title": TITLE,
            "generatedAt": "2026-07-23T00:00:00+08:00",
            "cards": [],
            "charts": [
                {
                    "id": "mobility_error_chart",
                    "title": "Median edge-mobility error by branch",
                    "description": "40 exact states; 360 edge samples per carrier and branch; absolute log10 error versus exported Sentaurus mobility",
                    "type": "bar",
                    "dataset": "mobility_summary",
                    "sourceId": "mobility_summary_query",
                    "encodings": {
                        "x": {"field": "branch", "type": "nominal"},
                        "y": {"field": "median_abs_dex", "type": "quantitative"},
                        "color": {"field": "carrier", "type": "nominal"},
                    },
                },
            ],
            "tables": [
                {
                    "id": "mobility_table",
                    "title": "Edge mobility comparison",
                    "description": "40 exact states; 360 edge samples per carrier and branch",
                    "dataset": "mobility_summary",
                    "sourceId": "mobility_summary_query",
                    "columns": [
                        {"field": "carrier", "label": "Carrier"},
                        {"field": "branch", "label": "Mobility branch"},
                        {"field": "n", "label": "N", "format": "number"},
                        {"field": "median_mobility", "label": "Median mobility (m2/V/s)", "format": "number"},
                        {"field": "median_abs_dex", "label": "Median abs error (dex)", "format": "number"},
                        {"field": "p95_abs_dex", "label": "P95 abs error (dex)", "format": "number"},
                        {"field": "median_relative_error", "label": "Median relative error", "format": "number"},
                    ],
                    "defaultSort": {"field": "carrier", "direction": "asc"},
                },
                {
                    "id": "qf_table",
                    "title": "Quasi-Fermi-gradient inversion",
                    "description": "Node, edge, and cell support with both projection orders",
                    "dataset": "qf_summary",
                    "sourceId": "qf_summary_query",
                    "columns": [
                        {"field": "support", "label": "Support/order"},
                        {"field": "carrier", "label": "Carrier"},
                        {"field": "branch", "label": "Mobility branch"},
                        {"field": "n_valid", "label": "N valid", "format": "number"},
                        {"field": "median_abs_dex", "label": "Median abs error (dex)", "format": "number"},
                        {"field": "p95_abs_dex", "label": "P95 abs error (dex)", "format": "number"},
                        {"field": "median_angle_deg", "label": "Median angle (deg)", "format": "number"},
                        {"field": "sign_agreement", "label": "Sign agreement", "format": "number"},
                    ],
                    "defaultSort": {"field": "support", "direction": "asc"},
                },
                {
                    "id": "orientation_table",
                    "title": "Sign and coordinate controls",
                    "description": "Eight fixed signed-axis transforms on locally inverted cell vectors",
                    "dataset": "orientation_summary",
                    "sourceId": "orientation_summary_query",
                    "columns": [
                        {"field": "carrier", "label": "Carrier"},
                        {"field": "transform", "label": "Transform"},
                        {"field": "n", "label": "N", "format": "number"},
                        {"field": "median_angle_deg", "label": "Median angle (deg)", "format": "number"},
                        {"field": "p95_angle_deg", "label": "P95 angle (deg)", "format": "number"},
                    ],
                    "defaultSort": {"field": "median_angle_deg", "direction": "asc"},
                },
            ],
            "blocks": [
                {"id": "title", "type": "markdown", "body": f"# {TITLE}"},
                {"id": "technical_summary", "type": "markdown", "body": summary, "sourceId": "diagnosis_output"},
                {"id": "mobility_finding", "type": "markdown", "body": "## Same-state Masetti is closer than the constant control, but tails remain large\n\nVela Masetti evaluated on the Sentaurus state has median edge-mobility errors of 0.102586 dex for electrons and 0.086174 dex for holes. Constant mobility is much farther away. The near-1 dex electron P95 shows that one global scale cannot replace the local high-field law.", "sourceId": "diagnosis_output"},
                {"id": "mobility_chart_block", "type": "chart", "chartId": "mobility_error_chart"},
                {"id": "mobility_block", "type": "table", "tableId": "mobility_table"},
                {"id": "qf_finding", "type": "markdown", "body": "## Projection order, not mobility, dominates the magnitude residual\n\nWith exported mobility, local inversion followed by projection gives median cell errors of 0.142764 dex (electron) and 0.149314 dex (hole). Projecting current, density, and mobility first and then dividing inflates the errors to 5.779672 and 5.687322 dex because division and averaging do not commute across the strongly nonuniform junction.", "sourceId": "diagnosis_output"},
                {"id": "qf_block", "type": "table", "tableId": "qf_table"},
                {"id": "orientation_finding", "type": "markdown", "body": "## The exported hole-QFP sign caused the old 90-degree aggregate\n\nIdentity gives median cell angles of 0.009142 degrees for electrons and 0.011362 degrees for holes. Negating either vector gives approximately 180 degrees; swapping axes gives 90 degrees. This supports the coordinate mapping and a negative current-to-QFP-gradient sign for both exported carrier potentials.", "sourceId": "diagnosis_output"},
                {"id": "orientation_block", "type": "table", "tableId": "orientation_table"},
                {"id": "scope", "type": "markdown", "body": "## Scope, data, and definitions\n\nThe analysis covers 40 exact states, two topologies, two carriers, 720 carrier-edge rows, and 320 carrier-cell rows. Edge support uses endpoint means and tangent current projection. Cell support uses equal P1 node means and the affine triangle QFP gradient. All source quantities are sealed SI conversions; no production formula was modified.", "sourceId": "diagnosis_output"},
                {"id": "method", "type": "markdown", "body": "## Methodology and robustness\n\nThe diagnostic reproduces Vela's Masetti low-field formula, endpoint-average net doping, and QFP-gradient high-field limiter. It compares exported Sentaurus mobility, same-state Vela Masetti, and the 1417/470.5 cm2/V/s constant control. A second deterministic root is byte-identical. The chart communicates the branch-level mobility comparison; tables retain exact branch/support lookup."},
                {"id": "limitations", "type": "markdown", "body": "## Limitations and uncertainty\n\nSentaurus current is node-exported, not an internal edge flux. The local J/(q n mu) identity is not a discrete Scharfetter-Gummel inverse. Eighty of 360 edge rows per carrier have zero reference tangent gradient and are excluded from magnitude summaries. Large P95 node errors are support-zero/outlier sensitive and should not be interpreted as a mobility fit target."},
                {"id": "next_steps", "type": "markdown", "body": "## Recommended next step\n\nDo not change production mobility or current formulas from this diagnostic. Export or reconstruct a Sentaurus edge flux/current on the same directed edges, then invert the discrete Scharfetter-Gummel relation using endpoint densities and QFP differences. That is the next test capable of separating current semantics from mobility magnitude."},
                {"id": "further_questions", "type": "markdown", "body": "## Further question\n\nCan the discrete edge-flux inversion close the remaining approximately 0.13-0.15 dex median magnitude gap without an empirical mobility scale?"},
            ],
            "sources": sources,
        },
        "snapshot": {
            "version": 1,
            "generatedAt": "2026-07-23T00:00:00+08:00",
            "status": "ready",
            "datasets": {
                "mobility_summary": mobility_rows,
                "qf_summary": qf_rows,
                "orientation_summary": orientation_rows,
            },
        },
        "sources": sources,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    artifact = build_artifact(args.report)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(artifact, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
