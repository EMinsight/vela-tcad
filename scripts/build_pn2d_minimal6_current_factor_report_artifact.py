#!/usr/bin/env python3
"""Build a portable report artifact for the Minimal6 current-factor audit."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


TITLE = "Minimal6 current influence-factor diagnosis"


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _number(row: dict[str, str], name: str) -> float:
    return float(row[name])


def build_artifact(root: Path) -> dict[str, object]:
    current = _read_csv(root / "current_factor_summary.csv")
    qfp = _read_csv(root / "qfp_mode_summary.csv")
    mobility = _read_csv(root / "mobility_residual_effect.csv")
    field = _read_csv(root / "driving_force_summary.csv")

    current_stages: list[dict[str, object]] = []
    current_table: list[dict[str, object]] = []
    for row in current:
        carrier = row["carrier"]
        before = _number(row, "self_consistent_median_error_dex")
        after = _number(row, "imported_state_median_error_dex")
        current_stages.extend(
            [
                {
                    "carrier": carrier,
                    "stage": "self-consistent Vela state",
                    "median_error_dex": before,
                    "paired_state_improvement_dex": _number(
                        row, "paired_state_improvement_median_dex"
                    ),
                    "sign_agreement": _number(
                        row, "self_consistent_sign_agreement_fraction"
                    ),
                    "active_edges": int(row["active_edge_count"]),
                },
                {
                    "carrier": carrier,
                    "stage": "imported Sentaurus state",
                    "median_error_dex": after,
                    "paired_state_improvement_dex": _number(
                        row, "paired_state_improvement_median_dex"
                    ),
                    "sign_agreement": _number(
                        row, "imported_state_sign_agreement_fraction"
                    ),
                    "active_edges": int(row["active_edge_count"]),
                },
            ]
        )
        current_table.append(
            {
                "carrier": carrier,
                "self_consistent_error_dex": before,
                "imported_state_error_dex": after,
                "paired_improvement_dex": _number(
                    row, "paired_state_improvement_median_dex"
                ),
                "current_weighted_improvement_dex": _number(
                    row, "current_weighted_state_improvement_dex"
                ),
                "self_consistent_sign_agreement": _number(
                    row, "self_consistent_sign_agreement_fraction"
                ),
                "imported_sign_agreement": _number(
                    row, "imported_state_sign_agreement_fraction"
                ),
            }
        )

    qfp_rows = [
        {
            "carrier": row["carrier"],
            "states": int(row["state_count"]),
            "common_mode_median_V": _number(
                row, "common_mode_abs_error_median_V"
            ),
            "differential_mode_median_V": _number(
                row, "differential_mode_abs_error_median_V"
            ),
            "sent_edge_delta_median_V": _number(
                row, "sentaurus_abs_edge_delta_median_V"
            ),
            "vela_edge_delta_median_V": _number(
                row, "vela_abs_edge_delta_median_V"
            ),
            "edge_sign_agreement": _number(
                row, "edge_qfp_sign_agreement_fraction"
            ),
        }
        for row in qfp
    ]
    mobility_rows = [
        {
            "carrier": row["carrier"],
            "samples": int(row["node_state_count"]),
            "sent_to_vela_residual_ratio": _number(
                row, "sentaurus_to_vela_residual_ratio_median"
            ),
            "residual_reduction_fraction": _number(
                row, "residual_reduction_fraction_median"
            ),
            "minimum_ratio": _number(row, "residual_ratio_minimum"),
            "maximum_ratio": _number(row, "residual_ratio_maximum"),
        }
        for row in mobility
    ]
    field_rows = [
        {
            "carrier": row["carrier"],
            "elements": int(row["element_count"]),
            "native_egrad_median_dex": _number(
                row, "native_egrad_mobility_error_median_dex"
            ),
            "native_egrad_p95_dex": _number(
                row, "native_egrad_mobility_error_p95_dex"
            ),
            "triangle_qfp_median_dex": _number(
                row, "triangle_qfp_mobility_error_median_dex"
            ),
            "triangle_qfp_p95_dex": _number(
                row, "triangle_qfp_mobility_error_p95_dex"
            ),
            "triangle_to_native_field_ratio": _number(
                row, "triangle_to_native_field_ratio_median"
            ),
        }
        for row in field
    ]

    sources = [
        {
            "id": "diagnosis_report",
            "label": "Current-factor diagnostic report",
            "path": (
                "build-release/pn2d-minimal6-current-factor-followup-"
                "20260724-a/report.md"
            ),
        },
        {
            "id": "current_summary",
            "label": "Paired edge-current factor summary",
            "path": (
                "build-release/pn2d-minimal6-current-factor-followup-"
                "20260724-a/current_factor_summary.csv"
            ),
            "query": {
                "engine": "python",
                "language": "python",
                "description": (
                    "Pair Phase F self-consistent edge-current errors with the "
                    "same 400 imported-state control edges."
                ),
                "tables_used": [
                    "directed_edge_current_comparison.csv",
                    "stage_edge_samples.csv",
                ],
                "filters": [
                    "mirror and sketch",
                    "bias -1 V through -20 V",
                    "valid nonzero carrier edges",
                ],
                "metric_definitions": {
                    "median_error_dex": (
                        "median abs(log10(abs(candidate/reference)))"
                    ),
                    "paired_state_improvement_dex": (
                        "self-consistent error minus imported-state error on "
                        "the same carrier edge"
                    ),
                },
            },
        },
        {
            "id": "qfp_summary",
            "label": "Internal-node QFP mode summary",
            "path": (
                "build-release/pn2d-minimal6-current-factor-followup-"
                "20260724-a/qfp_mode_summary.csv"
            ),
        },
        {
            "id": "mobility_summary",
            "label": "Fixed-state mobility residual effect",
            "path": (
                "build-release/pn2d-minimal6-current-factor-followup-"
                "20260724-a/mobility_residual_effect.csv"
            ),
        },
        {
            "id": "field_summary",
            "label": "Native and triangle driving-force comparison",
            "path": (
                "build-release/pn2d-minimal6-current-factor-followup-"
                "20260724-a/driving_force_summary.csv"
            ),
        },
    ]

    return {
        "surface": "report",
        "manifest": {
            "version": 1,
            "surface": "report",
            "title": TITLE,
            "generatedAt": "2026-07-24T00:00:00+08:00",
            "cards": [],
            "charts": [
                {
                    "id": "current_stage_chart",
                    "title": "Directed edge-current error by state branch",
                    "description": (
                        "Median absolute log10 error on 200 nonzero edges per "
                        "carrier; lower is better"
                    ),
                    "type": "bar",
                    "dataset": "current_stages",
                    "sourceId": "current_summary",
                    "source": {"label": "Paired edge-current factor summary", "path": "build-release/pn2d-minimal6-current-factor-followup-20260724-a/current_factor_summary.csv", "query": {"engine": "duckdb", "language": "sql", "description": "Reshape the paired carrier rows for the chart.", "sql": "SELECT carrier, 'self-consistent Vela state' AS stage, self_consistent_median_error_dex AS median_error_dex, paired_state_improvement_median_dex, self_consistent_sign_agreement_fraction AS sign_agreement, active_edge_count AS active_edges FROM read_csv_auto('current_factor_summary.csv', header = true) UNION ALL SELECT carrier, 'imported Sentaurus state' AS stage, imported_state_median_error_dex AS median_error_dex, paired_state_improvement_median_dex, imported_state_sign_agreement_fraction AS sign_agreement, active_edge_count AS active_edges FROM read_csv_auto('current_factor_summary.csv', header = true)", "tables_used": ["current_factor_summary.csv"]}},
                    "encodings": {
                        "x": {"field": "stage", "type": "nominal"},
                        "y": {
                            "field": "median_error_dex",
                            "type": "quantitative",
                        },
                        "color": {"field": "carrier", "type": "nominal"},
                    },
                }
            ],
            "tables": [
                {
                    "id": "current_table",
                    "title": "Paired state contribution",
                    "description": (
                        "All 40 states and five nonzero reference edges per "
                        "carrier"
                    ),
                    "dataset": "current_table",
                    "sourceId": "current_summary",
                    "source": {
                        "label": "Paired edge-current factor summary",
                        "path": "build-release/pn2d-minimal6-current-factor-followup-20260724-a/current_factor_summary.csv",
                        "query": {
                            "engine": "duckdb", "language": "sql", "description": "Select the paired current metrics shown in the table.", "sql": "SELECT carrier, self_consistent_median_error_dex AS self_consistent_error_dex, imported_state_median_error_dex AS imported_state_error_dex, paired_state_improvement_median_dex AS paired_improvement_dex, current_weighted_state_improvement_dex AS current_weighted_improvement_dex, self_consistent_sign_agreement_fraction AS self_consistent_sign_agreement, imported_state_sign_agreement_fraction AS imported_sign_agreement FROM read_csv_auto('current_factor_summary.csv', header = true)", "tables_used": ["current_factor_summary.csv"]
                        },
                    },
                    "columns": [
                        {"field": "carrier", "label": "Carrier"},
                        {
                            "field": "self_consistent_error_dex",
                            "label": "Self-consistent median (dex)",
                            "format": "number",
                        },
                        {
                            "field": "imported_state_error_dex",
                            "label": "Imported-state median (dex)",
                            "format": "number",
                        },
                        {
                            "field": "paired_improvement_dex",
                            "label": "Paired improvement (dex)",
                            "format": "number",
                        },
                        {
                            "field": "current_weighted_improvement_dex",
                            "label": "Current-weighted improvement (dex)",
                            "format": "number",
                        },
                        {
                            "field": "self_consistent_sign_agreement",
                            "label": "Self-consistent sign",
                            "format": "percent",
                        },
                        {
                            "field": "imported_sign_agreement",
                            "label": "Imported sign",
                            "format": "percent",
                        },
                    ],
                    "defaultSort": {"field": "carrier", "direction": "asc"},
                },
                {
                    "source": {
                        "label": "Internal-node QFP mode summary",
                        "path": "build-release/pn2d-minimal6-current-factor-followup-20260724-a/qfp_mode_summary.csv",
                        "query": {
                            "engine": "duckdb", "language": "sql", "description": "Select the QFP common- and differential-mode metrics.", "sql": "SELECT carrier, state_count AS states, common_mode_abs_error_median_V AS common_mode_median_V, differential_mode_abs_error_median_V AS differential_mode_median_V, sentaurus_abs_edge_delta_median_V AS sent_edge_delta_median_V, vela_abs_edge_delta_median_V AS vela_edge_delta_median_V, edge_qfp_sign_agreement_fraction AS edge_sign_agreement FROM read_csv_auto('qfp_mode_summary.csv', header = true)", "tables_used": ["qfp_mode_summary.csv"]
                        },
                    },
                    "id": "qfp_table",
                    "title": "QFP common and differential modes",
                    "description": (
                        "Internal nodes 1 and 5 over 40 exact states"
                    ),
                    "dataset": "qfp_modes",
                    "sourceId": "qfp_summary",
                    "columns": [
                        {"field": "carrier", "label": "Carrier"},
                        {"field": "states", "label": "States", "format": "number"},
                        {
                            "field": "common_mode_median_V",
                            "label": "Common mode (V)",
                            "format": "number",
                        },
                        {
                            "field": "differential_mode_median_V",
                            "label": "Differential mode (V)",
                            "format": "number",
                        },
                        {
                            "field": "sent_edge_delta_median_V",
                            "label": "Sent edge delta (V)",
                            "format": "number",
                        },
                        {
                            "field": "vela_edge_delta_median_V",
                            "label": "Vela edge delta (V)",
                            "format": "number",
                        },
                        {
                            "field": "edge_sign_agreement",
                            "label": "Edge sign",
                            "format": "percent",
                        },
                    ],
                    "defaultSort": {"field": "carrier", "direction": "asc"},
                },
                {
                    "id": "field_table",
                    "title": "Mobility replay under two QFP-field definitions",
                    "description": (
                        "160 native elements per carrier; same documented "
                        "mobility parameters"
                    ),
                    "dataset": "field_summary",
                    "sourceId": "field_summary",
                    "source": {
                        "label": "Native and triangle driving-force comparison",
                        "path": "build-release/pn2d-minimal6-current-factor-followup-20260724-a/driving_force_summary.csv",
                        "query": {
                            "engine": "duckdb", "language": "sql", "description": "Select the two mobility-replay field definitions.", "sql": "SELECT carrier, element_count AS elements, native_egrad_mobility_error_median_dex AS native_egrad_median_dex, native_egrad_mobility_error_p95_dex AS native_egrad_p95_dex, triangle_qfp_mobility_error_median_dex AS triangle_qfp_median_dex, triangle_qfp_mobility_error_p95_dex AS triangle_qfp_p95_dex, triangle_to_native_field_ratio_median AS triangle_to_native_field_ratio FROM read_csv_auto('driving_force_summary.csv', header = true)", "tables_used": ["driving_force_summary.csv"]
                        },
                    },
                    "columns": [
                        {"field": "carrier", "label": "Carrier"},
                        {
                            "field": "elements",
                            "label": "Elements",
                            "format": "number",
                        },
                        {
                            "field": "native_egrad_median_dex",
                            "label": "Native eGrad median (dex)",
                            "format": "number",
                        },
                        {
                            "field": "native_egrad_p95_dex",
                            "label": "Native eGrad P95 (dex)",
                            "format": "number",
                        },
                        {
                            "field": "triangle_qfp_median_dex",
                            "label": "Triangle QFP median (dex)",
                            "format": "number",
                        },
                        {
                            "field": "triangle_qfp_p95_dex",
                            "label": "Triangle QFP P95 (dex)",
                            "format": "number",
                        },
                        {
                            "field": "triangle_to_native_field_ratio",
                            "label": "Triangle/native field",
                            "format": "number",
                        },
                    ],
                    "defaultSort": {"field": "carrier", "direction": "asc"},
                },
            ],
            "blocks": [
                {"id": "title", "type": "markdown", "body": f"# {TITLE}"},
                {
                    "id": "technical_summary",
                    "type": "markdown",
                    "body": (
                        "## Technical summary\n\n"
                        "The self-consistent QFP state is the dominant directed-"
                        "current driver. Importing the Sentaurus potentials and "
                        "recomputing dependent Vela quantities removes 0.72 dex "
                        "of paired electron-edge error and 0.85 dex of paired "
                        "hole-edge error. The remaining fixed-state median is "
                        "about 0.06 dex and is attributable to mobility/operator "
                        "support rather than density or geometry."
                    ),
                    "sourceId": "diagnosis_report",
                },
                {
                    "id": "state_finding",
                    "type": "markdown",
                    "body": (
                        "## QFP state explains most of the current gap\n\n"
                        "The state replacement improves both unweighted and "
                        "absolute-current-weighted errors. All imported-state "
                        "nonzero edges have matching signs, whereas the "
                        "self-consistent branch has 80% sign agreement."
                    ),
                    "sourceId": "current_summary",
                },
                {
                    "id": "current_chart_block",
                    "type": "chart",
                    "chartId": "current_stage_chart",
                },
                {
                    "id": "current_table_block",
                    "type": "table",
                    "tableId": "current_table",
                },
                {
                    "id": "qfp_finding",
                    "type": "markdown",
                    "body": (
                        "## Common-mode QFP controls density; differential mode "
                        "controls the central-edge sign\n\n"
                        "The common-mode offsets are 0.0549 V for electrons and "
                        "0.0621 V for holes. The approximately 0.53 mV "
                        "differential mismatch is much smaller, but the physical "
                        "Sentaurus 1-5 drop is also about 0.5 mV, so Vela reverses "
                        "that local direction in every state."
                    ),
                    "sourceId": "qfp_summary",
                },
                {
                    "id": "qfp_table_block",
                    "type": "table",
                    "tableId": "qfp_table",
                },
                {
                    "id": "mobility_finding",
                    "type": "markdown",
                    "body": (
                        "## Mobility contributes but does not close the "
                        "continuity residual\n\n"
                        "On the identical imported state, coefficient-weighted "
                        "Sentaurus element mobility reduces the median residual "
                        "by 15.2% for electrons and 4.3% for holes. This is too "
                        "small to explain the final QFP displacement alone."
                    ),
                    "sourceId": "mobility_summary",
                },
                {
                    "id": "field_finding",
                    "type": "markdown",
                    "body": (
                        "## Exported electron eGradQuasiFermi is not a proven "
                        "high-field mobility drive\n\n"
                        "Using the exported native electron field gives 0.594 "
                        "dex median mobility error; using the affine node-QFP "
                        "gradient gives 0.0527 dex. Their field ratio is "
                        "7.633856879 across all 160 electron elements. This "
                        "supports an output-semantics or proprietary element "
                        "evaluation difference, not a production scale change."
                    ),
                    "sourceId": "field_summary",
                },
                {
                    "id": "field_table_block",
                    "type": "table",
                    "tableId": "field_table",
                },
                {
                    "id": "scope",
                    "type": "markdown",
                    "body": (
                        "## Scope, data, and metric definitions\n\n"
                        "The comparison covers mirror/sketch at -1 through -20 "
                        "V, 400 nonzero carrier edges, 80 internal-node states "
                        "per carrier, and 160 native elements per carrier. "
                        "Current error is abs(log10(abs(candidate/reference))). "
                        "All pairwise effects use the same topology, bias, "
                        "carrier, and unordered node pair."
                    ),
                },
                {
                    "id": "method",
                    "type": "markdown",
                    "body": (
                        "## Methodology and robustness\n\n"
                        "Phase F self-consistent edge errors are paired with the "
                        "Phase C imported-state control. QFP errors are split "
                        "into the mean of nodes 1/5 and their difference. "
                        "Mobility residual effects come from Phase E fixed-state "
                        "branches. Two generated roots are byte-identical."
                    ),
                },
                {
                    "id": "limitations",
                    "type": "markdown",
                    "body": (
                        "## Limitations and uncertainty\n\n"
                        "The Sentaurus edge reference is a terminal- and total-"
                        "KCL-closed box-operator reconstruction, not a native "
                        "directed-edge observation. The authorized Sentaurus VM "
                        "was unreachable, so a high-field-off native low-field "
                        "mobility probe remains pending."
                    ),
                },
                {
                    "id": "next",
                    "type": "markdown",
                    "body": (
                        "## Recommended next step\n\n"
                        "Regenerate one identical Sentaurus state with "
                        "HighFieldSaturation disabled. Pair its native low-field "
                        "element mobility with the current final mobility to "
                        "invert the actual high-field drive without assuming "
                        "Vela low-field interpolation. Do not change production "
                        "mobility, SG, Poisson, impact, or QFP formulas before "
                        "that control."
                    ),
                },
                {
                    "id": "question",
                    "type": "markdown",
                    "body": (
                        "## Further question\n\n"
                        "Does the native low-field control show that Sentaurus "
                        "internally drives mobility with the affine QFP gradient, "
                        "the exported eGrad field, or a third boundary-smoothed "
                        "field?"
                    ),
                },
            ],
            "sources": sources,
        },
        "snapshot": {
            "version": 1,
            "generatedAt": "2026-07-24T00:00:00+08:00",
            "status": "ready",
            "datasets": {
                "current_stages": current_stages,
                "current_table": current_table,
                "qfp_modes": qfp_rows,
                "mobility_effect": mobility_rows,
                "field_summary": field_rows,
            },
        },
        "sources": sources,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
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
