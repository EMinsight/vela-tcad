#!/usr/bin/env python3
"""Build a portable report artifact for the M2 BV knee diagnostic."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


TITLE = "PN2D M2 BV knee-region read-only error localization"


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def number(value: str) -> float:
    return float(value)


def build(root: Path) -> dict[str, object]:
    result = json.loads((root / "diagnostic.json").read_text(encoding="utf-8"))
    bias_rows = rows(root / "bias_summary.csv")
    edge_rows = rows(root / "edge_current_summary.csv")
    mapping_rows = rows(root / "source_mapping_summary.csv")
    tracking: list[dict[str, object]] = []
    stages: list[dict[str, object]] = []
    for row in bias_rows:
        voltage = number(row["reverse_bias_magnitude_V"])
        tracking.extend(
            [
                {"reverse_bias_V": voltage, "metric": "Terminal current", "error_dex": number(row["terminal_abs_log_error_dex"])},
                {"reverse_bias_V": voltage, "metric": "Integrated source", "error_dex": number(row["source_abs_log_error_dex"])},
                {"reverse_bias_V": voltage, "metric": "Electron density", "error_dex": abs(number(row["electron_density_log_ratio_dex"]))},
                {"reverse_bias_V": voltage, "metric": "Hole density", "error_dex": abs(number(row["hole_density_log_ratio_dex"]))},
            ]
        )
        for carrier, label in (("electron", "Electron"), ("hole", "Hole")):
            for field, stage in (
                ("current_log_ratio_dex", "Cell current"),
                ("drive_log_ratio_dex", "QFP drive"),
                ("mobility_log_ratio_dex", "Dominant-edge mobility"),
                ("alpha_log_ratio_dex", "Alpha counterfactual"),
            ):
                stages.append(
                    {
                        "reverse_bias_V": voltage,
                        "series": f"{label} {stage}",
                        "signed_log_ratio_dex": number(row[f"{carrier}_{field}"]),
                    }
                )
    edge = [
        {
            "reverse_bias_V": number(row["reverse_bias_magnitude_V"]),
            "native_projection_median_dex": number(row["native_projection_median_abs_error_dex"]),
            "replay_scaled_median_dex": number(row["operator_replay_x1e6_median_abs_error_dex"]),
            "replay_raw_median_dex": number(row["operator_replay_raw_median_abs_error_dex"]),
        }
        for row in edge_rows
    ]
    mapping = [
        {
            "reverse_bias_V": number(row["reverse_bias_magnitude_V"]),
            "vertex_overlap": number(row["vertex_source_overlap"]),
            "cell_overlap": number(row["cell_source_overlap"]),
            "source_measure_max_rel": number(row["source_measure_max_relative_error"]),
        }
        for row in mapping_rows
    ]
    sources = [
        {
            "id": "diagnostic",
            "label": "M2 knee machine-readable diagnostic",
            "path": "build-release/pn2d-bv-m2-knee-readonly-diagnostic-20260731/diagnostic.json",
        },
        {
            "id": "bias_summary",
            "label": "Bias-level process decomposition",
            "path": "build-release/pn2d-bv-m2-knee-readonly-diagnostic-20260731/bias_summary.csv",
        },
        {
            "id": "sentaurus_manifest",
            "label": "Sentaurus M2 process manifest",
            "path": "build-release/pn2d-task10-balanced-m2-sentaurus-process-v2-20260731/manifest.json",
        },
        {
            "id": "vela_manifest",
            "label": "Vela M2 process manifest",
            "path": "build-release/pn2d-bv-template-default-prospective-v2-default-20260731/M2/run-a/manifest.json",
        },
        {
            "id": "tracking_query",
            "label": "Knee error tracking query",
            "path": "build-release/pn2d-bv-m2-knee-readonly-diagnostic-20260731/bias_summary.csv",
            "query": {
                "engine": "duckdb",
                "language": "sql",
                "description": "Reshape terminal, source, and carrier-density errors over the frozen M2 knee lattice.",
                "sql": "SELECT reverse_bias_magnitude_V AS reverse_bias_V, 'Terminal current' AS metric, terminal_abs_log_error_dex AS error_dex FROM read_csv_auto('bias_summary.csv', header=true) UNION ALL SELECT reverse_bias_magnitude_V, 'Integrated source', source_abs_log_error_dex FROM read_csv_auto('bias_summary.csv', header=true) UNION ALL SELECT reverse_bias_magnitude_V, 'Electron density', abs(electron_density_log_ratio_dex) FROM read_csv_auto('bias_summary.csv', header=true) UNION ALL SELECT reverse_bias_magnitude_V, 'Hole density', abs(hole_density_log_ratio_dex) FROM read_csv_auto('bias_summary.csv', header=true)",
                "tables_used": ["bias_summary.csv"],
                "filters": ["-20 V <= bias_V <= -18 V", "branch = avalanche_on"],
            },
        },
        {
            "id": "stage_query",
            "label": "Local-stage ratio query",
            "path": "build-release/pn2d-bv-m2-knee-readonly-diagnostic-20260731/bias_summary.csv",
            "query": {
                "engine": "duckdb",
                "language": "sql",
                "description": "Reshape carrier current, drive, mobility, and alpha signed log ratios.",
                "sql": "SELECT reverse_bias_magnitude_V AS reverse_bias_V, 'Electron Cell current' AS series, electron_current_log_ratio_dex AS signed_log_ratio_dex FROM read_csv_auto('bias_summary.csv', header=true) UNION ALL SELECT reverse_bias_magnitude_V, 'Hole Cell current', hole_current_log_ratio_dex FROM read_csv_auto('bias_summary.csv', header=true) UNION ALL SELECT reverse_bias_magnitude_V, 'Electron QFP drive', electron_drive_log_ratio_dex FROM read_csv_auto('bias_summary.csv', header=true) UNION ALL SELECT reverse_bias_magnitude_V, 'Hole QFP drive', hole_drive_log_ratio_dex FROM read_csv_auto('bias_summary.csv', header=true) UNION ALL SELECT reverse_bias_magnitude_V, 'Electron Dominant-edge mobility', electron_mobility_log_ratio_dex FROM read_csv_auto('bias_summary.csv', header=true) UNION ALL SELECT reverse_bias_magnitude_V, 'Hole Dominant-edge mobility', hole_mobility_log_ratio_dex FROM read_csv_auto('bias_summary.csv', header=true) UNION ALL SELECT reverse_bias_magnitude_V, 'Electron Alpha counterfactual', electron_alpha_log_ratio_dex FROM read_csv_auto('bias_summary.csv', header=true) UNION ALL SELECT reverse_bias_magnitude_V, 'Hole Alpha counterfactual', hole_alpha_log_ratio_dex FROM read_csv_auto('bias_summary.csv', header=true)",
                "tables_used": ["bias_summary.csv"],
                "filters": ["Sentaurus active-cell source >= 0.1% of bias-local peak"],
            },
        },
        {
            "id": "mapping_query",
            "label": "Source-map overlap query",
            "path": "build-release/pn2d-bv-m2-knee-readonly-diagnostic-20260731/source_mapping_summary.csv",
            "query": {
                "engine": "duckdb",
                "language": "sql",
                "description": "Reshape normalized source-map overlap at element-vertex and cell support.",
                "sql": "SELECT reverse_bias_magnitude_V AS reverse_bias_V, 'Element vertex' AS support, vertex_source_overlap AS overlap FROM read_csv_auto('source_mapping_summary.csv', header=true) UNION ALL SELECT reverse_bias_magnitude_V, 'Cell', cell_source_overlap FROM read_csv_auto('source_mapping_summary.csv', header=true)",
                "tables_used": ["source_mapping_summary.csv"],
                "filters": ["branch = avalanche_on"],
            },
        },
        {
            "id": "edge_query",
            "label": "Edge-current control query",
            "path": "build-release/pn2d-bv-m2-knee-readonly-diagnostic-20260731/edge_current_summary.csv",
            "query": {
                "engine": "duckdb",
                "language": "sql",
                "description": "Select connectivity-aligned native and replay edge-current errors.",
                "sql": "SELECT reverse_bias_magnitude_V AS reverse_bias_V, native_projection_median_abs_error_dex AS native_projection_median_dex, operator_replay_x1e6_median_abs_error_dex AS replay_scaled_median_dex, operator_replay_raw_median_abs_error_dex AS replay_raw_median_dex FROM read_csv_auto('edge_current_summary.csv', header=true)",
                "tables_used": ["edge_current_summary.csv"],
                "filters": ["edges paired by unordered endpoint connectivity"],
            },
        },
    ]
    growth = result["knee_error_growth"]
    controls = result["controls"]
    correlations = result["correlations"]
    return {
        "surface": "report",
        "manifest": {
            "version": 1,
            "surface": "report",
            "title": TITLE,
            "generatedAt": "2026-07-31T00:00:00+08:00",
            "cards": [],
            "charts": [
                {
                    "id": "tracking_chart",
                    "title": "Error growth over the M2 knee window",
                    "description": "Absolute log10 error versus Sentaurus; no acceptance threshold was changed.",
                    "type": "line",
                    "dataset": "tracking",
                    "sourceId": "tracking_query",
                    "encodings": {
                        "x": {"field": "reverse_bias_V", "type": "quantitative", "title": "Reverse-bias magnitude (V)"},
                        "y": {"field": "error_dex", "type": "quantitative", "title": "Absolute error (dex)"},
                        "color": {"field": "metric", "type": "nominal"},
                    },
                },
                {
                    "id": "stage_chart",
                    "title": "Signed local-stage ratios",
                    "description": "log10(Vela/Sentaurus) on active Sentaurus source cells; alpha uses a source-integral counterfactual.",
                    "type": "line",
                    "dataset": "stages",
                    "sourceId": "stage_query",
                    "encodings": {
                        "x": {"field": "reverse_bias_V", "type": "quantitative", "title": "Reverse-bias magnitude (V)"},
                        "y": {"field": "signed_log_ratio_dex", "type": "quantitative", "title": "log10(Vela/Sentaurus)"},
                        "color": {"field": "series", "type": "nominal"},
                    },
                },
                {
                    "id": "mapping_chart",
                    "title": "Normalized source-map overlap",
                    "description": "One means identical spatial allocation after normalizing each total source.",
                    "type": "line",
                    "dataset": "mapping_long",
                    "sourceId": "mapping_query",
                    "encodings": {
                        "x": {"field": "reverse_bias_V", "type": "quantitative", "title": "Reverse-bias magnitude (V)"},
                        "y": {"field": "overlap", "type": "quantitative", "title": "Overlap"},
                        "color": {"field": "support", "type": "nominal"},
                    },
                },
            ],
            "tables": [
                {
                    "id": "edge_table",
                    "title": "Edge-current controls",
                    "description": "Edges aligned by unordered node connectivity, not local-edge index.",
                    "dataset": "edge",
                    "sourceId": "edge_query",
                    "columns": [
                        {"field": "reverse_bias_V", "label": "|V|", "format": "number"},
                        {"field": "native_projection_median_dex", "label": "Native projection median (dex)", "format": "number"},
                        {"field": "replay_scaled_median_dex", "label": "Replay x1e6 median (dex)", "format": "number"},
                        {"field": "replay_raw_median_dex", "label": "Raw replay median (dex)", "format": "number"},
                    ],
                    "defaultSort": {"field": "reverse_bias_V", "direction": "asc"},
                }
            ],
            "blocks": [
                {"id": "title", "type": "markdown", "body": f"# {TITLE}"},
                {
                    "id": "summary",
                    "type": "markdown",
                    "sourceId": "diagnostic",
                    "body": (
                        "## Technical summary\n\n"
                        "The growing M2 knee discrepancy is localized to the self-consistent carrier-density and SG/Laux current-amplitude loop. From 18 to 20 V reverse-bias magnitude, terminal-current error grows by "
                        f"{growth['terminal_error_dex']:.5f} dex and integrated-source error grows by {growth['integrated_source_error_dex']:.5f} dex. "
                        "Their bias-wise correlation is "
                        f"{correlations['terminal_vs_integrated_source_error']:.8f}. This is reproducible rather than continuation noise: the two independent Vela runs are byte-identical for IV, process probes, Newton attempts, and Newton history."
                    ),
                },
                {"id": "tracking", "type": "chart", "chartId": "tracking_chart"},
                {
                    "id": "findings",
                    "type": "markdown",
                    "body": (
                        "## Key findings\n\n"
                        "Both electron and hole cell-current deficits grow by about 0.053 dex, while source-weighted carrier-density deficits grow by 0.051 to 0.056 dex. Those changes match the 0.055 dex terminal-error growth. In contrast, the maximum QFP-drive mismatch is "
                        f"{controls['maximum_qfp_drive_abs_log_ratio_dex']:.4g} dex, the source-weighted alpha counterfactual is within {controls['maximum_source_weighted_alpha_abs_log_ratio_dex']:.4g} dex, and source geometry measures match within {controls['maximum_source_measure_relative_error']:.3g} relative. The carrier source-fraction difference never exceeds {controls['maximum_carrier_source_fraction_difference']:.4g}."
                    ),
                },
                {"id": "stages", "type": "chart", "chartId": "stage_chart"},
                {
                    "id": "mapping_text",
                    "type": "markdown",
                    "body": (
                        "## Source mapping\n\n"
                        "The source-measure geometry is not the cause: every nonzero carrier support closes to the same element-vertex measure within 9.47e-6 relative, and the hotspot cell and vertex are identical at every bias. The normalized source shape is not identical (minimum overlap 0.836 at vertex support and 0.902 after cell aggregation), but overlap improves as the terminal error worsens. The residual is therefore a support/interpolation difference in local generation, not a growing geometric-weight or hotspot relocation error."
                    ),
                },
                {"id": "mapping", "type": "chart", "chartId": "mapping_chart"},
                {
                    "id": "edge_text",
                    "type": "markdown",
                    "body": (
                        "## Edge-current evidence\n\n"
                        "After aligning edges by node connectivity, the median Vela-versus-native-Sentaurus edge-current projection error grows from about 0.18 to 0.24 dex, mirroring the knee error. The Sentaurus operator-replay export carries a consistent 1e6 scale signature: raw comparison is about 6 dex off, while applying the inferred um2-to-cm2 conversion reduces the median to 0.10-0.12 dex. This scale issue belongs to the diagnostic export and is not evidence of a production physics error."
                    ),
                },
                {"id": "edge", "type": "table", "tableId": "edge_table"},
                {
                    "id": "scope",
                    "type": "markdown",
                    "body": (
                        "## Scope, data, and definitions\n\n"
                        "Scope is restricted to the M2 avalanche-on exact lattice from -18 through -20 V. Sentaurus is the golden reference. Signed stage values are log10(Vela/Sentaurus); terminal and source curves use absolute log error. Active cells are selected prospectively from Sentaurus total local source at or above 0.1% of the bias-local peak. No solver, model, template default, or acceptance threshold was changed."
                    ),
                },
                {
                    "id": "method",
                    "type": "markdown",
                    "body": (
                        "## Methodology and robustness\n\n"
                        "Cell fields are paired by cell ID on the shared M2 mesh; edges are paired by unordered endpoint connectivity; element-vertex sources are paired by exact support key. Alpha is evaluated as a source-weighted counterfactual ratio, not a median over inactive nodes. Source maps are checked both per vertex and after cell aggregation. Four independent Vela artifacts are hashed across run A and run B."
                    ),
                },
                {
                    "id": "limitations",
                    "type": "markdown",
                    "body": (
                        "## Limitations\n\n"
                        "Self-consistent snapshots establish where the growing discrepancy resides but cannot order the feedback loop. In particular, they do not prove whether the carrier-density deficit or the difference between Sentaurus native current support and Vela SG/Laux support is the initiating perturbation. The operator-replay scale correction is an explicit inference from the stable 1e6 signature and is reported both raw and corrected."
                    ),
                },
                {
                    "id": "next",
                    "type": "markdown",
                    "body": (
                        "## Recommended next read-only control\n\n"
                        "Replay the exact M2 Sentaurus states at -18, -19.5, -19.7, and -20 V through Vela SG/Laux in postprocess-only mode. Compare the resulting source integral with both the Sentaurus source and the self-consistent Vela source. This will separate the frozen operator from self-consistent state feedback without changing production defaults or acceptance thresholds."
                    ),
                },
                {
                    "id": "question",
                    "type": "markdown",
                    "body": (
                        "## Further question\n\n"
                        "On the imported Sentaurus M2 state, does SG/Laux close the total source while retaining the native-current edge projection gap? If yes, the remaining production decision should target state/current feedback semantics rather than alpha, mobility, or source geometry."
                    ),
                },
            ],
            "sources": sources,
        },
        "snapshot": {
            "version": 1,
            "generatedAt": "2026-07-31T00:00:00+08:00",
            "status": "ready",
            "datasets": {
                "tracking": tracking,
                "stages": stages,
                "edge": edge,
                "mapping": mapping,
                "mapping_long": [
                    {"reverse_bias_V": row["reverse_bias_V"], "support": support, "overlap": row[field]}
                    for row in mapping
                    for support, field in (("Element vertex", "vertex_overlap"), ("Cell", "cell_overlap"))
                ],
            },
        },
        "sources": sources,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    artifact = build(args.root.resolve())
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="utf-8"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
