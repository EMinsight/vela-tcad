#!/usr/bin/env python3
"""Build the canonical Data Analytics artifact for the PN2D basis comparison."""

from __future__ import annotations

import argparse
import csv
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--comparison", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    comparison = read_csv(args.comparison)
    summary = json.loads(args.summary.read_text(encoding="utf-8"))
    by_bias: dict[int, dict[str, Any]] = {}
    for row in comparison:
        bias = int(round(abs(float(row["bias_V"]))))
        target = by_bias.setdefault(
            bias,
            {
                "reverse_bias_V": bias,
                "sentaurus_A_per_um": float(row["sentaurus_current_A_per_um"]),
            },
        )
        basis = row["basis"]
        target[f"{basis}_A_per_um"] = float(row["vela_current_A_per_um"])
        target[f"{basis}_over_sentaurus"] = float(row["vela_over_sentaurus"])
        target[f"{basis}_closure_ratio"] = float(
            row["global_closure_max_ratio"]
        )
    curve_rows = [by_bias[bias] for bias in sorted(by_bias)]
    chart_rows: list[dict[str, Any]] = []
    for row in curve_rows:
        bias = int(row["reverse_bias_V"])
        chart_rows.append(
            {
                "reverse_bias_V": bias,
                "series": "Sentaurus avalanche-off",
                "current_A_per_um": row["sentaurus_A_per_um"],
                "vela_over_sentaurus": 1.0,
                "closure_ratio": None,
                "is_reference": True,
            }
        )
        for basis in (
            "net_doping",
            "total_impurity",
            "cell_reconstructed_total_impurity",
        ):
            chart_rows.append(
                {
                    "reverse_bias_V": bias,
                    "series": basis,
                    "current_A_per_um": row[f"{basis}_A_per_um"],
                    "vela_over_sentaurus": row[
                        f"{basis}_over_sentaurus"
                    ],
                    "closure_ratio": row[f"{basis}_closure_ratio"],
                    "is_reference": False,
                }
            )

    summary_rows: list[dict[str, Any]] = []
    for basis, values in summary["bases"].items():
        pairwise = summary["pairwise"].get(f"{basis}_vs_net_doping", {})
        summary_rows.append(
            {
                "basis": basis,
                "converged_points": values["converged_points"],
                "log10_rmse": values["log10_ratio_rmse"],
                "median_abs_relative_error": values[
                    "median_absolute_relative_error"
                ],
                "max_closure_ratio": values["max_global_closure_ratio"],
                "ratio_minus_1V": values["ratio_at_minus_1V"],
                "ratio_minus_10V": values["ratio_at_minus_10V"],
                "ratio_minus_20V": values["ratio_at_minus_20V"],
                "max_relative_delta_vs_net": pairwise.get(
                    "max_absolute_relative_difference", 0.0
                ),
            }
        )

    generated_at = datetime.now(timezone.utc).isoformat()
    source = {
        "id": "basis_comparison",
        "label": "Vela/Sentaurus avalanche-off comparison",
        "path": "comparison.csv",
        "query": {
            "engine": "DuckDB",
            "language": "sql",
            "sql": (
                "SELECT abs(bias_V) AS reverse_bias_V, basis AS series, "
                "vela_current_A_per_um AS current_A_per_um, "
                "vela_over_sentaurus, global_closure_max_ratio AS "
                "closure_ratio, false AS is_reference "
                "FROM read_csv_auto('comparison.csv') "
                "UNION ALL "
                "SELECT DISTINCT abs(bias_V), 'Sentaurus avalanche-off', "
                "sentaurus_current_A_per_um, 1.0, NULL, true "
                "FROM read_csv_auto('comparison.csv')"
            ),
            "description": (
                "Three Vela runs differing only in mobility "
                "doping_concentration_basis, joined to the refreshed "
                "Sentaurus avalanche-off integer-bias reference."
            ),
            "tables_used": [
                "comparison.csv",
                "summary.json",
            ],
            "filters": [
                "PN2D coarse7x3",
                "impact ionization disabled",
                "reverse biases -1 V through -20 V",
                "SRH enabled",
            ],
            "metric_definitions": [
                "log10 RMSE = RMS(log10(|I_Vela|/|I_Sentaurus|)) over 20 nonzero biases",
                "median absolute relative error = median(abs(I_Vela/I_Sentaurus - 1))",
                "closure ratio = max(electron global closure ratio, hole global closure ratio)",
            ],
        },
    }
    title = "PN2D avalanche-off mobility doping-basis comparison"
    artifact = {
        "surface": "report",
        "manifest": {
            "version": 1,
            "surface": "report",
            "title": title,
            "description": (
                "Technical comparison of net, total, and cell-reconstructed "
                "mobility doping bases after the 2-D source-scale repair."
            ),
            "generatedAt": generated_at,
            "sources": [source],
            "charts": [
                {
                    "id": "reverse_current",
                    "title": "Avalanche-off reverse current",
                    "subtitle": (
                        "PN2D coarse7x3, SRH enabled, impact ionization disabled; "
                        "20 nonzero integer bias points"
                    ),
                    "intent": "comparison",
                    "question": (
                        "Does mobility doping basis explain the remaining "
                        "difference from Sentaurus?"
                    ),
                    "rationale": (
                        "A multi-series line chart shows voltage dependence "
                        "and whether any basis separates from the others."
                    ),
                    "comparisonContext": {
                        "baseline": "Sentaurus avalanche-off",
                        "grain": "integer reverse-bias point",
                        "unit": "A/um",
                    },
                    "type": "line",
                    "dataset": "curve",
                    "sourceId": "basis_comparison",
                    "encodings": {
                        "x": {
                            "field": "reverse_bias_V",
                            "type": "quantitative",
                            "label": "Reverse bias |V|",
                            "unit": "V",
                        },
                        "y": {
                            "field": "current_A_per_um",
                            "type": "quantitative",
                            "label": "|Anode current|",
                            "unit": "A/um",
                        },
                        "color": {
                            "field": "series",
                            "type": "nominal",
                            "label": "Simulation",
                        },
                        "tooltip": [
                            {
                                "field": "vela_over_sentaurus",
                                "type": "quantitative",
                                "label": "Vela/Sentaurus",
                            },
                            {
                                "field": "closure_ratio",
                                "type": "quantitative",
                                "label": "Closure ratio",
                            },
                        ],
                    },
                    "valueFormat": "compact",
                    "unit": "A/um",
                    "layout": "full",
                    "compatibleTypes": ["line"],
                    "maxRows": 20,
                }
            ],
            "tables": [
                {
                    "id": "basis_metrics",
                    "title": "Basis comparison metrics",
                    "subtitle": (
                        "Sentaurus error, numerical closure, and sensitivity "
                        "relative to net_doping"
                    ),
                    "dataset": "summary",
                    "sourceId": "basis_comparison",
                    "defaultSort": {
                        "field": "log10_rmse",
                        "direction": "asc",
                    },
                    "density": "dense",
                    "layout": "full",
                    "columns": [
                        {"field": "basis", "label": "Basis", "type": "text"},
                        {
                            "field": "converged_points",
                            "label": "Converged",
                            "format": "number",
                        },
                        {
                            "field": "log10_rmse",
                            "label": "log10 RMSE",
                            "format": "number",
                        },
                        {
                            "field": "median_abs_relative_error",
                            "label": "Median abs. error",
                            "format": "percent",
                        },
                        {
                            "field": "max_closure_ratio",
                            "label": "Max closure ratio",
                            "format": "percent",
                        },
                        {
                            "field": "ratio_minus_1V",
                            "label": "-1 V ratio",
                            "format": "number",
                        },
                        {
                            "field": "ratio_minus_10V",
                            "label": "-10 V ratio",
                            "format": "number",
                        },
                        {
                            "field": "ratio_minus_20V",
                            "label": "-20 V ratio",
                            "format": "number",
                        },
                        {
                            "field": "max_relative_delta_vs_net",
                            "label": "Max delta vs net",
                            "format": "percent",
                        },
                    ],
                }
            ],
            "blocks": [
                {"id": "title", "type": "markdown", "body": f"# {title}"},
                {
                    "id": "technical_summary",
                    "type": "markdown",
                    "sourceId": "basis_comparison",
                    "body": (
                        "## Technical summary\n\n"
                        "Changing `doping_concentration_basis` does not "
                        "materially change the repaired avalanche-off BV "
                        "curve. All three candidates converge at 21/21 "
                        "points. Their maximum pairwise current difference is "
                        "0.0113%, while the median absolute error against "
                        "Sentaurus remains about 36.18%. The remaining "
                        "voltage-shape mismatch is therefore not driven by "
                        "the mobility doping basis."
                    ),
                },
                {
                    "id": "findings",
                    "type": "markdown",
                    "sourceId": "basis_comparison",
                    "body": (
                        "## The three Vela curves are numerically indistinguishable\n\n"
                        "`cell_reconstructed_total_impurity` has the smallest "
                        "log10-ratio RMSE (0.221262), but the advantage over "
                        "`net_doping` (0.221264) is too small to be physically "
                        "meaningful. At -10 V its current differs from "
                        "`net_doping` by only 0.0113%; at -20 V the difference "
                        "is 0.000016%."
                    ),
                },
                {
                    "id": "chart",
                    "type": "chart",
                    "chartId": "reverse_current",
                    "layout": "full",
                },
                {
                    "id": "scope",
                    "type": "markdown",
                    "body": (
                        "## Scope and metric definitions\n\n"
                        "The population is the coarse7x3 PN2D device at 300 K. "
                        "The comparison uses integer biases from -1 V through "
                        "-20 V, SRH plus Old Slotboom BGN, Masetti-field "
                        "mobility, impact ionization disabled, and the repaired "
                        "2-D continuity-source scaling. Error metrics exclude "
                        "the zero-bias numerical floor."
                    ),
                },
                {
                    "id": "method",
                    "type": "markdown",
                    "body": (
                        "## Controlled comparison method\n\n"
                        "Each run starts from the same validated JSON. The only "
                        "changed physics field is "
                        "`solver.mobility.doping_concentration_basis`: "
                        "`net_doping`, `total_impurity`, or "
                        "`cell_reconstructed_total_impurity`. Convergence, "
                        "global carrier closure, terminal-current extraction, "
                        "bias points, and Sentaurus reference data are held "
                        "fixed."
                    ),
                },
                {
                    "id": "robustness",
                    "type": "markdown",
                    "sourceId": "basis_comparison",
                    "body": (
                        "## Numerical closure does not explain the overlap\n\n"
                        "The worst global carrier-closure ratio is 0.0132%, "
                        "well below the enforced 1% threshold. The three "
                        "curves' overlap is therefore a resolved physical "
                        "sensitivity result, not failed convergence or lost "
                        "source current."
                    ),
                },
                {
                    "id": "table",
                    "type": "table",
                    "tableId": "basis_metrics",
                    "layout": "full",
                },
                {
                    "id": "next_steps",
                    "type": "markdown",
                    "body": (
                        "## Recommended next steps\n\n"
                        "Keep `net_doping` as the BV default because no "
                        "alternative improves avalanche-off agreement. Retain "
                        "`cell_reconstructed_total_impurity` for forward IV, "
                        "where its benefit was previously measurable. For the "
                        "remaining BV shape difference, prioritize depletion-"
                        "region SRH spatial support, effective lifetime/trap "
                        "parameters, and mesh resolution."
                    ),
                },
                {
                    "id": "questions",
                    "type": "markdown",
                    "body": (
                        "## Further questions\n\n"
                        "The next discriminating experiments are: whether a "
                        "refined junction mesh restores Sentaurus-like SRH "
                        "growth from -1 V to -15 V; whether Sentaurus uses "
                        "field-, doping-, or temperature-dependent lifetime "
                        "terms absent from Vela; and how the corrected source "
                        "scale changes avalanche-on gain."
                    ),
                },
            ],
        },
        "snapshot": {
            "version": 1,
            "generatedAt": generated_at,
            "status": "ready",
            "datasets": {
                "curve": chart_rows,
                "summary": summary_rows,
            },
            "accessIssues": [],
        },
        "sources": [source],
    }
    args.output.write_text(
        json.dumps(artifact, indent=2) + "\n", encoding="utf-8"
    )
    print(args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
