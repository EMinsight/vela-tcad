#!/usr/bin/env python3
"""Build the portable report artifact for the M2 frozen-state experiment."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


TITLE = "PN2D M2 Sentaurus-state SG/Laux frozen replay"


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def build(root: Path) -> dict[str, object]:
    result = json.loads((root / "result.json").read_text(encoding="utf-8-sig"))
    source_rows = rows(root / "source_comparison.csv")
    totals = sorted(
        (row for row in source_rows if row["carrier"] == "total"),
        key=lambda row: abs(float(row["bias_V"])),
    )
    ratio_rows = []
    error_rows = []
    table_rows = []
    for row in totals:
        voltage = abs(float(row["bias_V"]))
        ratio_rows.extend(
            [
                {"reverse_bias_V": voltage, "series": "Sentaurus golden", "source_ratio": 1.0},
                {
                    "reverse_bias_V": voltage,
                    "series": "Vela self-consistent",
                    "source_ratio": float(row["self_consistent_to_sentaurus_ratio"]),
                },
                {
                    "reverse_bias_V": voltage,
                    "series": "Vela frozen Sentaurus state",
                    "source_ratio": float(row["frozen_to_sentaurus_ratio"]),
                },
            ]
        )
        error_rows.extend(
            [
                {
                    "reverse_bias_V": voltage,
                    "series": "Self-consistent",
                    "error_dex": float(row["self_consistent_abs_log10_error_dex"]),
                },
                {
                    "reverse_bias_V": voltage,
                    "series": "Frozen Sentaurus state",
                    "error_dex": float(row["frozen_abs_log10_error_dex"]),
                },
            ]
        )
        table_rows.append(
            {
                "reverse_bias_V": voltage,
                "sentaurus_A_per_um": float(row["sentaurus_source_A_per_um"]),
                "self_consistent_A_per_um": float(
                    row["vela_self_consistent_source_A_per_um"]
                ),
                "frozen_A_per_um": float(
                    row["vela_frozen_sentaurus_state_source_A_per_um"]
                ),
                "self_ratio": float(row["self_consistent_to_sentaurus_ratio"]),
                "frozen_ratio": float(row["frozen_to_sentaurus_ratio"]),
                "error_reduction_dex": float(row["frozen_error_reduction_dex"]),
            }
        )

    verdict = result["verdict"]
    sources = [
        {
            "id": "result",
            "label": "Frozen-state machine-readable result",
            "path": "build-release/pn2d-bv-m2-sentaurus-frozen-sg-laux-20260731/result.json",
        },
        {
            "id": "source_comparison",
            "label": "Carrier-resolved source comparison",
            "path": "build-release/pn2d-bv-m2-sentaurus-frozen-sg-laux-20260731/source_comparison.csv",
        },
        {
            "id": "determinism",
            "label": "Independent-run output hashes",
            "path": "build-release/pn2d-bv-m2-sentaurus-frozen-sg-laux-20260731/determinism.csv",
        },
        {
            "id": "sentaurus_manifest",
            "label": "Sentaurus M2 process manifest",
            "path": "build-release/pn2d-task10-balanced-m2-sentaurus-process-v2-20260731/manifest.json",
        },
        {
            "id": "vela_manifest",
            "label": "Vela M2 self-consistent process manifest",
            "path": "build-release/pn2d-bv-template-default-prospective-v2-default-20260731/M2/run-a/manifest.json",
        },
        {
            "id": "ratio_query",
            "label": "Total-source ratio transformation",
            "path": "build-release/pn2d-bv-m2-sentaurus-frozen-sg-laux-20260731/source_comparison.csv",
            "query": {
                "engine": "duckdb",
                "language": "sql",
                "description": "Compare self-consistent and frozen-state Vela total source against Sentaurus at four M2 knee biases.",
                "sql": "SELECT abs(bias_V) AS reverse_bias_V, carrier, sentaurus_source_A_per_um, vela_self_consistent_source_A_per_um, vela_frozen_sentaurus_state_source_A_per_um, self_consistent_to_sentaurus_ratio, frozen_to_sentaurus_ratio, self_consistent_abs_log10_error_dex, frozen_abs_log10_error_dex, frozen_error_reduction_dex FROM read_csv_auto('source_comparison.csv', header=true) WHERE carrier = 'total' ORDER BY reverse_bias_V",
                "tables_used": ["source_comparison.csv"],
                "filters": ["carrier = total", "bias in (-18, -19.5, -19.7, -20) V"],
                "metric_definitions": {
                    "source_ratio": "Vela integrated impact-ionization source divided by the Sentaurus source on the same M2 mesh and bias.",
                    "error_dex": "Absolute log10 of the source ratio.",
                },
            },
        },
    ]
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
                    "id": "ratio_chart",
                    "title": "Integrated source ratio on four M2 knee biases",
                    "description": "Ratio to Sentaurus; 1.0 is exact agreement on the shared mesh.",
                    "type": "bar",
                    "dataset": "source_ratios",
                    "sourceId": "ratio_query",
                    "encodings": {
                        "x": {"field": "reverse_bias_V", "type": "nominal", "title": "Reverse-bias magnitude (V)"},
                        "y": {"field": "source_ratio", "type": "quantitative", "title": "Source / Sentaurus"},
                        "color": {"field": "series", "type": "nominal"},
                    },
                },
                {
                    "id": "error_chart",
                    "title": "Absolute integrated-source error",
                    "description": "Four discrete knee-region checks; lower is better.",
                    "type": "bar",
                    "dataset": "source_errors",
                    "sourceId": "ratio_query",
                    "encodings": {
                        "x": {"field": "reverse_bias_V", "type": "nominal", "title": "Reverse-bias magnitude (V)"},
                        "y": {"field": "error_dex", "type": "quantitative", "title": "Absolute error (dex)"},
                        "color": {"field": "series", "type": "nominal"},
                    },
                },
            ],
            "tables": [
                {
                    "id": "source_table",
                    "title": "Total integrated source audit",
                    "description": "Exact values and ratios for the four predeclared frozen-state biases.",
                    "dataset": "source_table",
                    "sourceId": "ratio_query",
                    "columns": [
                        {"field": "reverse_bias_V", "label": "|V| (V)", "format": "number"},
                        {"field": "sentaurus_A_per_um", "label": "Sentaurus (A/um)", "format": "number"},
                        {"field": "self_consistent_A_per_um", "label": "Vela self-consistent (A/um)", "format": "number"},
                        {"field": "frozen_A_per_um", "label": "Vela frozen (A/um)", "format": "number"},
                        {"field": "self_ratio", "label": "Self / Sentaurus", "format": "number"},
                        {"field": "frozen_ratio", "label": "Frozen / Sentaurus", "format": "number"},
                        {"field": "error_reduction_dex", "label": "Error reduction (dex)", "format": "number"},
                    ],
                    "defaultSort": {"field": "reverse_bias_V", "direction": "asc"},
                }
            ],
            "blocks": [
                {"id": "title", "type": "markdown", "body": f"# {TITLE}"},
                {
                    "id": "summary",
                    "type": "markdown",
                    "sourceId": "result",
                    "body": (
                        "## Technical summary\n\n"
                        "The discriminating experiment classifies the M2 knee discrepancy as **state-feedback dominant**. With Sentaurus potential, electron/hole quasi-Fermi potentials, and carrier densities frozen, Vela SG/Laux reproduces the Sentaurus total impact-ionization source to a mean absolute error of "
                        f"{verdict['mean_frozen_abs_log10_error_dex']:.6f} dex. The corresponding self-consistent Vela error is {verdict['mean_self_consistent_abs_log10_error_dex']:.6f} dex, so freezing the golden state removes {verdict['mean_error_reduction_dex']:.6f} dex on average. The result does not authorize a solver or default-value change; it identifies where the next correction should be sought."
                    ),
                },
                {
                    "id": "finding_ratio",
                    "type": "markdown",
                    "sourceId": "source_comparison",
                    "body": (
                        "## The same SG/Laux operator closes on the Sentaurus state\n\n"
                        "At 18, 19.5, 19.7, and 20 V reverse-bias magnitude, frozen-state Vela/Sentaurus total-source ratios are 1.00248, 1.00238, 1.00240, and 1.00237. The carrier-resolved electron and hole ratios are also confined to 1.00224-1.00264. This makes a large frozen SG/Laux alpha-current-source mismatch incompatible with the observed result; the production operator chain is accurate to about 0.24% on the golden state."
                    ),
                },
                {"id": "ratio", "type": "chart", "chartId": "ratio_chart"},
                {
                    "id": "finding_error",
                    "type": "markdown",
                    "sourceId": "source_comparison",
                    "body": (
                        "## Self-consistent state formation creates the growing deficit\n\n"
                        "The self-consistent total-source ratio falls from 0.936 at -18 V to 0.825 at -20 V, while the frozen-state ratio remains essentially flat near 1.0024. The self-consistent source error grows monotonically from 0.02866 to 0.08354 dex; the frozen error stays near 0.00103 dex. The bias-dependent discrepancy therefore enters before or during the coupled state-feedback loop, not during read-only evaluation of SG/Laux on a fixed state."
                    ),
                },
                {"id": "error", "type": "chart", "chartId": "error_chart"},
                {
                    "id": "table_text",
                    "type": "markdown",
                    "body": (
                        "## Exact source values\n\n"
                        "The audit table keeps the A/um values needed to reproduce the ratios. Error reduction is the self-consistent absolute log error minus the frozen-state absolute log error; positive values show how much mismatch disappears when the golden state is substituted."
                    ),
                },
                {"id": "table", "type": "table", "tableId": "source_table"},
                {
                    "id": "scope",
                    "type": "markdown",
                    "body": (
                        "## Scope, state contract, and metric definitions\n\n"
                        "The experiment uses the shared M2 mesh at exactly -18, -19.5, -19.7, and -20 V on the avalanche-on branch. Sentaurus is the golden reference. Each state contains 115 common physical nodes and five imported fields: electrostatic potential, electron and hole quasi-Fermi potential, and electron and hole density. Seven duplicated Sentaurus contact-support vertices are excluded because their coordinates duplicate existing M2 contact nodes. Source error is |log10(Vela/Sentaurus)|."
                    ),
                },
                {
                    "id": "method",
                    "type": "markdown",
                    "sourceId": "result",
                    "body": (
                        "## Observation-only method and robustness checks\n\n"
                        "Vela evaluates van Overstraeten impact ionization with quasi-Fermi-gradient drive, current-density generation, the complete element-edge SG/GSS/Laux current vector, and element-vertex box-measure mapping. Coupling is forced to postprocess_only. The state round trip is bit-exact numerically, every process record has solver_coupled=0, residual-feedback arrays are empty, and qG closure is within 2.83e-15 relative. Two independent executions produce byte-identical node, edge, triangle, element, and process files at all four biases."
                    ),
                },
                {
                    "id": "limitations",
                    "type": "markdown",
                    "body": (
                        "## What this result does and does not establish\n\n"
                        "The test proves that the combined frozen SG/Laux mobility-current-alpha-source chain closes against Sentaurus on these four M2 states. It does not by itself identify which coupled state equation first diverges, and small compensating local errors could still exist beneath the 0.24% integrated closure. Spatial and carrier-resolved evidence from the prior read-only decomposition remains necessary when selecting a minimal correction."
                    ),
                },
                {
                    "id": "next",
                    "type": "markdown",
                    "body": (
                        "## Recommended next step\n\n"
                        "Keep SG/Laux and all acceptance thresholds unchanged. On the same four biases, perform one-at-a-time mixed-state substitutions for psi, QFP, and n/p inside the frozen audit, then evaluate the first coupled Newton update and carrier-row residual. This should distinguish whether the initiating error is density feedback, Poisson-QFP cross-coupling, or the update path that transports an initially small state mismatch into the 17.5% source deficit at -20 V."
                    ),
                },
                {
                    "id": "questions",
                    "type": "markdown",
                    "body": (
                        "## Further questions\n\n"
                        "Which single state family recovers most of the 0.0825 dex improvement at -20 V? Does the first coupled update move that field away from the Sentaurus state in the same direction as the final self-consistent deficit? Are local errors still compensated after matching the integrated source?"
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
                "source_ratios": ratio_rows,
                "source_errors": error_rows,
                "source_table": table_rows,
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
