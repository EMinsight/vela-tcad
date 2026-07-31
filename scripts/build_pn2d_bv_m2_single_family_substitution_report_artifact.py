#!/usr/bin/env python3
"""Build the portable report artifact for the M2 family substitutions."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


TITLE = "PN2D M2 single-family state substitution"


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def build(root: Path) -> dict[str, object]:
    result = json.loads((root / "result.json").read_text(encoding="utf-8-sig"))
    source = rows(root / "source_substitution.csv")
    newton = rows(root / "newton_first_update.csv")
    source_variants = {
        "sent_psi_only": "psi only",
        "sent_qfp_only": "QFP only",
        "sent_density_only": "n/p only",
    }
    recovery = [
        {
            "reverse_bias_V": abs(float(row["bias_V"])),
            "family": source_variants[row["variant"]],
            "source_ratio": float(row["source_to_sentaurus_ratio"]),
            "absolute_error_dex": float(row["abs_log10_error_dex"]),
            "error_reduction_dex": float(row["error_reduction_from_vela_dex"]),
            "recovery_fraction": float(
                row["fraction_of_all_sent_error_removal"]
            ),
        }
        for row in source
        if row["variant"] in source_variants
    ]
    projection_variants = {
        "feedback_density_only": "n/p feedback",
        "feedback_qfp_only": "QFP feedback",
    }
    projections = [
        {
            "reverse_bias_V": abs(float(row["bias_V"])),
            "intervention": projection_variants[row["variant"]],
            "qfp_target_projection": float(row["qfp_target_projection_fraction"]),
            "trial_residual_ratio": float(row["trial_combined_to_initial_ratio"]),
            "qfp_target_distance_ratio": float(
                row["qfp_trial_to_initial_distance_ratio"]
            ),
            "combined_target_distance_ratio": float(
                row["combined_trial_to_initial_distance_ratio"]
            ),
        }
        for row in newton
        if row["variant"] in projection_variants
    ]
    source_table = [
        {
            "reverse_bias_V": abs(float(row["bias_V"])),
            "variant": row["variant"],
            "source_ratio": float(row["source_to_sentaurus_ratio"]),
            "recovery_fraction": float(
                row["fraction_of_all_sent_error_removal"]
            ),
        }
        for row in source
    ]
    update_table = [
        {
            "reverse_bias_V": abs(float(row["bias_V"])),
            "variant": row["variant"],
            "trial_residual_ratio": float(row["trial_combined_to_initial_ratio"]),
            "combined_target_distance_ratio": float(
                row["combined_trial_to_initial_distance_ratio"]
            ),
            "qfp_target_projection": (
                float(row["qfp_target_projection_fraction"])
                if row["qfp_target_projection_fraction"]
                else None
            ),
        }
        for row in newton
        if row["variant"]
        in {"sent_qfp_only", "feedback_density_only", "feedback_qfp_only"}
    ]
    verdict = result["verdict"]
    sources = [
        {
            "id": "result",
            "label": "Machine-readable experiment result",
            "path": "build-release/pn2d-bv-m2-single-family-state-substitution-20260731/result.json",
        },
        {
            "id": "source",
            "label": "Fixed-source family substitutions",
            "path": "build-release/pn2d-bv-m2-single-family-state-substitution-20260731/source_substitution.csv",
            "query": {
                "engine": "duckdb",
                "language": "sql",
                "description": "Select the three one-family substitutions on the predeclared M2 bias lattice.",
                "sql": "SELECT abs(bias_V) AS reverse_bias_V, variant, source_to_sentaurus_ratio, abs_log10_error_dex, error_reduction_from_vela_dex, fraction_of_all_sent_error_removal FROM read_csv_auto('source_substitution.csv', header=true) WHERE variant IN ('sent_psi_only', 'sent_qfp_only', 'sent_density_only') ORDER BY reverse_bias_V, variant",
                "tables_used": ["source_substitution.csv"],
                "filters": ["variant is a single-family substitution", "bias in (-18, -19.5, -19.7, -20) V"],
                "metric_definitions": {
                    "source_ratio": "Vela integrated source divided by the Sentaurus integrated source.",
                    "recovery_fraction": "One-family error reduction divided by the full-Sentaurus-state error reduction.",
                },
            },
        },
        {
            "id": "newton",
            "label": "First coupled Newton updates",
            "path": "build-release/pn2d-bv-m2-single-family-state-substitution-20260731/newton_first_update.csv",
            "query": {
                "engine": "duckdb",
                "language": "sql",
                "description": "Select the density and QFP feedback interventions from the first-update audit.",
                "sql": "SELECT abs(bias_V) AS reverse_bias_V, variant, trial_combined_to_initial_ratio, qfp_trial_to_initial_distance_ratio, combined_trial_to_initial_distance_ratio, qfp_target_projection_fraction FROM read_csv_auto('newton_first_update.csv', header=true) WHERE variant IN ('feedback_density_only', 'feedback_qfp_only') ORDER BY reverse_bias_V, variant",
                "tables_used": ["newton_first_update.csv"],
                "filters": ["variant in density-only or QFP-only feedback", "first production update only"],
                "metric_definitions": {
                    "qfp_target_projection": "First QFP update projected onto the Vela-to-Sentaurus QFP target, normalized by squared target distance.",
                    "trial_residual_ratio": "Production trial residual norm divided by the intervention residual norm.",
                },
            },
        },
        {
            "id": "determinism",
            "label": "Independent-run output hashes",
            "path": "build-release/pn2d-bv-m2-single-family-state-substitution-20260731/determinism.csv",
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
                    "id": "recovery_chart",
                    "title": "Fraction of full-state source error removal",
                    "description": "Discrete family substitutions; negative values worsen agreement.",
                    "type": "bar",
                    "dataset": "recovery",
                    "sourceId": "source",
                    "encodings": {
                        "x": {"field": "reverse_bias_V", "type": "nominal", "title": "Reverse-bias magnitude (V)"},
                        "y": {"field": "recovery_fraction", "type": "quantitative", "title": "Recovered fraction"},
                        "color": {"field": "family", "type": "nominal"},
                    },
                },
                {
                    "id": "projection_chart",
                    "title": "First-update projection onto golden QFP direction",
                    "description": "Negative values mean the production update points away from Sentaurus QFP.",
                    "type": "bar",
                    "dataset": "projections",
                    "sourceId": "newton",
                    "encodings": {
                        "x": {"field": "reverse_bias_V", "type": "nominal", "title": "Reverse-bias magnitude (V)"},
                        "y": {"field": "qfp_target_projection", "type": "quantitative", "title": "Target projection fraction"},
                        "color": {"field": "intervention", "type": "nominal"},
                    },
                },
            ],
            "tables": [
                {
                    "id": "source_table",
                    "title": "Fixed-source substitutions",
                    "description": "Exact source ratios and normalized recovery fractions.",
                    "dataset": "source_table",
                    "sourceId": "source",
                    "columns": [
                        {"field": "reverse_bias_V", "label": "|V| (V)", "format": "number"},
                        {"field": "variant", "label": "Variant", "format": "text"},
                        {"field": "source_ratio", "label": "Source / Sentaurus", "format": "number"},
                        {"field": "recovery_fraction", "label": "Recovery fraction", "format": "number"},
                    ],
                    "defaultSort": {"field": "reverse_bias_V", "direction": "asc"},
                },
                {
                    "id": "update_table",
                    "title": "Selected first updates",
                    "description": "Residual and target-direction response for QFP and density interventions.",
                    "dataset": "update_table",
                    "sourceId": "newton",
                    "columns": [
                        {"field": "reverse_bias_V", "label": "|V| (V)", "format": "number"},
                        {"field": "variant", "label": "Variant", "format": "text"},
                        {"field": "trial_residual_ratio", "label": "Trial / initial residual", "format": "number"},
                        {"field": "combined_target_distance_ratio", "label": "Target distance ratio", "format": "number"},
                        {"field": "qfp_target_projection", "label": "QFP projection", "format": "number"},
                    ],
                    "defaultSort": {"field": "reverse_bias_V", "direction": "asc"},
                },
            ],
            "blocks": [
                {"id": "title", "type": "markdown", "body": f"# {TITLE}"},
                {
                    "id": "summary",
                    "type": "markdown",
                    "sourceId": "result",
                    "body": (
                        "## Answer-first finding\n\n"
                        "The typed outcome is **QFP dominant, with density feedback moving QFP away from Sentaurus**. At -20 V, QFP-only substitution recovers "
                        f"{100.0 * verdict['minus20_dominant_recovery_fraction']:.2f}% of the full-state source-error removal and wins three of four biases. The first density-feedback update has QFP target projection {verdict['minus20_density_feedback_qfp_projection_fraction']:.6f}."
                    ),
                },
                {
                    "id": "source_text",
                    "type": "markdown",
                    "sourceId": "source",
                    "body": (
                        "## Frozen source localization\n\n"
                        "QFP explanatory power grows into the knee region. Density-only substitution is source-neutral because this SG/Laux edge flux is reconstructed from psi, QFP, and intrinsic density; the configured Masetti plus field-dependent mobility depends on doping and field rather than imported n/p."
                    ),
                },
                {"id": "recovery", "type": "chart", "chartId": "recovery_chart"},
                {
                    "id": "newton_text",
                    "type": "markdown",
                    "sourceId": "newton",
                    "body": (
                        "## Dynamic response rejects the golden QFP direction\n\n"
                        "Starting from Sentaurus QFP with Vela psi increases the combined target distance at every bias. In the -19.5 to -20 V knee region, the first production update also increases the combined residual by 3.91-4.05 times. Both density and QFP feedback projections are negative."
                    ),
                },
                {"id": "projection", "type": "chart", "chartId": "projection_chart"},
                {"id": "source_values", "type": "table", "tableId": "source_table"},
                {"id": "update_values", "type": "table", "tableId": "update_table"},
                {
                    "id": "scope",
                    "type": "markdown",
                    "body": (
                        "## Scope and limitations\n\n"
                        "The experiment uses the shared M2 mesh and exactly -18, -19.5, -19.7, and -20 V. It is observation-only: SG/Laux, physics, defaults, continuation, and thresholds are unchanged. The evidence localizes the first discrepancy to carrier-QFP residual/Jacobian coupling, but does not yet prove a particular derivative or sign defect."
                    ),
                },
                {
                    "id": "next",
                    "type": "markdown",
                    "body": (
                        "## Next read-only check\n\n"
                        "Split electron and hole QFP, decompose transport, recombination, avalanche, and boundary residual terms, and finite-difference the carrier-QFP and Poisson-QFP Jacobian blocks on both baseline and mixed states before considering a production change."
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
                "recovery": recovery,
                "projections": projections,
                "source_table": source_table,
                "update_table": update_table,
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
