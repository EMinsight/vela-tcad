#!/usr/bin/env python3
"""Build the canonical technical-report artifact for the M2 carrier-QFP audit."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def write_json(path: Path, payload: object) -> None:
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    output = args.output or args.root / "report_artifact.json"

    result = json.loads((args.root / "result.json").read_text(encoding="utf-8"))
    verdict = result["verdict"]
    source_rows = read_rows(args.root / "source_carrier_substitution.csv")
    term_rows = read_rows(args.root / "carrier_term_decomposition.csv")
    update_rows = read_rows(args.root / "first_qfp_updates.csv")
    jacobian_rows = read_rows(args.root / "jacobian_fd_blocks.csv")
    sensitivity_rows = read_rows(args.root / "jacobian_fd_step_sensitivity.csv")

    recovery = [
        {
            "reverse_bias_V": abs(float(row["bias_V"])),
            "carrier_qfp": "electron QFP" if row["variant"] == "sent_phin_only" else "hole QFP",
            "source_to_sentaurus_ratio": float(row["source_to_sentaurus_ratio"]),
            "qfp_error_recovery_fraction": float(row["fraction_of_qfp_error_removal"]),
        }
        for row in source_rows
        if row["variant"] in {"sent_phin_only", "sent_phip_only"}
    ]
    terms = [
        {
            "carrier": row["carrier"],
            "term": "transport / SG flux" if row["term"] == "flux" else "avalanche source",
            "residual_delta_share": float(row["term_delta_share"]),
            "term_delta_l2": float(row["term_delta_l2"]),
        }
        for row in term_rows
        if float(row["bias_V"]) == -20.0
        and row["variant"] == "sent_qfp_only"
        and row["scope"] == "interior"
        and row["term"] in {"flux", "impact"}
    ]
    sensitivity_grouped: dict[tuple[float, str, float], list[float]] = {}
    for row in sensitivity_rows:
        key = (
            float(row["bias_V"]), row["state_variant"],
            float(row["finite_difference_step"]),
        )
        sensitivity_grouped.setdefault(key, []).extend(
            float(row[f"diff_{field}_norm"])
            for field in ("electron_phin", "electron_phip", "hole_phin", "hole_phip")
        )
    sensitivity = [
        {
            "finite_difference_step": step,
            "case": f"{abs(bias):g} V / {'mixed QFP' if state == 'sent_qfp_only' else 'baseline'}",
            "maximum_absolute_difference": max(values),
            "absolute_floor_threshold": 1.0e-13,
        }
        for (bias, state, step), values in sorted(sensitivity_grouped.items())
    ]

    relevant = {
        "poisson": ("rel_phin_column_diff", "rel_phip_column_diff"),
        "transport": (
            "rel_electron_phin_diff", "rel_electron_phip_diff",
            "rel_hole_phin_diff", "rel_hole_phip_diff",
        ),
        "srh_auger": (
            "rel_electron_phin_diff", "rel_electron_phip_diff",
            "rel_hole_phin_diff", "rel_hole_phip_diff",
        ),
        "sg_avalanche": (
            "rel_electron_phin_diff", "rel_electron_phip_diff",
            "rel_hole_phin_diff", "rel_hole_phip_diff",
        ),
        "dirichlet_or_gauge": (
            "rel_electron_phin_diff", "rel_electron_phip_diff",
            "rel_hole_phin_diff", "rel_hole_phip_diff",
        ),
    }
    jacobian_summary = []
    for block, fields in relevant.items():
        candidates = [
            (float(row[field]), row, field)
            for row in jacobian_rows if row["block"] == block for field in fields
        ]
        value, row, field = max(candidates)
        jacobian_summary.append({
            "block": block,
            "maximum_relative_difference": value,
            "bias_V": float(row["bias_V"]),
            "state": row["state_variant"],
            "subblock": field.removeprefix("rel_").removesuffix("_diff"),
            "formal_gate": "pass" if value <= float(result["finite_difference_threshold"]) else "fail",
        })
    updates = [
        {
            "carrier": row["carrier"],
            "initial_residual": float(row["initial_residual"]),
            "trial_residual": float(row["trial_residual"]),
            "trial_to_initial_residual": float(row["trial_residual"]) / float(row["initial_residual"]),
            "projection_toward_sentaurus": float(row["update_projection_on_vela_to_sentaurus_target"]),
            "trial_target_distance_ratio": float(row["trial_to_vela_target_distance_ratio"]),
        }
        for row in update_rows
        if float(row["bias_V"]) == -20.0 and row["variant"] == "sent_qfp_only"
    ]

    base = "build-release/pn2d-bv-m2-qfp-carrier-jacobian-verification-20260731"
    sources: list[dict[str, Any]] = [
        {"id": "result", "label": "Machine-readable verdict", "path": f"{base}/result.json"},
        {
            "id": "source", "label": "Carrier-resolved frozen source substitutions",
            "path": f"{base}/source_carrier_substitution.csv",
            "query": {
                "engine": "duckdb", "language": "sql",
                "description": "Select electron-QFP and hole-QFP substitutions on the frozen four-bias lattice.",
                "sql": "SELECT abs(bias_V) AS reverse_bias_V, CASE WHEN variant='sent_phin_only' THEN 'electron QFP' ELSE 'hole QFP' END AS carrier_qfp, source_to_sentaurus_ratio, fraction_of_qfp_error_removal AS qfp_error_recovery_fraction FROM read_csv_auto('source_carrier_substitution.csv', header=true) WHERE variant IN ('sent_phin_only','sent_phip_only') ORDER BY reverse_bias_V, carrier_qfp",
                "tables_used": ["source_carrier_substitution.csv"],
                "filters": ["bias in (-18, -19.5, -19.7, -20) V", "one carrier-QFP family replaced"],
                "metric_definitions": {"qfp_error_recovery_fraction": "One-carrier QFP log-source-error reduction divided by the joint-QFP reduction."},
            },
        },
        {
            "id": "terms", "label": "Carrier continuity residual decomposition",
            "path": f"{base}/carrier_term_decomposition.csv",
            "query": {
                "engine": "duckdb", "language": "sql",
                "description": "Select interior transport and avalanche residual shares for joint-QFP substitution at -20 V.",
                "sql": "SELECT carrier, CASE WHEN term='flux' THEN 'transport / SG flux' ELSE 'avalanche source' END AS term, term_delta_share AS residual_delta_share, term_delta_l2 FROM read_csv_auto('carrier_term_decomposition.csv', header=true) WHERE bias_V=-20 AND variant='sent_qfp_only' AND scope='interior' AND term IN ('flux','impact') ORDER BY carrier, term",
                "tables_used": ["carrier_term_decomposition.csv"],
                "filters": ["bias=-20 V", "variant=joint Sentaurus QFP", "scope=interior nodes"],
                "metric_definitions": {"residual_delta_share": "Term-delta L2 divided by the sum of all five carrier term-delta L2 norms."},
            },
        },
        {
            "id": "updates", "label": "First coupled Newton updates",
            "path": f"{base}/first_qfp_updates.csv",
            "query": {
                "engine": "duckdb", "language": "sql",
                "description": "Select electron and hole first updates from the joint Sentaurus-QFP state at -20 V.",
                "sql": "SELECT carrier, initial_residual, trial_residual, trial_residual/initial_residual AS trial_to_initial_residual, update_projection_on_vela_to_sentaurus_target AS projection_toward_sentaurus, trial_to_vela_target_distance_ratio AS trial_target_distance_ratio FROM read_csv_auto('first_qfp_updates.csv', header=true) WHERE bias_V=-20 AND variant='sent_qfp_only' ORDER BY carrier",
                "tables_used": ["first_qfp_updates.csv"],
                "filters": ["bias=-20 V", "variant=joint Sentaurus QFP", "first production Newton update"],
                "metric_definitions": {"projection_toward_sentaurus": "Carrier-QFP update projected onto the Vela-to-Sentaurus target and normalized by target norm squared."},
            },
        },
        {
            "id": "jacobian", "label": "Analytic versus finite-difference Jacobian blocks",
            "path": f"{base}/jacobian_fd_blocks.csv",
            "query": {
                "engine": "duckdb", "language": "sql",
                "description": "Take the maximum relevant carrier-QFP or Poisson-QFP relative difference for each physical block.",
                "sql": "WITH long AS (UNPIVOT read_csv_auto('jacobian_fd_blocks.csv', header=true) ON rel_phin_column_diff, rel_phip_column_diff, rel_electron_phin_diff, rel_electron_phip_diff, rel_hole_phin_diff, rel_hole_phip_diff INTO NAME subblock VALUE relative_difference), ranked AS (SELECT *, row_number() OVER (PARTITION BY block ORDER BY relative_difference DESC) AS rank FROM long) SELECT block, relative_difference AS maximum_relative_difference, bias_V, state_variant AS state, subblock, CASE WHEN relative_difference<=5e-5 THEN 'pass' ELSE 'fail' END AS formal_gate FROM ranked WHERE rank=1 ORDER BY maximum_relative_difference DESC",
                "tables_used": ["jacobian_fd_blocks.csv"],
                "filters": ["baseline and joint-QFP states", "bias in (-18, -19.5, -19.7, -20) V", "predeclared relative gate=5e-5"],
                "metric_definitions": {"maximum_relative_difference": "Frobenius norm of analytic-minus-finite-difference subblock divided by the larger analytic or finite-difference norm."},
            },
        },
        {
            "id": "sensitivity", "label": "SRH finite-difference step sensitivity",
            "path": f"{base}/jacobian_fd_step_sensitivity.csv",
            "query": {
                "engine": "duckdb", "language": "sql",
                "description": "Aggregate the maximum SRH/Auger absolute carrier-QFP subblock difference by bias, state, and perturbation step.",
                "sql": "SELECT finite_difference_step, concat(abs(bias_V),' V / ',CASE WHEN state_variant='sent_qfp_only' THEN 'mixed QFP' ELSE 'baseline' END) AS case, greatest(diff_electron_phin_norm,diff_electron_phip_norm,diff_hole_phin_norm,diff_hole_phip_norm) AS maximum_absolute_difference, 1e-13 AS absolute_floor_threshold FROM read_csv_auto('jacobian_fd_step_sensitivity.csv', header=true) ORDER BY case, finite_difference_step",
                "tables_used": ["jacobian_fd_step_sensitivity.csv"],
                "filters": ["block=SRH/Auger", "bias in (-19.5, -20) V", "seven predeclared perturbation steps"],
                "metric_definitions": {"maximum_absolute_difference": "Maximum Frobenius norm of analytic-minus-finite-difference across four carrier-QFP subblocks."},
            },
        },
        {"id": "determinism", "label": "Independent-run hashes", "path": f"{base}/determinism.csv"},
    ]
    charts = [
        {
            "id": "recovery_chart", "title": "Integrated-source recovery by QFP carrier",
            "description": "Four predeclared M2 reverse biases; negative recovery worsens agreement.",
            "type": "bar", "dataset": "recovery", "sourceId": "source",
            "encodings": {
                "x": {"field": "reverse_bias_V", "type": "nominal", "title": "Reverse-bias magnitude (V)"},
                "y": {"field": "qfp_error_recovery_fraction", "type": "quantitative", "title": "Fraction of joint-QFP error recovery"},
                "color": {"field": "carrier_qfp", "type": "nominal"},
            },
        },
        {
            "id": "term_chart", "title": "Carrier residual change by physical term",
            "description": "Interior L2 share at -20 V after joint Sentaurus-QFP substitution.",
            "type": "bar", "dataset": "terms", "sourceId": "terms",
            "encodings": {
                "x": {"field": "carrier", "type": "nominal", "title": "Continuity equation"},
                "y": {"field": "residual_delta_share", "type": "quantitative", "title": "Share of term-delta L2 sum"},
                "color": {"field": "term", "type": "nominal"},
            },
        },
        {
            "id": "sensitivity_chart", "title": "SRH/Auger Jacobian finite-difference sensitivity",
            "description": "Maximum absolute subblock difference at -19.5 and -20 V across seven perturbation steps.",
            "type": "line", "dataset": "sensitivity", "sourceId": "sensitivity",
            "encodings": {
                "x": {"field": "finite_difference_step", "type": "quantitative", "title": "Finite-difference step (V)"},
                "y": {"field": "maximum_absolute_difference", "type": "quantitative", "title": "Maximum absolute matrix-norm difference"},
                "color": {"field": "case", "type": "nominal"},
            },
        },
    ]
    tables = [
        {
            "id": "jacobian_table", "title": "Maximum finite-difference error by Jacobian block",
            "description": "Baseline and mixed-QFP states over -18, -19.5, -19.7, and -20 V.",
            "dataset": "jacobian_summary", "sourceId": "jacobian",
            "columns": [
                {"field": "block", "label": "Block", "format": "text"},
                {"field": "maximum_relative_difference", "label": "Maximum relative difference", "format": "number"},
                {"field": "bias_V", "label": "Bias (V)", "format": "number"},
                {"field": "state", "label": "State", "format": "text"},
                {"field": "subblock", "label": "Subblock", "format": "text"},
                {"field": "formal_gate", "label": "5e-5 gate", "format": "text"},
            ],
            "defaultSort": {"field": "maximum_relative_difference", "direction": "desc"},
        },
        {
            "id": "updates_table", "title": "First coupled update from the Sentaurus-QFP state at -20 V",
            "description": "Negative projection points away from the Sentaurus carrier QFP.",
            "dataset": "updates", "sourceId": "updates",
            "columns": [
                {"field": "carrier", "label": "Carrier", "format": "text"},
                {"field": "initial_residual", "label": "Initial residual", "format": "number"},
                {"field": "trial_to_initial_residual", "label": "Trial / initial residual", "format": "number"},
                {"field": "projection_toward_sentaurus", "label": "Projection toward Sentaurus", "format": "number"},
                {"field": "trial_target_distance_ratio", "label": "Target-distance ratio", "format": "number"},
            ],
            "defaultSort": {"field": "carrier", "direction": "asc"},
        },
    ]
    blocks = [
        {"id": "title", "type": "markdown", "body": "# M2 carrier-QFP residual and Jacobian verification"},
        {"id": "summary", "type": "markdown", "sourceId": "result", "body": (
            "## Technical summary\n\n"
            "The first self-consistent mismatch is localized to **carrier transport response to QFP**, not to the SG/Laux avalanche derivative or Poisson-QFP cross block. Hole QFP is the larger frozen-source driver at all four biases and recovers 99.15% of the joint-QFP improvement at -20 V. Electron and hole residual changes are 88.39% and 89.30% transport/SG-flux dominated. Both first coupled carrier updates point away from the Sentaurus QFP state."
        )},
        {"id": "carrier_finding", "type": "markdown", "sourceId": "source", "body": (
            "## Hole QFP is the larger frozen-source contributor\n\n"
            "Replacing only hole QFP has the larger absolute effect at all four biases. This is an attribution result under a frozen Vela state: it identifies which QFP family moves the integrated ionization source most, but does not by itself prove that the hole equation contains the defect."
        )},
        {"id": "recovery", "type": "chart", "chartId": "recovery_chart"},
        {"id": "term_finding", "type": "markdown", "sourceId": "terms", "body": (
            "## Transport dominates the carrier residual discrepancy\n\n"
            "At -20 V, transport/SG flux accounts for 88.39% of the electron residual-term change and 89.30% of the hole change; avalanche contributes 11.61% and 10.70%. Recombination is below 4.4e-8 of the term-delta sum, while boundary and gauge contributions are zero for the interior comparison."
        )},
        {"id": "terms", "type": "chart", "chartId": "term_chart"},
        {"id": "jacobian_finding", "type": "markdown", "sourceId": "jacobian", "body": (
            "## Dominant Jacobian blocks pass finite differences\n\n"
            "Poisson, transport, SG avalanche, and boundary/gauge blocks all pass the unchanged 5e-5 relative-error gate; the worst non-SRH value is 6.13e-8. The formal all-block outcome remains failed because SRH/Auger is evaluated at absolute matrix norms near 1e-15."
        )},
        {"id": "jacobian_values", "type": "table", "tableId": "jacobian_table"},
        {"id": "sensitivity_finding", "type": "markdown", "sourceId": "sensitivity", "body": (
            "## Step sensitivity identifies an SRH absolute noise floor\n\n"
            "Across 1e-5 to 1e-8 V perturbations, every SRH/Auger absolute discrepancy remains below 1e-13. The best check reaches 1.01e-21 absolute difference at a 3e-6 V step, while small steps amplify cancellation. The predeclared relative failure is retained, but it is not evidence of a dynamically relevant SRH derivative error."
        )},
        {"id": "sensitivity", "type": "chart", "chartId": "sensitivity_chart"},
        {"id": "update_finding", "type": "markdown", "sourceId": "updates", "body": (
            "## The first coupled update rejects both golden carrier QFPs\n\n"
            "At -20 V, electron and hole target projections are -0.9240 and -0.8892. Their trial residuals rise by about 4.19x and 3.68x, respectively, and both trial states move farther from the Sentaurus QFP target."
        )},
        {"id": "updates", "type": "table", "tableId": "updates_table"},
        {"id": "scope", "type": "markdown", "body": (
            "## Scope, definitions, and method\n\n"
            "The experiment uses the shared M2 mesh and fixed biases -18, -19.5, -19.7, and -20 V. Four frozen states isolate Vela baseline, electron QFP, hole QFP, and joint QFP. Carrier residuals are decomposed into transport, recombination, avalanche, gauge, and boundary terms. Analytic carrier-QFP and Poisson-QFP blocks are compared with double-symmetric finite differences on both baseline and joint-QFP states. Two complete runs are required to be byte-identical."
        )},
        {"id": "robustness", "type": "markdown", "sourceId": "determinism", "body": (
            "## Robustness and limitations\n\n"
            "All 148 audited outputs are byte-identical between independent runs, and term-sum closure is 1.32e-23. This is a local frozen-state and first-update diagnostic; it does not establish a unique code defect or validate a production correction. The SRH relative gate failure remains recorded because the gate was declared before the run."
        )},
        {"id": "next", "type": "markdown", "body": (
            "## Recommended next step\n\n"
            "Keep SG/Laux unchanged. The smallest discriminating experiment is an edge-level transport Jacobian audit at the hotspot support: separate mobility, Bernoulli/GSS coefficient, QFP driving force, and contact-row elimination contributions for electron and hole equations, then compare their analytic derivatives with finite differences before proposing any opt-in correction."
        )},
        {"id": "questions", "type": "markdown", "body": (
            "## Further questions\n\n"
            "Does the first wrong direction originate in the SG transport derivative itself, in row scaling/contact elimination, or only after the carrier blocks are coupled? Does the same hotspot and sign persist on M0 and under an independent perturbation basis?"
        )},
    ]
    manifest = {
        "version": 1, "surface": "report",
        "title": "M2 carrier-QFP residual and Jacobian verification",
        "generatedAt": "2026-08-01T00:14:00+08:00",
        "cards": [], "charts": charts, "tables": tables, "blocks": blocks,
        "sources": sources,
    }
    artifact = {
        "surface": "report", "manifest": manifest,
        "snapshot": {
            "version": 1, "generatedAt": "2026-08-01T00:14:00+08:00",
            "status": "ready",
            "datasets": {
                "recovery": recovery, "terms": terms, "sensitivity": sensitivity,
                "jacobian_summary": jacobian_summary, "updates": updates,
            },
        },
        "sources": sources,
    }
    write_json(output, artifact)


if __name__ == "__main__":
    main()
