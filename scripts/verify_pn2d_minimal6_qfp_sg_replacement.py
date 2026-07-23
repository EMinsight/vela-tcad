#!/usr/bin/env python3
"""Independent structural and numerical verification of the QFP replacement artifact."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def quantile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * fraction
    lower, upper = math.floor(position), math.ceil(position)
    if lower == upper:
        return ordered[lower]
    return ordered[lower] + (ordered[upper] - ordered[lower]) * (position - lower)


def close(actual: float, expected: float) -> None:
    if not math.isclose(actual, expected, rel_tol=2.0e-13, abs_tol=2.0e-15):
        raise ValueError(f"summary mismatch: {actual} != {expected}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    args = parser.parse_args()
    manifest = json.loads((args.root / "manifest.json").read_text(encoding="utf-8"))
    if manifest["status"] != "valid":
        raise ValueError("artifact manifest is not valid")
    if manifest["state_contract"]["state_count"] != 40:
        raise ValueError("manifest does not declare exactly 40 states")
    replay = manifest["baseline_cpp_replay"]
    if replay["status"] != "passed" or replay["sample_count"] != 720:
        raise ValueError("baseline C++ replay gate did not pass 720 samples")
    if replay["max_relative_error"] > replay["tolerance_max_relative_error"]:
        raise ValueError("baseline C++ replay exceeds its declared tolerance")

    outputs = manifest["outputs"]
    for name_key, hash_key in (
        ("edge_samples_csv", "edge_samples_sha256"),
        ("summary_csv", "summary_sha256"),
        ("report_md", "report_sha256"),
    ):
        path = args.root / outputs[name_key]
        if sha256(path) != outputs[hash_key]:
            raise ValueError(f"output hash mismatch for {path.name}")
    with (args.root / outputs["edge_samples_csv"]).open(
        newline="", encoding="utf-8"
    ) as handle:
        edges = list(csv.DictReader(handle))
    with (args.root / outputs["summary_csv"]).open(
        newline="", encoding="utf-8"
    ) as handle:
        summaries = list(csv.DictReader(handle))
    if len(edges) != 3600 or len(summaries) != 10:
        raise ValueError(f"expected 3600 edge rows and 10 summaries")
    states = {(row["topology"], float(row["bias_V"])) for row in edges}
    expected_states = {
        (topology, float(-bias))
        for topology in ("sketch", "mirror")
        for bias in range(1, 21)
    }
    if states != expected_states:
        raise ValueError("edge samples do not contain the exact 40 states")

    by_key = {
        (
            row["topology"],
            float(row["bias_V"]),
            row["carrier"],
            int(row["edge_id"]),
            row["variant"],
        ): row
        for row in edges
    }
    if len(by_key) != len(edges):
        raise ValueError("edge sample keys are not unique")
    for topology, bias in expected_states:
        for carrier in ("electron", "hole"):
            for edge_id in range(9):
                baseline = by_key[(topology, bias, carrier, edge_id, "baseline")]
                no_op = "hole_qfp" if carrier == "electron" else "electron_qfp"
                relevant = "electron_qfp" if carrier == "electron" else "hole_qfp"
                if (
                    by_key[(topology, bias, carrier, edge_id, no_op)][
                        "candidate_continuity_flux_per_m2_s"
                    ]
                    != baseline["candidate_continuity_flux_per_m2_s"]
                ):
                    raise ValueError("opposite-carrier QFP replacement is not an exact no-op")
                if (
                    by_key[(topology, bias, carrier, edge_id, relevant)][
                        "candidate_continuity_flux_per_m2_s"
                    ]
                    != by_key[(topology, bias, carrier, edge_id, "both_qfp")][
                        "candidate_continuity_flux_per_m2_s"
                    ]
                ):
                    raise ValueError("both-QFP branch differs from relevant single-QFP branch")

    summary_index = {(row["carrier"], row["variant"]): row for row in summaries}
    for carrier in ("electron", "hole"):
        baseline_rows = [
            row
            for row in edges
            if row["carrier"] == carrier
            and row["variant"] == "baseline"
            and row["affected_edge"] == "1"
        ]
        if len(baseline_rows) != 280:
            raise ValueError("each carrier must have 280 affected baseline edges")
        baseline_errors = {
            (row["topology"], row["bias_V"], row["edge_id"]): float(
                row["abs_log10_error"]
            )
            for row in baseline_rows
            if row["abs_log10_error"] != ""
        }
        for variant in (
            "baseline",
            "electron_qfp",
            "hole_qfp",
            "both_qfp",
            "strict_frozen_density",
        ):
            selected = [
                row
                for row in edges
                if row["carrier"] == carrier
                and row["variant"] == variant
                and row["affected_edge"] == "1"
            ]
            target = summary_index[(carrier, variant)]
            errors = [
                float(row["abs_log10_error"])
                for row in selected
                if row["abs_log10_error"] != ""
            ]
            residuals = [float(row["symmetric_relative_residual"]) for row in selected]
            signs = [
                float(row["sign_agreement"])
                for row in selected
                if row["sign_agreement"] != ""
            ]
            improvements = [
                baseline_errors[(row["topology"], row["bias_V"], row["edge_id"])]
                - float(row["abs_log10_error"])
                for row in selected
                if row["abs_log10_error"] != ""
                and (row["topology"], row["bias_V"], row["edge_id"])
                in baseline_errors
            ]
            close(float(target["median_abs_log10_error_dex"]), quantile(errors, 0.5))
            close(float(target["p95_abs_log10_error_dex"]), quantile(errors, 0.95))
            close(
                float(target["median_symmetric_relative_residual"]),
                quantile(residuals, 0.5),
            )
            close(float(target["sign_agreement_fraction"]), statistics.fmean(signs))
            close(
                float(target["median_paired_log_error_improvement_dex"]),
                quantile(improvements, 0.5),
            )

    electron = summary_index[("electron", "electron_qfp")]
    hole = summary_index[("hole", "hole_qfp")]
    if not (
        float(electron["median_abs_log10_error_dex"]) > 1.0
        and float(hole["median_abs_log10_error_dex"]) > 1.0
        and float(electron["median_symmetric_relative_residual"]) > 0.9
        and float(hole["median_symmetric_relative_residual"]) > 0.9
    ):
        raise ValueError("replacement unexpectedly satisfies the non-closure rejection rule")
    result = {
        "status": "passed",
        "state_count": len(states),
        "edge_row_count": len(edges),
        "summary_row_count": len(summaries),
        "scientific_classification": "sign_corrected_but_magnitude_not_closed",
    }
    verification_path = args.root / "independent_verification.json"
    verification_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
