#!/usr/bin/env python3
"""Independent verifier for the Minimal6 directed-edge SG inversion audit."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def quantile(values: list[float], fraction: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    position = (len(ordered) - 1) * fraction
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    return ordered[lower] + (ordered[upper] - ordered[lower]) * (
        position - lower
    )


def close(actual: str, expected: float | None, *, tolerance: float = 2.0e-12) -> bool:
    if expected is None:
        return actual == ""
    return math.isclose(float(actual), expected, rel_tol=tolerance, abs_tol=1.0e-14)


def verify(root: Path) -> dict[str, object]:
    manifest_path = root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    failures: list[str] = []
    for name, expected in manifest["outputs"].items():
        path = root / name
        if not path.is_file():
            failures.append(f"missing output {name}")
        elif sha256(path) != expected:
            failures.append(f"hash mismatch {name}")

    support = rows(root / "support_mapping_samples.csv")
    samples = rows(root / "sg_replacement_samples.csv")
    summaries = rows(root / "sg_replacement_summary.csv")
    contributions = rows(root / "replacement_contributions.csv")
    mobility = rows(root / "mobility_inversion_samples.csv")
    mobility_summaries = rows(root / "mobility_inversion_summary.csv")
    expected_counts = {
        "support": 720,
        "samples": 34560,
        "mobility": 4320,
    }
    if len(support) != expected_counts["support"]:
        failures.append("support sample count mismatch")
    if len(samples) != expected_counts["samples"]:
        failures.append("SG sample count mismatch")
    if len(mobility) != expected_counts["mobility"]:
        failures.append("mobility sample count mismatch")
    if any(
        float(row["p1_endpoint_identity_abs_difference_A_per_m2"]) != 0.0
        for row in support
    ):
        failures.append("P1 line and endpoint means differ")

    sample_groups: dict[tuple[str, ...], list[dict[str, str]]] = defaultdict(list)
    for row in samples:
        base = (
            row["support_mapping"],
            row["carrier"],
            row["formulation"],
            row["branch"],
        )
        sample_groups[base + (row["edge_scope"],)].append(row)
        sample_groups[base + ("all_edges",)].append(row)
        candidate = float(row["candidate_continuity_flux_per_m2_s"])
        reference = float(row["reference_continuity_flux_per_m2_s"])
        denominator = abs(candidate) + abs(reference)
        residual = 0.0 if denominator == 0.0 else abs(candidate - reference) / denominator
        if not math.isclose(
            residual,
            float(row["symmetric_relative_residual"]),
            rel_tol=2.0e-14,
            abs_tol=2.0e-15,
        ):
            failures.append("raw SG residual mismatch")
            break
    if len(sample_groups) != len(summaries):
        failures.append("SG summary group count mismatch")
    for row in summaries:
        key = (
            row["support_mapping"],
            row["carrier"],
            row["formulation"],
            row["branch"],
            row["edge_scope"],
        )
        selected = sample_groups.get(key, [])
        errors = [
            float(item["abs_log10_error_dex"])
            for item in selected
            if item["abs_log10_error_dex"] != ""
        ]
        signs = [
            float(item["sign_agreement"])
            for item in selected
            if item["sign_agreement"] != ""
        ]
        if int(row["sample_count"]) != len(selected):
            failures.append(f"SG summary count mismatch {key}")
        if not close(row["median_abs_log10_error_dex"], quantile(errors, 0.5)):
            failures.append(f"SG median mismatch {key}")
        if not close(row["p95_abs_log10_error_dex"], quantile(errors, 0.95)):
            failures.append(f"SG p95 mismatch {key}")
        expected_sign = sum(signs) / len(signs) if signs else None
        if not close(row["sign_agreement_fraction"], expected_sign):
            failures.append(f"SG sign fraction mismatch {key}")

    sample_index = {
        (
            row["topology"],
            row["bias_V"],
            row["carrier"],
            row["edge_id"],
            row["support_mapping"],
            row["formulation"],
            row["branch"],
        ): row
        for row in samples
    }
    paths = {
        "qf_sg": (
            ("vela_all", "sent_qf_only", "qf"),
            ("sent_qf_only", "sent_qf_and_mobility", "mobility"),
            ("sent_qf_and_mobility", "sent_all", "psi"),
        ),
        "density_sg": (
            ("vela_all", "sent_density_only", "density"),
            ("sent_density_only", "sent_density_and_mobility", "mobility"),
            ("sent_density_and_mobility", "sent_all", "psi"),
        ),
    }
    contribution_values: dict[tuple[str, ...], list[float]] = defaultdict(list)
    for row in samples:
        if row["branch"] != "vela_all":
            continue
        prefix = (
            row["topology"],
            row["bias_V"],
            row["carrier"],
            row["edge_id"],
            row["support_mapping"],
            row["formulation"],
        )
        for before, after, factor in paths[row["formulation"]]:
            left = sample_index[prefix + (before,)]
            right = sample_index[prefix + (after,)]
            if (
                left["abs_log10_error_dex"] == ""
                or right["abs_log10_error_dex"] == ""
            ):
                continue
            contribution_values[
                (
                    row["support_mapping"],
                    row["carrier"],
                    row["formulation"],
                    factor,
                )
            ].append(
                float(left["abs_log10_error_dex"])
                - float(right["abs_log10_error_dex"])
            )
    for row in contributions:
        key = (
            row["support_mapping"],
            row["carrier"],
            row["formulation"],
            row["replacement_step"],
        )
        values = contribution_values[key]
        if int(row["paired_sample_count"]) != len(values):
            failures.append(f"contribution count mismatch {key}")
        if not close(
            row["median_paired_error_reduction_dex"], quantile(values, 0.5)
        ):
            failures.append(f"contribution median mismatch {key}")

    mobility_groups: dict[tuple[str, ...], list[dict[str, str]]] = defaultdict(list)
    for row in mobility:
        key = (
            row["support_mapping"],
            row["carrier"],
            row["formulation"],
        )
        mobility_groups[key].append(row)
        reference = float(row["reference_continuity_flux_per_m2_s"])
        unit = float(row["unit_mobility_operator_flux_per_m2_s"])
        classification = row["classification"]
        if unit == 0.0:
            expected_classification = "zero_operator"
        elif reference / unit < 0.0:
            expected_classification = "sign_incompatible"
        else:
            expected_classification = "available"
        if classification != expected_classification:
            failures.append("mobility classification mismatch")
            break
        if classification == "available":
            required = reference / unit
            if not math.isclose(
                required,
                float(row["required_mobility_m2_per_Vs"]),
                rel_tol=2.0e-14,
                abs_tol=0.0,
            ):
                failures.append("required mobility mismatch")
                break
    for row in mobility_summaries:
        key = (
            row["support_mapping"],
            row["carrier"],
            row["formulation"],
        )
        selected = mobility_groups[key]
        available = [
            item for item in selected if item["classification"] == "available"
        ]
        dex = [
            abs(float(item["required_over_sentaurus_mobility_dex"]))
            for item in available
        ]
        if int(row["available_count"]) != len(available):
            failures.append(f"mobility available count mismatch {key}")
        if not close(
            row["median_abs_required_over_sentaurus_mobility_dex"],
            quantile(dex, 0.5),
        ):
            failures.append(f"mobility median mismatch {key}")

    policy = manifest["acceptance_policy"]
    if policy["native_edge_flux_available"] is not False:
        failures.append("native edge flux must remain unavailable")
    if policy["formula_change_authorized"] is not False:
        failures.append("formula change must remain unauthorized")
    result: dict[str, object] = {
        "schema_version": 1,
        "status": "passed" if not failures else "failed",
        "failure_count": len(failures),
        "failures": failures,
        "verified_counts": expected_counts,
        "verified_output_hash_count": len(manifest["outputs"]),
        "native_edge_flux_classification": manifest["support_audit"][
            "native_directed_edge_flux"
        ]["classification"],
        "formula_change_authorized": policy["formula_change_authorized"],
    }
    return result


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = args.root.resolve()
    result = verify(root)
    output = root / "independent_verification.json"
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, sort_keys=True))
    return 0 if result["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
