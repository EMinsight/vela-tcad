#!/usr/bin/env python3
"""Independently verify the corrected Sentaurus avalanche control matrix."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from typing import Any


EXPECTED_SOURCE_VARIANTS = {
    "base": (
        "implicit_default",
        "explicit_grad_qf",
        "explicit_electric_field",
        "grad_qf_aval_dens_grad_qf",
    ),
    "contact": (
        "grad_qf_use_qf_contacts",
        "grad_qf_use_qf_contacts_aval_dens_grad_qf",
    ),
}
STATIC_BUNDLE_FILES = ("pn2d_minimal6.tdr", "models.par")
EXPECTED_PAIRS = {
    "implicit_default": "explicit_grad_qf",
    "explicit_electric_field": "explicit_grad_qf",
    "grad_qf_aval_dens_grad_qf": "explicit_grad_qf",
    "grad_qf_use_qf_contacts": "explicit_grad_qf",
    "grad_qf_use_qf_contacts_aval_dens_grad_qf": (
        "grad_qf_use_qf_contacts"
    ),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-root", type=Path, required=True)
    parser.add_argument("--contact-root", type=Path, required=True)
    parser.add_argument("--comparison-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="ascii") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"empty CSV: {path}")
    return rows


def close(left: float, right: float) -> bool:
    return math.isclose(left, right, rel_tol=1.0e-15, abs_tol=1.0e-300)


def verify(
    base_root: Path,
    contact_root: Path,
    comparison_root: Path,
    output: Path,
) -> dict[str, Any]:
    roots = {
        "base": base_root.resolve(),
        "contact": contact_root.resolve(),
    }
    comparison_root = comparison_root.resolve()
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=True)

    manifest_path = comparison_root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="ascii"))
    if (
        manifest["status"]
        != "valid_sentaurus_avalanche_corrected_control_comparison"
    ):
        raise ValueError("comparison manifest status is invalid")
    observed_pairs = {
        item["candidate_variant"]: item["reference_variant"]
        for item in manifest["pair_specs"]
    }
    if observed_pairs != EXPECTED_PAIRS:
        raise ValueError("pair specification mismatch")

    source_manifests = {
        label: json.loads(
            (root / "manifest.json").read_text(encoding="ascii")
        )
        for label, root in roots.items()
    }
    releases = {
        item["sentaurus_release"] for item in source_manifests.values()
    }
    if releases != {manifest["sentaurus_release"]}:
        raise ValueError("source Sentaurus release mismatch")
    expected_biases = tuple(float(value) for value in manifest["biases_V"])
    for label, variants in EXPECTED_SOURCE_VARIANTS.items():
        source = source_manifests[label]
        if source["status"] != "passed":
            raise ValueError(f"{label} source manifest is not passed")
        if tuple(source["variants"]) != variants:
            raise ValueError(f"{label} source variant matrix mismatch")
        if tuple(float(value) for value in source["biases_V"]) != expected_biases:
            raise ValueError(f"{label} source bias matrix mismatch")
        if set(source["topologies"]) != {"mirror", "sketch"}:
            raise ValueError(f"{label} source topology matrix mismatch")

    recomputed_static_hashes = {}
    for topology in ("mirror", "sketch"):
        recomputed_static_hashes[topology] = {}
        for name in STATIC_BUNDLE_FILES:
            values = set()
            for label, variants in EXPECTED_SOURCE_VARIANTS.items():
                for variant in variants:
                    item = source_manifests[label]["topologies"][
                        topology
                    ][variant]
                    if item["status"] != "passed":
                        raise ValueError(
                            f"{label}/{topology}/{variant} is not passed"
                        )
                    values.add(item["bundle_sha256"][name])
            if len(values) != 1:
                raise ValueError(
                    f"{topology} static bundle mismatch for {name}"
                )
            recomputed_static_hashes[topology][name] = values.pop()
    if (
        recomputed_static_hashes
        != manifest["static_bundle_sha256_by_topology"]
    ):
        raise ValueError("reported static bundle hashes mismatch")

    for name, expected in manifest["input_sha256"].items():
        root_label, relative = name.split("/", 1)
        actual = sha256(roots[root_label] / Path(relative))
        if actual != expected:
            raise ValueError(f"input hash mismatch: {name}")
    for name, expected in manifest["output_sha256"].items():
        actual = sha256(comparison_root / name)
        if actual != expected:
            raise ValueError(f"output hash mismatch: {name}")

    states = load_csv(comparison_root / "state_summary.csv")
    quantities = load_csv(comparison_root / "quantity_comparison.csv")
    if len(states) != int(manifest["state_count"]):
        raise ValueError("state count mismatch")
    if len(quantities) != int(manifest["quantity_comparison_count"]):
        raise ValueError("quantity comparison count mismatch")

    biases = {float(value) for value in manifest["biases_V"]}
    expected_states = {
        (topology, bias, candidate, reference)
        for topology in ("mirror", "sketch")
        for bias in biases
        for candidate, reference in EXPECTED_PAIRS.items()
    }
    observed_states = {
        (
            row["topology"],
            float(row["bias_V"]),
            row["candidate_variant"],
            row["reference_variant"],
        )
        for row in states
    }
    if observed_states != expected_states:
        raise ValueError("state lattice mismatch")

    exact_by_variant = {}
    for candidate in EXPECTED_PAIRS:
        selected_states = [
            row for row in states if row["candidate_variant"] == candidate
        ]
        selected_quantities = [
            row
            for row in quantities
            if row["candidate_variant"] == candidate
        ]
        if any(
            row["reference_variant"] != EXPECTED_PAIRS[candidate]
            for row in selected_quantities
        ):
            raise ValueError(f"reference mismatch for {candidate}")
        exact = all(row["exact"] == "1" for row in selected_quantities)
        exact_by_variant[candidate] = exact
        if exact != all(
            row["all_parsed_values_exact"] == "1"
            for row in selected_states
        ):
            raise ValueError(f"state exact flag mismatch for {candidate}")

    if not exact_by_variant["implicit_default"]:
        raise ValueError("implicit default differs from explicit GradQF")
    if not exact_by_variant["explicit_electric_field"]:
        raise ValueError("ElectricField differs from contact-fallback GradQF")
    if exact_by_variant["grad_qf_use_qf_contacts"]:
        raise ValueError("forced QF-at-contacts control is not distinct")
    if exact_by_variant[
        "grad_qf_use_qf_contacts_aval_dens_grad_qf"
    ]:
        raise ValueError("forced-QF AvalDens control is not distinct")

    metric_names = [
        key for key in states[0] if key.startswith("max_")
    ]
    maxima = {}
    for candidate in EXPECTED_PAIRS:
        selected = [
            row for row in states if row["candidate_variant"] == candidate
        ]
        maxima[candidate] = {}
        for metric in metric_names:
            actual = max(abs(float(row[metric])) for row in selected)
            expected = float(
                manifest["maxima_by_variant"][candidate][metric]
            )
            if not close(actual, expected):
                raise ValueError(
                    f"maximum mismatch: {candidate}/{metric}"
                )
            maxima[candidate][metric] = actual

    result = {
        "schema_version": 1,
        "status": "independently_verified",
        "experiment": manifest["experiment"],
        "sentaurus_release": manifest["sentaurus_release"],
        "state_count": len(states),
        "quantity_comparison_count": len(quantities),
        "exact_match_by_variant": exact_by_variant,
        "static_bundle_sha256_by_topology": recomputed_static_hashes,
        "recomputed_maxima_by_variant": maxima,
        "input_sha256": {
            "manifest.json": sha256(manifest_path),
            "state_summary.csv": sha256(
                comparison_root / "state_summary.csv"
            ),
            "quantity_comparison.csv": sha256(
                comparison_root / "quantity_comparison.csv"
            ),
        },
    }
    (output / "verification.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n",
        encoding="ascii",
        newline="\n",
    )
    return result


def main() -> int:
    args = parse_args()
    result = verify(
        args.base_root,
        args.contact_root,
        args.comparison_root,
        args.output,
    )
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
