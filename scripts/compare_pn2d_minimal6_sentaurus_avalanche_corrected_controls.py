#!/usr/bin/env python3
"""Compare default, contact-forced, and AvalDens Sentaurus controls."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.compare_pn2d_minimal6_sentaurus_avalanche_drive_controls import (
    GROUP_IDENTITIES,
    compare_currentplot,
    compare_group,
    currentplot_rows,
    sha256,
    summary_row,
    write_csv,
)
from scripts.diagnose_pn2d_minimal6_element_avalanche_replay import parse_log


TOPOLOGIES = ("mirror", "sketch")
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
PAIR_SPECS = (
    (
        "implicit_default",
        "base",
        "explicit_grad_qf",
        "base",
    ),
    (
        "explicit_electric_field",
        "base",
        "explicit_grad_qf",
        "base",
    ),
    (
        "grad_qf_aval_dens_grad_qf",
        "base",
        "explicit_grad_qf",
        "base",
    ),
    (
        "grad_qf_use_qf_contacts",
        "contact",
        "explicit_grad_qf",
        "base",
    ),
    (
        "grad_qf_use_qf_contacts_aval_dens_grad_qf",
        "contact",
        "grad_qf_use_qf_contacts",
        "contact",
    ),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-root", type=Path, required=True)
    parser.add_argument("--contact-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def load_variant(
    root: Path,
    topology: str,
    variant: str,
    biases: tuple[float, ...],
) -> tuple[
    dict[str, list[dict[str, Any]]],
    list[dict[str, float]],
    tuple[Path, Path],
]:
    fetched = root / topology / variant / "fetched"
    log_path = fetched / f"run_{variant}.out"
    plt_path = fetched / f"runtime_element_avalanche_probe_{variant}.plt"
    return (
        parse_log(log_path, biases),
        currentplot_rows(plt_path, biases),
        (log_path, plt_path),
    )


def validate_source_manifests(
    manifests: dict[str, dict[str, Any]],
) -> tuple[str, dict[str, dict[str, str]]]:
    releases = {item["sentaurus_release"] for item in manifests.values()}
    if len(releases) != 1:
        raise ValueError("base and contact Sentaurus releases differ")
    for label, expected_variants in EXPECTED_SOURCE_VARIANTS.items():
        manifest = manifests[label]
        if manifest["status"] != "passed":
            raise ValueError(f"{label} source manifest is not passed")
        if tuple(manifest["variants"]) != expected_variants:
            raise ValueError(f"{label} source variant matrix mismatch")
        if set(manifest["topologies"]) != set(TOPOLOGIES):
            raise ValueError(f"{label} source topology matrix mismatch")

    static_hashes = {}
    for topology in TOPOLOGIES:
        static_hashes[topology] = {}
        for name in STATIC_BUNDLE_FILES:
            values = set()
            for label, variants in EXPECTED_SOURCE_VARIANTS.items():
                for variant in variants:
                    result = manifests[label]["topologies"][topology][variant]
                    if result["status"] != "passed":
                        raise ValueError(
                            f"{label}/{topology}/{variant} is not passed"
                        )
                    values.add(result["bundle_sha256"][name])
            if len(values) != 1:
                raise ValueError(
                    f"{topology} static bundle hash mismatch for {name}"
                )
            static_hashes[topology][name] = values.pop()
    return releases.pop(), static_hashes


def run(base_root: Path, contact_root: Path, output: Path) -> dict[str, Any]:
    roots = {
        "base": base_root.resolve(),
        "contact": contact_root.resolve(),
    }
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    source_manifests = {
        label: json.loads(
            (root / "manifest.json").read_text(encoding="ascii")
        )
        for label, root in roots.items()
    }
    sentaurus_release, static_bundle_hashes = validate_source_manifests(
        source_manifests
    )
    base_biases = tuple(
        float(value) for value in source_manifests["base"]["biases_V"]
    )
    contact_biases = tuple(
        float(value) for value in source_manifests["contact"]["biases_V"]
    )
    if base_biases != contact_biases:
        raise ValueError("base and contact bias matrices differ")

    cache = {}
    input_hashes = {}
    for label, root in roots.items():
        manifest_path = root / "manifest.json"
        input_hashes[f"{label}/manifest.json"] = sha256(manifest_path)
    all_rows = []
    summaries = []
    for topology in TOPOLOGIES:
        for (
            candidate_variant,
            candidate_root_label,
            reference_variant,
            reference_root_label,
        ) in PAIR_SPECS:
            candidate_key = (
                candidate_root_label,
                topology,
                candidate_variant,
            )
            reference_key = (
                reference_root_label,
                topology,
                reference_variant,
            )
            for key in (candidate_key, reference_key):
                if key not in cache:
                    root_label, key_topology, variant = key
                    cache[key] = load_variant(
                        roots[root_label],
                        key_topology,
                        variant,
                        base_biases,
                    )
                    for path in cache[key][2]:
                        relative = path.relative_to(roots[root_label])
                        input_hashes[
                            f"{root_label}/{relative.as_posix()}"
                        ] = sha256(path)
            candidate_groups, candidate_currentplot, _ = cache[candidate_key]
            reference_groups, reference_currentplot, _ = cache[reference_key]

            pair_rows = []
            for group in GROUP_IDENTITIES:
                pair_rows.extend(
                    compare_group(
                        topology=topology,
                        candidate_variant=candidate_variant,
                        group=group,
                        reference_rows=reference_groups[group],
                        candidate_rows=candidate_groups[group],
                    )
                )
            pair_rows.extend(
                compare_currentplot(
                    topology=topology,
                    candidate_variant=candidate_variant,
                    reference_rows=reference_currentplot,
                    candidate_rows=candidate_currentplot,
                )
            )
            for row in pair_rows:
                row["reference_variant"] = reference_variant
            all_rows.extend(pair_rows)

            for bias in base_biases:
                row = summary_row(
                    topology=topology,
                    bias=bias,
                    candidate_variant=candidate_variant,
                    long_rows=pair_rows,
                    reference_elements=reference_groups["elements"],
                    candidate_elements=candidate_groups["elements"],
                )
                row["reference_variant"] = reference_variant
                summaries.append(row)

    quantity_path = output / "quantity_comparison.csv"
    state_path = output / "state_summary.csv"
    write_csv(quantity_path, all_rows)
    write_csv(state_path, summaries)

    summary_by_variant = {
        variant: [
            row for row in summaries if row["candidate_variant"] == variant
        ]
        for variant, _, _, _ in PAIR_SPECS
    }
    metric_names = [
        key
        for key in summaries[0]
        if key.startswith("max_")
    ]
    maxima = {
        variant: {
            metric: max(float(row[metric]) for row in rows)
            for metric in metric_names
        }
        for variant, rows in summary_by_variant.items()
    }
    exact_by_variant = {
        variant: all(row["all_parsed_values_exact"] == 1 for row in rows)
        for variant, rows in summary_by_variant.items()
    }
    manifest = {
        "schema_version": 1,
        "status": "valid_sentaurus_avalanche_corrected_control_comparison",
        "experiment": "pn2d_minimal6_sentaurus_avalanche_corrected_controls",
        "sentaurus_release": sentaurus_release,
        "topologies": list(TOPOLOGIES),
        "static_bundle_sha256_by_topology": static_bundle_hashes,
        "biases_V": list(base_biases),
        "pair_specs": [
            {
                "candidate_variant": candidate,
                "candidate_root": candidate_root,
                "reference_variant": reference,
                "reference_root": reference_root,
            }
            for candidate, candidate_root, reference, reference_root
            in PAIR_SPECS
        ],
        "state_count": len(summaries),
        "quantity_comparison_count": len(all_rows),
        "exact_match_by_variant": exact_by_variant,
        "default_matches_explicit_grad_qf": exact_by_variant[
            "implicit_default"
        ],
        "electric_field_matches_contact_fallback_grad_qf": exact_by_variant[
            "explicit_electric_field"
        ],
        "forced_qf_contacts_distinct": not exact_by_variant[
            "grad_qf_use_qf_contacts"
        ],
        "forced_qf_aval_dens_distinct": not exact_by_variant[
            "grad_qf_use_qf_contacts_aval_dens_grad_qf"
        ],
        "maxima_by_variant": maxima,
        "input_sha256": input_hashes,
        "output_sha256": {
            quantity_path.name: sha256(quantity_path),
            state_path.name: sha256(state_path),
        },
    }
    manifest_path = output / "manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="ascii",
        newline="\n",
    )
    return manifest


def main() -> int:
    args = parse_args()
    result = run(args.base_root, args.contact_root, args.output)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
