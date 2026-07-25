#!/usr/bin/env python3
"""Compare low-field-mobility Sentaurus avalanche driving-force controls."""

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
REFERENCE = "lowfield_mobility_avalanche_electric_field"
CANDIDATE = "lowfield_mobility_avalanche_grad_qf"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-root", type=Path, required=True)
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


def run(raw_root: Path, output: Path) -> dict[str, Any]:
    raw_root = raw_root.resolve()
    output = output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    source_manifest_path = raw_root / "manifest.json"
    source_manifest = json.loads(
        source_manifest_path.read_text(encoding="ascii")
    )
    if source_manifest["status"] != "passed":
        raise ValueError("source manifest is not passed")
    biases = tuple(float(value) for value in source_manifest["biases_V"])
    input_hashes = {
        "manifest.json": sha256(source_manifest_path),
    }
    all_rows = []
    summaries = []
    for topology in TOPOLOGIES:
        reference_groups, reference_currentplot, reference_paths = (
            load_variant(raw_root, topology, REFERENCE, biases)
        )
        candidate_groups, candidate_currentplot, candidate_paths = (
            load_variant(raw_root, topology, CANDIDATE, biases)
        )
        for path in (*reference_paths, *candidate_paths):
            input_hashes[path.relative_to(raw_root).as_posix()] = sha256(path)

        pair_rows = []
        for group in GROUP_IDENTITIES:
            pair_rows.extend(
                compare_group(
                    topology=topology,
                    candidate_variant=CANDIDATE,
                    group=group,
                    reference_rows=reference_groups[group],
                    candidate_rows=candidate_groups[group],
                )
            )
        pair_rows.extend(
            compare_currentplot(
                topology=topology,
                candidate_variant=CANDIDATE,
                reference_rows=reference_currentplot,
                candidate_rows=candidate_currentplot,
            )
        )
        for row in pair_rows:
            row["reference_variant"] = REFERENCE
        all_rows.extend(pair_rows)
        for bias in biases:
            row = summary_row(
                topology=topology,
                bias=bias,
                candidate_variant=CANDIDATE,
                long_rows=pair_rows,
                reference_elements=reference_groups["elements"],
                candidate_elements=candidate_groups["elements"],
            )
            row["reference_variant"] = REFERENCE
            summaries.append(row)

    quantity_path = output / "quantity_comparison.csv"
    state_path = output / "state_summary.csv"
    write_csv(quantity_path, all_rows)
    write_csv(state_path, summaries)
    metric_names = [key for key in summaries[0] if key.startswith("max_")]
    manifest = {
        "schema_version": 1,
        "status": "valid_sentaurus_avalanche_mobility_isolation",
        "experiment": "pn2d_minimal6_sentaurus_avalanche_mobility_isolation",
        "sentaurus_release": source_manifest["sentaurus_release"],
        "reference_variant": REFERENCE,
        "candidate_variant": CANDIDATE,
        "topologies": list(TOPOLOGIES),
        "biases_V": list(biases),
        "state_count": len(summaries),
        "quantity_comparison_count": len(all_rows),
        "candidate_distinct": any(
            row["all_parsed_values_exact"] == 0 for row in summaries
        ),
        "maxima": {
            metric: max(float(row[metric]) for row in summaries)
            for metric in metric_names
        },
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
    result = run(args.raw_root, args.output)
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
