#!/usr/bin/env python3
"""Perturb one global-node field occurrence in a copied Sentaurus TDR."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
from pathlib import Path

import h5py


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def text_attribute(value: object) -> str:
    if isinstance(value, bytes):
        return value.decode("utf-8")
    return str(value)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--export-dir", type=Path, required=True)
    parser.add_argument("--field", required=True)
    parser.add_argument("--node", type=int, required=True)
    parser.add_argument("--delta", type=float, required=True)
    parser.add_argument("--manifest", type=Path)
    args = parser.parse_args()

    source_path = args.input.resolve()
    output_path = args.output.resolve()
    if source_path == output_path:
        raise ValueError("input and output TDR paths must differ")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source_path, output_path)

    field_manifest = json.loads(
        (args.export_dir / "field_manifest.json").read_text(encoding="utf-8")
    )
    occurrences: list[tuple[int, int]] = []
    for entry in field_manifest["fields"]:
        if entry["name"] != args.field or entry["mapping_status"] != "complete":
            continue
        with (args.export_dir / "fields" / entry["csv_file"]).open(
            newline="", encoding="utf-8"
        ) as stream:
            node_order = [int(row["node_id"]) for row in csv.DictReader(stream)]
        if args.node in node_order:
            occurrences.append((int(entry["region"]), node_order.index(args.node)))
    if not occurrences:
        raise ValueError(f"node {args.node} has no {args.field} occurrence")

    changes = []
    with h5py.File(output_path, "r+") as tdr:
        state = tdr["collection/geometry_0/state_0"]
        datasets = {}
        for key, group in state.items():
            if not key.startswith("dataset_") or "values" not in group:
                continue
            name = text_attribute(group.attrs.get("name", ""))
            region = int(group.attrs.get("region", -1))
            datasets[(name, region)] = group["values"]
        for region, offset in occurrences:
            values = datasets[(args.field, region)]
            before = float(values[offset])
            values[offset] = before + args.delta
            changes.append({
                "region": region,
                "offset": offset,
                "before": before,
                "after": before + args.delta,
            })
        tdr.flush()

    manifest = {
        "schema": "vela.sentaurus_tdr_node_perturbation.v1",
        "input": str(source_path),
        "input_sha256": digest(source_path),
        "output": str(output_path),
        "output_sha256": digest(output_path),
        "field": args.field,
        "node": args.node,
        "delta": args.delta,
        "changes": changes,
    }
    manifest_path = args.manifest or output_path.with_suffix(".manifest.json")
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
