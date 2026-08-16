#!/usr/bin/env python3
"""Copy a Sentaurus TDR and substitute one regionwise nodal field's values."""

from __future__ import annotations

import argparse
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
    parser.add_argument(
        "--source-field", action="append", required=True,
        help="source field to copy; repeat to substitute their componentwise sum",
    )
    parser.add_argument("--target-field", required=True)
    parser.add_argument("--manifest", type=Path)
    args = parser.parse_args()

    source_path = args.input.resolve()
    output_path = args.output.resolve()
    if source_path == output_path:
        raise ValueError("input and output TDR paths must differ")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source_path, output_path)

    substitutions = []
    with h5py.File(output_path, "r+") as tdr:
        state = tdr["collection/geometry_0/state_0"]
        by_name_region = {}
        for key, group in state.items():
            if not key.startswith("dataset_") or "values" not in group:
                continue
            name = text_attribute(group.attrs.get("name", ""))
            region = int(group.attrs.get("region", -1))
            by_name_region[(name, region)] = group
        source_region_sets = [
            {region for name, region in by_name_region if name == source_name}
            for source_name in args.source_field
        ]
        source_regions = source_region_sets[0]
        if any(regions != source_regions for regions in source_region_sets[1:]):
            raise ValueError(
                "source field region mismatch: "
                + ", ".join(
                    f"{name}={sorted(regions)}"
                    for name, regions in zip(args.source_field, source_region_sets)
                )
            )
        target_regions = {
            region for name, region in by_name_region if name == args.target_field
        }
        if source_regions != target_regions:
            raise ValueError(
                f"field region mismatch: source={sorted(source_regions)}, "
                f"target={sorted(target_regions)}"
            )
        for region in sorted(source_regions):
            sources = [
                by_name_region[(source_name, region)]["values"]
                for source_name in args.source_field
            ]
            target = by_name_region[(args.target_field, region)]["values"]
            if any(source.shape != target.shape for source in sources):
                raise ValueError(
                    f"region {region} shape mismatch: "
                    f"{[source.shape for source in sources]} vs {target.shape}"
                )
            target[...] = sources[0][...]
            for source in sources[1:]:
                target[...] += source[...]
            substitutions.append({"region": region, "values": int(target.size)})
        tdr.flush()

    manifest = {
        "schema": "vela.sentaurus_tdr_field_substitution.v1",
        "input": str(source_path),
        "input_sha256": digest(source_path),
        "output": str(output_path),
        "output_sha256": digest(output_path),
        "source_fields": args.source_field,
        "target_field": args.target_field,
        "substitutions": substitutions,
    }
    manifest_path = args.manifest or output_path.with_suffix(".manifest.json")
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
