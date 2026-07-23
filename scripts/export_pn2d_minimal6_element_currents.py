#!/usr/bin/env python3
"""Export the exact 40-state Minimal6 native element current vectors."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
from pathlib import Path


TOPOLOGIES = ("mirror", "sketch")
BIASES = tuple(float(-value) for value in range(1, 21))
CARRIERS = {
    "electron": "eCurrentDensity",
    "hole": "hCurrentDensity",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def cell_field(
    manifest: dict[str, object], name: str
) -> dict[str, object]:
    fields = manifest.get("fields")
    if not isinstance(fields, list):
        raise ValueError("field manifest lacks fields")
    matches = [
        field
        for field in fields
        if isinstance(field, dict)
        and field.get("name") == name
        and field.get("support_kind") == "cell"
    ]
    if len(matches) != 1:
        raise ValueError(f"expected one cell field {name}, got {len(matches)}")
    field = matches[0]
    expected = {
        "location_type": 3,
        "components": 2,
        "values": 4,
        "raw_value_count": 8,
        "global_node_mapping": "region_cell_order",
        "mapping_status": "complete",
        "unit": "A*cm^-2",
    }
    for key, value in expected.items():
        if field.get(key) != value:
            raise ValueError(
                f"{name} cell field {key}={field.get(key)!r}, expected {value!r}"
            )
    csv_file = field.get("csv_file")
    if not isinstance(csv_file, str) or not csv_file.endswith("_cells.csv"):
        raise ValueError(f"{name} cell field lacks a cell CSV")
    return field


def read_cell_vectors(path: Path) -> list[tuple[int, float, float]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    result = [
        (
            int(row["cell_id"]),
            float(row["component0"]),
            float(row["component1"]),
        )
        for row in rows
    ]
    if [row[0] for row in result] != [0, 1, 2, 3]:
        raise ValueError(f"{path} does not contain canonical cells 0..3")
    return result


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError("refusing to write an empty cell-current CSV")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def run(
    *, input_root: Path, importer: Path, output_root: Path
) -> dict[str, object]:
    input_root = input_root.resolve()
    importer = importer.resolve()
    output_root = output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    samples: list[dict[str, object]] = []
    input_hashes: dict[str, str] = {}
    states: list[dict[str, object]] = []

    for topology in TOPOLOGIES:
        for bias in BIASES:
            label = f"m{abs(int(bias))}V"
            filename = f"pn2d_minimal6_state_{label}.tdr"
            tdr = input_root / topology / label / filename
            if not tdr.is_file():
                raise ValueError(f"missing state TDR {tdr}")
            relative = tdr.relative_to(input_root).as_posix()
            digest = sha256(tdr)
            input_hashes[relative] = digest
            export = output_root / "states" / topology / label / "export"
            export.mkdir(parents=True, exist_ok=True)
            completed = subprocess.run(
                [
                    str(importer),
                    "--tdr",
                    str(tdr),
                    "--export-dir",
                    str(export),
                ],
                capture_output=True,
                text=True,
                check=False,
            )
            if completed.returncode != 0:
                raise RuntimeError(
                    f"import failed for {topology} {bias:g} V: "
                    f"{completed.stderr.strip()}"
                )
            manifest_path = export / "field_manifest.json"
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
            state_record: dict[str, object] = {
                "topology": topology,
                "bias_V": bias,
                "tdr": relative,
                "tdr_sha256": digest,
                "field_manifest": manifest_path.relative_to(output_root).as_posix(),
                "field_manifest_sha256": sha256(manifest_path),
                "carriers": {},
            }
            for carrier, field_name in CARRIERS.items():
                field = cell_field(manifest, field_name)
                vectors_path = export / "fields" / str(field["csv_file"])
                vectors = read_cell_vectors(vectors_path)
                state_record["carriers"][carrier] = {
                    "field_index": field["index"],
                    "csv": vectors_path.relative_to(output_root).as_posix(),
                    "csv_sha256": sha256(vectors_path),
                }
                for cell_id, x_cm, y_cm in vectors:
                    samples.append(
                        {
                            "topology": topology,
                            "bias_V": bias,
                            "carrier": carrier,
                            "cell_id": cell_id,
                            "current_x_A_per_cm2": format(x_cm, ".17g"),
                            "current_y_A_per_cm2": format(y_cm, ".17g"),
                            "current_x_A_per_m2": format(x_cm * 1.0e4, ".17g"),
                            "current_y_A_per_m2": format(y_cm * 1.0e4, ".17g"),
                            "source_tdr_sha256": digest,
                        }
                    )
            states.append(state_record)

    if len(states) != 40 or len(samples) != 320:
        raise ValueError("element-current export count differs from 40 x 2 x 4")
    samples_path = output_root / "element_current_vectors.csv"
    write_csv(samples_path, samples)
    manifest: dict[str, object] = {
        "schema_version": 1,
        "status": "valid",
        "experiment": "pn2d_minimal6_native_element_current_vectors",
        "state_count": len(states),
        "sample_count": len(samples),
        "support": {
            "kind": "cell",
            "location_type": 3,
            "cell_count_per_state": 4,
            "components": 2,
            "mapping": "region_cell_order",
        },
        "inputs": {
            "root": str(input_root),
            "tdr_sha256": input_hashes,
            "importer": str(importer),
            "importer_sha256": sha256(importer),
        },
        "states": states,
        "outputs": {
            "element_current_vectors_csv": samples_path.name,
            "element_current_vectors_sha256": sha256(samples_path),
        },
    }
    manifest_path = output_root / "manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--importer", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    manifest = run(
        input_root=args.input_root,
        importer=args.importer,
        output_root=args.output,
    )
    print(
        json.dumps(
            {
                "status": manifest["status"],
                "state_count": manifest["state_count"],
                "sample_count": manifest["sample_count"],
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
