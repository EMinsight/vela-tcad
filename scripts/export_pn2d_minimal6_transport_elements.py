#!/usr/bin/env python3
"""Export the exact 40-state Minimal6 native element transport fields."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
from pathlib import Path


TOPOLOGIES = ("mirror", "sketch")
BIASES = tuple(float(-value) for value in range(1, 21))
FIELD_SPECS = {
    "ElectricField": (2, "V*cm^-1"),
    "eGradQuasiFermi": (2, "V*cm^-1"),
    "hGradQuasiFermi": (2, "V*cm^-1"),
    "eMobility": (1, "cm^2*V^-1*s^-1"),
    "hMobility": (1, "cm^2*V^-1*s^-1"),
    "eCurrentDensity": (2, "A*cm^-2"),
    "hCurrentDensity": (2, "A*cm^-2"),
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
    components, unit = FIELD_SPECS[name]
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
        "components": components,
        "values": 4,
        "raw_value_count": 4 * components,
        "global_node_mapping": "region_cell_order",
        "mapping_status": "complete",
        "unit": unit,
    }
    for key, value in expected.items():
        if field.get(key) != value:
            raise ValueError(
                f"{name} cell field {key}={field.get(key)!r}, "
                f"expected {value!r}"
            )
    csv_file = field.get("csv_file")
    if not isinstance(csv_file, str) or not csv_file.endswith("_cells.csv"):
        raise ValueError(f"{name} cell field lacks a cell CSV")
    return field


def read_cell_values(
    path: Path, components: int
) -> list[tuple[int, tuple[float, ...]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    result = [
        (
            int(row["cell_id"]),
            tuple(float(row[f"component{index}"]) for index in range(components)),
        )
        for row in rows
    ]
    if [row[0] for row in result] != [0, 1, 2, 3]:
        raise ValueError(f"{path} does not contain canonical cells 0..3")
    return result


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError("refusing to write an empty transport CSV")
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
            cell_rows = {
                cell: {
                    "topology": topology,
                    "bias_V": bias,
                    "cell_id": cell,
                    "source_tdr_sha256": digest,
                }
                for cell in range(4)
            }
            state_fields: dict[str, object] = {}
            for name, (components, _) in FIELD_SPECS.items():
                field = cell_field(manifest, name)
                values_path = export / "fields" / str(field["csv_file"])
                values = read_cell_values(values_path, components)
                state_fields[name] = {
                    "field_index": field["index"],
                    "csv": values_path.relative_to(output_root).as_posix(),
                    "csv_sha256": sha256(values_path),
                }
                for cell, components_value in values:
                    for index, value in enumerate(components_value):
                        cell_rows[cell][f"{name}_component{index}_raw"] = format(
                            value, ".17g"
                        )
            for cell in range(4):
                row = cell_rows[cell]
                row.update(
                    {
                        "electric_field_x_V_per_m": format(
                            float(row["ElectricField_component0_raw"]) * 100.0,
                            ".17g",
                        ),
                        "electric_field_y_V_per_m": format(
                            float(row["ElectricField_component1_raw"]) * 100.0,
                            ".17g",
                        ),
                        "electron_grad_qf_x_V_per_m": format(
                            float(row["eGradQuasiFermi_component0_raw"])
                            * 100.0,
                            ".17g",
                        ),
                        "electron_grad_qf_y_V_per_m": format(
                            float(row["eGradQuasiFermi_component1_raw"])
                            * 100.0,
                            ".17g",
                        ),
                        "hole_grad_qf_x_V_per_m": format(
                            float(row["hGradQuasiFermi_component0_raw"])
                            * 100.0,
                            ".17g",
                        ),
                        "hole_grad_qf_y_V_per_m": format(
                            float(row["hGradQuasiFermi_component1_raw"])
                            * 100.0,
                            ".17g",
                        ),
                        "electron_mobility_m2_per_Vs": format(
                            float(row["eMobility_component0_raw"]) * 1.0e-4,
                            ".17g",
                        ),
                        "hole_mobility_m2_per_Vs": format(
                            float(row["hMobility_component0_raw"]) * 1.0e-4,
                            ".17g",
                        ),
                        "electron_current_x_A_per_m2": format(
                            float(row["eCurrentDensity_component0_raw"])
                            * 1.0e4,
                            ".17g",
                        ),
                        "electron_current_y_A_per_m2": format(
                            float(row["eCurrentDensity_component1_raw"])
                            * 1.0e4,
                            ".17g",
                        ),
                        "hole_current_x_A_per_m2": format(
                            float(row["hCurrentDensity_component0_raw"])
                            * 1.0e4,
                            ".17g",
                        ),
                        "hole_current_y_A_per_m2": format(
                            float(row["hCurrentDensity_component1_raw"])
                            * 1.0e4,
                            ".17g",
                        ),
                    }
                )
                samples.append(row)
            states.append(
                {
                    "topology": topology,
                    "bias_V": bias,
                    "tdr": relative,
                    "tdr_sha256": digest,
                    "field_manifest": (
                        manifest_path.relative_to(output_root).as_posix()
                    ),
                    "field_manifest_sha256": sha256(manifest_path),
                    "fields": state_fields,
                }
            )

    if len(states) != 40 or len(samples) != 160:
        raise ValueError("transport element export count differs from 40 x 4")
    samples_path = output_root / "transport_element_values.csv"
    write_csv(samples_path, samples)
    result: dict[str, object] = {
        "schema_version": 1,
        "status": "valid",
        "experiment": "pn2d_minimal6_native_element_transport_fields",
        "state_count": len(states),
        "sample_count": len(samples),
        "support": {
            "kind": "cell",
            "location_type": 3,
            "cell_count_per_state": 4,
            "mapping": "region_cell_order",
            "fields": FIELD_SPECS,
            "unavailable_native_cell_fields": [
                "ElectrostaticPotential",
                "eDensity",
                "hDensity",
                "eQuasiFermiPotential",
                "hQuasiFermiPotential",
            ],
        },
        "inputs": {
            "root": str(input_root),
            "tdr_sha256": input_hashes,
            "importer": str(importer),
            "importer_sha256": sha256(importer),
        },
        "states": states,
        "outputs": {
            "transport_element_values_csv": samples_path.name,
            "transport_element_values_sha256": sha256(samples_path),
        },
    }
    manifest_path = output_root / "manifest.json"
    manifest_path.write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-root", type=Path, required=True)
    parser.add_argument("--importer", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
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
