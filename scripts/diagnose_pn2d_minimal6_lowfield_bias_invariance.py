#!/usr/bin/env python3
"""Verify Minimal6 native low-field element mobility across bias."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError("refusing to write an empty CSV")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def load_cells(
    directory: Path, topology: str, bias: float
) -> tuple[list[dict[str, object]], dict[str, str]]:
    samples: list[dict[str, object]] = []
    hashes: dict[str, str] = {}
    for carrier, quantity in (
        ("electron", "eMobility"),
        ("hole", "hMobility"),
    ):
        path = directory / "fields" / f"{quantity}_region0_cells.csv"
        rows = read_csv(path)
        if [int(row["cell_id"]) for row in rows] != [0, 1, 2, 3]:
            raise ValueError(f"{path} lacks canonical cell ids")
        hashes[f"{topology}/{bias:g}/{quantity}"] = sha256(path)
        for row in rows:
            samples.append(
                {
                    "topology": topology,
                    "bias_V": bias,
                    "carrier": carrier,
                    "cell_id": int(row["cell_id"]),
                    "mobility_m2_per_Vs": float(row["component0"]) * 1.0e-4,
                    "source_sha256": sha256(path),
                }
            )
    return samples, hashes


def run(args: argparse.Namespace) -> dict[str, object]:
    control_root = args.control_root.resolve()
    manifest_path = control_root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("outputs_complete") is not True:
        raise ValueError("low-field control manifest is incomplete")
    states = manifest.get("states")
    if not isinstance(states, list) or len(states) != 4:
        raise ValueError("expected four -1/-10 V control states")

    samples: list[dict[str, object]] = []
    input_hashes: dict[str, str] = {
        "control_manifest": sha256(manifest_path)
    }
    for state in states:
        if state.get("status") != "passed":
            raise ValueError("all low-field control states must pass")
        state_samples, hashes = load_cells(
            Path(str(state["export_dir"])),
            str(state["topology_id"]),
            float(state["requested_bias_V"]),
        )
        samples.extend(state_samples)
        input_hashes.update(hashes)

    for topology, directory in (
        ("mirror", args.mirror_m20.resolve()),
        ("sketch", args.sketch_m20.resolve()),
    ):
        state_samples, hashes = load_cells(directory, topology, -20.0)
        samples.extend(state_samples)
        input_hashes.update(hashes)
    if len(samples) != 48:
        raise ValueError("expected 48 topology/bias/carrier/cell samples")

    reference = {
        (row["topology"], row["carrier"], row["cell_id"]): float(
            row["mobility_m2_per_Vs"]
        )
        for row in samples
        if float(row["bias_V"]) == -20.0
    }
    compared: list[dict[str, object]] = []
    for row in samples:
        key = (row["topology"], row["carrier"], row["cell_id"])
        value = float(row["mobility_m2_per_Vs"])
        baseline = reference[key]
        relative = abs(value - baseline) / max(abs(baseline), 1.0e-300)
        compared.append(
            {
                **row,
                "reference_bias_V": -20.0,
                "reference_mobility_m2_per_Vs": baseline,
                "relative_difference": relative,
            }
        )

    summary: list[dict[str, object]] = []
    for carrier in ("electron", "hole"):
        rows = [row for row in compared if row["carrier"] == carrier]
        relative = [float(row["relative_difference"]) for row in rows]
        summary.append(
            {
                "carrier": carrier,
                "sample_count": len(rows),
                "maximum_relative_difference": max(relative),
                "gate_threshold": 1.0e-10,
                "gate": "pass" if max(relative) <= 1.0e-10 else "fail",
            }
        )
    status = (
        "valid"
        if all(row["gate"] == "pass" for row in summary)
        else "bias_dependent_low_field_mobility"
    )

    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    write_csv(output / "lowfield_bias_samples.csv", compared)
    write_csv(output / "lowfield_bias_summary.csv", summary)
    report = [
        "# PN2D Minimal6 native low-field bias invariance",
        "",
        f"Status: `{status}`",
        "",
        "| Carrier | N | Maximum relative difference | Gate |",
        "|---|---:|---:|---|",
    ]
    for row in summary:
        report.append(
            f"| {row['carrier']} | {row['sample_count']} | "
            f"{float(row['maximum_relative_difference']):.6e} | "
            f"{row['gate']} |"
        )
    report.extend(
        [
            "",
            "The comparison covers mirror/sketch, -1/-10/-20 V, two "
            "carriers, and four native elements per state.",
            "",
        ]
    )
    (output / "report.md").write_text(
        "\n".join(report), encoding="utf-8", newline="\n"
    )
    outputs = {}
    for name in (
        "lowfield_bias_samples.csv",
        "lowfield_bias_summary.csv",
        "report.md",
    ):
        outputs[name] = sha256(output / name)
    result = {
        "schema_version": 1,
        "status": status,
        "experiment": "pn2d_minimal6_native_lowfield_bias_invariance",
        "control_state_count": 6,
        "sample_count": len(compared),
        "maximum_relative_difference": max(
            float(row["relative_difference"]) for row in compared
        ),
        "gate_threshold": 1.0e-10,
        "production_formula_modified": False,
        "inputs": input_hashes,
        "outputs": outputs,
    }
    (output / "manifest.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
        newline="\n",
    )
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--control-root", type=Path, required=True)
    parser.add_argument("--mirror-m20", type=Path, required=True)
    parser.add_argument("--sketch-m20", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(run(args), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
