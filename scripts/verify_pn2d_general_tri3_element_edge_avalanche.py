#!/usr/bin/env python3
"""Independently verify general-Tri3 Sentaurus avalanche evidence."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path
from typing import Any


SCHEMA = "pn2d_general_tri3_element_edge_avalanche/v1"
RELEASE = "O-2018.06-SP2"
BIASES = (-1.0, -10.0, -20.0)
VARIANTS = (
    "implicit_default",
    "explicit_grad_qf",
    "explicit_electric_field",
    "grad_qf_use_qf_contacts",
    "electric_field_use_qf_contacts",
    "grad_qf_aval_dens_grad_qf",
    "lowfield_mobility_avalanche_electric_field",
    "lowfield_mobility_avalanche_grad_qf",
)
PREFIX_GROUPS = {
    "AVAL_PROBE_BEGIN": "begin",
    "AVAL_PROBE_VERTEX": "vertex",
    "AVAL_PROBE_ELEMENT": "element",
    "AVAL_PROBE_MEASURE": "measure",
    "AVAL_PROBE_EDGE": "edge",
    "AVAL_PROBE_INTEGRAL": "integral",
    "AVAL_PROBE_END": "end",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-root", type=Path, required=True)
    parser.add_argument("--analysis-root", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="ascii") as stream:
        return list(csv.DictReader(stream))


def tokens(line: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for token in line.split()[1:]:
        if "=" in token:
            key, value = token.split("=", 1)
            result[key] = value
    return result


def verify_log(path: Path) -> dict[str, Any]:
    records: dict[str, list[dict[str, str]]] = {
        group: [] for group in PREFIX_GROUPS.values()
    }
    with path.open(encoding="ascii", errors="strict") as stream:
        for raw_line in stream:
            line = raw_line.strip()
            for prefix, group in PREFIX_GROUPS.items():
                if line.startswith(prefix + " "):
                    records[group].append(tokens(line))
                    break
    observed_biases = tuple(
        sorted(float(row["bias_V"]) for row in records["begin"])
    )
    require(observed_biases == tuple(sorted(BIASES)), f"{path}: bias mismatch")
    for bias in BIASES:
        begin = [
            row for row in records["begin"]
            if float(row["bias_V"]) == bias
        ]
        require(len(begin) == 1, f"{path}: begin count at {bias}")
        declared = begin[0]
        expected = {
            "vertex": int(declared["vertices"]),
            "element": int(declared["elements"]),
            "measure": int(declared["element_vertices"]),
            "edge": 3 * int(declared["elements"]),
            "integral": 1,
            "end": 1,
        }
        for group, count in expected.items():
            observed = sum(
                float(row["bias_V"]) == bias for row in records[group]
            )
            require(
                observed == count,
                f"{path}: {group} count {observed} != {count} at {bias}",
            )
    return {
        "biases": list(observed_biases),
        "vertices_per_state": int(records["begin"][0]["vertices"]),
        "elements_per_state": int(records["begin"][0]["elements"]),
    }


def verify_raw(raw_root: Path) -> dict[str, Any]:
    manifest_path = raw_root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="ascii"))
    require(manifest.get("schema") == SCHEMA, "raw schema mismatch")
    require(manifest.get("status") == "passed", "raw status mismatch")
    require(manifest.get("sentaurus_release") == RELEASE, "release mismatch")
    require(
        tuple(float(value) for value in manifest.get("exact_biases_V", ()))
        == BIASES,
        "raw exact bias mismatch",
    )
    require(tuple(manifest.get("variants", ())) == VARIANTS, "variant mismatch")
    cases = manifest.get("cases", {})
    require(len(cases) == 1, f"expected one raw case, got {tuple(cases)}")
    case_name = next(iter(cases))
    require(tuple(cases[case_name]) == tuple(sorted(VARIANTS)), "case variants")
    static = manifest["case_hashes"][case_name]
    require(set(static) >= {"tdr", "models.par"}, "static hashes incomplete")

    logs: dict[str, Any] = {}
    for variant in VARIANTS:
        result = cases[case_name][variant]
        require(result.get("status") == "passed", f"{variant}: not passed")
        require(
            tuple(float(value) for value in result["observed_biases_V"])
            == BIASES,
            f"{variant}: observed bias mismatch",
        )
        variant_root = raw_root / case_name / variant
        for name, expected in result["bundle_sha256"].items():
            path = variant_root / "bundle" / name
            require(path.is_file(), f"{variant}: missing bundle {name}")
            require(digest(path) == expected, f"{variant}: bundle hash {name}")
        require(
            result["bundle_sha256"]["pn2d_msh.tdr"] == static["tdr"],
            f"{variant}: TDR differs",
        )
        require(
            result["bundle_sha256"]["models.par"] == static["models.par"],
            f"{variant}: models.par differs",
        )
        for name, expected in result["output_sha256"].items():
            path = variant_root / "fetched" / name
            require(path.is_file(), f"{variant}: missing fetched {name}")
            require(digest(path) == expected, f"{variant}: output hash {name}")
        log = variant_root / "fetched" / f"run_{variant}.out"
        logs[variant] = verify_log(log)
    return {
        "case_name": case_name,
        "manifest_sha256": digest(manifest_path),
        "variant_count": len(cases[case_name]),
        "state_count": len(cases[case_name]) * len(BIASES),
        "log_contracts": logs,
    }


def angle_class(angles: list[float]) -> str:
    maximum = max(angles)
    if abs(maximum - 90.0) <= 1.0e-8:
        return "right"
    return "acute" if maximum < 90.0 else "obtuse"


def verify_analysis(
    analysis_root: Path,
    raw_summary: dict[str, Any],
) -> dict[str, Any]:
    manifest_path = analysis_root / "analysis_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="ascii"))
    require(manifest.get("schema") == SCHEMA, "analysis schema mismatch")
    require(manifest.get("status") == "valid", "analysis status mismatch")
    require(manifest.get("case_name") == raw_summary["case_name"], "case mismatch")
    require(
        manifest.get("source_manifest_sha256")
        == raw_summary["manifest_sha256"],
        "analysis/raw manifest hash mismatch",
    )
    for name, expected in manifest["outputs"].items():
        path = analysis_root / name
        require(path.is_file(), f"missing analysis output {name}")
        require(digest(path) == expected, f"analysis output hash {name}")

    geometry = read_csv(analysis_root / "mesh_classification.csv")
    require(len(geometry) == manifest["element_count"], "element row mismatch")
    classes: Counter[str] = Counter()
    has_acute_scalene = False
    has_obtuse_truncation = False
    for row in geometry:
        angles = [
            float(row["angle0_deg"]),
            float(row["angle1_deg"]),
            float(row["angle2_deg"]),
        ]
        require(abs(sum(angles) - 180.0) <= 1.0e-8, "angle sum mismatch")
        classification = angle_class(angles)
        require(row["angle_class"] == classification, "angle class mismatch")
        area = float(row["signed_area_um2"])
        require(area != 0.0, "degenerate element")
        require(
            row["orientation"] == ("ccw" if area > 0.0 else "cw"),
            "orientation mismatch",
        )
        positive = int(row["positive_support_count"])
        zero = int(row["zero_support_count"])
        negative = int(row["negative_support_count"])
        require(positive + zero + negative == 3, "support count mismatch")
        classes[classification] += 1
        has_acute_scalene |= (
            classification == "acute"
            and row["scalene"] == "1"
            and positive == 3
        )
        has_obtuse_truncation |= (
            classification == "obtuse" and (zero + negative) >= 1
        )

    counts = read_csv(analysis_root / "state_counts.csv")
    require(
        len(counts) == raw_summary["state_count"],
        "analysis state count mismatch",
    )
    comparisons = read_csv(analysis_root / "driver_control_summary.csv")
    require(len(comparisons) == 6 * len(BIASES), "driver comparison count")
    element_comparisons = read_csv(
        analysis_root / "driver_element_summary.csv"
    )
    require(
        len(element_comparisons) == 6 * len(BIASES) * len(geometry),
        "driver element comparison count",
    )
    class_summaries = read_csv(
        analysis_root / "driver_contact_class_summary.csv"
    )
    require(
        len(class_summaries) == 6 * len(BIASES) * 2,
        "driver contact/interior summary count",
    )
    require(
        {row["element_class"] for row in class_summaries}
        == {"contact", "interior"},
        "driver contact/interior classes incomplete",
    )
    case_name = raw_summary["case_name"]
    if case_name == "skewed_tri3":
        require(has_acute_scalene, "skewed mesh lacks acute scalene support")
    if case_name == "skewed_tri3_constrained":
        require(has_obtuse_truncation, "constrained mesh lacks obtuse support")
    return {
        "element_count": len(geometry),
        "angle_class_counts": dict(sorted(classes.items())),
        "has_acute_scalene_three_positive_supports": has_acute_scalene,
        "has_obtuse_truncation": has_obtuse_truncation,
        "scientific_role": (
            "diagnostic_only"
            if case_name == "skewed_tri3_constrained"
            else "device_physics_oracle"
        ),
    }


def main() -> int:
    args = parse_args()
    raw_root = args.raw_root.resolve()
    analysis_root = args.analysis_root.resolve()
    raw_summary = verify_raw(raw_root)
    analysis_summary = verify_analysis(analysis_root, raw_summary)
    result = {
        "schema": SCHEMA,
        "status": "verified",
        "raw": raw_summary,
        "analysis": analysis_summary,
    }
    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="ascii", newline="\n")
    print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
