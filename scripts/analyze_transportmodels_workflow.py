#!/usr/bin/env python3
"""Summarize a completed TransportModels DD or DG comparison by bias regime."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from statistics import median
from typing import Any


def read_curve(path: Path, current_column: str) -> list[dict[str, float]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"{path}: empty curve")
    return [{
        "bias_V": float(row["bias_V"]),
        "current": float(row[current_column]),
    } for row in rows]


def aligned_errors(reference: list[dict[str, float]],
                   candidate: list[dict[str, float]]) -> list[dict[str, float]]:
    if len(reference) != len(candidate):
        raise ValueError("reference and candidate point counts differ")
    result = []
    for ref, cand in zip(reference, candidate):
        if abs(ref["bias_V"] - cand["bias_V"]) > 1.0e-10:
            raise ValueError("reference and candidate bias lattices differ")
        ref_abs = abs(ref["current"])
        cand_abs = abs(cand["current"])
        result.append({
            "bias_V": ref["bias_V"],
            "reference_A_per_um": ref_abs,
            "vela_A_per_um": cand_abs,
            "relative_error": abs(cand_abs - ref_abs) / max(ref_abs, 1.0e-300),
            "absolute_log10_error_dex": abs(
                math.log10(max(cand_abs, 1.0e-300))
                - math.log10(max(ref_abs, 1.0e-300))),
        })
    return result


def region_metrics(rows: list[dict[str, float]]) -> dict[str, Any]:
    if not rows:
        raise ValueError("cannot summarize an empty region")
    relative = [row["relative_error"] for row in rows]
    dex = [row["absolute_log10_error_dex"] for row in rows]
    return {
        "points": len(rows),
        "bias_range_V": [rows[0]["bias_V"], rows[-1]["bias_V"]],
        "median_relative_error": median(relative),
        "max_relative_error": max(relative),
        "median_absolute_log10_error_dex": median(dex),
        "max_absolute_log10_error_dex": max(dex),
    }


def analyze(workflow_dir: Path, generated_dir: Path, branch: str,
            candidate_paths: dict[str, Path] | None = None) -> dict[str, Any]:
    candidate_paths = candidate_paths or {}
    manifest_path = workflow_dir / "workflow_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    result: dict[str, Any] = {
        "schema": "vela.transportmodels.bias_regime_analysis.v1",
        "branch": branch,
        "workflow_status": manifest.get("status"),
        "comparison_status": manifest.get("comparison_status"),
        "controlled_delta": manifest.get("controlled_delta", {}),
        "stages": [{
            "name": stage["name"],
            "status": stage.get("status"),
            "execution": stage.get("execution"),
            "config_sha256": stage.get("config_sha256"),
            "final_state_sha256": stage.get("final_state_sha256"),
        } for stage in manifest["stages"]],
    }
    curves: dict[str, Any] = {}
    for kind in ("idvg", "idvd"):
        candidate_path = candidate_paths.get(
            kind,
            workflow_dir / f"{branch}_{kind}_curve_comparison_candidate.csv",
        )
        reference = read_curve(
            generated_dir / "reference_curves" /
            f"transportmodels_sentaurus2022_{branch}_{kind}_reference.csv",
            "current_total",
        )
        candidate = read_curve(
            candidate_path,
            "current_total_A_per_um",
        )
        rows = aligned_errors(reference, candidate)
        curves[kind] = {
            "candidate_path": str(candidate_path.resolve()),
            "points": rows,
            "endpoint": rows[-1],
            "overall": region_metrics(rows),
        }
        if kind == "idvg":
            groups = {
                "off": [row for row in rows if row["bias_V"] <= -0.68],
                "transition": [row for row in rows
                               if -0.52 <= row["bias_V"] <= 0.12],
                "on": [row for row in rows if row["bias_V"] >= 0.28],
            }
            curves[kind]["regions"] = {
                name: region_metrics(group) for name, group in groups.items()
            }
        else:
            curves[kind]["nonzero_bias"] = region_metrics(
                [row for row in rows if abs(row["bias_V"]) > 1.0e-12])
    result["curves"] = curves
    return result


def markdown(result: dict[str, Any]) -> str:
    idvg = result["curves"]["idvg"]
    idvd = result["curves"]["idvd"]
    lines = [
        f"# TransportModels {result['branch'].upper()} bias-regime analysis",
        "",
        f"- workflow status: `{result['workflow_status']}`",
        f"- comparison trend status: `{result['comparison_status']}`",
        "",
        "| Curve / region | Points | Median error (dex) | Max error (dex) | Max relative error |",
        "|---|---:|---:|---:|---:|",
    ]
    entries = [
        ("Id-Vg off", idvg["regions"]["off"]),
        ("Id-Vg transition", idvg["regions"]["transition"]),
        ("Id-Vg on", idvg["regions"]["on"]),
        ("Id-Vd nonzero bias", idvd["nonzero_bias"]),
    ]
    for label, metrics in entries:
        lines.append(
            f"| {label} | {metrics['points']} | "
            f"{metrics['median_absolute_log10_error_dex']:.6g} | "
            f"{metrics['max_absolute_log10_error_dex']:.6g} | "
            f"{metrics['max_relative_error']:.6g} |")
    lines.extend([
        "",
        "## Endpoints",
        "",
        "| Curve | Bias (V) | Sentaurus (A/um) | Vela (A/um) | Relative error |",
        "|---|---:|---:|---:|---:|",
    ])
    for label, endpoint in (("Id-Vg", idvg["endpoint"]),
                            ("Id-Vd", idvd["endpoint"])):
        lines.append(
            f"| {label} | {endpoint['bias_V']:.6g} | "
            f"{endpoint['reference_A_per_um']:.12g} | "
            f"{endpoint['vela_A_per_um']:.12g} | "
            f"{endpoint['relative_error']:.6g} |")
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workflow-dir", type=Path, required=True)
    parser.add_argument("--generated-dir", type=Path, required=True)
    parser.add_argument("--branch", choices=("dd", "dg"), required=True)
    parser.add_argument("--idvg-candidate", type=Path)
    parser.add_argument("--idvd-candidate", type=Path)
    parser.add_argument("--output-json", type=Path, required=True)
    parser.add_argument("--output-md", type=Path, required=True)
    args = parser.parse_args()
    candidate_paths = {
        kind: path.resolve()
        for kind, path in (("idvg", args.idvg_candidate),
                           ("idvd", args.idvd_candidate))
        if path is not None
    }
    result = analyze(args.workflow_dir.resolve(), args.generated_dir.resolve(),
                     args.branch, candidate_paths)
    args.output_json.parent.mkdir(parents=True, exist_ok=True)
    args.output_md.parent.mkdir(parents=True, exist_ok=True)
    args.output_json.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    args.output_md.write_text(markdown(result), encoding="utf-8")
    print(json.dumps({
        "branch": args.branch,
        "workflow_status": result["workflow_status"],
        "output_json": str(args.output_json.resolve()),
        "output_md": str(args.output_md.resolve()),
    }))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
