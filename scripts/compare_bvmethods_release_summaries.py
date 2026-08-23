#!/usr/bin/env python3
"""Compare two outputs from analyze_sentaurus_bvmethods.py by method."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Sequence


def compare(baseline: dict, candidate: dict, max_relative_drift: float) -> dict:
    baseline_by_method = {row["method"]: row for row in baseline["results"]}
    candidate_by_method = {row["method"]: row for row in candidate["results"]}
    methods = sorted(set(baseline_by_method) | set(candidate_by_method))
    rows = []
    for method in methods:
        old = baseline_by_method.get(method)
        new = candidate_by_method.get(method)
        if old is None or new is None:
            rows.append({
                "method": method,
                "status": "missing",
                "baseline_present": old is not None,
                "candidate_present": new is not None,
            })
            continue
        old_bv = float(old["bv_V"])
        new_bv = float(new["bv_V"])
        absolute = new_bv - old_bv
        relative = abs(absolute) / max(abs(old_bv), 1.0e-300)
        rows.append({
            "method": method,
            "status": "pass" if relative <= max_relative_drift else "fail",
            "baseline_bv_V": old_bv,
            "candidate_bv_V": new_bv,
            "signed_bv_drift_V": absolute,
            "relative_bv_drift": relative,
            "baseline_rows": int(old["rows"]),
            "candidate_rows": int(new["rows"]),
            "criterion_match": old["criterion"] == new["criterion"],
        })
    comparable = [row for row in rows if "relative_bv_drift" in row]
    passed = (
        len(comparable) == len(methods)
        and all(row["status"] == "pass" and row["criterion_match"] for row in rows)
    )
    return {
        "schema": "vela.sentaurus_bvmethods_release_compare.v1",
        "status": "pass" if passed else "fail",
        "max_relative_bv_drift_gate": max_relative_drift,
        "maximum_relative_bv_drift": max(
            (row["relative_bv_drift"] for row in comparable), default=None),
        "methods": rows,
    }


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline", type=Path, required=True)
    parser.add_argument("--candidate", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--max-relative-bv-drift", type=float, default=0.01)
    args = parser.parse_args(argv)
    if args.max_relative_bv_drift < 0.0:
        parser.error("--max-relative-bv-drift must be non-negative")

    result = compare(
        json.loads(args.baseline.read_text()),
        json.loads(args.candidate.read_text()),
        args.max_relative_bv_drift,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
    return 0 if result["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
