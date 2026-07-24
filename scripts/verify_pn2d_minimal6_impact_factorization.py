#!/usr/bin/env python3
"""Independent arithmetic verifier for Minimal6 impact factorization."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path


OUTPUTS = (
    "local_edge_factorization.csv",
    "state_source_factorization.csv",
    "report.md",
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def close(actual: float, expected: float, tolerance: float = 1.0e-12) -> None:
    scale = max(abs(actual), abs(expected), 1.0e-300)
    if abs(actual - expected) / scale > tolerance:
        raise ValueError(
            f"arithmetic closure failed: actual={actual}, expected={expected}"
        )


def verify(root: Path) -> dict:
    root = root.resolve()
    manifest = json.loads((root / "manifest.json").read_text(encoding="utf-8"))
    if manifest.get("experiment") != "pn2d_minimal6_impact_factorization":
        raise ValueError("unexpected impact experiment")
    if manifest.get("status") != "valid":
        raise ValueError("impact factorization is not valid")
    for name in OUTPUTS:
        if sha256(root / name) != manifest["outputs"][name]:
            raise ValueError(f"hash mismatch for {name}")

    local = rows(root / "local_edge_factorization.csv")
    state = rows(root / "state_source_factorization.csv")
    if len(local) != 960 or len(state) != 40:
        raise ValueError("impact row-count contract failed")

    totals: dict[tuple[str, int], dict[str, float]] = defaultdict(
        lambda: defaultdict(float)
    )
    zero_count = 0
    for row in local:
        volume_m2 = float(row["volume_m2"])
        volume_cm2 = volume_m2 * 1.0e4
        alpha_cm_inv = float(row["candidate_alpha_per_m"]) / 100.0
        sent_alpha_cm_inv = (
            float(row["sentaurus_endpoint_alpha_per_m"]) / 100.0
        )
        candidate_flux = float(row["candidate_flux_proxy_per_cm2_s"])
        sent_flux = float(row["sentaurus_endpoint_flux_per_cm2_s"])
        candidate = float(row["candidate_source_per_cm_s"])
        sent_alpha = float(row["sentaurus_alpha_hybrid_source_per_cm_s"])
        sent_current = float(
            row["sentaurus_current_hybrid_source_per_cm_s"]
        )
        projected = float(row["projected_sentaurus_source_per_cm_s"])
        close(candidate, alpha_cm_inv * candidate_flux * volume_cm2)
        close(sent_alpha, sent_alpha_cm_inv * candidate_flux * volume_cm2)
        close(sent_current, alpha_cm_inv * sent_flux * volume_cm2)
        if int(row["geometric_zero"]):
            zero_count += 1
            if any(
                value != 0.0
                for value in (candidate, sent_alpha, sent_current, projected)
            ):
                raise ValueError("geometric zero has nonzero source")
        key = (row["topology"], int(float(row["bias_V"])))
        totals[key]["candidate"] += candidate
        totals[key]["sent_alpha"] += sent_alpha
        totals[key]["sent_current"] += sent_current
        totals[key]["projected"] += projected

    if zero_count != 320:
        raise ValueError(
            f"expected 320 carrier-local geometric zeros, found {zero_count}"
        )
    if manifest["contracts"]["zero_volume_count"] != 160:
        raise ValueError("manifest local-edge geometric-zero count must be 160")
    if manifest["contracts"]["zero_volume_nonzero_source_count"] != 0:
        raise ValueError("manifest reports a nonzero geometric-zero source")
    if manifest["contracts"]["maximum_geometry_relative_error"] > 1.0e-15:
        raise ValueError("triangle local-volume geometry does not close")
    if (
        manifest["contracts"]["maximum_node_mapping_relative_error"]
        > 1.0e-15
    ):
        raise ValueError("source-to-node mapping does not close")

    for row in state:
        key = (row["topology"], int(float(row["bias_V"])))
        aggregate = totals[key]
        close(
            aggregate["candidate"],
            float(row["vela_candidate_source_per_cm_s"]),
        )
        close(
            aggregate["sent_alpha"],
            float(row["sentaurus_alpha_hybrid_source_per_cm_s"]),
        )
        close(
            aggregate["sent_current"],
            float(row["sentaurus_current_hybrid_source_per_cm_s"]),
        )
        close(
            aggregate["projected"],
            float(row["sentaurus_projected_triangle_source_per_cm_s"]),
        )
        if float(row["node_mapping_relative_error"]) > 1.0e-15:
            raise ValueError("state node mapping exceeds tolerance")

    if manifest["summaries"]["nodal_reconstruction"]["p95"] > 1.0e-3:
        raise ValueError("Sentaurus alpha-current reconstruction does not close")
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("--compare-root", type=Path)
    args = parser.parse_args()
    first = verify(args.root)
    deterministic = None
    if args.compare_root is not None:
        second = verify(args.compare_root)
        for name in OUTPUTS:
            if (args.root / name).read_bytes() != (
                args.compare_root / name
            ).read_bytes():
                raise ValueError(f"A/B output mismatch for {name}")
        if first["contracts"] != second["contracts"]:
            raise ValueError("A/B impact contracts differ")
        if first["summaries"] != second["summaries"]:
            raise ValueError("A/B impact summaries differ")
        deterministic = True
    print(json.dumps(
        {
            "status": "passed",
            "outcome": first["outcome"],
            "deterministic_pair": deterministic,
            "verified_output_count": len(OUTPUTS),
        },
        sort_keys=True,
    ))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
