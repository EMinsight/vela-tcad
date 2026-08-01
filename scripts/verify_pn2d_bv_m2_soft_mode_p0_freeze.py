#!/usr/bin/env python3
"""Verify the frozen inputs for the PN2D M2 soft-mode P0 contract."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--repo-root", type=Path, default=Path(__file__).resolve().parents[1])
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    digest.update(path.read_bytes())
    return digest.hexdigest()


def main() -> int:
    args = parse_args()
    contract = json.loads(args.contract.read_text(encoding="utf-8"))
    rows = []
    for item in contract["frozen_inputs"]:
        path = args.repo_root / item["path"]
        actual = sha256(path) if path.is_file() else None
        rows.append({
            "path": item["path"],
            "expected_sha256": item["sha256"],
            "actual_sha256": actual,
            "passed": actual == item["sha256"],
        })
    passed = all(row["passed"] for row in rows)
    result = {
        "schema": "vela.pn2d_bv_m2_soft_mode_p0_freeze_verification.v1",
        "status": "passed" if passed else "failed",
        "contract": str(args.contract),
        "files": rows,
    }
    print(json.dumps(result, sort_keys=True))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
