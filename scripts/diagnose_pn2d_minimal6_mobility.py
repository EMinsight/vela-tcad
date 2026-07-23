#!/usr/bin/env python3
"""Run the deterministic Minimal6 same-support mobility diagnosis."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.pn2d_minimal6_diagnostics.mobility_diagnosis import (
    build_mobility_diagnosis,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--vela-root", required=True)
    parser.add_argument("--sentaurus-root", required=True)
    parser.add_argument("--supplemental-root", required=True)
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    report = build_mobility_diagnosis(
        args.vela_root, args.sentaurus_root, args.supplemental_root, args.output
    )
    print(json.dumps({
        "schema": report["schema"],
        "state_count": report["state_count"],
        "edge_sample_count": report["edge_sample_count"],
        "cell_sample_count": report["cell_sample_count"],
        "conclusions": report["conclusions"],
    }, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
