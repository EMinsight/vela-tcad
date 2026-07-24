#!/usr/bin/env python3
"""Run the deterministic Minimal6 mobility unit root-cause audit."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.pn2d_minimal6_diagnostics.mobility_unit_root_cause import (
    run_mobility_unit_root_cause,
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--self-consistent-root", required=True, type=Path)
    parser.add_argument("--inverse-inputs-root", required=True, type=Path)
    parser.add_argument("--sentaurus-element-csv", required=True, type=Path)
    parser.add_argument("--output-root", required=True, type=Path)
    args = parser.parse_args()
    result = run_mobility_unit_root_cause(
        self_consistent_root=args.self_consistent_root,
        inverse_inputs_root=args.inverse_inputs_root,
        sentaurus_element_csv=args.sentaurus_element_csv,
        output_root=args.output_root,
    )
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
