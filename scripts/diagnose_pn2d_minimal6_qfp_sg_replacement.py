#!/usr/bin/env python3
"""CLI for the Minimal6 internal-node QFP SG replacement experiment."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.pn2d_minimal6_diagnostics.qfp_sg_experiment import (
    run_qfp_replacement_experiment,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--observations", type=Path, required=True)
    parser.add_argument("--current-edges", type=Path, required=True)
    parser.add_argument("--inverse-inputs-root", type=Path, required=True)
    parser.add_argument("--operator-audit", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    manifest = run_qfp_replacement_experiment(
        observations_csv=args.observations,
        current_edges_csv=args.current_edges,
        inverse_inputs_root=args.inverse_inputs_root,
        operator_audit_executable=args.operator_audit,
        output_root=args.output,
    )
    print(json.dumps(manifest["baseline_cpp_replay"], sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
