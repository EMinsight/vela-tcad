#!/usr/bin/env python3
"""Generate the PN2D Minimal6 Phase F self-consistent comparison."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.pn2d_minimal6_diagnostics.phase_f_self_consistent import run_phase_f


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-sweep-root", type=Path, required=True)
    parser.add_argument("--sentaurus-sweep-root", type=Path, required=True)
    parser.add_argument("--inverse-inputs-root", type=Path, required=True)
    parser.add_argument("--phase-c-root", type=Path, required=True)
    parser.add_argument("--phase-d-root", type=Path, required=True)
    parser.add_argument("--operator-audit", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()
    manifest = run_phase_f(
        candidate_sweep_root=args.candidate_sweep_root,
        sentaurus_sweep_root=args.sentaurus_sweep_root,
        inverse_inputs_root=args.inverse_inputs_root,
        phase_c_root=args.phase_c_root,
        phase_d_root=args.phase_d_root,
        operator_audit=args.operator_audit,
        output_root=args.output_root,
    )
    print(json.dumps(manifest["outcome"], sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
