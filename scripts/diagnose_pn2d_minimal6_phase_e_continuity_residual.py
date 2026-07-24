#!/usr/bin/env python3
"""Run the PN2D Minimal6 Phase E continuity-residual audit."""

from __future__ import annotations

import argparse
import json

from pn2d_minimal6_diagnostics.phase_e_continuity_residual import run_phase_e


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--runner", required=True)
    parser.add_argument("--inverse-inputs-root", required=True)
    parser.add_argument("--self-consistent-root", required=True)
    parser.add_argument("--phase-d-root", required=True)
    parser.add_argument("--output-root", required=True)
    args = parser.parse_args()
    manifest = run_phase_e(
        runner=args.runner,
        inverse_inputs_root=args.inverse_inputs_root,
        self_consistent_root=args.self_consistent_root,
        phase_d_root=args.phase_d_root,
        output_root=args.output_root,
    )
    print(json.dumps(manifest["outcome"], sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

