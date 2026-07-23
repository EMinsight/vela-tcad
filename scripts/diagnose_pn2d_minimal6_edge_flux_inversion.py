#!/usr/bin/env python3
"""CLI for the Minimal6 directed-edge SG inversion audit."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.pn2d_minimal6_diagnostics.edge_flux_experiment import (
    run_edge_flux_experiment,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--observations", type=Path, required=True)
    parser.add_argument("--effective-mobility", type=Path, required=True)
    parser.add_argument("--inverse-inputs-root", type=Path, required=True)
    parser.add_argument("--sentaurus-export-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    manifest = run_edge_flux_experiment(
        observations_csv=args.observations,
        effective_mobility_csv=args.effective_mobility,
        inverse_inputs_root=args.inverse_inputs_root,
        sentaurus_export_root=args.sentaurus_export_root,
        output_root=args.output,
    )
    print(
        json.dumps(
            {
                "status": manifest["status"],
                "state_count": manifest["state_count"],
                "sg_replacement_sample_count": manifest[
                    "sg_replacement_sample_count"
                ],
                "formula_change_authorized": manifest["acceptance_policy"][
                    "formula_change_authorized"
                ],
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
