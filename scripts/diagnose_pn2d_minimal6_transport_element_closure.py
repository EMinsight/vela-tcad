#!/usr/bin/env python3
"""Run the native-element Minimal6 transport closure audit."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.pn2d_minimal6_diagnostics.transport_element_closure import (
    run_transport_element_closure,
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--transport-csv", required=True)
    parser.add_argument("--transport-manifest", required=True)
    parser.add_argument("--observations", required=True)
    parser.add_argument("--inverse-inputs-root", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    manifest = run_transport_element_closure(
        transport_csv=args.transport_csv,
        transport_manifest=args.transport_manifest,
        observations_csv=args.observations,
        inverse_inputs_root=args.inverse_inputs_root,
        output_root=args.output,
    )
    print(
        json.dumps(
            {
                "status": manifest["status"],
                "state_count": manifest["state_count"],
                "sample_count": manifest["sample_count"],
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
