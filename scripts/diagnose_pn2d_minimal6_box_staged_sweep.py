#!/usr/bin/env python3
"""CLI for the forty-state Minimal6 Sentaurus box staged replay."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.pn2d_minimal6_diagnostics.box_staged_sweep import (
    run_box_staged_sweep,
)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--observations", type=Path, required=True)
    parser.add_argument("--transport-elements", type=Path, required=True)
    parser.add_argument("--vela-mobility", type=Path, required=True)
    parser.add_argument("--vela-replay-root", type=Path, required=True)
    parser.add_argument("--mesh-root", type=Path, required=True)
    parser.add_argument("--sentaurus-state-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    manifest = run_box_staged_sweep(
        observations_csv=args.observations,
        transport_elements_csv=args.transport_elements,
        vela_mobility_csv=args.vela_mobility,
        vela_replay_root=args.vela_replay_root,
        mesh_root=args.mesh_root,
        sentaurus_state_root=args.sentaurus_state_root,
        output_root=args.output,
    )
    print(json.dumps(manifest["gates"], sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
