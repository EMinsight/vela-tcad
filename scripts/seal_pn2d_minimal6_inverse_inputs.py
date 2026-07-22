#!/usr/bin/env python3
"""Seal live Minimal6 sweep/export roots into canonical inverse-audit inputs."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys


REPO = Path(__file__).resolve().parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from scripts.pn2d_minimal6_diagnostics.inverse_input_sealer import (  # noqa: E402
    seal_inverse_input_roots,
)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create deterministic sealed canonical roots for the Minimal6 inverse audit."
    )
    parser.add_argument("--vela-sweep-root", type=Path, required=True)
    parser.add_argument("--sentaurus-sweep-root", type=Path, required=True)
    parser.add_argument("--supplemental-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--importer", type=Path, required=True)
    parser.add_argument("--vela-executable", type=Path, required=True)
    parser.add_argument("--phase-base", required=True)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    roots = seal_inverse_input_roots(
        args.vela_sweep_root, args.sentaurus_sweep_root,
        args.supplemental_root, args.output_root,
        importer=args.importer, vela_executable=args.vela_executable,
        phase_base=args.phase_base,
    )
    for name in ("vela", "sentaurus", "supplemental"):
        print(f"{name}_root={roots[name]}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
