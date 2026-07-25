#!/usr/bin/env python3
"""Replay the documented Sentaurus box operator on all 40 Minimal6 states."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


from scripts.diagnose_pn2d_minimal6_element_avalanche_replay import run


TARGET_BIASES = tuple(float(-magnitude) for magnitude in range(1, 21))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--raw-root",
        type=Path,
        default=Path(
            "build-release/"
            "pn2d-minimal6-element-avalanche-all20-runtime-20260725"
        ),
    )
    parser.add_argument(
        "--vela-factorization",
        type=Path,
        default=Path(
            "build-release/"
            "pn2d-minimal6-impact-factorization-final-20260724-b/"
            "state_source_factorization.csv"
        ),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(
            "build-release/"
            "pn2d-minimal6-element-avalanche-all20-runtime-20260725/"
            "analysis"
        ),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    manifest = run(
        args.raw_root,
        args.vela_factorization,
        args.output,
        target_biases=TARGET_BIASES,
        log_relative=Path("fetched/run_all20.out"),
        plt_relative=Path(
            "fetched/runtime_element_avalanche_probe_all20.plt"
        ),
        experiment="pn2d_minimal6_element_avalanche_all20_replay",
    )
    if manifest["scope"]["biases_V"] != list(TARGET_BIASES):
        raise RuntimeError("full 20-bias matrix was not preserved")
    print(json.dumps(manifest["rankings"], indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
