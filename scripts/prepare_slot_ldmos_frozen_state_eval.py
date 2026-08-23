#!/usr/bin/env python3
"""Prepare a no-Newton terminal-current evaluation of a saved DD state."""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle", type=Path, required=True)
    parser.add_argument("--source-config", required=True)
    parser.add_argument("--state", required=True)
    parser.add_argument("--bias", type=float, required=True)
    parser.add_argument("--name", required=True)
    args = parser.parse_args()
    bundle = args.bundle.resolve()
    deck = json.loads((bundle / args.source_config).read_text(encoding="utf-8"))
    output = f"outputs/ialmob_ablation/direct_bordered_20260822_v5/frozen_eval/{args.name}"
    (bundle / output).mkdir(parents=True, exist_ok=True)
    deck["simulation_type"] = "newton_solve_from_state"
    deck["state_file"] = args.state
    deck["solver"]["method"] = "newton"
    deck["solver"]["max_iter"] = 0
    for contact in deck["contacts"]:
        if contact["name"] == "drain":
            contact["bias"] = args.bias
    deck.pop("sweep", None)
    deck.pop("output_csv", None)
    config_name = f"simulation_frozen_eval_{args.name}.json"
    (bundle / config_name).write_text(
        json.dumps(deck, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(config_name)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
