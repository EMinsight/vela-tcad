#!/usr/bin/env python3
"""Run the paired Task 14 Sentaurus high-bias oracle variant matrix."""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.run_pn2d_high_bias_oracle_variant_vm import ORACLE_VARIANTS


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-prefix",
        type=Path,
        default=Path("build-release/pn2d-task14-high-bias-oracle-20260726"),
    )
    parser.add_argument(
        "--remote-prefix",
        default="/home/tcad/codex/pn2d-task14-high-bias-oracle-20260726",
    )
    return parser.parse_args()


def passed(path: Path) -> bool:
    manifest = path / "manifest.json"
    if not manifest.is_file():
        return False
    return json.loads(manifest.read_text(encoding="ascii")).get("status") == "passed"


def main() -> int:
    args = parse_args()
    progress_path = args.output_prefix.parent / (
        args.output_prefix.name + "-matrix-progress.json"
    )
    completed: list[dict[str, str]] = []
    for root in ("a", "b"):
        output_root = Path(str(args.output_prefix) + f"-{root}")
        remote_root = args.remote_prefix + f"-{root}"
        for variant in ORACLE_VARIANTS:
            destination = output_root / variant
            if not passed(destination):
                subprocess.run(
                    [
                        sys.executable,
                        "scripts/run_pn2d_high_bias_oracle_variant_vm.py",
                        "--variant",
                        variant,
                        "--output-root",
                        str(destination),
                        "--remote-root",
                        remote_root,
                    ],
                    check=True,
                )
            completed.append({"root": root, "variant": variant, "status": "passed"})
            progress_path.write_text(
                json.dumps(
                    {
                        "completed": completed,
                        "completed_count": len(completed),
                        "total_count": 2 * len(ORACLE_VARIANTS),
                    },
                    indent=2,
                    sort_keys=True,
                )
                + "\n",
                encoding="ascii",
                newline="\n",
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
