#!/usr/bin/env python3
"""Prepare reproducible Sentaurus Slot-LDMOS IALMob-off control decks."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = (
    REPO
    / "build-release/reference_tcad/slot_ldmos_sentaurus2022/run01/remote_payload"
)
DEFAULT_OUTPUT = (
    REPO
    / "build-release/reference_tcad/slot_ldmos_sentaurus2022/run01"
    / "ialmob_ablation/run01/payload"
)
CONTROL_DECKS = (
    "avalanche_off_60v",
    "bvds_external_resistor_final",
)
IALMOB_LINE = re.compile(
    r"^[ \t]*Enormal\s*\(\s*IALMob\s*\)[ \t]*\r?\n", re.MULTILINE
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_no_ialmob_deck(template: str, source_stem: str) -> str:
    """Remove exactly one IALMob selector and isolate every output prefix."""
    output_stem = f"{source_stem}_no_ialmob"
    text, count = IALMOB_LINE.subn("", template, count=1)
    if count != 1:
        raise ValueError(
            f"{source_stem}: expected exactly one Enormal(IALMob) line, found {count}"
        )
    text = text.replace(source_stem, output_stem)
    if "IALMob" in text or "Enormal" in text:
        raise ValueError(f"{source_stem}: IALMob/Enormal remained in control deck")
    if output_stem not in text:
        raise ValueError(f"{source_stem}: failed to isolate output names")
    return text


def prepare(source_dir: Path, output_dir: Path) -> dict[str, object]:
    source_dir = source_dir.resolve()
    output_dir = output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    cases: list[dict[str, object]] = []
    for stem in CONTROL_DECKS:
        source = source_dir / f"{stem}.cmd"
        if not source.is_file():
            raise FileNotFoundError(source)
        target = output_dir / f"{stem}_no_ialmob.cmd"
        target.write_text(
            build_no_ialmob_deck(source.read_text(encoding="utf-8"), stem),
            encoding="utf-8",
        )
        cases.append(
            {
                "case": f"{stem}_no_ialmob",
                "source_deck": str(source),
                "source_sha256": sha256(source),
                "control_deck": target.name,
                "control_sha256": sha256(target),
                "physics_delta": "remove Enormal(IALMob) only",
            }
        )

    parameter_source = source_dir / "pp2_des.par"
    if not parameter_source.is_file():
        raise FileNotFoundError(parameter_source)
    parameter_target = output_dir / parameter_source.name
    shutil.copy2(parameter_source, parameter_target)

    manifest: dict[str, object] = {
        "schema": "vela.slot_ldmos.sentaurus_ialmob_ablation.v1",
        "source_dir": str(source_dir),
        "sentaurus_release": "T-2022.03-SP2",
        "shared_grid": "n1_fps.tdr",
        "shared_parameters": parameter_target.name,
        "shared_parameters_sha256": sha256(parameter_target),
        "cases": cases,
    }
    (output_dir / "ialmob_ablation_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-dir", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    manifest = prepare(args.source_dir, args.output_dir)
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
