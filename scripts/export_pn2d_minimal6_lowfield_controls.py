#!/usr/bin/env python3
"""Export Minimal6 HighFieldSaturation-off native element controls."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
from typing import Sequence

from export_pn2d_minimal6_states import (
    DEFAULT_IMPORTER,
    DEFAULT_REMOTE_ROOT,
    _live_executor,
    _parse_csv_values,
    default_windows_openssh,
    prepare_exports,
    run_exports,
    write_manifest,
)


ELEMENT_PLOT_LINES = (
    "  Potential/Element",
    "  eDensity/Element",
    "  hDensity/Element",
    "  ElectricField/Element/Vector",
    "  eGradQuasiFermi/Element/Vector",
    "  hGradQuasiFermi/Element/Vector",
    "  eMobility/Element",
    "  hMobility/Element",
    "  eCurrentDensity/Element/Vector",
    "  hCurrentDensity/Element/Vector",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def transform_deck(path: Path) -> str:
    text = path.read_text(encoding="utf-8")
    high_field_line = "    HighFieldSaturation\n"
    if text.count(high_field_line) != 1:
        raise ValueError(
            f"{path} must contain exactly one HighFieldSaturation line"
        )
    text = text.replace(high_field_line, "")
    marker = "  TotalCurrent\n"
    if text.count(marker) != 1:
        raise ValueError(f"{path} lacks a unique TotalCurrent plot marker")
    if any(line in text for line in ELEMENT_PLOT_LINES):
        raise ValueError(f"{path} already contains native element plot fields")
    inserted = marker + "\n".join(ELEMENT_PLOT_LINES) + "\n"
    text = text.replace(marker, inserted)
    path.write_text(text, encoding="utf-8", newline="\n")
    return sha256(path)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topologies", default="mirror,sketch")
    parser.add_argument("--biases", default="-1,-10")
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--ssh-target", default="sentaurus")
    parser.add_argument("--ssh-bin", default=None)
    parser.add_argument("--scp-bin", default=None)
    parser.add_argument("--remote-root", default=DEFAULT_REMOTE_ROOT)
    parser.add_argument("--importer", type=Path, default=DEFAULT_IMPORTER)
    parser.add_argument("--dry-run", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    topologies = tuple(
        value.strip()
        for value in args.topologies.split(",")
        if value.strip()
    )
    biases = _parse_csv_values(args.biases)
    manifest = prepare_exports(
        topology_ids=topologies,
        biases=biases,
        run_id=args.run_id,
        output_dir=args.output_dir,
        ssh_target=args.ssh_target,
        remote_root=args.remote_root,
        importer=args.importer,
    )
    deck_hashes: dict[str, str] = {}
    for state in manifest["states"]:
        deck = (
            Path(str(state["bundle_dir"]))
            / str(state["deck_name"])
        )
        key = (
            f"{state['topology_id']}/"
            f"{state['bias_tag']}/{state['deck_name']}"
        )
        deck_hashes[key] = transform_deck(deck)
    manifest["experiment"] = (
        "pn2d_minimal6_native_lowfield_bias_invariance_control"
    )
    manifest["control_contract"] = {
        "disabled": ["HighFieldSaturation"],
        "retained": [
            "DopingDependence",
            "SRH",
            "Avalanche(VanOverstraeten)",
            "EffectiveIntrinsicDensity(OldSlotboom)",
        ],
        "native_element_fields": list(ELEMENT_PLOT_LINES),
        "temperature_K": 300.0,
    }
    manifest["control_deck_sha256"] = deck_hashes
    write_manifest(Path(str(manifest["manifest_path"])), manifest)

    if not args.dry_run:
        ssh_bin = args.ssh_bin or default_windows_openssh("ssh")
        scp_bin = args.scp_bin or default_windows_openssh("scp")
        importer = args.importer.resolve()
        run_exports(
            manifest,
            executor=lambda state: _live_executor(
                state,
                ssh_bin=ssh_bin,
                scp_bin=scp_bin,
                ssh_target=args.ssh_target,
                importer=importer,
            ),
        )
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    os.environ.setdefault("PYTHONUTF8", "1")
    raise SystemExit(main())
