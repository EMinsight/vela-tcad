#!/usr/bin/env python3
"""Prepare a minimal Sentaurus Device box-measure debug probe bundle."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--source-deck",
        default="pn2d_bv_process_avalanche_off.cmd",
    )
    parser.add_argument(
        "--box-method",
        choices=("average", "mix-average", "weighted-voronoi"),
        default="average",
        help="Sentaurus box discretization selected explicitly in Math.",
    )
    args = parser.parse_args()

    source = args.source_dir.resolve()
    output = args.output_dir.resolve()
    if output.exists():
        raise FileExistsError(f"probe output already exists: {output}")
    output.mkdir(parents=True)

    copied = []
    for name in (
        "models.par",
        "pn2d_msh.tdr",
        "runtime_element_avalanche_probe.tcl",
    ):
        src = source / name
        dst = output / name
        shutil.copy2(src, dst)
        copied.append({"name": name, "sha256": sha256(dst)})

    source_deck = source / args.source_deck
    text = source_deck.read_text(encoding="utf-8")
    math_marker = "Math {"
    solve_marker = "Solve {"
    if text.count(math_marker) != 1 or text.count(solve_marker) != 1:
        raise ValueError("expected exactly one Math and one Solve section")
    method_keyword = {
        "average": "AverageBoxMethod",
        "mix-average": "MixAverageBoxMethod",
        "weighted-voronoi": "WeightedVoronoiBox",
    }[args.box_method]
    text = text.replace(
        math_marker,
        math_marker
        + f"\n  {method_keyword}"
        + "\n  BoxMeasureFromFile(GrdNumbering)",
        1,
    )
    text, plot_count = re.subn(
        r"(?m)^Plot\s*\{",
        "Plot {"
        + "\n  BM_AngleElements"
        + "\n  BM_CoeffIntersectionNonDelaunayElements"
        + "\n  BM_ElementVolume"
        + "\n  BM_IntersectionNonDelaunayElements"
        + "\n  BM_VolumeIntersectionNonDelaunayElements",
        text,
    )
    if plot_count != 1:
        raise ValueError("expected exactly one top-level Plot section")
    text = text[: text.index(solve_marker)] + (
        "Solve {\n"
        "  Coupled(Iterations=1) { Poisson }\n"
        '  Plot(FilePrefix="box_measure_probe")\n'
        "}\n"
    )
    deck = output / "box_measure_probe.cmd"
    deck.write_text(text, encoding="utf-8")

    manifest = {
        "schema": "vela.sentaurus_box_measure_probe_bundle.v1",
        "source_dir": str(source),
        "source_deck": str(source_deck),
        "source_deck_sha256": sha256(source_deck),
        "probe_deck": str(deck),
        "probe_deck_sha256": sha256(deck),
        "math_keyword": "BoxMeasureFromFile(GrdNumbering)",
        "box_method": args.box_method,
        "box_method_keyword": method_keyword,
        "solve_scope": "single_poisson_initialization",
        "plot_file_prefix": "box_measure_probe",
        "copied_inputs": copied,
    }
    (output / "bundle_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
