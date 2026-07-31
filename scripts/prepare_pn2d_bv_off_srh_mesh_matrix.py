#!/usr/bin/env python3
"""Prepare sealed Sentaurus sources for the PN2D SRH mesh matrix."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
from pathlib import Path


LEVELS = {
    "M0": {"junction_x": 1.0 / 3.0, "junction_y": 0.25},
    "M1": {"junction_x": 1.0 / 6.0, "junction_y": 0.125},
    "M2": {"junction_x": 1.0 / 12.0, "junction_y": 0.0625},
}
# SDE applies a geometric tolerance when evaluating profile-window boundaries.
# 1e-9 um was below that tolerance in O-2018.06-SP2 and still double-counted
# x=XJ nodes.  1e-3 um is well above that tolerance while remaining 62.5x
# smaller than the finest requested junction spacing (0.0625 um).
JUNCTION_EPSILON_UM = 1.0e-3


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def replace_size(text: str, name: str, value: float) -> str:
    pattern = (
        rf'(\(sdedr:define-refinement-size\s+"{re.escape(name)}"\s+)'
        rf"[-+0-9.eE]+\s+[-+0-9.eE]+\s+[-+0-9.eE]+\s+[-+0-9.eE]+(\s*\))"
    )
    replacement = rf"\g<1>{value:g} {value:g} {value:g} {value:g}\g<2>"
    updated, count = re.subn(pattern, replacement, text, flags=re.MULTILINE)
    if count != 1:
        raise RuntimeError(f"expected one {name} refinement-size block, found {count}")
    return updated


def add_junction_refinement(text: str, spacing_x: float, spacing_y: float) -> str:
    marker = ";----------------------------------------------------------\n; Build mesh"
    if marker not in text:
        raise RuntimeError("missing build-mesh marker in sealed coarse7x3 source")
    block = f""";----------------------------------------------------------
; Junction-focused Task 4 refinement
;----------------------------------------------------------

(sdedr:define-refeval-window
  "Junction.Window"
  "Rectangle"
  (position 0.75 0.0 0.0)
  (position 1.25 H 0.0)
)

(sdedr:define-refinement-size
  "Junction.Mesh"
  {spacing_x:.17g} {spacing_y:.17g}
  {spacing_x:.17g} {spacing_y:.17g}
)

(sdedr:define-refinement-placement
  "Junction.Mesh.Place"
  "Junction.Mesh"
  "Junction.Window"
)

"""
    return text.replace(marker, block + marker, 1)


def make_dose_preserving_junction(text: str) -> str:
    """Assign the x=XJ nodes to exactly one constant-profile window.

    SDE rectangle windows include their boundaries.  The historical deck used
    XJ as both the P-window end and N-window start, so junction nodes received
    both active-species profiles.  Moving only the P-window end by a
    sub-mesh epsilon keeps the continuum junction at XJ while making nodal
    ownership unambiguous and the total-impurity dose mesh independent.
    """

    define_anchor = "(define XJ 1.0)      ; PN junction position"
    if text.count(define_anchor) != 1:
        raise RuntimeError("expected exactly one XJ definition")
    text = text.replace(
        define_anchor,
        define_anchor
        + f"\n(define XJ_P (- XJ {JUNCTION_EPSILON_UM:.17g}))"
        + " ; exclusive P-window end",
        1,
    )
    p_window_end = "(position XJ H 0.0)"
    if text.count(p_window_end) != 1:
        raise RuntimeError("expected exactly one P-window XJ endpoint")
    return text.replace(p_window_end, "(position XJ_P H 0.0)", 1)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sealed-source", type=Path, required=True)
    parser.add_argument("--out-root", type=Path, required=True)
    args = parser.parse_args()
    source = args.sealed_source.resolve()
    root = args.out_root.resolve()
    root.mkdir(parents=True, exist_ok=True)
    base_sde = make_dose_preserving_junction(
        (source / "pn2d_sde.cmd").read_text(encoding="utf-8")
    )
    rows = []
    for level, sizes in LEVELS.items():
        out = root / level
        out.mkdir(parents=True, exist_ok=True)
        sde = base_sde
        if level != "M0":
            sde = add_junction_refinement(
                sde, sizes["junction_x"], sizes["junction_y"]
            )
        (out / "pn2d_sde.cmd").write_text(sde, encoding="utf-8")
        shutil.copyfile(source / "pn2d_bv_sdevice.cmd", out / "pn2d_bv_sdevice.cmd")
        shutil.copyfile(source / "models.par", out / "models.par")
        rows.append(
            {
                "level": level,
                "global_spacing_x_um": 1.0 / 3.0,
                "global_spacing_y_um": 0.25,
                "junction_spacing_x_um": sizes["junction_x"],
                "junction_spacing_y_um": sizes["junction_y"],
                "relative_junction_refinement_vs_M0": (
                    LEVELS["M0"]["junction_x"] / sizes["junction_x"]
                ),
                "source_dir": str(out),
                "hashes": {
                    name: digest(out / name)
                    for name in ("pn2d_sde.cmd", "pn2d_bv_sdevice.cmd", "models.par")
                },
            }
        )
    manifest = {
        "schema": "vela.pn2d_bv_off_srh_mesh_matrix_sources.v1",
        "geometry_um": {"length": 2.0, "height": 0.5, "junction_x": 1.0},
        "doping_cm3": {"p": 1.0e17, "n": 1.0e17},
        "junction_profile_contract": {
            "mode": "single_owner_submesh_epsilon",
            "junction_owner": "N.Window",
            "p_window_end_um": 1.0 - JUNCTION_EPSILON_UM,
            "n_window_start_um": 1.0,
            "epsilon_um": JUNCTION_EPSILON_UM,
            "expected_nodes_in_gap": 0,
            "double_counted_junction_nodes_allowed": False,
        },
        "levels": rows,
    }
    (root / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
