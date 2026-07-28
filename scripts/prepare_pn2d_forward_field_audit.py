#!/usr/bin/env python3
"""Prepare matched coarse7x3 Sentaurus and Vela forward-field audit inputs."""

from __future__ import annotations

import json
import shutil
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "build-release" / "pn2d-forward-field-audit-20260727"
SENTAURUS_SOURCE = ROOT / "reference_tcad" / "pn2d_sentaurus2018" / "source"
COARSE_SOURCE = (
    ROOT
    / "build-release"
    / "reference_tcad"
    / "pn2d_sentaurus2018_coarse7x3"
    / "sentaurus_vm_runs"
    / "coarse7x3_vector_bv_20260627"
    / "source"
)
VELA_SOURCE = (
    ROOT
    / "build-release"
    / "pn2d-forward-iv-0v20v-20260727"
    / "simulation_vela_coarse7x3_forward_iv_legacy_lowfield_0v20v.json"
)


def prepare_sentaurus() -> Path:
    bundle = OUT / "bundle"
    bundle.mkdir(parents=True, exist_ok=True)
    for name in ("pn2d_msh.tdr", "models.par"):
        shutil.copy2(COARSE_SOURCE / name, bundle / name)

    text = (SENTAURUS_SOURCE / "pn2d_iv_sdevice.cmd").read_text(encoding="ascii")
    text = text.replace('Voltage=10.0', 'Voltage=20.0', 1)
    old_plot = (
        '    # 0.05 V spacing over the 0-10 V normalized sweep for per-bias TDR comparison.\n'
        '    Plot(FilePrefix="pn2d_iv_multibias" Time=(Range=(0 1) Intervals=200) NoOverWrite)'
    )
    new_plot = (
        '    # Exact field snapshots at 0, 1, 2, 5, 10, 15, and 20 V.\n'
        '    Plot(FilePrefix="pn2d_forward_fields" '
        'Time=(0;0.05;0.1;0.25;0.5;0.75;1.0) NoOverWrite)'
    )
    if old_plot not in text:
        raise RuntimeError("Sentaurus multibias Plot anchor was not found")
    text = text.replace(old_plot, new_plot, 1)
    deck = bundle / "pn2d_forward_fields_sdevice.cmd"
    deck.write_text(text, encoding="ascii", newline="\n")
    return deck


def prepare_vela() -> Path:
    config = json.loads(VELA_SOURCE.read_text(encoding="utf-8-sig"))
    config["output_csv"] = str((OUT / "vela_forward_fields_0v20v.csv").resolve())
    sweep = config["sweep"]
    sweep["write_vtk"] = True
    sweep["vtk_prefix"] = str((OUT / "vela_forward_fields").resolve())
    sweep["write_state_file"] = str((OUT / "vela_forward_fields_last_state.csv").resolve())
    diagnostics = sweep.get("diagnostics", {})
    if "newton_history" in diagnostics:
        diagnostics["newton_history"]["csv_file"] = str(
            (OUT / "vela_forward_fields_newton_history.csv").resolve()
        )
    deck = OUT / "simulation_vela_forward_fields_0v20v.json"
    deck.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    return deck


def main() -> None:
    sentaurus = prepare_sentaurus()
    vela = prepare_vela()
    print(json.dumps({"sentaurus": str(sentaurus), "vela": str(vela)}, indent=2))


if __name__ == "__main__":
    main()
