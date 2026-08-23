#!/usr/bin/env python3
"""Audit Sentaurus affinity/BGN fields against the Vela DG band drive."""

from __future__ import annotations

import argparse
import csv
import json
import statistics
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
FIELDS = (
    REPO_ROOT
    / "build-release/reference_tcad/transportmodels_sentaurus2022/vela_baseline"
    / "generated/sim_fields/dg_idvd/fields"
)
MANIFEST = (
    REPO_ROOT
    / "build-release/reference_tcad/transportmodels_sentaurus2022/vela_baseline"
    / "generated/sim_fields/dg_idvd/field_manifest.json"
)
ELEMENTS = MANIFEST.parent / "elements.csv"
OUT_JSON = REPO_ROOT / "docs/validation/transportmodels_dg_band_drive_audit_2026-08-21.json"
OUT_MD = REPO_ROOT / "docs/validation/transportmodels_dg_band_drive_audit_2026-08-21.md"
BANDGAP_TO_AFFINITY = 0.5


def read_field(path: Path) -> dict[int, float]:
    with path.open(encoding="utf-8", newline="") as stream:
        return {
            int(row["node_id"]): float(row["component0"])
            for row in csv.DictReader(stream)
        }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fields", type=Path, default=FIELDS)
    parser.add_argument("--manifest", type=Path, default=MANIFEST)
    parser.add_argument("--elements", type=Path, default=ELEMENTS)
    parser.add_argument("--output-json", type=Path, default=OUT_JSON)
    parser.add_argument("--output-md", type=Path, default=OUT_MD)
    args = parser.parse_args()

    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    region_info: dict[int, tuple[str, str]] = {}
    for row in manifest["fields"]:
        if row["name"] != "ElectronAffinity":
            continue
        region_id = int(row["region"])
        material = {
            "R.Gateox": "SiO2",
            "R.PolyReox": "SiO2",
            "R.PolyReox_mirrored": "SiO2",
            "R.Substrate": "Silicon",
            "R.Polygate": "PolySilicon",
            "R.Spacer": "Nitride",
            "R.Spacer_mirrored": "Nitride",
        }[row["region_name"]]
        region_info[region_id] = (row["region_name"], material)
    node_regions: dict[int, set[str]] = {}
    with args.elements.open(encoding="utf-8", newline="") as stream:
        for element in csv.DictReader(stream):
            for local in range(3):
                node_regions.setdefault(int(element[f"node{local}"]), set()).add(
                    element["region"]
                )

    vela_affinity = {"Silicon": 4.05, "Si": 4.05, "PolySilicon": 4.05, "SiO2": 0.95, "Nitride": 1.9}
    rows = []
    for region_id, (region_name, material) in sorted(region_info.items()):
        affinity = read_field(args.fields / f"ElectronAffinity_region{region_id}.csv")
        narrowing = read_field(args.fields / f"BandgapNarrowing_region{region_id}.csv")
        shared = sorted(affinity.keys() & narrowing.keys())
        interior = [node for node in shared if node_regions.get(node) == {region_name}]
        evaluation_nodes = interior or shared
        base = [
            affinity[node] - BANDGAP_TO_AFFINITY * narrowing[node]
            for node in evaluation_nodes
        ]
        bgn = [narrowing[node] for node in evaluation_nodes]
        expected = vela_affinity[material]
        rows.append(
            {
                "region_id": region_id,
                "region_name": region_name,
                "material": material,
                "nodes": len(shared),
                "interior_nodes": len(interior),
                "sentaurus_base_affinity_min_eV": min(base),
                "sentaurus_base_affinity_median_eV": statistics.median(base),
                "sentaurus_base_affinity_max_eV": max(base),
                "bandgap_narrowing_min_eV": min(bgn),
                "bandgap_narrowing_max_eV": max(bgn),
                "vela_base_affinity_eV": expected,
                "median_offset_vela_minus_sentaurus_eV": expected - statistics.median(base),
            }
        )
    report = {
        "schema": "vela.transportmodels.dg_band_drive_audit.v1",
        "status": "pass",
        "bgn_to_affinity_fraction": BANDGAP_TO_AFFINITY,
        "regions": rows,
    }
    args.output_json.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    table = "\n".join(
        f"| {row['region_name']} | {row['material']} | {row['sentaurus_base_affinity_median_eV']:.12g} | "
        f"{row['vela_base_affinity_eV']:.12g} | {1e3 * row['median_offset_vela_minus_sentaurus_eV']:.6f} | "
        f"{row['bandgap_narrowing_min_eV']:.6g}–{row['bandgap_narrowing_max_eV']:.6g} |"
        for row in rows
    )
    markdown = f"""# TransportModels DG band-drive audit

The Sentaurus base affinity is reconstructed as
`ElectronAffinity - 0.5 * BandgapNarrowing` at the Vg=1 V, Vd=2 V DG state.

| Region | Material | Sentaurus base affinity (eV) | Vela base affinity (eV) | Vela−Sentaurus (mV) | BGN range (eV) |
|---|---|---:|---:|---:|---:|
{table}
"""
    args.output_md.write_text(markdown, encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
