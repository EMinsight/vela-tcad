#!/usr/bin/env python3
"""Run the sealed PN2D Minimal6 Phase D mobility/support diagnostic."""

from __future__ import annotations

import argparse
import json

from pn2d_minimal6_diagnostics.phase_d_mobility_support import run_phase_d


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--models-par", required=True)
    parser.add_argument("--sdevice-cmd", required=True)
    parser.add_argument("--vela-deck", required=True)
    parser.add_argument("--mobility-source", required=True)
    parser.add_argument("--inverse-inputs-root", required=True)
    parser.add_argument("--mapped-transport-csv", required=True)
    parser.add_argument("--cell-mapping-csv", required=True)
    parser.add_argument("--geometry-csv", required=True)
    parser.add_argument("--local-mobility-csv", required=True)
    parser.add_argument("--global-mobility-csv", required=True)
    parser.add_argument("--stage-edge-csv", required=True)
    parser.add_argument("--output-root", required=True)
    args = parser.parse_args()
    manifest = run_phase_d(
        models_par=args.models_par,
        sdevice_cmd=args.sdevice_cmd,
        vela_deck=args.vela_deck,
        mobility_source=args.mobility_source,
        inverse_inputs_root=args.inverse_inputs_root,
        mapped_transport_csv=args.mapped_transport_csv,
        cell_mapping_csv=args.cell_mapping_csv,
        geometry_csv=args.geometry_csv,
        local_mobility_csv=args.local_mobility_csv,
        global_mobility_csv=args.global_mobility_csv,
        stage_edge_csv=args.stage_edge_csv,
        output_root=args.output_root,
    )
    print(json.dumps(manifest["outcome"], sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

