#!/usr/bin/env python3
"""Build a full Vela restart state with Sentaurus semiconductor fields overlaid.

Sentaurus exports the carrier state only for the silicon region, whereas Vela's
restart CSV requires every mesh node.  Preserve the baseline values outside the
transport region and replace psi, quasi-Fermi potentials, and carrier densities
on every node present in the Sentaurus region export.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


STATE_COLUMNS = (
    "node_id",
    "psi",
    "phin",
    "phip",
    "electrons_m3",
    "holes_m3",
)


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def scalar_field(root: Path, name: str, region: int) -> dict[int, float]:
    path = root / f"{name}_region{region}.csv"
    values = rows(path)
    if not values or "component0" not in values[0]:
        raise RuntimeError(f"missing scalar field values: {path}")
    return {int(row["node_id"]): float(row["component0"]) for row in values}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-state", type=Path, required=True)
    parser.add_argument("--sentaurus-fields", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--audit-output", type=Path)
    parser.add_argument("--baseline-audit-output", type=Path)
    parser.add_argument("--simulation-config", type=Path)
    parser.add_argument("--audit-config-output", type=Path)
    parser.add_argument("--probe-output-dir", type=Path)
    parser.add_argument("--bias", type=float)
    parser.add_argument("--region", type=int, default=3)
    parser.add_argument(
        "--preserve-eparallel",
        action="store_true",
        help="Keep the production eparallel avalanche drive in generated probes.",
    )
    parser.add_argument(
        "--btbt-source-integration",
        choices=("semiconductor_cell_lumped", "transport_node_lumped"),
        help="Override the E2 spatial recovery used by generated probes.",
    )
    parser.add_argument(
        "--node-volume-policy",
        choices=("barycentric", "mixed_voronoi"),
        help="Override the box-method nodal control-volume policy.",
    )
    parser.add_argument(
        "--high-field-gradient-discretization",
        choices=("edge_projection", "transport_cell_vector"),
        help="Override the high-field quasi-Fermi gradient recovery.",
    )
    args = parser.parse_args()

    baseline_rows = rows(args.baseline_state)
    if not baseline_rows or tuple(baseline_rows[0]) != STATE_COLUMNS:
        raise RuntimeError(
            "baseline header must be exactly " + ",".join(STATE_COLUMNS)
        )
    baseline = {int(row["node_id"]): row for row in baseline_rows}
    if len(baseline) != len(baseline_rows):
        raise RuntimeError("baseline state contains duplicate node IDs")
    original_baseline = {node: dict(row) for node, row in baseline.items()}

    fields = {
        "psi": scalar_field(
            args.sentaurus_fields, "ElectrostaticPotential", args.region
        ),
        "phin": scalar_field(
            args.sentaurus_fields, "eQuasiFermiPotential", args.region
        ),
        "phip": scalar_field(
            args.sentaurus_fields, "hQuasiFermiPotential", args.region
        ),
        "electrons_m3": scalar_field(
            args.sentaurus_fields, "eDensity", args.region
        ),
        "holes_m3": scalar_field(
            args.sentaurus_fields, "hDensity", args.region
        ),
    }
    imported_nodes = set(fields["psi"])
    if any(set(field) != imported_nodes for field in fields.values()):
        raise RuntimeError("Sentaurus state fields do not use the same node set")
    missing = imported_nodes - set(baseline)
    if missing:
        raise RuntimeError(f"Sentaurus export has unknown node IDs: {sorted(missing)[:8]}")

    for node_id in imported_nodes:
        row = baseline[node_id]
        row["psi"] = format(fields["psi"][node_id], ".17g")
        row["phin"] = format(fields["phin"][node_id], ".17g")
        row["phip"] = format(fields["phip"][node_id], ".17g")
        row["electrons_m3"] = format(fields["electrons_m3"][node_id] * 1.0e6, ".17g")
        row["holes_m3"] = format(fields["holes_m3"][node_id] * 1.0e6, ".17g")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=STATE_COLUMNS)
        writer.writeheader()
        for node_id in sorted(baseline):
            writer.writerow(baseline[node_id])

    if args.audit_output is not None:
        args.audit_output.parent.mkdir(parents=True, exist_ok=True)
        audit_columns = ("node_id", "psi_V", "phin_V", "phip_V", "n_m3", "p_m3")
        with args.audit_output.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=audit_columns)
            writer.writeheader()
            for node_id in sorted(baseline):
                row = baseline[node_id]
                writer.writerow({
                    "node_id": node_id,
                    "psi_V": row["psi"],
                    "phin_V": row["phin"],
                    "phip_V": row["phip"],
                    "n_m3": row["electrons_m3"],
                    "p_m3": row["holes_m3"],
                })

    if args.baseline_audit_output is not None:
        args.baseline_audit_output.parent.mkdir(parents=True, exist_ok=True)
        audit_columns = ("node_id", "psi_V", "phin_V", "phip_V", "n_m3", "p_m3")
        with args.baseline_audit_output.open(
            "w", newline="", encoding="utf-8"
        ) as stream:
            writer = csv.DictWriter(stream, fieldnames=audit_columns)
            writer.writeheader()
            for node_id in sorted(original_baseline):
                row = original_baseline[node_id]
                writer.writerow({
                    "node_id": node_id,
                    "psi_V": row["psi"],
                    "phin_V": row["phin"],
                    "phip_V": row["phip"],
                    "n_m3": row["electrons_m3"],
                    "p_m3": row["holes_m3"],
                })

    if (args.simulation_config is None) != (args.audit_config_output is None):
        raise RuntimeError(
            "--simulation-config and --audit-config-output must be provided together"
        )
    if args.simulation_config is not None:
        config = json.loads(args.simulation_config.read_text(encoding="utf-8"))
        impact = config.get("solver", {}).get("impact_ionization", {})
        if not args.preserve_eparallel and impact.get("driving_force") == "eparallel":
            impact["driving_force"] = "effective_field_parallel_j"
        if args.btbt_source_integration is not None:
            config.setdefault("solver", {}).setdefault("band_to_band", {})[
                "source_integration"
            ] = args.btbt_source_integration
        if args.node_volume_policy is not None:
            config.setdefault("mesh_geometry", {})[
                "node_volume_policy"
            ] = args.node_volume_policy
        if args.high_field_gradient_discretization is not None:
            config.setdefault("solver", {}).setdefault("mobility", {})[
                "high_field_gradient_discretization"
            ] = args.high_field_gradient_discretization
        args.audit_config_output.parent.mkdir(parents=True, exist_ok=True)
        args.audit_config_output.write_text(
            json.dumps(config, indent=2), encoding="utf-8"
        )
        if (args.probe_output_dir is None) != (args.bias is None):
            raise RuntimeError(
                "--probe-output-dir and --bias must be provided together"
            )
        if args.probe_output_dir is not None:
            args.probe_output_dir.mkdir(parents=True, exist_ok=True)
            contact_name = config.get("sweep", {}).get("contact", "drain")
            for contact in config.get("contacts", []):
                if contact.get("name") == contact_name:
                    contact["bias"] = args.bias

            sentaurus_psi_vela_qf = {
                node_id: dict(row) for node_id, row in original_baseline.items()
            }
            vela_psi_sentaurus_qf = {
                node_id: dict(row) for node_id, row in original_baseline.items()
            }
            for node_id in imported_nodes:
                sentaurus_psi_vela_qf[node_id]["psi"] = baseline[node_id]["psi"]
                vela_psi_sentaurus_qf[node_id]["phin"] = baseline[node_id]["phin"]
                vela_psi_sentaurus_qf[node_id]["phip"] = baseline[node_id]["phip"]

            probe_states = (
                ("vela", original_baseline),
                ("sentaurus_overlay", baseline),
                ("sentaurus_psi_vela_qf", sentaurus_psi_vela_qf),
                ("vela_psi_sentaurus_qf", vela_psi_sentaurus_qf),
            )
            for label, state_source in probe_states:
                fields_dir = args.probe_output_dir / f"{label}_fields"
                fields_dir.mkdir(parents=True, exist_ok=True)
                for field_name, column in (
                    ("ElectrostaticPotential", "psi"),
                    ("eQuasiFermiPotential", "phin"),
                    ("hQuasiFermiPotential", "phip"),
                ):
                    with (fields_dir / f"{field_name}_region0.csv").open(
                        "w", newline="", encoding="utf-8"
                    ) as stream:
                        writer = csv.writer(stream)
                        writer.writerow(("node_id", "component0"))
                        for node_id in sorted(state_source):
                            writer.writerow((node_id, state_source[node_id][column]))
                probe = json.loads(json.dumps(config))
                probe["simulation_type"] = "newton_carrier_term_probe"
                probe["state_fields_dir"] = str(fields_dir.resolve())
                probe["output_csv"] = str(
                    (args.probe_output_dir / f"{label}_carrier_terms.csv").resolve()
                )
                probe.pop("sweep", None)
                (args.probe_output_dir / f"{label}_carrier_terms.json").write_text(
                    json.dumps(probe, indent=2), encoding="utf-8"
                )
                srh_only_probe = json.loads(json.dumps(probe))
                srh_only_probe.setdefault("solver", {}).setdefault(
                    "band_to_band", {}
                )["model"] = "none"
                srh_only_probe["output_csv"] = str(
                    (args.probe_output_dir / f"{label}_srh_only_terms.csv").resolve()
                )
                (args.probe_output_dir / f"{label}_srh_only_terms.json").write_text(
                    json.dumps(srh_only_probe, indent=2), encoding="utf-8"
                )
                btbt_only_probe = json.loads(json.dumps(probe))
                btbt_only_probe.setdefault("solver", {})["recombination"] = ["none"]
                btbt_only_probe["output_csv"] = str(
                    (args.probe_output_dir / f"{label}_btbt_only_terms.csv").resolve()
                )
                (args.probe_output_dir / f"{label}_btbt_only_terms.json").write_text(
                    json.dumps(btbt_only_probe, indent=2), encoding="utf-8"
                )
                residual_probe = json.loads(json.dumps(probe))
                residual_probe["simulation_type"] = "newton_residual_probe"
                residual_probe["output_csv"] = str(
                    (args.probe_output_dir / f"{label}_residual.csv").resolve()
                )
                (args.probe_output_dir / f"{label}_residual.json").write_text(
                    json.dumps(residual_probe, indent=2), encoding="utf-8"
                )
                step_probe = json.loads(json.dumps(probe))
                step_probe["simulation_type"] = "newton_step_probe"
                step_probe["output_csv"] = str(
                    (args.probe_output_dir / f"{label}_step.csv").resolve()
                )
                (args.probe_output_dir / f"{label}_step.json").write_text(
                    json.dumps(step_probe, indent=2), encoding="utf-8"
                )
                mobility_probe = json.loads(json.dumps(probe))
                mobility_probe["simulation_type"] = "edge_mobility_probe"
                mobility_probe["output_csv"] = str(
                    (args.probe_output_dir / f"{label}_edge_mobility.csv").resolve()
                )
                (args.probe_output_dir / f"{label}_edge_mobility.json").write_text(
                    json.dumps(mobility_probe, indent=2), encoding="utf-8"
                )
                sg_probe = json.loads(json.dumps(probe))
                sg_probe["simulation_type"] = "sg_edge_flux_probe"
                sg_probe["output_csv"] = str(
                    (args.probe_output_dir / f"{label}_sg_edges.csv").resolve()
                )
                (args.probe_output_dir / f"{label}_sg_edges.json").write_text(
                    json.dumps(sg_probe, indent=2), encoding="utf-8"
                )

    print(
        f"wrote {len(baseline)} nodes to {args.output}; "
        f"overlaid {len(imported_nodes)} Sentaurus region-{args.region} nodes"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
