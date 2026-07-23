"""Run the Minimal6 SG audit against native Sentaurus element currents."""

from __future__ import annotations

import csv
import json
import math
from collections import defaultdict
from pathlib import Path

from .edge_flux_experiment import (
    EXPECTED_STATES,
    _contribution_rows,
    _edge_scope,
    _error_fields,
    _load_effective_mobility,
    _load_mesh,
    _load_observations,
    _mobility_summary,
    _node_state,
    _operator_state,
    _scalar_supports,
    _sha256,
    _summary_rows,
    _value,
    _write_csv,
)
from .edge_flux_inversion import (
    canonical_edges,
    continuity_flux_from_current,
    required_positive_mobility,
    staged_sg_flux,
)


SUPPORT = "native_element_vector_adjacent_mean_tangent"


def _load_cell_currents(
    path: Path,
) -> dict[tuple[str, float, str, int], tuple[float, float]]:
    result: dict[tuple[str, float, str, int], tuple[float, float]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            key = (
                row["topology"],
                float(row["bias_V"]),
                row["carrier"],
                int(row["cell_id"]),
            )
            if key in result:
                raise ValueError(f"duplicate native element current {key}")
            result[key] = (
                float(row["current_x_A_per_m2"]),
                float(row["current_y_A_per_m2"]),
            )
    expected = {
        (topology, bias, carrier, cell)
        for topology, bias in EXPECTED_STATES
        for carrier in ("electron", "hole")
        for cell in range(4)
    }
    if set(result) != expected:
        raise ValueError("native element currents differ from 40 x 2 x 4")
    return result


def _edge_cells(
    triangles: tuple[tuple[int, int, int], ...],
) -> dict[tuple[int, int], tuple[int, ...]]:
    adjacent: dict[tuple[int, int], list[int]] = defaultdict(list)
    for cell, triangle in enumerate(triangles):
        for edge in canonical_edges((triangle,)):
            adjacent[edge].append(cell)
    return {
        edge: tuple(adjacent[edge])
        for edge in canonical_edges(triangles)
    }


def _tangent(
    coordinates: dict[int, tuple[float, float]], edge: tuple[int, int]
) -> tuple[float, float]:
    x0, y0 = coordinates[edge[0]]
    x1, y1 = coordinates[edge[1]]
    length = math.hypot(x1 - x0, y1 - y0)
    return (x1 - x0) / length, (y1 - y0) / length


def _native_edge_current(
    currents: dict[tuple[str, float, str, int], tuple[float, float]],
    topology: str,
    bias: float,
    carrier: str,
    cells: tuple[int, ...],
    tangent: tuple[float, float],
) -> tuple[float, tuple[float, ...]]:
    projections = tuple(
        currents[(topology, bias, carrier, cell)][0] * tangent[0]
        + currents[(topology, bias, carrier, cell)][1] * tangent[1]
        for cell in cells
    )
    return sum(projections) / len(projections), projections


def _report(
    summary: list[dict[str, object]],
    contributions: list[dict[str, object]],
    mobility_summary: list[dict[str, object]],
) -> str:
    selected = [
        row
        for row in summary
        if row["edge_scope"] == "all_edges"
        and row["branch"] in {"vela_all", "sent_all"}
    ]
    lines = [
        "# PN2D Minimal6 native-element current SG inversion",
        "",
        "## Support",
        "",
        "- Sentaurus source: native four-triangle e/h current vectors.",
        "- Directed-edge reference: adjacent native element vectors averaged and projected on the canonical edge tangent.",
        "- A native Sentaurus directed-edge flux is still not present.",
        "",
        "## Baseline and full replacement",
        "",
        "| Carrier | Formulation | Branch | Median error (dex) | p95 error (dex) | Sign agreement |",
        "|---|---|---|---:|---:|---:|",
    ]
    for row in selected:
        lines.append(
            "| {carrier} | {formulation} | {branch} | "
            "{median_abs_log10_error_dex:.6g} | "
            "{p95_abs_log10_error_dex:.6g} | "
            "{sign_agreement_fraction:.6g} |".format(**row)
        )
    lines.extend(
        [
            "",
            "## Sequential replacement contribution",
            "",
            "| Carrier | Formulation | Step | Median paired reduction (dex) |",
            "|---|---|---|---:|",
        ]
    )
    for row in contributions:
        lines.append(
            "| {carrier} | {formulation} | {replacement_step} | "
            "{median_paired_error_reduction_dex:.6g} |".format(**row)
        )
    lines.extend(
        [
            "",
            "## Required mobility",
            "",
            "| Carrier | Formulation | Available | Sign incompatible | Median abs gap to Sentaurus mobility (dex) |",
            "|---|---|---:|---:|---:|",
        ]
    )
    for row in mobility_summary:
        lines.append(
            "| {carrier} | {formulation} | {available_count} | "
            "{sign_incompatible_count} | "
            "{median_abs_required_over_sentaurus_mobility_dex:.6g} |".format(
                **row
            )
        )
    return "\n".join(lines) + "\n"


def run_native_cell_sg_experiment(
    *,
    observations_csv: str | Path,
    effective_mobility_csv: str | Path,
    inverse_inputs_root: str | Path,
    element_currents_csv: str | Path,
    element_currents_manifest: str | Path,
    output_root: str | Path,
) -> dict[str, object]:
    observations_path = Path(observations_csv).resolve()
    mobility_path = Path(effective_mobility_csv).resolve()
    inverse_root = Path(inverse_inputs_root).resolve()
    currents_path = Path(element_currents_csv).resolve()
    currents_manifest_path = Path(element_currents_manifest).resolve()
    output = Path(output_root).resolve()
    output.mkdir(parents=True, exist_ok=True)

    observations = _load_observations(observations_path)
    effective = _load_effective_mobility(mobility_path)
    currents = _load_cell_currents(currents_path)
    currents_manifest = json.loads(
        currents_manifest_path.read_text(encoding="utf-8")
    )
    if (
        currents_manifest.get("status") != "valid"
        or currents_manifest.get("state_count") != 40
        or currents_manifest.get("sample_count") != 320
    ):
        raise ValueError("native element-current manifest is not valid")

    support_rows: list[dict[str, object]] = []
    samples: list[dict[str, object]] = []
    mobility_rows: list[dict[str, object]] = []
    for topology, bias in EXPECTED_STATES:
        coordinates, triangles = _load_mesh(inverse_root, topology)
        adjacent = _edge_cells(triangles)
        edges = canonical_edges(triangles)
        for carrier in ("electron", "hole"):
            vela = _node_state(observations, "vela", topology, bias, carrier)
            sentaurus = _node_state(
                observations, "sentaurus", topology, bias, carrier
            )
            mobility_quantity = (
                "eMobility" if carrier == "electron" else "hMobility"
            )
            sent_mobility = _scalar_supports(
                triangles,
                {
                    node: _value(
                        observations,
                        "sentaurus",
                        topology,
                        bias,
                        node,
                        mobility_quantity,
                    )
                    for node in range(6)
                },
            )
            for edge_id, edge in enumerate(edges):
                tangent = _tangent(coordinates, edge)
                current, projections = _native_edge_current(
                    currents,
                    topology,
                    bias,
                    carrier,
                    adjacent[edge],
                    tangent,
                )
                reference = continuity_flux_from_current(carrier, current)
                length = math.dist(coordinates[edge[0]], coordinates[edge[1]])
                native_mobility = effective[
                    (topology, bias, carrier, edge[0], edge[1])
                ]["vela"]
                mapped_sent_mobility = sent_mobility[edge][
                    "adjacent_cell_mean_tangent"
                ]
                support_rows.append(
                    {
                        "topology": topology,
                        "bias_V": bias,
                        "carrier": carrier,
                        "edge_id": edge_id,
                        "node0": edge[0],
                        "node1": edge[1],
                        "adjacent_cell_ids": ";".join(
                            str(cell) for cell in adjacent[edge]
                        ),
                        "adjacent_cell_count": len(adjacent[edge]),
                        "cell_tangent_currents_A_per_m2": ";".join(
                            format(value, ".17g") for value in projections
                        ),
                        "native_element_mean_tangent_A_per_m2": current,
                        "internal_cell_projection_spread_A_per_m2": (
                            0.0
                            if len(projections) == 1
                            else abs(projections[0] - projections[1])
                        ),
                    }
                )
                for formulation in ("qf_sg", "density_sg"):
                    primary_name = "qf" if formulation == "qf_sg" else "density"
                    branches = (
                        ("vela_all", False, False, False),
                        (f"sent_{primary_name}_only", True, False, False),
                        ("sent_psi_only", False, True, False),
                        ("sent_mobility_only", False, False, True),
                        (f"sent_{primary_name}_and_psi", True, True, False),
                        (
                            f"sent_{primary_name}_and_mobility",
                            True,
                            False,
                            True,
                        ),
                        ("sent_psi_and_mobility", False, True, True),
                        ("sent_all", True, True, True),
                    )
                    for branch, primary_flag, psi_flag, mobility_flag in branches:
                        mobility = (
                            mapped_sent_mobility
                            if mobility_flag
                            else native_mobility
                        )
                        state = _operator_state(
                            formulation,
                            edge,
                            length,
                            vela,
                            sentaurus,
                            replace_primary=primary_flag,
                            replace_psi=psi_flag,
                            mobility=mobility,
                        )
                        candidate = staged_sg_flux(formulation, carrier, state)
                        row: dict[str, object] = {
                            "topology": topology,
                            "bias_V": bias,
                            "carrier": carrier,
                            "edge_id": edge_id,
                            "node0": edge[0],
                            "node1": edge[1],
                            "edge_scope": _edge_scope(edge),
                            "support_mapping": SUPPORT,
                            "support_classification": (
                                "native_cell_to_edge_reconstruction"
                            ),
                            "eligible_for_formula_acceptance": 0,
                            "formulation": formulation,
                            "branch": branch,
                            "sent_primary": int(primary_flag),
                            "sent_psi": int(psi_flag),
                            "sent_mobility": int(mobility_flag),
                            "mobility_m2_per_Vs": mobility,
                            "candidate_continuity_flux_per_m2_s": candidate,
                            "reference_continuity_flux_per_m2_s": reference,
                        }
                        row.update(_error_fields(candidate, reference))
                        samples.append(row)

                    unit_state = _operator_state(
                        formulation,
                        edge,
                        length,
                        vela,
                        sentaurus,
                        replace_primary=True,
                        replace_psi=True,
                        mobility=1.0,
                    )
                    unit_flux = staged_sg_flux(
                        formulation, carrier, unit_state
                    )
                    required = required_positive_mobility(
                        reference_flux=reference,
                        unit_mobility_flux=unit_flux,
                    )
                    required_value = required["mobility_m2_per_Vs"]
                    mobility_rows.append(
                        {
                            "topology": topology,
                            "bias_V": bias,
                            "carrier": carrier,
                            "edge_id": edge_id,
                            "node0": edge[0],
                            "node1": edge[1],
                            "edge_scope": _edge_scope(edge),
                            "support_mapping": SUPPORT,
                            "support_classification": (
                                "native_cell_to_edge_reconstruction"
                            ),
                            "formulation": formulation,
                            "classification": required["classification"],
                            "reference_continuity_flux_per_m2_s": reference,
                            "unit_mobility_operator_flux_per_m2_s": unit_flux,
                            "required_mobility_m2_per_Vs": (
                                "" if required_value is None else required_value
                            ),
                            "sentaurus_mobility_m2_per_Vs": mapped_sent_mobility,
                            "vela_production_mobility_m2_per_Vs": native_mobility,
                            "required_over_sentaurus_mobility_dex": (
                                ""
                                if required_value is None
                                else math.log10(
                                    float(required_value) / mapped_sent_mobility
                                )
                            ),
                            "required_over_vela_mobility_dex": (
                                ""
                                if required_value is None
                                else math.log10(
                                    float(required_value) / native_mobility
                                )
                            ),
                        }
                    )

    if len(support_rows) != 720:
        raise ValueError("native support sample count differs from 720")
    if len(samples) != 11520 or len(mobility_rows) != 1440:
        raise ValueError("native SG audit emitted an unexpected sample count")
    summary = _summary_rows(samples)
    contributions = _contribution_rows(samples)
    mobility_summary = _mobility_summary(mobility_rows)
    tables = {
        "native_element_support_samples.csv": support_rows,
        "native_cell_sg_replacement_samples.csv": samples,
        "native_cell_sg_replacement_summary.csv": summary,
        "native_cell_replacement_contributions.csv": contributions,
        "native_cell_mobility_inversion_samples.csv": mobility_rows,
        "native_cell_mobility_inversion_summary.csv": mobility_summary,
    }
    hashes: dict[str, str] = {}
    for name, rows in tables.items():
        path = output / name
        _write_csv(path, rows)
        hashes[name] = _sha256(path)
    report = output / "report.md"
    report.write_text(
        _report(summary, contributions, mobility_summary), encoding="utf-8"
    )
    hashes[report.name] = _sha256(report)
    manifest: dict[str, object] = {
        "schema_version": 1,
        "status": "valid",
        "experiment": "minimal6_native_element_current_sg_inversion",
        "state_count": 40,
        "support_sample_count": len(support_rows),
        "sg_replacement_sample_count": len(samples),
        "mobility_inversion_sample_count": len(mobility_rows),
        "support_audit": {
            "native_element_current_vectors_available": True,
            "native_directed_edge_flux_available": False,
            "edge_mapping": SUPPORT,
        },
        "acceptance_policy": {
            "formula_change_authorized": False,
            "reason": (
                "Sentaurus element vectors are native, but the shared "
                "directed-edge reference still requires cell averaging and "
                "tangent projection."
            ),
        },
        "inputs": {
            "observations_csv": str(observations_path),
            "observations_sha256": _sha256(observations_path),
            "effective_mobility_csv": str(mobility_path),
            "effective_mobility_sha256": _sha256(mobility_path),
            "element_currents_csv": str(currents_path),
            "element_currents_sha256": _sha256(currents_path),
            "element_currents_manifest": str(currents_manifest_path),
            "element_currents_manifest_sha256": _sha256(
                currents_manifest_path
            ),
            "inverse_inputs_root": str(inverse_root),
        },
        "outputs": hashes,
    }
    manifest_path = output / "manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest
