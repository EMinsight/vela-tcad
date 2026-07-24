"""Run the Minimal6 exact-support SG replacement and inversion audit."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import statistics
from collections import defaultdict
from collections.abc import Mapping, Sequence
from pathlib import Path

from .edge_flux_inversion import (
    canonical_edges,
    continuity_flux_from_current,
    edge_current_supports,
    required_positive_mobility,
    staged_sg_flux,
)
from .qfp_sg_replacement import (
    absolute_log10_error,
    symmetric_relative_residual,
)


THERMAL_VOLTAGE_300K_V = 1.380649e-23 * 300.0 / 1.602176634e-19


def _effective_intrinsic_density_m3(
    carrier: str,
    state: Mapping[int, Mapping[str, float]],
    node: int,
) -> float:
    values = state[node]
    psi = float(values["psi_V"])
    qf = float(values["qf_V"])
    density = float(values["density_m3"])
    if carrier == "electron":
        exponent = (psi - qf) / THERMAL_VOLTAGE_300K_V
    elif carrier == "hole":
        exponent = (qf - psi) / THERMAL_VOLTAGE_300K_V
    else:
        raise ValueError(f"unsupported carrier {carrier!r}")
    if not math.isfinite(density) or density <= 0.0:
        raise ValueError("carrier density must be finite and positive")
    ni = math.exp(math.log(density) - exponent)
    if not math.isfinite(ni) or ni <= 0.0:
        raise ValueError("effective intrinsic density must be finite and positive")
    return ni
EXPECTED_STATES = tuple(
    (topology, float(-bias))
    for topology in ("mirror", "sketch")
    for bias in range(1, 21)
)
SUPPORTS = (
    "endpoint_mean_tangent",
    "adjacent_cell_mean_tangent",
    "qfp_aligned_endpoint_magnitude_control",
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty CSV {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _quantile(values: Sequence[float], fraction: float) -> float | None:
    if not values:
        return None
    ordered = sorted(float(value) for value in values)
    position = (len(ordered) - 1) * fraction
    lower, upper = math.floor(position), math.ceil(position)
    if lower == upper:
        return ordered[lower]
    return ordered[lower] + (ordered[upper] - ordered[lower]) * (
        position - lower
    )


def _load_observations(
    path: Path,
) -> dict[tuple[str, str, float, int, str, str], float]:
    index: dict[tuple[str, str, float, int, str, str], float] = {}
    states: set[tuple[str, float]] = set()
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row["support_kind"] != "node" or row["status"] != "valid":
                continue
            key = (
                row["solver"],
                row["topology"],
                float(row["bias_V"]),
                int(row["support_id"]),
                row["quantity"],
                row["component"],
            )
            if key in index:
                raise ValueError(f"duplicate observation {key}")
            index[key] = float(row["value_si"])
            if (
                row["solver"] == "sentaurus"
                and row["quantity"] == "ElectrostaticPotential"
            ):
                states.add((row["topology"], float(row["bias_V"])))
    if states != set(EXPECTED_STATES):
        raise ValueError("observations differ from the exact 40-state contract")
    return index


def _value(
    index: Mapping[tuple[str, str, float, int, str, str], float],
    solver: str,
    topology: str,
    bias: float,
    node: int,
    quantity: str,
    component: str = "component0",
) -> float:
    key = (solver, topology, bias, node, quantity, component)
    if key not in index:
        raise ValueError(f"missing observation {key}")
    return float(index[key])


def _load_mesh(
    inverse_root: Path, topology: str
) -> tuple[dict[int, tuple[float, float]], tuple[tuple[int, int, int], ...]]:
    path = inverse_root / "vela" / "source" / "topologies" / topology / "mesh.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    scale = 1.0e-6 if payload.get("coordinate_unit") == "um" else 1.0
    coordinates = {
        int(row["id"]): (float(row["x"]) * scale, float(row["y"]) * scale)
        for row in payload["nodes"]
    }
    triangles = tuple(
        tuple(int(node) for node in row["node_ids"])
        for row in payload["triangles"]
    )
    if len(coordinates) != 6 or len(triangles) != 4:
        raise ValueError(f"{topology} is not a Minimal6 mesh")
    if len(canonical_edges(triangles)) != 9:
        raise ValueError(f"{topology} does not have 9 canonical edges")
    return coordinates, triangles


def _load_effective_mobility(
    path: Path,
) -> dict[tuple[str, float, str, int, int], dict[str, float]]:
    result: dict[tuple[str, float, str, int, int], dict[str, float]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            key = (
                row["topology"],
                float(row["bias_V"]),
                row["carrier"],
                int(row["node0"]),
                int(row["node1"]),
            )
            if key in result:
                raise ValueError(f"duplicate effective mobility {key}")
            result[key] = {
                "vela": float(
                    row["vela_masetti_native_state_mobility_m2_per_Vs"]
                ),
                "sentaurus_endpoint": float(
                    row["sentaurus_exported_mobility_m2_per_Vs"]
                ),
            }
    if len(result) != 720:
        raise ValueError(f"effective mobility input expected 720 rows, got {len(result)}")
    return result


def _scalar_supports(
    triangles: Sequence[Sequence[int]],
    node_values: Mapping[int, float],
) -> dict[tuple[int, int], dict[str, float]]:
    adjacent: dict[tuple[int, int], list[float]] = defaultdict(list)
    for triangle in triangles:
        cell_value = sum(node_values[int(node)] for node in triangle) / 3.0
        for edge in canonical_edges((triangle,)):
            adjacent[edge].append(cell_value)
    return {
        edge: {
            "endpoint_mean_tangent": (
                node_values[edge[0]] + node_values[edge[1]]
            )
            * 0.5,
            "adjacent_cell_mean_tangent": statistics.fmean(adjacent[edge]),
            "qfp_aligned_endpoint_magnitude_control": (
                node_values[edge[0]] + node_values[edge[1]]
            )
            * 0.5,
        }
        for edge in canonical_edges(triangles)
    }


def _node_state(
    index: Mapping[tuple[str, str, float, int, str, str], float],
    solver: str,
    topology: str,
    bias: float,
    carrier: str,
) -> dict[int, dict[str, float]]:
    qf = "eQuasiFermiPotential" if carrier == "electron" else "hQuasiFermiPotential"
    density = "eDensity" if carrier == "electron" else "hDensity"
    return {
        node: {
            "psi_V": _value(
                index, solver, topology, bias, node, "ElectrostaticPotential"
            ),
            "qf_V": _value(index, solver, topology, bias, node, qf),
            "density_m3": _value(index, solver, topology, bias, node, density),
        }
        for node in range(6)
    }


def _sentaurus_supports(
    index: Mapping[tuple[str, str, float, int, str, str], float],
    topology: str,
    bias: float,
    carrier: str,
    coordinates: Mapping[int, tuple[float, float]],
    triangles: Sequence[Sequence[int]],
    sent_state: Mapping[int, Mapping[str, float]],
) -> tuple[
    dict[tuple[int, int], dict[str, float]],
    dict[tuple[int, int], dict[str, float]],
]:
    current_quantity = (
        "eCurrentDensity" if carrier == "electron" else "hCurrentDensity"
    )
    mobility_quantity = "eMobility" if carrier == "electron" else "hMobility"
    vectors = {
        node: (
            _value(
                index,
                "sentaurus",
                topology,
                bias,
                node,
                current_quantity,
                "component0",
            ),
            _value(
                index,
                "sentaurus",
                topology,
                bias,
                node,
                current_quantity,
                "component1",
            ),
        )
        for node in range(6)
    }
    current = edge_current_supports(coordinates, triangles, vectors)
    mobility = _scalar_supports(
        triangles,
        {
            node: _value(
                index,
                "sentaurus",
                topology,
                bias,
                node,
                mobility_quantity,
            )
            for node in range(6)
        },
    )
    for edge in canonical_edges(triangles):
        length = math.dist(coordinates[edge[0]], coordinates[edge[1]])
        unit_state = {
            "psi0_V": sent_state[edge[0]]["psi_V"],
            "psi1_V": sent_state[edge[1]]["psi_V"],
            "qf0_V": sent_state[edge[0]]["qf_V"],
            "qf1_V": sent_state[edge[1]]["qf_V"],
            "density0_m3": sent_state[edge[0]]["density_m3"],
            "density1_m3": sent_state[edge[1]]["density_m3"],
            "ni0_m3": _effective_intrinsic_density_m3(carrier, sent_state, edge[0]),
            "ni1_m3": _effective_intrinsic_density_m3(carrier, sent_state, edge[1]),
            "mobility_m2_per_Vs": 1.0,
            "length_m": length,
            "thermal_voltage_V": THERMAL_VOLTAGE_300K_V,
        }
        unit_flux = staged_sg_flux("qf_sg", carrier, unit_state)
        current_sign = -1.0 if carrier == "electron" else 1.0
        if unit_flux != 0.0:
            current_sign *= math.copysign(1.0, unit_flux)
        else:
            current_sign = math.copysign(
                1.0, current[edge]["endpoint_mean_tangent"]
            )
        current[edge]["qfp_aligned_endpoint_magnitude_control"] = (
            current_sign * current[edge]["endpoint_mean_magnitude"]
        )
    return current, mobility


def _operator_state(
    formulation: str,
    carrier: str,
    edge: tuple[int, int],
    length: float,
    vela: Mapping[int, Mapping[str, float]],
    sentaurus: Mapping[int, Mapping[str, float]],
    *,
    replace_primary: bool,
    replace_psi: bool,
    mobility: float,
) -> dict[str, float]:
    psi = sentaurus if replace_psi else vela
    primary = sentaurus if replace_primary else vela
    return {
        "psi0_V": psi[edge[0]]["psi_V"],
        "psi1_V": psi[edge[1]]["psi_V"],
        "qf0_V": primary[edge[0]]["qf_V"],
        "qf1_V": primary[edge[1]]["qf_V"],
        "density0_m3": primary[edge[0]]["density_m3"],
        "density1_m3": primary[edge[1]]["density_m3"],
        "ni0_m3": _effective_intrinsic_density_m3(carrier, vela, edge[0]),
        "ni1_m3": _effective_intrinsic_density_m3(carrier, vela, edge[1]),
        "mobility_m2_per_Vs": mobility,
        "length_m": length,
        "thermal_voltage_V": THERMAL_VOLTAGE_300K_V,
    }


def _error_fields(candidate: float, reference: float) -> dict[str, object]:
    if candidate == 0.0 and reference == 0.0:
        classification = "exact_zero"
    elif candidate == 0.0 or reference == 0.0:
        classification = "zero_mismatch"
    else:
        classification = "available"
    log_error = absolute_log10_error(candidate, reference)
    sign = (
        None
        if candidate == 0.0 or reference == 0.0
        else float(math.copysign(1.0, candidate) == math.copysign(1.0, reference))
    )
    return {
        "error_classification": classification,
        "abs_log10_error_dex": "" if log_error is None else log_error,
        "symmetric_relative_residual": symmetric_relative_residual(
            candidate, reference
        ),
        "sign_agreement": "" if sign is None else sign,
    }


def _edge_scope(edge: tuple[int, int]) -> str:
    internal = {1, 5}
    count = int(edge[0] in internal) + int(edge[1] in internal)
    if count == 2:
        return "internal_edge"
    if count == 1:
        return "boundary_to_internal"
    return "boundary_edge"


def _support_audit(export_root: Path) -> dict[str, object]:
    checked = 0
    vector_fields = 0
    value_counts: set[int] = set()
    mappings: set[str] = set()
    for topology, bias in EXPECTED_STATES:
        label = f"m{abs(int(bias))}V"
        path = (
            export_root
            / "states"
            / topology
            / label
            / "export"
            / "field_manifest.json"
        )
        payload = json.loads(path.read_text(encoding="utf-8"))
        for carrier_field in ("eCurrentDensity", "hCurrentDensity"):
            matches = [
                row
                for row in payload["fields"]
                if row["name"] == carrier_field and int(row["components"]) == 2
            ]
            if len(matches) != 1:
                raise ValueError(f"{path} lacks one vector {carrier_field}")
            row = matches[0]
            vector_fields += 1
            value_counts.add(int(row["values"]))
            mappings.add(str(row["global_node_mapping"]))
        checked += 1
    if checked != 40 or value_counts != {6}:
        raise ValueError("Sentaurus current support audit did not find 40 node datasets")
    return {
        "state_count": checked,
        "carrier_vector_dataset_count": vector_fields,
        "value_counts": sorted(value_counts),
        "global_node_mappings": sorted(mappings),
        "native_directed_edge_flux": {
            "classification": "unavailable",
            "reason": (
                "TDR contains six global-node current vectors per carrier, "
                "not nine canonical directed-edge fluxes"
            ),
        },
        "p1_line_mean_identity": (
            "For linear nodal interpolation, the edge line mean is exactly "
            "the endpoint mean and is not an independent support."
        ),
    }


def _summary_rows(samples: list[dict[str, object]]) -> list[dict[str, object]]:
    groups: dict[tuple[str, ...], list[dict[str, object]]] = defaultdict(list)
    for row in samples:
        key = (
            str(row["support_mapping"]),
            str(row["carrier"]),
            str(row["formulation"]),
            str(row["branch"]),
            str(row["edge_scope"]),
        )
        groups[key].append(row)
        groups[key[:-1] + ("all_edges",)].append(row)
    output: list[dict[str, object]] = []
    for key, rows in sorted(groups.items()):
        errors = [
            float(row["abs_log10_error_dex"])
            for row in rows
            if row["abs_log10_error_dex"] != ""
        ]
        residuals = [float(row["symmetric_relative_residual"]) for row in rows]
        signs = [
            float(row["sign_agreement"])
            for row in rows
            if row["sign_agreement"] != ""
        ]
        output.append(
            {
                "support_mapping": key[0],
                "carrier": key[1],
                "formulation": key[2],
                "branch": key[3],
                "edge_scope": key[4],
                "sample_count": len(rows),
                "log_error_count": len(errors),
                "median_abs_log10_error_dex": _quantile(errors, 0.5),
                "p95_abs_log10_error_dex": _quantile(errors, 0.95),
                "median_symmetric_relative_residual": _quantile(residuals, 0.5),
                "sign_agreement_fraction": (
                    statistics.fmean(signs) if signs else None
                ),
            }
        )
    return output


def _contribution_rows(samples: list[dict[str, object]]) -> list[dict[str, object]]:
    index = {
        (
            row["topology"],
            row["bias_V"],
            row["carrier"],
            row["edge_id"],
            row["support_mapping"],
            row["formulation"],
            row["branch"],
        ): row
        for row in samples
    }
    paths = {
        "qf_sg": (
            ("vela_all", "sent_qf_only", "qf"),
            ("sent_qf_only", "sent_qf_and_mobility", "mobility"),
            ("sent_qf_and_mobility", "sent_all", "psi"),
        ),
        "density_sg": (
            ("vela_all", "sent_density_only", "density"),
            ("sent_density_only", "sent_density_and_mobility", "mobility"),
            ("sent_density_and_mobility", "sent_all", "psi"),
        ),
    }
    groups: dict[tuple[str, str, str, str], list[float]] = defaultdict(list)
    for row in samples:
        formulation = str(row["formulation"])
        if row["branch"] != "vela_all":
            continue
        base = (
            row["topology"],
            row["bias_V"],
            row["carrier"],
            row["edge_id"],
            row["support_mapping"],
            formulation,
        )
        for before, after, factor in paths[formulation]:
            left = index[base + (before,)]
            right = index[base + (after,)]
            if (
                left["abs_log10_error_dex"] == ""
                or right["abs_log10_error_dex"] == ""
            ):
                continue
            groups[
                (
                    str(row["support_mapping"]),
                    str(row["carrier"]),
                    formulation,
                    factor,
                )
            ].append(
                float(left["abs_log10_error_dex"])
                - float(right["abs_log10_error_dex"])
            )
    return [
        {
            "support_mapping": key[0],
            "carrier": key[1],
            "formulation": key[2],
            "replacement_step": key[3],
            "paired_sample_count": len(values),
            "median_paired_error_reduction_dex": _quantile(values, 0.5),
            "p05_paired_error_reduction_dex": _quantile(values, 0.05),
            "p95_paired_error_reduction_dex": _quantile(values, 0.95),
        }
        for key, values in sorted(groups.items())
    ]


def _mobility_summary(
    rows: list[dict[str, object]],
) -> list[dict[str, object]]:
    groups: dict[tuple[str, str, str], list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        groups[
            (
                str(row["support_mapping"]),
                str(row["carrier"]),
                str(row["formulation"]),
            )
        ].append(row)
    output: list[dict[str, object]] = []
    for key, selected in sorted(groups.items()):
        available = [
            row for row in selected if row["classification"] == "available"
        ]
        sent_dex = [
            float(row["required_over_sentaurus_mobility_dex"])
            for row in available
        ]
        vela_dex = [
            float(row["required_over_vela_mobility_dex"]) for row in available
        ]
        output.append(
            {
                "support_mapping": key[0],
                "carrier": key[1],
                "formulation": key[2],
                "sample_count": len(selected),
                "available_count": len(available),
                "sign_incompatible_count": sum(
                    row["classification"] == "sign_incompatible"
                    for row in selected
                ),
                "zero_operator_count": sum(
                    row["classification"] == "zero_operator" for row in selected
                ),
                "median_abs_required_over_sentaurus_mobility_dex": _quantile(
                    [abs(value) for value in sent_dex], 0.5
                ),
                "p95_abs_required_over_sentaurus_mobility_dex": _quantile(
                    [abs(value) for value in sent_dex], 0.95
                ),
                "median_required_over_sentaurus_mobility_dex": _quantile(
                    sent_dex, 0.5
                ),
                "median_required_over_vela_mobility_dex": _quantile(
                    vela_dex, 0.5
                ),
            }
        )
    return output


def _report(
    support_audit: Mapping[str, object],
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
        "# PN2D Minimal6 directed-edge SG inversion audit",
        "",
        "## Observable support",
        "",
        f"- TDR states checked: {support_audit['state_count']}.",
        "- Native Sentaurus directed-edge flux: unavailable.",
        "- Available data are six global-node current vectors per carrier.",
        "- Endpoint and adjacent-cell mappings are deterministic reconstructions.",
        "- QFP-aligned magnitude is a localization control and is ineligible for formula acceptance.",
        "",
        "## Baseline and full replacement",
        "",
        "| Support | Carrier | Formulation | Branch | Median error (dex) | p95 error (dex) | Sign agreement |",
        "|---|---|---|---|---:|---:|---:|",
    ]
    for row in selected:
        lines.append(
            "| {support_mapping} | {carrier} | {formulation} | {branch} | "
            "{median_abs_log10_error_dex:.6g} | {p95_abs_log10_error_dex:.6g} | "
            "{sign_agreement_fraction:.6g} |".format(**row)
        )
    lines.extend(
        [
            "",
            "## Sequential replacement contributions",
            "",
            "| Support | Carrier | Formulation | Step | Median paired reduction (dex) |",
            "|---|---|---|---|---:|",
        ]
    )
    for row in contributions:
        lines.append(
            "| {support_mapping} | {carrier} | {formulation} | "
            "{replacement_step} | {median_paired_error_reduction_dex:.6g} |".format(
                **row
            )
        )
    lines.extend(
        [
            "",
            "## Required effective mobility",
            "",
            "| Support | Carrier | Formulation | Available | Sign incompatible | Median abs gap to Sentaurus mobility (dex) |",
            "|---|---|---|---:|---:|---:|",
        ]
    )
    for row in mobility_summary:
        lines.append(
            "| {support_mapping} | {carrier} | {formulation} | "
            "{available_count} | {sign_incompatible_count} | "
            "{median_abs_required_over_sentaurus_mobility_dex:.6g} |".format(**row)
        )
    return "\n".join(lines) + "\n"


def run_edge_flux_experiment(
    *,
    observations_csv: str | Path,
    effective_mobility_csv: str | Path,
    inverse_inputs_root: str | Path,
    sentaurus_export_root: str | Path,
    output_root: str | Path,
) -> dict[str, object]:
    observations_path = Path(observations_csv).resolve()
    mobility_path = Path(effective_mobility_csv).resolve()
    inverse_root = Path(inverse_inputs_root).resolve()
    export_root = Path(sentaurus_export_root).resolve()
    output = Path(output_root).resolve()
    output.mkdir(parents=True, exist_ok=True)

    observations = _load_observations(observations_path)
    effective = _load_effective_mobility(mobility_path)
    support_audit = _support_audit(export_root)
    samples: list[dict[str, object]] = []
    mobility_rows: list[dict[str, object]] = []
    support_rows: list[dict[str, object]] = []

    for topology, bias in EXPECTED_STATES:
        coordinates, triangles = _load_mesh(inverse_root, topology)
        edges = canonical_edges(triangles)
        for carrier in ("electron", "hole"):
            vela = _node_state(observations, "vela", topology, bias, carrier)
            sentaurus = _node_state(
                observations, "sentaurus", topology, bias, carrier
            )
            current_supports, sent_mobility_supports = _sentaurus_supports(
                observations,
                topology,
                bias,
                carrier,
                coordinates,
                triangles,
                sentaurus,
            )
            for edge_id, edge in enumerate(edges):
                length = math.dist(coordinates[edge[0]], coordinates[edge[1]])
                native = effective[
                    (topology, bias, carrier, edge[0], edge[1])
                ]["vela"]
                endpoint_input = effective[
                    (topology, bias, carrier, edge[0], edge[1])
                ]["sentaurus_endpoint"]
                endpoint_recomputed = sent_mobility_supports[edge][
                    "endpoint_mean_tangent"
                ]
                if not math.isclose(
                    endpoint_input, endpoint_recomputed, rel_tol=5.0e-13, abs_tol=0.0
                ):
                    raise ValueError(
                        f"Sentaurus mobility projection drift at {topology} "
                        f"{bias:g} {carrier} edge {edge}"
                    )
                support_rows.append(
                    {
                        "topology": topology,
                        "bias_V": bias,
                        "carrier": carrier,
                        "edge_id": edge_id,
                        "node0": edge[0],
                        "node1": edge[1],
                        "endpoint_mean_tangent_A_per_m2": current_supports[edge][
                            "endpoint_mean_tangent"
                        ],
                        "p1_line_mean_tangent_A_per_m2": current_supports[edge][
                            "p1_line_mean_tangent"
                        ],
                        "adjacent_cell_mean_tangent_A_per_m2": current_supports[
                            edge
                        ]["adjacent_cell_mean_tangent"],
                        "qfp_aligned_endpoint_magnitude_control_A_per_m2": current_supports[
                            edge
                        ]["qfp_aligned_endpoint_magnitude_control"],
                        "p1_endpoint_identity_abs_difference_A_per_m2": abs(
                            current_supports[edge]["p1_line_mean_tangent"]
                            - current_supports[edge]["endpoint_mean_tangent"]
                        ),
                    }
                )
                for support in SUPPORTS:
                    current = current_supports[edge][support]
                    reference = continuity_flux_from_current(carrier, current)
                    sent_mobility = sent_mobility_supports[edge][support]
                    for formulation in ("qf_sg", "density_sg"):
                        primary_name = (
                            "qf" if formulation == "qf_sg" else "density"
                        )
                        branches = (
                            ("vela_all", False, False, False),
                            (f"sent_{primary_name}_only", True, False, False),
                            ("sent_psi_only", False, True, False),
                            ("sent_mobility_only", False, False, True),
                            (
                                f"sent_{primary_name}_and_psi",
                                True,
                                True,
                                False,
                            ),
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
                            mobility = sent_mobility if mobility_flag else native
                            state = _operator_state(
                                formulation,
                                carrier,
                                edge,
                                length,
                                vela,
                                sentaurus,
                                replace_primary=primary_flag,
                                replace_psi=psi_flag,
                                mobility=mobility,
                            )
                            candidate = staged_sg_flux(
                                formulation, carrier, state
                            )
                            row: dict[str, object] = {
                                "topology": topology,
                                "bias_V": bias,
                                "carrier": carrier,
                                "edge_id": edge_id,
                                "node0": edge[0],
                                "node1": edge[1],
                                "edge_scope": _edge_scope(edge),
                                "support_mapping": support,
                                "support_classification": (
                                    "localization_control"
                                    if support.endswith("_control")
                                    else "deterministic_reconstruction"
                                ),
                                "eligible_for_formula_acceptance": int(
                                    not support.endswith("_control")
                                ),
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

                        full_unit = _operator_state(
                            formulation,
                            carrier,
                            edge,
                            length,
                            vela,
                            sentaurus,
                            replace_primary=True,
                            replace_psi=True,
                            mobility=1.0,
                        )
                        unit_flux = staged_sg_flux(
                            formulation, carrier, full_unit
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
                                "support_mapping": support,
                                "support_classification": (
                                    "localization_control"
                                    if support.endswith("_control")
                                    else "deterministic_reconstruction"
                                ),
                                "formulation": formulation,
                                "classification": required["classification"],
                                "reference_continuity_flux_per_m2_s": reference,
                                "unit_mobility_operator_flux_per_m2_s": unit_flux,
                                "required_mobility_m2_per_Vs": (
                                    "" if required_value is None else required_value
                                ),
                                "sentaurus_mobility_m2_per_Vs": sent_mobility,
                                "vela_production_mobility_m2_per_Vs": native,
                                "required_over_sentaurus_mobility_dex": (
                                    ""
                                    if required_value is None
                                    else math.log10(
                                        float(required_value) / sent_mobility
                                    )
                                ),
                                "required_over_vela_mobility_dex": (
                                    ""
                                    if required_value is None
                                    else math.log10(float(required_value) / native)
                                ),
                            }
                        )

    if len(support_rows) != 720 or len(samples) != 34560:
        raise ValueError("edge flux experiment emitted an unexpected sample count")
    if len(mobility_rows) != 4320:
        raise ValueError("mobility inversion emitted an unexpected sample count")
    if max(
        float(row["p1_endpoint_identity_abs_difference_A_per_m2"])
        for row in support_rows
    ) != 0.0:
        raise ValueError("P1 line mean did not equal the endpoint mean")

    summary = _summary_rows(samples)
    contributions = _contribution_rows(samples)
    mobility_summary = _mobility_summary(mobility_rows)
    paths = {
        "support_mapping_samples.csv": support_rows,
        "sg_replacement_samples.csv": samples,
        "sg_replacement_summary.csv": summary,
        "replacement_contributions.csv": contributions,
        "mobility_inversion_samples.csv": mobility_rows,
        "mobility_inversion_summary.csv": mobility_summary,
    }
    hashes: dict[str, str] = {}
    for name, rows in paths.items():
        path = output / name
        _write_csv(path, rows)
        hashes[name] = _sha256(path)
    report_path = output / "report.md"
    report_path.write_text(
        _report(support_audit, summary, contributions, mobility_summary),
        encoding="utf-8",
    )
    hashes[report_path.name] = _sha256(report_path)
    manifest: dict[str, object] = {
        "schema_version": 1,
        "status": "valid",
        "experiment": "minimal6_directed_edge_sg_inversion",
        "state_count": 40,
        "carrier_edge_count": 720,
        "support_mapping_sample_count": len(support_rows),
        "sg_replacement_sample_count": len(samples),
        "mobility_inversion_sample_count": len(mobility_rows),
        "support_audit": support_audit,
        "acceptance_policy": {
            "native_edge_flux_required_for_formula_acceptance": True,
            "native_edge_flux_available": False,
            "formula_change_authorized": False,
            "reason": (
                "All current references are reconstructed from node vectors; "
                "the QFP-aligned magnitude branch is a localization control."
            ),
        },
        "inputs": {
            "observations_csv": str(observations_path),
            "observations_sha256": _sha256(observations_path),
            "effective_mobility_csv": str(mobility_path),
            "effective_mobility_sha256": _sha256(mobility_path),
            "inverse_inputs_root": str(inverse_root),
            "sentaurus_export_root": str(export_root),
        },
        "outputs": hashes,
    }
    manifest_path = output / "manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest
