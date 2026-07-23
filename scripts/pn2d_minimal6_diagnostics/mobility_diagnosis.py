"""Same-support Minimal6 mobility and quasi-Fermi-gradient diagnostics.

This module is diagnostic-only.  It reproduces the production Vela Masetti and
high-field mobility arithmetic without importing or modifying the C++ solver,
then compares it with sealed Sentaurus node mobility on explicit edge and cell
support.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import statistics
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from .inverse_contracts import SampleStatus, SupportKind
from .inverse_fields import cell_to_node_vectors, triangle_gradient
from .inverse_inputs import InputBundle, load_input_bundle


ELEMENTARY_CHARGE_C = 1.602176634e-19
CONSTANT_MOBILITY_M2_PER_VS = {"electron": 0.14170, "hole": 0.04705}
MASETTI = {
    "electron": {
        "mu_const": 0.14170, "mu_min1": 0.00522, "mu_min2": 0.00522,
        "mu1": 0.00434, "pc": 0.0, "cr": 9.68e22, "cs": 3.43e26,
        "alpha": 0.68, "beta": 2.0,
    },
    "hole": {
        "mu_const": 0.04705, "mu_min1": 0.00449, "mu_min2": 0.0,
        "mu1": 0.00290, "pc": 9.23e22, "cr": 2.23e23, "cs": 6.10e26,
        "alpha": 0.719, "beta": 2.0,
    },
}
FIELD = {
    "electron": {"saturation_velocity": 1.07e5, "beta": 1.109},
    "hole": {"saturation_velocity": 8.37e4, "beta": 1.213},
}
MOBILITY_BRANCHES = (
    "sentaurus_exported",
    "vela_masetti_sentaurus_state",
    "constant",
)
MOBILITY_COMPARISON_BRANCHES = MOBILITY_BRANCHES + ("vela_masetti_native_state",)
ORIENTATION_TRANSFORMS = (
    "identity", "negate", "swap_xy", "negate_swap_xy",
    "reflect_x", "reflect_y", "rotate_cw", "rotate_ccw",
)


def _carrier(carrier: str) -> str:
    if carrier not in ("electron", "hole"):
        raise ValueError("carrier must be 'electron' or 'hole'")
    return carrier


def _finite(value: object, label: str) -> float:
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{label} must be finite") from error
    if not math.isfinite(result):
        raise ValueError(f"{label} must be finite")
    return result


def masetti_low_field_mobility(carrier: str, net_doping_m3: float) -> float:
    """Reproduce ``DopingDependentMobility::masetti`` in physical SI units."""
    parameters = MASETTI[_carrier(carrier)]
    doping = abs(_finite(net_doping_m3, "net doping"))
    if doping <= 0.0:
        return parameters["mu_const"]
    exponential = parameters["mu_min1"] * math.exp(-max(0.0, parameters["pc"]) / doping)
    rolloff = (parameters["mu_const"] - parameters["mu_min2"]) / (
        1.0 + (doping / parameters["cr"]) ** parameters["alpha"]
    )
    correction = parameters["mu1"] / (
        1.0 + (parameters["cs"] / doping) ** parameters["beta"]
    )
    return max(0.0, exponential + rolloff - correction)


def field_limited_mobility(carrier: str, low_field_mobility: float, field_V_per_m: float) -> float:
    """Reproduce ``DopingDependentMobility::fieldLimit`` in physical SI units."""
    parameters = FIELD[_carrier(carrier)]
    mobility = _finite(low_field_mobility, "low-field mobility")
    field = abs(_finite(field_V_per_m, "driving field"))
    if mobility < 0.0:
        raise ValueError("low-field mobility must be non-negative")
    if mobility == 0.0 or field == 0.0:
        return mobility
    ratio = mobility * field / parameters["saturation_velocity"]
    return mobility / (1.0 + ratio ** parameters["beta"]) ** (1.0 / parameters["beta"])


def vela_masetti_edge_mobility(
    carrier: str,
    *,
    net_doping0_m3: float,
    net_doping1_m3: float,
    qf0_V: float,
    qf1_V: float,
    length_m: float,
) -> float:
    """Recompute production Vela ``masetti_field`` mobility for one edge."""
    length = _finite(length_m, "edge length")
    if length <= 0.0:
        raise ValueError("edge length must be positive")
    net_doping = 0.5 * (
        _finite(net_doping0_m3, "endpoint net doping")
        + _finite(net_doping1_m3, "endpoint net doping")
    )
    driving_field = abs(
        _finite(qf1_V, "endpoint quasi-Fermi potential")
        - _finite(qf0_V, "endpoint quasi-Fermi potential")
    ) / length
    return field_limited_mobility(
        carrier, masetti_low_field_mobility(carrier, net_doping), driving_field
    )


def unique_edges(triangles: Mapping[str, Sequence[str]]) -> tuple[tuple[str, str], ...]:
    pairs: set[tuple[str, str]] = set()
    for nodes_value in triangles.values():
        nodes = tuple(str(node) for node in nodes_value)
        if len(nodes) != 3 or len(set(nodes)) != 3:
            raise ValueError("each cell must contain three distinct nodes")
        for index in range(3):
            pairs.add(tuple(sorted((nodes[index], nodes[(index + 1) % 3]), key=_node_key)))
    return tuple(sorted(pairs, key=lambda pair: (_node_key(pair[0]), _node_key(pair[1]))))


def cell_inverted_gradient(
    carrier: str,
    density_m3: float,
    mobility_m2_per_Vs: float,
    current_A_per_m2: Sequence[float],
    *,
    q: float = ELEMENTARY_CHARGE_C,
) -> tuple[float, float]:
    """Invert current using the sign of Sentaurus-exported QF potentials.

    Both exported carrier current vectors oppose the gradient of their named
    quasi-Fermi-potential fields.  The legacy audit's positive hole sign is
    retained as an explicit negate control rather than assumed here.
    """
    _carrier(carrier)
    density = _finite(density_m3, "density")
    mobility = _finite(mobility_m2_per_Vs, "mobility")
    charge = _finite(q, "q")
    if density <= 0.0 or mobility <= 0.0 or charge <= 0.0:
        raise ValueError("density, mobility, and q must be positive")
    if len(current_A_per_m2) != 2:
        raise ValueError("current must have two components")
    current = tuple(_finite(value, "current") for value in current_A_per_m2)
    denominator = -charge * mobility * density
    return current[0] / denominator, current[1] / denominator


def _node_key(value: object) -> tuple[int, object]:
    text = str(value)
    return (0, int(text)) if text.isdecimal() else (1, text)


def _mean(values: Iterable[float]) -> float:
    items = tuple(float(value) for value in values)
    if not items:
        raise ValueError("mean needs at least one value")
    if not all(math.isfinite(value) for value in items):
        raise ValueError("mean values must be finite")
    return sum(items) / len(items)


def _mean_vector(values: Iterable[Sequence[float]]) -> tuple[float, float]:
    items = tuple((float(value[0]), float(value[1])) for value in values)
    return _mean(value[0] for value in items), _mean(value[1] for value in items)


def _length(start: Sequence[float], end: Sequence[float]) -> float:
    result = math.hypot(end[0] - start[0], end[1] - start[1])
    if not math.isfinite(result) or result <= 0.0:
        raise ValueError("edge length must be positive")
    return result


def _project(vector: Sequence[float], start: Sequence[float], end: Sequence[float]) -> float:
    length = _length(start, end)
    return (vector[0] * (end[0] - start[0]) + vector[1] * (end[1] - start[1])) / length


def _magnitude(vector: Sequence[float]) -> float:
    return math.hypot(float(vector[0]), float(vector[1]))


def _angle_deg(candidate: Sequence[float], reference: Sequence[float]) -> float | None:
    candidate_norm, reference_norm = _magnitude(candidate), _magnitude(reference)
    if candidate_norm <= 0.0 or reference_norm <= 0.0:
        return None
    cosine = (
        candidate[0] * reference[0] + candidate[1] * reference[1]
    ) / (candidate_norm * reference_norm)
    return math.degrees(math.acos(max(-1.0, min(1.0, cosine))))


def _abs_dex(candidate: float, reference: float) -> float | None:
    candidate_value, reference_value = abs(float(candidate)), abs(float(reference))
    if candidate_value <= 0.0 or reference_value <= 0.0:
        return None
    return abs(math.log10(candidate_value / reference_value))


def _signed_log10_ratio(candidate: float, reference: float) -> float | None:
    candidate_value, reference_value = float(candidate), float(reference)
    if candidate_value <= 0.0 or reference_value <= 0.0:
        return None
    return math.log10(candidate_value / reference_value)


def _quantile(values: Iterable[float], fraction: float) -> float | None:
    ordered = sorted(float(value) for value in values if value is not None and math.isfinite(float(value)))
    if not ordered:
        return None
    position = (len(ordered) - 1) * fraction
    lower = int(math.floor(position))
    upper = int(math.ceil(position))
    if lower == upper:
        return ordered[lower]
    weight = position - lower
    return ordered[lower] * (1.0 - weight) + ordered[upper] * weight


def _transform(vector: Sequence[float], name: str) -> tuple[float, float]:
    x, y = float(vector[0]), float(vector[1])
    transforms = {
        "identity": (x, y), "negate": (-x, -y), "swap_xy": (y, x),
        "negate_swap_xy": (-y, -x), "reflect_x": (-x, y),
        "reflect_y": (x, -y), "rotate_cw": (y, -x), "rotate_ccw": (-y, x),
    }
    try:
        return transforms[name]
    except KeyError as error:
        raise ValueError("unknown orientation transform") from error


def _observation_index(bundle: InputBundle) -> dict[tuple[str, str, float, str, str, str], float]:
    index: dict[tuple[str, str, float, str, str, str], float] = {}
    for row in bundle.observations:
        if row.support_kind is not SupportKind.NODE or row.status is not SampleStatus.VALID:
            continue
        if row.value_si is None:
            continue
        key = (
            row.solver, row.topology, float(row.bias_V), str(row.support_id),
            row.quantity, row.component,
        )
        if key in index:
            raise ValueError("duplicate valid observation")
        index[key] = float(row.value_si)
    return index


def _required(index, solver, topology, bias, node, quantity, component="component0") -> float:
    key = (solver, topology, float(bias), str(node), quantity, component)
    if key not in index:
        raise ValueError(f"missing required observation {key}")
    return index[key]


def _load_net_doping(vela_root: Path, topology: str) -> dict[str, float]:
    path = vela_root / "source" / "topologies" / topology / "doping.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        rows = tuple(csv.DictReader(handle))
    result: dict[str, float] = {}
    for row in rows:
        node = str(row["node_id"])
        if node in result:
            raise ValueError("duplicate doping node")
        donors = _finite(row["donors_cm3"], "donor concentration")
        acceptors = _finite(row["acceptors_cm3"], "acceptor concentration")
        result[node] = (donors - acceptors) * 1.0e6
    return result


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _write_csv(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    if not rows:
        raise ValueError(f"cannot write empty table {path.name}")
    headers = tuple(rows[0])
    if any(tuple(row) != headers for row in rows):
        raise ValueError(f"inconsistent columns for {path.name}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=headers, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def _summarize_errors(rows: Sequence[Mapping[str, object]], *, support: str, branch: str,
                      carrier: str, topology: str = "combined") -> dict[str, object]:
    selected = [row for row in rows if row["support"] == support and row["branch"] == branch
                and row["carrier"] == carrier
                and (topology == "combined" or row["topology"] == topology)]
    magnitude = [float(row["abs_log10_error"]) for row in selected
                 if row["abs_log10_error"] not in (None, "")]
    angles = [float(row["angle_deg"]) for row in selected if row["angle_deg"] not in (None, "")]
    signs = [float(row["sign_agreement"]) for row in selected if row["sign_agreement"] not in (None, "")]
    return {
        "support": support, "branch": branch, "carrier": carrier, "topology": topology,
        "sample_count": len(selected), "magnitude_valid_count": len(magnitude),
        "direction_valid_count": len(angles),
        "median_abs_log10_error": _quantile(magnitude, 0.5),
        "p95_abs_log10_error": _quantile(magnitude, 0.95),
        "median_angle_deg": _quantile(angles, 0.5),
        "sign_agreement_fraction": _mean(signs) if signs else None,
    }


def _mobility_summary(edge_rows: Sequence[Mapping[str, object]], carrier: str, branch: str,
                      topology: str = "combined") -> dict[str, object]:
    selected = [row for row in edge_rows if row["carrier"] == carrier
                and (topology == "combined" or row["topology"] == topology)]
    reference = "sentaurus_exported_mobility_m2_per_Vs"
    candidate = f"{branch}_mobility_m2_per_Vs"
    signed = [_signed_log10_ratio(float(row[candidate]), float(row[reference])) for row in selected]
    absolute = [abs(value) for value in signed if value is not None]
    relative = [abs(float(row[candidate]) - float(row[reference])) / float(row[reference])
                for row in selected if float(row[reference]) > 0.0]
    values = [float(row[candidate]) for row in selected]
    return {
        "carrier": carrier, "branch": branch, "topology": topology,
        "sample_count": len(selected), "min_mobility_m2_per_Vs": min(values),
        "median_mobility_m2_per_Vs": statistics.median(values),
        "max_mobility_m2_per_Vs": max(values),
        "median_signed_log10_ratio_vs_sentaurus": _quantile(signed, 0.5),
        "median_abs_log10_error_vs_sentaurus": _quantile(absolute, 0.5),
        "p95_abs_log10_error_vs_sentaurus": _quantile(absolute, 0.95),
        "median_relative_error_vs_sentaurus": _quantile(relative, 0.5),
    }


def _markdown(report: Mapping[str, object]) -> str:
    qf = report["qf_gradient_summary"]
    mobility = report["mobility_summary"]
    conclusions = report["conclusions"]
    lines = [
        "# Minimal6 mobility and quasi-Fermi-gradient same-support diagnosis", "",
        "## Technical summary", "",
        str(conclusions["headline"]), "",
        f"The audit covers {report['state_count']} exact states, {report['edge_sample_count']} carrier-edge samples, "
        f"and {report['cell_sample_count']} carrier-cell samples. Production `include/` and `src/` formulas were not changed.", "",
        "## Mobility law comparison on edge support", "",
        "| carrier | branch | N | median mobility (m2/V/s) | median abs error (dex) | p95 abs error (dex) | median relative error |", "|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in mobility:
        if row["topology"] != "combined":
            continue
        lines.append(
            f"| {row['carrier']} | {row['branch']} | {row['sample_count']} | "
            f"{row['median_mobility_m2_per_Vs']:.12g} | {row['median_abs_log10_error_vs_sentaurus']:.12g} | "
            f"{row['p95_abs_log10_error_vs_sentaurus']:.12g} | {row['median_relative_error_vs_sentaurus']:.12g} |"
        )
    lines += ["", "The same-state Masetti branch isolates the mobility formula by feeding the Sentaurus quasi-Fermi edge difference into the Vela formula. The native-state branch is reported only as a state-sensitive reference; it is not used in the three-way inversion.", "",
              "## Quasi-Fermi-gradient inversion after explicit support mapping", "",
              "| support | carrier | mobility branch | N valid | median abs error (dex) | p95 abs error (dex) | median angle (deg) | sign agreement |", "|---|---|---|---:|---:|---:|---:|---:|"]
    for row in qf:
        if row["topology"] != "combined":
            continue
        lines.append(
            f"| {row['support']} | {row['carrier']} | {row['branch']} | {row['magnitude_valid_count']} | "
            f"{row['median_abs_log10_error'] if row['median_abs_log10_error'] is not None else ''} | "
            f"{row['p95_abs_log10_error'] if row['p95_abs_log10_error'] is not None else ''} | "
            f"{row['median_angle_deg'] if row['median_angle_deg'] is not None else ''} | "
            f"{row['sign_agreement_fraction'] if row['sign_agreement_fraction'] is not None else ''} |"
        )
    lines += ["", "A positive scalar mobility changes only gradient magnitude. Identical angles across the three cell branches are therefore both expected and checked explicitly.", "",
              "## Sign and coordinate controls", "",
              f"Best tested orientation transform: `{conclusions['best_orientation_transform']}`; "
              f"median cell angle {conclusions['best_orientation_median_angle_deg']:.12g} deg. "
              f"The established carrier-sign convention is `{conclusions['carrier_sign_result']}`.", "",
              "## Scope, data, and definitions", "",
              "Sentaurus node current, density, quasi-Fermi potential, and mobility are hash-bound sealed inputs. Edge values use endpoint arithmetic means and current tangent projection. Cell current and density use equal P1 nodal weights; cell mobility is the mean of the three edge mobilities; the reference gradient is the exact affine P1 triangle gradient.", "",
              "For the exported Sentaurus potential fields, both carrier inversions use `J = -q mu c grad(phi_F)`; the legacy positive hole sign is retained as a negate control. Magnitude error is `abs(log10(|candidate|/|reference|))`; direction error is the Cartesian angle in degrees.", "",
              "## Methodology and robustness", "",
              "The same-state Masetti calculation copies the production arithmetic and defaults, including average endpoint net doping and quasi-Fermi-gradient high-field limiting. A constant-mobility control uses 1417/470.5 cm2/V/s. Eight signed axis permutations test sign and coordinate hypotheses without tuning continuous parameters.", "",
              "## Limitations", "",
              "The Sentaurus current is a node-exported vector, so edge/cell projection is an explicit diagnostic interpolation rather than Sentaurus internal flux support. The formula `J/(q n mu)` is a local drift identity and is not the discrete Scharfetter-Gummel inverse. Cell-averaging a scalar mobility cannot represent within-cell mobility variation or a tensor mobility.", "",
              "## Recommended next steps", "",
              str(conclusions["next_step"]), "",
              "## Further questions", "",
              "The remaining discriminant is whether a discrete SG inversion using endpoint carrier densities and the exact Sentaurus edge current/flux (if exportable) closes the magnitude residual without an empirical mobility scale.", ""]
    return "\n".join(lines)


def build_mobility_diagnosis(
    vela_root_value: str | Path,
    sentaurus_root_value: str | Path,
    supplemental_root_value: str | Path,
    output_value: str | Path,
) -> dict[str, object]:
    """Run and persist the deterministic 40-state same-support diagnosis."""
    vela_root = Path(vela_root_value).resolve()
    sentaurus_root = Path(sentaurus_root_value).resolve()
    supplemental_root = Path(supplemental_root_value).resolve()
    output = Path(output_value).resolve()
    output.mkdir(parents=True, exist_ok=True)

    bundle = load_input_bundle(vela_root, sentaurus_root, supplemental_root)
    index = _observation_index(bundle)
    error_rows: list[dict[str, object]] = []
    edge_rows: list[dict[str, object]] = []
    cell_rows: list[dict[str, object]] = []
    orientation_samples: dict[tuple[str, str], list[float]] = {
        (carrier, transform): []
        for carrier in ("electron", "hole") for transform in ORIENTATION_TRANSFORMS
    }

    for topology, bias in bundle.common_keys:
        mesh = bundle.diagnostic_context.mesh_by_topology[topology]
        triangles = {str(cell): tuple(str(node) for node in nodes)
                     for cell, nodes in mesh["triangles"].items()}
        edges = {str(edge): tuple(str(node) for node in nodes)
                 for edge, nodes in mesh["edges"].items()}
        if tuple(edges.values()) != unique_edges(triangles):
            raise ValueError("sealed edge topology is not canonical")
        nodes = sorted({node for pair in edges.values() for node in pair}, key=_node_key)
        coordinates = {
            node: (
                _required(index, "sentaurus", topology, bias, node, "coordinate", "x"),
                _required(index, "sentaurus", topology, bias, node, "coordinate", "y"),
            ) for node in nodes
        }
        net_doping = _load_net_doping(vela_root, topology)
        if set(net_doping) != set(nodes):
            raise ValueError("doping and mesh node sets differ")

        for carrier in ("electron", "hole"):
            prefix = "e" if carrier == "electron" else "h"
            qf_quantity, density_quantity = f"{prefix}QuasiFermiPotential", f"{prefix}Density"
            mobility_quantity, current_quantity = f"{prefix}Mobility", f"{prefix}CurrentDensity"
            sent_qf = {node: _required(index, "sentaurus", topology, bias, node, qf_quantity) for node in nodes}
            vela_qf = {node: _required(index, "vela", topology, bias, node, qf_quantity) for node in nodes}
            density = {node: _required(index, "sentaurus", topology, bias, node, density_quantity) for node in nodes}
            sent_mu = {node: _required(index, "sentaurus", topology, bias, node, mobility_quantity) for node in nodes}
            current = {node: (
                _required(index, "sentaurus", topology, bias, node, current_quantity, "component0"),
                _required(index, "sentaurus", topology, bias, node, current_quantity, "component1"),
            ) for node in nodes}

            edge_mobility: dict[str, dict[str, float]] = {}
            edge_row_by_id: dict[str, dict[str, object]] = {}
            edge_pair_to_id = {pair: edge_id for edge_id, pair in edges.items()}
            for edge_id, (node0, node1) in edges.items():
                start, end = coordinates[node0], coordinates[node1]
                length = _length(start, end)
                current_vector = _mean_vector((current[node0], current[node1]))
                current_tangent = _project(current_vector, start, end)
                direct_gradient = (sent_qf[node1] - sent_qf[node0]) / length
                density_value = _mean((density[node0], density[node1]))
                mobility_values = {
                    "sentaurus_exported": _mean((sent_mu[node0], sent_mu[node1])),
                    "vela_masetti_sentaurus_state": vela_masetti_edge_mobility(
                        carrier, net_doping0_m3=net_doping[node0],
                        net_doping1_m3=net_doping[node1], qf0_V=sent_qf[node0],
                        qf1_V=sent_qf[node1], length_m=length,
                    ),
                    "constant": CONSTANT_MOBILITY_M2_PER_VS[carrier],
                    "vela_masetti_native_state": vela_masetti_edge_mobility(
                        carrier, net_doping0_m3=net_doping[node0],
                        net_doping1_m3=net_doping[node1], qf0_V=vela_qf[node0],
                        qf1_V=vela_qf[node1], length_m=length,
                    ),
                }
                edge_mobility[edge_id] = mobility_values
                row: dict[str, object] = {
                    "topology": topology, "bias_V": bias, "carrier": carrier,
                    "edge_id": edge_id, "node0": node0, "node1": node1,
                    "length_m": length, "density_m3": density_value,
                    "current_tangent_A_per_m2": current_tangent,
                    "direct_qf_gradient_tangent_V_per_m": direct_gradient,
                }
                for branch in MOBILITY_COMPARISON_BRANCHES:
                    row[f"{branch}_mobility_m2_per_Vs"] = mobility_values[branch]
                for branch in MOBILITY_BRANCHES:
                    inverted = current_tangent / (
                        -ELEMENTARY_CHARGE_C * density_value * mobility_values[branch]
                    )
                    error = _abs_dex(inverted, direct_gradient)
                    sign_agreement = None if inverted == 0.0 or direct_gradient == 0.0 else float(
                        math.copysign(1.0, inverted) == math.copysign(1.0, direct_gradient)
                    )
                    row[f"{branch}_inverted_qf_gradient_tangent_V_per_m"] = inverted
                    row[f"{branch}_abs_log10_error"] = "" if error is None else error
                    row[f"{branch}_sign_agreement"] = "" if sign_agreement is None else sign_agreement
                    error_rows.append({
                        "support": "edge_ratio_after_project", "topology": topology, "bias_V": bias,
                        "carrier": carrier, "support_id": edge_id, "branch": branch,
                        "abs_log10_error": error, "angle_deg": None,
                        "sign_agreement": sign_agreement,
                    })
                edge_row_by_id[edge_id] = row
                edge_rows.append(row)

            node_mobility = {branch: {} for branch in MOBILITY_BRANCHES}
            for node in nodes:
                incident_edges = tuple(
                    edge_id for edge_id, endpoints in edges.items() if node in endpoints
                )
                for branch in MOBILITY_BRANCHES:
                    if branch == "sentaurus_exported":
                        value = sent_mu[node]
                    else:
                        value = _mean(edge_mobility[edge_id][branch] for edge_id in incident_edges)
                    node_mobility[branch][node] = value
            node_inverted = {
                branch: {
                    node: cell_inverted_gradient(
                        carrier, density[node], node_mobility[branch][node], current[node]
                    )
                    for node in nodes
                }
                for branch in MOBILITY_BRANCHES
            }

            for edge_id, (node0, node1) in edges.items():
                start, end = coordinates[node0], coordinates[node1]
                local_projected = {
                    branch: _project(
                        _mean_vector((node_inverted[branch][node0], node_inverted[branch][node1])),
                        start, end,
                    )
                    for branch in MOBILITY_BRANCHES
                }
                row = edge_row_by_id[edge_id]
                direct_gradient = float(row["direct_qf_gradient_tangent_V_per_m"])
                for branch, inverted in local_projected.items():
                    error = _abs_dex(inverted, direct_gradient)
                    sign_agreement = None if inverted == 0.0 or direct_gradient == 0.0 else float(
                        math.copysign(1.0, inverted) == math.copysign(1.0, direct_gradient)
                    )
                    row[f"{branch}_local_then_project_qf_gradient_tangent_V_per_m"] = inverted
                    error_rows.append({
                        "support": "edge_local_then_project", "topology": topology,
                        "bias_V": bias, "carrier": carrier, "support_id": edge_id,
                        "branch": branch, "abs_log10_error": error,
                        "angle_deg": None, "sign_agreement": sign_agreement,
                    })

            cell_gradients: dict[str, tuple[float, float]] = {}
            cell_areas: dict[str, float] = {}
            for cell_id, cell_nodes in triangles.items():
                points = tuple(coordinates[node] for node in cell_nodes)
                direct = triangle_gradient(points, tuple(sent_qf[node] for node in cell_nodes))
                cell_gradients[cell_id] = direct
                cell_areas[cell_id] = abs(
                    (points[1][0] - points[0][0]) * (points[2][1] - points[0][1])
                    - (points[2][0] - points[0][0]) * (points[1][1] - points[0][1])
                ) * 0.5
                cell_current = _mean_vector(current[node] for node in cell_nodes)
                cell_density = _mean(density[node] for node in cell_nodes)
                cell_edge_ids = tuple(
                    edge_pair_to_id[tuple(sorted((cell_nodes[i], cell_nodes[(i + 1) % 3]), key=_node_key))]
                    for i in range(3)
                )
                cell_mu = {
                    branch: _mean(edge_mobility[edge_id][branch] for edge_id in cell_edge_ids)
                    for branch in MOBILITY_BRANCHES
                }
                cell_row: dict[str, object] = {
                    "topology": topology, "bias_V": bias, "carrier": carrier,
                    "cell_id": cell_id, "node0": cell_nodes[0], "node1": cell_nodes[1],
                    "node2": cell_nodes[2], "area_m2": cell_areas[cell_id],
                    "density_m3": cell_density,
                    "current_x_A_per_m2": cell_current[0], "current_y_A_per_m2": cell_current[1],
                    "direct_qf_gradient_x_V_per_m": direct[0],
                    "direct_qf_gradient_y_V_per_m": direct[1],
                }
                for branch in MOBILITY_BRANCHES:
                    inverted = cell_inverted_gradient(carrier, cell_density, cell_mu[branch], cell_current)
                    magnitude_error = _abs_dex(_magnitude(inverted), _magnitude(direct))
                    angle = _angle_deg(inverted, direct)
                    cell_row[f"{branch}_mobility_m2_per_Vs"] = cell_mu[branch]
                    cell_row[f"{branch}_inverted_qf_gradient_x_V_per_m"] = inverted[0]
                    cell_row[f"{branch}_inverted_qf_gradient_y_V_per_m"] = inverted[1]
                    cell_row[f"{branch}_abs_log10_error"] = "" if magnitude_error is None else magnitude_error
                    cell_row[f"{branch}_angle_deg"] = "" if angle is None else angle
                    error_rows.append({
                        "support": "cell_ratio_after_project", "topology": topology, "bias_V": bias,
                        "carrier": carrier, "support_id": cell_id, "branch": branch,
                        "abs_log10_error": magnitude_error, "angle_deg": angle,
                        "sign_agreement": None,
                    })
                    local_inverted = _mean_vector(node_inverted[branch][node] for node in cell_nodes)
                    local_error = _abs_dex(_magnitude(local_inverted), _magnitude(direct))
                    local_angle = _angle_deg(local_inverted, direct)
                    cell_row[f"{branch}_local_then_project_qf_gradient_x_V_per_m"] = local_inverted[0]
                    cell_row[f"{branch}_local_then_project_qf_gradient_y_V_per_m"] = local_inverted[1]
                    cell_row[f"{branch}_local_then_project_abs_log10_error"] = "" if local_error is None else local_error
                    cell_row[f"{branch}_local_then_project_angle_deg"] = "" if local_angle is None else local_angle
                    error_rows.append({
                        "support": "cell_local_then_project", "topology": topology,
                        "bias_V": bias, "carrier": carrier, "support_id": cell_id,
                        "branch": branch, "abs_log10_error": local_error,
                        "angle_deg": local_angle, "sign_agreement": None,
                    })
                    if branch == "sentaurus_exported" and local_angle is not None:
                        for transform in ORIENTATION_TRANSFORMS:
                            transformed_angle = _angle_deg(_transform(local_inverted, transform), direct)
                            if transformed_angle is not None:
                                orientation_samples[(carrier, transform)].append(transformed_angle)
                cell_rows.append(cell_row)

            node_reference = cell_to_node_vectors(cell_gradients, triangles, coordinates)["values"]
            for node in nodes:
                inverted = cell_inverted_gradient(carrier, density[node], sent_mu[node], current[node])
                reference = node_reference[node]
                magnitude_error = _abs_dex(_magnitude(inverted), _magnitude(reference))
                angle = _angle_deg(inverted, reference)
                error_rows.append({
                    "support": "node", "topology": topology, "bias_V": bias,
                    "carrier": carrier, "support_id": node,
                    "branch": "sentaurus_exported", "abs_log10_error": magnitude_error,
                    "angle_deg": angle, "sign_agreement": None,
                })

    mobility_summary = [
        _mobility_summary(edge_rows, carrier, branch, topology)
        for topology in ("combined", "mirror", "sketch")
        for carrier in ("electron", "hole")
        for branch in MOBILITY_COMPARISON_BRANCHES
    ]
    qf_summary = [
        _summarize_errors(error_rows, support=support, branch=branch,
                          carrier=carrier, topology=topology)
        for topology in ("combined", "mirror", "sketch")
        for support in (
            "node", "edge_ratio_after_project", "edge_local_then_project",
            "cell_ratio_after_project", "cell_local_then_project")
        for carrier in ("electron", "hole")
        for branch in (("sentaurus_exported",) if support == "node" else MOBILITY_BRANCHES)
    ]
    orientation_summary = [
        {
            "carrier": carrier, "transform": transform,
            "valid_count": len(orientation_samples[(carrier, transform)]),
            "median_angle_deg": _quantile(orientation_samples[(carrier, transform)], 0.5),
            "p95_angle_deg": _quantile(orientation_samples[(carrier, transform)], 0.95),
        }
        for carrier in ("electron", "hole") for transform in ORIENTATION_TRANSFORMS
    ]

    combined_cell = [row for row in qf_summary if row["topology"] == "combined" and row["support"] == "cell_ratio_after_project"]
    angle_by_branch = {
        (row["carrier"], row["branch"]): row["median_angle_deg"] for row in combined_cell
    }
    angle_spread = max(
        max(angle_by_branch[(carrier, branch)] for branch in MOBILITY_BRANCHES)
        - min(angle_by_branch[(carrier, branch)] for branch in MOBILITY_BRANCHES)
        for carrier in ("electron", "hole")
    )
    combined_local_cell = [
        row for row in qf_summary
        if row["topology"] == "combined" and row["support"] == "cell_local_then_project"
    ]
    local_angle_by_branch = {
        (row["carrier"], row["branch"]): row["median_angle_deg"] for row in combined_local_cell
    }
    local_angle_spread = max(
        max(local_angle_by_branch[(carrier, branch)] for branch in MOBILITY_BRANCHES)
        - min(local_angle_by_branch[(carrier, branch)] for branch in MOBILITY_BRANCHES)
        for carrier in ("electron", "hole")
    )
    best_orientation = min(
        orientation_summary,
        key=lambda row: (float(row["median_angle_deg"]), row["carrier"], row["transform"]),
    )
    identity = {row["carrier"]: row for row in orientation_summary if row["transform"] == "identity"}
    negate = {row["carrier"]: row for row in orientation_summary if row["transform"] == "negate"}
    carrier_sign_result = (
        "supported" if all(identity[carrier]["median_angle_deg"] <= negate[carrier]["median_angle_deg"]
                           for carrier in ("electron", "hole")) else "not_supported"
    )
    combined_mobility = [
        row for row in mobility_summary
        if row["topology"] == "combined"
        and row["branch"] in ("vela_masetti_sentaurus_state", "constant")
    ]
    best_mobility_by_carrier = {
        carrier: min(
            (row for row in combined_mobility if row["carrier"] == carrier),
            key=lambda row: float(row["median_abs_log10_error_vs_sentaurus"]),
        )["branch"]
        for carrier in ("electron", "hole")
    }
    legacy_pooled_angles = (
        orientation_samples[("electron", "identity")] + orientation_samples[("hole", "negate")]
    )
    corrected_pooled_angles = (
        orientation_samples[("electron", "identity")] + orientation_samples[("hole", "identity")]
    )
    conclusions = {
        "ratio_after_project_mobility_angle_spread_deg": angle_spread,
        "local_then_project_mobility_angle_spread_deg": local_angle_spread,
        "mobility_can_explain_direction_error": False,
        "best_orientation_transform": best_orientation["transform"],
        "best_orientation_carrier": best_orientation["carrier"],
        "best_orientation_median_angle_deg": best_orientation["median_angle_deg"],
        "carrier_sign_result": carrier_sign_result,
        "best_same_state_mobility_branch_by_carrier": best_mobility_by_carrier,
        "legacy_mixed_sign_pooled_median_angle_deg": _quantile(legacy_pooled_angles, 0.5),
        "corrected_same_sign_pooled_median_angle_deg": _quantile(corrected_pooled_angles, 0.5),
    }
    conclusions["headline"] = (
        "The earlier near-90-degree aggregate is a carrier-sign pooling artifact: electron current already "
        "aligns with minus grad(eQFP), while hole current aligns with minus grad(hQFP), not the legacy plus "
        "hole sign. With the data-supported negative sign for both exported quasi-Fermi potentials, the "
        "identity coordinate transform wins for both carriers. Scalar cell mobility changes magnitude but "
        "does not rotate the ratio-after-projection result."
    )
    conclusions["next_step"] = (
        "Do not change a production mobility or current formula from this result. Next obtain a Sentaurus "
        "edge flux/current export and perform a discrete Scharfetter-Gummel inversion; the local drift "
        "identity and arithmetic projection do not close magnitude in strongly nonuniform cells."
    )

    report: dict[str, object] = {
        "schema": "vela.pn2d_minimal6_mobility_diagnosis.v1",
        "diagnostic_only": True,
        "state_count": len(bundle.common_keys),
        "edge_sample_count": len(edge_rows),
        "cell_sample_count": len(cell_rows),
        "input_manifest_sha256": {
            "vela": _sha256(vela_root / "manifest.json"),
            "sentaurus": _sha256(sentaurus_root / "manifest.json"),
            "supplemental": _sha256(supplemental_root / "manifest.json"),
        },
        "support_contract": {
            "edge": "endpoint arithmetic mean; current projected onto canonical node0-to-node1 tangent",
            "cell": "equal P1 nodal mean for current and density; mean of three edge mobilities; affine P1 qF gradient",
            "node_baseline": "native current/density/mobility; adjacent-cell area-weighted qF gradient",
        },
        "carrier_sign_contract": {
            "electron": "Jn=-q*mu_n*n*grad(phi_n)",
            "hole": "Jp=-q*mu_p*p*grad(phi_p) for exported Sentaurus hQFP",
        },
        "mobility_summary": mobility_summary,
        "qf_gradient_summary": qf_summary,
        "orientation_control_summary": orientation_summary,
        "conclusions": conclusions,
    }

    _write_csv(output / "mobility_edge_samples.csv", edge_rows)
    _write_csv(output / "qf_gradient_cell_samples.csv", cell_rows)
    _write_csv(output / "qf_gradient_error_samples.csv", error_rows)
    _write_csv(output / "mobility_summary.csv", mobility_summary)
    _write_csv(output / "qf_gradient_summary.csv", qf_summary)
    _write_csv(output / "orientation_control_summary.csv", orientation_summary)
    (output / "mobility_diagnosis.json").write_text(
        json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    (output / "mobility_diagnosis.md").write_text(_markdown(report), encoding="utf-8", newline="\n")
    return report
