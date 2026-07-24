"""Forty-state Minimal6 Vela-to-Sentaurus box-current staged replay."""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
import statistics
from pathlib import Path


Q_SENT_C = 1.6021918e-19
Q_SI_C = 1.602176634e-19
VT_V = 1.380649e-23 * 300.0 / Q_SI_C
STAGES = (
    "vela_baseline",
    "sentaurus_qfp",
    "sentaurus_qfp_density",
    "sentaurus_qfp_density_element_mobility",
    "sentaurus_qfp_density_element_mobility_geometry",
)
CONTROL_STAGE = "sentaurus_qfp_recomputed_density_control"
QUANTITIES = {
    "psi_V": "ElectrostaticPotential",
    "phin_V": "eQuasiFermiPotential",
    "phip_V": "hQuasiFermiPotential",
    "n_m3": "eDensity",
    "p_m3": "hDensity",
}


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


def _bernoulli(value: float) -> float:
    if abs(value) < 1.0e-10:
        return 1.0 - value * 0.5 + value * value / 12.0
    if value > 500.0:
        return value * math.exp(-value)
    if value < -500.0:
        return -value
    return value / math.expm1(value)


def _limited_exp(value: float) -> float:
    return math.exp(max(-500.0, min(500.0, value)))


def _element_current_A_per_um(
    carrier: str,
    state: dict[int, dict[str, float]],
    start: int,
    end: int,
    mobility_cm2_V_s: float,
    kappa: float,
) -> float:
    first = state[start]
    second = state[end]
    if carrier == "electron":
        qf_delta = (first["phin_V"] - second["phin_V"]) / VT_V
        x = math.log(first["n_m3"] / second["n_m3"]) + qf_delta
        bracket = -first["n_m3"] * _bernoulli(x) * math.expm1(qf_delta)
    elif carrier == "hole":
        qf_delta = (first["phip_V"] - second["phip_V"]) / VT_V
        x = math.log(first["p_m3"] / second["p_m3"]) - qf_delta
        bracket = first["p_m3"] * _bernoulli(x) * math.expm1(-qf_delta)
    else:
        raise ValueError(f"unsupported carrier {carrier!r}")
    # m^-3 to cm^-3 contributes 1e-6; A/cm to A/um contributes 1e-4.
    return Q_SENT_C * VT_V * kappa * mobility_cm2_V_s * bracket * 1.0e-10


def _effective_ni(state: dict[str, float], carrier: str) -> float:
    if carrier == "electron":
        density = state["n_m3"]
        exponent = (state["psi_V"] - state["phin_V"]) / VT_V
    else:
        density = state["p_m3"]
        exponent = (state["phip_V"] - state["psi_V"]) / VT_V
    return math.exp(math.log(density) - exponent)


def _load_observations(
    path: Path,
) -> dict[tuple[str, str, float], dict[int, dict[str, float]]]:
    flat: dict[tuple[str, str, float, int, str], float] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if (
                row["support_kind"] != "node"
                or row["status"] != "valid"
                or row["component"] != "component0"
            ):
                continue
            quantity = row["quantity"]
            if quantity not in QUANTITIES.values():
                continue
            key = (
                row["solver"],
                row["topology"],
                float(row["bias_V"]),
                int(row["support_id"]),
                quantity,
            )
            if key in flat:
                raise ValueError(f"duplicate observation {key}")
            flat[key] = float(row["value_si"])
    states: dict[tuple[str, str, float], dict[int, dict[str, float]]] = {}
    for solver in ("vela", "sentaurus"):
        for topology in ("mirror", "sketch"):
            for magnitude in range(1, 21):
                bias = -float(magnitude)
                nodes: dict[int, dict[str, float]] = {}
                for node in range(6):
                    nodes[node] = {}
                    for field, quantity in QUANTITIES.items():
                        key = (solver, topology, bias, node, quantity)
                        if key not in flat:
                            raise ValueError(f"missing observation {key}")
                        nodes[node][field] = flat[key]
                states[(solver, topology, bias)] = nodes
    if len(states) != 80:
        raise ValueError("observation contract must contain 80 solver states")
    return states


def _load_transport(
    path: Path,
) -> dict[tuple[str, float, int], dict[str, float]]:
    output: dict[tuple[str, float, int], dict[str, float]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            key = (row["topology"], float(row["bias_V"]), int(row["cell_id"]))
            if key in output:
                raise ValueError(f"duplicate transport element {key}")
            output[key] = {
                "electron_cm2_V_s": float(row["eMobility_component0_raw"]),
                "hole_cm2_V_s": float(row["hMobility_component0_raw"]),
            }
    if len(output) != 160:
        raise ValueError("transport contract must contain 160 element rows")
    return output


def _load_vela_mobility(
    path: Path,
) -> dict[tuple[str, float, str, tuple[int, int]], float | None]:
    output: dict[tuple[str, float, str, tuple[int, int]], float | None] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            pair = tuple(sorted((int(row["node0"]), int(row["node1"]))))
            key = (row["topology"], float(row["bias_V"]), row["carrier"], pair)
            raw = row["baseline_mobility_m2_per_Vs"].strip()
            output[key] = None if not raw else float(raw)
    if len(output) != 720:
        raise ValueError("Vela mobility contract must contain 720 carrier-edge rows")
    return output


def _load_mesh(
    path: Path,
) -> tuple[
    list[dict[str, object]],
    dict[tuple[int, int], float],
    dict[tuple[int, int], int],
    dict[str, set[int]],
]:
    mesh = json.loads(path.read_text(encoding="utf-8"))
    scale = 1.0e-6 if mesh["coordinate_unit"] == "um" else 1.0
    coordinates = {
        int(node["id"]): (float(node["x"]) * scale, float(node["y"]) * scale)
        for node in mesh["nodes"]
    }
    local_rows: list[dict[str, object]] = []
    kappa_by_local: dict[tuple[int, int], float] = {}
    pairs: set[tuple[int, int]] = set()
    for triangle in mesh["triangles"]:
        element = int(triangle["id"])
        nodes = [int(value) for value in triangle["node_ids"]]
        definitions = (
            (nodes[0], nodes[1], nodes[2]),
            (nodes[0], nodes[2], nodes[1]),
            (nodes[1], nodes[2], nodes[0]),
        )
        for local_edge, (raw_start, raw_end, opposite) in enumerate(definitions):
            start, end = sorted((raw_start, raw_end))
            sx, sy = coordinates[start]
            ex, ey = coordinates[end]
            ox, oy = coordinates[opposite]
            ax, ay = sx - ox, sy - oy
            bx, by = ex - ox, ey - oy
            cross = abs(ax * by - ay * bx)
            dot = ax * bx + ay * by
            kappa = 0.5 * dot / cross
            if abs(kappa) < 1.0e-15:
                kappa = 0.0
            pair = (start, end)
            pairs.add(pair)
            kappa_by_local[(element, local_edge)] = kappa
            local_rows.append(
                {
                    "element": element,
                    "local_edge": local_edge,
                    "node0": start,
                    "node1": end,
                    "opposite_node": opposite,
                    "pair": pair,
                    "kappa": kappa,
                }
            )
    pair_to_edge = {pair: edge_id for edge_id, pair in enumerate(sorted(pairs))}
    if len(local_rows) != 12 or len(pair_to_edge) != 9:
        raise ValueError("mesh must contain 12 local and 9 global edges")
    contacts = {
        contact["name"]: {int(node) for node in contact["node_ids"]}
        for contact in mesh["contacts"]
    }
    return local_rows, kappa_by_local, pair_to_edge, contacts


def _load_vela_edges(path: Path) -> dict[tuple[int, int], dict[str, float]]:
    output: dict[tuple[int, int], dict[str, float]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            pair = tuple(sorted((int(row["node0"]), int(row["node1"]))))
            output[pair] = {
                "vela_edge_id": int(row["edge_id"]),
                "length_m": float(row["length_m"]),
                "electron_flux_m2_s": float(
                    row["electron_raw_signed_flux_per_m2_s"]
                ),
                "hole_flux_m2_s": float(row["hole_raw_signed_flux_per_m2_s"]),
            }
    if len(output) != 9:
        raise ValueError(f"Vela edge file {path} must contain 9 rows")
    return output


def _parse_final_plt(path: Path) -> dict[str, float]:
    text = path.read_text(encoding="utf-8", errors="replace")
    info, data = text.split("Data {", 1)
    dataset_block = info.split("datasets", 1)[1].split("]", 1)[0]
    names = re.findall(r'"([^"]+)"', dataset_block)
    values = [
        float(value)
        for value in re.findall(
            r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[Ee][-+]?\d+)?", data
        )
    ]
    if not names or len(values) % len(names) != 0:
        raise ValueError(f"cannot parse final PLT record {path}")
    final = values[-len(names) :]
    return dict(zip(names, final))


def _contact_current(
    edge_values: dict[tuple[int, int], float],
    contact_nodes: set[int],
) -> float:
    outward = 0.0
    for (start, end), value in edge_values.items():
        start_in = start in contact_nodes
        end_in = end in contact_nodes
        if start_in == end_in:
            continue
        outward += value if start_in else -value
    return -outward


def _node_divergence(
    edge_values: dict[tuple[int, int], float], node: int
) -> float:
    result = 0.0
    for (start, end), value in edge_values.items():
        if start == node:
            result += value
        elif end == node:
            result -= value
    return result


def _metric(values: list[float]) -> dict[str, float | int | None]:
    return {
        "count": len(values),
        "median": None if not values else statistics.median(values),
        "p95": None if not values else _quantile(values, 0.95),
        "maximum": None if not values else max(values),
    }


def _quantile(values: list[float], fraction: float) -> float:
    ordered = sorted(values)
    position = (len(ordered) - 1) * fraction
    lower = math.floor(position)
    upper = math.ceil(position)
    if lower == upper:
        return ordered[lower]
    return ordered[lower] + (ordered[upper] - ordered[lower]) * (
        position - lower
    )


def run_box_staged_sweep(
    *,
    observations_csv: str | Path,
    transport_elements_csv: str | Path,
    vela_mobility_csv: str | Path,
    vela_replay_root: str | Path,
    mesh_root: str | Path,
    sentaurus_state_root: str | Path,
    output_root: str | Path,
) -> dict[str, object]:
    observations_path = Path(observations_csv).resolve()
    transport_path = Path(transport_elements_csv).resolve()
    vela_mobility_path = Path(vela_mobility_csv).resolve()
    vela_root = Path(vela_replay_root).resolve()
    meshes = Path(mesh_root).resolve()
    sentaurus_root = Path(sentaurus_state_root).resolve()
    output = Path(output_root).resolve()
    output.mkdir(parents=True, exist_ok=True)

    observations = _load_observations(observations_path)
    transport = _load_transport(transport_path)
    vela_mobility = _load_vela_mobility(vela_mobility_path)
    mesh_data = {
        topology: _load_mesh(meshes / topology / "mesh.json")
        for topology in ("mirror", "sketch")
    }

    sample_rows: list[dict[str, object]] = []
    density_rows: list[dict[str, object]] = []
    mobility_rows: list[dict[str, object]] = []
    geometry_rows: list[dict[str, object]] = []
    terminal_rows: list[dict[str, object]] = []
    kcl_rows: list[dict[str, object]] = []
    baseline_rows: list[dict[str, object]] = []
    state_rows: list[dict[str, object]] = []

    for topology in ("mirror", "sketch"):
        local_rows, kappa_by_local, pair_to_edge, contacts = mesh_data[topology]
        for local in local_rows:
            geometry_rows.append(
                {
                    "topology": topology,
                    "element": local["element"],
                    "local_edge": local["local_edge"],
                    "node0": local["node0"],
                    "node1": local["node1"],
                    "kappa": local["kappa"],
                    "status": "geometric_zero"
                    if float(local["kappa"]) == 0.0
                    else "valid",
                }
            )

        for magnitude in range(1, 21):
            bias = -float(magnitude)
            label = f"m{magnitude}V"
            vela_state = observations[("vela", topology, bias)]
            sent_state = observations[("sentaurus", topology, bias)]
            edge_path = vela_root / topology / label / "edges.csv"
            state_path = vela_root / topology / label / "state.csv"
            if not edge_path.is_file() or not state_path.is_file():
                raise FileNotFoundError(f"missing Vela replay state {topology} {bias:g}")
            vela_edges = _load_vela_edges(edge_path)

            ni: dict[int, float] = {}
            recomputed: dict[int, dict[str, float]] = {}
            for node in range(6):
                ni_e = _effective_ni(vela_state[node], "electron")
                ni_h = _effective_ni(vela_state[node], "hole")
                ni[node] = math.sqrt(ni_e * ni_h)
                recomputed[node] = dict(sent_state[node])
                recomputed[node]["n_m3"] = ni[node] * _limited_exp(
                    (sent_state[node]["psi_V"] - sent_state[node]["phin_V"])
                    / VT_V
                )
                recomputed[node]["p_m3"] = ni[node] * _limited_exp(
                    (sent_state[node]["phip_V"] - sent_state[node]["psi_V"])
                    / VT_V
                )
                density_rows.append(
                    {
                        "topology": topology,
                        "bias_V": bias,
                        "node": node,
                        "vela_ni_e_m3": ni_e,
                        "vela_ni_h_m3": ni_h,
                        "vela_ni_used_m3": ni[node],
                        "recomputed_n_m3": recomputed[node]["n_m3"],
                        "sentaurus_n_m3": sent_state[node]["n_m3"],
                        "n_abs_dex": abs(
                            math.log10(
                                recomputed[node]["n_m3"]
                                / sent_state[node]["n_m3"]
                            )
                        ),
                        "recomputed_p_m3": recomputed[node]["p_m3"],
                        "sentaurus_p_m3": sent_state[node]["p_m3"],
                        "p_abs_dex": abs(
                            math.log10(
                                recomputed[node]["p_m3"]
                                / sent_state[node]["p_m3"]
                            )
                        ),
                    }
                )

            stage_states = {
                "vela_baseline": vela_state,
                "sentaurus_qfp": {
                    node: {
                        **vela_state[node],
                        "phin_V": sent_state[node]["phin_V"],
                        "phip_V": sent_state[node]["phip_V"],
                    }
                    for node in range(6)
                },
                "sentaurus_qfp_density": {
                    node: {
                        **vela_state[node],
                        "phin_V": sent_state[node]["phin_V"],
                        "phip_V": sent_state[node]["phip_V"],
                        "n_m3": sent_state[node]["n_m3"],
                        "p_m3": sent_state[node]["p_m3"],
                    }
                    for node in range(6)
                },
                "sentaurus_qfp_density_element_mobility": sent_state,
                "sentaurus_qfp_density_element_mobility_geometry": sent_state,
                CONTROL_STAGE: recomputed,
            }

            reference: dict[tuple[str, tuple[int, int]], float] = {
                (carrier, pair): 0.0
                for carrier in ("electron", "hole")
                for pair in pair_to_edge
            }
            candidates: dict[tuple[str, str, tuple[int, int]], float] = {
                (stage, carrier, pair): 0.0
                for stage in (*STAGES, CONTROL_STAGE)
                for carrier in ("electron", "hole")
                for pair in pair_to_edge
            }

            for local in local_rows:
                element = int(local["element"])
                local_edge = int(local["local_edge"])
                pair = local["pair"]
                start, end = pair
                kappa = kappa_by_local[(element, local_edge)]
                for carrier in ("electron", "hole"):
                    sent_mu = transport[(topology, bias, element)][
                        f"{carrier}_cm2_V_s"
                    ]
                    final_value = _element_current_A_per_um(
                        carrier, sent_state, start, end, sent_mu, kappa
                    )
                    reference[(carrier, pair)] += final_value
                    for stage, state in stage_states.items():
                        use_sent_mu = stage in (
                            "sentaurus_qfp_density_element_mobility",
                            "sentaurus_qfp_density_element_mobility_geometry",
                        )
                        if use_sent_mu:
                            mobility_cm2 = sent_mu
                        else:
                            value = vela_mobility[
                                (topology, bias, carrier, pair)
                            ]
                            mobility_cm2 = 0.0 if value is None else value * 1.0e4
                        candidates[(stage, carrier, pair)] += (
                            _element_current_A_per_um(
                                carrier,
                                state,
                                start,
                                end,
                                mobility_cm2,
                                kappa,
                            )
                        )

            for pair in pair_to_edge:
                adjacent = [
                    local for local in local_rows if local["pair"] == pair
                ]
                kappa_sum = sum(float(local["kappa"]) for local in adjacent)
                for carrier in ("electron", "hole"):
                    value = vela_mobility[(topology, bias, carrier, pair)]
                    weighted_sent = (
                        None
                        if kappa_sum == 0.0
                        else sum(
                            float(local["kappa"])
                            * transport[
                                (topology, bias, int(local["element"]))
                            ][f"{carrier}_cm2_V_s"]
                            for local in adjacent
                        )
                        / kappa_sum
                    )
                    mobility_rows.append(
                        {
                            "topology": topology,
                            "bias_V": bias,
                            "carrier": carrier,
                            "edge_id": pair_to_edge[pair],
                            "node0": pair[0],
                            "node1": pair[1],
                            "vela_mobility_m2_V_s": ""
                            if value is None
                            else value,
                            "sentaurus_kappa_weighted_mobility_m2_V_s": ""
                            if weighted_sent is None
                            else weighted_sent * 1.0e-4,
                            "abs_log10_ratio_dex": ""
                            if value is None or weighted_sent is None
                            else abs(math.log10(value * 1.0e4 / weighted_sent)),
                        }
                    )

            state_error_index: dict[tuple[str, str], list[float]] = {}
            state_sign_index: dict[tuple[str, str], list[float]] = {}
            for stage in (*STAGES, CONTROL_STAGE):
                for carrier in ("electron", "hole"):
                    errors: list[float] = []
                    signs: list[float] = []
                    for pair, edge_id in pair_to_edge.items():
                        ref = reference[(carrier, pair)]
                        candidate = candidates[(stage, carrier, pair)]
                        if ref == 0.0 and candidate == 0.0:
                            status = "exact_zero"
                            error: float | str = ""
                            sign: float | str = ""
                        elif ref == 0.0:
                            status = "reference_zero_candidate_nonzero"
                            error = ""
                            sign = ""
                        elif candidate == 0.0:
                            status = "candidate_zero"
                            error = ""
                            sign = 0.0
                            signs.append(0.0)
                        else:
                            status = "valid"
                            error = abs(
                                math.log10(abs(candidate) / abs(ref))
                            )
                            sign = float(
                                math.copysign(1.0, candidate)
                                == math.copysign(1.0, ref)
                            )
                            errors.append(float(error))
                            signs.append(float(sign))
                        sample_rows.append(
                            {
                                "topology": topology,
                                "bias_V": bias,
                                "stage": stage,
                                "carrier": carrier,
                                "edge_id": edge_id,
                                "node0": pair[0],
                                "node1": pair[1],
                                "reference_A_per_um": ref,
                                "candidate_A_per_um": candidate,
                                "absolute_log10_error_dex": error,
                                "sign_agreement": sign,
                                "status": status,
                            }
                        )
                    state_error_index[(stage, carrier)] = errors
                    state_sign_index[(stage, carrier)] = signs
                    metric = _metric(errors)
                    state_rows.append(
                        {
                            "topology": topology,
                            "bias_V": bias,
                            "stage": stage,
                            "carrier": carrier,
                            "valid_count": metric["count"],
                            "median_abs_dex": metric["median"],
                            "maximum_abs_dex": metric["maximum"],
                            "sign_agreement_fraction": ""
                            if not signs
                            else statistics.mean(signs),
                        }
                    )

            kappa_by_pair = {
                pair: sum(
                    float(local["kappa"])
                    for local in local_rows
                    if local["pair"] == pair
                )
                for pair in pair_to_edge
            }
            for pair, edge in vela_edges.items():
                dual_length = kappa_by_pair[pair] * edge["length_m"]
                for carrier, charge_sign in (
                    ("electron", 1.0),
                    ("hole", -1.0),
                ):
                    production = (
                        charge_sign
                        * Q_SENT_C
                        * edge[f"{carrier}_flux_m2_s"]
                        * dual_length
                        * 1.0e-6
                    )
                    replay = candidates[
                        ("vela_baseline", carrier, pair)
                    ]
                    scale = max(abs(production), abs(replay), 1.0e-300)
                    baseline_rows.append(
                        {
                            "topology": topology,
                            "bias_V": bias,
                            "carrier": carrier,
                            "edge_id": pair_to_edge[pair],
                            "node0": pair[0],
                            "node1": pair[1],
                            "production_A_per_um": production,
                            "replay_A_per_um": replay,
                            "relative_difference": abs(production - replay)
                            / scale,
                        }
                    )

            plt_path = (
                sentaurus_root
                / topology
                / label
                / f"pn2d_minimal6_state_{label}.plt"
            )
            terminal = _parse_final_plt(plt_path)
            reference_edges = {
                carrier: {
                    pair: reference[(carrier, pair)] for pair in pair_to_edge
                }
                for carrier in ("electron", "hole")
            }
            for contact in ("Anode", "Cathode"):
                for carrier, field in (
                    ("electron", "eCurrent"),
                    ("hole", "hCurrent"),
                ):
                    predicted = _contact_current(
                        reference_edges[carrier], contacts[contact]
                    )
                    observed = terminal[f"{contact} {field}"]
                    terminal_rows.append(
                        {
                            "topology": topology,
                            "bias_V": bias,
                            "contact": contact,
                            "carrier": carrier,
                            "predicted_A_per_um": predicted,
                            "sentaurus_A_per_um": observed,
                            "relative_error": abs(predicted - observed)
                            / max(abs(observed), 1.0e-300),
                            "absolute_log10_error_dex": ""
                            if predicted == 0.0 or observed == 0.0
                            else abs(
                                math.log10(
                                    abs(predicted) / abs(observed)
                                )
                            ),
                        }
                    )
            total_reference = {
                pair: reference[("electron", pair)]
                + reference[("hole", pair)]
                for pair in pair_to_edge
            }
            terminal_scale = max(
                abs(terminal["Anode TotalCurrent"]),
                abs(terminal["Cathode TotalCurrent"]),
                1.0e-300,
            )
            for node in (1, 5):
                residual = _node_divergence(total_reference, node)
                kcl_rows.append(
                    {
                        "topology": topology,
                        "bias_V": bias,
                        "node": node,
                        "total_current_divergence_A_per_um": residual,
                        "relative_to_terminal_total": abs(residual)
                        / terminal_scale,
                    }
                )

    summary_rows: list[dict[str, object]] = []
    for scope in ("all", "mirror", "sketch"):
        for stage in (*STAGES, CONTROL_STAGE):
            for carrier in ("electron", "hole"):
                selected = [
                    row
                    for row in sample_rows
                    if row["stage"] == stage
                    and row["carrier"] == carrier
                    and row["status"] == "valid"
                    and (scope == "all" or row["topology"] == scope)
                ]
                errors = [
                    float(row["absolute_log10_error_dex"])
                    for row in selected
                ]
                signs = [float(row["sign_agreement"]) for row in selected]
                metric = _metric(errors)
                summary_rows.append(
                    {
                        "scope": scope,
                        "stage": stage,
                        "carrier": carrier,
                        "valid_count": metric["count"],
                        "median_abs_dex": metric["median"],
                        "p95_abs_dex": metric["p95"],
                        "maximum_abs_dex": metric["maximum"],
                        "sign_agreement_fraction": statistics.mean(signs),
                    }
                )

    sample_index = {
        (
            row["topology"],
            float(row["bias_V"]),
            row["stage"],
            row["carrier"],
            int(row["edge_id"]),
        ): row
        for row in sample_rows
    }
    contribution_rows: list[dict[str, object]] = []
    for scope in ("all", "mirror", "sketch"):
        for carrier in ("electron", "hole"):
            for previous, current in zip(STAGES, STAGES[1:]):
                paired: list[float] = []
                for topology in ("mirror", "sketch"):
                    if scope != "all" and topology != scope:
                        continue
                    for magnitude in range(1, 21):
                        bias = -float(magnitude)
                        for edge_id in range(9):
                            before = sample_index[
                                (
                                    topology,
                                    bias,
                                    previous,
                                    carrier,
                                    edge_id,
                                )
                            ]
                            after = sample_index[
                                (
                                    topology,
                                    bias,
                                    current,
                                    carrier,
                                    edge_id,
                                )
                            ]
                            if (
                                before["status"] == "valid"
                                and after["status"] == "valid"
                            ):
                                paired.append(
                                    float(
                                        before[
                                            "absolute_log10_error_dex"
                                        ]
                                    )
                                    - float(
                                        after[
                                            "absolute_log10_error_dex"
                                        ]
                                    )
                                )
                metric = _metric(paired)
                contribution_rows.append(
                    {
                        "scope": scope,
                        "carrier": carrier,
                        "previous_stage": previous,
                        "current_stage": current,
                        "paired_count": metric["count"],
                        "median_error_reduction_dex": metric["median"],
                        "p95_error_reduction_dex": metric["p95"],
                        "minimum_error_reduction_dex": min(paired),
                        "maximum_error_reduction_dex": metric["maximum"],
                    }
                )

    outputs = {
        "stage_edge_samples.csv": sample_rows,
        "stage_summary.csv": summary_rows,
        "paired_contributions.csv": contribution_rows,
        "state_summary.csv": state_rows,
        "density_recompute_control.csv": density_rows,
        "mobility_comparison.csv": mobility_rows,
        "geometry_coefficients.csv": geometry_rows,
        "terminal_closure.csv": terminal_rows,
        "total_current_kcl.csv": kcl_rows,
        "baseline_operator_crosscheck.csv": baseline_rows,
    }
    for name, rows in outputs.items():
        _write_csv(output / name, rows)

    final_rows = [
        row
        for row in sample_rows
        if row["stage"] == STAGES[-1] and row["status"] == "valid"
    ]
    final_relative = [
        abs(float(row["candidate_A_per_um"]) - float(row["reference_A_per_um"]))
        / abs(float(row["reference_A_per_um"]))
        for row in final_rows
    ]
    manifest: dict[str, object] = {
        "schema_version": 1,
        "status": "valid",
        "experiment": "minimal6_sentaurus_box_staged_replacement_40_state",
        "state_contract": {
            "state_count": 40,
            "topologies": ["mirror", "sketch"],
            "biases_V": [-float(value) for value in range(1, 21)],
            "nodes_per_state": 6,
            "global_edges_per_state": 9,
            "elements_per_state": 4,
            "element_local_edges_per_state": 12,
        },
        "stage_order": list(STAGES),
        "control_stage": CONTROL_STAGE,
        "current_convention": {
            "electron": "Sentaurus qfp-plus box branch",
            "hole": "Sentaurus qfp-minus box branch",
            "sentaurus_elementary_charge_C": Q_SENT_C,
            "thermal_voltage_V": VT_V,
        },
        "gates": {
            "final_valid_carrier_edge_count": len(final_rows),
            "final_max_relative_error": max(final_relative),
            "terminal_closure_max_relative_error": max(
                float(row["relative_error"]) for row in terminal_rows
            ),
            "total_current_kcl_max_relative_error": max(
                float(row["relative_to_terminal_total"]) for row in kcl_rows
            ),
            "baseline_operator_max_relative_difference": max(
                float(row["relative_difference"]) for row in baseline_rows
            ),
            "recomputed_density_max_abs_dex": max(
                max(float(row["n_abs_dex"]), float(row["p_abs_dex"]))
                for row in density_rows
            ),
            "geometric_zero_local_edge_count": sum(
                row["status"] == "geometric_zero" for row in geometry_rows
            ),
        },
        "limitations": [
            "the directed-edge reference is a documented Sentaurus box-operator reconstruction",
            "native Sentaurus directed-edge current remains unavailable",
            "mirror ReadCoefficient was observed directly; sketch coefficients follow the same documented cotangent geometry",
        ],
        "inputs": {
            "observations_csv": str(observations_path),
            "observations_sha256": _sha256(observations_path),
            "transport_elements_csv": str(transport_path),
            "transport_elements_sha256": _sha256(transport_path),
            "vela_mobility_csv": str(vela_mobility_path),
            "vela_mobility_sha256": _sha256(vela_mobility_path),
            "vela_replay_root": str(vela_root),
            "mesh_root": str(meshes),
            "sentaurus_state_root": str(sentaurus_root),
        },
        "outputs": {},
    }
    for name in outputs:
        manifest["outputs"][name] = _sha256(output / name)
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    return manifest
