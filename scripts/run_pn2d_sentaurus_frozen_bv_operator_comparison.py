#!/usr/bin/env python3
"""Replay frozen Sentaurus states through Vela's production BV operator.

The generated audit configuration always uses impact-ionization
``coupling_mode=postprocess_only``.  The imported state is therefore observed
but never advanced by a continuity solve.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
import subprocess
from collections import defaultdict
from pathlib import Path


Q_C = 1.602176634e-19


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def write_rows(path: Path, values: list[dict[str, object]]) -> None:
    if not values:
        raise RuntimeError(f"refusing to write empty CSV: {path}")
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(values[0]))
        writer.writeheader()
        writer.writerows(values)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def field(case: Path, name: str) -> dict[int, tuple[float, ...]]:
    path = case / "fields" / f"{name}_region0.csv"
    values: dict[int, tuple[float, ...]] = {}
    for row in read_rows(path):
        components = tuple(
            float(row[key])
            for key in ("component0", "component1", "component2")
            if key in row and row[key] != ""
        )
        values[int(row["node_id"])] = components
    return values


def scalar(values: dict[int, tuple[float, ...]], node: int) -> float:
    return values[node][0]


def endpoint_mean(
    values: dict[int, tuple[float, ...]], node0: int, node1: int
) -> tuple[float, ...]:
    a = values[node0]
    b = values[node1]
    return tuple(0.5 * (x + y) for x, y in zip(a, b))


def norm(values: tuple[float, ...]) -> float:
    return math.sqrt(sum(value * value for value in values))


def logarithmic_mean(a: float, b: float) -> float:
    if a <= 0.0 or b <= 0.0:
        return 0.0
    if abs(a - b) <= 1.0e-12 * max(a, b):
        return 0.5 * (a + b)
    return (a - b) / math.log(a / b)


def has_nonzero_joined_values(value: str) -> bool:
    return any(float(token) != 0.0 for token in value.split(";") if token)


def relative_l2(pairs: list[tuple[float, float]]) -> float:
    denominator = sum(reference * reference for _, reference in pairs)
    numerator = sum((candidate - reference) ** 2 for candidate, reference in pairs)
    if denominator == 0.0:
        return 0.0 if numerator == 0.0 else math.inf
    return math.sqrt(numerator / denominator)


def median_relative_error(pairs: list[tuple[float, float]]) -> float:
    peak = max((abs(reference) for _, reference in pairs), default=0.0)
    floor = max(peak * 1.0e-12, 1.0e-300)
    selected = [
        abs(candidate - reference) / abs(reference)
        for candidate, reference in pairs
        if abs(reference) > floor
    ]
    return statistics.median(selected) if selected else 0.0


def peak_ratio(pairs: list[tuple[float, float]]) -> float:
    reference = max((abs(value) for _, value in pairs), default=0.0)
    candidate = max((abs(value) for value, _ in pairs), default=0.0)
    if reference == 0.0:
        return 1.0 if candidate == 0.0 else math.inf
    return candidate / reference


def stage_row(
    bias: float,
    stage_index: int,
    stage: str,
    quantity: str,
    pairs: list[tuple[float, float]],
    note: str,
) -> dict[str, object]:
    return {
        "bias_V": bias,
        "stage_index": stage_index,
        "stage": stage,
        "quantity": quantity,
        "record_count": len(pairs),
        "relative_l2": relative_l2(pairs),
        "median_relative_error": median_relative_error(pairs),
        "peak_ratio": peak_ratio(pairs),
        "note": note,
    }


def triangle_gradient(
    coordinates_um: dict[int, tuple[float, float]],
    nodes: list[int],
    values: dict[int, float],
) -> tuple[float, float]:
    n0, n1, n2 = nodes
    x0, y0 = coordinates_um[n0]
    x1, y1 = coordinates_um[n1]
    x2, y2 = coordinates_um[n2]
    determinant = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)
    if determinant == 0.0:
        raise RuntimeError(f"degenerate triangle {nodes}")
    v0, v1, v2 = values[n0], values[n1], values[n2]
    gx_per_um = (
        v0 * (y1 - y2) + v1 * (y2 - y0) + v2 * (y0 - y1)
    ) / determinant
    gy_per_um = (
        v0 * (x2 - x1) + v1 * (x0 - x2) + v2 * (x1 - x0)
    ) / determinant
    return gx_per_um * 1.0e6, gy_per_um * 1.0e6


def parse_case(value: str) -> tuple[float, Path]:
    label, separator, raw_path = value.partition("=")
    if not separator:
        raise argparse.ArgumentTypeError("--case must be BIAS=PATH")
    try:
        bias = float(label)
    except ValueError as error:
        raise argparse.ArgumentTypeError(f"invalid case bias: {label}") from error
    path = Path(raw_path)
    if not path.is_dir():
        raise argparse.ArgumentTypeError(f"case directory does not exist: {path}")
    return bias, path


def verify_coordinates(
    mesh: dict, case: Path
) -> tuple[dict[int, tuple[float, float]], float]:
    coordinates = {
        int(node["id"]): (float(node["x"]), float(node["y"]))
        for node in mesh["nodes"]
    }
    sentaurus = {
        int(row["id"]): (float(row["x_um"]), float(row["y_um"]))
        for row in read_rows(case / "nodes.csv")
    }
    if coordinates.keys() != sentaurus.keys():
        raise RuntimeError(
            f"{case}: Vela/Sentaurus node ID sets do not match "
            f"({len(coordinates)} versus {len(sentaurus)})"
        )
    maximum = max(
        math.hypot(
            coordinates[node][0] - sentaurus[node][0],
            coordinates[node][1] - sentaurus[node][1],
        )
        for node in coordinates
    )
    if maximum > 1.0e-12:
        raise RuntimeError(
            f"{case}: maximum coordinate mismatch is {maximum:.17g} um"
        )
    return coordinates, maximum


def write_state(case: Path, path: Path) -> dict[str, dict[int, float]]:
    sources = {
        "psi_V": field(case, "ElectrostaticPotential"),
        "phin_V": field(case, "eQuasiFermiPotential"),
        "phip_V": field(case, "hQuasiFermiPotential"),
        "n_m3": field(case, "eDensity"),
        "p_m3": field(case, "hDensity"),
    }
    node_ids = set(sources["psi_V"])
    if any(set(values) != node_ids for values in sources.values()):
        raise RuntimeError(f"{case}: imported state fields use different node sets")
    state = {
        name: {
            node: scalar(values, node) * (1.0e6 if name in {"n_m3", "p_m3"} else 1.0)
            for node in node_ids
        }
        for name, values in sources.items()
    }
    with path.open("w", newline="", encoding="utf-8") as stream:
        columns = ["node_id", "psi_V", "phin_V", "phip_V", "n_m3", "p_m3"]
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        for node in sorted(node_ids):
            writer.writerow(
                {"node_id": node, **{name: state[name][node] for name in columns[1:]}}
            )
    return state


def mesh_integral_per_m_s(
    mesh: dict,
    coordinates_um: dict[int, tuple[float, float]],
    generation_cm3_s: dict[int, tuple[float, ...]],
) -> float:
    total = 0.0
    for triangle in mesh["triangles"]:
        nodes = [int(node) for node in triangle["node_ids"]]
        a, b, c = (coordinates_um[node] for node in nodes)
        area_m2 = abs(
            (b[0] - a[0]) * (c[1] - a[1])
            - (b[1] - a[1]) * (c[0] - a[0])
        ) * 0.5e-12
        average_m3_s = (
            sum(scalar(generation_cm3_s, node) for node in nodes) / 3.0
        ) * 1.0e6
        total += average_m3_s * area_m2
    return total


def run_case(
    bias: float,
    case: Path,
    args: argparse.Namespace,
    mesh: dict,
    baseline: dict,
) -> tuple[
    list[dict[str, object]],
    list[dict[str, object]],
    list[dict[str, object]],
    list[dict[str, object]],
    list[dict[str, object]],
    dict[str, object],
]:
    tag = ("m" if bias < 0 else "p") + str(abs(bias)).replace(".", "p")
    output = args.out_dir / tag
    output.mkdir(parents=True, exist_ok=True)
    coordinates, coordinate_error = verify_coordinates(mesh, case)
    state_path = output / "sentaurus_frozen_state.csv"
    imported_state = write_state(case, state_path)

    audit_config = json.loads(json.dumps(baseline))
    impact = audit_config["solver"]["impact_ionization"]
    impact["coupling_mode"] = "postprocess_only"
    if impact.get("model", "none") == "none":
        raise RuntimeError("baseline impact-ionization model must not be none")
    config_path = output / "audit_postprocess_only.json"
    config_path.write_text(
        json.dumps(audit_config, indent=2), encoding="utf-8"
    )

    paths = {
        "node": output / "vela_node_state.csv",
        "edge": output / "vela_edge_audit.csv",
        "triangle": output / "vela_triangle_audit.csv",
        "element": output / "vela_element_edge_gss_laux.csv",
        "process": output / "vela_bv_process_probe.csv",
    }
    command = [
        str(args.audit.resolve()),
        "--mesh", str(args.mesh.resolve()),
        "--doping", str(args.doping.resolve()),
        "--state", str(state_path.resolve()),
        "--config", str(config_path.resolve()),
        "--node-out", str(paths["node"].resolve()),
        "--edge-out", str(paths["edge"].resolve()),
        "--triangle-out", str(paths["triangle"].resolve()),
        "--element-out", str(paths["element"].resolve()),
        "--process-out", str(paths["process"].resolve()),
        "--scope", "general_tri3",
    ]
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    (output / "command.json").write_text(
        json.dumps(command, indent=2), encoding="utf-8"
    )
    (output / "stdout.log").write_text(completed.stdout, encoding="utf-8")
    (output / "stderr.log").write_text(completed.stderr, encoding="utf-8")
    if completed.returncode:
        raise RuntimeError(
            f"{bias:g} V fixed-state audit failed: {completed.stderr[-2000:]}"
        )

    sg_vector_config = json.loads(json.dumps(audit_config))
    sg_vector_impact = sg_vector_config["solver"]["impact_ionization"]
    sg_vector_impact["current_approximation"] = "element_edge_sg_gss_laux"
    sg_vector_impact["source_mapping_mode"] = "element_vertex_box_measure"
    sg_vector_impact["cell_reconstructed_midpoint_density"] = "bernoulli"
    sg_vector_config_path = output / "audit_sg_vector_postprocess_only.json"
    sg_vector_config_path.write_text(
        json.dumps(sg_vector_config, indent=2), encoding="utf-8"
    )
    sg_vector_paths = {
        "sg_vector_node": output / "sg_vector_node_state.csv",
        "sg_vector_edge": output / "sg_vector_edge_audit.csv",
        "sg_vector_triangle": output / "sg_vector_triangle_audit.csv",
        "sg_vector_element": output / "sg_vector_element_edge_gss_laux.csv",
        "sg_vector_process": output / "sg_vector_bv_process_probe.csv",
    }
    sg_vector_command = [
        str(args.audit.resolve()),
        "--mesh", str(args.mesh.resolve()),
        "--doping", str(args.doping.resolve()),
        "--state", str(state_path.resolve()),
        "--config", str(sg_vector_config_path.resolve()),
        "--node-out", str(sg_vector_paths["sg_vector_node"].resolve()),
        "--edge-out", str(sg_vector_paths["sg_vector_edge"].resolve()),
        "--triangle-out", str(sg_vector_paths["sg_vector_triangle"].resolve()),
        "--element-out", str(sg_vector_paths["sg_vector_element"].resolve()),
        "--process-out", str(sg_vector_paths["sg_vector_process"].resolve()),
        "--scope", "general_tri3",
    ]
    sg_vector_completed = subprocess.run(
        sg_vector_command, text=True, capture_output=True, check=False
    )
    (output / "sg_vector_command.json").write_text(
        json.dumps(sg_vector_command, indent=2), encoding="utf-8"
    )
    (output / "sg_vector_stdout.log").write_text(
        sg_vector_completed.stdout, encoding="utf-8"
    )
    (output / "sg_vector_stderr.log").write_text(
        sg_vector_completed.stderr, encoding="utf-8"
    )
    if sg_vector_completed.returncode:
        raise RuntimeError(
            f"{bias:g} V SG-vector audit failed: "
            f"{sg_vector_completed.stderr[-2000:]}"
        )
    paths.update(sg_vector_paths)

    summaries: list[dict[str, object]] = []
    details: list[dict[str, object]] = []
    roundtrip = read_rows(paths["node"])
    for quantity in ("psi_V", "phin_V", "phip_V", "n_m3", "p_m3"):
        pairs = [
            (float(row[quantity]), imported_state[quantity][int(row["node_id"])])
            for row in roundtrip
        ]
        summaries.append(
            stage_row(
                bias, 1, "state_import", quantity, pairs,
                "Vela fixed-state CSV versus the exact imported Sentaurus state",
            )
        )

    psi = imported_state["psi_V"]
    phin = imported_state["phin_V"]
    phip = imported_state["phip_V"]
    triangles = {
        int(triangle["id"]): [int(node) for node in triangle["node_ids"]]
        for triangle in mesh["triangles"]
    }
    gradients = {
        ("psi", cell): triangle_gradient(coordinates, nodes, psi)
        for cell, nodes in triangles.items()
    }
    gradients.update({
        ("electron", cell): triangle_gradient(coordinates, nodes, phin)
        for cell, nodes in triangles.items()
    })
    gradients.update({
        ("hole", cell): triangle_gradient(coordinates, nodes, phip)
        for cell, nodes in triangles.items()
    })

    sentaurus_fields = {
        "electric": field(case, "ElectricField"),
        "electron_current": field(case, "eCurrentDensity"),
        "hole_current": field(case, "hCurrentDensity"),
        "electron_alpha": field(case, "eAlphaAvalanche"),
        "hole_alpha": field(case, "hAlphaAvalanche"),
        "electron_mobility": field(case, "eMobility"),
        "hole_mobility": field(case, "hMobility"),
        "generation": field(case, "ImpactIonization"),
    }
    sentaurus_export_closure_pairs = []
    for node in sorted(sentaurus_fields["generation"]):
        electron_current_A_cm2 = norm(
            sentaurus_fields["electron_current"][node]
        )
        hole_current_A_cm2 = norm(
            sentaurus_fields["hole_current"][node]
        )
        reconstructed_cm3_s = (
            scalar(sentaurus_fields["electron_alpha"], node)
            * electron_current_A_cm2
            + scalar(sentaurus_fields["hole_alpha"], node)
            * hole_current_A_cm2
        ) / Q_C
        sentaurus_export_closure_pairs.append((
            reconstructed_cm3_s,
            scalar(sentaurus_fields["generation"], node),
        ))
    summaries.append(
        stage_row(
            bias, 0, "sentaurus_export_closure",
            "alpha_times_current_over_q",
            sentaurus_export_closure_pairs,
            "Sentaurus nodal alpha and vector-current magnitudes versus "
            "exported ImpactIonization",
        )
    )

    process = read_rows(paths["process"])
    edge_audit = {
        int(row["edge_id"]): row for row in read_rows(paths["edge"])
    }
    source_closure_errors = []
    qg_closure_errors = []
    for row in process:
        measure = float(row["source_measure_m2"])
        generation = float(row["generation_rate_per_m3_s"])
        source = float(row["source_integral_per_m_s"])
        qg = float(row["qG_contribution_A_per_m"])
        reconstructed_source = measure * generation
        source_scale = max(abs(source), abs(reconstructed_source), 1.0e-300)
        source_closure_errors.append(
            abs(source - reconstructed_source) / source_scale
        )
        reconstructed_qg = Q_C * source
        qg_scale = max(abs(qg), abs(reconstructed_qg), 1.0e-300)
        qg_closure_errors.append(abs(qg - reconstructed_qg) / qg_scale)
    if max(source_closure_errors, default=0.0) > 1.0e-12:
        raise RuntimeError(f"{bias:g} V process generation/source closure failed")
    if max(qg_closure_errors, default=0.0) > 1.0e-12:
        raise RuntimeError(f"{bias:g} V process source/qG closure failed")

    sg_vector_process = read_rows(paths["sg_vector_process"])
    if any(int(row["solver_coupled"]) != 0 for row in sg_vector_process):
        raise RuntimeError(f"{bias:g} V SG-vector control is not postprocess-only")
    sg_vector_details: list[dict[str, object]] = []
    sg_vector_pairs: dict[str, list[tuple[float, float]]] = defaultdict(list)
    for row in sg_vector_process:
        if row["support_kind"] != "element_vertex_gss_laux":
            continue
        carrier = row["carrier"]
        node = int(row["node0"])
        vela_current_A_m2 = Q_C * math.hypot(
            float(row["current_vector_x_per_m2_s"]),
            float(row["current_vector_y_per_m2_s"]),
        )
        sentaurus_current_A_m2 = (
            norm(sentaurus_fields[f"{carrier}_current"][node]) * 1.0e4
        )
        sg_vector_pairs[carrier].append(
            (vela_current_A_m2, sentaurus_current_A_m2)
        )
        sg_vector_details.append({
            "bias_V": bias,
            "carrier": carrier,
            "cell_id": int(row["cell_id"]),
            "node_id": node,
            "source_measure_m2": float(row["source_measure_m2"]),
            "vela_sg_vector_current_A_per_m2": vela_current_A_m2,
            "sentaurus_node_current_A_per_m2": sentaurus_current_A_m2,
            "vela_alpha_per_m": float(row["alpha_per_m"]),
            "vela_generation_rate_per_m3_s": float(
                row["generation_rate_per_m3_s"]
            ),
            "vela_source_integral_per_m_s": float(
                row["source_integral_per_m_s"]
            ),
        })
    for carrier, pairs in sorted(sg_vector_pairs.items()):
        summaries.append(
            stage_row(
                bias, 5, "actual_sg_vector_current",
                carrier, pairs,
                "element-edge GSS/Laux reconstructed vector versus "
                "Sentaurus nodal vector-current magnitude",
            )
        )
    sg_vector_total_source = sum(
        float(row["source_integral_per_m_s"]) for row in sg_vector_process
    )

    metric_pairs: dict[str, list[tuple[float, float]]] = defaultdict(list)
    generation_groups: dict[
        tuple[str, int, int, int, int, int], list[dict[str, str]]
    ] = defaultdict(list)
    external_by_record: dict[
        tuple[str, int, int, int, int, int, str], dict[str, float]
    ] = {}
    factorization_rows: list[dict[str, object]] = []
    factorization_pairs: dict[
        tuple[str, str], list[tuple[float, float]]
    ] = defaultdict(list)
    factorization_source_totals: dict[str, float] = defaultdict(float)

    for row in process:
        carrier = row["carrier"]
        node0, node1 = int(row["node0"]), int(row["node1"])
        cell_id = int(row["cell_id"])
        x0, y0 = coordinates[node0]
        x1, y1 = coordinates[node1]
        edge_length_um = math.hypot(x1 - x0, y1 - y0)
        if edge_length_um == 0.0:
            raise RuntimeError(f"zero-length process support {node0}->{node1}")
        edge_unit = ((x1 - x0) / edge_length_um, (y1 - y0) / edge_length_um)
        edge_field = -(
            psi[node1] - psi[node0]
        ) / (edge_length_um * 1.0e-6)
        electric_expected = tuple(edge_field * value for value in edge_unit)
        qf_expected = gradients[(carrier, cell_id)]
        vela_e = (
            float(row["electric_field_x_V_per_m"]),
            float(row["electric_field_y_V_per_m"]),
        )
        vela_qf = (
            float(row["qf_gradient_x_V_per_m"]),
            float(row["qf_gradient_y_V_per_m"]),
        )
        sent_e = tuple(
            value * 1.0e2
            for value in endpoint_mean(
                sentaurus_fields["electric"], node0, node1
            )
        )
        sent_e_edge_projection = abs(
            sent_e[0] * edge_unit[0] + sent_e[1] * edge_unit[1]
        )
        sent_current = tuple(
            value * 1.0e4
            for value in endpoint_mean(
                sentaurus_fields[f"{carrier}_current"], node0, node1
            )
        )
        sent_alpha = endpoint_mean(
            sentaurus_fields[f"{carrier}_alpha"], node0, node1
        )[0] * 1.0e2
        sent_mobility = endpoint_mean(
            sentaurus_fields[f"{carrier}_mobility"], node0, node1
        )[0] * 1.0e-4
        vela_current_A_m2 = Q_C * math.hypot(
            float(row["current_vector_x_per_m2_s"]),
            float(row["current_vector_y_per_m2_s"]),
        )
        sent_current_A_m2 = norm(sent_current)
        sent_current_edge_projection_A_m2 = abs(
            sent_current[0] * edge_unit[0]
            + sent_current[1] * edge_unit[1]
        )
        vela_alpha_per_m = float(row["alpha_per_m"])
        external_generation_vela_alpha = (
            vela_alpha_per_m * sent_current_A_m2 / Q_C
        )
        external_generation_sentaurus_alpha = (
            sent_alpha * sent_current_A_m2 / Q_C
        )
        mobility = float(row["final_mobility_m2_per_V_s"])
        edge_qf_drive = float(row["high_field_drive_V_per_m"])
        midpoint = float(row["midpoint_density_m3"])
        density0 = float(row["density0_m3"])
        density1 = float(row["density1_m3"])
        arithmetic_midpoint = 0.5 * (density0 + density1)
        geometric_midpoint = math.sqrt(max(density0, 0.0) * max(density1, 0.0))
        log_midpoint = logarithmic_mean(density0, density1)
        reconstructed_proxy_A_m2 = (
            Q_C * mobility * midpoint * edge_qf_drive
        )
        arithmetic_proxy_A_m2 = (
            Q_C * mobility * arithmetic_midpoint * edge_qf_drive
        )
        geometric_proxy_A_m2 = (
            Q_C * mobility * geometric_midpoint * edge_qf_drive
        )
        log_proxy_A_m2 = Q_C * mobility * log_midpoint * edge_qf_drive
        sentaurus_mobility_proxy_A_m2 = (
            Q_C * sent_mobility * midpoint * edge_qf_drive
        )
        raw_flux_column = f"{carrier}_raw_signed_flux_per_m2_s"
        edge_record = edge_audit[int(row["edge_id"])]
        raw_sg_current_A_m2 = Q_C * abs(float(edge_record[raw_flux_column]))
        edge_midpoint_density = float(
            edge_record[f"{carrier}_midpoint_density_m3"]
        )
        edge_midpoint_proxy_A_m2 = (
            Q_C * mobility * edge_midpoint_density * edge_qf_drive
        )
        low_endpoint_density = min(density0, density1)
        high_endpoint_density = max(density0, density1)
        low_endpoint_proxy_A_m2 = (
            Q_C * mobility * low_endpoint_density * edge_qf_drive
        )
        high_endpoint_proxy_A_m2 = (
            Q_C * mobility * high_endpoint_density * edge_qf_drive
        )
        implied_midpoint_from_sentaurus = (
            sent_current_A_m2 / (Q_C * mobility * edge_qf_drive)
            if mobility > 0.0 and edge_qf_drive > 0.0 else 0.0
        )
        implied_drive_from_sentaurus = (
            sent_current_A_m2 / (Q_C * mobility * midpoint)
            if mobility > 0.0 and midpoint > 0.0 else 0.0
        )
        values = {
            "state_to_electric_field": (norm(vela_e), norm(electric_expected)),
            f"state_to_qf_gradient_{carrier}": (norm(vela_qf), norm(qf_expected)),
            "sentaurus_electric_field": (
                norm(vela_e), sent_e_edge_projection
            ),
            f"sentaurus_mobility_{carrier}": (
                float(row["final_mobility_m2_per_V_s"]), sent_mobility
            ),
            f"sentaurus_current_{carrier}": (
                vela_current_A_m2, sent_current_A_m2
            ),
            f"sentaurus_alpha_{carrier}": (
                float(row["alpha_per_m"]), sent_alpha
            ),
        }
        for metric, pair in values.items():
            metric_pairs[metric].append(pair)
        details.append({
            "bias_V": bias,
            "support_kind": row["support_kind"],
            "carrier": carrier,
            "cell_id": cell_id,
            "local_edge": int(row["local_edge"]),
            "edge_id": int(row["edge_id"]),
            "node0": node0,
            "node1": node1,
            "vela_electric_field_V_per_m": norm(vela_e),
            "sentaurus_electric_field_endpoint_mean_V_per_m": norm(sent_e),
            "sentaurus_electric_field_edge_projection_V_per_m": (
                sent_e_edge_projection
            ),
            "vela_qf_gradient_V_per_m": norm(vela_qf),
            "vela_final_mobility_m2_per_V_s": values[
                f"sentaurus_mobility_{carrier}"
            ][0],
            "sentaurus_mobility_endpoint_mean_m2_per_V_s": sent_mobility,
            "vela_current_magnitude_A_per_m2": vela_current_A_m2,
            "sentaurus_current_endpoint_mean_magnitude_A_per_m2": (
                sent_current_A_m2
            ),
            "vela_alpha_per_m": vela_alpha_per_m,
            "sentaurus_alpha_endpoint_mean_per_m": sent_alpha,
            "vela_generation_rate_per_m3_s": float(
                row["generation_rate_per_m3_s"]
            ),
            "vela_source_integral_per_m_s": float(
                row["source_integral_per_m_s"]
            ),
            "external_current_vela_alpha_generation_rate_per_m3_s": (
                external_generation_vela_alpha
            ),
            "external_current_sentaurus_alpha_generation_rate_per_m3_s": (
                external_generation_sentaurus_alpha
            ),
            "solver_coupled": int(row["solver_coupled"]),
        })
        group_key = (
            row["support_kind"], cell_id, int(row["local_edge"]),
            int(row["edge_id"]), node0, node1,
        )
        generation_groups[group_key].append(row)
        external_by_record[(*group_key, carrier)] = {
            "sentaurus_current_A_m2": sent_current_A_m2,
            "vela_alpha_per_m": vela_alpha_per_m,
            "sentaurus_alpha_per_m": sent_alpha,
            "vela_alpha_generation_m3_s": external_generation_vela_alpha,
            "sentaurus_alpha_generation_m3_s": (
                external_generation_sentaurus_alpha
            ),
        }
        positive_source_support = (
            float(row["source_measure_m2"]) > 0.0
        )
        factorization_rows.append({
            "bias_V": bias,
            "support_kind": row["support_kind"],
            "carrier": carrier,
            "cell_id": cell_id,
            "local_edge": int(row["local_edge"]),
            "edge_id": int(row["edge_id"]),
            "node0": node0,
            "node1": node1,
            "positive_source_support": int(positive_source_support),
            "mobility_m2_per_V_s": mobility,
            "sentaurus_mobility_endpoint_mean_m2_per_V_s": sent_mobility,
            "density0_m3": density0,
            "density1_m3": density1,
            "gss_midpoint_density_m3": midpoint,
            "arithmetic_midpoint_density_m3": arithmetic_midpoint,
            "geometric_midpoint_density_m3": geometric_midpoint,
            "logarithmic_midpoint_density_m3": log_midpoint,
            "vela_edge_audit_midpoint_density_m3": edge_midpoint_density,
            "low_endpoint_density_m3": low_endpoint_density,
            "high_endpoint_density_m3": high_endpoint_density,
            "sentaurus_current_implied_midpoint_density_m3": (
                implied_midpoint_from_sentaurus
            ),
            "edge_qf_drive_V_per_m": edge_qf_drive,
            "cell_qf_drive_V_per_m": float(row["impact_field_V_per_m"]),
            "sentaurus_current_implied_qf_drive_V_per_m": (
                implied_drive_from_sentaurus
            ),
            "production_proxy_current_A_per_m2": vela_current_A_m2,
            "factor_product_reconstructed_proxy_current_A_per_m2": (
                reconstructed_proxy_A_m2
            ),
            "sentaurus_mobility_proxy_current_A_per_m2": (
                sentaurus_mobility_proxy_A_m2
            ),
            "arithmetic_midpoint_proxy_current_A_per_m2": (
                arithmetic_proxy_A_m2
            ),
            "geometric_midpoint_proxy_current_A_per_m2": (
                geometric_proxy_A_m2
            ),
            "logarithmic_midpoint_proxy_current_A_per_m2": (
                log_proxy_A_m2
            ),
            "vela_edge_midpoint_proxy_current_A_per_m2": (
                edge_midpoint_proxy_A_m2
            ),
            "low_endpoint_proxy_current_A_per_m2": (
                low_endpoint_proxy_A_m2
            ),
            "high_endpoint_proxy_current_A_per_m2": (
                high_endpoint_proxy_A_m2
            ),
            "vela_raw_sg_transport_current_A_per_m2": raw_sg_current_A_m2,
            "sentaurus_current_magnitude_A_per_m2": sent_current_A_m2,
            "sentaurus_current_edge_projection_A_per_m2": (
                sent_current_edge_projection_A_m2
            ),
            "proxy_over_raw_sg": (
                vela_current_A_m2 / raw_sg_current_A_m2
                if raw_sg_current_A_m2 > 0.0 else math.inf
            ),
            "gss_midpoint_over_arithmetic": (
                midpoint / arithmetic_midpoint
                if arithmetic_midpoint > 0.0 else math.inf
            ),
            "gss_midpoint_over_geometric": (
                midpoint / geometric_midpoint
                if geometric_midpoint > 0.0 else math.inf
            ),
            "gss_midpoint_over_logarithmic": (
                midpoint / log_midpoint
                if log_midpoint > 0.0 else math.inf
            ),
            "triangle_gss_midpoint_over_edge_midpoint": (
                midpoint / edge_midpoint_density
                if edge_midpoint_density > 0.0 else math.inf
            ),
        })
        if positive_source_support:
            pairs = {
                "production_proxy_vs_sentaurus_magnitude": (
                    vela_current_A_m2, sent_current_A_m2
                ),
                "factor_product_closure": (
                    reconstructed_proxy_A_m2, vela_current_A_m2
                ),
                "sentaurus_mobility_proxy_vs_sentaurus_magnitude": (
                    sentaurus_mobility_proxy_A_m2, sent_current_A_m2
                ),
                "arithmetic_midpoint_proxy_vs_sentaurus_magnitude": (
                    arithmetic_proxy_A_m2, sent_current_A_m2
                ),
                "geometric_midpoint_proxy_vs_sentaurus_magnitude": (
                    geometric_proxy_A_m2, sent_current_A_m2
                ),
                "logarithmic_midpoint_proxy_vs_sentaurus_magnitude": (
                    log_proxy_A_m2, sent_current_A_m2
                ),
                "edge_midpoint_proxy_vs_sentaurus_magnitude": (
                    edge_midpoint_proxy_A_m2, sent_current_A_m2
                ),
                "low_endpoint_proxy_vs_sentaurus_magnitude": (
                    low_endpoint_proxy_A_m2, sent_current_A_m2
                ),
                "high_endpoint_proxy_vs_sentaurus_magnitude": (
                    high_endpoint_proxy_A_m2, sent_current_A_m2
                ),
                "raw_sg_vs_sentaurus_edge_projection": (
                    raw_sg_current_A_m2,
                    sent_current_edge_projection_A_m2,
                ),
                "production_proxy_vs_raw_sg": (
                    vela_current_A_m2, raw_sg_current_A_m2
                ),
            }
            for quantity, pair in pairs.items():
                factorization_pairs[(carrier, quantity)].append(pair)
            source_measure = float(row["source_measure_m2"])
            current_candidates = {
                "production_proxy": vela_current_A_m2,
                "sentaurus_mobility_proxy": sentaurus_mobility_proxy_A_m2,
                "arithmetic_midpoint_proxy": arithmetic_proxy_A_m2,
                "geometric_midpoint_proxy": geometric_proxy_A_m2,
                "logarithmic_midpoint_proxy": log_proxy_A_m2,
                "edge_midpoint_proxy": edge_midpoint_proxy_A_m2,
                "low_endpoint_proxy": low_endpoint_proxy_A_m2,
                "high_endpoint_proxy": high_endpoint_proxy_A_m2,
                "raw_sg_transport": raw_sg_current_A_m2,
                "sentaurus_current": sent_current_A_m2,
            }
            for candidate, current_A_m2 in current_candidates.items():
                factorization_source_totals[candidate] += (
                    vela_alpha_per_m * current_A_m2 / Q_C * source_measure
                )

    operator_specs = [
        ("state_to_electric_field", 2, "state_to_driver", "electric_field"),
        ("state_to_qf_gradient_electron", 2, "state_to_driver", "electron_qf_gradient"),
        ("state_to_qf_gradient_hole", 2, "state_to_driver", "hole_qf_gradient"),
        (
            "sentaurus_electric_field", 3, "auxiliary_field",
            "electric_field_edge_projection",
        ),
        ("sentaurus_mobility_electron", 4, "mobility", "electron_mobility"),
        ("sentaurus_mobility_hole", 4, "mobility", "hole_mobility"),
        ("sentaurus_current_electron", 5, "current", "electron_current"),
        ("sentaurus_current_hole", 5, "current", "hole_current"),
        ("sentaurus_alpha_electron", 6, "ionization_coefficient", "electron_alpha"),
        ("sentaurus_alpha_hole", 6, "ionization_coefficient", "hole_alpha"),
    ]
    for key, index, stage, quantity in operator_specs:
        note = (
            "independent reconstruction from imported state"
            if stage == "state_to_driver"
            else (
                "Sentaurus nodal endpoint mean projected to the Vela edge; "
                "the active impact-ionization drive is quasi-Fermi gradient"
                if stage == "auxiliary_field"
                else "Sentaurus nodal endpoint mean projected to Vela production support"
            )
        )
        summaries.append(
            stage_row(bias, index, stage, quantity, metric_pairs[key], note)
        )

    local_generation_pairs = []
    external_current_local_pairs = []
    sentaurus_model_local_pairs = []
    external_support_rows: list[dict[str, object]] = []
    external_current_total = 0.0
    sentaurus_model_total = 0.0
    sentaurus_geometry_total = 0.0
    for key, records in generation_groups.items():
        _, _, _, _, node0, node1 = key
        positive = [row for row in records if float(row["source_measure_m2"]) > 0.0]
        if not positive:
            continue
        source_measure = float(positive[0]["source_measure_m2"])
        if any(
            abs(float(row["source_measure_m2"]) - source_measure)
            > max(source_measure, 1.0) * 1.0e-14
            for row in positive
        ):
            raise RuntimeError(f"inconsistent carrier source measure for {key}")
        vela_generation = sum(
            float(row["generation_rate_per_m3_s"]) for row in positive
        )
        sentaurus_generation = endpoint_mean(
            sentaurus_fields["generation"], node0, node1
        )[0] * 1.0e6
        local_generation_pairs.append((vela_generation, sentaurus_generation))
        external_generation = sum(
            external_by_record[(*key, row["carrier"])][
                "vela_alpha_generation_m3_s"
            ]
            for row in positive
        )
        sentaurus_model_generation = sum(
            external_by_record[(*key, row["carrier"])][
                "sentaurus_alpha_generation_m3_s"
            ]
            for row in positive
        )
        external_source = external_generation * source_measure
        sentaurus_model_source = sentaurus_model_generation * source_measure
        sentaurus_geometry_source = sentaurus_generation * source_measure
        external_current_total += external_source
        sentaurus_model_total += sentaurus_model_source
        sentaurus_geometry_total += sentaurus_geometry_source
        external_current_local_pairs.append(
            (external_generation, sentaurus_generation)
        )
        sentaurus_model_local_pairs.append(
            (sentaurus_model_generation, sentaurus_generation)
        )
        external_support_rows.append({
            "bias_V": bias,
            "support_kind": key[0],
            "cell_id": key[1],
            "local_edge": key[2],
            "edge_id": key[3],
            "node0": node0,
            "node1": node1,
            "source_measure_m2": source_measure,
            "vela_baseline_generation_rate_per_m3_s": vela_generation,
            "vela_alpha_sentaurus_current_generation_rate_per_m3_s": (
                external_generation
            ),
            "sentaurus_alpha_current_generation_rate_per_m3_s": (
                sentaurus_model_generation
            ),
            "sentaurus_impact_ionization_endpoint_mean_per_m3_s": (
                sentaurus_generation
            ),
            "vela_baseline_source_integral_per_m_s": (
                vela_generation * source_measure
            ),
            "vela_alpha_sentaurus_current_source_integral_per_m_s": (
                external_source
            ),
            "sentaurus_alpha_current_vela_geometry_source_integral_per_m_s": (
                sentaurus_model_source
            ),
            "sentaurus_generation_vela_geometry_source_integral_per_m_s": (
                sentaurus_geometry_source
            ),
        })
    summaries.append(
        stage_row(
            bias, 7, "local_generation", "total_impact_ionization",
            local_generation_pairs,
            "carrier-summed Vela generation versus Sentaurus nodal endpoint mean",
        )
    )
    summaries.append(
        stage_row(
            bias, 7, "local_generation_external_current",
            "vela_alpha_sentaurus_current",
            external_current_local_pairs,
            "Vela alpha and production source support with only Jn/Jp "
            "replaced by Sentaurus vector-current endpoint means",
        )
    )
    summaries.append(
        stage_row(
            bias, 7, "local_generation_sentaurus_model",
            "sentaurus_alpha_current_vela_geometry",
            sentaurus_model_local_pairs,
            "Sentaurus alpha/current projected to Vela production support",
        )
    )

    vela_total = sum(float(row["source_integral_per_m_s"]) for row in process)
    sentaurus_total = mesh_integral_per_m_s(
        mesh, coordinates, sentaurus_fields["generation"]
    )
    summaries.append(
        stage_row(
            bias, 8, "integrated_source", "total_impact_ionization",
            [(vela_total, sentaurus_total)],
            "Vela production support sum versus P1 Sentaurus nodal area integral",
        )
    )
    summaries.append(
        stage_row(
            bias, 8, "integrated_source_external_current",
            "vela_alpha_sentaurus_current",
            [(external_current_total, sentaurus_total)],
            "Only current magnitude is substituted; Vela alpha and source "
            "geometry remain unchanged",
        )
    )
    summaries.append(
        stage_row(
            bias, 8, "integrated_source_sentaurus_model",
            "sentaurus_alpha_current_vela_geometry",
            [(sentaurus_model_total, sentaurus_total)],
            "Sentaurus alpha/current with Vela source geometry",
        )
    )
    summaries.append(
        stage_row(
            bias, 8, "integrated_source_geometry_control",
            "sentaurus_generation_vela_geometry",
            [(sentaurus_geometry_total, sentaurus_total)],
            "Sentaurus ImpactIonization projected directly to Vela source geometry",
        )
    )
    summaries.append(
        stage_row(
            bias, 8, "integrated_source_actual_sg_vector",
            "element_edge_gss_laux",
            [(sg_vector_total_source, sentaurus_total)],
            "Vela alpha and element-vertex geometry driven by reconstructed "
            "actual SG transport vectors",
        )
    )
    factorization_notes = {
        "production_proxy_vs_sentaurus_magnitude": (
            "production q*mu*n_mid*edge-QFP-drive versus Sentaurus |J|"
        ),
        "factor_product_closure": (
            "independent q*mu*n_mid*edge-QFP-drive reconstruction"
        ),
        "sentaurus_mobility_proxy_vs_sentaurus_magnitude": (
            "only mobility replaced by Sentaurus endpoint mean"
        ),
        "arithmetic_midpoint_proxy_vs_sentaurus_magnitude": (
            "only GSS midpoint replaced by arithmetic endpoint mean"
        ),
        "geometric_midpoint_proxy_vs_sentaurus_magnitude": (
            "only GSS midpoint replaced by geometric endpoint mean"
        ),
        "logarithmic_midpoint_proxy_vs_sentaurus_magnitude": (
            "only GSS midpoint replaced by logarithmic endpoint mean"
        ),
        "edge_midpoint_proxy_vs_sentaurus_magnitude": (
            "triangle proxy rebuilt with the fixed-state edge-audit midpoint"
        ),
        "low_endpoint_proxy_vs_sentaurus_magnitude": (
            "triangle proxy rebuilt with the lower endpoint density"
        ),
        "high_endpoint_proxy_vs_sentaurus_magnitude": (
            "triangle proxy rebuilt with the higher endpoint density"
        ),
        "raw_sg_vs_sentaurus_edge_projection": (
            "Vela raw SG transport current versus Sentaurus edge projection"
        ),
        "production_proxy_vs_raw_sg": (
            "production avalanche current proxy versus Vela raw SG transport"
        ),
    }
    for (carrier, quantity), pairs in sorted(factorization_pairs.items()):
        summaries.append(
            stage_row(
                bias, 5, "current_proxy_factorization",
                f"{carrier}_{quantity}", pairs,
                factorization_notes[quantity],
            )
        )
    for candidate, total in sorted(factorization_source_totals.items()):
        summaries.append(
            stage_row(
                bias, 8, "integrated_source_current_factorization",
                candidate, [(total, sentaurus_total)],
                "Vela alpha and production geometry with the named "
                "current-factorization candidate",
            )
        )

    if any(int(row["solver_coupled"]) != 0 for row in process):
        raise RuntimeError(f"{bias:g} V process output is not postprocess-only")
    if any(
        has_nonzero_joined_values(
            row["electron_residual_contributions_per_m_s"]
        )
        or has_nonzero_joined_values(
            row["hole_residual_contributions_per_m_s"]
        )
        for row in process
    ):
        raise RuntimeError(f"{bias:g} V process output has continuity residual feedback")

    finite_hotspots = [
        row for row in factorization_rows
        if row["positive_source_support"]
        and math.isfinite(float(row["proxy_over_raw_sg"]))
    ]
    hotspot = max(
        finite_hotspots,
        key=lambda row: float(row["proxy_over_raw_sg"]),
        default=None,
    )
    positive_source_integrals = sorted(
        (
            max(float(row["source_integral_per_m_s"]), 0.0)
            for row in process
        ),
        reverse=True,
    )
    top8_source_integral = sum(positive_source_integrals[:8])
    positive_source_integral = sum(positive_source_integrals)

    manifest = {
        "bias_V": bias,
        "sentaurus_case": str(case.resolve()),
        "maximum_coordinate_mismatch_um": coordinate_error,
        "state_node_count": len(imported_state["psi_V"]),
        "process_record_count": len(process),
        "impact_ionization_coupling_mode": "postprocess_only",
        "solver_coupled_record_count": sum(
            int(row["solver_coupled"]) for row in process
        ),
        "nonempty_residual_feedback_record_count": sum(
            has_nonzero_joined_values(
                row["electron_residual_contributions_per_m_s"]
            )
            or has_nonzero_joined_values(
                row["hole_residual_contributions_per_m_s"]
            )
            for row in process
        ),
        "maximum_generation_source_closure_relative_error": max(
            source_closure_errors, default=0.0
        ),
        "maximum_source_qG_closure_relative_error": max(
            qg_closure_errors, default=0.0
        ),
        "vela_total_source_integral_per_m_s": vela_total,
        "top8_positive_source_integral_per_m_s": top8_source_integral,
        "top8_positive_source_share": (
            top8_source_integral / positive_source_integral
            if positive_source_integral > 0.0
            else 0.0
        ),
        "vela_alpha_sentaurus_current_total_source_integral_per_m_s": (
            external_current_total
        ),
        "sentaurus_alpha_current_vela_geometry_total_source_integral_per_m_s": (
            sentaurus_model_total
        ),
        "sentaurus_generation_vela_geometry_total_source_integral_per_m_s": (
            sentaurus_geometry_total
        ),
        "sentaurus_total_source_integral_per_m_s": sentaurus_total,
        "sentaurus_export_alpha_current_closure_relative_l2": relative_l2(
            sentaurus_export_closure_pairs
        ),
        "current_proxy_factorization_hotspot": hotspot,
        "current_factorization_total_source_integrals_per_m_s": dict(
            sorted(factorization_source_totals.items())
        ),
        "sg_vector_total_source_integral_per_m_s": sg_vector_total_source,
        "sg_vector_process_record_count": len(sg_vector_process),
        "inputs": {
            "mesh": {"path": str(args.mesh.resolve()), "sha256": sha256(args.mesh)},
            "doping": {
                "path": str(args.doping.resolve()), "sha256": sha256(args.doping)
            },
            "baseline_config": {
                "path": str(args.baseline_config.resolve()),
                "sha256": sha256(args.baseline_config),
            },
            "sentaurus_field_manifest": {
                "path": str((case / "field_manifest.json").resolve()),
                "sha256": sha256(case / "field_manifest.json"),
            },
        },
        "outputs": {
            name: {"path": str(path.resolve()), "sha256": sha256(path)}
            for name, path in paths.items()
        },
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2), encoding="utf-8"
    )
    return (
        summaries,
        details,
        external_support_rows,
        factorization_rows,
        sg_vector_details,
        manifest,
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--doping", type=Path, required=True)
    parser.add_argument("--baseline-config", type=Path, required=True)
    parser.add_argument(
        "--case", action="append", type=parse_case, required=True,
        help="repeatable BIAS=PATH Sentaurus imported-state case",
    )
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    mesh = json.loads(args.mesh.read_text(encoding="utf-8-sig"))
    baseline = json.loads(args.baseline_config.read_text(encoding="utf-8-sig"))
    all_summaries: list[dict[str, object]] = []
    all_details: list[dict[str, object]] = []
    all_external_support: list[dict[str, object]] = []
    all_factorization: list[dict[str, object]] = []
    all_sg_vector: list[dict[str, object]] = []
    manifests = []
    for bias, case in args.case:
        (
            summaries,
            details,
            external_support,
            factorization,
            sg_vector,
            manifest,
        ) = run_case(
            bias, case, args, mesh, baseline
        )
        all_summaries.extend(summaries)
        all_details.extend(details)
        all_external_support.extend(external_support)
        all_factorization.extend(factorization)
        all_sg_vector.extend(sg_vector)
        manifests.append(manifest)

    stage_summary_path = args.out_dir / "stage_summary.csv"
    support_comparison_path = args.out_dir / "support_comparison.csv"
    external_substitution_path = (
        args.out_dir / "external_current_substitution.csv"
    )
    factorization_path = args.out_dir / "current_proxy_factorization.csv"
    sg_vector_path = args.out_dir / "sg_vector_current_control.csv"
    write_rows(stage_summary_path, all_summaries)
    write_rows(support_comparison_path, all_details)
    write_rows(
        external_substitution_path,
        all_external_support,
    )
    write_rows(factorization_path, all_factorization)
    write_rows(sg_vector_path, all_sg_vector)
    first_departures = {}
    for bias, _ in args.case:
        candidates = [
            row for row in all_summaries
            if float(row["bias_V"]) == bias
            and int(row["stage_index"]) >= 4
            and float(row["relative_l2"]) > 0.05
        ]
        first_departures[str(bias)] = (
            min(candidates, key=lambda row: int(row["stage_index"]))
            if candidates else None
        )
    result = {
        "schema": "vela.sentaurus_frozen_bv_operator_comparison.v1",
        "observation_only": True,
        "state_advanced": False,
        "continuity_feedback_enabled": False,
        "external_current_usage": "diagnostic_alpha_source_postprocess_only",
        "first_departure_threshold_relative_l2": 0.05,
        "first_departures": first_departures,
        "first_departure_scope": (
            "production mobility/current/alpha/source stages; auxiliary "
            "electric-field projection is excluded because the configured "
            "impact-ionization drive is quasi-Fermi gradient"
        ),
        "cases": manifests,
        "artifacts": {
            "stage_summary": {
                "path": str(stage_summary_path.resolve()),
                "sha256": sha256(stage_summary_path),
            },
            "support_comparison": {
                "path": str(support_comparison_path.resolve()),
                "sha256": sha256(support_comparison_path),
            },
            "external_current_substitution": {
                "path": str(external_substitution_path.resolve()),
                "sha256": sha256(external_substitution_path),
            },
            "current_proxy_factorization": {
                "path": str(factorization_path.resolve()),
                "sha256": sha256(factorization_path),
            },
            "sg_vector_current_control": {
                "path": str(sg_vector_path.resolve()),
                "sha256": sha256(sg_vector_path),
            },
        },
    }
    (args.out_dir / "result.json").write_text(
        json.dumps(result, indent=2), encoding="utf-8"
    )

    lines = [
        "# Sentaurus frozen-state Vela BV operator comparison",
        "",
        "The Sentaurus state was imported without modification. Vela evaluated "
        "the production impact-ionization chain with `postprocess_only`; no "
        "continuity residual or Jacobian feedback was enabled.",
        "",
        "| Bias (V) | Stage | Quantity | Relative L2 | Median relative error | Peak ratio |",
        "|---:|---|---|---:|---:|---:|",
    ]
    for row in all_summaries:
        lines.append(
            f"| {float(row['bias_V']):g} | {row['stage']} | {row['quantity']} | "
            f"{float(row['relative_l2']):.6g} | "
            f"{float(row['median_relative_error']):.6g} | "
            f"{float(row['peak_ratio']):.6g} |"
        )
    lines.extend([
        "",
        "Sentaurus nodal quantities are projected to each Vela production "
        "support by endpoint averaging. Integrated Sentaurus impact "
        "ionization uses a linear-triangle nodal area integral. These "
        "projection definitions are part of the comparison contract and "
        "should not be interpreted as an IIC path integral.",
        "",
    ])
    (args.out_dir / "report.md").write_text("\n".join(lines), encoding="utf-8")


if __name__ == "__main__":
    main()
