#!/usr/bin/env python3
"""Build a contract-valid Vela PN2D BV process manifest from exact-lattice runs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import shutil
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

try:
    from scripts.pn2d_bv_process_contract import validate_process_run
except ModuleNotFoundError:
    from pn2d_bv_process_contract import validate_process_run


BRANCHES = ("avalanche_off", "iic_postprocess", "avalanche_on")
ELEMENTARY_CHARGE_C = 1.602176634e-19
INTERNAL_LINE_SOURCE_TO_A_PER_UM = 1.0e-12
INTERNAL_PARTICLE_FLUX_TO_PER_CM2_S = 1.0e6
M3_TO_CM3 = 1.0e-6
UM_INV_TO_CM_INV = 1.0e4
EXACT_BIAS_TOLERANCE_V = 1.0e-10


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def bias_token(value: float) -> str:
    prefix = "m" if value < 0.0 else ""
    fixed = f"{abs(value):.6f}".replace(".", "p")
    return prefix + fixed


def relative(root: Path, path: Path) -> str:
    return path.resolve().relative_to(root.resolve()).as_posix()


def artifact(root: Path, path: Path) -> dict[str, str]:
    return {"path": relative(root, path), "sha256": sha256(path)}


def write_jsonl(path: Path, records: Iterable[dict[str, Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for record in records:
            handle.write(
                json.dumps(
                    record,
                    sort_keys=True,
                    separators=(",", ":"),
                    allow_nan=False,
                )
                + "\n"
            )


def load_mesh(path: Path) -> tuple[dict[int, tuple[float, float]], dict[int, list[int]]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    nodes = {
        int(node["id"]): (float(node["x"]), float(node["y"]))
        for node in payload["nodes"]
    }
    cells = {
        int(cell["id"]): [int(node) for node in cell["node_ids"]]
        for cell in payload["triangles"]
    }
    return nodes, cells


def triangle_gradient(
    node_ids: list[int],
    nodes: dict[int, tuple[float, float]],
    values: dict[int, float],
) -> tuple[float, float]:
    first, second, third = node_ids
    x0, y0 = nodes[first]
    x1, y1 = nodes[second]
    x2, y2 = nodes[third]
    area2 = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)
    if area2 == 0.0:
        raise ValueError(f"degenerate cell {node_ids}")
    gx = (
        values[first] * (y1 - y2)
        + values[second] * (y2 - y0)
        + values[third] * (y0 - y1)
    ) / area2
    gy = (
        values[first] * (x2 - x1)
        + values[second] * (x0 - x2)
        + values[third] * (x1 - x0)
    ) / area2
    return gx, gy


def nodal_control_areas_um2(
    nodes: dict[int, tuple[float, float]],
    cells: dict[int, list[int]],
) -> dict[int, float]:
    result = {node: 0.0 for node in nodes}
    for node_ids in cells.values():
        x0, y0 = nodes[node_ids[0]]
        x1, y1 = nodes[node_ids[1]]
        x2, y2 = nodes[node_ids[2]]
        area_um2 = 0.5 * abs(
            (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0)
        )
        for node in node_ids:
            result[node] += area_um2 / 3.0
    return result


def source(
    root: Path,
    path: Path,
    dataset: str,
    index: int,
) -> dict[str, Any]:
    return {
        "file": relative(root, path),
        "dataset": dataset,
        "index": index,
    }


def field(
    *,
    branch: str,
    requested: float,
    actual: float,
    support_kind: str,
    support_key: str,
    provenance: str,
    carrier: str,
    quantity: str,
    components: list[str],
    unit: str,
    values: list[float],
    source_record: dict[str, Any],
    coordinates_um: list[float] | None = None,
    connectivity: list[int] | None = None,
) -> dict[str, Any]:
    centering = {
        "physical_node": "vertex",
        "cell": "cell",
        "element_local_edge": "element_edge",
        "element_local_vertex": "element_vertex",
        "contact": "contact",
    }[support_kind]
    result: dict[str, Any] = {
        "branch": branch,
        "requested_bias_V": requested,
        "actual_bias_V": actual,
        "support_kind": support_kind,
        "support_key": support_key,
        "centering": centering,
        "provenance": provenance,
        "carrier": carrier,
        "quantity": quantity,
        "components": components,
        "unit": unit,
        "values": values,
        "source": source_record,
    }
    if coordinates_um is not None:
        result["coordinates_um"] = coordinates_um
    if connectivity is not None:
        result["connectivity"] = connectivity
    return result


def state_fields(
    root: Path,
    path: Path,
    branch: str,
    requested: float,
    actual: float,
    nodes: dict[int, tuple[float, float]],
) -> tuple[list[dict[str, Any]], dict[str, dict[int, float]]]:
    result: list[dict[str, Any]] = []
    state: dict[str, dict[int, float]] = {
        "psi": {},
        "phin": {},
        "phip": {},
        "electron": {},
        "hole": {},
    }
    for index, row in enumerate(rows(path)):
        node = int(row["node_id"])
        coordinates = list(nodes[node])
        state["psi"][node] = float(row["psi"])
        state["phin"][node] = float(row["phin"])
        state["phip"][node] = float(row["phip"])
        state["electron"][node] = float(row["electrons_m3"]) * M3_TO_CM3
        state["hole"][node] = float(row["holes_m3"]) * M3_TO_CM3
        record_source = source(root, path, "state", index)
        common = {
            "branch": branch,
            "requested": requested,
            "actual": actual,
            "support_kind": "physical_node",
            "support_key": f"node:{node}",
            "provenance": "solver_used",
            "source_record": record_source,
            "coordinates_um": coordinates,
        }
        result.append(
            field(
                **common,
                carrier="none",
                quantity="potential",
                components=["scalar"],
                unit="V",
                values=[state["psi"][node]],
            )
        )
        for carrier, qf_name in (("electron", "phin"), ("hole", "phip")):
            result.append(
                field(
                    **common,
                    carrier=carrier,
                    quantity="quasi_fermi",
                    components=["scalar"],
                    unit="V",
                    values=[state[qf_name][node]],
                )
            )
            result.append(
                field(
                    **common,
                    carrier=carrier,
                    quantity="density",
                    components=["scalar"],
                    unit="cm^-3",
                    values=[state[carrier][node]],
                )
            )
    return result, state


def probe_fields(
    root: Path,
    path: Path,
    branch: str,
    requested: float,
    actual: float,
    nodes: dict[int, tuple[float, float]],
    cells: dict[int, list[int]],
    state: dict[str, dict[int, float]],
    provenance: str,
) -> tuple[list[dict[str, Any]], dict[str, float]]:
    all_rows = rows(path)
    selected: list[tuple[int, dict[str, str]]] = [
        (index, row)
        for index, row in enumerate(all_rows)
        if abs(float(row["bias_V"]) - actual) <= EXACT_BIAS_TOLERANCE_V
    ]
    if not selected:
        raise ValueError(f"{branch}: no process rows at {actual:.17g} V")
    result: list[dict[str, Any]] = []
    grouped: dict[tuple[int, str], list[tuple[int, dict[str, str]]]] = defaultdict(list)
    for index, row in selected:
        grouped[(int(row["cell_id"]), row["carrier"])].append((index, row))

    for cell_id, node_ids in cells.items():
        psi_gradient = triangle_gradient(node_ids, nodes, state["psi"])
        result.append(
            field(
                branch=branch,
                requested=requested,
                actual=actual,
                support_kind="cell",
                support_key=f"cell:{cell_id}",
                provenance=provenance,
                carrier="none",
                quantity="electric_field",
                components=["x", "y"],
                unit="V/cm",
                values=[
                    -psi_gradient[0] * UM_INV_TO_CM_INV,
                    -psi_gradient[1] * UM_INV_TO_CM_INV,
                ],
                connectivity=sorted(node_ids),
                source_record=source(
                    root, path, "bv_process_probe", grouped[(cell_id, "electron")][0][0]
                ),
            )
        )
        for carrier, qf_name in (("electron", "phin"), ("hole", "phip")):
            carrier_rows = grouped[(cell_id, carrier)]
            qf_gradient = triangle_gradient(node_ids, nodes, state[qf_name])
            mobility = sum(float(row["final_mobility"]) for _, row in carrier_rows) / len(
                carrier_rows
            )
            current_x = sum(
                float(row["current_vector_x"]) for _, row in carrier_rows
            ) / len(carrier_rows)
            current_y = sum(
                float(row["current_vector_y"]) for _, row in carrier_rows
            ) / len(carrier_rows)
            record_source = source(
                root, path, "bv_process_probe", carrier_rows[0][0]
            )
            common = {
                "branch": branch,
                "requested": requested,
                "actual": actual,
                "support_kind": "cell",
                "support_key": f"cell:{cell_id}",
                "provenance": provenance,
                "carrier": carrier,
                "connectivity": sorted(node_ids),
                "source_record": record_source,
            }
            result.append(
                field(
                    **common,
                    quantity="quasi_fermi_gradient",
                    components=["x", "y"],
                    unit="V/cm",
                    values=[
                        qf_gradient[0] * UM_INV_TO_CM_INV,
                        qf_gradient[1] * UM_INV_TO_CM_INV,
                    ],
                )
            )
            result.append(
                field(
                    **common,
                    quantity="mobility",
                    components=["scalar"],
                    unit="cm^2/(V s)",
                    values=[mobility],
                )
            )
            result.append(
                field(
                    **common,
                    quantity="current_density",
                    components=["x", "y"],
                    unit="A/cm^2",
                    values=[
                        ELEMENTARY_CHARGE_C
                        * current_x
                        * INTERNAL_PARTICLE_FLUX_TO_PER_CM2_S,
                        ELEMENTARY_CHARGE_C
                        * current_y
                        * INTERNAL_PARTICLE_FLUX_TO_PER_CM2_S,
                    ],
                )
            )

    control_areas = nodal_control_areas_um2(nodes, cells)
    alpha_values: dict[tuple[int, str], list[float]] = defaultdict(list)
    generation_integrals: dict[tuple[int, str], float] = defaultdict(float)
    vertex_sources: dict[tuple[int, int, str], float] = defaultdict(float)
    vertex_source_index: dict[tuple[int, int, str], int] = {}
    total_qg = {"electron": 0.0, "hole": 0.0}
    for index, row in selected:
        cell_id = int(row["cell_id"])
        carrier = row["carrier"]
        node0 = int(row["node0"])
        node1 = int(row["node1"])
        alpha = float(row["alpha"])
        alpha_values[(node0, carrier)].append(alpha)
        alpha_values[(node1, carrier)].append(alpha)
        source_integral = float(row["source_integral"])
        total_qg[carrier] += (
            ELEMENTARY_CHARGE_C
            * source_integral
            * INTERNAL_LINE_SOURCE_TO_A_PER_UM
        )
        scatter_nodes = [
            int(value) for value in row["scatter_nodes"].split(";") if value
        ]
        weights = [
            float(value) for value in row["source_weights"].split(";") if value
        ]
        for node, weight in zip(scatter_nodes, weights):
            generation_integrals[(node, carrier)] += source_integral * weight
            local_vertex = cells[cell_id].index(node)
            key = (cell_id, local_vertex, carrier)
            vertex_sources[key] += (
                ELEMENTARY_CHARGE_C
                * source_integral
                * weight
                * INTERNAL_LINE_SOURCE_TO_A_PER_UM
            )
            vertex_source_index.setdefault(key, index)

    for node, coordinates in nodes.items():
        for carrier in ("electron", "hole"):
            values = alpha_values[(node, carrier)]
            alpha = sum(values) / len(values) if values else 0.0
            generation = generation_integrals[(node, carrier)] / control_areas[node]
            record_source = source(
                root,
                path,
                "bv_process_probe",
                next(
                    (
                        index
                        for index, row in selected
                        if row["carrier"] == carrier
                        and node in (int(row["node0"]), int(row["node1"]))
                    ),
                    selected[0][0],
                ),
            )
            common = {
                "branch": branch,
                "requested": requested,
                "actual": actual,
                "support_kind": "physical_node",
                "support_key": f"node:{node}",
                "provenance": "reconstructed",
                "carrier": carrier,
                "coordinates_um": list(coordinates),
                "source_record": record_source,
            }
            result.append(
                field(
                    **common,
                    quantity="avalanche_alpha",
                    components=["scalar"],
                    unit="cm^-1",
                    values=[alpha],
                )
            )
            result.append(
                field(
                    **common,
                    quantity="avalanche_generation",
                    components=["scalar"],
                    unit="cm^-3 s^-1",
                    values=[generation],
                )
            )
        result.append(
            field(
                branch=branch,
                requested=requested,
                actual=actual,
                support_kind="physical_node",
                support_key=f"node:{node}",
                provenance="reconstructed",
                carrier="total",
                quantity="avalanche_generation",
                components=["scalar"],
                unit="cm^-3 s^-1",
                values=[
                    (
                        generation_integrals[(node, "electron")]
                        + generation_integrals[(node, "hole")]
                    )
                    / control_areas[node]
                ],
                coordinates_um=list(coordinates),
                source_record=source(
                    root, path, "bv_process_probe", selected[0][0]
                ),
            )
        )

    for cell_id, node_ids in cells.items():
        for local_vertex, node in enumerate(node_ids):
            carrier_values = {}
            for carrier in ("electron", "hole"):
                key = (cell_id, local_vertex, carrier)
                carrier_values[carrier] = vertex_sources[key]
                result.append(
                    field(
                        branch=branch,
                        requested=requested,
                        actual=actual,
                        support_kind="element_local_vertex",
                        support_key=f"cell:{cell_id}/local_vertex:{local_vertex}",
                        provenance=provenance,
                        carrier=carrier,
                        quantity="integrated_source",
                        components=["scalar"],
                        unit="A/um",
                        values=[carrier_values[carrier]],
                        connectivity=[node],
                        source_record=source(
                            root,
                            path,
                            "bv_process_probe",
                            vertex_source_index.get(key, selected[0][0]),
                        ),
                    )
                )
            result.append(
                field(
                    branch=branch,
                    requested=requested,
                    actual=actual,
                    support_kind="element_local_vertex",
                    support_key=f"cell:{cell_id}/local_vertex:{local_vertex}",
                    provenance=provenance,
                    carrier="total",
                    quantity="integrated_source",
                    components=["scalar"],
                    unit="A/um",
                    values=[
                        carrier_values["electron"] + carrier_values["hole"]
                    ],
                    connectivity=[node],
                    source_record=source(
                        root, path, "bv_process_probe", selected[0][0]
                    ),
                )
            )
    total_qg["total"] = total_qg["electron"] + total_qg["hole"]
    return result, total_qg


def aggregate(
    branch: str,
    requested: float,
    actual: float,
    carrier: str,
    quantity: str,
    unit: str,
    value: float,
    provenance: str,
    source_record: dict[str, Any],
) -> dict[str, Any]:
    return {
        "branch": branch,
        "requested_bias_V": requested,
        "actual_bias_V": actual,
        "carrier": carrier,
        "quantity": quantity,
        "unit": unit,
        "value": value,
        "provenance": provenance,
        "source": source_record,
    }


def exact_row(
    sweep_rows: list[dict[str, str]], requested: float
) -> tuple[int, dict[str, str]]:
    matches = [
        (index, row)
        for index, row in enumerate(sweep_rows)
        if abs(float(row["bias_V"]) - requested) <= EXACT_BIAS_TOLERANCE_V
    ]
    if len(matches) != 1:
        raise ValueError(
            f"expected one sweep row at {requested:.17g} V, got {len(matches)}"
        )
    index, row = matches[0]
    if row["converged"] != "1":
        raise ValueError(f"sweep row did not converge at {requested:.17g} V")
    return index, row


def build_manifest(
    root: Path,
    base_config: Path,
    biases: list[float],
    *,
    run_id: str,
) -> dict[str, Any]:
    root = root.resolve()
    base = json.loads(base_config.read_text(encoding="utf-8-sig"))
    mesh_path = Path(base["mesh_file"]).resolve()
    nodes, cells = load_mesh(mesh_path)

    off_probe = root / "avalanche_off" / "process_probe.csv"
    iic_probe = root / "iic_postprocess" / "process_probe.csv"
    if not iic_probe.is_file():
        raise FileNotFoundError(iic_probe)
    off_states = root / "avalanche_off" / "states"
    iic_states = root / "iic_postprocess" / "states"
    for bias in biases:
        name = f"state_bias_{bias_token(bias)}.csv"
        if sha256(off_states / name) != sha256(iic_states / name):
            raise ValueError(f"off/IIC state drift at {bias:.17g} V")
    shutil.copyfile(iic_probe, off_probe)

    field_records: list[dict[str, Any]] = []
    aggregate_records: list[dict[str, Any]] = []
    newton_records: list[dict[str, Any]] = []
    branch_records: list[dict[str, Any]] = []
    input_hashes: dict[str, str] = {}
    output_files: set[Path] = set()

    for branch in BRANCHES:
        case_dir = root / branch
        config_path = case_dir / "simulation.json"
        sweep_path = case_dir / "iv.csv"
        probe_path = case_dir / "process_probe.csv"
        attempts_path = case_dir / "newton_attempts.csv"
        input_hashes[relative(root, config_path)] = sha256(config_path)
        output_files.update({sweep_path, probe_path, attempts_path})
        sweep_rows = rows(sweep_path)
        provenance = "postprocessed" if branch == "avalanche_off" else "solver_used"
        bias_records = []
        for requested in biases:
            sweep_index, sweep_row = exact_row(sweep_rows, requested)
            actual = float(sweep_row["bias_V"])
            state_path = (
                case_dir / "states" / f"state_bias_{bias_token(requested)}.csv"
            )
            output_files.add(state_path)
            state_result, state = state_fields(
                root, state_path, branch, requested, actual, nodes
            )
            process_result, qg = probe_fields(
                root,
                probe_path,
                branch,
                requested,
                actual,
                nodes,
                cells,
                state,
                provenance,
            )
            field_records.extend(state_result)
            field_records.extend(process_result)
            qg_reintegrated = {
                carrier: sum(
                    float(record["values"][0])
                    for record in process_result
                    if record["quantity"] == "integrated_source"
                    and record["carrier"] == carrier
                )
                for carrier in ("electron", "hole", "total")
            }
            sweep_source = source(root, sweep_path, "sweep", sweep_index)
            for carrier in ("electron", "hole", "total"):
                for source_provenance, value in (
                    (provenance, qg[carrier]),
                    ("operator_replay", qg_reintegrated[carrier]),
                ):
                    aggregate_records.append(
                        aggregate(
                            branch,
                            requested,
                            actual,
                            carrier,
                            "integrated_source",
                            "A/um",
                            value,
                            source_provenance,
                            source(root, probe_path, "bv_process_probe", 0),
                        )
                    )
            terminal_columns = {
                "electron": "current_electron_A_per_um",
                "hole": "current_hole_A_per_um",
                "total": "current_total_A_per_um",
            }
            for carrier, column in terminal_columns.items():
                aggregate_records.append(
                    aggregate(
                        branch,
                        requested,
                        actual,
                        carrier,
                        "terminal_current",
                        "A/um",
                        float(sweep_row[column]),
                        "solver_used",
                        sweep_source,
                    )
                )
            bias_records.append(
                {
                    "requested_bias_V": requested,
                    "actual_bias_V": actual,
                    "snapshot_tdr": artifact(root, state_path),
                    "currentplot": artifact(root, sweep_path),
                    "process_record": artifact(root, probe_path),
                }
            )
        branch_records.append(
            {
                "branch": branch,
                "requested_biases_V": biases,
                "bias_records": bias_records,
            }
        )
        for index, row in enumerate(rows(attempts_path)):
            requested = float(row["requested_target_bias_V"])
            actual = float(row["actual_target_bias_V"])
            if not any(
                abs(requested - bias) <= EXACT_BIAS_TOLERANCE_V
                for bias in biases
            ):
                continue
            if abs(requested - actual) > EXACT_BIAS_TOLERANCE_V:
                continue
            status = row["status"]
            newton_records.append(
                {
                    "branch": branch,
                    "attempt_id": f"{branch}:{row['attempt_id']}",
                    "requested_bias_V": requested,
                    "actual_bias_V": actual,
                    "status": status if status in {"accepted", "rejected"} else "failed",
                    "reason": row["reason"] or status,
                    "source": source(
                        root, attempts_path, "newton_attempts", index
                    ),
                }
            )

    field_path = root / "field_records.jsonl"
    aggregate_path = root / "aggregate_records.jsonl"
    newton_path = root / "newton_attempt_records.jsonl"
    write_jsonl(field_path, field_records)
    write_jsonl(aggregate_path, aggregate_records)
    write_jsonl(newton_path, newton_records)
    output_files.update({field_path, aggregate_path, newton_path})
    normalized_output_hashes = {
        relative(root, path): sha256(path)
        for path in sorted(output_files)
    }
    manifest = {
        "schema": "vela.pn2d_bv_process_run.v1",
        "status": "passed",
        "outcome": "complete_exact_lattice_process_manifest",
        "run_id": run_id,
        "simulator": "vela",
        "release": "current-code-opt-in-newton-80",
        "missing_value_policy": "reject",
        "input_hashes": input_hashes,
        "normalized_output_hashes": normalized_output_hashes,
        "branch_records": branch_records,
        "field_records": field_records,
        "aggregate_records": aggregate_records,
        "newton_attempt_records": newton_records,
    }
    validate_process_run(manifest, base_dir=root)
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--base-config", type=Path, required=True)
    parser.add_argument("--sentaurus-manifest", type=Path, required=True)
    parser.add_argument("--run-id")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = args.run_root.resolve()
    sentaurus = json.loads(
        args.sentaurus_manifest.read_text(encoding="utf-8")
    )
    records = {
        row["branch"]: [float(value) for value in row["requested_biases_V"]]
        for row in sentaurus["branch_records"]
        if row["branch"] in BRANCHES
    }
    biases = records[BRANCHES[0]]
    if any(records.get(branch) != biases for branch in BRANCHES):
        raise ValueError("Sentaurus branches do not share one exact lattice")
    manifest = build_manifest(
        root,
        args.base_config.resolve(),
        biases,
        run_id=args.run_id or root.name,
    )
    path = root / "manifest.json"
    path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(
        json.dumps(
            {
                "manifest": str(path),
                "sha256": sha256(path),
                "field_records": len(manifest["field_records"]),
                "aggregate_records": len(manifest["aggregate_records"]),
                "newton_attempt_records": len(
                    manifest["newton_attempt_records"]
                ),
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
