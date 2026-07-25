#!/usr/bin/env python3
"""Replay the opt-in element-edge GSS/Laux operator on six Minimal6 states."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
import subprocess
from pathlib import Path

Q_C = 1.602176634e-19


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--binary",
        type=Path,
        default=Path("build-release/pn2d_minimal6_operator_audit.exe"),
    )
    parser.add_argument(
        "--state-root",
        type=Path,
        default=Path(
            "build-release/"
            "pn2d-minimal6-inverse-sentaurus-fields-20260717-d/"
            "minimal6_inverse_fields_20260717_d"
        ),
    )
    parser.add_argument(
        "--config-root",
        type=Path,
        default=Path(
            "build-release/"
            "pn2d-minimal6-qfp-sg-bootstrap-unitfix-20260723-b/"
            "baseline_replay"
        ),
    )
    parser.add_argument(
        "--density-control",
        type=Path,
        default=Path(
            "build-release/"
            "pn2d-minimal6-sentaurus-box-staged-sweep-20260724-a/"
            "density_recompute_control.csv"
        ),
    )
    parser.add_argument(
        "--reference-root",
        type=Path,
        default=Path(
            "build-release/"
            "pn2d-minimal6-element-avalanche-replay-20260725/"
            "analysis"
        ),
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=Path(
            "build-release/"
            "pn2d-minimal6-element-edge-gss-laux-fixed-state-20260725"
        ),
    )
    parser.add_argument(
        "--magnitudes",
        type=int,
        nargs="+",
        default=(1, 10, 20),
    )
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="ascii") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        raise RuntimeError(f"refusing to write empty CSV: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="ascii") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def abs_log_error(candidate: float, reference: float) -> float | None:
    if candidate <= 0.0 or reference <= 0.0:
        return None
    return abs(math.log10(candidate / reference))


def summary(values: list[float | None]) -> dict[str, float | int | None]:
    finite = [value for value in values if value is not None and math.isfinite(value)]
    return {
        "finite_count": len(finite),
        "median_absolute_error_dex": statistics.median(finite) if finite else None,
        "max_absolute_error_dex": max(finite) if finite else None,
    }


def configure(base: dict[str, object]) -> dict[str, object]:
    configured = json.loads(json.dumps(base))
    solver = configured["solver"]
    assert isinstance(solver, dict)
    impact = solver["impact_ionization"]
    assert isinstance(impact, dict)
    impact.update(
        {
            "generation": "current_density",
            "driving_force": "electric_field",
            "current_approximation": "element_edge_sg_gss_laux",
            "quasi_fermi_gradient_discretization": "edge_difference",
            "source_mapping_mode": "element_vertex_box_measure",
        }
    )
    return configured


def main() -> int:
    args = parse_args()
    magnitudes = tuple(args.magnitudes)
    if (
        len(set(magnitudes)) != len(magnitudes)
        or any(item < 1 or item > 20 for item in magnitudes)
    ):
        raise ValueError("magnitudes must be unique integers in [1, 20]")
    binary = args.binary.resolve()
    state_root = args.state_root.resolve()
    reference_root = args.reference_root.resolve()
    config_root = args.config_root.resolve()
    density_control_path = args.density_control.resolve()
    output_root = args.output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    reference_edge_rows = read_csv(reference_root / "element_edges.csv")
    ref_edges = {
        (
            row["topology"],
            int(round(float(row["bias_V"]))),
            int(row["element"]),
            min(int(row["start"]), int(row["end"])),
            max(int(row["start"]), int(row["end"])),
        ): row
        for row in reference_edge_rows
    }
    reference_cell_nodes: dict[tuple[str, int, int], set[int]] = {}
    for row in reference_edge_rows:
        key = (
            row["topology"],
            int(round(float(row["bias_V"]))),
            int(row["element"]),
        )
        reference_cell_nodes.setdefault(key, set()).update(
            (int(row["start"]), int(row["end"]))
        )
    reference_element_by_nodes = {
        (topology, bias, tuple(sorted(nodes))): element
        for (topology, bias, element), nodes in reference_cell_nodes.items()
    }
    ref_vectors = {
        (
            row["topology"],
            int(round(float(row["bias_V"]))),
            int(row["element"]),
            row["carrier"],
        ): row
        for row in read_csv(reference_root / "element_reconstructions.csv")
        if row["candidate"] == "gss_laux_edge_volume_weighted"
    }
    ref_vertices = {
        (
            row["topology"],
            int(round(float(row["bias_V"]))),
            int(row["element"]),
            int(row["vertex"]),
        ): row
        for row in read_csv(reference_root / "element_vertex_measures.csv")
    }
    recomputed_density = {
        (
            row["topology"],
            int(round(float(row["bias_V"]))),
            int(row["node"]),
        ): (float(row["recomputed_n_m3"]), float(row["recomputed_p_m3"]))
        for row in read_csv(density_control_path)
    }

    detail_rows: list[dict[str, object]] = []
    generated: list[Path] = []
    for topology in ("mirror", "sketch"):
        for magnitude in magnitudes:
            bias = -magnitude
            state_dir = state_root / "states" / topology / f"m{magnitude}V" / "export"
            config_dir = config_root / topology / f"m{magnitude}V"
            out_dir = output_root / topology / f"m{magnitude}V"
            out_dir.mkdir(parents=True, exist_ok=True)
            config_path = out_dir / "audit_config.json"
            with (config_dir / "audit_config.json").open(
                encoding="ascii"
            ) as stream:
                config = configure(json.load(stream))
            with config_path.open("w", encoding="ascii", newline="\n") as stream:
                json.dump(config, stream, indent=2, sort_keys=True)
                stream.write("\n")

            state_rows = read_csv(state_dir / "state.csv")
            for state_row in state_rows:
                node = int(state_row["node_id"])
                density_n, density_p = recomputed_density[
                    (topology, bias, node)
                ]
                state_row["n_m3"] = format(density_n, ".17g")
                state_row["p_m3"] = format(density_p, ".17g")
            state_path = out_dir / "state_recomputed_density.csv"
            write_csv(state_path, state_rows)
            element_path = out_dir / "element_edge_gss_laux.csv"
            command = [
                str(binary),
                "--mesh",
                str(config["mesh_file"]),
                "--doping",
                str(config["node_doping_file"]),
                "--state",
                str(state_path),
                "--config",
                str(config_path),
                "--node-out",
                str(out_dir / "nodes.csv"),
                "--edge-out",
                str(out_dir / "edges.csv"),
                "--triangle-out",
                str(out_dir / "triangles.csv"),
                "--element-out",
                str(element_path),
            ]
            subprocess.run(command, check=True)
            generated.extend(
                [
                    config_path,
                    state_path,
                    out_dir / "nodes.csv",
                    out_dir / "edges.csv",
                    out_dir / "triangles.csv",
                    element_path,
                ]
            )

            element_rows = read_csv(element_path)
            cpp_cell_nodes: dict[int, set[int]] = {}
            for element_row in element_rows:
                cpp_cell_nodes.setdefault(
                    int(element_row["cell_id"]), set()
                ).update(
                    (int(element_row["node_id"]), int(element_row["next_node_id"]))
                )
            for row in element_rows:
                cell = int(row["cell_id"])
                local = int(row["local_index"])
                node = int(row["node_id"])
                next_node = int(row["next_node_id"])
                reference_cell = reference_element_by_nodes[
                    (topology, bias, tuple(sorted(cpp_cell_nodes[cell])))
                ]
                ref_edge = ref_edges[
                    (
                        topology, bias, reference_cell,
                        min(node, next_node), max(node, next_node),
                    )
                ]
                ref_vertex = ref_vertices[
                    (topology, bias, reference_cell, node)
                ]
                output: dict[str, object] = {
                    "topology": topology,
                    "bias_V": bias,
                    "cell_id": cell,
                    "local_index": local,
                    "sentaurus_element_id": reference_cell,
                    "node_id": node,
                    "edge_id": int(row["edge_id"]),
                    "edge_partial_volume_m2": float(
                        row["edge_partial_volume_m2"]
                    ),
                    "vertex_measure_m2": float(row["vertex_measure_m2"]),
                }
                for carrier, short in (("electron", "n"), ("hole", "p")):
                    cpp_edge = (
                        abs(float(row[f"{carrier}_signed_edge_flux_per_m2_s"]))
                        * Q_C
                        / 1.0e4
                    )
                    ref_edge_current = abs(float(ref_edge[f"sg_j{short}_A_cm2"]))
                    cpp_vector = (
                        float(row[f"{carrier}_current_magnitude_per_m2_s"])
                        * Q_C
                        / 1.0e4
                    )
                    ref_vector = ref_vectors[
                        (topology, bias, reference_cell, carrier)
                    ]
                    ref_vector_current = float(ref_vector["magnitude_A_cm2"])
                    cpp_alpha = float(row[f"{carrier}_alpha_per_m"]) / 100.0
                    ref_alpha = float(ref_vector["alpha_cm_inv"])
                    cpp_qg = (
                        float(row[f"{carrier}_source_integral_per_m_s"])
                        * Q_C
                        * 1.0e-6
                    )
                    ref_qg = float(ref_vertex[f"qg_{short}_A_um"])
                    output.update(
                        {
                            f"{carrier}_cpp_edge_A_cm2": cpp_edge,
                            f"{carrier}_sentaurus_edge_replay_A_cm2": (
                                ref_edge_current
                            ),
                            f"{carrier}_edge_error_dex": abs_log_error(
                                cpp_edge, ref_edge_current
                            ),
                            f"{carrier}_cpp_vector_A_cm2": cpp_vector,
                            f"{carrier}_sentaurus_vector_replay_A_cm2": (
                                ref_vector_current
                            ),
                            f"{carrier}_vector_error_dex": abs_log_error(
                                cpp_vector, ref_vector_current
                            ),
                            f"{carrier}_cpp_alpha_cm_inv": cpp_alpha,
                            f"{carrier}_sentaurus_alpha_cm_inv": ref_alpha,
                            f"{carrier}_alpha_error_dex": abs_log_error(
                                cpp_alpha, ref_alpha
                            ),
                            f"{carrier}_cpp_qg_A_um": cpp_qg,
                            f"{carrier}_sentaurus_qg_A_um": ref_qg,
                            f"{carrier}_qg_error_dex": abs_log_error(
                                cpp_qg, ref_qg
                            ),
                        }
                    )
                detail_rows.append(output)

    detail_path = output_root / "fixed_state_comparison.csv"
    write_csv(detail_path, detail_rows)
    generated.append(detail_path)

    metrics: dict[str, object] = {}
    for carrier in ("electron", "hole"):
        for quantity in ("edge", "vector", "alpha", "qg"):
            metrics[f"{carrier}_{quantity}"] = summary(
                [
                    row[f"{carrier}_{quantity}_error_dex"]
                    for row in detail_rows
                ]
            )
    manifest = {
        "schema_version": 1,
        "status": "valid_diagnostic_replay",
        "scope": {
            "topologies": ["mirror", "sketch"],
            "biases_V": [-magnitude for magnitude in magnitudes],
            "state_semantics": (
                "Sentaurus psi/phin/phip imported; Vela n/p, mobility, SG, "
                "GSS/Laux vector, alpha, and source recomputed"
            ),
            "native_sentaurus_edge_current_observed": False,
            "reference_semantics": "documented Sentaurus box-operator replay",
            "production_default_changed": False,
        },
        "row_count": len(detail_rows),
        "metrics": metrics,
        "output_sha256": {
            str(path.relative_to(output_root)).replace("\\", "/"): sha256(path)
            for path in generated
        },
    }
    manifest_path = output_root / "manifest.json"
    with manifest_path.open("w", encoding="ascii", newline="\n") as stream:
        json.dump(manifest, stream, indent=2, sort_keys=True)
        stream.write("\n")
    print(json.dumps(manifest["metrics"], indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
