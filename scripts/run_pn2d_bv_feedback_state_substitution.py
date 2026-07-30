#!/usr/bin/env python3
"""Run density/QFP one-stage feedback substitutions at exact PN2D knee states."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
import sys
from pathlib import Path
from typing import Any

if not __package__:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.diagnose_pn2d_bv_predictor_first_step_audit import (
    make_probe_config,
    state_csv_to_fields,
)


EXACT_TOLERANCE_V = 1.0e-10
COORDINATE_TOLERANCE_UM = 1.0e-12
REQUIRED_QUANTITIES = (
    ("quasi_fermi", "electron", "eQuasiFermiPotential", "V", 1.0),
    ("quasi_fermi", "hole", "hQuasiFermiPotential", "V", 1.0),
    ("density", "electron", "eDensity_m3", "cm^-3", 1.0e6),
    ("density", "hole", "hDensity_m3", "cm^-3", 1.0e6),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sentaurus-chain", type=Path, required=True)
    parser.add_argument("--vela-manifest", type=Path, required=True)
    parser.add_argument("--base-config", type=Path, required=True)
    parser.add_argument("--runner", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--biases", nargs="+", type=float, default=[-19.7, -19.8])
    parser.add_argument("--branch", default="avalanche_on")
    parser.add_argument("--resume", action="store_true")
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def bias_token(bias: float) -> str:
    return f"{abs(bias):.6f}".replace(".", "p").join(
        ("m" if bias < 0 else "p", "")
    )


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def branch_record(manifest: dict[str, Any], branch: str) -> dict[str, Any]:
    matches = [
        item for item in manifest["branch_records"]
        if str(item["branch"]) == branch
    ]
    if len(matches) != 1:
        raise ValueError(f"expected one {branch!r} branch record")
    return matches[0]


def exact_state_record(
    manifest: dict[str, Any],
    branch: str,
    bias: float,
) -> dict[str, Any]:
    matches = [
        record for record in branch_record(manifest, branch)["bias_records"]
        if abs(float(record["requested_bias_V"]) - bias) <= EXACT_TOLERANCE_V
    ]
    if len(matches) != 1:
        raise ValueError(f"{branch} {bias:g} V: expected one exact state")
    record = matches[0]
    if abs(float(record["actual_bias_V"]) - bias) > EXACT_TOLERANCE_V:
        raise ValueError(f"{branch} {bias:g} V: actual bias is not exact")
    return record


def sentaurus_node_records(
    chain: dict[str, Any],
    branch: str,
    bias: float,
    quantity: str,
    carrier: str,
    unit: str,
    expected_node_count: int,
) -> dict[int, dict[str, Any]]:
    result: dict[int, dict[str, Any]] = {}
    for record in chain["records"]:
        if (
            str(record["branch"]) != branch
            or abs(float(record["bias_V"]) - bias) > EXACT_TOLERANCE_V
            or str(record["quantity"]) != quantity
            or str(record["carrier"]) != carrier
            or str(record["support_kind"]) != "physical_node"
            or str(record["provenance"]) != "native"
        ):
            continue
        if str(record["unit"]) != unit:
            raise ValueError(
                f"{branch} {bias:g} V {quantity}/{carrier}: "
                f"expected {unit}, got {record['unit']}"
            )
        key = str(record["support_key"])
        if not key.startswith("node:"):
            raise ValueError(f"unexpected Sentaurus physical-node key {key!r}")
        node = int(key.split(":", 1)[1])
        if node in result:
            raise ValueError(f"duplicate Sentaurus node {node}")
        result[node] = record
    canonical = {
        node: result[node]
        for node in range(expected_node_count)
        if node in result
    }
    if sorted(canonical) != list(range(expected_node_count)):
        raise ValueError(
            f"{branch} {bias:g} V {quantity}/{carrier}: noncanonical node support"
        )
    canonical_coordinates = {
        tuple(float(value) for value in record["coordinates_um"][:2])
        for record in canonical.values()
    }
    for node, record in result.items():
        if node < expected_node_count:
            continue
        coordinates = tuple(
            float(value) for value in record["coordinates_um"][:2]
        )
        if coordinates not in canonical_coordinates:
            raise ValueError(
                f"{branch} {bias:g} V {quantity}/{carrier}: "
                f"unmapped extra support node {node}"
            )
    return canonical


def write_scalar_field(path: Path, values: list[float]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["node_id", "component0"])
        for node, value in enumerate(values):
            if not math.isfinite(value):
                raise ValueError(f"{path}: nonfinite node {node}")
            writer.writerow([node, f"{value:.17g}"])


def mesh_coordinates(base_config: Path) -> list[tuple[float, float]]:
    config = json.loads(base_config.read_text(encoding="utf-8-sig"))
    mesh_path = Path(config["mesh_file"])
    if not mesh_path.is_absolute():
        mesh_path = (base_config.parent / mesh_path).resolve()
    mesh = json.loads(mesh_path.read_text(encoding="utf-8"))
    nodes = sorted(mesh["nodes"], key=lambda item: int(item["id"]))
    if [int(node["id"]) for node in nodes] != list(range(len(nodes))):
        raise ValueError(f"{mesh_path}: noncanonical node order")
    return [(float(node["x"]), float(node["y"])) for node in nodes]


def build_replacement_fields(
    chain: dict[str, Any],
    branch: str,
    bias: float,
    baseline_rows: list[dict[str, str]],
    baseline_coordinates: list[tuple[float, float]],
    output: Path,
) -> dict[str, str]:
    hashes: dict[str, str] = {}
    node_count = len(baseline_rows)
    for quantity, carrier, field_name, unit, factor in REQUIRED_QUANTITIES:
        records = sentaurus_node_records(
            chain, branch, bias, quantity, carrier, unit, node_count
        )
        if len(records) != node_count:
            raise ValueError(
                f"{branch} {bias:g} V {quantity}/{carrier}: "
                "Sentaurus/Vela node-count mismatch"
            )
        values: list[float] = []
        for node, baseline in enumerate(baseline_rows):
            record = records[node]
            coordinates = [float(value) for value in record["coordinates_um"][:2]]
            error = math.hypot(
                coordinates[0] - baseline_coordinates[node][0],
                coordinates[1] - baseline_coordinates[node][1],
            )
            if error > COORDINATE_TOLERANCE_UM:
                raise ValueError(
                    f"{branch} {bias:g} V node {node}: coordinate mismatch {error:g} um"
                )
            value = float(record["values"][0]) * factor
            if quantity == "density" and value <= 0.0:
                raise ValueError(
                    f"{branch} {bias:g} V {carrier} density is not positive"
                )
            values.append(value)
        path = output / f"{field_name}_region0.csv"
        write_scalar_field(path, values)
        hashes[str(path)] = sha256(path)
    return hashes


def run_case(
    *,
    chain: dict[str, Any],
    manifest: dict[str, Any],
    manifest_root: Path,
    branch: str,
    bias: float,
    base_config: Path,
    runner: Path,
    output: Path,
    resume: bool,
) -> dict[str, Any]:
    record = exact_state_record(manifest, branch, bias)
    state_path = manifest_root / str(record["snapshot_tdr"]["path"])
    baseline_rows = read_csv(state_path)
    if [int(row["node_id"]) for row in baseline_rows] != list(range(len(baseline_rows))):
        raise ValueError(f"{state_path}: noncanonical node order")
    coordinates = mesh_coordinates(base_config)
    if len(coordinates) != len(baseline_rows):
        raise ValueError(f"{state_path}: state/mesh node-count mismatch")

    case_root = output / branch / bias_token(bias)
    baseline_fields = case_root / "baseline_state_fields"
    replacement_fields = case_root / "feedback_state_fields"
    output_csv = case_root / "feedback_state_substitution.csv"
    config_path = case_root / "feedback_state_substitution.json"
    status_path = case_root / "status.json"
    case_root.mkdir(parents=True, exist_ok=True)
    state_csv_to_fields(state_path, baseline_fields)
    artifact_hashes = build_replacement_fields(
        chain, branch, bias, baseline_rows, coordinates, replacement_fields
    )
    artifact_hashes[str(state_path)] = sha256(state_path)

    config = make_probe_config(
        base_config,
        output_csv,
        baseline_fields,
        "newton_feedback_substitution_probe",
        bias,
        "Anode",
        "Cathode",
    )
    config["feedback_state_fields_dir"] = str(replacement_fields.resolve())
    config_path.write_text(
        json.dumps(config, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    if not (
        resume
        and output_csv.is_file()
        and status_path.is_file()
    ):
        completed = subprocess.run(
            [str(runner), "--config", str(config_path)],
            cwd=case_root,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        if completed.returncode != 0:
            raise RuntimeError(
                f"{branch} {bias:g} V: "
                f"{completed.stderr.strip() or completed.stdout.strip()}"
            )
        status = json.loads(completed.stdout)
        status_path.write_text(
            json.dumps(status, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
    artifact_hashes[str(config_path)] = sha256(config_path)
    artifact_hashes[str(output_csv)] = sha256(output_csv)
    artifact_hashes[str(status_path)] = sha256(status_path)
    return {
        "branch": branch,
        "bias_V": bias,
        "baseline_state": str(state_path),
        "config": str(config_path),
        "output_csv": str(output_csv),
        "status": str(status_path),
        "artifact_hashes": artifact_hashes,
    }


def main() -> int:
    args = parse_args()
    sentaurus_chain_path = args.sentaurus_chain.resolve()
    vela_manifest_path = args.vela_manifest.resolve()
    base_config = args.base_config.resolve()
    runner = args.runner.resolve()
    output = args.output_root.resolve()
    output.mkdir(parents=True, exist_ok=True)
    chain = json.loads(sentaurus_chain_path.read_text(encoding="utf-8"))
    manifest = json.loads(vela_manifest_path.read_text(encoding="utf-8"))
    cases = [
        run_case(
            chain=chain,
            manifest=manifest,
            manifest_root=vela_manifest_path.parent,
            branch=args.branch,
            bias=float(bias),
            base_config=base_config,
            runner=runner,
            output=output,
            resume=args.resume,
        )
        for bias in args.biases
    ]
    execution = {
        "schema": "vela.pn2d_bv_feedback_state_substitution_execution.v1",
        "status": "passed",
        "outcome": "feedback_state_substitution_observations_available",
        "sentaurus_chain": str(sentaurus_chain_path),
        "sentaurus_chain_sha256": sha256(sentaurus_chain_path),
        "vela_manifest": str(vela_manifest_path),
        "vela_manifest_sha256": sha256(vela_manifest_path),
        "base_config": str(base_config),
        "base_config_sha256": sha256(base_config),
        "runner": str(runner),
        "runner_sha256": sha256(runner),
        "contract": {
            "baseline": "converged_vela_avalanche_on_exact_state",
            "replacement": "sentaurus_avalanche_on_exact_state",
            "variants": [
                "baseline",
                "electron_density_only",
                "hole_density_only",
                "density_only",
                "electron_qfp_only",
                "hole_qfp_only",
                "qfp_only",
                "density_qfp",
            ],
            "jacobian": "single_production_baseline_jacobian",
            "boundary_rows": "baseline_preserved",
            "production_defaults_changed": False,
        },
        "cases": cases,
    }
    execution_path = output / "execution.json"
    execution_path.write_text(
        json.dumps(execution, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(execution, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
