#!/usr/bin/env python3
"""Augment PN2D process-chain inputs with fixed-transition Newton observations."""

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

from scripts.analyze_pn2d_bv_process_chain import CHAIN_SCHEMA, from_process_run


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def normalized(values: list[float]) -> list[float]:
    norm = math.sqrt(sum(value * value for value in values))
    if not math.isfinite(norm):
        raise ValueError("nonfinite Newton observation norm")
    if norm == 0.0:
        return [0.0 for _ in values]
    return [value / norm for value in values]


def process_input(path: Path) -> dict[str, Any]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    return from_process_run(payload)


def record(
    *,
    simulator: str,
    branch: str,
    bias: float,
    stage: str,
    quantity: str,
    carrier: str,
    node: int,
    value: float,
    coordinates: list[float],
    provenance: str,
) -> dict[str, Any]:
    return {
        "simulator": simulator,
        "branch": branch,
        "bias_V": bias,
        "stage": stage,
        "quantity": quantity,
        "carrier": carrier,
        "support_kind": "physical_node",
        "support_key": f"v{node}",
        "values": [value],
        "unit": "1",
        "coordinates_um": coordinates,
        "provenance": provenance,
    }


def field_map(payload: dict[str, Any]) -> dict[str, dict[str, Any]]:
    return {str(field["name"]): field for field in payload["fields"]}


def sentaurus_records(
    execution_path: Path,
    importer: Path,
) -> tuple[list[dict[str, Any]], dict[str, str]]:
    execution = json.loads(execution_path.read_text(encoding="utf-8"))
    result: list[dict[str, Any]] = []
    hashes: dict[str, str] = {str(execution_path): sha256(execution_path)}
    field_names = {
        "residual_signature": (
            ("none", "PoissonRhs"),
            ("electron", "eContinuityRhs"),
            ("hole", "hContinuityRhs"),
        ),
        "newton_update_signature": (
            ("none", "NewtonStepElectrostaticPotentialUpdate"),
            ("electron", "NewtonStepEDensityUpdate"),
            ("hole", "NewtonStepHDensityUpdate"),
        ),
    }
    for case in execution["cases"]:
        branch = str(case["branch"])
        bias = float(case["target_bias_V"])
        tdr = Path(case["first_update_tdr"])
        case_root = tdr.parent.parent
        values_path = case_root / "field_values.json"
        completed = subprocess.run(
            [
                str(importer),
                "--tdr",
                str(tdr),
                "--field-values-json",
                str(values_path),
            ],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        if completed.returncode != 0:
            raise RuntimeError(
                f"sentaurus_import failed for {tdr}: "
                f"{completed.stderr.strip() or completed.stdout.strip()}"
            )
        payload = json.loads(values_path.read_text(encoding="utf-8"))
        fields = field_map(payload)
        vertices = payload["geometry"]["vertices"]
        hashes[str(tdr)] = sha256(tdr)
        hashes[str(values_path)] = sha256(values_path)
        for quantity, definitions in field_names.items():
            stage = (
                "residual_jacobian"
                if quantity == "residual_signature"
                else "newton_update"
            )
            for carrier, name in definitions:
                if name not in fields:
                    raise ValueError(f"{tdr}: missing {name}")
                values = [float(value) for value in fields[name]["raw_values"]]
                signature = normalized(values)
                if len(signature) != len(vertices):
                    raise ValueError(f"{tdr}: {name} support mismatch")
                for node, value in enumerate(signature):
                    result.append(
                        record(
                            simulator="sentaurus",
                            branch=branch,
                            bias=bias,
                            stage=stage,
                            quantity=quantity,
                            carrier=carrier,
                            node=node,
                            value=value,
                            coordinates=[float(x) for x in vertices[node]],
                            provenance="native",
                        )
                    )
    return result, hashes


def csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def vela_records(
    execution_path: Path,
) -> tuple[list[dict[str, Any]], dict[str, str]]:
    execution = json.loads(execution_path.read_text(encoding="utf-8"))
    result: list[dict[str, Any]] = []
    hashes: dict[str, str] = {str(execution_path): sha256(execution_path)}
    for case in execution["cases"]:
        branch = str(case["branch"])
        bias = float(case["target_bias_V"])
        probe_path = Path(case["csv"])
        source_state_path = Path(case["source_state"])
        probe = csv_rows(probe_path)
        source_state = {
            int(row["node_id"]): row for row in csv_rows(source_state_path)
        }
        if len(probe) != len(source_state):
            raise ValueError(f"{probe_path}: state support mismatch")
        hashes[str(probe_path)] = sha256(probe_path)
        hashes[str(source_state_path)] = sha256(source_state_path)
        residual_columns = (
            ("none", "psi_residual"),
            ("electron", "phin_residual"),
            ("hole", "phip_residual"),
        )
        update_values = {
            "none": [float(row["delta_psi_V"]) for row in probe],
            "electron": [
                (
                    float(row["trial_electron_density_m3"])
                    - float(source_state[int(row["node_id"])]["electrons_m3"])
                )
                * 1.0e-6
                for row in probe
            ],
            "hole": [
                (
                    float(row["trial_hole_density_m3"])
                    - float(source_state[int(row["node_id"])]["holes_m3"])
                )
                * 1.0e-6
                for row in probe
            ],
        }
        for carrier, column in residual_columns:
            signature = normalized([float(row[column]) for row in probe])
            for index, (row, value) in enumerate(zip(probe, signature)):
                node = int(row["node_id"])
                if node != index:
                    raise ValueError(f"{probe_path}: noncanonical node order")
                result.append(
                    record(
                        simulator="vela",
                        branch=branch,
                        bias=bias,
                        stage="residual_jacobian",
                        quantity="residual_signature",
                        carrier=carrier,
                        node=node,
                        value=value,
                        coordinates=[float(row["x"]), float(row["y"])],
                        provenance="solver_used",
                    )
                )
        for carrier, values in update_values.items():
            signature = normalized(values)
            for index, (row, value) in enumerate(zip(probe, signature)):
                node = int(row["node_id"])
                result.append(
                    record(
                        simulator="vela",
                        branch=branch,
                        bias=bias,
                        stage="newton_update",
                        quantity="newton_update_signature",
                        carrier=carrier,
                        node=node,
                        value=value,
                        coordinates=[float(row["x"]), float(row["y"])],
                        provenance="solver_used",
                    )
                )
    return result, hashes


def write_chain(
    output: Path,
    base: dict[str, Any],
    observations: list[dict[str, Any]],
    *,
    process_manifest: Path,
    execution: Path,
    hashes: dict[str, str],
) -> None:
    payload = {
        **base,
        "schema": CHAIN_SCHEMA,
        "records": [*base["records"], *observations],
        "newton_observation_contract": {
            "transition": "previous_exact_accepted_state_to_target_exact_bias",
            "iteration": 1,
            "residual_scaling": "per_equation_l2_spatial_signature",
            "update_scaling": "per_variable_l2_spatial_signature",
            "jacobian_observation": "implicit_inverse_action_delta_x_equals_minus_J_inverse_R",
            "full_jacobian_matrix_available": False,
            "process_manifest": str(process_manifest),
            "process_manifest_sha256": sha256(process_manifest),
            "execution": str(execution),
            "artifact_hashes": hashes,
        },
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sentaurus-manifest", type=Path, required=True)
    parser.add_argument("--vela-manifest", type=Path, required=True)
    parser.add_argument("--sentaurus-execution", type=Path, required=True)
    parser.add_argument("--vela-execution", type=Path, required=True)
    parser.add_argument("--sentaurus-import", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    sentaurus_manifest = args.sentaurus_manifest.resolve()
    vela_manifest = args.vela_manifest.resolve()
    sentaurus_execution = args.sentaurus_execution.resolve()
    vela_execution = args.vela_execution.resolve()
    sentaurus, sentaurus_hashes = sentaurus_records(
        sentaurus_execution, args.sentaurus_import.resolve()
    )
    vela, vela_hashes = vela_records(vela_execution)
    output = args.output_root.resolve()
    write_chain(
        output / "sentaurus_chain_input.json",
        process_input(sentaurus_manifest),
        sentaurus,
        process_manifest=sentaurus_manifest,
        execution=sentaurus_execution,
        hashes=sentaurus_hashes,
    )
    write_chain(
        output / "vela_chain_input.json",
        process_input(vela_manifest),
        vela,
        process_manifest=vela_manifest,
        execution=vela_execution,
        hashes=vela_hashes,
    )
    summary = {
        "schema": "vela.pn2d_bv_newton_chain_inputs.v1",
        "status": "passed",
        "sentaurus_observation_records": len(sentaurus),
        "vela_observation_records": len(vela),
        "sentaurus_chain_input": str(output / "sentaurus_chain_input.json"),
        "vela_chain_input": str(output / "vela_chain_input.json"),
    }
    (output / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
