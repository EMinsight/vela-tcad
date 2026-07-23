#!/usr/bin/env python3
"""Independently verify Minimal6 native-element transport closure outputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


Q_C = 1.602176634e-19


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def close(actual: float, expected: float, tolerance: float = 2.0e-12) -> bool:
    return abs(actual - expected) <= tolerance * max(
        1.0, abs(actual), abs(expected)
    )


def verify(
    *,
    transport_csv: Path,
    transport_manifest: Path,
    closure_root: Path,
    output: Path,
) -> dict[str, object]:
    transport_csv = transport_csv.resolve()
    transport_manifest = transport_manifest.resolve()
    closure_root = closure_root.resolve()
    transport_meta = json.loads(
        transport_manifest.read_text(encoding="utf-8")
    )
    closure_meta_path = closure_root / "manifest.json"
    closure_meta = json.loads(closure_meta_path.read_text(encoding="utf-8"))
    if transport_meta.get("status") != "valid":
        raise ValueError("transport manifest is not valid")
    if closure_meta.get("status") != "valid":
        raise ValueError("closure manifest is not valid")
    expected_transport_hash = transport_meta["outputs"][
        "transport_element_values_sha256"
    ]
    if sha256(transport_csv) != expected_transport_hash:
        raise ValueError("transport CSV hash mismatch")
    for name, digest in closure_meta["outputs"].items():
        if sha256(closure_root / name) != digest:
            raise ValueError(f"closure output hash mismatch: {name}")

    with transport_csv.open(newline="", encoding="utf-8") as handle:
        transport_rows = list(csv.DictReader(handle))
    with (closure_root / "transport_element_closure_samples.csv").open(
        newline="", encoding="utf-8"
    ) as handle:
        closure_rows = list(csv.DictReader(handle))
    transport = {
        (row["topology"], float(row["bias_V"]), int(row["cell_id"])): row
        for row in transport_rows
    }
    closure = {
        (
            row["topology"],
            float(row["bias_V"]),
            int(row["cell_id"]),
            row["carrier"],
        ): row
        for row in closure_rows
    }
    if len(transport_rows) != 160 or len(transport) != 160:
        raise ValueError("transport cardinality is not 160")
    if len(closure_rows) != 320 or len(closure) != 320:
        raise ValueError("closure cardinality is not 320")

    max_density_error = 0.0
    max_angle_error = 0.0
    max_residual_error = 0.0
    valid_count = 0
    for key, row in transport.items():
        for carrier in ("electron", "hole"):
            result = closure[key + (carrier,)]
            prefix = "electron" if carrier == "electron" else "hole"
            jx = float(row[f"{prefix}_current_x_A_per_m2"])
            jy = float(row[f"{prefix}_current_y_A_per_m2"])
            gx = float(row[f"{prefix}_grad_qf_x_V_per_m"])
            gy = float(row[f"{prefix}_grad_qf_y_V_per_m"])
            mobility = float(row[f"{prefix}_mobility_m2_per_Vs"])
            grad2 = gx * gx + gy * gy
            current2 = jx * jx + jy * jy
            if mobility <= 0.0 or grad2 == 0.0 or current2 == 0.0:
                if result["status"] != "degenerate":
                    raise ValueError(f"wrong degenerate classification {key}")
                continue
            dot = jx * gx + jy * gy
            density = dot / (Q_C * mobility * grad2)
            expected_status = "valid" if density > 0.0 else "sign_incompatible"
            if result["status"] != expected_status:
                raise ValueError(f"wrong sign classification {key}")
            cosine = max(-1.0, min(1.0, dot / math.sqrt(current2 * grad2)))
            angle = math.degrees(math.acos(cosine))
            projected_x = Q_C * mobility * density * gx
            projected_y = Q_C * mobility * density * gy
            residual = math.hypot(jx - projected_x, jy - projected_y)
            residual /= math.sqrt(current2)
            angle_error = abs(
                float(result["current_grad_angle_deg"]) - angle
            )
            residual_error = abs(
                float(result["orthogonal_current_residual"]) - residual
            )
            max_angle_error = max(max_angle_error, angle_error)
            max_residual_error = max(max_residual_error, residual_error)
            if expected_status == "valid":
                density_error = abs(
                    float(result["effective_density_m3"]) - density
                )
                density_error /= max(abs(density), 1.0)
                max_density_error = max(max_density_error, density_error)
                valid_count += 1
                if not close(
                    float(result["effective_density_m3"]), density
                ):
                    raise ValueError(f"density reconstruction mismatch {key}")
            if not close(
                float(result["current_grad_angle_deg"]), angle
            ) or not close(
                float(result["orthogonal_current_residual"]), residual
            ):
                raise ValueError(f"vector reconstruction mismatch {key}")

    verification = {
        "schema_version": 1,
        "status": "valid",
        "experiment": "minimal6_transport_element_closure_independent_verify",
        "transport_sample_count": len(transport),
        "closure_sample_count": len(closure),
        "valid_density_count": valid_count,
        "maximum_relative_density_reconstruction_error": max_density_error,
        "maximum_angle_reconstruction_error_deg": max_angle_error,
        "maximum_orthogonal_residual_reconstruction_error": (
            max_residual_error
        ),
        "inputs": {
            "transport_csv_sha256": sha256(transport_csv),
            "transport_manifest_sha256": sha256(transport_manifest),
            "closure_manifest_sha256": sha256(closure_meta_path),
        },
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(verification, indent=2) + "\n", encoding="utf-8"
    )
    return verification


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--transport-csv", type=Path, required=True)
    parser.add_argument("--transport-manifest", type=Path, required=True)
    parser.add_argument("--closure-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = verify(
        transport_csv=args.transport_csv,
        transport_manifest=args.transport_manifest,
        closure_root=args.closure_root,
        output=args.output,
    )
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
