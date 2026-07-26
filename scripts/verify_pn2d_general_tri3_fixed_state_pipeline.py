#!/usr/bin/env python3
"""Independently verify deterministic general-Tri3 Tasks 3-5 outputs."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any


RAW_SCHEMA = "pn2d_general_tri3_element_edge_avalanche/v1"
STAGES = {
    "imported_state": "pn2d_general_tri3_imported_state/v1",
    "element_edge_current": "pn2d_general_tri3_element_edge_current/v1",
    "current_closure": (
        "pn2d_general_tri3_element_edge_current_closure/v1"
    ),
}
SOURCE_SCHEMA = "pn2d_general_tri3_element_edge_source/v1"
SOURCE_METHODS = ("gss", "charon", "genius", "active")
EXACT_BIASES_V = [-1.0, -10.0, -20.0]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root-a", type=Path, required=True)
    parser.add_argument("--root-b", type=Path, required=True)
    parser.add_argument(
        "--scientific-role",
        choices=("device_physics_oracle", "diagnostic_only"),
        required=True,
    )
    parser.add_argument("--skip-source", action="store_true")
    return parser.parse_args()


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def load(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="ascii"))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def verify_manifest_outputs(
    root: Path,
    manifest: dict[str, Any],
) -> None:
    for name, expected in manifest["outputs"].items():
        path = root / name
        require(path.is_file(), f"missing output: {path}")
        require(digest(path) == expected, f"output hash mismatch: {path}")


def verify_stage_pair(
    root_a: Path,
    root_b: Path,
    stage: str,
    schema: str,
) -> dict[str, Any]:
    manifests = []
    for root in (root_a, root_b):
        stage_root = root / stage
        manifest = load(stage_root / "analysis_manifest.json")
        require(manifest["schema"] == schema, f"{stage} schema mismatch")
        require(manifest["status"] == "valid", f"{stage} is not valid")
        require(
            manifest["exact_biases_V"] == EXACT_BIASES_V,
            f"{stage} bias contract mismatch",
        )
        verify_manifest_outputs(stage_root, manifest)
        manifests.append(manifest)
    require(
        manifests[0]["case_name"] == manifests[1]["case_name"],
        f"{stage} case mismatch",
    )
    require(
        manifests[0]["outputs"] == manifests[1]["outputs"],
        f"{stage} A/B outputs differ",
    )
    return manifests[0]


def main() -> int:
    args = parse_args()
    roots = (args.root_a.resolve(), args.root_b.resolve())
    raw_manifests = [load(root / "manifest.json") for root in roots]
    for manifest in raw_manifests:
        require(manifest["schema"] == RAW_SCHEMA, "raw schema mismatch")
        require(manifest["status"] == "passed", "raw root is not passed")
        require(
            manifest["exact_biases_V"] == EXACT_BIASES_V,
            "raw bias contract mismatch",
        )
    require(
        sorted(raw_manifests[0]["cases"])
        == sorted(raw_manifests[1]["cases"]),
        "raw A/B case mismatch",
    )

    verified: dict[str, Any] = {}
    for stage, schema in STAGES.items():
        verified[stage] = verify_stage_pair(
            roots[0],
            roots[1],
            stage,
            schema,
        )

    imported = verified["imported_state"]
    for carrier in ("electron", "hole"):
        require(
            imported["density"][carrier]["maximum"] <= 1.0e-4,
            f"{carrier} density replay exceeds fixed-state gate",
        )
    for field in (
        "electric_field",
        "electron_qfp_gradient",
        "hole_qfp_gradient",
    ):
        require(
            imported["vectors"][field]["active_relative"]["maximum"]
            <= 1.0e-9,
            f"{field} replay exceeds fixed-state gate",
        )

    closure = verified["current_closure"]
    require(
        closure["native_element_current_observation"]
        == "insufficient_native_observation_undocumented_element_vector",
        "native element-current observation must remain typed insufficient",
    )
    require(
        closure["geometry_support"][
            "maximum_sentaurus_area_relative_error"
        ]
        <= 1.0e-12,
        "Sentaurus geometry support does not close",
    )
    if args.scientific_role == "device_physics_oracle":
        require(
            closure["geometry_support"][
                "maximum_vela_truncated_area_relative_error"
            ]
            <= 1.0e-12,
            "Vela geometry support does not close on scientific oracle",
        )

    source_methods: dict[str, Any] = {}
    if not args.skip_source:
        for method in SOURCE_METHODS:
            stage = f"element_edge_source_{method}"
            manifest = verify_stage_pair(
                roots[0],
                roots[1],
                stage,
                SOURCE_SCHEMA,
            )
            require(
                manifest["source_identity"][
                    "maximum_readmeasure_currentplot_relative_error"
                ]
                <= 1.0e-10,
                f"{method} ReadMeasure/CurrentPlot identity failed",
            )
            source_methods[method] = {
                "current_vector_method": manifest[
                    "current_vector_method"
                ],
                "source_identity": manifest["source_identity"],
            }

    result = {
        "schema": "pn2d_general_tri3_fixed_state_verification/v1",
        "status": "passed",
        "case_name": verified["imported_state"]["case_name"],
        "scientific_role": args.scientific_role,
        "a_b_outputs": "exact",
        "verified_stages": list(STAGES),
        "source_methods": source_methods,
    }
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
