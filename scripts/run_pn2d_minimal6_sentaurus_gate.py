#!/usr/bin/env python3
"""Build and run the PN2D minimal6 explicit-grid Sentaurus gate."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import shutil
import sys
from datetime import datetime
from pathlib import Path
from typing import Sequence


REPO = Path(__file__).resolve().parents[1]
if __package__ in (None, ""):
    sys.path.insert(0, str(REPO))

from scripts.pn2d_minimal6_topology import (  # noqa: E402
    Topology,
    canonical_edges,
    canonical_triangle,
    load_topology,
    validate_dfise_roundtrip,
    validate_topology,
    write_dfise_doping,
    write_dfise_grid,
)
from scripts.run_sentaurus_vm_reference import (  # noqa: E402
    default_windows_openssh,
    run_checked,
    write_manifest,
)


REFERENCE_ROOT = REPO / "reference_tcad" / "pn2d_sentaurus2018_minimal6"
SOURCE_DIR = REFERENCE_ROOT / "source"
REFERENCE_DESCRIPTOR = REFERENCE_ROOT / "pn2d_sentaurus2018_minimal6_reference.json"
TOPOLOGY_FIXTURE = SOURCE_DIR / "minimal6_topologies.json"
MODELS_SOURCE = SOURCE_DIR / "models.par"
DECK_SOURCE = SOURCE_DIR / "pn2d_minimal6_gate_sdevice.cmd"
DEFAULT_OUTPUT_DIR = (
    REPO
    / "build-release"
    / "reference_tcad"
    / "pn2d_sentaurus2018_minimal6"
    / "sentaurus_gate_runs"
)
DEFAULT_REMOTE_ROOT = "~/sentaurus_runs/vela_oracle"
DEFAULT_IMPORTER = REPO / "build-release" / (
    "sentaurus_import.exe" if os.name == "nt" else "sentaurus_import"
)
STAGED_FILES = (
    "pn2d_minimal6.grd",
    "pn2d_minimal6.dat",
    "pn2d_minimal6_gate_sdevice.cmd",
    "models.par",
)
RETURNED_FILES = (
    "pn2d_minimal6_gate_des.tdr",
    "pn2d_minimal6_gate_des.log",
    "pn2d_minimal6.tdr",
    "pn2d_minimal6.grd",
    "pn2d_minimal6.dat",
    "run_tdx_dfise_to_tdr.out",
    "run_pn2d_minimal6_gate.out",
)
COORDINATE_TOLERANCE_UM = 1.0e-12
DOPING_RELATIVE_TOLERANCE = 1.0e-15
_SAFE_REMOTE_COMPONENT = re.compile(r"^[A-Za-z0-9_./~:-]+$")
_SAFE_RUN_ID = re.compile(r"^[A-Za-z0-9_.-]+$")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _load_descriptor() -> dict[str, object]:
    descriptor = json.loads(REFERENCE_DESCRIPTOR.read_text(encoding="utf-8"))
    if descriptor.get("schema") != "vela.pn2d_minimal6_reference.v1":
        raise ValueError("unexpected PN2D minimal6 reference descriptor schema")
    if descriptor.get("diagnostic_only") is not True:
        raise ValueError("PN2D minimal6 reference must remain diagnostic-only")
    if descriptor.get("coordinate_tolerance_um") != COORDINATE_TOLERANCE_UM:
        raise ValueError("unexpected PN2D minimal6 coordinate tolerance")
    return descriptor


def _validate_remote_component(value: str, name: str) -> str:
    if not value or not _SAFE_REMOTE_COMPONENT.fullmatch(value):
        raise ValueError(f"{name} contains unsupported remote-shell characters")
    return value.rstrip("/")


def _remote_dir(remote_root: str, run_id: str, topology_id: str) -> str:
    root = _validate_remote_component(remote_root, "remote root")
    if not _SAFE_RUN_ID.fullmatch(run_id):
        raise ValueError("run ID contains unsupported characters")
    return f"{root}/{run_id}/{topology_id}"


def build_gate_bundle(topology_id: str, output_dir: Path) -> dict[str, object]:
    descriptor = _load_descriptor()
    approved = tuple(str(value) for value in descriptor["topologies"])
    if topology_id not in approved:
        raise ValueError(f"unsupported minimal6 topology: {topology_id}")
    for required in (TOPOLOGY_FIXTURE, MODELS_SOURCE, DECK_SOURCE):
        if not required.is_file():
            raise FileNotFoundError(f"missing minimal6 gate source file: {required}")

    topology = load_topology(TOPOLOGY_FIXTURE, topology_id)
    summary = validate_topology(topology)
    output_dir.mkdir(parents=True, exist_ok=True)
    grid = output_dir / STAGED_FILES[0]
    doping = output_dir / STAGED_FILES[1]
    deck = output_dir / STAGED_FILES[2]
    models = output_dir / STAGED_FILES[3]
    write_dfise_grid(topology, grid)
    write_dfise_doping(topology, doping)
    shutil.copyfile(DECK_SOURCE, deck)
    shutil.copyfile(MODELS_SOURCE, models)
    roundtrip = validate_dfise_roundtrip(topology, grid, doping)
    if not roundtrip["passed"]:
        raise ValueError(f"generated {topology_id} DF-ISE bundle failed roundtrip")

    staged_paths = [grid, doping, deck, models]
    return {
        "topology_id": topology_id,
        "bundle_dir": str(output_dir.resolve()),
        "staged_files": [path.name for path in staged_paths],
        "file_sha256": {path.name: sha256_file(path) for path in staged_paths},
        "topology_fixture_sha256": sha256_file(TOPOLOGY_FIXTURE),
        "topology_contract": {
            "nodes": summary.nodes,
            "triangles": summary.triangles,
            "edges": summary.edges,
            "contact_edges": {
                name: list(edge) for name, edge in summary.contact_edges.items()
            },
            "triangle_connectivity": [list(triangle) for triangle in topology.triangles],
        },
    }


def prepare_gate(
    *,
    topology_ids: Sequence[str],
    run_id: str,
    output_dir: Path,
    ssh_target: str,
    remote_root: str = DEFAULT_REMOTE_ROOT,
) -> dict[str, object]:
    descriptor = _load_descriptor()
    requested = list(topology_ids)
    if not requested or len(set(requested)) != len(requested):
        raise ValueError("topology IDs must be non-empty and unique")
    approved = [str(value) for value in descriptor["topologies"]]
    invalid = [value for value in requested if value not in approved]
    if invalid:
        raise ValueError(f"unsupported topology IDs: {', '.join(invalid)}")
    if not _SAFE_RUN_ID.fullmatch(run_id):
        raise ValueError("run ID contains unsupported characters")

    run_root = (output_dir / run_id).resolve()
    runs: list[dict[str, object]] = []
    for topology_id in requested:
        bundle_dir = run_root / "topologies" / topology_id / "source"
        artifacts_dir = run_root / "topologies" / topology_id / "artifacts"
        run = build_gate_bundle(topology_id, bundle_dir)
        remote_dir = _remote_dir(remote_root, run_id, topology_id)
        run.update(
            {
                "remote_dir": remote_dir,
                "artifacts_dir": str(artifacts_dir),
                "remote_commands": [
                    f"cd {remote_dir} && "
                    "tdx -d pn2d_minimal6.grd pn2d_minimal6.dat "
                    "pn2d_minimal6.tdr > run_tdx_dfise_to_tdr.out 2>&1",
                    f"cd {remote_dir} && "
                    "sdevice pn2d_minimal6_gate_sdevice.cmd "
                    "> run_pn2d_minimal6_gate.out 2>&1",
                ],
                "status": "prepared",
            }
        )
        runs.append(run)

    manifest: dict[str, object] = {
        "schema": "vela.pn2d_minimal6_sentaurus_gate.v1",
        "diagnostic_only": True,
        "run_id": run_id,
        "dry_run": True,
        "ssh_target": ssh_target,
        "remote_root": remote_root,
        "coordinate_tolerance_um": COORDINATE_TOLERANCE_UM,
        "doping_relative_tolerance": DOPING_RELATIVE_TOLERANCE,
        "reference_descriptor": str(REFERENCE_DESCRIPTOR.resolve()),
        "reference_descriptor_sha256": sha256_file(REFERENCE_DESCRIPTOR),
        "run_root": str(run_root),
        "runs": runs,
        "passed": None,
    }
    manifest_path = run_root / "manifest.json"
    manifest["manifest_path"] = str(manifest_path)
    write_manifest(manifest_path, manifest)
    return manifest


def build_live_argv(
    run: dict[str, object],
    *,
    ssh_bin: str,
    scp_bin: str,
    ssh_target: str,
) -> list[list[str]]:
    bundle_dir = Path(str(run["bundle_dir"]))
    artifacts_dir = Path(str(run["artifacts_dir"]))
    remote_dir = str(run["remote_dir"])
    commands: list[list[str]] = [
        [ssh_bin, ssh_target, f"mkdir -p {remote_dir}"],
    ]
    for name in run["staged_files"]:
        commands.append(
            [scp_bin, str(bundle_dir / str(name)), f"{ssh_target}:{remote_dir}/"]
        )
    for remote_command in run["remote_commands"]:
        commands.append([ssh_bin, ssh_target, str(remote_command)])
    for name in RETURNED_FILES:
        commands.append(
            [
                scp_bin,
                f"{ssh_target}:{remote_dir}/{name}",
                str(artifacts_dir) + os.sep,
            ]
        )
    return commands


def _read_csv(path: Path, required: set[str]) -> list[dict[str, str]]:
    if not path.is_file():
        raise ValueError(f"missing neutral export file: {path.name}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fields = set(reader.fieldnames or [])
        if not required.issubset(fields):
            missing = ", ".join(sorted(required - fields))
            raise ValueError(f"{path.name} missing columns: {missing}")
        return list(reader)


def _canonical_node_mapping(
    topology: Topology,
    rows: list[dict[str, str]],
) -> dict[int, int]:
    if len(rows) != 6:
        raise ValueError(f"expected 6 nodes, found {len(rows)}")
    mapping: dict[int, int] = {}
    used_labels: set[int] = set()
    for row in rows:
        source_id = int(row["id"])
        if source_id in mapping:
            raise ValueError(f"duplicate neutral node ID: {source_id}")
        point = (float(row["x_um"]), float(row["y_um"]))
        matches = [
            label
            for label, expected in topology.nodes.items()
            if abs(point[0] - expected[0]) < COORDINATE_TOLERANCE_UM
            and abs(point[1] - expected[1]) < COORDINATE_TOLERANCE_UM
        ]
        if len(matches) != 1:
            raise ValueError(
                f"node {source_id} has no unique exact canonical coordinate match"
            )
        label = matches[0]
        if label in used_labels:
            raise ValueError(f"duplicate canonical coordinate for node {label}")
        mapping[source_id] = label
        used_labels.add(label)
    if used_labels != set(topology.nodes):
        raise ValueError("canonical coordinate mapping is incomplete")
    return mapping


def _triangle_area2(topology: Topology, triangle: tuple[int, int, int]) -> float:
    (ax, ay), (bx, by), (cx, cy) = (
        topology.nodes[node_id] for node_id in triangle
    )
    return (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)


def validate_returned_tdr(
    topology: Topology,
    neutral_export: Path,
) -> dict[str, object]:
    validate_topology(topology)
    node_rows = _read_csv(neutral_export / "nodes.csv", {"id", "x_um", "y_um"})
    node_mapping = _canonical_node_mapping(topology, node_rows)

    element_rows = _read_csv(
        neutral_export / "elements.csv",
        {"id", "node0", "node1", "node2", "region", "material"},
    )
    if len(element_rows) != 4:
        raise ValueError(f"expected 4 triangles, found {len(element_rows)}")
    triangles: list[tuple[int, int, int]] = []
    for row in element_rows:
        if row["region"] != "R.Si" or row["material"] != "Si":
            raise ValueError("triangle region/material mismatch")
        try:
            triangle = tuple(
                node_mapping[int(row[name])] for name in ("node0", "node1", "node2")
            )
        except KeyError as error:
            raise ValueError("triangle references an unknown node") from error
        if _triangle_area2(topology, triangle) <= 0.0:
            raise ValueError("triangle orientation mismatch")
        triangles.append(canonical_triangle(triangle))
    expected_triangles = {canonical_triangle(value) for value in topology.triangles}
    if len(set(triangles)) != 4 or set(triangles) != expected_triangles:
        raise ValueError("triangle connectivity mismatch")
    edges = canonical_edges(triangles)
    if len(edges) != 9:
        raise ValueError(f"expected 9 unique edges, found {len(edges)}")

    contact_rows = _read_csv(
        neutral_export / "contacts.csv", {"name", "node_ids", "region"}
    )
    contacts: dict[str, tuple[int, int]] = {}
    for row in contact_rows:
        if row["name"] in contacts or row["region"] != "R.Si":
            raise ValueError("contact edges mismatch")
        try:
            labels = tuple(
                sorted(node_mapping[int(value)] for value in row["node_ids"].split(";"))
            )
        except (KeyError, ValueError) as error:
            raise ValueError("contact edges reference an unknown node") from error
        if len(labels) != 2:
            raise ValueError("contact edges must contain exactly two nodes")
        contacts[row["name"]] = labels
    expected_contacts = {
        name: tuple(sorted(edge)) for name, edge in topology.contacts.items()
    }
    if contacts != expected_contacts:
        raise ValueError("contact edges mismatch")

    doping_rows = _read_csv(
        neutral_export / "doping.csv",
        {"node_id", "donors_cm3", "acceptors_cm3"},
    )
    if len(doping_rows) != 6:
        raise ValueError(f"expected 6 doping rows, found {len(doping_rows)}")
    doping: dict[int, tuple[float, float]] = {}
    for row in doping_rows:
        try:
            label = node_mapping[int(row["node_id"])]
        except KeyError as error:
            raise ValueError("doping references an unknown node") from error
        if label in doping:
            raise ValueError("duplicate doping node")
        doping[label] = (float(row["donors_cm3"]), float(row["acceptors_cm3"]))
    expected_doping = {
        label: (topology.donors_cm3[label], topology.acceptors_cm3[label])
        for label in topology.nodes
    }
    if set(doping) != set(expected_doping) or any(
        not math.isclose(
            actual,
            expected,
            rel_tol=DOPING_RELATIVE_TOLERANCE,
            abs_tol=0.0,
        )
        for label in expected_doping
        for actual, expected in zip(doping[label], expected_doping[label])
    ):
        raise ValueError("doping mismatch")

    metadata_path = neutral_export / "metadata.json"
    if not metadata_path.is_file():
        raise ValueError("missing neutral export file: metadata.json")
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if int(metadata.get("vertex_count", -1)) != 6:
        raise ValueError("metadata vertex count mismatch")

    return {
        "schema": "vela.pn2d_minimal6_topology_gate.v1",
        "topology_id": topology.topology_id,
        "passed": True,
        "coordinate_tolerance_um": COORDINATE_TOLERANCE_UM,
        "doping_relative_tolerance": DOPING_RELATIVE_TOLERANCE,
        "node_count": len(node_rows),
        "triangle_count": len(triangles),
        "edge_count": len(edges),
        "triangles": [list(value) for value in sorted(triangles)],
        "contact_edges": {name: list(edge) for name, edge in contacts.items()},
        "doping_matches": True,
    }


def run_live(
    manifest: dict[str, object],
    *,
    ssh_bin: str,
    scp_bin: str,
    importer: Path,
) -> None:
    if not importer.is_file():
        raise FileNotFoundError(f"Sentaurus importer is not built: {importer}")
    ssh_target = str(manifest["ssh_target"])
    for run in manifest["runs"]:
        artifacts_dir = Path(str(run["artifacts_dir"]))
        artifacts_dir.mkdir(parents=True, exist_ok=True)
        argv_commands = build_live_argv(
            run,
            ssh_bin=ssh_bin,
            scp_bin=scp_bin,
            ssh_target=ssh_target,
        )
        run["local_argv"] = argv_commands
        execution_count = (
            1 + len(run["staged_files"]) + len(run["remote_commands"])
        )
        try:
            for argv in argv_commands[:execution_count]:
                run_checked(argv)
        except Exception as execution_error:
            recovery_errors: list[str] = []
            for argv in argv_commands[execution_count:]:
                try:
                    run_checked(argv)
                except Exception as recovery_error:
                    recovery_errors.append(str(recovery_error))
            run["status"] = "failed"
            run["execution_error"] = str(execution_error)
            run["artifact_recovery_errors"] = recovery_errors
            raise
        for argv in argv_commands[execution_count:]:
            run_checked(argv)

        neutral_export = artifacts_dir / "neutral_export"
        tdr = artifacts_dir / "pn2d_minimal6_gate_des.tdr"
        run_checked(
            [
                str(importer),
                "--tdr",
                str(tdr),
                "--export-dir",
                str(neutral_export),
                "--compensated-doping-policy",
                "reported",
            ]
        )
        topology = load_topology(TOPOLOGY_FIXTURE, str(run["topology_id"]))
        gate_report = validate_returned_tdr(topology, neutral_export)
        gate_path = artifacts_dir / "topology_gate.json"
        write_manifest(gate_path, gate_report)
        run["topology_gate"] = str(gate_path)
        run["returned_file_sha256"] = {
            name: sha256_file(artifacts_dir / name) for name in RETURNED_FILES
        }
        run["status"] = "passed"


def _parse_topologies(raw: str) -> tuple[str, ...]:
    values = tuple(value.strip() for value in raw.split(",") if value.strip())
    if not values:
        raise ValueError("--topologies must contain at least one topology")
    return values


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topologies", default="sketch,mirror")
    parser.add_argument("--ssh-target", default="sentaurus")
    parser.add_argument("--ssh-bin", default=None)
    parser.add_argument("--scp-bin", default=None)
    parser.add_argument("--remote-root", default=DEFAULT_REMOTE_ROOT)
    parser.add_argument("--run-id", default=None)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--importer", type=Path, default=DEFAULT_IMPORTER)
    parser.add_argument("--dry-run", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    manifest: dict[str, object] | None = None
    try:
        run_id = args.run_id or datetime.now().strftime("minimal6_gate_%Y%m%d_%H%M%S")
        manifest = prepare_gate(
            topology_ids=_parse_topologies(args.topologies),
            run_id=run_id,
            output_dir=args.output_dir,
            ssh_target=args.ssh_target,
            remote_root=args.remote_root,
        )
        manifest["dry_run"] = bool(args.dry_run)
        if not args.dry_run:
            run_live(
                manifest,
                ssh_bin=args.ssh_bin or default_windows_openssh("ssh"),
                scp_bin=args.scp_bin or default_windows_openssh("scp"),
                importer=args.importer.resolve(),
            )
            manifest["passed"] = True
        write_manifest(Path(str(manifest["manifest_path"])), manifest)
        print(json.dumps(manifest, indent=2))
        return 0
    except Exception as error:  # noqa: BLE001 - preserve a failed live manifest.
        if manifest is not None:
            manifest["passed"] = False
            manifest["error"] = str(error)
            write_manifest(Path(str(manifest["manifest_path"])), manifest)
        print(str(error), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
