#!/usr/bin/env python3
"""Seal the forty-state Minimal6 box staged replay evidence."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def copy(source: Path, destination: Path) -> None:
    if not source.is_file():
        raise FileNotFoundError(source)
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, destination)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, required=True)
    parser.add_argument("--root", type=Path, required=True)
    args = parser.parse_args()
    repo = args.repo.resolve()
    root = args.root.resolve()
    raw = root / "raw"

    fixed_inputs = {
        "inputs/observations_node.csv": repo
        / "build-release/pn2d-minimal6-physics-inverse-audit-unitfix-20260723-b/observations_node.csv",
        "inputs/transport_element_values.csv": repo
        / "build-release/pn2d-minimal6-transport-elements-20260723-b/export/transport_element_values.csv",
        "inputs/vela_self_consistent_edge_samples.csv": repo
        / "build-release/pn2d-minimal6-self-consistent-replacement-velocityfix-20260724-a/self_consistent_edge_samples.csv",
        "inputs/mirror_mesh.json": repo
        / "build-release/pn2d-minimal6-inverse-inputs-unitfix-20260723-b/vela/source/topologies/mirror/mesh.json",
        "inputs/sketch_mesh.json": repo
        / "build-release/pn2d-minimal6-inverse-inputs-unitfix-20260723-b/vela/source/topologies/sketch/mesh.json",
        "inputs/single_state_closure_summary.json": repo
        / "build-release/pn2d-minimal6-sentaurus-box-current-replay-20260724-a/closure_summary.json",
        "scripts/box_staged_sweep.py": repo
        / "scripts/pn2d_minimal6_diagnostics/box_staged_sweep.py",
        "scripts/diagnose_pn2d_minimal6_box_staged_sweep.py": repo
        / "scripts/diagnose_pn2d_minimal6_box_staged_sweep.py",
        "scripts/remap_pn2d_minimal6_transport_cells.py": repo
        / "scripts/remap_pn2d_minimal6_transport_cells.py",
        "scripts/verify_pn2d_minimal6_box_staged_sweep.py": repo
        / "scripts/verify_pn2d_minimal6_box_staged_sweep.py",
        "scripts/seal_pn2d_minimal6_box_staged_sweep.py": Path(__file__).resolve(),
        "reports/pn2d_minimal6_sentaurus_box_staged_sweep_2026-07-24.md": repo
        / "docs/validation/pn2d_minimal6_sentaurus_box_staged_sweep_2026-07-24.md",
    }
    for relative, source in fixed_inputs.items():
        copy(source, raw / relative)

    sentaurus_root = (
        repo
        / "build-release/pn2d-minimal6-transport-elements-20260723-b/codex_pn2d_minimal6_transport_elements_20260723_b"
    )
    vela_root = (
        repo
        / "build-release/pn2d-minimal6-self-consistent-replacement-velocityfix-20260724-a/baseline_replay"
    )
    for topology in ("mirror", "sketch"):
        for magnitude in range(1, 21):
            label = f"m{magnitude}V"
            copy(
                sentaurus_root
                / topology
                / label
                / f"pn2d_minimal6_state_{label}.plt",
                raw
                / "terminal"
                / topology
                / label
                / f"pn2d_minimal6_state_{label}.plt",
            )
            copy(
                vela_root / topology / label / "edges.csv",
                raw / "vela_edges" / topology / label / "edges.csv",
            )

    files = []
    for path in sorted(root.rglob("*")):
        if not path.is_file() or path.name == "evidence_manifest.json":
            continue
        files.append(
            {
                "path": path.relative_to(root).as_posix(),
                "sha256": sha256(path),
                "size_bytes": path.stat().st_size,
            }
        )
    manifest = {
        "schema_version": 1,
        "status": "sealed",
        "artifact": root.name,
        "state_count": 40,
        "file_count": len(files),
        "files": files,
    }
    (root / "evidence_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({"status": "sealed", "file_count": len(files)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
