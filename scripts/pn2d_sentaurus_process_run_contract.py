#!/usr/bin/env python3
"""Manifest contract for exact Sentaurus PN2D process-probe runs."""

from __future__ import annotations

import hashlib
import json
import re
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from scripts.pn2d_high_bias_process_contract import SENTAURUS_RELEASE


SCHEMA_ID = "vela.pn2d_sentaurus_process_run.v1"
REQUIRED_FIELDS = frozenset(
    {
        "schema",
        "status",
        "experiment",
        "sentaurus_release",
        "exact_biases_V",
        "observed_biases_V",
        "variant",
        "remote_root",
        "bundle_sha256",
        "output_sha256",
    }
)
HASH_RE = re.compile(r"^[0-9a-f]{64}$")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def file_hashes(directory: Path) -> dict[str, str]:
    files = sorted(path for path in directory.iterdir() if path.is_file())
    if not files:
        raise ValueError(f"empty evidence directory: {directory}")
    return {path.name: sha256(path) for path in files}


def build_run_manifest(
    *,
    status: str,
    experiment: str,
    variant: str,
    exact_biases: Sequence[float],
    observed_biases: Sequence[float],
    remote_root: str,
    bundle: Path,
    fetched: Path,
) -> dict[str, Any]:
    return {
        "schema": SCHEMA_ID,
        "status": status,
        "experiment": experiment,
        "sentaurus_release": SENTAURUS_RELEASE,
        "exact_biases_V": [float(value) for value in exact_biases],
        "observed_biases_V": [float(value) for value in observed_biases],
        "variant": variant,
        "remote_root": remote_root,
        "bundle_sha256": file_hashes(bundle),
        "output_sha256": file_hashes(fetched),
    }


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def _validate_hash_map(value: Any, label: str) -> None:
    _require(isinstance(value, Mapping) and bool(value), f"missing {label}")
    for name, digest in value.items():
        _require(isinstance(name, str) and bool(name), f"invalid {label} name")
        _require(
            isinstance(digest, str) and HASH_RE.fullmatch(digest) is not None,
            f"invalid {label} SHA-256 for {name}",
        )


def validate_run_manifest(
    manifest: Mapping[str, Any],
    *,
    experiment: str,
    variant: str,
    exact_biases: Sequence[float],
) -> None:
    expected_biases = tuple(float(value) for value in exact_biases)
    _require(manifest.get("schema") == SCHEMA_ID, "run manifest schema mismatch")
    _require(
        set(manifest) == REQUIRED_FIELDS,
        "run manifest field set mismatch",
    )
    _require(manifest.get("status") == "passed", "run manifest did not pass")
    _require(manifest.get("experiment") == experiment, "experiment mismatch")
    _require(
        manifest.get("sentaurus_release") == SENTAURUS_RELEASE,
        "Sentaurus release mismatch",
    )
    _require(manifest.get("variant") == variant, "variant mismatch")
    _require(
        tuple(float(value) for value in manifest.get("exact_biases_V", ()))
        == expected_biases,
        "declared exact lattice mismatch",
    )
    _require(
        tuple(float(value) for value in manifest.get("observed_biases_V", ()))
        == expected_biases,
        "observed exact lattice mismatch",
    )
    _require(
        isinstance(manifest.get("remote_root"), str)
        and bool(manifest["remote_root"]),
        "missing remote root",
    )
    _validate_hash_map(manifest.get("bundle_sha256"), "bundle hashes")
    _validate_hash_map(manifest.get("output_sha256"), "output hashes")


def validate_case(
    case: Path,
    *,
    experiment: str,
    variant: str,
    exact_biases: Sequence[float],
) -> dict[str, Any]:
    manifest_path = case / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="ascii"))
    validate_run_manifest(
        manifest,
        experiment=experiment,
        variant=variant,
        exact_biases=exact_biases,
    )
    for directory_name, manifest_key in (
        ("bundle", "bundle_sha256"),
        ("fetched", "output_sha256"),
    ):
        actual = file_hashes(case / directory_name)
        _require(
            actual == dict(manifest[manifest_key]),
            f"{case}: {directory_name} hash closure mismatch",
        )
    return manifest
