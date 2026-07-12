#!/usr/bin/env python3
"""Reproduce the current-HEAD PN2D compensated-junction SG replay artifacts."""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence


REPO = Path(__file__).resolve().parents[1]
CASE_NAME = "pn2d_sentaurus2018_coarse7x3"
DEFAULT_SOURCE_DIR = (
    REPO
    / "build-release"
    / "reference_tcad"
    / CASE_NAME
    / "sentaurus_vm_runs"
    / "pn2d_coarse7x3_vector_current_20260630"
    / "source"
)
DEFAULT_REFERENCE_CONFIG = (
    REPO
    / "reference_tcad"
    / CASE_NAME
    / f"{CASE_NAME}_reference.json"
)
DEFAULT_OUT_DIR = (
    REPO
    / "build-release"
    / "reference_tcad"
    / CASE_NAME
    / "reports"
    / "pn2d_bv_compensated_sg_replay"
)
DEFAULT_RUNNER = REPO / "build-release" / (
    "vela_example_runner.exe" if os.name == "nt" else "vela_example_runner"
)
DEFAULT_IMPORTER = REPO / "build-release" / (
    "sentaurus_import.exe" if os.name == "nt" else "sentaurus_import"
)
DEFAULT_IMPORT_SCRIPT = REPO / "scripts" / "sentaurus_import.py"
DEFAULT_DIAGNOSTIC_SCRIPT = (
    REPO / "scripts" / "diagnose_pn2d_bv_compensated_source_proxy.py"
)
DEFAULT_BIASES = [-12.0, -19.0, -20.0]
BIAS_TO_MULTIBIAS_INDEX = {
    -12.0: 240,
    -19.0: 380,
    -20.0: 400,
}
REQUIRED_SCALAR_FIELD_SPECS = {
    "ElectrostaticPotential": "V",
    "eQuasiFermiPotential": "V",
    "eDensity": "cm^-3",
    "eMobility": "cm^2*V^-1*s^-1",
    "eAlphaAvalanche": "cm^-1",
}
ARTIFACT_MANIFEST_NAME = "artifact_manifest.json"
MANIFEST_SCHEMA = "vela.pn2d_compensated_sg_replay.artifact_manifest.v3"
IMPACT_IONIZATION_BASE = {
    "model": "van_overstraeten",
    "driving_force": "quasi_fermi_gradient",
    "generation": "current_density",
    "current_magnitude_mode": "edge_scalar_abs",
    "cell_reconstructed_midpoint_density": "bernoulli",
    "quasi_fermi_gradient_discretization": "edge_difference",
    "source_volume_policy": "genius_truncated",
    "source_volume_factor": 0.0,
    "source_geometry_scale": 1.0,
    "edge_source_partition": "symmetric",
}
LEGACY_DOPING_REPLAY = {
    "strategy": "legacy_p_side_unresolved_compensated",
    "source_policy": "dominant_signed_region",
    "resolution_source": "signed_aggregate_zero",
}


def multibias_index_for_bias(bias: float) -> int:
    """Return the audited coarse sweep index for one supported replay bias."""
    normalized = float(bias)
    try:
        return BIAS_TO_MULTIBIAS_INDEX[normalized]
    except KeyError as exc:
        supported = ", ".join(f"{item:g}" for item in DEFAULT_BIASES)
        raise ValueError(
            f"unsupported replay bias {bias:g} V; expected one of {supported}"
        ) from exc


def required_tdr_path(source_dir: Path, bias: float) -> Path:
    """Resolve a required multibias TDR and fail before any export can run."""
    index = multibias_index_for_bias(bias)
    path = Path(source_dir) / f"pn2d_bv_multibias_{index:04d}_des.tdr"
    if not path.is_file():
        raise FileNotFoundError(
            f"missing multibias TDR for {bias:g} V at index {index:04d}: {path}"
        )
    return path


def _load_manifest(
    manifest: Path | Mapping[str, Any],
) -> Mapping[str, Any]:
    if isinstance(manifest, Path):
        return json.loads(manifest.read_text(encoding="utf-8-sig"))
    return manifest


def _require_manifest_field(
    fields: Sequence[Any],
    name: str,
    *,
    components: int,
    unit: str,
) -> dict[str, Any]:
    matches = [
        dict(field)
        for field in fields
        if (
            isinstance(field, Mapping)
            and field.get("name") == name
            and field.get("components") == components
        )
    ]
    if len(matches) != 1:
        raise ValueError(
            f"field_manifest requires exactly one {name} components={components} entry"
        )
    field = matches[0]
    if field.get("unit") != unit:
        raise ValueError(f"{name} unit must be {unit}")
    if int(field.get("region", -1)) != 0:
        raise ValueError(f"{name} components={components} must be region0")
    if field.get("mapping_status") != "complete":
        raise ValueError(f"{name} mapping_status must be complete")
    if field.get("global_node_mapping") != "global_vertex_order":
        raise ValueError(f"{name} must use global_vertex_order mapping")
    return field


def validate_vector_field_manifest(
    manifest: Path | Mapping[str, Any],
) -> dict[str, Any]:
    """Require strict scalar replay fields and one vector electron current."""
    data = _load_manifest(manifest)
    fields = data.get("fields", [])
    if not isinstance(fields, list):
        raise ValueError("field_manifest fields must be a list")
    for name, unit in REQUIRED_SCALAR_FIELD_SPECS.items():
        _require_manifest_field(
            fields,
            name,
            components=1,
            unit=unit,
        )
    return _require_manifest_field(
        fields,
        "eCurrentDensity",
        components=2,
        unit="A*cm^-2",
    )


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _manifest_path(out_dir: Path, path: Path) -> str:
    root = out_dir.resolve()
    resolved = path.resolve()
    try:
        return resolved.relative_to(root).as_posix()
    except ValueError:
        return resolved.as_posix()


def _artifact_records(
    out_dir: Path,
    paths: Iterable[Path],
) -> list[dict[str, Any]]:
    manifest_path = (out_dir / ARTIFACT_MANIFEST_NAME).resolve()
    unique: dict[str, Path] = {}
    for raw_path in paths:
        path = Path(raw_path)
        resolved = path.resolve()
        if resolved == manifest_path:
            continue
        unique[_manifest_path(out_dir, path)] = path

    records = []
    for display_path, path in sorted(unique.items()):
        records.append({
            "path": display_path,
            "sha256": sha256_file(path),
            "size_bytes": path.stat().st_size,
        })
    return records


def _sorted_mapping(value: Mapping[str, Any] | None) -> dict[str, Any]:
    if value is None:
        return {}
    return {str(key): value[key] for key in sorted(value, key=str)}


def build_artifact_manifest(
    *,
    out_dir: Path,
    git_head: str,
    dirty: bool,
    parameters: Mapping[str, Any] | None = None,
    commands: Sequence[Mapping[str, Any]] = (),
    artifact_paths: Iterable[Path] = (),
    tdrs: Sequence[Mapping[str, Any]] = (),
    generated_decks: Iterable[Path] = (),
    model_configs: Iterable[Path] = (),
    mesh_files: Iterable[Path] = (),
    doping_files: Iterable[Path] = (),
    material_files: Iterable[Path] = (),
    edge_mapping: Any = None,
) -> dict[str, Any]:
    """Build schema-v3 data with stable ordering and no self hash."""
    normalized_commands = [
        {
            "argv": [str(item) for item in command.get("argv", [])],
            "cwd": str(command.get("cwd", "")),
            "returncode": command.get("returncode"),
        }
        for command in commands
    ]
    normalized_tdrs = sorted(
        (dict(item) for item in tdrs),
        key=lambda item: (
            float(item.get("bias_V", 0.0)),
            str(item.get("path", "")),
        ),
    )
    return {
        "schema": MANIFEST_SCHEMA,
        "schema_version": 3,
        "git_head": str(git_head),
        "dirty": bool(dirty),
        "parameters": _sorted_mapping(parameters),
        "commands": normalized_commands,
        "tdrs": normalized_tdrs,
        "generated": {
            "decks": _artifact_records(out_dir, generated_decks),
            "model_configs": _artifact_records(out_dir, model_configs),
        },
        "inputs": {
            "mesh": _artifact_records(out_dir, mesh_files),
            "doping": _artifact_records(out_dir, doping_files),
            "materials": _artifact_records(out_dir, material_files),
        },
        "edge_mapping": [] if edge_mapping is None else edge_mapping,
        "artifacts": _artifact_records(out_dir, artifact_paths),
    }


def variant_specs(out_dir: Path) -> dict[str, dict[str, Any]]:
    variants_root = Path(out_dir) / "variants"
    doping_strategies = (
        ("legacy", "dominant_signed_region"),
        ("reported", "reported"),
    )
    current_variants = (
        ("density_gradient", "density_gradient"),
        ("gss_midpoint", "cell_reconstructed"),
        ("triangle_gss_gradqf", "cell_reconstructed"),
    )
    specs: dict[str, dict[str, Any]] = {}
    for doping_strategy, policy in doping_strategies:
        for current_variant, current_approximation in current_variants:
            name = f"{doping_strategy}_{current_variant}"
            impact_ionization = {
                **IMPACT_IONIZATION_BASE,
                "current_approximation": current_approximation,
            }
            diagnostic_csv_kind = "sg_avalanche_edges"
            diagnostic_csv_name = f"sg_avalanche_edges_{current_variant}.csv"
            if current_variant == "triangle_gss_gradqf":
                impact_ionization.update({
                    "cell_reconstructed_midpoint_density": "gss_logistic",
                    "quasi_fermi_gradient_discretization": "cell_gradient",
                    "source_mapping_mode": "triangle_gss_gradqf_truncated",
                })
                diagnostic_csv_kind = "triangle_gss_sources"
                diagnostic_csv_name = (
                    f"triangle_gss_sources_{current_variant}.csv"
                )
            specs[name] = {
                "name": name,
                "implementation": "current_head",
                "doping_strategy": doping_strategy,
                "current_variant": current_variant,
                "current_approximation": current_approximation,
                "impact_ionization": impact_ionization,
                "compensated_doping_policy": policy,
                "root": variants_root / name,
                "reference_config": variants_root / name / "reference_config.json",
                "imported_dir": variants_root / name / "imported",
                "run_dir": variants_root / name / "run",
                "deck_name": f"simulation_pn2d_bv_{current_variant}.json",
                "output_csv_name": f"pn2d_bv_{current_variant}.csv",
                "vtk_subdir": f"vtk/{current_variant}",
                "diagnostics_suffix": f"_{current_variant}",
                "diagnostic_csv_kind": diagnostic_csv_kind,
                "diagnostic_csv_name": diagnostic_csv_name,
            }
    return specs


def export_command(
    *,
    importer: Path,
    tdr: Path,
    export_dir: Path,
) -> list[str]:
    return [
        str(importer),
        "--tdr",
        str(tdr),
        "--inventory-json",
        str(export_dir / "tdr_inventory.json"),
        "--export-dir",
        str(export_dir),
        "--compensated-doping-policy",
        "reported",
    ]


def reference_import_command(
    *,
    import_script: Path,
    reference_config: Path,
    source_dir: Path,
    imported_dir: Path,
    importer: Path,
    runner: Path,
) -> list[str]:
    return [
        sys.executable,
        str(import_script),
        "reference",
        "--config",
        str(reference_config),
        "--source-dir",
        str(source_dir),
        "--output-dir",
        str(imported_dir),
        "--tdr-importer",
        str(importer),
        "--runner",
        str(runner),
        "--skip-vela-run",
    ]


def vela_run_command(*, runner: Path, deck: Path) -> list[str]:
    return [str(runner), "--config", str(deck)]


def diagnostic_command(
    *,
    diagnostic_script: Path,
    out_dir: Path,
) -> list[str]:
    return [
        sys.executable,
        str(diagnostic_script),
        "--variants-root",
        str(out_dir / "variants"),
        "--sentaurus-root",
        str(out_dir / "sentaurus_exports"),
        "--out-dir",
        str(out_dir / "report"),
    ]


def run_recorded_command(
    argv: Sequence[str],
    cwd: Path,
    commands: list[dict[str, Any]],
    *,
    execute: bool = True,
    command_runner: Any = None,
    check: bool = True,
) -> int | None:
    normalized = [str(item) for item in argv]
    record = {
        "argv": normalized,
        "cwd": str(cwd),
        "returncode": None,
    }
    if not execute:
        commands.append(record)
        return None

    runner = subprocess.run if command_runner is None else command_runner
    completed = runner(normalized, cwd=cwd, check=False)
    record["returncode"] = int(completed.returncode)
    commands.append(record)
    if check and completed.returncode != 0:
        raise subprocess.CalledProcessError(completed.returncode, normalized)
    return int(completed.returncode)

def write_json(path: Path, data: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(data, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def git_state(repo: Path = REPO) -> tuple[str, bool]:
    git = shutil.which("git")
    if git is None:
        preferred = Path("D:/msys64/usr/bin/git.exe")
        if preferred.is_file():
            git = str(preferred)
    if git is None:
        return "unknown", True

    head = subprocess.run(
        [git, "-C", str(repo), "rev-parse", "HEAD"],
        capture_output=True,
        text=True,
        check=False,
    )
    status = subprocess.run(
        [git, "-C", str(repo), "status", "--porcelain"],
        capture_output=True,
        text=True,
        check=False,
    )
    if head.returncode != 0 or status.returncode != 0:
        return "unknown", True
    return head.stdout.strip(), bool(status.stdout.strip())


def write_variant_reference_configs(
    reference_config: Path,
    out_dir: Path,
) -> dict[str, dict[str, Any]]:
    base = json.loads(reference_config.read_text(encoding="utf-8-sig"))
    specs = variant_specs(out_dir)
    for spec in specs.values():
        config = json.loads(json.dumps(base))
        config.setdefault("tdr_doping", {})[
            "compensated_node_policy"
        ] = spec["compensated_doping_policy"]
        write_json(spec["reference_config"], config)
        spec["root"].mkdir(parents=True, exist_ok=True)
        spec["imported_dir"].mkdir(parents=True, exist_ok=True)
        spec["run_dir"].mkdir(parents=True, exist_ok=True)
    return specs


def apply_variant_doping_strategy(spec: Mapping[str, Any]) -> dict[str, Any]:
    """Materialize the historical p-side junction convention per variant."""
    imported_dir = Path(spec["imported_dir"])
    doping_path = imported_dir / "vela" / "doping.csv"
    metadata_path = imported_dir / "doping_metadata.json"
    if not doping_path.is_file() or not metadata_path.is_file():
        raise FileNotFoundError(
            "variant doping strategy requires Vela doping and importer metadata: "
            f"{doping_path}, {metadata_path}"
        )

    metadata = json.loads(metadata_path.read_text(encoding="utf-8-sig"))
    compensated = metadata.get("compensated_nodes", {})
    transformed_nodes: list[int] = []
    if spec["doping_strategy"] == "legacy":
        transformed_nodes = sorted(
            int(node["node_id"])
            for node in compensated.get("nodes", [])
            if (
                not bool(node.get("resolved"))
                and node.get("resolution_source")
                == LEGACY_DOPING_REPLAY["resolution_source"]
            )
        )
        transformed = set(transformed_nodes)
        with doping_path.open(newline="", encoding="utf-8-sig") as handle:
            reader = csv.DictReader(handle)
            fieldnames = list(reader.fieldnames or [])
            rows = [dict(row) for row in reader]
        required = {"node_id", "donors_cm3", "acceptors_cm3"}
        if not required.issubset(fieldnames):
            raise ValueError(f"invalid Vela doping CSV columns: {doping_path}")
        seen: set[int] = set()
        for row in rows:
            node_id = int(row["node_id"])
            if node_id not in transformed:
                continue
            donors = float(row["donors_cm3"])
            acceptors = float(row["acceptors_cm3"])
            scale = max(abs(donors), abs(acceptors), 1.0)
            if donors == 0.0 and acceptors > 0.0:
                seen.add(node_id)
                continue
            if (
                donors <= 0.0
                or acceptors <= 0.0
                or abs(donors - acceptors) > 1.0e-6 * scale
            ):
                raise ValueError(
                    f"legacy replay node {node_id} is not compensated in {doping_path}"
                )
            row["donors_cm3"] = "0"
            row["acceptors_cm3"] = f"{max(donors, acceptors):.17g}"
            seen.add(node_id)
        if seen != transformed:
            raise ValueError(
                f"legacy replay nodes missing from {doping_path}: "
                f"{sorted(transformed - seen)}"
            )
        with doping_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=fieldnames,
                lineterminator="\n",
            )
            writer.writeheader()
            writer.writerows(rows)

    record = {
        "schema": "vela.pn2d_compensated_doping_strategy.v1",
        "doping_strategy": spec["doping_strategy"],
        "import_policy": spec["compensated_doping_policy"],
        "legacy_replay": dict(LEGACY_DOPING_REPLAY),
        "transformed_node_ids": transformed_nodes,
    }
    write_json(imported_dir / "vela" / "doping_strategy.json", record)
    return record


def validate_doping_strategy_matrix(
    specs: Mapping[str, Mapping[str, Any]],
) -> dict[str, Any]:
    """Require one mesh and exactly two distinct doping strategy snapshots."""
    mesh_hashes: dict[str, str] = {}
    doping_hashes: dict[str, list[str]] = {}
    for name, spec in specs.items():
        vela_dir = Path(spec["imported_dir"]) / "vela"
        mesh_hashes[name] = sha256_file(vela_dir / "mesh.json")
        doping_hashes.setdefault(str(spec["doping_strategy"]), []).append(
            sha256_file(vela_dir / "doping.csv")
        )
    if len(set(mesh_hashes.values())) != 1:
        raise ValueError("2x3 doping comparison requires identical meshes")
    normalized = {
        strategy: sorted(set(hashes))
        for strategy, hashes in doping_hashes.items()
    }
    if any(len(hashes) != 1 for hashes in normalized.values()):
        raise ValueError("current variants within one doping strategy must match")
    if (
        len(normalized) != 2
        or len({hashes[0] for hashes in normalized.values()}) != 2
    ):
        raise ValueError("2x3 doping strategy matrix collapsed to identical inputs")
    return {
        "mesh_sha256": next(iter(mesh_hashes.values())),
        "doping_sha256_by_strategy": {
            strategy: hashes[0] for strategy, hashes in sorted(normalized.items())
        },
    }


def _output_files(out_dir: Path) -> list[Path]:
    if not out_dir.exists():
        return []
    return sorted(
        (
            path
            for path in out_dir.rglob("*")
            if path.is_file() and path.name != ARTIFACT_MANIFEST_NAME
        ),
        key=lambda path: path.relative_to(out_dir).as_posix(),
    )


def _tdr_records(source_dir: Path, biases: Sequence[float]) -> list[dict[str, Any]]:
    records = []
    for bias in biases:
        path = required_tdr_path(source_dir, bias)
        records.append({
            "bias_V": float(bias),
            "index": multibias_index_for_bias(bias),
            "path": path.resolve().as_posix(),
            "sha256": sha256_file(path),
            "size_bytes": path.stat().st_size,
        })
    return records


def input_file_signature(path: Path) -> dict[str, Any]:
    resolved = Path(path).resolve()
    if not resolved.is_file():
        return {
            "path": str(resolved), "exists": False,
            "size_bytes": None, "sha256": None,
        }
    return {
        "path": str(resolved), "exists": True,
        "size_bytes": resolved.stat().st_size, "sha256": sha256_file(resolved),
    }


def build_reuse_signature(
    *,
    biases: Sequence[float],
    reference_config: Path,
    source_dir: Path,
    specs: Mapping[str, Mapping[str, Any]],
    runner: Path,
    importer: Path,
    import_script: Path,
    diagnostic_script: Path,
) -> dict[str, Any]:
    return {
        "biases_V": [float(value) for value in biases],
        "implementation": "current_head",
        "reference_config": str(reference_config),
        "reference_config_sha256": sha256_file(reference_config),
        "source_dir": str(source_dir),
        "tools": {
            "runner": input_file_signature(runner),
            "importer": input_file_signature(importer),
            "import_script": input_file_signature(import_script),
            "diagnostic_script": input_file_signature(diagnostic_script),
        },
        "variants": {
            name: {
                "compensated_doping_policy": spec["compensated_doping_policy"],
                "doping_strategy": spec["doping_strategy"],
                "current_variant": spec["current_variant"],
                "current_approximation": spec["current_approximation"],
                "impact_ionization": dict(spec["impact_ionization"]),
                "doping_strategy_replay": (
                    dict(LEGACY_DOPING_REPLAY)
                    if spec["doping_strategy"] == "legacy" else None
                ),
            }
            for name, spec in specs.items()
        },
    }


def _manifest_record_is_current(out_dir: Path, record: Mapping[str, Any]) -> bool:
    raw_path = Path(str(record.get("path", "")))
    path = raw_path if raw_path.is_absolute() else out_dir / raw_path
    if not path.is_file():
        return False
    try:
        return (
            int(record.get("size_bytes", -1)) == path.stat().st_size
            and str(record.get("sha256", "")) == sha256_file(path)
        )
    except (OSError, TypeError, ValueError):
        return False


def reuse_manifest_matches(
    manifest_path: Path,
    *,
    out_dir: Path,
    git_head: str,
    dirty: bool,
    tdrs: Sequence[Mapping[str, Any]],
    signature: Mapping[str, Any],
) -> bool:
    """Accept reuse only for a clean, hash-identical recorded run."""
    if dirty or not manifest_path.is_file():
        return False
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8-sig"))
        if (
            manifest.get("schema") != MANIFEST_SCHEMA
            or manifest.get("schema_version") != 3
        ):
            return False
        if manifest.get("git_head") != git_head or bool(manifest.get("dirty")):
            return False
        if not isinstance(manifest.get("commands"), list) or not manifest["commands"]:
            return False
        parameters = manifest.get("parameters", {})
        if any(parameters.get(key) != value for key, value in signature.items()):
            return False
        recorded_tdrs = sorted(
            (dict(item) for item in manifest.get("tdrs", [])),
            key=lambda item: (float(item["bias_V"]), str(item["path"])),
        )
        expected_tdrs = sorted(
            (dict(item) for item in tdrs),
            key=lambda item: (float(item["bias_V"]), str(item["path"])),
        )
        if recorded_tdrs != expected_tdrs:
            return False
        records = list(manifest.get("artifacts", []))
        generated = manifest.get("generated", {})
        inputs = manifest.get("inputs", {})
        for key in ("decks", "model_configs"):
            records.extend(generated.get(key, []))
        for key in ("mesh", "doping", "materials"):
            records.extend(inputs.get(key, []))
        return bool(records) and all(
            _manifest_record_is_current(out_dir, record) for record in records
        )
    except (KeyError, OSError, TypeError, ValueError, json.JSONDecodeError):
        return False

def full_bias_points() -> list[float]:
    return [round(-0.05 * index, 12) for index in range(401)]


def _load_previous_full20_module() -> Any:
    module_path = REPO / "scripts" / "run_pn2d_coarse7x3_previous_full20_compare.py"
    spec = importlib.util.spec_from_file_location(
        "run_pn2d_coarse7x3_previous_full20_compare_reproducer",
        module_path,
    )
    if spec is None or spec.loader is None:
        raise ImportError(f"cannot load previous-full20 helper: {module_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_variant_deck(spec: Mapping[str, Any]) -> Path:
    base_config = Path(spec["imported_dir"]) / "vela" / "simulation_bv.json"
    if not base_config.is_file():
        raise FileNotFoundError(
            f"reference import did not generate Vela BV config: {base_config}"
        )
    module = _load_previous_full20_module()
    deck_path = module.write_previous_full20_config(
        base_config=base_config,
        out_dir=Path(spec["run_dir"]),
        output_csv_name=str(spec["output_csv_name"]),
        bias_points=full_bias_points(),
        config_name=str(spec["deck_name"]),
        vtk_subdir=str(spec["vtk_subdir"]),
        diagnostics_suffix=str(spec["diagnostics_suffix"]),
        current_approximation=str(spec["current_approximation"]),
    )
    deck = json.loads(deck_path.read_text(encoding="utf-8-sig"))
    deck.setdefault("solver", {})["impact_ionization"] = dict(
        spec["impact_ionization"]
    )
    diagnostics = deck.setdefault("sweep", {}).setdefault("diagnostics", {})
    diagnostics[str(spec["diagnostic_csv_kind"])] = {
        "enabled": True,
        "csv_file": str(Path(spec["run_dir"]) / spec["diagnostic_csv_name"]),
    }
    write_json(deck_path, deck)
    return deck_path


def _vela_outputs_complete(deck_path: Path) -> bool:
    if not deck_path.is_file():
        return False
    deck = json.loads(deck_path.read_text(encoding="utf-8-sig"))
    run_dir = deck_path.parent
    output_csv = Path(str(deck.get("output_csv", "")))
    if not output_csv.is_absolute():
        output_csv = run_dir / output_csv
    diagnostics = deck.get("sweep", {}).get("diagnostics", {})
    impact = deck.get("solver", {}).get("impact_ionization", {})
    diagnostic_kind = (
        "triangle_gss_sources"
        if impact.get("source_mapping_mode") == "triangle_gss_gradqf_truncated"
        else "sg_avalanche_edges"
    )
    source_raw = diagnostics.get(diagnostic_kind, {}).get("csv_file", "")
    source_csv = Path(str(source_raw))
    if not source_csv.is_absolute():
        source_csv = run_dir / source_csv
    vtk_prefix = Path(str(deck.get("sweep", {}).get("vtk_prefix", "")))
    if not vtk_prefix.is_absolute():
        vtk_prefix = run_dir / vtk_prefix
    vtk_matches = (
        list(vtk_prefix.parent.glob(vtk_prefix.name + "*.vtk"))
        if vtk_prefix.name
        else []
    )
    if not output_csv.is_file() or not source_csv.is_file() or not vtk_matches:
        return False
    with output_csv.open("r", newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    return any(
        abs(float(row.get("bias_V", "nan")) + 20.0) <= 1.0e-9
        and str(row.get("converged", "")).strip().lower()
            in {"1", "true", "yes"}
        for row in rows
    )


def vela_run_failed(
    returncode: int | None,
    deck_path: Path,
    *,
    prepare_only: bool,
) -> bool:
    if prepare_only:
        return False
    return returncode != 0 or not _vela_outputs_complete(deck_path)


def standard_report_paths(report_dir: Path) -> tuple[Path, Path, Path]:
    return (
        report_dir / "compensated_sg_replay.csv",
        report_dir / "compensated_sg_replay.json",
        report_dir / "compensated_sg_replay_report.md",
    )



def read_edge_mapping(report_dir: Path) -> list[dict[str, str]]:
    path = report_dir / "compensated_sg_replay.csv"
    if not path.is_file():
        return []
    keys = ("variant", "bias_V", "y_um", "side", "edge_id")
    rows: list[dict[str, str]] = []
    with path.open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            if all(row.get(key) not in (None, "") for key in keys):
                rows.append({key: str(row[key]) for key in keys})
    return sorted(
        rows,
        key=lambda row: (
            row["variant"],
            float(row["bias_V"]),
            float(row["y_um"]),
            row["side"],
            int(row["edge_id"]),
        ),
    )


def _variant_input_files(
    specs: Mapping[str, Mapping[str, Any]],
    deck_paths: Sequence[Path],
) -> tuple[list[Path], list[Path], list[Path]]:
    mesh_files: list[Path] = []
    doping_files: list[Path] = []
    material_files: list[Path] = []
    for spec, deck_path in zip(specs.values(), deck_paths):
        vela_dir = Path(spec["imported_dir"]) / "vela"
        mesh = vela_dir / "mesh.json"
        doping = vela_dir / "doping.csv"
        if mesh.is_file():
            mesh_files.append(mesh)
        if doping.is_file():
            doping_files.append(doping)
        if deck_path.is_file():
            deck = json.loads(deck_path.read_text(encoding="utf-8-sig"))
            material = Path(str(deck.get("materials_file", "")))
            if material and not material.is_absolute():
                material = deck_path.parent / material
            if material.is_file():
                material_files.append(material)
    return mesh_files, doping_files, material_files


def run_reproduction(
    args: argparse.Namespace,
    command_runner: Any = None,
) -> dict[str, Any]:
    source_dir = Path(args.source_dir).resolve()
    reference_config = Path(args.reference_config).resolve()
    out_dir = Path(args.out_dir).resolve()
    runner = Path(args.runner).resolve()
    importer = Path(args.importer).resolve()
    import_script = Path(args.import_script).resolve()
    diagnostic_script = Path(args.diagnostic_script).resolve()
    biases = [float(item) for item in args.biases]

    # Validate every audited source snapshot before creating or reusing outputs.
    tdrs = _tdr_records(source_dir, biases)
    head, dirty = git_state(REPO)

    sentaurus_root = out_dir / "sentaurus_exports"
    report_dir = out_dir / "report"
    sentaurus_root.mkdir(parents=True, exist_ok=True)
    report_dir.mkdir(parents=True, exist_ok=True)
    base_specs = variant_specs(out_dir)
    reuse_signature = build_reuse_signature(
        biases=biases,
        reference_config=reference_config,
        source_dir=source_dir,
        specs=base_specs,
        runner=runner,
        importer=importer,
        import_script=import_script,
        diagnostic_script=diagnostic_script,
    )
    reuse_allowed = bool(args.reuse_existing) and reuse_manifest_matches(
        out_dir / ARTIFACT_MANIFEST_NAME, out_dir=out_dir, git_head=head,
        dirty=dirty, tdrs=tdrs, signature=reuse_signature,
    )
    specs = write_variant_reference_configs(reference_config, out_dir)
    commands: list[dict[str, Any]] = []

    for tdr in tdrs:
        export_dir = sentaurus_root / f"sentaurus_{float(tdr['bias_V']):g}v"
        manifest_path = export_dir / "field_manifest.json"
        if args.skip_export:
            if not manifest_path.is_file():
                raise FileNotFoundError(
                    f"--skip-export requires existing field manifest: {manifest_path}"
                )
            validate_vector_field_manifest(manifest_path)
            continue
        if reuse_allowed and manifest_path.is_file():
            validate_vector_field_manifest(manifest_path)
            continue

        export_dir.mkdir(parents=True, exist_ok=True)
        run_recorded_command(
            export_command(
                importer=importer,
                tdr=Path(tdr["path"]),
                export_dir=export_dir,
            ),
            REPO,
            commands,
            execute=not args.prepare_only,
            command_runner=command_runner,
        )
        if not args.prepare_only:
            if not manifest_path.is_file():
                raise FileNotFoundError(
                    f"Sentaurus export did not write field manifest: {manifest_path}"
                )
            validate_vector_field_manifest(manifest_path)

    for spec in specs.values():
        base_config = spec["imported_dir"] / "vela" / "simulation_bv.json"
        if not (reuse_allowed and base_config.is_file()):
            run_recorded_command(
                reference_import_command(
                    import_script=import_script,
                    reference_config=spec["reference_config"],
                    source_dir=source_dir,
                    imported_dir=spec["imported_dir"],
                    importer=importer,
                    runner=runner,
                ),
                REPO,
                commands,
                command_runner=command_runner,
            )
        apply_variant_doping_strategy(spec)

    doping_matrix = validate_doping_strategy_matrix(specs)
    deck_paths: list[Path] = []
    failed_variants: list[str] = []
    for spec in specs.values():
        deck = write_variant_deck(spec)
        deck_paths.append(deck)

        if args.skip_vela_run:
            if not args.prepare_only and not _vela_outputs_complete(deck):
                raise FileNotFoundError(
                    "--skip-vela-run requires complete terminal, SG, and VTK outputs "
                    f"for {deck}"
                )
            continue
        if reuse_allowed and _vela_outputs_complete(deck):
            continue
        run_recorded_command(
            vela_run_command(runner=runner, deck=deck),
            spec["run_dir"],
            commands,
            execute=not args.prepare_only,
            command_runner=command_runner,
            check=False,
        )
        if vela_run_failed(
            commands[-1]["returncode"],
            deck,
            prepare_only=args.prepare_only,
        ):
            failed_variants.append(str(spec["name"]))

    standard_reports = standard_report_paths(report_dir)
    diagnostic_summary = standard_reports[1]
    outputs_complete = args.prepare_only or (
        not failed_variants
        and all(_vela_outputs_complete(deck) for deck in deck_paths)
    )
    if not (reuse_allowed and diagnostic_summary.is_file()):
        run_recorded_command(
            diagnostic_command(
                diagnostic_script=diagnostic_script,
                out_dir=out_dir,
            ),
            REPO,
            commands,
            execute=not args.prepare_only and outputs_complete,
            command_runner=command_runner,
        )
    if not args.prepare_only and outputs_complete:
        missing_reports = [path for path in standard_reports if not path.is_file()]
        if missing_reports:
            raise FileNotFoundError(
                "SG replay diagnostic missing standard artifacts: "
                + ", ".join(str(path) for path in missing_reports)
            )

    parameters = {
        **reuse_signature,
        "doping_matrix": doping_matrix,
        "failed_variants": failed_variants,
        "outputs_complete": bool(outputs_complete),
        "importer": str(importer),
        "prepare_only": bool(args.prepare_only),
        "reuse_existing": bool(args.reuse_existing),
        "reuse_accepted": reuse_allowed,
        "runner": str(runner),
        "skip_export": bool(args.skip_export),
        "skip_vela_run": bool(args.skip_vela_run),
    }
    model_configs = [
        spec["reference_config"] for spec in specs.values()
    ]
    existing_decks = [path for path in deck_paths if path.is_file()]
    mesh_files, doping_files, material_files = _variant_input_files(
        specs,
        deck_paths,
    )
    generation_commands: Sequence[Mapping[str, Any]] = commands
    if reuse_allowed:
        previous_manifest = json.loads(
            (out_dir / ARTIFACT_MANIFEST_NAME).read_text(encoding="utf-8-sig")
        )
        generation_commands = previous_manifest["commands"]
    manifest = build_artifact_manifest(
        out_dir=out_dir,
        git_head=head,
        dirty=dirty,
        parameters=parameters,
        commands=generation_commands,
        artifact_paths=_output_files(out_dir),
        tdrs=tdrs,
        generated_decks=existing_decks,
        model_configs=model_configs,
        mesh_files=mesh_files,
        doping_files=doping_files,
        material_files=material_files,
        edge_mapping=read_edge_mapping(report_dir),
    )
    manifest["invocation"] = {
        "reuse_requested": bool(args.reuse_existing),
        "reuse_accepted": reuse_allowed,
        "commands": [
            {
                "argv": [str(item) for item in command.get("argv", [])],
                "cwd": str(command.get("cwd", "")),
                "returncode": command.get("returncode"),
            }
            for command in commands
        ],
    }
    write_json(out_dir / ARTIFACT_MANIFEST_NAME, manifest)
    return manifest

def parse_biases(raw: str) -> list[float]:
    biases = [float(item.strip()) for item in raw.split(",") if item.strip()]
    if not biases:
        raise argparse.ArgumentTypeError("bias list must not be empty")
    try:
        for bias in biases:
            multibias_index_for_bias(bias)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(str(exc)) from exc
    return biases


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=DEFAULT_SOURCE_DIR)
    parser.add_argument(
        "--reference-config",
        type=Path,
        default=DEFAULT_REFERENCE_CONFIG,
    )
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--runner", type=Path, default=DEFAULT_RUNNER)
    parser.add_argument(
        "--importer",
        "--tdr-importer",
        dest="importer",
        type=Path,
        default=DEFAULT_IMPORTER,
    )
    parser.add_argument(
        "--import-script",
        type=Path,
        default=DEFAULT_IMPORT_SCRIPT,
    )
    parser.add_argument(
        "--diagnostic-script",
        type=Path,
        default=DEFAULT_DIAGNOSTIC_SCRIPT,
    )
    parser.add_argument(
        "--biases",
        type=parse_biases,
        default=list(DEFAULT_BIASES),
    )
    parser.add_argument("--prepare-only", action="store_true")
    parser.add_argument("--skip-export", action="store_true")
    parser.add_argument("--skip-vela-run", action="store_true")
    parser.add_argument("--reuse-existing", action="store_true")
    return parser


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    return build_parser().parse_args(argv)

def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    manifest = run_reproduction(args)
    print(json.dumps({
        "artifact_manifest": str(
            Path(args.out_dir).resolve() / ARTIFACT_MANIFEST_NAME
        ),
        "schema": manifest["schema"],
        "failed_variants": manifest.get("parameters", {}).get("failed_variants", []),
    }, indent=2))
    return 1 if manifest.get("parameters", {}).get("failed_variants") else 0


if __name__ == "__main__":
    raise SystemExit(main())