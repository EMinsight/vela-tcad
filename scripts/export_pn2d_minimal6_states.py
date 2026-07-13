#!/usr/bin/env python3
"""Export exact-bias Sentaurus states on both PN2D minimal6 topologies."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import sys
from datetime import datetime
from pathlib import Path
from typing import Callable, Sequence


REPO = Path(__file__).resolve().parents[1]
if __package__ in (None, ""):
    sys.path.insert(0, str(REPO))

from scripts.pn2d_minimal6_topology import load_topology  # noqa: E402
from scripts.run_pn2d_minimal6_sentaurus_gate import (  # noqa: E402
    DEFAULT_REMOTE_ROOT,
    MODELS_SOURCE,
    TOPOLOGY_FIXTURE,
    build_gate_bundle,
)
from scripts.run_sentaurus_vm_reference import (  # noqa: E402
    default_windows_openssh,
    run_checked,
    write_manifest,
)
from scripts.sentaurus_import import parse_quoted_list, parse_values_block  # noqa: E402


SCHEMA = "vela.pn2d_minimal6_states.v1"
REQUIRED_TOPOLOGIES = ("sketch", "mirror")
REQUIRED_BIASES = (0.0, -12.0, -19.0)
BIAS_TOLERANCE_V = 1.0e-12
COORDINATE_TOLERANCE_UM = 1.0e-12
SOURCE_DECK = (
    REPO / "reference_tcad" / "pn2d_sentaurus2018_minimal6" / "source"
    / "pn2d_minimal6_state_sdevice.cmd"
)
DEFAULT_OUTPUT_DIR = (
    REPO / "build-release" / "reference_tcad" / "pn2d_sentaurus2018_minimal6"
    / "state_exports"
)
DEFAULT_IMPORTER = REPO / "build-release" / (
    "sentaurus_import.exe" if os.name == "nt" else "sentaurus_import"
)
_SAFE_RUN_ID = re.compile(r"^[A-Za-z0-9_.-]+$")
_SAFE_REMOTE_COMPONENT = re.compile(r"^[A-Za-z0-9_./~:-]+$")

_FIELD_CONTRACT = {
    "ElectrostaticPotential": (1, "V"),
    "eQuasiFermiPotential": (1, "V"),
    "hQuasiFermiPotential": (1, "V"),
    "eDensity": (1, "cm^-3"),
    "hDensity": (1, "cm^-3"),
    "ElectricField": (2, "V*cm^-1"),
    "eCurrentDensity": (2, "A*cm^-2"),
    "hCurrentDensity": (2, "A*cm^-2"),
    "eMobility": (1, "cm^2*V^-1*s^-1"),
    "hMobility": (1, "cm^2*V^-1*s^-1"),
    "eAlphaAvalanche": (1, "cm^-1"),
    "hAlphaAvalanche": (1, "cm^-1"),
}


def _bias_tag(bias_V: float) -> str:
    sign = "m" if bias_V < 0.0 else "p"
    magnitude = format(abs(bias_V), ".17g").replace(".", "p")
    return f"{sign}{magnitude}V"


def validate_final_bias(requested_bias_V: float, actual_bias_V: float) -> float:
    requested = float(requested_bias_V)
    actual = float(actual_bias_V)
    if not math.isfinite(actual):
        raise ValueError(f"final Anode contact voltage is not finite: {actual_bias_V}")
    if abs(actual - requested) > BIAS_TOLERANCE_V:
        raise ValueError(
            f"final Anode contact voltage {actual:.17g} V does not match requested "
            f"{requested:.17g} V within 1e-12 V"
        )
    return actual


def validate_state_matrix(states: Sequence[dict[str, object]]) -> list[tuple[str, float]]:
    required = {(topology, bias) for topology in REQUIRED_TOPOLOGIES for bias in REQUIRED_BIASES}
    matrix: list[tuple[str, float]] = []
    for state in states:
        topology = str(state.get("topology_id", ""))
        requested = float(state.get("requested_bias_V", math.nan))
        if state.get("status") != "passed":
            raise ValueError(f"state {topology} at {requested:g} V is not passed")
        try:
            actual = validate_final_bias(requested, float(state.get("actual_bias_V", math.nan)))
        except ValueError as error:
            raise ValueError(f"{topology} at {requested:g} V: {error}") from error
        key = (topology, requested)
        if key in matrix:
            raise ValueError(f"duplicate state {topology} at {requested:g} V")
        matrix.append(key)
        if actual != actual:  # Defensive; validate_final_bias already rejects NaN.
            raise ValueError("unreachable non-finite bias")
    if set(matrix) != required or len(matrix) != len(required):
        missing = sorted(required - set(matrix))
        extra = sorted(set(matrix) - required)
        raise ValueError(f"exact six-state matrix mismatch; missing={missing}, extra={extra}")
    return matrix


def validate_field_manifest(manifest: dict[str, object]) -> dict[str, dict[str, object]]:
    fields = manifest.get("fields")
    if not isinstance(fields, list):
        raise ValueError("field_manifest.json must contain a fields list")
    selected: dict[str, dict[str, object]] = {}
    for name, (components, unit) in _FIELD_CONTRACT.items():
        candidates = [field for field in fields if isinstance(field, dict) and field.get("name") == name]
        expected = {
            "region": 0,
            "components": components,
            "unit": unit,
            "mapping_status": "complete",
            "global_node_mapping": "global_vertex_order",
        }
        match = next(
            (field for field in candidates if all(field.get(key) == value for key, value in expected.items())),
            None,
        )
        if match is None:
            contract = ", ".join(f"{key}={value}" for key, value in expected.items())
            raise ValueError(f"{name} must satisfy {contract}; got {candidates}")
        selected[name] = match
    return selected


def canonical_minimal6_coordinates() -> dict[int, tuple[float, float]]:
    topology = load_topology(TOPOLOGY_FIXTURE, "sketch")
    return {
        label - 1: tuple(topology.nodes[label])
        for label in sorted(topology.nodes)
    }


def _read_csv(path: Path, required: set[str]) -> list[dict[str, str]]:
    if not path.is_file():
        raise ValueError(f"missing neutral export file: {path.name}")
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        columns = set(reader.fieldnames or [])
        if not required.issubset(columns):
            raise ValueError(f"{path.name} missing columns: {sorted(required - columns)}")
        return list(reader)


def _canonical_source_ids(
    export_dir: Path, coordinates: dict[int, tuple[float, float]]
) -> dict[int, int]:
    rows = _read_csv(export_dir / "nodes.csv", {"id", "x_um", "y_um"})
    if len(rows) != len(coordinates):
        raise ValueError(f"expected {len(coordinates)} canonical nodes, found {len(rows)}")
    source_for_canonical: dict[int, int] = {}
    used_sources: set[int] = set()
    for canonical, expected in coordinates.items():
        matches = [
            int(row["id"])
            for row in rows
            if abs(float(row["x_um"]) - expected[0]) < COORDINATE_TOLERANCE_UM
            and abs(float(row["y_um"]) - expected[1]) < COORDINATE_TOLERANCE_UM
        ]
        if len(matches) != 1:
            raise ValueError(
                f"canonical node {canonical} requires one exact coordinate mapping; got {matches}"
            )
        if matches[0] in used_sources:
            raise ValueError("canonical node mapping reuses a source node")
        source_for_canonical[canonical] = matches[0]
        used_sources.add(matches[0])
    return source_for_canonical


def _read_scalar_field(export_dir: Path, name: str) -> dict[int, float]:
    path = export_dir / "fields" / f"{name}_region0.csv"
    if not path.is_file():
        role = {"eQuasiFermiPotential": "phin", "hQuasiFermiPotential": "phip"}.get(name)
        suffix = f" ({role})" if role else ""
        raise ValueError(f"missing required field {name}{suffix}; density-derived QF is forbidden")
    rows = _read_csv(path, {"node_id", "component0"})
    values: dict[int, float] = {}
    for row in rows:
        node_id = int(row["node_id"])
        value = float(row["component0"])
        if node_id in values or not math.isfinite(value):
            raise ValueError(f"{name} has duplicate or non-finite node {node_id}")
        values[node_id] = value
    return values


def write_state_csv(
    export_dir: Path, coordinates: dict[int, tuple[float, float]] | None = None
) -> Path:
    export_dir = Path(export_dir)
    manifest_path = export_dir / "field_manifest.json"
    if not manifest_path.is_file():
        raise ValueError("missing neutral export file: field_manifest.json")
    validate_field_manifest(json.loads(manifest_path.read_text(encoding="utf-8")))
    raw = {
        "psi_V": _read_scalar_field(export_dir, "ElectrostaticPotential"),
        "phin_V": _read_scalar_field(export_dir, "eQuasiFermiPotential"),
        "phip_V": _read_scalar_field(export_dir, "hQuasiFermiPotential"),
        "n_m3": _read_scalar_field(export_dir, "eDensity"),
        "p_m3": _read_scalar_field(export_dir, "hDensity"),
    }
    source_ids = _canonical_source_ids(export_dir, coordinates or canonical_minimal6_coordinates())
    output = export_dir / "state.csv"
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["node_id", "psi_V", "phin_V", "phip_V", "n_m3", "p_m3"])
        for canonical in sorted(source_ids):
            source = source_ids[canonical]
            try:
                row = [raw[name][source] for name in ("psi_V", "phin_V", "phip_V", "n_m3", "p_m3")]
            except KeyError as error:
                raise ValueError(f"state field mapping is missing source node {source}") from error
            row[3] *= 1.0e6
            row[4] *= 1.0e6
            writer.writerow([canonical, *(format(value, ".17g") for value in row)])
    return output


def _render_deck(target_bias_V: float, tag: str, destination: Path) -> None:
    source = SOURCE_DECK.read_text(encoding="utf-8")
    rendered = source.replace("__TARGET_BIAS_V__", format(target_bias_V, ".17g"))
    rendered = rendered.replace("__BIAS_TAG__", tag)
    if "__TARGET_BIAS_V__" in rendered or "__BIAS_TAG__" in rendered:
        raise ValueError("state deck placeholder replacement is incomplete")
    destination.write_text(rendered, encoding="utf-8", newline="\n")


def _validated_requested_matrix(
    topology_ids: Sequence[str], biases: Sequence[float]
) -> tuple[tuple[str, ...], tuple[float, ...]]:
    topologies = tuple(topology_ids)
    bias_values = tuple(float(value) for value in biases)
    if topologies != REQUIRED_TOPOLOGIES:
        raise ValueError(f"topologies must be exactly {REQUIRED_TOPOLOGIES}")
    if bias_values != REQUIRED_BIASES:
        raise ValueError(f"biases must be exactly {REQUIRED_BIASES}; nearest substitutes are forbidden")
    return topologies, bias_values


def prepare_exports(
    *, topology_ids: Sequence[str], biases: Sequence[float], run_id: str,
    output_dir: Path, ssh_target: str, remote_root: str = DEFAULT_REMOTE_ROOT,
    importer: Path = DEFAULT_IMPORTER,
) -> dict[str, object]:
    topologies, bias_values = _validated_requested_matrix(topology_ids, biases)
    if not _SAFE_RUN_ID.fullmatch(run_id):
        raise ValueError("run ID contains unsupported characters")
    if not _SAFE_REMOTE_COMPONENT.fullmatch(remote_root):
        raise ValueError("remote root contains unsupported characters")
    if not SOURCE_DECK.is_file() or not MODELS_SOURCE.is_file():
        raise FileNotFoundError("missing minimal6 state deck or model parameters")
    run_root = (Path(output_dir) / run_id).resolve()
    states: list[dict[str, object]] = []
    for topology_id in topologies:
        for bias in bias_values:
            tag = _bias_tag(bias)
            state_root = run_root / "states" / topology_id / tag
            bundle = state_root / "source"
            artifacts = state_root / "artifacts"
            neutral = state_root / "export"
            gate_bundle = build_gate_bundle(topology_id, bundle)
            gate_deck = bundle / "pn2d_minimal6_gate_sdevice.cmd"
            gate_deck.unlink()
            deck_name = f"pn2d_minimal6_state_{tag}_sdevice.cmd"
            _render_deck(bias, tag, bundle / deck_name)
            staged = ["pn2d_minimal6.grd", "pn2d_minimal6.dat", "models.par", deck_name]
            remote_dir = f"{remote_root.rstrip('/')}/{run_id}/{topology_id}/{tag}"
            plot_stem = f"pn2d_minimal6_state_{tag}"
            final_tdr_name = f"{plot_stem}.tdr"
            current_plt_name = f"{plot_stem}.plt"
            log_name = f"{plot_stem}_des.log"
            stdout_name = f"run_{plot_stem}.out"
            returned_files = [
                final_tdr_name,
                current_plt_name,
                log_name,
                "pn2d_minimal6.tdr",
                "pn2d_minimal6.grd",
                "pn2d_minimal6.dat",
                "run_tdx_dfise_to_tdr.out",
                stdout_name,
            ]
            states.append({
                "topology_id": topology_id,
                "requested_bias_V": bias,
                "bias_tag": tag,
                "bundle_dir": str(bundle),
                "artifacts_dir": str(artifacts),
                "export_dir": str(neutral),
                "remote_dir": remote_dir,
                "deck_name": deck_name,
                "staged_files": staged,
                "topology_contract": gate_bundle["topology_contract"],
                "remote_commands": [
                    f"cd {remote_dir} && tdx -d pn2d_minimal6.grd pn2d_minimal6.dat "
                    "pn2d_minimal6.tdr > run_tdx_dfise_to_tdr.out 2>&1",
                    f"cd {remote_dir} && sdevice {deck_name} > {stdout_name} 2>&1",
                ],
                "final_tdr_name": final_tdr_name,
                "current_plt_name": current_plt_name,
                "log_name": log_name,
                "stdout_name": stdout_name,
                "returned_files": returned_files,
                "status": "prepared",
            })
    manifest: dict[str, object] = {
        "schema": SCHEMA,
        "run_id": run_id,
        "ssh_target": ssh_target,
        "remote_root": remote_root,
        "importer": str(Path(importer).resolve()),
        "bias_tolerance_V": BIAS_TOLERANCE_V,
        "outputs_complete": False,
        "states": states,
        "manifest_path": str(run_root / "manifest.json"),
    }
    write_manifest(Path(str(manifest["manifest_path"])), manifest)
    return manifest


def _parse_final_anode_bias(plt_path: Path) -> float:
    text = plt_path.read_text(errors="replace")
    datasets = parse_quoted_list(text, "datasets")
    if "Anode OuterVoltage" not in datasets:
        raise ValueError(f"{plt_path.name} lacks Anode OuterVoltage")
    rows = parse_values_block(text, len(datasets))
    if not rows:
        raise ValueError(f"{plt_path.name} has no bias rows")
    return float(rows[-1][datasets.index("Anode OuterVoltage")])


def _live_executor(
    state: dict[str, object], *, ssh_bin: str, scp_bin: str,
    ssh_target: str, importer: Path,
) -> dict[str, object]:
    if not importer.is_file():
        raise FileNotFoundError(f"Sentaurus importer is not built: {importer}")
    bundle = Path(str(state["bundle_dir"]))
    artifacts = Path(str(state["artifacts_dir"]))
    neutral = Path(str(state["export_dir"]))
    artifacts.mkdir(parents=True, exist_ok=True)
    remote_dir = str(state["remote_dir"])
    returned = [str(name) for name in state["returned_files"]]

    def return_argv(name: str) -> list[str]:
        return [
            scp_bin,
            f"{ssh_target}:{remote_dir}/{name}",
            str(artifacts) + os.sep,
        ]

    try:
        run_checked([ssh_bin, ssh_target, f"mkdir -p {remote_dir}"])
        for name in state["staged_files"]:
            run_checked([scp_bin, str(bundle / str(name)), f"{ssh_target}:{remote_dir}/"])
        for command in state["remote_commands"]:
            run_checked([ssh_bin, ssh_target, str(command)])
    except Exception:
        recovery_errors: list[str] = []
        for name in returned:
            try:
                run_checked(return_argv(name))
            except Exception as recovery_error:
                recovery_errors.append(f"{name}: {recovery_error}")
        state["artifact_recovery_errors"] = recovery_errors
        raise
    for name in returned:
        run_checked(return_argv(name))
    final_tdr = artifacts / str(state["final_tdr_name"])
    # The current C++ CLI performs the export-neutral operation when --export-dir is present.
    run_checked([
        str(importer), "--tdr", str(final_tdr), "--export-dir", str(neutral),
        "--compensated-doping-policy", "reported",
    ])
    return {
        "actual_bias_V": _parse_final_anode_bias(artifacts / str(state["current_plt_name"])),
        "export_dir": str(neutral),
    }


def run_exports(
    manifest: dict[str, object],
    *, executor: Callable[[dict[str, object]], dict[str, object] | None],
) -> None:
    manifest_path = Path(str(manifest["manifest_path"]))
    manifest["outputs_complete"] = False
    write_manifest(manifest_path, manifest)
    for state in manifest["states"]:
        try:
            result = executor(state) or {}
            if "actual_bias_V" not in result:
                raise ValueError("executor result is missing actual_bias_V")
            try:
                actual_bias = float(result["actual_bias_V"])
            except (TypeError, ValueError) as error:
                raise ValueError(
                    f"executor result has invalid actual_bias_V: {result['actual_bias_V']!r}"
                ) from error
            state["actual_bias_V"] = validate_final_bias(
                float(state["requested_bias_V"]), actual_bias
            )
            export_dir = Path(str(result.get("export_dir", state["export_dir"])))
            if not export_dir.is_dir():
                raise ValueError(f"missing neutral export directory: {export_dir}")
            state["state_csv"] = str(write_state_csv(export_dir))
            state["field_manifest"] = str(export_dir / "field_manifest.json")
            state["status"] = "passed"
            write_manifest(manifest_path, manifest)
        except Exception as error:
            state["status"] = "failed"
            state["error"] = str(error)
            manifest["error"] = str(error)
            write_manifest(manifest_path, manifest)
            raise
    validate_state_matrix(manifest["states"])
    manifest["outputs_complete"] = True
    manifest.pop("error", None)
    write_manifest(manifest_path, manifest)


def _parse_csv_values(raw: str) -> tuple[float, ...]:
    values = tuple(float(value.strip()) for value in raw.split(",") if value.strip())
    if not values:
        raise ValueError("bias list is empty")
    return values


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--topologies", default="sketch,mirror")
    parser.add_argument("--biases", default="0,-12,-19")
    parser.add_argument("--run-id", default=None)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--ssh-target", default="sentaurus")
    parser.add_argument("--ssh-bin", default=None)
    parser.add_argument("--scp-bin", default=None)
    parser.add_argument("--remote-root", default=DEFAULT_REMOTE_ROOT)
    parser.add_argument("--importer", type=Path, default=DEFAULT_IMPORTER)
    parser.add_argument("--dry-run", action="store_true")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    manifest: dict[str, object] | None = None
    try:
        run_id = args.run_id or datetime.now().strftime("minimal6_states_%Y%m%d_%H%M%S")
        manifest = prepare_exports(
            topology_ids=tuple(value.strip() for value in args.topologies.split(",") if value.strip()),
            biases=_parse_csv_values(args.biases),
            run_id=run_id,
            output_dir=args.output_dir,
            ssh_target=args.ssh_target,
            remote_root=args.remote_root,
            importer=args.importer,
        )
        if not args.dry_run:
            ssh_bin = args.ssh_bin or default_windows_openssh("ssh")
            scp_bin = args.scp_bin or default_windows_openssh("scp")
            run_exports(
                manifest,
                executor=lambda state: _live_executor(
                    state, ssh_bin=ssh_bin, scp_bin=scp_bin,
                    ssh_target=args.ssh_target, importer=args.importer.resolve(),
                ),
            )
        print(json.dumps(manifest, indent=2))
        return 0
    except Exception as error:  # noqa: BLE001 - partial manifest is the contract.
        if manifest is not None:
            manifest["outputs_complete"] = False
            manifest["error"] = str(error)
            write_manifest(Path(str(manifest["manifest_path"])), manifest)
        print(str(error), file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
