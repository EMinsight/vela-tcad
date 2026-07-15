#!/usr/bin/env python3
"""Create immutable, segmented Minimal6 Vela/Sentaurus diagnostic-sweep packages."""
from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import math
import shutil
import subprocess
from pathlib import Path
from typing import Any

REPO = Path(__file__).resolve().parents[1]
DISCLAIMER = "minimal6 diagnostic sweep; not a physical BV curve"
SCHEMA = "vela.pn2d_minimal6_sweep_manifest.v1"


def integer_targets() -> tuple[float, ...]:
    return tuple(float(-value) for value in range(21))


def classify_branch(sentaurus_current_A_per_um: float | None, vela_current_A_per_um: float | None) -> str:
    if sentaurus_current_A_per_um is None or vela_current_A_per_um is None or sentaurus_current_A_per_um == 0.0:
        return "unidentified"
    ratio = abs(vela_current_A_per_um / sentaurus_current_A_per_um)
    if 0.1 <= ratio <= 10.0:
        return "multiplication_like"
    if ratio <= 1.0e-3:
        return "leakage_like"
    return "unidentified"


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _strip_allowed(deck: dict[str, Any]) -> dict[str, Any]:
    value = copy.deepcopy(deck)
    for key in ("mesh_file", "node_doping_file"):
        value.pop(key, None)
    sweep = value.get("sweep", {})
    for key in ("start", "stop", "initial_state_file", "write_state_every_point_prefix", "csv_file"):
        sweep.pop(key, None)
    for cfg in sweep.get("diagnostics", {}).values():
        if isinstance(cfg, dict):
            cfg.pop("csv_file", None)
            cfg.pop("summary_file", None)
    return value


def validate_segment_deck(template: dict[str, Any], generated: dict[str, Any]) -> None:
    if _strip_allowed(template) != _strip_allowed(generated):
        raise ValueError("generated segment changed immutable physics or solver configuration")
    for key in ("mesh_file", "node_doping_file"):
        if not isinstance(generated.get(key), str) or not generated[key]:
            raise ValueError(f"segment deck lacks topology input {key}")
    sweep = generated.get("sweep")
    if not isinstance(sweep, dict) or not all(math.isfinite(float(sweep[key])) for key in ("start", "stop")):
        raise ValueError("segment deck lacks finite start/stop")
    if sweep["stop"] != sweep["start"] - 1.0:
        raise ValueError("diagnostic segments must end at the next exact integer volt")
    if not isinstance(sweep.get("write_state_every_point_prefix"), str) or not sweep["write_state_every_point_prefix"]:
        raise ValueError("segment deck lacks unique state prefix")


def bias_token(bias: float) -> str:
    return ("m" if bias < 0.0 else "") + f"{abs(bias):.6f}".replace(".", "p")


def segment_state_path(root: Path, topology: str, segment_start: float, target_bias: float) -> Path:
    return root / "vela" / topology / "states" / f"segment_{abs(int(segment_start)):02d}_bias_{bias_token(target_bias)}.csv"

def make_segment_deck(template: dict[str, Any], *, topology: str, segment_start: float, root: Path, restart: Path | None) -> dict[str, Any]:
    deck = copy.deepcopy(template)
    deck["mesh_file"] = str(root / "inputs" / topology / "mesh.json")
    deck["node_doping_file"] = str(root / "inputs" / topology / "doping.csv")
    sweep = deck["sweep"]
    sweep.update({"start": segment_start, "stop": segment_start - 1.0,
                  "csv_file": str(root / "vela" / topology / f"segment_{abs(int(segment_start)):02d}.csv"),
                  "write_state_every_point_prefix": str(root / "vela" / topology / "states" / f"segment_{abs(int(segment_start)):02d}")})
    diagnostics = sweep.get("diagnostics", {})
    for name, cfg in diagnostics.items():
        if isinstance(cfg, dict) and "csv_file" in cfg:
            cfg["csv_file"] = str(root / "vela" / topology / "diagnostics" / f"segment_{abs(int(segment_start)):02d}_{name}.csv")
    if restart is None:
        sweep.pop("initial_state_file", None)
    else:
        sweep["initial_state_file"] = str(restart)
    validate_segment_deck(template, deck)
    return deck


def read_vela_endpoint(curve_csv: Path, terminal_csv: Path, target_bias_V: float) -> dict[str, float]:
    """Read exact endpoint observables from Vela's primary and terminal diagnostics."""
    with curve_csv.open(newline="", encoding="utf-8") as handle:
        curve_rows = [row for row in csv.DictReader(handle) if abs(float(row["bias_V"]) - target_bias_V) <= 1.0e-12]
    if len(curve_rows) != 1:
        raise ValueError("Vela curve lacks exactly one target-bias endpoint")
    with terminal_csv.open(newline="", encoding="utf-8") as handle:
        terminal_rows = [row for row in csv.DictReader(handle) if abs(float(row["bias_V"]) - target_bias_V) <= 1.0e-12]
    contacts = {row["contact"]: row for row in terminal_rows}
    if set(contacts) != {"Anode", "Cathode"}:
        raise ValueError("Vela terminal diagnostics lack both contacts at the endpoint")
    values = {
        "anode_current_A_per_um": float(contacts["Anode"]["I_sgflux_A_per_um"]),
        "cathode_current_A_per_um": float(contacts["Cathode"]["I_sgflux_A_per_um"]),
        "max_field_V_per_m": float(curve_rows[0]["max_electric_field_V_per_m"]),
        "reconstructed_source_integral_s_inv_per_cm": float(contacts["Anode"]["sg_avalanche_source_integral_total"]),
    }
    if not all(math.isfinite(value) for value in values.values()):
        raise ValueError("Vela endpoint contains non-finite observables")
    return values

def copy_topology_inputs(authoritative_state_root: Path, destination_root: Path, topologies: tuple[str, ...] = ("sketch", "mirror")) -> dict[str, dict[str, str]]:
    """Copy only canonical 0 V mesh/doping inputs into the independent sweep root."""
    hashes: dict[str, dict[str, str]] = {}
    for topology in topologies:
        source = authoritative_state_root / "states" / topology / "p0V" / "export"
        target = destination_root / "inputs" / topology
        target.mkdir(parents=True, exist_ok=True)
        hashes[topology] = {}
        for name in ("mesh.json", "doping.csv"):
            member = source / name
            if not member.is_file():
                raise FileNotFoundError(f"authoritative topology input is missing: {member}")
            copied = target / name
            shutil.copyfile(member, copied)
            hashes[topology][name] = _sha(copied)
    return hashes

def record_transition(
    manifest: dict[str, Any], *, solver: str, topology: str, start_bias_V: float,
    target_bias_V: float, exit_code: int, actual_bias_V: float | None,
    state_path: Path | None, observables: dict[str, float] | None, diagnostics: dict[str, str] | None = None,
    incomplete_reason: str | None = None,
) -> dict[str, Any]:
    """Append immutable evidence for one segment; rejected rows carry no observables."""
    exact = actual_bias_V is not None and abs(actual_bias_V - target_bias_V) <= 1.0e-12
    accepted = exit_code == 0 and exact and state_path is not None and state_path.is_file() and observables is not None
    row: dict[str, Any] = {
        "solver": solver, "topology": topology, "start_bias_V": start_bias_V,
        "target_bias_V": target_bias_V, "actual_bias_V": actual_bias_V,
        "exit_code": exit_code, "status": "accepted" if accepted else "rejected",
        "state_path": None if state_path is None else str(state_path),
        "state_sha256": _sha(state_path) if accepted else None,
        "observables": dict(observables) if accepted else None,
        "stdout": "" if diagnostics is None else diagnostics.get("stdout", ""),
        "stderr": "" if diagnostics is None else diagnostics.get("stderr", ""),
    }
    if accepted:
        required = {"anode_current_A_per_um", "cathode_current_A_per_um", "max_field_V_per_m", "native_source_integral_s_inv_per_cm", "reconstructed_source_integral_s_inv_per_cm"}
        if observables is None or set(observables) != required or not all(math.isfinite(float(value)) for value in observables.values()):
            raise ValueError("accepted checkpoint lacks complete finite observables")
        manifest.setdefault("accepted_checkpoints", []).append(row)
    else:
        if incomplete_reason:
            row["incomplete_reason"] = incomplete_reason
        if manifest.get("failed_transition") is None:
            manifest["failed_transition"] = row
    return row


def run_vela_subprocess_segment(root: Path, executable: Path, segment: dict[str, Any]) -> dict[str, Any]:
    """Run one generated Vela deck and return evidence without inventing a native source."""
    deck = root / str(segment["deck"])
    completed = subprocess.run([str(executable), "--config", str(deck)], cwd=root, text=True, capture_output=True)
    target = float(segment["target_bias_V"])
    topology = str(segment["topology"])
    start = float(segment["start_bias_V"])
    state = segment_state_path(root, topology, start, target)
    diagnostics = root / "vela" / topology / "diagnostics"
    terminal = diagnostics / f"segment_{abs(int(start)):02d}_terminal_current_method_compare.csv"
    curve = root / "vela" / topology / f"segment_{abs(int(start)):02d}.csv"
    if completed.returncode != 0 or not state.is_file() or not terminal.is_file() or not curve.is_file():
        return {"exit_code": completed.returncode, "actual_bias_V": None, "state_path": None, "observables": None,
                "stdout": completed.stdout, "stderr": completed.stderr}
    endpoint = read_vela_endpoint(curve, terminal, target)
    # Vela has a reconstructed source only.  Native Sentaurus source is deliberately
    # absent here and causes the common package gate to reject rather than relabel it.
    return {"exit_code": completed.returncode, "actual_bias_V": target, "state_path": state, "observables": None,
            "partial_observables": endpoint, "stdout": completed.stdout, "stderr": completed.stderr,
            "incomplete_reason": "Vela has no native Sentaurus source integral"}

def execute_segments(manifest: dict[str, Any], root: Path, runner) -> None:
    """Run pending segments in order and permanently retain the first failed one."""
    for segment in manifest.get("segments", []):
        if segment.get("status") != "pending":
            continue
        result = runner(segment)
        row = record_transition(
            manifest, solver=str(segment["solver"]), topology=str(segment["topology"]),
            start_bias_V=float(segment["start_bias_V"]), target_bias_V=float(segment["target_bias_V"]),
            exit_code=int(result.get("exit_code", 1)), actual_bias_V=result.get("actual_bias_V"),
            state_path=result.get("state_path"), observables=result.get("observables"),
            diagnostics={"stdout": str(result.get("stdout", "")), "stderr": str(result.get("stderr", ""))},
            incomplete_reason=None if result.get("incomplete_reason") is None else str(result["incomplete_reason"]),
        )
        segment["status"] = row["status"]
        if row["status"] != "accepted":
            break

def validate_sweep_manifest(manifest: dict[str, Any]) -> None:
    if manifest.get("schema") != SCHEMA:
        raise ValueError("invalid sweep manifest schema")
    if tuple(manifest.get("targets_V", ())) != integer_targets() and tuple(manifest.get("targets_V", ())) != (0.0, -1.0):
        raise ValueError("sweep manifest has non-canonical targets")
    if manifest.get("interpolation", "forbidden") != "forbidden":
        raise ValueError("sweep manifest must forbid interpolation")
    for row in manifest.get("accepted_checkpoints", []):
        if row.get("status") != "accepted" or abs(float(row["actual_bias_V"]) - float(row["target_bias_V"])) > 1.0e-12:
            raise ValueError("accepted checkpoint is not exact")
        state = Path(str(row["state_path"]))
        if not state.is_file() or _sha(state) != row.get("state_sha256"):
            raise ValueError("accepted checkpoint state is missing or hash-tampered")
        if not isinstance(row.get("observables"), dict):
            raise ValueError("accepted checkpoint lacks observables")
    failed = manifest.get("failed_transition")
    if failed is not None and (failed.get("status") != "rejected" or failed.get("observables") is not None):
        raise ValueError("failed transition must preserve no fabricated observables")

def _write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, sort_keys=True, indent=2) + "\n", encoding="utf-8")


def write_sentaurus_decks(root: Path, deck_template: Path, topologies: tuple[str, ...] = ("sketch", "mirror")) -> list[dict[str, Any]]:
    """Create independent exact-bias Sentaurus decks; no authoritative manifest is touched."""
    source = deck_template.read_text(encoding="utf-8")
    rows: list[dict[str, Any]] = []
    for topology in topologies:
        for bias in integer_targets():
            tag = f"{topology}_{bias_token(bias)}"
            path = root / "sentaurus" / topology / "decks" / f"{tag}.cmd"
            payload = source.replace("__BIAS_TAG__", tag).replace("__TARGET_BIAS_V__", f"{bias:.1f}")
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(payload, encoding="utf-8")
            rows.append({"solver": "sentaurus", "topology": topology, "target_bias_V": bias,
                         "deck": str(path.relative_to(root)), "deck_sha256": _sha(path), "status": "pending",
                         "checkpoint_tdr": str((root / "sentaurus" / topology / "checkpoints" / f"{tag}.tdr").relative_to(root))})
    return rows

def initialise_package(root: Path, template_path: Path, topologies: tuple[str, ...] = ("sketch", "mirror"), authoritative_state_root: Path | None = None) -> Path:
    template = json.loads(template_path.read_text(encoding="utf-8"))
    input_hashes = {} if authoritative_state_root is None else copy_topology_inputs(authoritative_state_root, root, topologies)
    targets = integer_targets()
    segments: list[dict[str, Any]] = []
    for topology in topologies:
        restart: Path | None = None
        for start in targets[:-1]:
            deck = make_segment_deck(template, topology=topology, segment_start=start, root=root, restart=restart)
            deck_path = root / "vela" / topology / "decks" / f"segment_{abs(int(start)):02d}.json"
            _write_json(deck_path, deck)
            segments.append({"solver": "vela", "topology": topology, "start_bias_V": start, "target_bias_V": start - 1.0,
                             "deck": str(deck_path.relative_to(root)), "deck_sha256": _sha(deck_path), "status": "pending"})
            restart = segment_state_path(root, topology, start, start - 1.0)
    sentaurus_template = REPO / "reference_tcad" / "pn2d_sentaurus2018_minimal6" / "source" / "pn2d_minimal6_sweep_sdevice.cmd"
    sentaurus_segments = write_sentaurus_decks(root, sentaurus_template, topologies)
    manifest = {"schema": SCHEMA, "diagnostic_disclaimer": DISCLAIMER, "targets_V": list(targets),
                "template": {"path": str(template_path), "sha256": _sha(template_path)}, "topology_input_sha256": input_hashes, "segments": segments, "sentaurus_segments": sentaurus_segments,
                "accepted_checkpoints": [], "failed_transition": None, "interpolation": "forbidden",
                "branch_threshold_version": "v1: multiplication=[0.1,10], leakage<=1e-3"}
    manifest_path = root / "sweep_manifest.json"
    _write_json(manifest_path, manifest)
    return manifest_path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--template", type=Path, default=REPO / "reference_tcad" / "pn2d_sentaurus2018_minimal6" / "vela" / "pn2d_minimal6_sweep_template.json")
    parser.add_argument("--state-root", type=Path)
    parser.add_argument("--vela-runner", type=Path)
    args = parser.parse_args()
    args.out_dir = args.out_dir.resolve()
    manifest_path = initialise_package(args.out_dir, args.template, authoritative_state_root=args.state_root)
    if args.vela_runner is not None:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        execute_segments(manifest, args.out_dir, lambda segment: run_vela_subprocess_segment(args.out_dir, args.vela_runner, segment))
        validate_sweep_manifest(manifest)
        _write_json(manifest_path, manifest)
    return 0

if __name__ == "__main__":
    raise SystemExit(main())