#!/usr/bin/env python3
"""Run the aligned 21-point PN2D avalanche-off spatial SRH baseline."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import subprocess
from pathlib import Path
from typing import Any, Sequence


REPO = Path(__file__).resolve().parents[1]
DEFAULT_RUNNER = REPO / "build-release" / (
    "vela_example_runner.exe" if os.name == "nt" else "vela_example_runner"
)
DEFAULT_IMPORTER = REPO / "build-release" / (
    "sentaurus_import.exe" if os.name == "nt" else "sentaurus_import"
)
DEFAULT_SENT_RUN = (
    REPO
    / "build-release"
    / "reference_tcad"
    / "pn2d_sentaurus2018"
    / "sentaurus_vm_runs"
    / "pn2d_bv_sentaurus_aligned_off_20260728"
)
DEFAULT_SENT_CURVE = (
    REPO
    / "build-release"
    / "pn2d-reverse-sentaurus-scheduled-20260728"
    / "sentaurus_on_off_refresh_integer.csv"
)
DEFAULT_BASE_CONFIG = (
    REPO
    / "build-release"
    / "pn2d-bv-off-rootcause-20260728"
    / "source_scale_fix_full20"
    / "cases"
    / "source_aware_rows_qf_reference"
    / "source_aware_rows_qf_reference.json"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run_checked(command: Sequence[str], cwd: Path, log_prefix: Path) -> None:
    completed = subprocess.run(
        list(command), cwd=cwd, text=True, capture_output=True, check=False
    )
    log_prefix.with_suffix(".stdout.log").write_text(
        completed.stdout, encoding="utf-8"
    )
    log_prefix.with_suffix(".stderr.log").write_text(
        completed.stderr, encoding="utf-8"
    )
    if completed.returncode:
        raise RuntimeError(
            f"command failed ({completed.returncode}): {' '.join(command)}\n"
            f"{completed.stderr[-2000:]}"
        )


def git_metadata() -> dict[str, Any]:
    commit = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=REPO,
        text=True,
        capture_output=True,
        check=True,
    ).stdout.strip()
    dirty = bool(
        subprocess.run(
            ["git", "status", "--porcelain"],
            cwd=REPO,
            text=True,
            capture_output=True,
            check=True,
        ).stdout.strip()
    )
    return {"commit": commit, "dirty": dirty}


def configure(base: dict[str, Any], out_dir: Path) -> tuple[dict[str, Any], Path]:
    cfg = json.loads(json.dumps(base))
    solver = cfg["solver"]
    basis = solver["mobility"].get("doping_concentration_basis")
    impact = solver["impact_ionization"].get("model")
    if basis != "net_doping":
        raise RuntimeError(f"BV baseline must use net_doping, got {basis!r}")
    if impact != "none":
        raise RuntimeError(f"avalanche-off baseline requires model none, got {impact!r}")
    if "srh" not in solver.get("recombination", []):
        raise RuntimeError("SRH must remain enabled")
    if solver.get("bandgap_narrowing") != "old_slotboom":
        raise RuntimeError("Old Slotboom BGN must remain enabled")
    # This is a diagnostic acceptance gate, not a physical-model or Newton
    # residual-tolerance change.  It prevents an otherwise converged point
    # from being accepted when the Task 1 source/contact contract is missed.
    closure = solver.setdefault("global_continuity_closure", {})
    closure["mode"] = "enforce"
    closure["tolerance"] = 1.0e-5

    cfg["output_csv"] = str((out_dir / "vela_bv_off.csv").resolve())
    sweep = cfg["sweep"]
    sweep["bias_points"] = [float(-value) for value in range(21)]
    sweep["write_vtk"] = True
    sweep["vtk_prefix"] = str((out_dir / "vtk" / "dc_sweep").resolve())
    diagnostics = sweep.setdefault("diagnostics", {})
    diagnostics["newton_history"] = {
        "enabled": True,
        "csv_file": str((out_dir / "newton_history.csv").resolve()),
    }
    diagnostics["continuity_balance"] = {
        "enabled": True,
        "contacts": ["Anode", "Cathode"],
        "csv_file": str((out_dir / "continuity_balance.csv").resolve()),
    }
    diagnostics["terminal_balance"] = {
        "enabled": True,
        "contacts": ["Anode", "Cathode"],
        "csv_file": str((out_dir / "terminal_balance.csv").resolve()),
    }
    diagnostics["terminal_current_method_compare"] = {
        "enabled": True,
        "contacts": ["Anode", "Cathode"],
        "csv_file": str((out_dir / "terminal_current_method.csv").resolve()),
    }
    path = out_dir / "simulation.json"
    return cfg, path


def sentaurus_exact_points(path: Path) -> set[int]:
    result: set[int] = set()
    with path.open(newline="", encoding="utf-8-sig") as stream:
        for row in csv.DictReader(stream):
            bias = abs(float(row["bias_V"]))
            integer = int(round(bias))
            if abs(bias - integer) <= 1.0e-8:
                result.add(integer)
    return result


def import_sentaurus(
    importer: Path,
    source: Path,
    out_root: Path,
    skip_existing: bool,
    intervals: int,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    out_root.mkdir(parents=True, exist_ok=True)
    for reverse_bias in range(21):
        index = int(round(reverse_bias * intervals / 20.0))
        tdr = source / f"pn2d_bv_multibias_{index:04d}_des.tdr"
        token = "0" if reverse_bias == 0 else f"-{reverse_bias}"
        export_dir = out_root / f"sentaurus_{token}v"
        if not tdr.is_file():
            raise FileNotFoundError(tdr)
        command = [
            str(importer),
            "--tdr",
            str(tdr),
            "--export-dir",
            str(export_dir),
            "--compensated-doping-policy",
            "reported",
        ]
        if not (skip_existing and (export_dir / "field_manifest.json").is_file()):
            run_checked(command, REPO, out_root / f"import_{reverse_bias:02d}")
        rows.append(
            {
                "bias_V": -reverse_bias,
                "tdr": str(tdr),
                "tdr_sha256": sha256(tdr),
                "export_dir": str(export_dir),
                "command": command,
            }
        )
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-config", type=Path, default=DEFAULT_BASE_CONFIG)
    parser.add_argument("--runner", type=Path, default=DEFAULT_RUNNER)
    parser.add_argument("--tdr-importer", type=Path, default=DEFAULT_IMPORTER)
    parser.add_argument(
        "--sentaurus-source", type=Path, default=DEFAULT_SENT_RUN / "source"
    )
    parser.add_argument("--sentaurus-curve", type=Path, default=DEFAULT_SENT_CURVE)
    parser.add_argument("--sentaurus-intervals", type=int, default=200)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--skip-vela-run", action="store_true")
    parser.add_argument("--skip-existing-imports", action="store_true")
    parser.add_argument("--skip-report", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    out_dir = args.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "vtk").mkdir(exist_ok=True)
    base = json.loads(args.base_config.read_text(encoding="utf-8-sig"))
    cfg, config_path = configure(base, out_dir)
    config_path.write_text(json.dumps(cfg, indent=2) + "\n", encoding="utf-8")

    if not args.skip_vela_run:
        run_checked(
            [str(args.runner.resolve()), "--config", str(config_path)],
            out_dir,
            out_dir / "vela_runner",
        )
    exact_sentaurus = sentaurus_exact_points(args.sentaurus_curve)
    if exact_sentaurus not in (set(range(1, 21)), set(range(21))):
        raise RuntimeError(
            "Sentaurus curve lacks exact nonzero integer-bias coverage: "
            f"{sorted(exact_sentaurus)}"
        )
    exports = import_sentaurus(
        args.tdr_importer.resolve(),
        args.sentaurus_source.resolve(),
        out_dir / "sentaurus_exports",
        args.skip_existing_imports,
        args.sentaurus_intervals,
    )
    manifest = {
        "schema": "vela.pn2d_bv_off_srh_spatial_audit.v1",
        "git": git_metadata(),
        "frozen_configuration": {
            "bv_doping_concentration_basis": "net_doping",
            "forward_iv_doping_concentration_basis": "cell_reconstructed_total_impurity",
            "impact_ionization": "none",
            "recombination": ["srh"],
            "bandgap_narrowing": "old_slotboom",
            "diagnostic_global_continuity_closure_tolerance": 1.0e-5,
            "bias_points_V": [float(-value) for value in range(21)],
        },
        "hashes": {
            "base_config": sha256(args.base_config),
            "generated_config": sha256(config_path),
            "mesh": sha256(Path(cfg["mesh_file"])),
            "doping": sha256(Path(cfg["node_doping_file"])),
            "materials": sha256(Path(cfg["materials_file"])),
            "sentaurus_curve": sha256(args.sentaurus_curve),
            "sentaurus_multibias_intervals": args.sentaurus_intervals,
        },
        "sentaurus_exports": exports,
    }
    (out_dir / "run_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    if not args.skip_report:
        report = REPO / "scripts" / "build_pn2d_bv_off_srh_spatial_report.py"
        command = [
            os.fspath(Path(os.sys.executable)),
            str(report),
            "--mesh",
            cfg["mesh_file"],
            "--vela-curve",
            cfg["output_csv"],
            "--sentaurus-curve",
            str(args.sentaurus_curve.resolve()),
            "--terminal-balance",
            str(out_dir / "terminal_balance.csv"),
            "--vtk-dir",
            str(out_dir / "vtk"),
            "--sentaurus-export-root",
            str(out_dir / "sentaurus_exports"),
            "--out-dir",
            str(out_dir / "report"),
        ]
        run_checked(command, REPO, out_dir / "report_builder")
    print(json.dumps({"out_dir": str(out_dir), "manifest": str(out_dir / "run_manifest.json")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
