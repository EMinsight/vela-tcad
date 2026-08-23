#!/usr/bin/env python3
"""Export the Sentaurus 2022 DD TransportModels state at Vg=-1 V, Vd=1.1 V."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import tarfile
from pathlib import Path
from typing import Sequence


REPO_ROOT = Path(__file__).resolve().parents[1]
SOURCE = REPO_ROOT / "build-release/reference_tcad/transportmodels_sentaurus2022/run02/full_raw"
OUTPUT_ROOT = (
    REPO_ROOT
    / "build-release/reference_tcad/transportmodels_sentaurus2022/sentaurus_vm_runs"
    / "dd_deep_off_oracle_20260823"
)
REMOTE_ROOT = "~/sentaurus_runs/vela_oracle/transportmodels_dd_deep_off_oracle_20260823"
IMPORTER = REPO_ROOT / "build-release/sentaurus_import.exe"
MANIFEST = OUTPUT_ROOT / "dd_deep_off_oracle_manifest.json"


def executable(name: str) -> str:
    if os.name == "nt":
        candidate = (
            Path(os.environ.get("SystemRoot", r"C:\Windows"))
            / "System32/OpenSSH"
            / f"{name}.exe"
        )
        if candidate.is_file():
            return str(candidate)
    return shutil.which(name) or name


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def run(argv: Sequence[str], *, capture: bool = False) -> str:
    completed = subprocess.run(
        list(argv),
        cwd=REPO_ROOT,
        check=True,
        text=True,
        stdout=subprocess.PIPE if capture else None,
        stderr=subprocess.STDOUT if capture else None,
    )
    return completed.stdout or ""


def prepare_bundle() -> Path:
    bundle = OUTPUT_ROOT / "bundle"
    bundle.mkdir(parents=True, exist_ok=True)
    for name in ("n1_msh.tdr", "pp6_des.par"):
        shutil.copy2(SOURCE / name, bundle / name)

    deck = (SOURCE / "pp6_des.cmd").read_text(encoding="utf-8")
    drain_solve = ") { Coupled { Poisson Electron Hole   } }"
    if deck.count(drain_solve) != 1:
        raise RuntimeError("Expected exactly one DD drain-ramp solve block")
    deck = deck.replace(
        drain_solve,
        ") { Coupled { Poisson Electron Hole   }\n"
        "       Plot(FilePrefix=\"dd_deep_off\" NoOverWrite Time=(1))\n"
        "     }",
        1,
    )
    gate_marker = "   *-  gate voltage sweep"
    if gate_marker not in deck:
        raise RuntimeError("Expected DD gate-sweep marker not found")
    deck = deck.split(gate_marker, 1)[0].rstrip() + "\n}\n"
    (bundle / "dd_deep_off_des.cmd").write_text(deck, encoding="utf-8")
    return bundle


def find_state_tdr(bundle: Path) -> Path:
    candidates = sorted(bundle.glob("dd_deep_off*_des.tdr"))
    if len(candidates) != 1:
        raise RuntimeError(
            f"Expected one DD deep-off TDR, found {len(candidates)}: {candidates}"
        )
    return candidates[0]


def export_state(raw_bundle: Path) -> dict[str, str]:
    state_tdr = find_state_tdr(raw_bundle)
    export_dir = OUTPUT_ROOT / "export"
    if export_dir.exists():
        shutil.rmtree(export_dir)
    run([str(IMPORTER), "--tdr", str(state_tdr), "--export-dir", str(export_dir)])
    return {
        "tdr": str(state_tdr.resolve()),
        "tdr_sha256": sha256(state_tdr),
        "export_dir": str(export_dir.resolve()),
        "field_manifest": str((export_dir / "field_manifest.json").resolve()),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--live", action="store_true")
    parser.add_argument("--export-only", action="store_true")
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--ssh-target", default="sentaurus")
    parser.add_argument("--ssh-bin", default=executable("ssh"))
    parser.add_argument("--scp-bin", default=executable("scp"))
    args = parser.parse_args()

    if args.check:
        manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
        state = manifest["state"]
        assert sha256(Path(state["tdr"])) == state["tdr_sha256"]
        assert Path(state["field_manifest"]).is_file()
        print("TransportModels Sentaurus DD deep-off oracle check: PASS")
        return 0

    bundle = prepare_bundle()
    raw = OUTPUT_ROOT / "raw"
    banner_file = OUTPUT_ROOT / "sentaurus_banner.txt"
    if args.live and not args.export_only:
        banner = run(
            [args.ssh_bin, args.ssh_target, "sdevice -h 2>&1 | sed -n '1,5p'"],
            capture=True,
        ).strip()
        if "T-2022.03-SP2" not in banner:
            raise RuntimeError(f"Unexpected Sentaurus release:\n{banner}")
        run([args.ssh_bin, args.ssh_target, f"mkdir -p {REMOTE_ROOT}"])
        run([args.scp_bin, "-r", str(bundle), f"{args.ssh_target}:{REMOTE_ROOT}/"])
        run(
            [
                args.ssh_bin,
                args.ssh_target,
                f"cd {REMOTE_ROOT}/bundle && "
                "sdevice dd_deep_off_des.cmd > run_sdevice.out 2>&1",
            ]
        )
        archive_name = "dd_deep_off_oracle_results.tgz"
        run(
            [
                args.ssh_bin,
                args.ssh_target,
                f"cd {REMOTE_ROOT} && tar -czf {archive_name} bundle",
            ]
        )
        raw.mkdir(parents=True, exist_ok=True)
        archive = raw / archive_name
        run(
            [
                args.scp_bin,
                f"{args.ssh_target}:{REMOTE_ROOT}/{archive_name}",
                str(archive),
            ]
        )
        with tarfile.open(archive, "r:gz") as stream:
            stream.extractall(raw, filter="data")
        banner_file.write_text(banner + "\n", encoding="utf-8")
    elif not (raw / "bundle").is_dir():
        print("Prepared DD deep-off bundle. Use --live to execute it.")
        return 0

    state = export_state(raw / "bundle")
    manifest = {
        "schema": "vela.transportmodels.sentaurus_dd_deep_off_oracle.v1",
        "as_of": "2026-08-23",
        "status": "complete",
        "sentaurus_banner": (
            banner_file.read_text(encoding="utf-8").strip()
            if banner_file.exists()
            else None
        ),
        "model": "DD",
        "gate_bias_V": -1.0,
        "drain_bias_V": 1.1,
        "source_bias_V": 0.0,
        "substrate_bias_V": 0.0,
        "state": state,
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
