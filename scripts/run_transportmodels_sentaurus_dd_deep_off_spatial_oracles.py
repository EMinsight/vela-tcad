#!/usr/bin/env python3
"""Export the three Sentaurus DD deep-off Id-Vg spatial states.

The deck is derived mechanically from the official TransportModels pp6 DD
example.  Physics, Math, and the bias sequence are unchanged.  Explicit Plot
snapshots are added at Vg=-1.00, -0.84, and -0.68 V so each immutable state can
be replayed through Vela's production formulas.
"""

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
    / "dd_deep_off_spatial_oracles_20260824"
)
REMOTE_ROOT = "~/sentaurus_runs/vela_oracle/transportmodels_dd_deep_off_spatial_oracles_20260824"
IMPORTER = REPO_ROOT / "build-release/sentaurus_import.exe"
MANIFEST = OUTPUT_ROOT / "dd_deep_off_spatial_oracles_manifest.json"

TIMES = (0.00, 0.05, 0.10)
BIASES = (-1.00, -0.84, -0.68)


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
    attempts = 3 if Path(argv[0]).name.lower() in {"ssh.exe", "scp.exe", "ssh", "scp"} else 1
    completed: subprocess.CompletedProcess[str] | None = None
    for _ in range(attempts):
        completed = subprocess.run(
            list(argv),
            cwd=REPO_ROOT,
            check=False,
            text=True,
            stdout=subprocess.PIPE if capture else None,
            stderr=subprocess.STDOUT if capture else None,
        )
        if completed.returncode == 0:
            return completed.stdout or ""
        if completed.returncode != 255:
            break
    assert completed is not None
    raise subprocess.CalledProcessError(
        completed.returncode, list(argv), output=completed.stdout
    )


def prepare_bundle() -> Path:
    bundle = OUTPUT_ROOT / "bundle"
    bundle.mkdir(parents=True, exist_ok=True)
    for name in ("n1_msh.tdr", "pp6_des.par"):
        shutil.copy2(SOURCE / name, bundle / name)

    deck = (SOURCE / "pp6_des.cmd").read_text(encoding="utf-8")
    gate_marker = "CurrentPlot(Time=(Range=(0 1) Intervals=20))"
    gate_replacement = (
        "CurrentPlot(Time=(Range=(0 1) Intervals=20))\n"
        "       Plot(FilePrefix=\"deepoff_sweep\" NoOverWrite Time=(0.00; 0.05; 0.10))"
    )
    if gate_marker not in deck:
        raise RuntimeError("Expected pp6 gate-sweep CurrentPlot marker not found")
    deck = deck.replace(gate_marker, gate_replacement, 1)
    (bundle / "dd_deep_off_spatial_des.cmd").write_text(deck, encoding="utf-8")
    return bundle


def slug(value: float) -> str:
    sign = "m" if value < 0 else "p"
    return f"vg_{sign}{abs(value):.2f}".replace(".", "p")


def export_one(tdr: Path, export_dir: Path) -> dict[str, str]:
    run([str(IMPORTER), "--tdr", str(tdr), "--export-dir", str(export_dir)])
    return {
        "tdr": str(tdr.resolve()),
        "tdr_sha256": sha256(tdr),
        "export_dir": str(export_dir.resolve()),
        "field_manifest": str((export_dir / "field_manifest.json").resolve()),
    }


def export_tdrs(raw_bundle: Path) -> list[dict[str, object]]:
    exports = OUTPUT_ROOT / "exports"
    exports.mkdir(parents=True, exist_ok=True)
    tdrs = sorted(raw_bundle.glob("deepoff_sweep*.tdr"))
    if len(tdrs) != len(BIASES):
        raise RuntimeError(
            f"Expected {len(BIASES)} DD deep-off TDRs, found {len(tdrs)} "
            f"({[path.name for path in tdrs]})"
        )
    states: list[dict[str, object]] = []
    for time, bias, tdr in zip(TIMES, BIASES, tdrs):
        entry: dict[str, object] = {"gate_bias_V": bias, "time": time}
        entry.update(export_one(tdr, exports / slug(bias)))
        states.append(entry)
    return states


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
        for artifact in manifest["dd_deep_off_states"]:
            if sha256(Path(artifact["tdr"])) != artifact["tdr_sha256"]:
                raise RuntimeError(f"Hash mismatch: {artifact['tdr']}")
            if not Path(artifact["field_manifest"]).is_file():
                raise RuntimeError(f"Missing field manifest: {artifact['field_manifest']}")
        print("TransportModels DD deep-off Sentaurus spatial oracles check: PASS")
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
        print("running DD deep-off Id-Vg spatial oracle", flush=True)
        run(
            [
                args.ssh_bin,
                args.ssh_target,
                f"cd {REMOTE_ROOT}/bundle && sdevice dd_deep_off_spatial_des.cmd > run_deepoff.out 2>&1",
            ]
        )
        archive_name = "dd_deep_off_spatial_oracles_results.tgz"
        run([args.ssh_bin, args.ssh_target, f"cd {REMOTE_ROOT} && tar -czf {archive_name} bundle"])
        raw.mkdir(parents=True, exist_ok=True)
        archive = raw / archive_name
        run([args.scp_bin, f"{args.ssh_target}:{REMOTE_ROOT}/{archive_name}", str(archive)])
        with tarfile.open(archive, "r:gz") as stream:
            stream.extractall(raw, filter="data")
        banner_file.write_text(banner + "\n", encoding="utf-8")
    elif not (raw / "bundle").is_dir():
        print("Prepared DD deep-off spatial-oracle bundle. Use --live to execute it.")
        return 0

    states = export_tdrs(raw / "bundle")
    manifest = {
        "schema": "vela.transportmodels.sentaurus_dd_deep_off_spatial_oracles.v1",
        "as_of": "2026-08-24",
        "status": "complete",
        "sentaurus_banner": banner_file.read_text(encoding="utf-8").strip()
        if banner_file.exists()
        else None,
        "fixed_drain_bias_V": 1.1,
        "dd_deep_off_states": states,
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "complete", "deep_off_states": len(states)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
