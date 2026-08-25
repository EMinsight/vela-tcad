#!/usr/bin/env python3
"""Export matched DD/DG Sentaurus Id-Vg states in three operating regimes."""

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
    / "three_regime_spatial_oracles_20260824"
)
REMOTE_ROOT = "~/sentaurus_runs/vela_oracle/transportmodels_three_regime_spatial_20260824"
IMPORTER = REPO_ROOT / "build-release/sentaurus_import.exe"
MANIFEST = OUTPUT_ROOT / "three_regime_spatial_oracles_manifest.json"

# The official sweep starts at -1 V and ends at 2.2 V.  Normalized
# quasistationary time is therefore (Vg + 1) / 3.2.
STATES = (
    ("deep_off", -1.00, 0.00),
    ("threshold", 0.12, 0.35),
    ("on", 0.92, 0.60),
)
MODES = {
    "dd": ("pp6_des.cmd", "pp6_des.par", "dd_regime", "dd_three_regime_des.cmd"),
    "dg": ("pp7_des.cmd", "pp7_des.par", "dg_regime", "dg_three_regime_des.cmd"),
}


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
            list(argv), cwd=REPO_ROOT, check=False, text=True,
            stdout=subprocess.PIPE if capture else None,
            stderr=subprocess.STDOUT if capture else None,
        )
        if completed.returncode == 0:
            return completed.stdout or ""
        if completed.returncode != 255:
            break
    assert completed is not None
    raise subprocess.CalledProcessError(
        completed.returncode, list(argv), output=completed.stdout)


def prepare_bundle() -> Path:
    bundle = OUTPUT_ROOT / "bundle"
    bundle.mkdir(parents=True, exist_ok=True)
    shutil.copy2(SOURCE / "n1_msh.tdr", bundle / "n1_msh.tdr")
    marker = "CurrentPlot(Time=(Range=(0 1) Intervals=20))"
    times = "; ".join(f"{time:.2f}" for _, _, time in STATES)
    for mode, (deck_name, parameter_name, prefix, output_name) in MODES.items():
        shutil.copy2(SOURCE / parameter_name, bundle / parameter_name)
        deck = (SOURCE / deck_name).read_text(encoding="utf-8")
        replacement = (
            marker + "\n"
            f"       Plot(FilePrefix=\"{prefix}\" NoOverWrite Time=({times}))"
        )
        if marker not in deck:
            raise RuntimeError(f"Expected gate-sweep marker not found in {deck_name}")
        (bundle / output_name).write_text(
            deck.replace(marker, replacement, 1), encoding="utf-8")
    return bundle


def slug(value: float) -> str:
    return ("m" if value < 0 else "p") + f"{abs(value):.2f}".replace(".", "p")


def export_states(raw_bundle: Path) -> list[dict[str, object]]:
    artifacts: list[dict[str, object]] = []
    for mode, (_, _, prefix, _) in MODES.items():
        tdrs = sorted(raw_bundle.glob(f"{prefix}*.tdr"))
        if len(tdrs) != len(STATES):
            raise RuntimeError(
                f"Expected {len(STATES)} {mode} TDRs, found {len(tdrs)}: "
                f"{[path.name for path in tdrs]}")
        for (regime, bias, time), tdr in zip(STATES, tdrs):
            export_dir = OUTPUT_ROOT / "exports" / mode / f"vg_{slug(bias)}"
            run([str(IMPORTER), "--tdr", str(tdr), "--export-dir", str(export_dir)])
            artifacts.append({
                "mode": mode,
                "regime": regime,
                "gate_bias_V": bias,
                "drain_bias_V": 1.1,
                "quasistationary_time": time,
                "tdr": str(tdr.resolve()),
                "tdr_sha256": sha256(tdr),
                "export_dir": str(export_dir.resolve()),
                "field_manifest": str((export_dir / "field_manifest.json").resolve()),
            })
    return artifacts


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--live", action="store_true")
    parser.add_argument("--export-only", action="store_true")
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--ssh-target", default="sentaurus")
    parser.add_argument("--ssh-bin", default=executable("ssh"))
    parser.add_argument("--scp-bin", default=executable("scp"))
    args = parser.parse_args()

    if args.check:
        manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
        for state in manifest["states"]:
            if sha256(Path(state["tdr"])) != state["tdr_sha256"]:
                raise RuntimeError(f"Hash mismatch: {state['tdr']}")
            if not Path(state["field_manifest"]).is_file():
                raise RuntimeError(f"Missing field manifest: {state['field_manifest']}")
        print("TransportModels matched three-regime Sentaurus oracles: PASS")
        return 0

    bundle = prepare_bundle()
    raw = OUTPUT_ROOT / "raw"
    banner_path = OUTPUT_ROOT / "sentaurus_banner.txt"
    if args.live and not args.export_only:
        banner = run(
            [args.ssh_bin, args.ssh_target, "sdevice -h 2>&1 | sed -n '1,5p'"],
            capture=True,
        ).strip()
        if "T-2022.03-SP2" not in banner:
            raise RuntimeError(f"Unexpected Sentaurus release:\n{banner}")
        run([args.ssh_bin, args.ssh_target, f"mkdir -p {REMOTE_ROOT}"])
        run([args.scp_bin, "-r", str(bundle), f"{args.ssh_target}:{REMOTE_ROOT}/"])
        for mode, (_, _, _, deck_name) in MODES.items():
            print(f"running {mode.upper()} three-regime spatial oracle", flush=True)
            run([
                args.ssh_bin, args.ssh_target,
                f"cd {REMOTE_ROOT}/bundle && sdevice {deck_name} > run_{mode}.out 2>&1",
            ])
        archive_name = "three_regime_spatial_oracles_results.tgz"
        run([args.ssh_bin, args.ssh_target,
             f"cd {REMOTE_ROOT} && tar -czf {archive_name} bundle"])
        raw.mkdir(parents=True, exist_ok=True)
        archive = raw / archive_name
        run([args.scp_bin, f"{args.ssh_target}:{REMOTE_ROOT}/{archive_name}", str(archive)])
        with tarfile.open(archive, "r:gz") as stream:
            stream.extractall(raw, filter="data")
        banner_path.write_text(banner + "\n", encoding="utf-8")
    elif not (raw / "bundle").is_dir():
        print("Prepared bundle. Use --live to run Sentaurus.")
        return 0

    states = export_states(raw / "bundle")
    manifest = {
        "schema": "vela.transportmodels.sentaurus_three_regime_spatial.v1",
        "as_of": "2026-08-24",
        "status": "complete",
        "sentaurus_banner": banner_path.read_text(encoding="utf-8").strip()
        if banner_path.exists() else None,
        "states": states,
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "complete", "states": len(states)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
