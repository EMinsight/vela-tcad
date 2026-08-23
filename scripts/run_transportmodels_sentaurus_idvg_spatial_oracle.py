#!/usr/bin/env python3
"""Export five-bias Sentaurus DG Id-Vg spatial TDR oracle states."""

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
    / "idvg_spatial_oracle_20260821"
)
REMOTE_ROOT = "~/sentaurus_runs/vela_oracle/transportmodels_idvg_spatial_oracle_20260821"
IMPORTER = REPO_ROOT / "build-release/sentaurus_import.exe"
MANIFEST = OUTPUT_ROOT / "spatial_oracle_manifest.json"
TIMES = (0.25, 0.30, 0.35, 0.40, 0.625)
GATE_BIASES = (-0.20, -0.04, 0.12, 0.28, 1.00)


def executable(name: str) -> str:
    if os.name == "nt":
        candidate = Path(os.environ.get("SystemRoot", r"C:\Windows")) / "System32/OpenSSH" / f"{name}.exe"
        if candidate.is_file():
            return str(candidate)
    return shutil.which(name) or name


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def run(argv: Sequence[str], *, capture: bool = False) -> str:
    completed = subprocess.run(
        list(argv), cwd=REPO_ROOT, check=True, text=True,
        stdout=subprocess.PIPE if capture else None,
        stderr=subprocess.STDOUT if capture else None,
    )
    return completed.stdout or ""


def prepare_bundle() -> Path:
    bundle = OUTPUT_ROOT / "bundle"
    bundle.mkdir(parents=True, exist_ok=True)
    for name in ("n1_msh.tdr", "pp7_des.par"):
        shutil.copy2(SOURCE / name, bundle / name)
    deck = (SOURCE / "pp7_des.cmd").read_text(encoding="utf-8")
    current_plot = "CurrentPlot(Time=(Range=(0 1) Intervals=20))"
    replacement = (
        "CurrentPlot(Time=(Range=(0 1) Intervals=40))\n"
        "       Plot(FilePrefix=\"spatial\" NoOverWrite "
        "Time=(0.25; 0.30; 0.35; 0.40; 0.625))"
    )
    if current_plot not in deck:
        raise RuntimeError("Expected CurrentPlot line not found")
    deck = deck.replace(current_plot, replacement, 1)
    (bundle / "idvg_spatial_des.cmd").write_text(deck, encoding="utf-8")
    return bundle


def export_tdrs(raw_bundle: Path) -> list[dict[str, str]]:
    exports = OUTPUT_ROOT / "exports"
    exports.mkdir(parents=True, exist_ok=True)
    tdrs = sorted(raw_bundle.glob("spatial*.tdr"))
    if len(tdrs) != len(GATE_BIASES):
        raise RuntimeError(f"Expected {len(GATE_BIASES)} spatial TDRs, found {len(tdrs)}")
    artifacts = []
    for bias, tdr in zip(GATE_BIASES, tdrs):
        slug = ("m" if bias < 0 else "p") + f"{abs(bias):.2f}".replace(".", "p")
        export = exports / f"vg_{slug}"
        run([str(IMPORTER), "--tdr", str(tdr), "--export-dir", str(export)])
        artifacts.append(
            {
                "gate_bias_V": str(bias),
                "time": str(TIMES[len(artifacts)]),
                "tdr": str(tdr.resolve()),
                "tdr_sha256": sha256(tdr),
                "export_dir": str(export.resolve()),
                "field_manifest": str((export / "field_manifest.json").resolve()),
            }
        )
    return artifacts


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
        for artifact in manifest["states"]:
            assert sha256(Path(artifact["tdr"])) == artifact["tdr_sha256"]
            assert Path(artifact["field_manifest"]).is_file()
        print("TransportModels Sentaurus spatial oracle check: PASS")
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
        print("[1/1] running default DG five-state spatial oracle", flush=True)
        run(
            [
                args.ssh_bin,
                args.ssh_target,
                f"cd {REMOTE_ROOT}/bundle && sdevice idvg_spatial_des.cmd > run_sdevice.out 2>&1",
            ]
        )
        archive_name = "idvg_spatial_oracle_results.tgz"
        run([args.ssh_bin, args.ssh_target, f"cd {REMOTE_ROOT} && tar -czf {archive_name} bundle"])
        raw.mkdir(parents=True, exist_ok=True)
        archive = raw / archive_name
        run([args.scp_bin, f"{args.ssh_target}:{REMOTE_ROOT}/{archive_name}", str(archive)])
        with tarfile.open(archive, "r:gz") as stream:
            stream.extractall(raw, filter="data")
        banner_file.write_text(banner + "\n", encoding="utf-8")
    elif not (raw / "bundle").is_dir():
        print("Prepared spatial-oracle bundle. Use --live to execute it.")
        return 0

    states = export_tdrs(raw / "bundle")
    manifest = {
        "schema": "vela.transportmodels.sentaurus_idvg_spatial_oracle.v1",
        "as_of": "2026-08-21",
        "status": "complete",
        "sentaurus_banner": banner_file.read_text(encoding="utf-8").strip() if banner_file.exists() else None,
        "fixed_drain_bias_V": 1.1,
        "gate_biases_V": list(GATE_BIASES),
        "quasistationary_times": list(TIMES),
        "states": states,
    }
    MANIFEST.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": manifest["status"], "states": states}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
