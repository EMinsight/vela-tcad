#!/usr/bin/env python3
"""Run the frozen silicon Schottky oracle on an SSH-accessible Sentaurus VM.

Dry-run is the default.  The source fixture deliberately excludes image-force
lowering, tunnelling, series resistance, high-field transport, AC, and heat.
"""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import tarfile
from datetime import datetime
from pathlib import Path
from typing import Sequence

from run_singledevice_sentaurus_vm import (
    capture,
    executable,
    file_hashes,
    parse_plt,
    run,
)


REPO = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = REPO / "reference_tcad" / "schottky_charon_sentaurus2018" / "source"
DEFAULT_OUTPUT = (
    REPO / "build-release" / "reference_tcad" / "schottky_sentaurus2022"
    / "sentaurus_vm_runs"
)
DEFAULT_REMOTE_ROOT = "~/sentaurus_runs/vela_oracle_2022"
REQUIRED_FILES = ("schottky_n_sde.cmd", "schottky_n_des.cmd")


def prepare_bundle(source: Path, bundle: Path) -> None:
    missing = [name for name in REQUIRED_FILES if not (source / name).is_file()]
    if missing:
        raise FileNotFoundError(f"missing Schottky source files: {', '.join(missing)}")
    bundle.mkdir(parents=True, exist_ok=True)
    for name in REQUIRED_FILES:
        shutil.copy2(source / name, bundle / name)


def remote_commands(remote_dir: str) -> list[str]:
    return [
        f"cd {remote_dir} && sde -e -l schottky_n_sde.cmd > run_sde.out 2>&1",
        f"cd {remote_dir} && sdevice schottky_n_des.cmd > run_sdevice.out 2>&1",
        (
            f"cd {remote_dir} && tar -czf schottky_results.tgz "
            "*.cmd *.tdr *.plt *.log *.out"
        ),
    ]


def extract_curve(path: Path, output: Path) -> None:
    datasets, rows = parse_plt(path)
    bias = datasets.index("anode OuterVoltage")
    current = datasets.index("anode TotalCurrent")
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["bias_V", "current_total_A_per_um"])
        writer.writerows((row[bias], abs(row[current])) for row in rows)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--local-output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--remote-root", default=DEFAULT_REMOTE_ROOT)
    parser.add_argument("--run-id", default=None)
    parser.add_argument("--ssh-target", default="sentaurus")
    parser.add_argument("--ssh-bin", default=executable("ssh"))
    parser.add_argument("--scp-bin", default=executable("scp"))
    parser.add_argument("--sentaurus-version", default=None)
    parser.add_argument("--live", action="store_true")
    args = parser.parse_args(argv)

    run_id = args.run_id or datetime.now().strftime("schottky_vm_%Y%m%d_%H%M%S")
    local = args.local_output_dir.resolve() / run_id
    bundle = local / "bundle"
    raw = local / "raw"
    remote = f"{args.remote_root.rstrip('/')}/{run_id}"
    prepare_bundle(args.source_dir.resolve(), bundle)
    commands = remote_commands(remote)
    sentaurus_banner = None

    if args.live:
        raw.mkdir(parents=True, exist_ok=True)
        sentaurus_banner = capture([
            args.ssh_bin, args.ssh_target, "sdevice -h 2>&1 | sed -n '1,5p'",
        ]).strip()
        if args.sentaurus_version and args.sentaurus_version not in sentaurus_banner:
            raise RuntimeError(
                "live SDevice banner does not contain expected release "
                f"{args.sentaurus_version!r}:\n{sentaurus_banner}"
            )
        run([args.ssh_bin, args.ssh_target, f"mkdir -p {remote}"])
        run([
            args.scp_bin,
            *(str(path) for path in sorted(bundle.iterdir())),
            f"{args.ssh_target}:{remote}/",
        ])
        for command in commands:
            run([args.ssh_bin, args.ssh_target, command])
        archive = raw / "schottky_results.tgz"
        run([args.scp_bin, f"{args.ssh_target}:{remote}/{archive.name}", str(archive)])
        with tarfile.open(archive, "r:gz") as stream:
            stream.extractall(raw, filter="data")
        extract_curve(
            raw / "forward_schottky_n_iv_des.plt",
            local / "schottky_forward.csv",
        )

    manifest = {
        "schema": "vela.schottky_sentaurus_vm_run.v1",
        "run_id": run_id,
        "dry_run": not args.live,
        "source_dir": str(args.source_dir.resolve()),
        "local_run_dir": str(local),
        "remote_dir": remote,
        "ssh_target": args.ssh_target,
        "expected_sentaurus_version": args.sentaurus_version,
        "sentaurus_banner": sentaurus_banner,
        "required_files": list(REQUIRED_FILES),
        "source_sha256": file_hashes(args.source_dir.resolve(), REQUIRED_FILES),
        "bundle_sha256": file_hashes(bundle, REQUIRED_FILES),
        "commands": commands,
    }
    local.mkdir(parents=True, exist_ok=True)
    (local / "sentaurus_vm_run_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
