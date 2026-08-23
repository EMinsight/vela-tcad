#!/usr/bin/env python3
"""Run the frozen Sentaurus SingleDevice oracle on an SSH-accessible VM.

Dry-run is the default.  A live run requires an SSH configuration or agent that
can authenticate without embedding credentials in this repository.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
import tarfile
from datetime import datetime
from pathlib import Path
from typing import Sequence


REPO = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = REPO / "reference_tcad" / "singledevice_sentaurus2018" / "source"
DEFAULT_OUTPUT = (
    REPO / "build-release" / "reference_tcad" / "singledevice_sentaurus2018"
    / "sentaurus_vm_runs"
)
DEFAULT_REMOTE_ROOT = "~/sentaurus_runs/vela_oracle"
REQUIRED_FILES = (
    "singledevice_sde.cmd", "singledevice_sdevice.cmd", "Silicon.par", "sdevice.par"
)


def executable(name: str) -> str:
    if os.name == "nt":
        path = Path(os.environ.get("SystemRoot", r"C:\Windows")) / "System32" / "OpenSSH" / f"{name}.exe"
        if path.is_file():
            return str(path)
    return shutil.which(name) or name


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def file_hashes(root: Path, names: Sequence[str]) -> dict[str, str]:
    return {name: sha256_file(root / name) for name in names}


def preprocess(text: str, variables: dict[str, str]) -> str:
    for name, value in variables.items():
        text = text.replace(f"@{name}@", value)
    return text


def prepare_bundle(source: Path, bundle: Path) -> None:
    missing = [name for name in REQUIRED_FILES if not (source / name).is_file()]
    if missing:
        raise FileNotFoundError(f"missing SingleDevice source files: {', '.join(missing)}")
    bundle.mkdir(parents=True, exist_ok=True)
    (bundle / "singledevice_sde.cmd").write_text(
        preprocess((source / "singledevice_sde.cmd").read_text(), {"node": "2"})
    )
    variables = {
        "previous": "1",
        "tdr": "n2_msh.tdr",
        "tdrdat": "singledevice_des.tdr",
        "parameter": "sdevice.par",
        "plot": "singledevice.plt",
        "log": "singledevice.log",
        "node": "2",
    }
    (bundle / "singledevice_sdevice.cmd").write_text(
        preprocess((source / "singledevice_sdevice.cmd").read_text(), variables)
    )
    shutil.copy2(source / "Silicon.par", bundle / "Silicon.par")
    # Current Applications Library parameter files may use the Workbench
    # ``@pwd@`` token.  The isolated bundle keeps the include beside the deck,
    # so a relative current-directory replacement is both portable and exact.
    (bundle / "sdevice.par").write_text(
        preprocess((source / "sdevice.par").read_text(), {"pwd": "."})
    )


def remote_commands(remote_dir: str) -> list[str]:
    return [
        f"cd {remote_dir} && sde -e -l singledevice_sde.cmd > run_sde.out 2>&1",
        f"cd {remote_dir} && sdevice singledevice_sdevice.cmd > run_sdevice.out 2>&1",
        (
            f"cd {remote_dir} && tar -czf singledevice_results.tgz "
            "*.cmd *.par *.tdr *.plt *.log *.out"
        ),
    ]


def parse_plt(path: Path) -> tuple[list[str], list[list[float]]]:
    # Reuse the repository's tested DF-ISE parser without duplicating syntax.
    from sentaurus_import import parse_quoted_list, parse_values_block

    text = path.read_text(errors="ignore")
    datasets = parse_quoted_list(text, "datasets")
    if not datasets:
        raise ValueError(f"{path} has no datasets block")
    return datasets, parse_values_block(text, len(datasets))


def extract_curve(path: Path, output: Path) -> None:
    datasets, rows = parse_plt(path)
    gate = datasets.index("gate OuterVoltage")
    drain = datasets.index("drain TotalCurrent")
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["bias_V", "current_total"])
        writer.writerows((row[gate], row[drain]) for row in rows)


def run(argv: Sequence[str], *, cwd: Path | None = None) -> None:
    subprocess.run(list(argv), cwd=cwd, check=True)


def capture(argv: Sequence[str], *, cwd: Path | None = None) -> str:
    completed = subprocess.run(
        list(argv), cwd=cwd, check=True, text=True,
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
    )
    return completed.stdout


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--local-output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--remote-root", default=DEFAULT_REMOTE_ROOT)
    parser.add_argument("--run-id", default=None)
    parser.add_argument("--ssh-target", default="sentaurus")
    parser.add_argument("--ssh-bin", default=executable("ssh"))
    parser.add_argument("--scp-bin", default=executable("scp"))
    parser.add_argument(
        "--sentaurus-version", default=None,
        help="expected live SDevice release, for example T-2022.03-SP2",
    )
    parser.add_argument(
        "--live", action="store_true",
        help="perform the SSH/SCP run; without this flag only prepare and print the manifest",
    )
    args = parser.parse_args(argv)

    run_id = args.run_id or datetime.now().strftime("singledevice_vm_%Y%m%d_%H%M%S")
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
            args.ssh_bin, args.ssh_target,
            "sdevice -h 2>&1 | sed -n '1,5p'",
        ]).strip()
        if args.sentaurus_version and args.sentaurus_version not in sentaurus_banner:
            raise RuntimeError(
                "live SDevice banner does not contain expected release "
                f"{args.sentaurus_version!r}:\n{sentaurus_banner}"
            )
        run([args.ssh_bin, args.ssh_target, f"mkdir -p {remote}"])
        run([args.scp_bin, *(str(path) for path in sorted(bundle.iterdir())), f"{args.ssh_target}:{remote}/"])
        for command in commands:
            run([args.ssh_bin, args.ssh_target, command])
        archive = raw / "singledevice_results.tgz"
        run([args.scp_bin, f"{args.ssh_target}:{remote}/{archive.name}", str(archive)])
        with tarfile.open(archive, "r:gz") as stream:
            stream.extractall(raw, filter="data")
        extract_curve(raw / "IdVgsLin_singledevice.plt", local / "singledevice_idvg_lin.csv")
        extract_curve(raw / "IdVgsSat_singledevice.plt", local / "singledevice_idvg_sat.csv")

    manifest = {
        "schema": "vela.singledevice_sentaurus_vm_run.v1",
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
    (local / "sentaurus_vm_run_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
