#!/usr/bin/env python3
"""Run standalone Sentaurus BVmethods IIC states at millivolt drain biases."""

from __future__ import annotations

import argparse
import re
import subprocess
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = (
    REPO
    / "build-release/reference_tcad/bvmethods_sentaurus2018/run01/full_raw"
)
DEFAULT_OUTPUT = (
    REPO
    / "build-release/reference_tcad/bvmethods_sentaurus2018/run01"
    / "sentaurus_low_bias_20260802"
)


def bias_tag(value: float) -> str:
    return f"{value:.6f}".replace("-", "m").replace(".", "p")


def run(command: list[str], *, cwd: Path | None = None) -> None:
    print("::", " ".join(command), flush=True)
    subprocess.run(command, cwd=cwd, check=True)


def ssh_host_options(host_name: str | None) -> list[str]:
    if not host_name:
        return []
    return ["-o", f"HostName={host_name}"]


def build_command(template: str, bias: float, tag: str) -> str:
    text = template
    text = text.replace('Plot      = "n4_des.tdr"', f'Plot      = "{tag}_des.tdr"')
    text = text.replace('Current   = "n4_des.plt"', f'Current   = "{tag}_des.plt"')
    text = text.replace('Output    = "n4_des.log"', f'Output    = "{tag}_des.log"')
    text = text.replace(
        "  eDensity hDensity",
        "  eDensity hDensity eQuasiFermi hQuasiFermi eMobility hMobility",
    )
    if bias == 0.0:
        text = re.sub(
            r"\n\s*Quasistationary\(\s*\n"
            r"\s*InitialStep=0\.0001 Increment=1\.41\s*\n"
            r"\s*MinStep=1e-07 MaxStep=0\.025\s*\n"
            r"\s*Goal\{ Name=\"drain\" Voltage=100\. \}\s*\n"
            r"\s*\) \{ Coupled \{ Poisson Electron Hole \} \}\s*",
            "\n",
            text,
            count=1,
        )
    else:
        text = text.replace(
            'Goal{ Name="drain" Voltage=100. }',
            f'Goal{{ Name="drain" Voltage={bias:.17g} }}',
        )
    return text


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--biases", default="0,0.001,0.002,0.005,0.01")
    parser.add_argument("--ssh-target", default="sentaurus")
    parser.add_argument(
        "--host-name",
        default=None,
        help="optional HostName override; defaults to the --ssh-target SSH config",
    )
    parser.add_argument(
        "--remote-dir",
        default="/home/tcad/sentaurus_runs/vela_oracle/bvmethods_low_bias_20260802",
    )
    parser.add_argument(
        "--ssh-bin",
        type=Path,
        default=Path(r"C:\Windows\System32\OpenSSH\ssh.exe"),
    )
    parser.add_argument(
        "--scp-bin",
        type=Path,
        default=Path(r"C:\Windows\System32\OpenSSH\scp.exe"),
    )
    args = parser.parse_args()

    source = args.source_dir.resolve()
    output = args.out_dir.resolve()
    bundle = output / "bundle"
    raw = output / "raw"
    bundle.mkdir(parents=True, exist_ok=True)
    raw.mkdir(parents=True, exist_ok=True)

    biases = [float(item) for item in args.biases.split(",") if item.strip()]
    template = (source / "pp4_des.cmd").read_text(encoding="utf-8")
    commands: list[Path] = []
    for bias in biases:
        tag = f"lowbias_{bias_tag(bias)}"
        path = bundle / f"{tag}.cmd"
        path.write_text(build_command(template, bias, tag), encoding="utf-8")
        commands.append(path)

    for name in ("n1_msh.tdr", "pp4_des.par"):
        (bundle / name).write_bytes((source / name).read_bytes())

    ssh_common = ssh_host_options(args.host_name)
    run(
        [
            str(args.ssh_bin),
            *ssh_common,
            args.ssh_target,
            f"mkdir -p {args.remote_dir}",
        ]
    )
    run(
        [
            str(args.scp_bin),
            *ssh_common,
            *(str(path) for path in sorted(bundle.iterdir())),
            f"{args.ssh_target}:{args.remote_dir}/",
        ]
    )

    command_names = " ".join(path.name for path in commands)
    remote_command = (
        f"cd {args.remote_dir} && "
        f"for cmd in {command_names}; do sdevice \"$cmd\" || exit $?; done && "
        "tar -czf lowbias_results.tgz lowbias_*_des.tdr lowbias_*_des.plt "
        "lowbias_*_des.log lowbias_*.cmd pp4_des.par"
    )
    run([str(args.ssh_bin), *ssh_common, args.ssh_target, remote_command])
    run(
        [
            str(args.scp_bin),
            *ssh_common,
            f"{args.ssh_target}:{args.remote_dir}/lowbias_results.tgz",
            str(raw),
        ]
    )
    print(f":: archive {raw / 'lowbias_results.tgz'}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
