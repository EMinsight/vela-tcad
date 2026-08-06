#!/usr/bin/env python3
"""Run one exact-checkpoint Sentaurus BVmethods IIC branch on the local VM."""

from __future__ import annotations

import argparse
import re
import subprocess
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = REPO / "build-release/reference_tcad/bvmethods_sentaurus2018/run01/full_raw"
DEFAULT_OUTPUT = (
    REPO
    / "build-release/reference_tcad/bvmethods_sentaurus2018/run01"
    / "sentaurus_iic_multibias_exact_20260803"
)


def bias_tag(value: float) -> str:
    return f"iic_v{value:.6f}".replace("-", "m").replace(".", "p")


def run(command: list[str], *, cwd: Path | None = None) -> None:
    print("::", " ".join(command), flush=True)
    subprocess.run(command, cwd=cwd, check=True)


def exact_solve_block(biases: list[float]) -> str:
    lines = [
        "Solve {",
        "  Coupled(Iterations=100) { Poisson }",
        "  Coupled { Poisson Electron Hole }",
    ]
    previous = 0.0
    for bias in biases:
        delta = bias - previous
        # Quasistationary steps are normalized to the voltage span of this
        # segment.  Express controls as absolute 10 mV initial, 10 uV minimum,
        # and 100 mV maximum voltage steps.
        initial = min(0.01 / delta, 1.0)
        minimum = min(1.0e-5 / delta, initial)
        maximum = min(0.1 / delta, 1.0)
        lines.extend(
            [
                "  Quasistationary(",
                f"    InitialStep={initial:.17g} Increment=1.2",
                f"    MinStep={minimum:.17g} MaxStep={maximum:.17g}",
                f'    Goal {{ Name="drain" Voltage={bias:.17g} }}',
                "  ) {",
                "    Coupled { Poisson Electron Hole }",
                "    CurrentPlot(Time=(1))",
                f'    Plot(FilePrefix="{bias_tag(bias)}" Time=(1) NoOverWrite)',
                "  }",
            ]
        )
        previous = bias
    lines.append("}")
    return "\n".join(lines)


def build_deck(template: str, biases: list[float]) -> str:
    text = template
    text = text.replace('Plot      = "n4_des.tdr"', 'Plot      = "iic_multibias_des.tdr"')
    text = text.replace('Current   = "n4_des.plt"', 'Current   = "iic_multibias_des.plt"')
    text = text.replace('Output    = "n4_des.log"', 'Output    = "iic_multibias_des.log"')
    text = text.replace(
        "  eDensity hDensity",
        "  eDensity hDensity eQuasiFermi hQuasiFermi eMobility hMobility\n"
        "  eEparallel hEparallel",
    )
    text = text.replace(
        "  ComputeIonizationIntegrals",
        "  ComputeIonizationIntegrals(WriteAll)",
    )
    text, count = re.subn(r"\s*BreakAtIonIntegral\s*\(\s*3\s+1\.\s*\)\s*", "\n", text)
    if count != 1:
        raise ValueError("expected one BreakAtIonIntegral(3 1.) control")
    text, count = re.subn(r"Solve\s*\{[\s\S]*?\}\s*$", exact_solve_block(biases) + "\n", text)
    if count != 1:
        raise ValueError("expected one terminal Solve block")
    return text


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--biases", default="1,2,4,5,6,6.32,6.34,6.36,6.37,6.38,6.39,6.4")
    parser.add_argument("--ssh-target", default="sentaurus")
    parser.add_argument("--host-name", default="192.168.119.130")
    parser.add_argument(
        "--remote-dir",
        default="/home/tcad/sentaurus_runs/vela_oracle/bvmethods_iic_multibias_exact_20260803",
    )
    parser.add_argument(
        "--ssh-bin", type=Path,
        default=Path(r"C:\Windows\System32\OpenSSH\ssh.exe"),
    )
    parser.add_argument(
        "--scp-bin", type=Path,
        default=Path(r"C:\Windows\System32\OpenSSH\scp.exe"),
    )
    args = parser.parse_args()

    biases = [float(item) for item in args.biases.split(",") if item.strip()]
    if not biases or biases != sorted(set(biases)) or biases[0] <= 0.0:
        raise ValueError("biases must be unique, positive, and ascending")

    source = args.source_dir.resolve()
    output = args.out_dir.resolve()
    bundle = output / "bundle"
    raw = output / "raw"
    bundle.mkdir(parents=True, exist_ok=True)
    raw.mkdir(parents=True, exist_ok=True)

    deck = build_deck((source / "pp4_des.cmd").read_text(encoding="utf-8"), biases)
    deck_path = bundle / "iic_multibias.cmd"
    deck_path.write_text(deck, encoding="utf-8")
    for name in ("n1_msh.tdr", "pp4_des.par"):
        (bundle / name).write_bytes((source / name).read_bytes())

    ssh_common = ["-o", f"HostName={args.host_name}"]
    run([str(args.ssh_bin), *ssh_common, args.ssh_target, f"mkdir -p {args.remote_dir}"])
    run([
        str(args.scp_bin), *ssh_common,
        *(str(path) for path in sorted(bundle.iterdir())),
        f"{args.ssh_target}:{args.remote_dir}/",
    ])
    remote_command = (
        f"cd {args.remote_dir} && sdevice iic_multibias.cmd && "
        "tar -czf iic_multibias_results.tgz iic_v*_des.tdr "
        "iic_multibias_des.tdr iic_multibias_des.plt iic_multibias_des.log "
        "iic_multibias.cmd pp4_des.par"
    )
    run([str(args.ssh_bin), *ssh_common, args.ssh_target, remote_command])
    run([
        str(args.scp_bin), *ssh_common,
        f"{args.ssh_target}:{args.remote_dir}/iic_multibias_results.tgz",
        str(raw),
    ])
    print(f":: archive {raw / 'iic_multibias_results.tgz'}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
