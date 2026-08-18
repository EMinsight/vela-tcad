#!/usr/bin/env python3
"""Run fixed-state Sentaurus controls that identify MeanIonIntegral semantics.

The electron-only and hole-only variants preserve the converged electrostatic
and carrier solution but set one van Overstraeten-de Man coefficient to zero.
For a one-carrier coefficient alpha, Eqs. (469)-(470) give

    I = 1 - exp(-A),  A = integral(alpha ds).

Consequently an independently integrated mean coefficient gives A/2, whereas
an arithmetic mean of the final carrier integrals gives I/2.  The two controls
therefore distinguish the definitions without fitting an avalanche model.
"""

from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
from pathlib import Path, PurePosixPath


REPO = Path(__file__).resolve().parents[1]
DEFAULT_SOURCE = (
    REPO
    / "build-release/reference_tcad/bvmethods_sentaurus2018/run01/full_raw"
)
DEFAULT_OUTPUT = (
    REPO
    / "build-release/reference_tcad/bvmethods_sentaurus2018/run01"
    / "sentaurus_mean_ionization_controls_20260805"
)

VAN_OVERSTRAETEN = {
    "baseline": None,
    "electron_only": """
vanOverstraetendeMan {
  a(low)       = 7.0300e+05 , 0.0000e+00
  a(high)      = 7.0300e+05 , 0.0000e+00
  b(low)       = 1.2310e+06 , 2.0360e+06
  b(high)      = 1.2310e+06 , 1.6930e+06
  E0           = 4.0000e+05 , 4.0000e+05
  hbarOmega    = 0.063 , 0.063
}
""",
    "electron_half": """
vanOverstraetendeMan {
  a(low)       = 3.5150e+05 , 0.0000e+00
  a(high)      = 3.5150e+05 , 0.0000e+00
  b(low)       = 1.2310e+06 , 2.0360e+06
  b(high)      = 1.2310e+06 , 1.6930e+06
  E0           = 4.0000e+05 , 4.0000e+05
  hbarOmega    = 0.063 , 0.063
}
""",
    "hole_only": """
vanOverstraetendeMan {
  a(low)       = 0.0000e+00 , 1.5820e+06
  a(high)      = 0.0000e+00 , 6.7100e+05
  b(low)       = 1.2310e+06 , 2.0360e+06
  b(high)      = 1.2310e+06 , 1.6930e+06
  E0           = 4.0000e+05 , 4.0000e+05
  hbarOmega    = 0.063 , 0.063
}
""",
    "hole_half": """
vanOverstraetendeMan {
  a(low)       = 0.0000e+00 , 7.9100e+05
  a(high)      = 0.0000e+00 , 3.3550e+05
  b(low)       = 1.2310e+06 , 2.0360e+06
  b(high)      = 1.2310e+06 , 1.6930e+06
  E0           = 4.0000e+05 , 4.0000e+05
  hbarOmega    = 0.063 , 0.063
}
""",
}


def run(command: list[str]) -> None:
    subprocess.run(command, check=True)


def ssh_host_options(host_name: str | None) -> list[str]:
    if not host_name:
        return []
    return ["-o", f"HostName={host_name}"]


def safe_remote_root(value: str) -> str:
    if re.fullmatch(r"/[A-Za-z0-9._/-]+", value) is None:
        raise ValueError("remote root must be a safe absolute POSIX path")
    path = PurePosixPath(value)
    if str(path) != value or any(part in {".", ".."} for part in path.parts):
        raise ValueError("remote root must be normalized")
    return value


def deck(variant: str, drive: str) -> str:
    avalanche_drive = {
        "eparallel": "Eparallel",
        "electric_field": "ElectricField",
    }[drive]
    return f'''File {{
  Grid      = "n1_msh.tdr"
  Plot      = "mean_control_{variant}.tdr"
  Parameter = "mean_control_{variant}.par"
  Current   = "mean_control_{variant}.plt"
  Output    = "mean_control_{variant}.log"
}}

Electrode {{
  {{ Name="drain" Voltage=0.0 }}
  {{ Name="source" Voltage=0.0 }}
  {{ Name="gate" Voltage=0.0 Barrier=-0.55 }}
  {{ Name="substrate" Voltage=0.0 }}
}}

Physics {{
  EffectiveIntrinsicDensity(OldSlotboom)
  Mobility(DopingDep HighFieldsaturation(GradQuasiFermi) Enormal)
  Recombination(SRH(DopingDep) Band2Band(E2) Avalanche({avalanche_drive}))
  Fermi
}}

Plot {{
  Potential ElectricField/Vector eEparallel hEparallel
  eIonIntegral hIonIntegral MeanIonIntegral
  eAlphaAvalanche hAlphaAvalanche
}}

Math {{
  Iterations=20
  Notdamped=100
  RelErrControl
  ComputeIonizationIntegrals(WriteAll)
  AvalPostProcessing
}}

Solve {{
  Load(FilePrefix="source_n4_des")
  Coupled(Iterations=1) {{ Poisson Electron Hole }}
  Plot(FilePrefix="mean_control_{variant}" NoOverWrite)
}}
'''


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--ssh-target", default="sentaurus")
    parser.add_argument(
        "--host-name",
        default=None,
        help="optional HostName override; defaults to the --ssh-target SSH config",
    )
    parser.add_argument(
        "--drive",
        choices=("eparallel", "electric_field"),
        default="eparallel",
    )
    parser.add_argument(
        "--variants",
        default=",".join(VAN_OVERSTRAETEN),
        help="comma-separated subset of baseline,electron_only,electron_half,"
        "hole_only,hole_half",
    )
    parser.add_argument(
        "--remote-root",
        default=(
            "/home/tcad/sentaurus_runs/vela_oracle/"
            "bvmethods_mean_ionization_controls_20260805"
        ),
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
    output = args.output_dir.resolve()
    remote_root = safe_remote_root(args.remote_root)
    variants = [item.strip() for item in args.variants.split(",") if item.strip()]
    if not variants or len(set(variants)) != len(variants):
        raise ValueError("variants must be nonempty and unique")
    unknown = sorted(set(variants) - set(VAN_OVERSTRAETEN))
    if unknown:
        raise ValueError(f"unknown variants: {', '.join(unknown)}")
    bundle = output / "bundle"
    raw = output / "raw"
    if bundle.exists():
        shutil.rmtree(bundle)
    if raw.exists():
        shutil.rmtree(raw)
    bundle.mkdir(parents=True)
    raw.mkdir(parents=True)

    shutil.copy2(source / "n1_msh.tdr", bundle / "n1_msh.tdr")
    # Sentaurus appends ``_des`` to a Load FilePrefix before opening the TDR.
    shutil.copy2(source / "n4_des.tdr", bundle / "source_n4_des_des.tdr")
    base_parameter = (source / "pp4_des.par").read_text(encoding="utf-8")
    for variant in variants:
        override = VAN_OVERSTRAETEN[variant]
        (bundle / f"mean_control_{variant}.cmd").write_text(
            deck(variant, args.drive), encoding="utf-8"
        )
        parameter = base_parameter.rstrip() + "\n"
        if override is not None:
            parameter += "\n" + override.strip() + "\n"
        (bundle / f"mean_control_{variant}.par").write_text(
            parameter, encoding="utf-8"
        )

    ssh_common = ssh_host_options(args.host_name)
    run([str(args.ssh_bin), *ssh_common, args.ssh_target, f"mkdir -p {remote_root}"])
    run(
        [
            str(args.scp_bin),
            *ssh_common,
            *(str(path) for path in sorted(bundle.iterdir())),
            f"{args.ssh_target}:{remote_root}/",
        ]
    )
    remote_commands = " && ".join(
        f"sdevice mean_control_{variant}.cmd"
        for variant in variants
    )
    archive_files = " ".join(
        f"mean_control_{variant}.*" for variant in variants
    )
    run(
        [
            str(args.ssh_bin),
            *ssh_common,
            args.ssh_target,
            (
                f"cd {remote_root} && {remote_commands} && "
                f"tar -czf mean_controls.tgz {archive_files}"
            ),
        ]
    )
    run(
        [
            str(args.scp_bin),
            *ssh_common,
            f"{args.ssh_target}:{remote_root}/mean_controls.tgz",
            str(raw),
        ]
    )
    manifest = {
        "source_state": str((source / "n4_des.tdr").resolve()),
        "variants": variants,
        "drive": args.drive,
        "remote_root": remote_root,
        "archive": str((raw / "mean_controls.tgz").resolve()),
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(raw / "mean_controls.tgz")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
