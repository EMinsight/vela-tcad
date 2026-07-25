#!/usr/bin/env python3
"""Run explicit Minimal6 Sentaurus avalanche driving-force controls."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import subprocess
from pathlib import Path, PurePosixPath


VARIANTS = {
    "implicit_default": {
        "avalanche": "Avalanche(VanOverstraeten)",
        "aval_dens_grad_qf": False,
    },
    "explicit_grad_qf": {
        "avalanche": "Avalanche(VanOverstraeten GradQuasiFermi)",
        "aval_dens_grad_qf": False,
    },
    "explicit_electric_field": {
        "avalanche": "Avalanche(VanOverstraeten ElectricField)",
        "aval_dens_grad_qf": False,
    },
    "grad_qf_aval_dens_grad_qf": {
        "avalanche": "Avalanche(VanOverstraeten GradQuasiFermi)",
        "aval_dens_grad_qf": True,
    },
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ssh-target", default="sentaurus")
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
    parser.add_argument(
        "--transport-root",
        type=Path,
        default=Path(
            "build-release/pn2d-minimal6-transport-elements-20260723-b/"
            "codex_pn2d_minimal6_transport_elements_20260723_b"
        ),
    )
    parser.add_argument(
        "--template-root",
        type=Path,
        default=Path(
            "build-release/"
            "pn2d-minimal6-element-avalanche-replay-20260725"
        ),
    )
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--remote-root", required=True)
    parser.add_argument(
        "--biases",
        nargs="+",
        type=int,
        default=[-1, -10, -20],
    )
    return parser.parse_args()


def run(command: list[str]) -> None:
    subprocess.run(command, check=True)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_biases(biases: tuple[int, ...]) -> None:
    if not biases:
        raise ValueError("at least one bias is required")
    if len(set(biases)) != len(biases):
        raise ValueError("biases must be unique")
    if any(bias >= 0 for bias in biases):
        raise ValueError("all avalanche control biases must be negative")
    if tuple(sorted(biases, reverse=True)) != biases:
        raise ValueError("biases must be ordered from low to high magnitude")


def validate_remote_root(value: str) -> str:
    if re.fullmatch(r"/[A-Za-z0-9._/-]+", value) is None:
        raise ValueError("remote root must be a safe absolute POSIX path")
    path = PurePosixPath(value)
    if (
        not path.is_absolute()
        or str(path) != value
        or any(part in {".", ".."} for part in path.parts)
    ):
        raise ValueError("remote root must be a normalized absolute POSIX path")
    return value


def make_solve_block(biases: tuple[int, ...]) -> str:
    validate_biases(biases)
    lines = [
        "Solve {",
        "  Coupled(Iterations=100) { Poisson }",
        "  Coupled(Iterations=100) { Poisson Electron Hole }",
    ]
    for index, bias in enumerate(biases):
        lines.extend(
            [
                "  Quasistationary(",
                (
                    "    InitialStep=1e-4 MinStep=1e-10 MaxStep=0.05"
                    if index == 0
                    else
                    "    InitialStep=1e-3 MinStep=1e-10 MaxStep=0.05"
                ),
                "    Increment=1.2 Decrement=2.0",
                f'    Goal {{ Name="Anode" Voltage={bias} }}',
                "  ) { Coupled { Poisson Electron Hole } }",
            ]
        )
    lines.append("}")
    return "\n".join(lines)


def make_variant_deck(
    template: str,
    variant: str,
    biases: tuple[int, ...],
) -> str:
    if variant not in VARIANTS:
        raise ValueError(f"unknown avalanche drive variant: {variant}")
    validate_biases(biases)
    default_output = "runtime_element_avalanche_probe_default"
    if default_output not in template:
        raise ValueError("default output stem not found in template")
    result = template.replace(
        default_output,
        f"runtime_element_avalanche_probe_{variant}",
    )
    default_avalanche = "Avalanche(VanOverstraeten)"
    if result.count(default_avalanche) != 1:
        raise ValueError("default Avalanche selector must occur exactly once")
    result = result.replace(
        default_avalanche,
        str(VARIANTS[variant]["avalanche"]),
        1,
    )
    if VARIANTS[variant]["aval_dens_grad_qf"]:
        result, math_count = re.subn(
            r"Math\s*\{",
            "Math {\n  AvalDensGradQF",
            result,
            count=1,
        )
        if math_count != 1:
            raise ValueError("Math block was not found exactly once")
    result, solve_count = re.subn(
        r"Solve\s*\{.*\}\s*$",
        make_solve_block(biases),
        result,
        count=1,
        flags=re.S,
    )
    if solve_count != 1:
        raise ValueError("Solve block was not replaced exactly once")
    return result.rstrip() + "\n"


def make_tcl(template: str, biases: tuple[int, ...]) -> str:
    validate_biases(biases)
    expected = "foreach candidate {-1 -10 -20} {"
    replacement = (
        "foreach candidate {"
        + " ".join(str(bias) for bias in biases)
        + "} {"
    )
    if expected not in template:
        raise ValueError("Tcl target list was not found")
    return template.replace(expected, replacement, 1)


def sentaurus_release(ssh_bin: Path, ssh_target: str) -> str:
    command = (
        'resolved=$(readlink -f "$(command -v sdevice)") && '
        'printf "path=%s\\n" "$resolved" && sdevice --version 2>&1'
    )
    completed = subprocess.run(
        [str(ssh_bin), ssh_target, command],
        check=False,
        capture_output=True,
        text=True,
    )
    output = completed.stdout + completed.stderr
    match = re.search(r"Version\s+([^\s*]+)", output)
    if match is None:
        raise RuntimeError("Sentaurus release was not found")
    return match.group(1)


def write_manifest(path: Path, manifest: dict[str, object]) -> None:
    path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="ascii",
        newline="\n",
    )


def main() -> int:
    args = parse_args()
    biases = tuple(args.biases)
    validate_biases(biases)
    remote_root = validate_remote_root(args.remote_root)
    transport_root = args.transport_root.resolve()
    template_root = args.template_root.resolve()
    output_root = args.output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    manifest_path = output_root / "manifest.json"

    tcl_template = (
        template_root / "runtime_element_avalanche_probe.tcl"
    ).read_text(encoding="ascii")
    deck_template = (
        template_root / "runtime_element_avalanche_probe_default.cmd"
    ).read_text(encoding="ascii")
    tcl_text = make_tcl(tcl_template, biases)

    manifest: dict[str, object] = {
        "schema_version": 1,
        "status": "running",
        "experiment": "pn2d_minimal6_sentaurus_avalanche_drive_controls",
        "sentaurus_release": sentaurus_release(
            args.ssh_bin,
            args.ssh_target,
        ),
        "ssh_target": args.ssh_target,
        "remote_root": remote_root,
        "biases_V": list(biases),
        "variants": list(VARIANTS),
        "topologies": {},
    }
    write_manifest(manifest_path, manifest)

    for topology in ("mirror", "sketch"):
        topology_results: dict[str, object] = {}
        source = transport_root / topology / "m1V"
        for variant in VARIANTS:
            local = output_root / topology / variant
            bundle = local / "bundle"
            fetched = local / "fetched"
            bundle.mkdir(parents=True, exist_ok=True)
            fetched.mkdir(parents=True, exist_ok=True)
            shutil.copy2(source / "pn2d_minimal6.tdr", bundle)
            shutil.copy2(source / "models.par", bundle)
            tcl_name = "runtime_element_avalanche_probe.tcl"
            (bundle / tcl_name).write_text(
                tcl_text,
                encoding="ascii",
                newline="\n",
            )
            deck_name = (
                f"runtime_element_avalanche_probe_{variant}.cmd"
            )
            deck_text = make_variant_deck(
                deck_template,
                variant,
                biases,
            )
            (bundle / deck_name).write_text(
                deck_text,
                encoding="ascii",
                newline="\n",
            )

            remote = f"{remote_root}/{topology}/{variant}"
            run(
                [
                    str(args.ssh_bin),
                    args.ssh_target,
                    f"mkdir -p {remote}",
                ]
            )
            for name in (
                "pn2d_minimal6.tdr",
                "models.par",
                tcl_name,
                deck_name,
            ):
                run(
                    [
                        str(args.scp_bin),
                        str(bundle / name),
                        f"{args.ssh_target}:{remote}/",
                    ]
                )
            run_name = f"run_{variant}.out"
            run(
                [
                    str(args.ssh_bin),
                    args.ssh_target,
                    (
                        f"cd {remote} && sdevice {deck_name} "
                        f"> {run_name} 2>&1"
                    ),
                ]
            )
            stem = f"runtime_element_avalanche_probe_{variant}"
            fetched_names = (
                run_name,
                f"{stem}.plt",
                f"{stem}_des.log",
            )
            for name in fetched_names:
                run(
                    [
                        str(args.scp_bin),
                        f"{args.ssh_target}:{remote}/{name}",
                        str(fetched / name),
                    ]
                )
            text = (fetched / run_name).read_text(
                encoding="ascii",
                errors="strict",
            )
            observed = sorted(
                {
                    int(value)
                    for value in re.findall(
                        r"AVAL_PROBE_BEGIN bias_V=(-?\d+)",
                        text,
                    )
                },
                reverse=True,
            )
            status = "passed" if tuple(observed) == biases else "failed"
            topology_results[variant] = {
                "status": status,
                "observed_biases_V": observed,
                "bundle_sha256": {
                    name: sha256(bundle / name)
                    for name in (
                        "pn2d_minimal6.tdr",
                        "models.par",
                        tcl_name,
                        deck_name,
                    )
                },
                "output_sha256": {
                    name: sha256(fetched / name)
                    for name in fetched_names
                },
            }
            manifest["topologies"][topology] = topology_results
            write_manifest(manifest_path, manifest)
            if status != "passed":
                raise RuntimeError(
                    f"{topology}/{variant}: expected {biases}, got "
                    f"{tuple(observed)}"
                )

    manifest["status"] = "passed"
    write_manifest(manifest_path, manifest)
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
