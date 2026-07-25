#!/usr/bin/env python3
"""Run the Minimal6 element-avalanche runtime probe at all 20 biases."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import subprocess
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
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
    parser.add_argument(
        "--output-root",
        type=Path,
        default=Path(
            "build-release/"
            "pn2d-minimal6-element-avalanche-all20-runtime-20260725"
        ),
    )
    parser.add_argument(
        "--remote-root",
        default="/home/tcad/codex_pn2d_minimal6_element_avalanche_all20_20260725",
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


def make_tcl(template: str) -> str:
    targets = " ".join(f"-{value}" for value in range(1, 21))
    expected = "foreach candidate {-1 -10 -20} {"
    replacement = f"foreach candidate {{{targets}}} {{"
    if expected not in template:
        raise RuntimeError("Tcl template target list was not found")
    return template.replace(expected, replacement, 1)


def make_solve_block() -> str:
    lines = [
        "Solve {",
        "  Coupled(Iterations=100) { Poisson }",
        "  Coupled(Iterations=100) { Poisson Electron Hole }",
    ]
    for magnitude in range(1, 21):
        lines.extend(
            [
                "  Quasistationary(",
                "    InitialStep=1e-3 MinStep=1e-10 MaxStep=0.05",
                "    Increment=1.2 Decrement=2.0",
                f'    Goal {{ Name="Anode" Voltage=-{magnitude} }}',
                "  ) { Coupled { Poisson Electron Hole } }",
            ]
        )
    lines.append("}")
    return "\n".join(lines)


def make_deck(template: str) -> str:
    renamed = template.replace(
        "runtime_element_avalanche_probe_default",
        "runtime_element_avalanche_probe_all20",
    )
    result, count = re.subn(
        r"Solve\s*\{.*\}\s*$",
        make_solve_block(),
        renamed,
        count=1,
        flags=re.S,
    )
    if count != 1:
        raise RuntimeError("SDevice Solve block was not replaced exactly once")
    return result.rstrip() + "\n"


def main() -> int:
    args = parse_args()
    transport_root = args.transport_root.resolve()
    template_root = args.template_root.resolve()
    output_root = args.output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)

    tcl_template = (
        template_root / "runtime_element_avalanche_probe.tcl"
    ).read_text(encoding="ascii")
    deck_template = (
        template_root / "runtime_element_avalanche_probe_default.cmd"
    ).read_text(encoding="ascii")
    tcl_text = make_tcl(tcl_template)
    deck_text = make_deck(deck_template)

    manifest: dict[str, object] = {
        "schema_version": 1,
        "status": "running",
        "ssh_target": args.ssh_target,
        "remote_root": args.remote_root,
        "topologies": {},
    }
    manifest_path = output_root / "manifest.json"

    for topology in ("mirror", "sketch"):
        local = output_root / topology
        bundle = local / "bundle"
        fetched = local / "fetched"
        bundle.mkdir(parents=True, exist_ok=True)
        fetched.mkdir(parents=True, exist_ok=True)
        source = transport_root / topology / "m1V"
        shutil.copy2(source / "pn2d_minimal6.tdr", bundle)
        shutil.copy2(source / "models.par", bundle)
        (bundle / "runtime_element_avalanche_probe.tcl").write_text(
            tcl_text, encoding="ascii", newline="\n"
        )
        deck_name = "runtime_element_avalanche_probe_all20.cmd"
        (bundle / deck_name).write_text(
            deck_text, encoding="ascii", newline="\n"
        )

        remote = f"{args.remote_root}/{topology}"
        run([str(args.ssh_bin), args.ssh_target, f"mkdir -p {remote}"])
        for name in (
            "pn2d_minimal6.tdr",
            "models.par",
            "runtime_element_avalanche_probe.tcl",
            deck_name,
        ):
            run(
                [
                    str(args.scp_bin),
                    str(bundle / name),
                    f"{args.ssh_target}:{remote}/",
                ]
            )
        run(
            [
                str(args.ssh_bin),
                args.ssh_target,
                (
                    f"cd {remote} && "
                    f"sdevice {deck_name} > run_all20.out 2>&1"
                ),
            ]
        )
        fetched_names = (
            "run_all20.out",
            "runtime_element_avalanche_probe_all20.plt",
            "runtime_element_avalanche_probe_all20_des.log",
        )
        for name in fetched_names:
            run(
                [
                    str(args.scp_bin),
                    f"{args.ssh_target}:{remote}/{name}",
                    str(fetched / name),
                ]
            )
        text = (fetched / "run_all20.out").read_text(
            encoding="ascii", errors="strict"
        )
        observed = sorted(
            {
                int(value)
                for value in re.findall(
                    r"AVAL_PROBE_BEGIN bias_V=(-?\d+)", text
                )
            }
        )
        expected = list(range(-20, 0))
        topology_result = {
            "status": "passed" if observed == expected else "failed",
            "observed_biases_V": observed,
            "output_sha256": {
                name: sha256(fetched / name) for name in fetched_names
            },
        }
        manifest["topologies"][topology] = topology_result
        manifest_path.write_text(
            json.dumps(manifest, indent=2, sort_keys=True) + "\n",
            encoding="ascii",
            newline="\n",
        )
        if observed != expected:
            raise RuntimeError(
                f"{topology}: expected all integer biases -20..-1, got "
                f"{observed}"
            )

    manifest["status"] = "passed"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="ascii",
        newline="\n",
    )
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
