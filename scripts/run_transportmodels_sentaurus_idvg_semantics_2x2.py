#!/usr/bin/env python3
"""Run the TransportModels DG Id-Vg DirectQC/NoFermi 2x2 oracle on the VM."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
import tarfile
from pathlib import Path
from typing import Any, Sequence


REPO_ROOT = Path(__file__).resolve().parents[1]
SOURCE = REPO_ROOT / "build-release/reference_tcad/transportmodels_sentaurus2022/run02/full_raw"
OUTPUT_ROOT = (
    REPO_ROOT
    / "build-release/reference_tcad/transportmodels_sentaurus2022/sentaurus_vm_runs"
    / "idvg_semantics_2x2_20260821"
)
REPORT_JSON = REPO_ROOT / "docs/validation/transportmodels_sentaurus_idvg_semantics_2x2_2026-08-21.json"
REPORT_MD = REPO_ROOT / "docs/validation/transportmodels_sentaurus_idvg_semantics_2x2_2026-08-21.md"
REMOTE_ROOT = "~/sentaurus_runs/vela_oracle/transportmodels_idvg_semantics_2x2_20260821"
KEY_BIASES = (-0.20, -0.04, 0.12, 0.28, 1.00)
VARIANTS = (
    {"name": "default_default", "direct_qc": False, "no_fermi": False},
    {"name": "directqc_default", "direct_qc": True, "no_fermi": False},
    {"name": "default_nofermi", "direct_qc": False, "no_fermi": True},
    {"name": "directqc_nofermi", "direct_qc": True, "no_fermi": True},
)


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
    base_deck = (SOURCE / "pp7_des.cmd").read_text(encoding="utf-8")
    required = (SOURCE / "n1_msh.tdr", SOURCE / "pp7_des.par")
    for source in required:
        if not source.exists():
            raise FileNotFoundError(source)
    for variant in VARIANTS:
        variant_dir = bundle / str(variant["name"])
        variant_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(required[0], variant_dir / required[0].name)
        shutil.copy2(required[1], variant_dir / required[1].name)
        deck = base_deck
        # The official 20-interval output contains 0.92 V and 1.08 V but not
        # the planned 1.00 V oracle point.  Forty intervals preserve the same
        # -1.0 -> 2.2 V continuation path while adding an exact 1.00 V sample.
        current_plot = "CurrentPlot(Time=(Range=(0 1) Intervals=20))"
        if current_plot not in deck:
            raise RuntimeError("Expected CurrentPlot interval specification not found")
        deck = deck.replace(
            current_plot,
            "CurrentPlot(Time=(Range=(0 1) Intervals=40))",
            1,
        )
        if variant["no_fermi"]:
            old = "EffectiveIntrinsicDensity( OldSlotboom )"
            new = "EffectiveIntrinsicDensity( OldSlotboom NoFermi )"
            if old not in deck:
                raise RuntimeError("Expected EffectiveIntrinsicDensity line not found")
            deck = deck.replace(old, new, 1)
        if variant["direct_qc"]:
            marker = "Math {\n"
            if marker not in deck:
                raise RuntimeError("Expected Math block not found")
            deck = deck.replace(marker, marker + "   DirectQuantumCorrection\n", 1)
        (variant_dir / "idvg_2x2_des.cmd").write_text(deck, encoding="utf-8")
    return bundle


def parse_plt(path: Path) -> tuple[list[str], list[list[float]]]:
    sys.path.insert(0, str(REPO_ROOT / "scripts"))
    from analyze_sentaurus_bvmethods import read_plt

    return read_plt(path)


def extract_curve(path: Path, output: Path) -> list[dict[str, float]]:
    datasets, values = parse_plt(path)
    gate_index = datasets.index("gate OuterVoltage")
    drain_index = datasets.index("drain TotalCurrent")
    rows = [
        {"bias_V": float(row[gate_index]), "current_A_per_um": abs(float(row[drain_index]))}
        for row in values
    ]
    # The current file can include solver-internal rows.  Retain the last value
    # for each requested CurrentPlot bias.
    unique = {round(row["bias_V"], 10): row for row in rows}
    rows = [unique[key] for key in sorted(unique)]
    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["bias_V", "current_A_per_um"])
        writer.writeheader()
        writer.writerows(rows)
    return rows


def key_rows(rows: list[dict[str, float]]) -> dict[float, float]:
    result: dict[float, float] = {}
    for target in KEY_BIASES:
        closest = min(rows, key=lambda row: abs(row["bias_V"] - target))
        if abs(closest["bias_V"] - target) > 1.0e-6:
            raise RuntimeError(f"Missing requested Vg={target:g} V point")
        result[target] = closest["current_A_per_um"]
    return result


def make_plot(curves: dict[str, list[dict[str, float]]]) -> tuple[Path, Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    labels = {
        "default_default": "Default DG + default Fermi-BGN",
        "directqc_default": "DirectQC + default Fermi-BGN",
        "default_nofermi": "Default DG + NoFermi",
        "directqc_nofermi": "DirectQC + NoFermi",
    }
    markers = ("o-", "s--", "^-.", "d:")
    fig, axes = plt.subplots(1, 2, figsize=(12.2, 5.1))
    default = key_rows(curves["default_default"])
    for marker, variant in zip(markers, VARIANTS):
        name = str(variant["name"])
        selected = key_rows(curves[name])
        axes[0].semilogy(KEY_BIASES, [selected[bias] for bias in KEY_BIASES], marker, label=labels[name])
        axes[1].plot(
            KEY_BIASES,
            [math.log10(max(selected[bias], 1.0e-30) / max(default[bias], 1.0e-30)) for bias in KEY_BIASES],
            marker,
            label=labels[name],
        )
    axes[0].set_xlabel("Gate voltage Vg (V)")
    axes[0].set_ylabel("Drain current Id (A/µm)")
    axes[0].set_title("Five-point Sentaurus DG oracle")
    axes[1].set_xlabel("Gate voltage Vg (V)")
    axes[1].set_ylabel("log10(Id / default Id)")
    axes[1].set_title("Semantic-model contribution")
    for axis in axes:
        axis.grid(True, which="both", alpha=0.25)
        axis.legend(fontsize=8)
    fig.tight_layout()
    png = OUTPUT_ROOT / "sentaurus_idvg_semantics_2x2.png"
    svg = OUTPUT_ROOT / "sentaurus_idvg_semantics_2x2.svg"
    fig.savefig(png, dpi=180)
    fig.savefig(svg)
    plt.close(fig)
    return png, svg


def analyze(raw: Path, banner: str | None) -> dict[str, Any]:
    curves: dict[str, list[dict[str, float]]] = {}
    artifacts: dict[str, Any] = {}
    for variant in VARIANTS:
        name = str(variant["name"])
        variant_raw = raw / "bundle" / name
        plot = variant_raw / "IdVgs_n7_des.plt"
        log = variant_raw / "n7_des.log"
        if not plot.exists() or not log.exists():
            raise FileNotFoundError(f"Missing Sentaurus output for {name}")
        curve_csv = OUTPUT_ROOT / f"{name}_full_curve.csv"
        curves[name] = extract_curve(plot, curve_csv)
        artifacts[name] = {
            "deck": str((variant_raw / "idvg_2x2_des.cmd").resolve()),
            "plot": str(plot.resolve()),
            "log": str(log.resolve()),
            "curve_csv": str(curve_csv.resolve()),
            "curve_sha256": sha256(curve_csv),
            "plot_sha256": sha256(plot),
            "log_sha256": sha256(log),
        }
    keys = {name: key_rows(rows) for name, rows in curves.items()}
    baseline = keys["default_default"]
    table = []
    for bias in KEY_BIASES:
        row: dict[str, Any] = {"gate_bias_V": bias}
        for variant in VARIANTS:
            name = str(variant["name"])
            current = keys[name][bias]
            row[name + "_A_per_um"] = current
            row[name + "_delta_vs_default_dex"] = math.log10(
                max(current, 1.0e-30) / max(baseline[bias], 1.0e-30)
            )
        row["directqc_effect_default_bgn_dex"] = math.log10(
            keys["directqc_default"][bias] / baseline[bias]
        )
        row["nofermi_effect_default_dg_dex"] = math.log10(
            keys["default_nofermi"][bias] / baseline[bias]
        )
        row["interaction_dex"] = (
            math.log10(keys["directqc_nofermi"][bias] / baseline[bias])
            - row["directqc_effect_default_bgn_dex"]
            - row["nofermi_effect_default_dg_dex"]
        )
        table.append(row)
    png, svg = make_plot(curves)
    return {
        "schema": "vela.transportmodels.sentaurus_idvg_semantics_2x2.v1",
        "as_of": "2026-08-21",
        "status": "complete",
        "sentaurus_banner": banner,
        "design": {
            "direct_quantum_correction": [False, True],
            "fermi_bgn_correction": ["default", "NoFermi"],
            "gate_biases_V": list(KEY_BIASES),
            "fixed_drain_bias_V": 1.1,
            "current_plot_intervals": 40,
        },
        "five_point_comparison": table,
        "artifacts": {**artifacts, "png": str(png.resolve()), "svg": str(svg.resolve())},
    }


def markdown(report: dict[str, Any]) -> str:
    lines = [
        "# Sentaurus TransportModels DG Id-Vg semantic 2x2 oracle",
        "",
        f"Status: **{report['status']}**. Release: `{report['sentaurus_banner']}`.",
        "",
        "The four decks differ only by `Math { DirectQuantumCorrection }` and "
        "`EffectiveIntrinsicDensity(OldSlotboom NoFermi)`.",
        "",
        "| Vg (V) | DirectQC effect, default BGN (dex) | NoFermi effect, default DG (dex) | Interaction (dex) |",
        "|---:|---:|---:|---:|",
    ]
    for row in report["five_point_comparison"]:
        lines.append(
            f"| {row['gate_bias_V']:.2f} | {row['directqc_effect_default_bgn_dex']:.6g} | "
            f"{row['nofermi_effect_default_dg_dex']:.6g} | {row['interaction_dex']:.6g} |"
        )
    lines.extend(["", f"Figure: `{report['artifacts']['png']}`", ""])
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--live", action="store_true")
    parser.add_argument("--report-only", action="store_true")
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--ssh-target", default="sentaurus")
    parser.add_argument("--ssh-bin", default=executable("ssh"))
    parser.add_argument("--scp-bin", default=executable("scp"))
    args = parser.parse_args()
    if args.check:
        report = json.loads(REPORT_JSON.read_text(encoding="utf-8"))
        for variant in VARIANTS:
            artifact = report["artifacts"][str(variant["name"])]
            assert sha256(Path(artifact["curve_csv"])) == artifact["curve_sha256"]
            assert sha256(Path(artifact["plot"])) == artifact["plot_sha256"]
            assert sha256(Path(artifact["log"])) == artifact["log_sha256"]
        print("TransportModels Sentaurus Id-Vg semantic 2x2 check: PASS")
        return 0

    bundle = prepare_bundle()
    raw = OUTPUT_ROOT / "raw"
    banner = None
    if args.live and not args.report_only:
        banner = run(
            [args.ssh_bin, args.ssh_target, "sdevice -h 2>&1 | sed -n '1,5p'"],
            capture=True,
        ).strip()
        if "T-2022.03-SP2" not in banner:
            raise RuntimeError(f"Unexpected Sentaurus release:\n{banner}")
        run([args.ssh_bin, args.ssh_target, f"mkdir -p {REMOTE_ROOT}"])
        run([args.scp_bin, "-r", str(bundle), f"{args.ssh_target}:{REMOTE_ROOT}/"])
        for index, variant in enumerate(VARIANTS, start=1):
            name = str(variant["name"])
            print(f"[{index}/4] running {name}", flush=True)
            command = (
                f"cd {REMOTE_ROOT}/bundle/{name} && "
                "sdevice idvg_2x2_des.cmd > run_sdevice.out 2>&1"
            )
            run([args.ssh_bin, args.ssh_target, command])
        archive_name = "idvg_semantics_2x2_results.tgz"
        run(
            [
                args.ssh_bin,
                args.ssh_target,
                f"cd {REMOTE_ROOT} && tar -czf {archive_name} bundle",
            ]
        )
        raw.mkdir(parents=True, exist_ok=True)
        archive = raw / archive_name
        run([args.scp_bin, f"{args.ssh_target}:{REMOTE_ROOT}/{archive_name}", str(archive)])
        with tarfile.open(archive, "r:gz") as stream:
            stream.extractall(raw, filter="data")
        (OUTPUT_ROOT / "sentaurus_banner.txt").write_text(banner + "\n", encoding="utf-8")
    elif (OUTPUT_ROOT / "sentaurus_banner.txt").exists():
        banner = (OUTPUT_ROOT / "sentaurus_banner.txt").read_text(encoding="utf-8").strip()
    else:
        print("Prepared the 2x2 bundle. Use --live to execute it on the Sentaurus VM.")
        return 0

    report = analyze(raw, banner)
    REPORT_JSON.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    REPORT_MD.write_text(markdown(report), encoding="utf-8")
    print(json.dumps({"status": report["status"], "five_point_comparison": report["five_point_comparison"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
