#!/usr/bin/env python3
"""Run the five-point Vela DD/DG Fermi-BGN correction A/B matrix."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
from pathlib import Path
from typing import Any


REPO_ROOT = Path(__file__).resolve().parents[1]
BASELINE = REPO_ROOT / "build-release/reference_tcad/transportmodels_sentaurus2022/vela_baseline"
OUTPUT_ROOT = BASELINE / "vela_fermi_bgn_ab_2026-08-21"
PHASE7_SCRIPT = REPO_ROOT / "scripts/run_transportmodels_dg_phase7_regression.py"
SENTAURUS_ORACLE = REPO_ROOT / "docs/validation/transportmodels_sentaurus_idvg_semantics_2x2_2026-08-21.json"
REPORT_JSON = REPO_ROOT / "docs/validation/transportmodels_vela_fermi_bgn_ab_2026-08-21.json"
REPORT_MD = REPO_ROOT / "docs/validation/transportmodels_vela_fermi_bgn_ab_2026-08-21.md"
POINTS = (-0.20, -0.04, 0.12, 0.28, 1.00)


def load_phase7_module():
    spec = importlib.util.spec_from_file_location("transportmodels_phase7", PHASE7_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load {PHASE7_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def sha256(path: Path) -> str:
    import hashlib

    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def load_curve(path: Path) -> dict[float, dict[str, Any]]:
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    return {
        round(float(row["bias_V"]), 10): {
            "current_A_per_um": abs(float(row["current_total_A_per_um"])),
            "converged": row["converged"] == "1",
            "iterations": int(row["iterations"]),
        }
        for row in rows
    }


def make_plot(results: dict[str, dict[float, dict[str, Any]]], sentaurus: dict[float, dict[str, float]]) -> tuple[Path, Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(12.0, 5.1))
    style = {
        "dd_off": ("o--", "Vela DD, Fermi-BGN off"),
        "dd_on": ("o-", "Vela DD, Fermi-BGN on"),
        "dg_off": ("s--", "Vela DG, Fermi-BGN off"),
        "dg_on": ("s-", "Vela DG, Fermi-BGN on"),
    }
    for name, (marker, label) in style.items():
        axes[0].semilogy(POINTS, [results[name][round(bias, 10)]["current_A_per_um"] for bias in POINTS], marker, label=label)
    axes[0].semilogy(POINTS, [sentaurus[bias]["default"] for bias in POINTS], "k^-", label="Sentaurus DG default")
    axes[0].semilogy(POINTS, [sentaurus[bias]["nofermi"] for bias in POINTS], "k^--", label="Sentaurus DG NoFermi")
    for branch, color in (("dd", "tab:blue"), ("dg", "tab:orange")):
        off = results[branch + "_off"]
        on = results[branch + "_on"]
        axes[1].plot(
            POINTS,
            [math.log10(on[round(bias, 10)]["current_A_per_um"] / off[round(bias, 10)]["current_A_per_um"]) for bias in POINTS],
            "o-",
            color=color,
            label=f"Vela {branch.upper()}: on/off",
        )
    axes[1].plot(
        POINTS,
        [math.log10(sentaurus[bias]["default"] / sentaurus[bias]["nofermi"]) for bias in POINTS],
        "k^--",
        label="Sentaurus DG: default/NoFermi",
    )
    axes[0].set_xlabel("Gate voltage Vg (V)")
    axes[0].set_ylabel("Drain current Id (A/µm)")
    axes[0].set_title("Five-point current comparison")
    axes[1].set_xlabel("Gate voltage Vg (V)")
    axes[1].set_ylabel("log10(Id on / Id off)")
    axes[1].set_title("Fermi-BGN correction effect")
    for axis in axes:
        axis.grid(True, which="both", alpha=0.25)
        axis.legend(fontsize=8)
    fig.tight_layout()
    png = OUTPUT_ROOT / "vela_fermi_bgn_ab.png"
    svg = OUTPUT_ROOT / "vela_fermi_bgn_ab.svg"
    fig.savefig(png, dpi=180)
    fig.savefig(svg)
    plt.close(fig)
    return png, svg


def markdown(report: dict[str, Any]) -> str:
    lines = [
        "# TransportModels Vela Fermi-BGN five-point A/B",
        "",
        f"Status: **{report['status']}**.",
        "",
        "| Vg (V) | Vela DD on/off (dex) | Vela DG on/off (dex) | Sentaurus default/NoFermi (dex) | DG-on vs Sentaurus default (dex) |",
        "|---:|---:|---:|---:|---:|",
    ]
    for row in report["comparison"]:
        lines.append(
            f"| {row['gate_bias_V']:.2f} | {row['vela_dd_on_off_dex']:.6g} | "
            f"{row['vela_dg_on_off_dex']:.6g} | {row['sentaurus_default_nofermi_dex']:.6g} | "
            f"{row['vela_dg_on_vs_sentaurus_default_dex']:.6g} |"
        )
    lines.extend(["", f"Figure: `{report['artifacts']['png']}`", ""])
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--report-only", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.check:
        report = json.loads(REPORT_JSON.read_text(encoding="utf-8"))
        for artifact in report["variants"].values():
            assert sha256(Path(artifact["curve_csv"])) == artifact["curve_sha256"]
            assert artifact["completed_points"] == len(POINTS)
            assert artifact["all_converged"]
        print("TransportModels Vela Fermi-BGN A/B check: PASS")
        return 0

    phase7 = load_phase7_module()
    phase7.OUTPUT_ROOT = OUTPUT_ROOT
    variants = (
        {
            "name": "dd_off", "branch": "dd", "correction": False,
            "initial_state": BASELINE / "dd_phase7_shared_baseline_2026-08-21/idvg/state_bias_m0p200000.csv",
        },
        {
            "name": "dd_on", "branch": "dd", "correction": True,
            "initial_state": BASELINE / "dd_phase7_shared_baseline_2026-08-21/idvg/state_bias_m0p200000.csv",
        },
        {
            "name": "dg_off", "branch": "dg", "correction": False,
            "initial_state": BASELINE / "dg_phase7_regression_2026-08-21/idvg/state_bias_m0p200000.csv",
        },
        {
            "name": "dg_on", "branch": "dg", "correction": True,
            "initial_state": BASELINE / "dg_phase7_regression_2026-08-21/idvg/state_bias_m0p200000.csv",
        },
    )
    curves = [
        {
            **variant,
            "label": f"{variant['branch'].upper()} Fermi-BGN {'on' if variant['correction'] else 'off'}",
            "config": phase7.GENERATED / "simulation_dg_idvg.json",
            "reference": phase7.REFERENCE / "transportmodels_sentaurus2022_dg_idvg_reference.csv",
            "contact": "gate", "current_contact": "drain", "points": list(POINTS),
        }
        for variant in variants
    ]
    original_make_config = phase7.make_config

    def make_ab_config(curve: dict[str, Any]) -> Path:
        path = original_make_config(curve)
        config = json.loads(path.read_text(encoding="utf-8"))
        config["_comment"] = curve["label"] + " five-point controlled A/B"
        config["solver"]["electron_quantum_potential"]["enabled"] = curve["branch"] == "dg"
        config["solver"]["bandgap_narrowing"] = {
            "model": "old_slotboom",
            "fermi_statistics_correction": bool(curve["correction"]),
        }
        path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
        return path

    phase7.make_config = make_ab_config
    execution = []
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    if not args.report_only:
        for curve in curves:
            execution.append(phase7.execute_curve(curve))

    results: dict[str, dict[float, dict[str, Any]]] = {}
    artifacts: dict[str, Any] = {}
    for curve in curves:
        candidate = OUTPUT_ROOT / str(curve["name"]) / "curve_combined.csv"
        results[str(curve["name"])] = load_curve(candidate)
        artifacts[str(curve["name"])] = {
            "branch": curve["branch"],
            "fermi_statistics_correction": curve["correction"],
            "config": str((OUTPUT_ROOT / str(curve["name"]) / "config.json").resolve()),
            "curve_csv": str(candidate.resolve()),
            "curve_sha256": sha256(candidate),
            "completed_points": len(results[str(curve["name"])]),
            "all_converged": all(row["converged"] for row in results[str(curve["name"])].values()),
        }
    oracle = json.loads(SENTAURUS_ORACLE.read_text(encoding="utf-8"))
    sentaurus = {
        float(row["gate_bias_V"]): {
            "default": float(row["default_default_A_per_um"]),
            "nofermi": float(row["default_nofermi_A_per_um"]),
        }
        for row in oracle["five_point_comparison"]
    }
    comparison = []
    for bias in POINTS:
        key = round(bias, 10)
        dd_off = results["dd_off"][key]["current_A_per_um"]
        dd_on = results["dd_on"][key]["current_A_per_um"]
        dg_off = results["dg_off"][key]["current_A_per_um"]
        dg_on = results["dg_on"][key]["current_A_per_um"]
        comparison.append(
            {
                "gate_bias_V": bias,
                "vela_dd_off_A_per_um": dd_off,
                "vela_dd_on_A_per_um": dd_on,
                "vela_dg_off_A_per_um": dg_off,
                "vela_dg_on_A_per_um": dg_on,
                "sentaurus_default_A_per_um": sentaurus[bias]["default"],
                "sentaurus_nofermi_A_per_um": sentaurus[bias]["nofermi"],
                "vela_dd_on_off_dex": math.log10(dd_on / dd_off),
                "vela_dg_on_off_dex": math.log10(dg_on / dg_off),
                "sentaurus_default_nofermi_dex": math.log10(sentaurus[bias]["default"] / sentaurus[bias]["nofermi"]),
                "vela_dg_on_vs_sentaurus_default_dex": math.log10(dg_on / sentaurus[bias]["default"]),
                "vela_dg_off_vs_sentaurus_nofermi_dex": math.log10(dg_off / sentaurus[bias]["nofermi"]),
            }
        )
    complete = all(
        artifact["completed_points"] == len(POINTS) and artifact["all_converged"]
        for artifact in artifacts.values()
    )
    png, svg = make_plot(results, sentaurus)
    report = {
        "schema": "vela.transportmodels.vela_fermi_bgn_ab.v1",
        "as_of": "2026-08-21",
        "status": "complete" if complete else "partial",
        "execution": execution,
        "design": {"gate_biases_V": list(POINTS), "fixed_drain_bias_V": 1.1},
        "variants": artifacts,
        "comparison": comparison,
        "sentaurus_oracle": str(SENTAURUS_ORACLE.resolve()),
        "artifacts": {"png": str(png.resolve()), "svg": str(svg.resolve())},
    }
    REPORT_JSON.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    REPORT_MD.write_text(markdown(report), encoding="utf-8")
    print(json.dumps({"status": report["status"], "comparison": comparison}, indent=2))
    return 0 if complete else 1


if __name__ == "__main__":
    raise SystemExit(main())
