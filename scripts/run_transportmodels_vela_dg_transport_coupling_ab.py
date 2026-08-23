#!/usr/bin/env python3
"""Run the five-point Vela DirectQC/Sentaurus-exponential DG coupling A/B."""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import math
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "build-release/reference_tcad/transportmodels_sentaurus2022/vela_baseline"
OUT = BASE / "dg_transport_coupling_ab_2026-08-21"
PHASE7 = ROOT / "scripts/run_transportmodels_dg_phase7_regression.py"
ORACLE = ROOT / "docs/validation/transportmodels_sentaurus_idvg_semantics_2x2_2026-08-21.json"
REPORT_JSON = ROOT / "docs/validation/transportmodels_dg_transport_coupling_ab_2026-08-21.json"
REPORT_MD = ROOT / "docs/validation/transportmodels_dg_transport_coupling_ab_2026-08-21.md"
POINTS = (-0.20, -0.04, 0.12, 0.28, 1.00)


def phase7_module():
    spec = importlib.util.spec_from_file_location("transportmodels_phase7", PHASE7)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Unable to load {PHASE7}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def curve(path: Path) -> dict[float, dict[str, Any]]:
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


def bias_slug(value: float) -> str:
    sign = "m" if value < 0.0 else ""
    return sign + f"{abs(value):.6f}".replace(".", "p")


def write_report(report: dict[str, Any]) -> None:
    lines = [
        "# TransportModels DG transport-coupling A/B",
        "",
        f"Status: **{report['status']}**.",
        "",
        "| Vg (V) | Vela exponential/direct (dex) | Sentaurus default/DirectQC (dex) | residual effect error (dex) |",
        "|---:|---:|---:|---:|",
    ]
    for row in report["comparison"]:
        vela_effect = row["vela_exponential_direct_dex"]
        effect_error = row["effect_error_dex"]
        vela_text = f"{vela_effect:.6g}" if vela_effect is not None else "failed"
        error_text = f"{effect_error:.6g}" if effect_error is not None else "failed"
        lines.append(
            f"| {row['gate_bias_V']:.2f} | "
            f"{vela_text} | "
            f"{row['sentaurus_default_directqc_dex']:.6g} | "
            f"{error_text} |"
        )
    lines.extend([
        "",
        "The implementation intentionally leaves `direct_band_edge` as the compatibility default.",
        "",
    ])
    REPORT_JSON.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    REPORT_MD.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--report-only", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.check:
        report = json.loads(REPORT_JSON.read_text(encoding="utf-8"))
        assert report["status"] == "complete"
        for item in report["variants"].values():
            path = Path(item["curve_csv"])
            assert digest(path) == item["curve_sha256"]
            assert item["all_converged"] and item["completed_points"] == len(POINTS)
        print("TransportModels DG transport-coupling A/B check: PASS")
        return 0

    phase7 = phase7_module()
    phase7.OUTPUT_ROOT = OUT
    variants = ("direct_band_edge", "sentaurus_exponential")
    common = {
        "config": phase7.GENERATED / "simulation_dg_idvg.json",
        "reference": phase7.REFERENCE / "transportmodels_sentaurus2022_dg_idvg_reference.csv",
        "contact": "gate",
        "current_contact": "drain",
    }
    direct_dir = OUT / "direct_band_edge"
    curves = [
        {
            **common,
            "name": f"frozen2_exponential_{bias_slug(bias)}",
            "label": f"DG frozen-Q sentaurus_exponential Vg={bias:g}",
            "coupling": "sentaurus_exponential",
            "weight": 1.0,
            "points": [bias],
            "initial_state": direct_dir / f"state_bias_{bias_slug(bias)}.csv",
        }
        for bias in POINTS
    ]
    original_make_config = phase7.make_config

    def make_config(spec: dict[str, Any]) -> Path:
        path = original_make_config(spec)
        config = json.loads(path.read_text(encoding="utf-8"))
        config["_comment"] = spec["label"] + " controlled five-point A/B"
        config["solver"]["electron_quantum_potential"]["enabled"] = True
        config["solver"]["electron_quantum_potential"]["transport_coupling"] = spec["coupling"]
        config["solver"]["electron_quantum_potential"]["transport_coupling_weight"] = spec["weight"]
        config["solver"]["electron_quantum_potential"]["coupling_mode"] = "frozen"
        config["solver"]["bandgap_narrowing"] = {
            "model": "old_slotboom",
            "fermi_statistics_correction": True,
        }
        path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
        return path

    phase7.make_config = make_config
    OUT.mkdir(parents=True, exist_ok=True)
    execution = []
    if not args.report_only:
        execution = [phase7.execute_curve(item) for item in curves]

    exponential_dir = OUT / "sentaurus_exponential"
    exponential_dir.mkdir(parents=True, exist_ok=True)
    exponential_csv = exponential_dir / "curve_combined.csv"
    with exponential_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["bias_V", "current_total_A_per_um", "converged", "iterations"])
        for spec in curves:
            item = next(iter(curve(OUT / spec["name"] / "curve_combined.csv").items()))
            bias, row = item
            writer.writerow([bias, row["current_A_per_um"], int(row["converged"]), row["iterations"]])

    results = {name: curve(OUT / name / "curve_combined.csv") for name in variants}
    oracle = json.loads(ORACLE.read_text(encoding="utf-8"))
    sentaurus = {round(float(row["gate_bias_V"]), 10): row for row in oracle["five_point_comparison"]}
    comparison = []
    for bias in POINTS:
        key = round(bias, 10)
        direct = results["direct_band_edge"][key]["current_A_per_um"]
        exponential = results["sentaurus_exponential"][key]["current_A_per_um"]
        sent_effect = math.log10(
            float(sentaurus[key]["default_default_A_per_um"]) /
            float(sentaurus[key]["directqc_default_A_per_um"])
        )
        vela_effect = math.log10(exponential / direct) if exponential > 0.0 else None
        comparison.append({
            "gate_bias_V": bias,
            "vela_direct_A_per_um": direct,
            "vela_exponential_A_per_um": exponential,
            "sentaurus_default_A_per_um": float(sentaurus[key]["default_default_A_per_um"]),
            "sentaurus_directqc_A_per_um": float(sentaurus[key]["directqc_default_A_per_um"]),
            "vela_exponential_direct_dex": vela_effect,
            "sentaurus_default_directqc_dex": sent_effect,
            "effect_error_dex": (
                vela_effect - sent_effect if vela_effect is not None else None
            ),
        })
    artifacts = {}
    for name in variants:
        path = OUT / name / "curve_combined.csv"
        artifacts[name] = {
            "config": (
                str((OUT / name / "config.json").resolve())
                if name == "direct_band_edge"
                else [str((OUT / spec["name"] / "config.json").resolve()) for spec in curves]
            ),
            "curve_csv": str(path.resolve()),
            "curve_sha256": digest(path),
            "completed_points": len(results[name]),
            "all_converged": all(row["converged"] for row in results[name].values()),
        }
    complete = all(item["all_converged"] and item["completed_points"] == len(POINTS) for item in artifacts.values())
    report = {
        "schema": "vela.transportmodels.dg_transport_coupling_ab.v1",
        "as_of": "2026-08-21",
        "status": "complete" if complete else "partial",
        "design": {"gate_biases_V": list(POINTS), "fixed_drain_bias_V": 1.1},
        "execution": execution,
        "variants": artifacts,
        "comparison": comparison,
        "sentaurus_oracle": str(ORACLE.resolve()),
    }
    write_report(report)
    print(json.dumps({"status": report["status"], "comparison": comparison}, indent=2))
    return 0 if complete else 1


if __name__ == "__main__":
    raise SystemExit(main())
