#!/usr/bin/env python3
"""Compare corrected Vela DG states with the new Sentaurus spatial oracles."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any, Iterable


REPO = Path(__file__).resolve().parents[1]
REF = REPO / "build-release/reference_tcad/transportmodels_sentaurus2022"
SENT_MANIFEST = REF / "sentaurus_vm_runs/remaining_spatial_oracles_20260823/remaining_spatial_oracles_manifest.json"
VELA = REF / "vela_baseline/dd_dg_srh_corrected_cold_regression_2026-08-23/runs/dg"
PROFILE_REPORT = REPO / "docs/validation/transportmodels_idvg_spatial_oracle_2026-08-21.json"
OUTPUT = REF / "reports/transportmodels_remaining_spatial_state_compare_20260823"
REPORT_JSON = REPO / "docs/validation/transportmodels_remaining_spatial_state_compare_2026-08-23.json"
REPORT_MD = REPO / "docs/validation/transportmodels_remaining_spatial_state_compare_2026-08-23.md"


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def percentile(values: Iterable[float], fraction: float) -> float:
    ordered = sorted(float(value) for value in values if math.isfinite(float(value)))
    if not ordered:
        return math.nan
    position = (len(ordered) - 1) * fraction
    lo, hi = math.floor(position), math.ceil(position)
    if lo == hi:
        return ordered[lo]
    return ordered[lo] * (hi - position) + ordered[hi] * (position - lo)


def stats(values: Iterable[float]) -> dict[str, float]:
    data = list(values)
    return {
        "count": len(data),
        "median": percentile(data, 0.5),
        "p95": percentile(data, 0.95),
        "maximum": max(data, default=math.nan),
    }


def field(export_dir: Path, name: str, region: int = 3) -> dict[int, float]:
    return {
        int(row["node_id"]): float(row["component0"])
        for row in read_csv(export_dir / "fields" / f"{name}_region{region}.csv")
    }


def vela_slug(value: float) -> str:
    return f"{value:.6f}".replace("-", "m").replace(".", "p")


def vela_state(group: str, bias: float) -> Path:
    if group == "dg_idvd":
        return VELA / f"dg_idvd_curve_state_bias_{vela_slug(bias)}.csv"
    if math.isclose(bias, -1.0):
        return VELA / "dg_idvg_final_bias_relax_final_state.csv"
    return VELA / f"dg_idvg_curve_state_bias_{vela_slug(bias)}.csv"


def compare_case(case: dict[str, Any], profiles: dict[str, set[int]]) -> dict[str, Any]:
    export_dir = Path(case["export_dir"])
    state_path = vela_state(case["group"], case["bias_V"])
    vela_rows = {int(row["node_id"]): row for row in read_csv(state_path)}
    sent = {
        "psi": field(export_dir, "ElectrostaticPotential"),
        "phin": field(export_dir, "eQuasiFermiPotential"),
        "qn": field(export_dir, "eQuantumPotential"),
        "n": field(export_dir, "eDensity"),
        "mobility": field(export_dir, "eMobility"),
        "enormal": field(export_dir, "eEnormal"),
        "eparallel": field(export_dir, "eEparallel"),
    }
    silicon = set(sent["n"])
    zones = {"all_substrate": silicon}
    zones.update({name: nodes & silicon for name, nodes in profiles.items()})
    zone_results: dict[str, Any] = {}
    for name, nodes in zones.items():
        zone_results[name] = {
            "psi_abs_error_mV": stats(
                1.0e3 * abs(float(vela_rows[node]["psi"]) - sent["psi"][node]) for node in nodes
            ),
            "phin_abs_error_mV": stats(
                1.0e3 * abs(float(vela_rows[node]["phin"]) - sent["phin"][node]) for node in nodes
            ),
            "qn_abs_error_mV": stats(
                1.0e3 * abs(float(vela_rows[node]["electron_quantum_potential_V"]) - sent["qn"][node])
                for node in nodes
            ),
            "electron_density_abs_error_dex": stats(
                abs(
                    math.log10(max(float(vela_rows[node]["electrons_m3"]) / 1.0e6, 1.0))
                    - math.log10(max(sent["n"][node], 1.0))
                )
                for node in nodes
            ),
        }
    return {
        "group": case["group"],
        "bias_V": case["bias_V"],
        "gate_bias_V": case["gate_bias_V"],
        "drain_bias_V": case["drain_bias_V"],
        "zones": zone_results,
        "vela_state": str(state_path.resolve()),
        "sentaurus_export": str(export_dir.resolve()),
    }


def markdown(report: dict[str, Any]) -> str:
    lines = [
        "# TransportModels remaining DG spatial-state comparison",
        "",
        "Corrected self-consistent Vela states are compared node-for-node with the new Sentaurus 2022 TDR snapshots.",
        "",
        "| Group | Bias | Qn p95 all/drain (mV) | n p95 all/drain (dex) | phin p95 all/drain (mV) |",
        "|---|---:|---:|---:|---:|",
    ]
    for row in report["cases"]:
        all_zone = row["zones"]["all_substrate"]
        drain = row["zones"]["drain_end"]
        lines.append(
            f"| {row['group']} | {row['bias_V']:.2f} V | "
            f"{all_zone['qn_abs_error_mV']['p95']:.4g} / {drain['qn_abs_error_mV']['p95']:.4g} | "
            f"{all_zone['electron_density_abs_error_dex']['p95']:.4g} / {drain['electron_density_abs_error_dex']['p95']:.4g} | "
            f"{all_zone['phin_abs_error_mV']['p95']:.4g} / {drain['phin_abs_error_mV']['p95']:.4g} |"
        )
    lines.extend(["", f"Raw report: `{REPORT_JSON.resolve()}`", ""])
    return "\n".join(lines)


def main() -> int:
    manifest = json.loads(SENT_MANIFEST.read_text(encoding="utf-8"))
    profile_doc = json.loads(PROFILE_REPORT.read_text(encoding="utf-8"))
    profiles = {
        name: {int(node) for node in spec["node_ids"]}
        for name, spec in profile_doc["profiles"].items()
    }
    cases: list[dict[str, Any]] = []
    for state in manifest["idvd_states"]:
        bias = float(state["drain_bias_V"])
        cases.append(
            {
                "group": "dg_idvd",
                "bias_V": bias,
                "gate_bias_V": 1.0,
                "drain_bias_V": bias,
                "export_dir": state["export_dir"],
            }
        )
    for state in manifest["dg_deep_off_states"]:
        bias = float(state["gate_bias_V"])
        cases.append(
            {
                "group": "dg_idvg_deep_off",
                "bias_V": bias,
                "gate_bias_V": bias,
                "drain_bias_V": 1.1,
                "export_dir": state["export_dir"],
            }
        )
    OUTPUT.mkdir(parents=True, exist_ok=True)
    results = [compare_case(case, profiles) for case in cases]
    report = {
        "schema": "vela.transportmodels.remaining_spatial_state_compare.v1",
        "as_of": "2026-08-23",
        "status": "complete",
        "case_count": len(results),
        "cases": results,
    }
    REPORT_JSON.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    REPORT_MD.write_text(markdown(report), encoding="utf-8")
    print(json.dumps({"status": "complete", "cases": len(results), "report": str(REPORT_JSON)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
