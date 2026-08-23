#!/usr/bin/env python3
"""Compare hard-closed DG deep-off states with Sentaurus spatial oracles."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any, Iterable


REPO = Path(__file__).resolve().parents[1]
REF = REPO / "build-release/reference_tcad/transportmodels_sentaurus2022"
SENT_MANIFEST = REF / "sentaurus_vm_runs/remaining_spatial_oracles_20260823/remaining_spatial_oracles_manifest.json"
PROFILE_REPORT = REPO / "docs/validation/transportmodels_idvg_spatial_oracle_2026-08-21.json"
BASE_RUN = REF / "vela_baseline/dg_quantum_contract_regression_2026-08-23/runs/dg"
STRICT_RUN = REF / "reports/transportmodels_dg_deep_off_strict_20260823/scaled_filter"
OUTPUT = REF / "reports/transportmodels_dg_deep_off_strict_20260823"
REPORT = REPO / "docs/validation/transportmodels_dg_deep_off_strict_spatial_2026-08-23.json"


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


def sent_field(export_dir: Path, name: str, region: int = 3) -> dict[int, float]:
    return {
        int(row["node_id"]): float(row["component0"])
        for row in read_csv(export_dir / "fields" / f"{name}_region{region}.csv")
    }


def baseline_state(bias: float) -> Path:
    if math.isclose(bias, -1.0):
        return BASE_RUN / "dg_idvg_final_bias_relax_final_state.csv"
    slug = f"{abs(bias):.6f}".replace(".", "p")
    return BASE_RUN / f"dg_idvg_curve_state_bias_m{slug}.csv"


def strict_state(bias: float) -> Path:
    slug = "m" + f"{abs(bias):.2f}".replace(".", "p")
    return STRICT_RUN / slug / "final_state.csv"


def compare_state(
    state_path: Path,
    sent: dict[str, dict[int, float]],
    zones: dict[str, set[int]],
) -> dict[str, Any]:
    vela = {int(row["node_id"]): row for row in read_csv(state_path)}
    result: dict[str, Any] = {}
    for zone, nodes in zones.items():
        result[zone] = {
            "psi_abs_error_mV": stats(
                1.0e3 * abs(float(vela[node]["psi"]) - sent["psi"][node])
                for node in nodes
            ),
            "phin_abs_error_mV": stats(
                1.0e3 * abs(float(vela[node]["phin"]) - sent["phin"][node])
                for node in nodes
            ),
            "qn_abs_error_mV": stats(
                1.0e3
                * abs(
                    float(vela[node]["electron_quantum_potential_V"])
                    - sent["qn"][node]
                )
                for node in nodes
            ),
            "electron_density_abs_error_dex": stats(
                abs(
                    math.log10(max(float(vela[node]["electrons_m3"]) / 1.0e6, 1.0))
                    - math.log10(max(sent["n"][node], 1.0))
                )
                for node in nodes
            ),
        }
    return {"state": str(state_path.resolve()), "zones": result}


def main() -> int:
    manifest = json.loads(SENT_MANIFEST.read_text(encoding="utf-8"))
    profiles_doc = json.loads(PROFILE_REPORT.read_text(encoding="utf-8"))
    profiles = {
        name: {int(node) for node in spec["node_ids"]}
        for name, spec in profiles_doc["profiles"].items()
    }
    cases: list[dict[str, Any]] = []
    flat: list[dict[str, Any]] = []
    for case in manifest["dg_deep_off_states"]:
        bias = float(case["gate_bias_V"])
        export_dir = Path(case["export_dir"])
        sent = {
            "psi": sent_field(export_dir, "ElectrostaticPotential"),
            "phin": sent_field(export_dir, "eQuasiFermiPotential"),
            "qn": sent_field(export_dir, "eQuantumPotential"),
            "n": sent_field(export_dir, "eDensity"),
        }
        silicon = set(sent["n"])
        zones = {"all_substrate": silicon}
        zones.update({name: nodes & silicon for name, nodes in profiles.items()})
        baseline = compare_state(baseline_state(bias), sent, zones)
        strict = compare_state(strict_state(bias), sent, zones)
        cases.append(
            {
                "bias_V": bias,
                "sentaurus_export": str(export_dir.resolve()),
                "baseline": baseline,
                "strict_closed": strict,
            }
        )
        for zone in zones:
            row: dict[str, Any] = {"bias_V": bias, "zone": zone}
            for variant_name, variant in (("baseline", baseline), ("strict", strict)):
                metrics = variant["zones"][zone]
                row[f"{variant_name}_psi_p95_mV"] = metrics["psi_abs_error_mV"]["p95"]
                row[f"{variant_name}_phin_p95_mV"] = metrics["phin_abs_error_mV"]["p95"]
                row[f"{variant_name}_qn_p95_mV"] = metrics["qn_abs_error_mV"]["p95"]
                row[f"{variant_name}_n_p95_dex"] = metrics[
                    "electron_density_abs_error_dex"
                ]["p95"]
            flat.append(row)

    OUTPUT.mkdir(parents=True, exist_ok=True)
    csv_path = OUTPUT / "deep_off_strict_spatial_comparison.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(flat[0]))
        writer.writeheader()
        writer.writerows(flat)
    report = {
        "schema": "vela.transportmodels.dg_deep_off_strict_spatial.v1",
        "status": "complete",
        "cases": cases,
        "artifacts": {"comparison_csv": str(csv_path.resolve())},
    }
    REPORT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "complete", "cases": len(cases), "report": str(REPORT.resolve())}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
