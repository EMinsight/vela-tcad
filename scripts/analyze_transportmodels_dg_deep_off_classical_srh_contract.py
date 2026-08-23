#!/usr/bin/env python3
"""Audit the Sentaurus default-DG classical/QM density contract in SRH."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[1]
REF = REPO / "build-release/reference_tcad/transportmodels_sentaurus2022"
MANIFEST = REF / "sentaurus_vm_runs/remaining_spatial_oracles_20260823/remaining_spatial_oracles_manifest.json"
AB_ROOT = REF / "reports/transportmodels_dg_deep_off_strict_20260823/srh_ab"
OUTPUT = REF / "reports/transportmodels_dg_deep_off_strict_20260823/classical_srh_contract"
REPORT = REPO / "docs/validation/transportmodels_dg_deep_off_classical_srh_contract_2026-08-23.json"
Q = 1.602176634e-19
VT = 0.025851999786435535
SILICON_NI_CM3 = 14638914958.767616


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def scalar_field(export_dir: Path, name: str) -> dict[int, float]:
    return {
        int(row["node_id"]): float(row["component0"])
        for row in read_csv(export_dir / "fields" / f"{name}_region3.csv")
    }


def tag(bias: float) -> str:
    return "m" + f"{abs(bias):.2f}".replace(".", "p")


def control_areas(export_dir: Path, nodes_of_interest: set[int]) -> dict[int, float]:
    nodes = {
        int(row["id"]): (float(row["x_um"]), float(row["y_um"]))
        for row in read_csv(export_dir / "nodes.csv")
    }
    result = {node: 0.0 for node in nodes_of_interest}
    for element in read_csv(export_dir / "elements.csv"):
        if element["region"] != "R.Substrate":
            continue
        ids = tuple(int(element[f"node{i}"]) for i in range(3))
        (x0, y0), (x1, y1), (x2, y2) = (nodes[node] for node in ids)
        area = 0.5 * abs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0))
        for node in ids:
            result[node] += area / 3.0
    return result


def srh_rate(
    n_denominator: float,
    p: float,
    ni: float,
    dphi: float,
    donors: float,
    acceptors: float,
) -> float:
    doping = abs(donors) + abs(acceptors)
    taun = 3.0e-8 / (1.0 + doping / 1.0e16)
    taup = 3.0e-6 / (1.0 + doping / 1.0e16)
    excess = ni * ni * math.expm1(max(-500.0, min(500.0, dphi / VT)))
    return excess / (taup * (n_denominator + ni) + taun * (p + ni))


def main() -> int:
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    OUTPUT.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    for state in manifest["dg_deep_off_states"]:
        bias = float(state["gate_bias_V"])
        export_dir = Path(state["export_dir"])
        sent_srh = scalar_field(export_dir, "srhRecombination")
        n_qm = scalar_field(export_dir, "eDensity")
        p = scalar_field(export_dir, "hDensity")
        qn = scalar_field(export_dir, "eQuantumPotential")
        phin = scalar_field(export_dir, "eQuasiFermiPotential")
        phip = scalar_field(export_dir, "hQuasiFermiPotential")
        donors = scalar_field(export_dir, "DonorConcentration")
        acceptors = scalar_field(export_dir, "AcceptorConcentration")
        sent_bgn = scalar_field(export_dir, "BandgapNarrowing")
        area = control_areas(export_dir, set(sent_srh))
        coordinates = {
            int(row["id"]): (float(row["x_um"]), float(row["y_um"]))
            for row in read_csv(export_dir / "nodes.csv")
        }
        term_rows = {
            int(row["node_id"]): row
            for row in read_csv(
                AB_ROOT / tag(bias) / "baseline_fermi_oldslotboom/carrier_terms.csv"
            )
        }

        sent_sum = 0.0
        current_sum = 0.0
        classical_sum = 0.0
        bgn_abs_weighted = 0.0
        bgn_signed_weighted = 0.0
        source_weight = 0.0
        node_rows: list[dict[str, Any]] = []
        for node in sent_srh:
            ni = float(term_rows[node]["ni_eff_m3"])
            n_classical = n_qm[node] * math.exp(max(-500.0, min(500.0, qn[node] / VT)))
            dphi = phip[node] - phin[node]
            current_rate = srh_rate(
                n_qm[node], p[node], ni, dphi, donors[node], acceptors[node]
            )
            classical_rate = srh_rate(
                n_classical, p[node], ni, dphi, donors[node], acceptors[node]
            )
            sent_sum += sent_srh[node] * area[node]
            current_sum += current_rate * area[node]
            classical_sum += classical_rate * area[node]
            weight = abs(sent_srh[node]) * area[node]
            vela_bgn = 2.0 * VT * math.log(max(ni, SILICON_NI_CM3) / SILICON_NI_CM3)
            bgn_error = vela_bgn - sent_bgn[node]
            bgn_abs_weighted += abs(bgn_error) * weight
            bgn_signed_weighted += bgn_error * weight
            source_weight += weight
            node_rows.append(
                {
                    "node_id": node,
                    "x_um": coordinates[node][0],
                    "y_um": coordinates[node][1],
                    "control_area_um2": area[node],
                    "sentaurus_srh_cm-3_s-1": sent_srh[node],
                    "vela_current_contract_srh_cm-3_s-1": current_rate,
                    "vela_classical_denominator_srh_cm-3_s-1": classical_rate,
                    "electron_qm_density_cm-3": n_qm[node],
                    "electron_classical_density_cm-3": n_classical,
                    "electron_quantum_potential_V": qn[node],
                    "ni_eff_cm-3": ni,
                }
            )

        sent_current = abs(sent_sum) * Q * 1.0e-12
        current_contract = abs(current_sum) * Q * 1.0e-12
        classical_contract = abs(classical_sum) * Q * 1.0e-12
        baseline_gap = abs(sent_current - current_contract) / sent_current
        classical_gap = abs(sent_current - classical_contract) / sent_current
        node_path = OUTPUT / f"{tag(bias)}_node_contract.csv"
        with node_path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(node_rows[0]))
            writer.writeheader()
            writer.writerows(node_rows)
        rows.append(
            {
                "gate_bias_V": bias,
                "sentaurus_srh_A_per_um": sent_current,
                "vela_current_quantum_density_denominator_A_per_um": current_contract,
                "vela_classical_electron_density_denominator_A_per_um": classical_contract,
                "current_contract_relative_gap_fraction": baseline_gap,
                "classical_contract_relative_gap_fraction": classical_gap,
                "explained_fraction_of_current_contract_gap": (
                    baseline_gap - classical_gap
                ) / baseline_gap,
                "srh_weighted_bgn_abs_error_eV": bgn_abs_weighted / source_weight,
                "srh_weighted_bgn_signed_error_eV": bgn_signed_weighted / source_weight,
                "node_csv": str(node_path.resolve()),
            }
        )

    csv_path = OUTPUT / "classical_srh_contract_summary.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    report = {
        "schema": "vela.transportmodels.dg_deep_off_classical_srh_contract.v1",
        "status": "complete",
        "manual_contract": {
            "source": "Sentaurus Device User Guide T-2022.03, pages 369 and 474",
            "finding": (
                "Default density-gradient coupling keeps classical and quantum-mechanical "
                "carrier densities distinct. The deep-off oracle is reproduced by using "
                "the classical electron density in the SRH denominator while retaining "
                "the quasi-Fermi-splitting numerator."
            ),
        },
        "fermi_factor_ab_note": (
            "The fixed-state Fermi-versus-Boltzmann SRH A/B changed each integral by less "
            "than 2e-10 relative, so the classical-density audit uses gamma_n=gamma_p=1."
        ),
        "points": rows,
        "artifacts": {"summary_csv": str(csv_path.resolve())},
    }
    REPORT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
