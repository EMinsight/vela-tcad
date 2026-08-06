#!/usr/bin/env python3
"""Audit BVmethods NMOS OldSlotboom and Fermi-Dirac density mapping."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


REPO = Path(__file__).resolve().parents[1]
RUN_ROOT = REPO / "build-release/reference_tcad/bvmethods_sentaurus2018/run01"
DEFAULT_STATE = (
    RUN_ROOT
    / "vela_validation/btbt_e2_iic_sentaurus_ni_fixed6p4_20260805/"
      "states/accepted_state_bias_6p400000.csv"
)
DEFAULT_SENT = (
    RUN_ROOT
    / "sentaurus_iic_multibias_exact_extended_20260803/imported/iic_v6p400000"
)
DEFAULT_DOPING = RUN_ROOT / "mesh_import_sorted_node_order/doping.csv"
DEFAULT_OUT = (
    RUN_ROOT
    / "vela_validation/btbt_e2_iic_sentaurus_ni_fixed6p4_20260805/bgn_fd_audit"
)

KB_J_K = 1.380649e-23
Q_C = 1.602176634e-19
TEMPERATURE_K = 300.0
VT_V = KB_J_K * TEMPERATURE_K / Q_C
SQRT_PI = 1.772453850905516
BEDNARCZYK_COEFFICIENT = 0.75 * SQRT_PI


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def scalar_field(root: Path, name: str) -> dict[int, float]:
    return {
        int(row["node_id"]): float(row["component0"])
        for row in read_rows(root / "fields" / f"{name}_region3.csv")
    }


def fermi_half(eta: float) -> float:
    if eta < -40.0:
        return math.exp(eta)
    exponential = math.exp(max(-eta, -700.0))
    shifted = eta + 1.0
    gaussian = math.exp(-0.17 * shifted * shifted)
    value = eta**4 + 50.0 + 33.6 * eta * (1.0 - 0.68 * gaussian)
    return 1.0 / (
        exponential + BEDNARCZYK_COEFFICIENT * value ** -0.375
    )


def inverse_fermi_half(value: float) -> float:
    lower, upper = -500.0, 100.0
    for _ in range(160):
        midpoint = 0.5 * (lower + upper)
        if fermi_half(midpoint) < value:
            lower = midpoint
        else:
            upper = midpoint
    return 0.5 * (lower + upper)


def old_slotboom_delta_eg(total_impurity_cm3: float) -> float:
    x = math.log(total_impurity_cm3 / 1.0e17)
    return max(0.0, 9.0e-3 * (x + math.sqrt(x * x + 0.5)))


def effective_ni(material_ni_cm3: float, delta_eg_ev: float) -> float:
    return material_ni_cm3 * math.exp(delta_eg_ev / (2.0 * VT_V))


def fermi_bgn_correction(
    donors_cm3: float, acceptors_cm3: float, nc_cm3: float, nv_cm3: float
) -> float:
    def species(concentration: float, density_of_states: float) -> float:
        if concentration <= 0.0:
            return 0.0
        ratio = concentration / density_of_states
        return max(0.0, inverse_fermi_half(ratio) - math.log(ratio))

    return VT_V * (species(donors_cm3, nc_cm3) + species(acceptors_cm3, nv_cm3))


def electron_density(
    ni_eff_cm3: float, nc_cm3: float, psi: float, phin: float
) -> tuple[float, float]:
    eta = (psi - phin) / VT_V + math.log(ni_eff_cm3 / nc_cm3)
    return nc_cm3 * fermi_half(eta), eta


def inferred_ni_eff(
    density_cm3: float, nc_cm3: float, psi: float, phin: float
) -> float:
    eta = inverse_fermi_half(density_cm3 / nc_cm3)
    return nc_cm3 * math.exp(eta - (psi - phin) / VT_V)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vela-state", type=Path, default=DEFAULT_STATE)
    parser.add_argument("--sentaurus-root", type=Path, default=DEFAULT_SENT)
    parser.add_argument("--doping", type=Path, default=DEFAULT_DOPING)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--nodes", nargs="+", type=int, default=[344, 1335])
    parser.add_argument("--default-ni-cm3", type=float, default=1.0e10)
    parser.add_argument(
        "--sentaurus-ni-cm3", type=float, default=1.4638914958767616e10
    )
    parser.add_argument("--nc-cm3", type=float, default=2.8e19)
    parser.add_argument("--nv-cm3", type=float, default=1.04e19)
    args = parser.parse_args()

    vela = {int(row["node_id"]): row for row in read_rows(args.vela_state)}
    doping = {int(row["node_id"]): row for row in read_rows(args.doping)}
    sent_psi = scalar_field(args.sentaurus_root, "ElectrostaticPotential")
    sent_phin = scalar_field(args.sentaurus_root, "eQuasiFermiPotential")
    sent_n = scalar_field(args.sentaurus_root, "eDensity")

    results: list[dict[str, float | int]] = []
    for node in args.nodes:
        donors = float(doping[node]["donors_cm3"])
        acceptors = float(doping[node]["acceptors_cm3"])
        total_impurity = donors + acceptors
        delta_eg = old_slotboom_delta_eg(total_impurity)
        default_ni_eff = effective_ni(args.default_ni_cm3, delta_eg)
        corrected_ni_eff = effective_ni(args.sentaurus_ni_cm3, delta_eg)
        fermi_correction = fermi_bgn_correction(
            donors, acceptors, args.nc_cm3, args.nv_cm3
        )
        corrected_fermi_ni_eff = effective_ni(
            args.sentaurus_ni_cm3, delta_eg + fermi_correction
        )
        sent_density_cm3 = sent_n[node]
        inferred = inferred_ni_eff(
            sent_density_cm3,
            args.nc_cm3,
            sent_psi[node],
            sent_phin[node],
        )
        inferred_delta_eg = 2.0 * VT_V * math.log(
            inferred / args.sentaurus_ni_cm3
        )
        predicted_default, eta_default = electron_density(
            default_ni_eff,
            args.nc_cm3,
            sent_psi[node],
            sent_phin[node],
        )
        predicted_corrected, eta_corrected = electron_density(
            corrected_ni_eff,
            args.nc_cm3,
            sent_psi[node],
            sent_phin[node],
        )
        predicted_corrected_fermi, eta_corrected_fermi = electron_density(
            corrected_fermi_ni_eff,
            args.nc_cm3,
            sent_psi[node],
            sent_phin[node],
        )
        vela_density_cm3 = float(vela[node]["electrons_m3"]) / 1.0e6
        results.append({
            "node_id": node,
            "donors_cm3": donors,
            "acceptors_cm3": acceptors,
            "total_impurity_cm3": total_impurity,
            "old_slotboom_delta_Eg_eV": delta_eg,
            "fermi_statistics_delta_Eg_correction_eV": fermi_correction,
            "corrected_total_delta_Eg_eV": delta_eg + fermi_correction,
            "default_ni_eff_cm3": default_ni_eff,
            "sentaurus_material_ni_eff_cm3": corrected_ni_eff,
            "sentaurus_material_fermi_corrected_ni_eff_cm3": corrected_fermi_ni_eff,
            "sentaurus_inferred_ni_eff_cm3": inferred,
            "corrected_over_inferred_ni_eff": corrected_ni_eff / inferred,
            "sentaurus_inferred_delta_Eg_eV": inferred_delta_eg,
            "inferred_minus_old_slotboom_delta_Eg_eV": inferred_delta_eg - delta_eg,
            "inferred_over_old_slotboom_delta_Eg": inferred_delta_eg / delta_eg,
            "sentaurus_psi_V": sent_psi[node],
            "sentaurus_phin_V": sent_phin[node],
            "sentaurus_eDensity_cm3": sent_density_cm3,
            "vela_eDensity_cm3": vela_density_cm3,
            "vela_over_sentaurus_density": vela_density_cm3 / sent_density_cm3,
            "predicted_from_sent_state_default_ni_cm3": predicted_default,
            "predicted_from_sent_state_corrected_ni_cm3": predicted_corrected,
            "default_prediction_over_sentaurus": predicted_default / sent_density_cm3,
            "corrected_prediction_over_sentaurus": predicted_corrected / sent_density_cm3,
            "corrected_fermi_prediction_over_sentaurus": (
                predicted_corrected_fermi / sent_density_cm3
            ),
            "default_reduced_fermi_eta": eta_default,
            "corrected_reduced_fermi_eta": eta_corrected,
            "corrected_fermi_reduced_eta": eta_corrected_fermi,
            "default_band_offset_proxy_V": VT_V * math.log(args.nc_cm3 / default_ni_eff),
            "corrected_band_offset_proxy_V": VT_V * math.log(args.nc_cm3 / corrected_ni_eff),
            "sentaurus_inferred_band_offset_proxy_V": VT_V * math.log(args.nc_cm3 / inferred),
        })

    args.out_dir.mkdir(parents=True, exist_ok=True)
    fields = list(results[0])
    with (args.out_dir / "focus_nodes.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(results)

    summary = {
        "temperature_K": TEMPERATURE_K,
        "thermal_voltage_V": VT_V,
        "nc_cm3": args.nc_cm3,
        "nv_cm3": args.nv_cm3,
        "default_material_ni_cm3": args.default_ni_cm3,
        "sentaurus_material_ni_cm3": args.sentaurus_ni_cm3,
        "nodes": results,
    }
    (args.out_dir / "summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(args.out_dir / "summary.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
