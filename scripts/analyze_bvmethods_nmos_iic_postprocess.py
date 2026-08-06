#!/usr/bin/env python3
"""Close the BVmethods NMOS postprocessed IIC validation ledger."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any


Q_C = 1.602176634e-19
TCAD_SOURCE_INTEGRAL_TO_PER_M_S = 1.0e-6
PER_M_DEPTH_TO_PER_UM_DEPTH = 1.0e-6
TCAD_PARTICLE_FLUX_TO_PER_M2_S = 1.0e4


def avalanche_current_from_native_source(source_integral: float) -> float:
    """Convert alpha[cm^-1]*flux[cm^-2/s]*area[um^2] to A/um."""
    return (
        Q_C
        * source_integral
        * TCAD_SOURCE_INTEGRAL_TO_PER_M_S
        * PER_M_DEPTH_TO_PER_UM_DEPTH
    )


def current_density_from_native_particle_flux(flux: float) -> float:
    """Convert the unit_scaling particle-flux dump to conventional A/m^2."""
    return Q_C * TCAD_PARTICLE_FLUX_TO_PER_M2_S * flux
REPO = Path(__file__).resolve().parents[1]
RUN_ROOT = REPO / "build-release/reference_tcad/bvmethods_sentaurus2018/run01"
DEFAULT_VELA = RUN_ROOT / "vela_validation/iic_postprocess_20260803"
DEFAULT_SENT_CURVE = RUN_ROOT / "analysis/curves/ABA_coupled.csv"
DEFAULT_SENT_SPATIAL = RUN_ROOT / "imported/aba_coupled"
DEFAULT_LOW_BIAS = (
    RUN_ROOT
    / "vela_validation/fermi_dirac_20260802/low_bias_strict/postprocess_only"
)
DEFAULT_OUTPUT = DEFAULT_VELA / "analysis"


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def f(row: dict[str, str], key: str, default: float = math.nan) -> float:
    value = row.get(key, "")
    return float(value) if value not in (None, "") else default


def write_csv(path: Path, records: list[dict[str, Any]]) -> None:
    if not records:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def interpolate(points: list[dict[str, str]], bias: float, key: str) -> float:
    ordered = sorted(points, key=lambda row: f(row, "inner_voltage_V"))
    for left, right in zip(ordered, ordered[1:]):
        x0 = f(left, "inner_voltage_V")
        x1 = f(right, "inner_voltage_V")
        if x0 <= bias <= x1:
            y0 = f(left, key)
            y1 = f(right, key)
            if y0 != 0.0 and y1 != 0.0 and y0 * y1 > 0.0:
                sign = math.copysign(1.0, y0)
                log_y = math.log10(abs(y0)) + (bias - x0) / (x1 - x0) * (
                    math.log10(abs(y1)) - math.log10(abs(y0))
                )
                return sign * 10.0**log_y
            return y0 + (bias - x0) / (x1 - x0) * (y1 - y0)
    return math.nan


def probe_bias(path: Path) -> float:
    record = rows(path / "avalanche_summary.csv")[0]
    return f(record, "bias_V")


def load_probe(path: Path, sent_curve: list[dict[str, str]]) -> dict[str, Any]:
    avalanche = rows(path / "avalanche_summary.csv")[0]
    sweep = rows(path / "sweep.csv")[-1]
    methods = rows(path / "terminal_current_method_compare.csv")
    contacts = {row["contact"]: f(row, "I_sgflux_A_per_um") for row in methods}
    bias = f(avalanche, "bias_V")
    drain = contacts["drain"]
    source = contacts.get("source", 0.0)
    substrate = contacts.get("substrate", 0.0)
    source_integral = f(avalanche, "sum_edge_source_integral")
    iava = avalanche_current_from_native_source(source_integral)
    sent_id = interpolate(sent_curve, bias, "drain_total_current_A_per_um")
    sent_iava = interpolate(sent_curve, bias, "avalanche_current_A_per_um")
    sent_phi_e = interpolate(sent_curve, bias, "phi_electron")
    sent_phi_h = interpolate(sent_curve, bias, "phi_hole")
    return {
        "bias_V": bias,
        "vela_Id_A_per_um": drain,
        "sentaurus_Id_A_per_um_log_interp": sent_id,
        "abs_vela_over_sentaurus_Id": abs(drain / sent_id) if sent_id else math.nan,
        "vela_max_electric_field_V_per_m": f(sweep, "max_electric_field_V_per_m"),
        "vela_max_electron_alpha_m_inv": f(avalanche, "max_electron_alpha_m_inv"),
        "vela_max_hole_alpha_m_inv": f(avalanche, "max_hole_alpha_m_inv"),
        "vela_source_integral_per_s_per_m_depth": source_integral,
        "vela_avalanche_current_A_per_um": iava,
        "sentaurus_avalanche_current_A_per_um_log_interp": sent_iava,
        "abs_vela_over_sentaurus_avalanche_current": (
            abs(iava / sent_iava) if sent_iava else math.nan
        ),
        "vela_Iava_over_abs_Id": iava / abs(drain) if drain else math.inf,
        "sentaurus_Iava_over_abs_Id": abs(sent_iava / sent_id) if sent_id else math.nan,
        "sentaurus_phi_electron": sent_phi_e,
        "sentaurus_phi_hole": sent_phi_h,
        "vela_terminal_kcl_relative": (
            abs(drain + source + substrate)
            / max(abs(drain), abs(source), abs(substrate), 1.0e-300)
        ),
        "vela_sg_vs_residual_drain_relative": abs(
            drain - f(next(row for row in methods if row["contact"] == "drain"),
                      "I_residual_A_per_um")
        ) / max(abs(drain), 1.0e-300),
        "qf_bounds_violations": int(float(sweep.get("qf_bounds_violations", "0") or 0)),
        "newton_convergence_reason": sweep.get("newton_convergence_reason", ""),
    }


def load_low_bias(path: Path, sent_curve: list[dict[str, str]]) -> list[dict[str, Any]]:
    avalanche = {round(f(row, "bias_V"), 12): row for row in rows(path / "avalanche_summary.csv")}
    records: list[dict[str, Any]] = []
    for sweep in rows(path / "sweep.csv"):
        bias = f(sweep, "bias_V")
        if bias <= 0.0:
            continue
        source = f(avalanche[round(bias, 12)], "sum_edge_source_integral")
        iava = avalanche_current_from_native_source(source)
        drain = f(sweep, "current_total_A_per_um")
        sent_id = interpolate(sent_curve, bias, "drain_total_current_A_per_um")
        sent_iava = interpolate(sent_curve, bias, "avalanche_current_A_per_um")
        records.append({
            "bias_V": bias,
            "vela_Id_A_per_um": drain,
            "sentaurus_Id_A_per_um_log_interp": sent_id,
            "vela_avalanche_current_A_per_um": iava,
            "sentaurus_avalanche_current_A_per_um_log_interp": sent_iava,
            "vela_Iava_over_abs_Id": iava / abs(drain),
            "sentaurus_Iava_over_abs_Id": abs(sent_iava / sent_id),
        })
    return records


def top_edges(probe_paths: list[Path], top_n: int) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    for path in sorted(probe_paths, key=probe_bias):
        edge_rows = rows(path / "sg_avalanche_edges.csv")
        total = sum(f(row, "edge_source_integral", 0.0) for row in edge_rows)
        ranked = sorted(edge_rows, key=lambda row: abs(f(row, "edge_source_integral", 0.0)), reverse=True)
        for rank, row in enumerate(ranked[:top_n], start=1):
            edge_source = f(row, "edge_source_integral", 0.0)
            output.append({
                "bias_V": f(row, "bias_V"),
                "rank": rank,
                "edge_id": row["edge_id"],
                "node0": row["node0"],
                "node1": row["node1"],
                "x_mid_um": 0.5 * (f(row, "x0_um") + f(row, "x1_um")),
                "y_mid_um": 0.5 * (f(row, "y0_um") + f(row, "y1_um")),
                "edge_class": row["edge_class"],
                "electric_field_V_per_m": f(row, "electric_field_V_per_m"),
                "electron_impact_field_V_per_m": f(row, "electron_impact_field_V_per_m"),
                "hole_impact_field_V_per_m": f(row, "hole_impact_field_V_per_m"),
                "electron_alpha_m_inv": f(row, "electron_alpha_m_inv"),
                "hole_alpha_m_inv": f(row, "hole_alpha_m_inv"),
                "electron_source_integral": f(row, "electron_source_integral"),
                "hole_source_integral": f(row, "hole_source_integral"),
                "edge_source_integral": edge_source,
                "edge_fraction_of_signed_total": edge_source / total if total else math.nan,
                "qG_contribution_A_per_um": avalanche_current_from_native_source(edge_source),
            })
    return output


def field_records(path: Path, name: str) -> list[dict[str, str]]:
    return rows(path / "fields" / f"{name}_region3.csv")


def sentaurus_spatial_peaks(path: Path) -> dict[str, Any]:
    nodes = {int(row["id"]): row for row in rows(path / "nodes.csv")}

    def scalar_peak(name: str, scale: float) -> dict[str, Any]:
        peak = max(field_records(path, name), key=lambda row: abs(f(row, "component0")))
        node_id = int(peak["node_id"])
        return {
            "value": abs(f(peak, "component0")) * scale,
            "node_id": node_id,
            "x_um": f(nodes[node_id], "x_um"),
            "y_um": f(nodes[node_id], "y_um"),
        }

    def vector_peak(name: str, scale: float) -> dict[str, Any]:
        data = field_records(path, name)
        peak = max(data, key=lambda row: math.hypot(f(row, "component0"), f(row, "component1")))
        node_id = int(peak["node_id"])
        return {
            "value": math.hypot(f(peak, "component0"), f(peak, "component1")) * scale,
            "node_id": node_id,
            "x_um": f(nodes[node_id], "x_um"),
            "y_um": f(nodes[node_id], "y_um"),
        }

    return {
        "electric_field_V_per_m": vector_peak("ElectricField", 100.0),
        "electron_alpha_m_inv": scalar_peak("eAlphaAvalanche", 100.0),
        "hole_alpha_m_inv": scalar_peak("hAlphaAvalanche", 100.0),
        "electron_current_density_A_per_m2": vector_peak("eCurrentDensity", 1.0e4),
        "hole_current_density_A_per_m2": vector_peak("hCurrentDensity", 1.0e4),
    }


def vela_spatial_peaks(probe_path: Path) -> dict[str, Any]:
    data = rows(probe_path / "sg_avalanche_edges.csv")

    def peak(key: str) -> dict[str, Any]:
        row = max(data, key=lambda item: abs(f(item, key)))
        return {
            "value": abs(f(row, key)),
            "edge_id": int(row["edge_id"]),
            "x_um": 0.5 * (f(row, "x0_um") + f(row, "x1_um")),
            "y_um": 0.5 * (f(row, "y0_um") + f(row, "y1_um")),
        }

    return {
        "electric_field_V_per_m": peak("electric_field_V_per_m"),
        "electron_alpha_m_inv": peak("electron_alpha_m_inv"),
        "hole_alpha_m_inv": peak("hole_alpha_m_inv"),
        "electron_current_density_A_per_m2": peak("electron_flux_proxy"),
        "hole_current_density_A_per_m2": peak("hole_flux_proxy"),
    }


def spatial_compare(vela: dict[str, Any], sent: dict[str, Any]) -> list[dict[str, Any]]:
    output = []
    for quantity in vela:
        v = vela[quantity]
        s = sent[quantity]
        vela_value = v["value"]
        if quantity.endswith("current_density_A_per_m2"):
            # The SG edge dump stores the carrier particle-flux magnitude.
            vela_value = current_density_from_native_particle_flux(vela_value)
        output.append({
            "quantity": quantity,
            "vela_peak": vela_value,
            "sentaurus_peak": s["value"],
            "abs_vela_over_sentaurus": vela_value / s["value"] if s["value"] else math.nan,
            "vela_support_id": v["edge_id"],
            "vela_x_um": v["x_um"],
            "vela_y_um": v["y_um"],
            "sentaurus_support_id": s["node_id"],
            "sentaurus_x_um": s["x_um"],
            "sentaurus_y_um": s["y_um"],
        })
    return output


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vela-root", type=Path, default=DEFAULT_VELA)
    parser.add_argument("--sentaurus-curve", type=Path, default=DEFAULT_SENT_CURVE)
    parser.add_argument("--sentaurus-spatial", type=Path, default=DEFAULT_SENT_SPATIAL)
    parser.add_argument("--low-bias-root", type=Path, default=DEFAULT_LOW_BIAS)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--top-n", type=int, default=10)
    args = parser.parse_args()

    sent_curve = rows(args.sentaurus_curve)
    probe_paths = [
        path / "postprocess_only"
        for path in (args.vela_root / "probes").iterdir()
        if (path / "postprocess_only/avalanche_summary.csv").exists()
    ]
    scalar = [load_probe(path, sent_curve) for path in sorted(probe_paths, key=probe_bias)]
    low_bias = load_low_bias(args.low_bias_root, sent_curve)
    hotspots = top_edges(probe_paths, args.top_n)
    near_iic = min(probe_paths, key=lambda path: abs(probe_bias(path) - 6.377494277837012))
    spatial = spatial_compare(
        vela_spatial_peaks(near_iic), sentaurus_spatial_peaks(args.sentaurus_spatial)
    )

    write_csv(args.out_dir / "iic_scalar_compare.csv", scalar)
    write_csv(args.out_dir / "iic_low_bias_compare.csv", low_bias)
    write_csv(args.out_dir / "iic_hotspot_top_edges.csv", hotspots)
    write_csv(args.out_dir / "iic_spatial_peak_compare_6p38V.csv", spatial)

    all_positive = low_bias + scalar
    summary = {
        "sentaurus_iic_BV_V": 6.377494277837012,
        "vela_bias_coverage_V": [record["bias_V"] for record in scalar],
        "minimum_tested_positive_bias_V": min(record["bias_V"] for record in all_positive),
        "minimum_vela_Iava_over_abs_Id": min(record["vela_Iava_over_abs_Id"] for record in all_positive),
        "vela_iic_crossing_bracketed": any(
            left["vela_Iava_over_abs_Id"] < 1.0 <= right["vela_Iava_over_abs_Id"]
            for left, right in zip(
                sorted(all_positive, key=lambda row: row["bias_V"]),
                sorted(all_positive, key=lambda row: row["bias_V"])[1:],
            )
        ),
        "near_iic_probe_V": probe_bias(near_iic),
        "near_iic_scalar": min(scalar, key=lambda row: abs(row["bias_V"] - 6.377494277837012)),
        "spatial_peak_compare": spatial,
        "sentaurus_multistate_spatial_status": "unavailable_vm_unreachable",
        "interpretation": (
            "Vela postprocessed Iava/|Id| is already above unity at the lowest "
            "positive tested bias, so no physical IIC crossing is bracketed."
        ),
    }
    args.out_dir.mkdir(parents=True, exist_ok=True)
    (args.out_dir / "iic_summary.json").write_text(
        json.dumps(summary, indent=2, allow_nan=True) + "\n", encoding="utf-8"
    )

    near = summary["near_iic_scalar"]
    report = [
        "# BVmethods NMOS postprocess-only IIC closure (2026-08-03)",
        "",
        "## Outcome",
        "",
        f"- Sentaurus IIC reference: `{summary['sentaurus_iic_BV_V']:.9f} V`.",
        f"- Vela positive-bias coverage starts at `{summary['minimum_tested_positive_bias_V']:.3g} V` and its minimum tested `Iava/|Id|` is `{summary['minimum_vela_Iava_over_abs_Id']:.6e}`.",
        "- Vela does not bracket an IIC crossing: the postprocessed avalanche current already exceeds the terminal current in the millivolt range.",
        f"- At `{near['bias_V']:.2f} V`, Vela gives `Id={near['vela_Id_A_per_um']:.6e} A/um`, `Iava={near['vela_avalanche_current_A_per_um']:.6e} A/um`, and `Iava/|Id|={near['vela_Iava_over_abs_Id']:.6e}`.",
        f"- At the same voltage, log-interpolated Sentaurus gives `Id={near['sentaurus_Id_A_per_um_log_interp']:.6e} A/um`, `Iava={near['sentaurus_avalanche_current_A_per_um_log_interp']:.6e} A/um`, and `Iava/|Id|={near['sentaurus_Iava_over_abs_Id']:.6e}`.",
        "",
        "## Numerical closure",
        "",
        f"- The Vela drain SG-flux and residual current agree to `{near['vela_sg_vs_residual_drain_relative']:.3e}` relative error.",
        f"- Drain/source/substrate KCL mismatch at the near-IIC probe is `{near['vela_terminal_kcl_relative']:.3e}`; this point used a relaxed continuation tolerance and has `{near['qf_bounds_violations']}` QF-bound warnings.",
        "",
        "## Spatial comparison",
        "",
        "The local Sentaurus comparison uses the imported final `aba_coupled` TDR. Vela values are edge peaks at 6.38 V; Sentaurus values are node peaks at its final IIC state, so support locations are reported and ratios are diagnostic rather than a same-support norm.",
        "",
        "| quantity | Vela peak | Sentaurus peak | Vela/Sentaurus |",
        "|---|---:|---:|---:|",
    ]
    for row in spatial:
        report.append(
            f"| {row['quantity']} | {row['vela_peak']:.6e} | {row['sentaurus_peak']:.6e} | {row['abs_vela_over_sentaurus']:.6e} |"
        )
    report += [
        "",
        "## Remaining oracle gap",
        "",
        "The local Sentaurus VM did not respond, so 1/2/4/5/6 V multistate TDR fields could not be regenerated. Those checkpoints currently have terminal-current and integrated-avalanche comparisons only; the exact 6.377 V neighborhood has the local field/alpha comparison above.",
        "",
        "## Decision",
        "",
        "Do not use the current Vela postprocessed source to declare BV. First remove the low-bias avalanche floor and the high-bias state/current discrepancy, then rerun this ledger and require an actual `Iava/|Id|=1` bracket near 6.377 V.",
    ]
    (args.out_dir / "iic_summary.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    print(args.out_dir / "iic_summary.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
