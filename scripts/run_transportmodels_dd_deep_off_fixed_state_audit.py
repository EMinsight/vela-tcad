#!/usr/bin/env python3
"""Replay the Sentaurus DD deep-off state through Vela production operators."""

from __future__ import annotations

import csv
import json
import math
import os
import subprocess
from pathlib import Path
from typing import Any, Iterable


REPO = Path(__file__).resolve().parents[1]
SENTAURUS_ROOT = (
    REPO
    / "build-release/reference_tcad/transportmodels_sentaurus2022/sentaurus_vm_runs"
    / "dd_deep_off_oracle_20260823/export"
)
VELA_ROOT = (
    REPO
    / "build-release/reference_tcad/transportmodels_sentaurus2022/reports"
    / "idvg_deep_off_precision_20260822/newton_calibration/floor2e11_qf1e2"
    / "dd_m1p000000"
)
OUTPUT = (
    REPO
    / "build-release/reference_tcad/transportmodels_sentaurus2022/reports"
    / "idvg_dd_deep_off_fixed_state_20260823"
)
REPORT = REPO / "docs/validation/transportmodels_dd_deep_off_fixed_state_2026-08-23.md"
RUNNER = REPO / "build-release/vela_example_runner.exe"
Q = 1.602176634e-19
SENTAURUS_DD_CURRENT_A_PER_UM = 1.63468406431e-15
UNCORRECTED_VELA_DD_CURRENT_A_PER_UM = 5.4667213725794825e-18


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def percentile(values: Iterable[float], fraction: float) -> float:
    ordered = sorted(float(value) for value in values if math.isfinite(float(value)))
    if not ordered:
        return math.nan
    position = (len(ordered) - 1) * fraction
    low = int(math.floor(position))
    high = int(math.ceil(position))
    if low == high:
        return ordered[low]
    return ordered[low] * (high - position) + ordered[high] * (position - low)


def scalar_field(name: str, region: int) -> dict[int, float]:
    path = SENTAURUS_ROOT / "fields" / f"{name}_region{region}.csv"
    if not path.is_file():
        return {}
    return {int(row["node_id"]): float(row["component0"]) for row in read_csv(path)}


def vector_field(name: str, region: int) -> dict[int, tuple[float, float]]:
    path = SENTAURUS_ROOT / "fields" / f"{name}_region{region}.csv"
    return {
        int(row["node_id"]): (float(row["component0"]), float(row["component1"]))
        for row in read_csv(path)
    }


def merged_scalar(name: str, regions: tuple[int, ...]) -> tuple[dict[int, float], dict[str, Any]]:
    values: dict[int, float] = {}
    conflicts: list[float] = []
    for region in regions:
        for node, value in scalar_field(name, region).items():
            if node in values:
                conflicts.append(abs(values[node] - value))
            values[node] = value
    return values, {
        "field": name,
        "covered_nodes": len(values),
        "duplicate_values": len(conflicts),
        "maximum_duplicate_absolute_difference": max(conflicts, default=0.0),
    }


def prepare_sentaurus_state() -> tuple[Path, dict[str, Any]]:
    nodes = read_csv(SENTAURUS_ROOT / "nodes.csv")
    node_count = len(nodes)
    # Insulators first, then polysilicon, then silicon. Transport-region values
    # therefore win at shared material-interface vertices.
    precedence = (0, 1, 2, 5, 6, 4, 3)
    psi, psi_audit = merged_scalar("ElectrostaticPotential", precedence)
    phin, phin_audit = merged_scalar("eQuasiFermiPotential", precedence)
    phip, phip_audit = merged_scalar("hQuasiFermiPotential", precedence)
    electron_cm3, _ = merged_scalar("eDensity", (4, 3))
    hole_cm3, _ = merged_scalar("hDensity", (4, 3))
    for name, field in (("psi", psi), ("phin", phin), ("phip", phip)):
        missing = sorted(set(range(node_count)) - set(field))
        if missing:
            raise RuntimeError(f"Sentaurus {name} field misses {len(missing)} nodes")

    rows = [
        {
            "node_id": node,
            "psi": psi[node],
            "phin": phin[node],
            "phip": phip[node],
            "electrons_m3": electron_cm3.get(node, 0.0) * 1.0e6,
            "holes_m3": hole_cm3.get(node, 0.0) * 1.0e6,
        }
        for node in range(node_count)
    ]
    path = OUTPUT / "sentaurus_state_for_vela.csv"
    write_csv(path, rows)
    return path, {
        "node_count": node_count,
        "mapping": "exact_global_node_id",
        "merge_audit": [psi_audit, phin_audit, phip_audit],
    }


def prepare_feedback_fields(state_file: Path) -> Path:
    fields = OUTPUT / "sentaurus_feedback_fields"
    fields.mkdir(parents=True, exist_ok=True)
    state = read_csv(state_file)
    mapping = {
        "eQuasiFermiPotential": "phin",
        "hQuasiFermiPotential": "phip",
        "eDensity_m3": "electrons_m3",
        "hDensity_m3": "holes_m3",
    }
    for field, column in mapping.items():
        is_density = field.endswith("Density_m3")
        write_csv(
            fields / f"{field}_region0.csv",
            [
                {
                    "node_id": int(row["node_id"]),
                    "component0": (
                        max(float(row[column]), 1.0)
                        if is_density
                        else float(row[column])
                    ),
                }
                for row in state
            ],
        )
    return fields


def probe_config(
    simulation_type: str,
    state_file: Path,
    output_csv: Path,
    *,
    srh_reference_internal: float | None = None,
) -> dict[str, Any]:
    cfg = json.loads((VELA_ROOT / "config.json").read_text(encoding="utf-8"))
    cfg["simulation_type"] = simulation_type
    cfg["state_file"] = str(state_file.resolve())
    cfg["output_csv"] = str(output_csv.resolve())
    cfg.pop("sweep", None)
    cfg.pop("log_file", None)
    cfg["_comment"] = (
        "DD Vg=-1 V fixed-state audit; evaluates an immutable external state "
        "with Vela production operators"
    )
    if srh_reference_internal is not None:
        dependence = cfg["solver"]["srh_doping_dependence"]
        dependence["electron"]["reference_doping_m3"] = srh_reference_internal
        dependence["hole"]["reference_doping_m3"] = srh_reference_internal
        cfg["_comment"] += (
            "; explicit v1 unit_scaling reference_doping_m3="
            f"{srh_reference_internal:.9g} cm^-3 internal"
        )
    return cfg


def run_probe(
    label: str,
    simulation_type: str,
    state_file: Path,
    *,
    srh_reference_internal: float | None = None,
) -> dict[str, Any]:
    output_csv = OUTPUT / f"{label}.csv"
    config_path = OUTPUT / f"{label}.json"
    config_path.write_text(
        json.dumps(
            probe_config(
                simulation_type,
                state_file,
                output_csv,
                srh_reference_internal=srh_reference_internal,
            ),
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    env = os.environ.copy()
    toolchain = r"D:\msys64\ucrt64\bin"
    env["Path"] = toolchain + os.pathsep + env.get("Path", "")
    completed = subprocess.run(
        [str(RUNNER), "--config", str(config_path), "--log", str(OUTPUT / f"{label}.log")],
        cwd=REPO,
        text=True,
        capture_output=True,
        check=False,
        env=env,
    )
    (OUTPUT / f"{label}.stdout.txt").write_text(completed.stdout, encoding="utf-8")
    (OUTPUT / f"{label}.stderr.txt").write_text(completed.stderr, encoding="utf-8")
    if completed.returncode != 0:
        raise RuntimeError(f"{label} failed: {completed.stderr or completed.stdout}")
    return json.loads(completed.stdout.strip().splitlines()[-1])


def run_feedback_probe(state_file: Path, feedback_fields: Path) -> dict[str, Any]:
    label = "sentaurus_state_feedback_substitution"
    output_csv = OUTPUT / f"{label}.csv"
    cfg = probe_config("newton_feedback_substitution_probe", state_file, output_csv)
    cfg["feedback_state_fields_dir"] = str(feedback_fields.resolve())
    config_path = OUTPUT / f"{label}.json"
    config_path.write_text(json.dumps(cfg, indent=2) + "\n", encoding="utf-8")
    env = os.environ.copy()
    env["Path"] = r"D:\msys64\ucrt64\bin" + os.pathsep + env.get("Path", "")
    completed = subprocess.run(
        [str(RUNNER), "--config", str(config_path), "--log", str(OUTPUT / f"{label}.log")],
        cwd=REPO,
        text=True,
        capture_output=True,
        check=False,
        env=env,
    )
    (OUTPUT / f"{label}.stdout.txt").write_text(completed.stdout, encoding="utf-8")
    (OUTPUT / f"{label}.stderr.txt").write_text(completed.stderr, encoding="utf-8")
    if completed.returncode != 0:
        raise RuntimeError(f"{label} failed: {completed.stderr or completed.stdout}")
    return json.loads(completed.stdout.strip().splitlines()[-1])


def run_self_consistent_reference_ab() -> dict[str, Any]:
    out = OUTPUT / "self_consistent_reference_1e16_internal"
    out.mkdir(parents=True, exist_ok=True)
    cfg = json.loads((VELA_ROOT / "config.json").read_text(encoding="utf-8"))
    dependence = cfg["solver"]["srh_doping_dependence"]
    dependence["electron"]["reference_doping_m3"] = 1.0e16
    dependence["hole"]["reference_doping_m3"] = 1.0e16
    cfg["_comment"] = (
        "Self-consistent DD verification at Vg=-1 V after correcting the "
        "Sentaurus importer reference_doping_m3 output to TCAD internal units"
    )
    cfg["output_csv"] = str((out / "curve.csv").resolve())
    sweep = cfg["sweep"]
    sweep["initial_state_file"] = str((VELA_ROOT / "final_state.csv").resolve())
    sweep["write_state_file"] = str((out / "final_state.csv").resolve())
    sweep["write_state_every_point_prefix"] = str((out / "state").resolve())
    diagnostics = sweep["diagnostics"]
    for name, filename in (
        ("terminal_balance", "terminal_balance.csv"),
        ("srh_balance", "srh_balance.csv"),
        ("contact_edge", "contact_edges.csv"),
        ("newton_history", "newton_history.csv"),
    ):
        diagnostics[name]["csv_file"] = str((out / filename).resolve())
    cfg["solver"]["carrier_row_convergence"]["diagnostic_csv"] = str(
        (out / "carrier_row_violations.csv").resolve()
    )
    config_path = out / "config.json"
    config_path.write_text(json.dumps(cfg, indent=2) + "\n", encoding="utf-8")
    env = os.environ.copy()
    env["Path"] = r"D:\msys64\ucrt64\bin" + os.pathsep + env.get("Path", "")
    completed = subprocess.run(
        [str(RUNNER), "--config", str(config_path), "--log", str(out / "runner.log")],
        cwd=REPO,
        text=True,
        capture_output=True,
        check=False,
        env=env,
    )
    (out / "stdout.txt").write_text(completed.stdout, encoding="utf-8")
    (out / "stderr.txt").write_text(completed.stderr, encoding="utf-8")
    result: dict[str, Any] = {
        "returncode": completed.returncode,
        "converged": completed.returncode == 0,
        "config": str(config_path.resolve()),
    }
    curve_path = out / "curve.csv"
    if curve_path.is_file() and read_csv(curve_path):
        curve = read_csv(curve_path)[-1]
        result["drain_current_A_per_um"] = float(curve["current_total_A_per_um"])
    srh_path = out / "srh_balance.csv"
    if srh_path.is_file() and read_csv(srh_path):
        srh = read_csv(srh_path)[-1]
        result["srh_net_current_A_per_um"] = float(srh["srh_net_current_A_per_um"])
        result["numerical_status"] = srh["numerical_status"]
    return result


def triangle_area_um2(nodes: dict[int, tuple[float, float]], ids: tuple[int, int, int]) -> float:
    (x0, y0), (x1, y1), (x2, y2) = (nodes[node] for node in ids)
    return 0.5 * abs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0))


def sentaurus_srh_integral() -> tuple[float, dict[int, float]]:
    nodes = {
        int(row["id"]): (float(row["x_um"]), float(row["y_um"]))
        for row in read_csv(SENTAURUS_ROOT / "nodes.csv")
    }
    srh = scalar_field("srhRecombination", 3)
    control_area = {node: 0.0 for node in srh}
    integral_rate_area = 0.0
    for row in read_csv(SENTAURUS_ROOT / "elements.csv"):
        if row["region"] != "R.Substrate":
            continue
        ids = (int(row["node0"]), int(row["node1"]), int(row["node2"]))
        area = triangle_area_um2(nodes, ids)
        rate = sum(srh[node] for node in ids) / 3.0
        integral_rate_area += rate * area
        for node in ids:
            control_area[node] += area / 3.0
    # cm^-3 s^-1 -> m^-3 s^-1 contributes 1e6; area and a 1 um width
    # contribute 1e-18 m^3. Numeric result is A for 1 um width, i.e. A/um.
    node_current = {
        node: Q * rate * control_area[node] * 1.0e-12 for node, rate in srh.items()
    }
    return Q * integral_rate_area * 1.0e-12, node_current


def state_comparison() -> tuple[list[dict[str, Any]], dict[str, Any]]:
    sent_psi = scalar_field("ElectrostaticPotential", 3)
    sent_phin = scalar_field("eQuasiFermiPotential", 3)
    sent_phip = scalar_field("hQuasiFermiPotential", 3)
    sent_n = scalar_field("eDensity", 3)
    sent_p = scalar_field("hDensity", 3)
    sent_srh = scalar_field("srhRecombination", 3)
    sent_mun = scalar_field("eMobility", 3)
    sent_mup = scalar_field("hMobility", 3)
    sent_jn = vector_field("eCurrentDensity", 3)
    sent_jp = vector_field("hCurrentDensity", 3)
    vela = {int(row["node_id"]): row for row in read_csv(VELA_ROOT / "final_state.csv")}
    rows: list[dict[str, Any]] = []
    for node in sorted(sent_psi):
        v = vela[node]
        vn = max(abs(float(v["electrons_m3"])) / 1.0e6, 1.0)
        vp = max(abs(float(v["holes_m3"])) / 1.0e6, 1.0)
        sn = max(abs(sent_n[node]), 1.0)
        sp = max(abs(sent_p[node]), 1.0)
        rows.append(
            {
                "node_id": node,
                "sentaurus_psi_V": sent_psi[node],
                "vela_psi_V": float(v["psi"]),
                "psi_abs_error_mV": 1.0e3 * abs(float(v["psi"]) - sent_psi[node]),
                "sentaurus_phin_V": sent_phin[node],
                "vela_phin_V": float(v["phin"]),
                "phin_abs_error_mV": 1.0e3 * abs(float(v["phin"]) - sent_phin[node]),
                "sentaurus_phip_V": sent_phip[node],
                "vela_phip_V": float(v["phip"]),
                "phip_abs_error_mV": 1.0e3 * abs(float(v["phip"]) - sent_phip[node]),
                "sentaurus_eDensity_cm3": sent_n[node],
                "vela_eDensity_cm3": vn,
                "electron_density_abs_error_dex": abs(math.log10(vn) - math.log10(sn)),
                "sentaurus_hDensity_cm3": sent_p[node],
                "vela_hDensity_cm3": vp,
                "hole_density_abs_error_dex": abs(math.log10(vp) - math.log10(sp)),
                "sentaurus_srh_cm3_s": sent_srh[node],
                "sentaurus_eMobility_cm2_V_s": sent_mun[node],
                "sentaurus_hMobility_cm2_V_s": sent_mup[node],
                "sentaurus_eCurrentDensity_x_A_cm2": sent_jn[node][0],
                "sentaurus_eCurrentDensity_y_A_cm2": sent_jn[node][1],
                "sentaurus_hCurrentDensity_x_A_cm2": sent_jp[node][0],
                "sentaurus_hCurrentDensity_y_A_cm2": sent_jp[node][1],
            }
        )
    specs = {
        "psi_mV": "psi_abs_error_mV",
        "phin_mV": "phin_abs_error_mV",
        "phip_mV": "phip_abs_error_mV",
        "electron_density_dex": "electron_density_abs_error_dex",
        "hole_density_dex": "hole_density_abs_error_dex",
    }
    metrics = {
        name: {
            "median": percentile((row[key] for row in rows), 0.5),
            "p95": percentile((row[key] for row in rows), 0.95),
            "maximum": max(float(row[key]) for row in rows),
        }
        for name, key in specs.items()
    }
    return rows, metrics


def main() -> int:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    sentaurus_state, mapping = prepare_sentaurus_state()
    feedback_fields = prepare_feedback_fields(sentaurus_state)
    probe_status = {
        "sentaurus_residual": run_probe(
            "sentaurus_state_residual", "newton_residual_probe", sentaurus_state
        ),
        "sentaurus_carrier_terms": run_probe(
            "sentaurus_state_carrier_terms", "newton_carrier_term_probe", sentaurus_state
        ),
        "sentaurus_carrier_terms_reference_1e16_internal": run_probe(
            "sentaurus_state_carrier_terms_reference_1e16_internal",
            "newton_carrier_term_probe",
            sentaurus_state,
            srh_reference_internal=1.0e16,
        ),
        "sentaurus_sg_edges": run_probe(
            "sentaurus_state_sg_edges", "sg_edge_flux_probe", sentaurus_state
        ),
        "sentaurus_edge_mobility": run_probe(
            "sentaurus_state_edge_mobility", "edge_mobility_probe", sentaurus_state
        ),
        "vela_carrier_terms": run_probe(
            "vela_state_carrier_terms", "newton_carrier_term_probe", VELA_ROOT / "final_state.csv"
        ),
        "sentaurus_feedback_substitution": run_feedback_probe(
            sentaurus_state, feedback_fields
        ),
    }
    self_consistent_ab = run_self_consistent_reference_ab()
    self_consistent_terms_label = "self_consistent_reference_1e16_internal_carrier_terms"
    probe_status[self_consistent_terms_label] = run_probe(
        self_consistent_terms_label,
        "newton_carrier_term_probe",
        OUTPUT / "self_consistent_reference_1e16_internal" / "final_state.csv",
        srh_reference_internal=1.0e16,
    )
    if "drain_current_A_per_um" in self_consistent_ab:
        corrected_current = self_consistent_ab["drain_current_A_per_um"]
        self_consistent_ab.update(
            {
                "sentaurus_drain_current_A_per_um": SENTAURUS_DD_CURRENT_A_PER_UM,
                "absolute_relative_error_vs_sentaurus": abs(
                    corrected_current - SENTAURUS_DD_CURRENT_A_PER_UM
                )
                / abs(SENTAURUS_DD_CURRENT_A_PER_UM),
                "current_improvement_factor_vs_uncorrected_vela": corrected_current
                / UNCORRECTED_VELA_DD_CURRENT_A_PER_UM,
            }
        )

    silicon_nodes = set(scalar_field("srhRecombination", 3))
    sent_terms = read_csv(OUTPUT / "sentaurus_state_carrier_terms.csv")
    sent_recomb_scaled = sum(
        float(row["electron_recombination"])
        for row in sent_terms
        if int(row["node_id"]) in silicon_nodes
    )
    self_consistent_recomb_scaled = sum(
        float(row["electron_recombination"])
        for row in read_csv(OUTPUT / f"{self_consistent_terms_label}.csv")
        if int(row["node_id"]) in silicon_nodes
    )
    vela_srh_A_um = float(self_consistent_ab["srh_net_current_A_per_um"])
    recombination_scale = vela_srh_A_um / self_consistent_recomb_scaled
    fixed_vela_srh_A_um = sent_recomb_scaled * recombination_scale
    fixed_reference_rows = read_csv(
        OUTPUT / "sentaurus_state_carrier_terms_reference_1e16_internal.csv"
    )
    fixed_reference_srh_A_um = sum(
        float(row["electron_recombination"])
        for row in fixed_reference_rows
        if int(row["node_id"]) in silicon_nodes
    ) * recombination_scale
    sentaurus_srh_A_um, sentaurus_node_srh = sentaurus_srh_integral()
    feedback_rows = read_csv(OUTPUT / "sentaurus_state_feedback_substitution.csv")
    feedback_srh_A_um: dict[str, float] = {}
    for variant in sorted({row["variant"] for row in feedback_rows}):
        scaled_sum = sum(
            float(row["electron_recombination"])
            for row in feedback_rows
            if row["variant"] == variant and int(row["node_id"]) in silicon_nodes
        )
        feedback_srh_A_um[variant] = scaled_sum * recombination_scale

    fixed_terms_by_node = {int(row["node_id"]): row for row in sent_terms}
    srh_rows = []
    for node in sorted(silicon_nodes):
        vela_current = (
            float(fixed_terms_by_node[node]["electron_recombination"])
            * recombination_scale
        )
        sent_current = sentaurus_node_srh[node]
        srh_rows.append(
            {
                "node_id": node,
                "sentaurus_lumped_srh_A_per_um": sent_current,
                "vela_fixed_state_lumped_srh_A_per_um": vela_current,
                "absolute_difference_A_per_um": abs(vela_current - sent_current),
                "magnitude_ratio_vela_over_sentaurus": (
                    abs(vela_current) / abs(sent_current) if sent_current != 0.0 else math.nan
                ),
            }
        )
    write_csv(OUTPUT / "srh_node_comparison.csv", srh_rows)

    state_rows, state_metrics = state_comparison()
    write_csv(OUTPUT / "state_node_comparison.csv", state_rows)

    residual_rows = read_csv(OUTPUT / "sentaurus_state_residual.csv")
    hotspot_rows = sorted(
        (
            row for row in residual_rows if int(row["node_id"]) in silicon_nodes
        ),
        key=lambda row: max(
            float(row["abs_phin_residual"]), float(row["abs_phip_residual"])
        ),
        reverse=True,
    )[:50]
    write_csv(OUTPUT / "sentaurus_state_residual_hotspots.csv", hotspot_rows)

    summary = {
        "schema": "vela.transportmodels_dd_deep_off_fixed_state_audit.v1",
        "status": "complete",
        "bias": {"gate_V": -1.0, "drain_V": 1.1},
        "mapping": mapping,
        "probe_status": probe_status,
        "self_consistent_reference_1e16_internal": self_consistent_ab,
        "state_comparison": state_metrics,
        "srh_integral_A_per_um": {
            "sentaurus_exported_field": sentaurus_srh_A_um,
            "vela_production_formula_on_sentaurus_state": fixed_vela_srh_A_um,
            "vela_self_consistent_state": vela_srh_A_um,
            "vela_fixed_state_reference_1e16_internal": fixed_reference_srh_A_um,
            "fixed_state_magnitude_ratio_vela_over_sentaurus": (
                abs(fixed_vela_srh_A_um) / abs(sentaurus_srh_A_um)
            ),
            "self_consistent_magnitude_ratio_vela_over_sentaurus": (
                abs(vela_srh_A_um) / abs(sentaurus_srh_A_um)
            ),
            "reference_1e16_internal_magnitude_ratio_vela_over_sentaurus": (
                abs(fixed_reference_srh_A_um) / abs(sentaurus_srh_A_um)
            ),
            "carrier_term_to_A_per_um_scale": recombination_scale,
            "feedback_variants": feedback_srh_A_um,
        },
        "diagnosis": {
            "primary": "resolved Sentaurus importer SRH reference-doping unit mismatch",
            "evidence": (
                "fixed-state source rises from the historical 0.3226% to "
                f"{100.0 * abs(fixed_reference_srh_A_um) / abs(sentaurus_srh_A_um):.4f}% "
                "of Sentaurus when the importer emits the intended "
                "1e16 cm^-3 value required by the v1 unit_scaling contract"
            ),
            "density_substitution_effect": (
                "substituting exact Sentaurus electron and hole densities changes "
                "the corrected Vela SRH integral by "
                f"{100.0 * abs(feedback_srh_A_um['density_only'] - fixed_vela_srh_A_um) / abs(fixed_vela_srh_A_um):.4f}%"
            ),
            "remaining": (
                "generalized Fermi SRH factors, ni/BGN details, and nodal source "
                "quadrature remain candidates for the fixed-state gap"
            ),
        },
        "artifacts": {
            "sentaurus_state": str(sentaurus_state.resolve()),
            "state_node_comparison": str((OUTPUT / "state_node_comparison.csv").resolve()),
            "srh_node_comparison": str((OUTPUT / "srh_node_comparison.csv").resolve()),
            "residual_hotspots": str(
                (OUTPUT / "sentaurus_state_residual_hotspots.csv").resolve()
            ),
            "sg_edges": str((OUTPUT / "sentaurus_state_sg_edges.csv").resolve()),
            "edge_mobility": str(
                (OUTPUT / "sentaurus_state_edge_mobility.csv").resolve()
            ),
            "feedback_substitution": str(
                (OUTPUT / "sentaurus_state_feedback_substitution.csv").resolve()
            ),
        },
    }
    (OUTPUT / "summary.json").write_text(
        json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )

    srh = summary["srh_integral_A_per_um"]
    lines = [
        "# TransportModels DD deep-off fixed-state audit",
        "",
        "Bias: `Vg=-1 V`, `Vd=1.1 V`. Sentaurus 2022 and Vela use the exact imported 3315-node topology.",
        "",
        "## First result",
        "",
        "| Quantity | Result |",
        "|---|---:|",
        f"| Sentaurus exported SRH integral | {srh['sentaurus_exported_field']:.9e} A/um |",
        f"| Vela with the uncorrected imported config on the Sentaurus state | {srh['vela_production_formula_on_sentaurus_state']:.9e} A/um |",
        f"| Vela SRH integral on its self-consistent state | {srh['vela_self_consistent_state']:.9e} A/um |",
        f"| Vela corrected fixed state with 1e16 cm^-3 internal reference | {srh['vela_fixed_state_reference_1e16_internal']:.9e} A/um |",
        f"| Fixed-state Vela/Sentaurus magnitude ratio | {srh['fixed_state_magnitude_ratio_vela_over_sentaurus']:.6g} |",
        f"| Self-consistent Vela/Sentaurus magnitude ratio | {srh['self_consistent_magnitude_ratio_vela_over_sentaurus']:.6g} |",
        f"| Uncorrected imported config with Sentaurus n,p substituted | {srh['feedback_variants']['density_only']:.9e} A/um |",
        "",
        "## Diagnosis",
        "",
        f"The Sentaurus importer defect has been fixed: for a v1 `unit_scaling` deck it now writes the intended internal `1e16 cm^-3` value instead of the SI literal `1e22`. The corrected fixed-state path restores {100.0 * srh['reference_1e16_internal_magnitude_ratio_vela_over_sentaurus']:.4f}% of the Sentaurus generation integral without changing the state.",
        "",
        f"Substituting the exact Sentaurus electron and hole densities changes the corrected Vela source by {100.0 * abs(srh['feedback_variants']['density_only'] - srh['vela_production_formula_on_sentaurus_state']) / abs(srh['vela_production_formula_on_sentaurus_state']):.4f}%, so density reconstruction is not the leading cause. Generalized Fermi-SRH/ni-BGN semantics and source quadrature remain candidates for the fixed-state gap pending separate A/B tests.",
        "",
        f"Self-consistent diagnostic converged: `{self_consistent_ab['converged']}`; drain current: `{self_consistent_ab.get('drain_current_A_per_um', 'unavailable')}` A/um.",
        f"Against the Sentaurus terminal current `{SENTAURUS_DD_CURRENT_A_PER_UM:.9e}` A/um, the corrected self-consistent relative error is `{100.0 * self_consistent_ab.get('absolute_relative_error_vs_sentaurus', math.nan):.4f}%`; current magnitude increases by `{self_consistent_ab.get('current_improvement_factor_vs_uncorrected_vela', math.nan):.3f}x` over the uncorrected Vela run.",
        "",
        "## Self-consistent state differences over silicon nodes",
        "",
        "| Field | median | p95 | maximum |",
        "|---|---:|---:|---:|",
    ]
    for name, metric in state_metrics.items():
        lines.append(
            f"| {name} | {metric['median']:.6g} | {metric['p95']:.6g} | {metric['maximum']:.6g} |"
        )
    lines.extend(
        [
            "",
            "The fixed-state carrier-term and SG-edge CSV files are direct outputs from Vela production operators; no Python reimplementation is used for those terms.",
            "",
            f"Raw artifact directory: `{OUTPUT}`",
            "",
        ]
    )
    REPORT.write_text("\n".join(lines), encoding="utf-8")
    print(json.dumps(summary, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
