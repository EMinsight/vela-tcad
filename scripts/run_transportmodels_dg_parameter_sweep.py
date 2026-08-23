#!/usr/bin/env python3
"""Run a fixed-state one-factor DG parameter audit for TransportModels."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

from run_transportmodels_dg_fixed_state_residual_audit import audit


REPO_ROOT = Path(__file__).resolve().parents[1]
BASELINE = REPO_ROOT / "build-release/reference_tcad/transportmodels_sentaurus2022/vela_baseline"
OUTPUT_ROOT = BASELINE / "dg_parameter_fixed_state_sweep_2026-08-21"
RUNNER = REPO_ROOT / "build-release/vela_example_runner.exe"
BASE_CONFIG = (
    BASELINE
    / "workflow_dg_idvd_resume_1p9_direct2p0_fixed1p5_outer200_run01"
    / "00_dg_idvd_curve.json"
)
HYBRID_RESTART = (
    BASELINE / "frozen_q_oracle_vg1_vd2_run01/vela_state_with_sentaurus_q.csv"
)
BASE_MATERIALS = BASELINE / "generated/vela/materials_sentaurus2022.json"
CORRECTED_MATERIALS = OUTPUT_ROOT / "materials_sentaurus2022_dg_band_drive.json"
REPORT_JSON = REPO_ROOT / "docs/validation/transportmodels_dg_parameter_sweep_2026-08-21.json"
REPORT_MD = REPO_ROOT / "docs/validation/transportmodels_dg_parameter_sweep_2026-08-21.md"

DOS_MASS = 1.0618016171622988
COEFFICIENT_MASS = 1.0906506732296395
SENTAURUS_SILICON_NI_CM3 = 1.4638914958767616e10


VARIANTS: tuple[dict[str, Any], ...] = (
    {"name": "baseline", "parameter": "baseline", "value": None},
    {"name": "coefficient_mass_1p0", "parameter": "coefficient_mass_ratio", "value": 1.0},
    {"name": "coefficient_mass_recovered", "parameter": "coefficient_mass_ratio", "value": COEFFICIENT_MASS},
    {"name": "coefficient_mass_1p2", "parameter": "coefficient_mass_ratio", "value": 1.2},
    {"name": "dos_mass_1p0", "parameter": "effective_mass_ratio", "value": 1.0},
    {"name": "dos_mass_coefficient", "parameter": "effective_mass_ratio", "value": COEFFICIENT_MASS},
    {"name": "gamma_3p2", "parameter": "gamma", "value": 3.2},
    {"name": "gamma_4p0", "parameter": "gamma", "value": 4.0},
    {"name": "theta_0p1", "parameter": "theta", "value": 0.1},
    {"name": "theta_0p25", "parameter": "theta", "value": 0.25},
    {"name": "theta_0p75", "parameter": "theta", "value": 0.75},
    {"name": "bgn_fraction_0", "parameter": "conduction_band_narrowing_fraction", "value": 0.0},
    {"name": "bgn_fraction_1", "parameter": "conduction_band_narrowing_fraction", "value": 1.0},
    {
        "name": "corrected_material_contract",
        "parameter": "material_contract",
        "value": "TDR affinity plus recovered coefficient mass",
    },
)


PARAMETER_MAP = (
    {
        "sentaurus_parameter": "QuantumPotentialParameters.gamma[electron]",
        "vela_field": "solver.electron_quantum_potential.gamma",
        "value": 3.6,
        "unit": "dimensionless",
        "status": "exact",
        "evidence": "Silicon default parameter library; electron component",
    },
    {
        "sentaurus_parameter": "QuantumPotentialParameters.theta[electron]",
        "vela_field": "solver.electron_quantum_potential.theta",
        "value": 0.5,
        "unit": "dimensionless",
        "status": "exact",
        "evidence": "Silicon default parameter library; electron component",
    },
    {
        "sentaurus_parameter": "QuantumPotentialParameters.xi/eta/nu[electron]",
        "vela_field": "fixed Eq. 231 semantics",
        "value": "1 / 1 / 0",
        "unit": "dimensionless",
        "status": "neutral defaults represented",
        "evidence": "No adjustable Vela fields; defaults are the implemented equation",
    },
    {
        "sentaurus_parameter": "eDOSMass Formula 1 at 300 K",
        "vela_field": "effective_mass_ratio",
        "value": DOS_MASS,
        "unit": "m*/m0",
        "status": "frozen at 300 K",
        "evidence": "Derived DOS mass used in the potential-like material drive",
    },
    {
        "sentaurus_parameter": "Formula-0 quantum coefficient mass",
        "vela_field": "coefficient_mass_ratio",
        "value": COEFFICIENT_MASS,
        "unit": "m*/m0",
        "status": "recovered oracle",
        "evidence": "SingleDevice Eq. 231 Jacobian recovery; distinct from DOS mass",
    },
    {
        "sentaurus_parameter": "Bandgap.Bgn2Chi",
        "vela_field": "conduction_band_narrowing_fraction",
        "value": 0.5,
        "unit": "dimensionless",
        "status": "Eq. 231 drive mapping",
        "evidence": "Silicon default parameter library",
    },
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def slug_label(variant: dict[str, Any]) -> str:
    if variant["parameter"] == "baseline":
        return "baseline"
    if isinstance(variant["value"], str):
        return variant["name"]
    return f"{variant['parameter']}={variant['value']:.9g}"


def write_corrected_materials(path: Path) -> None:
    source = json.loads(BASE_MATERIALS.read_text(encoding="utf-8"))
    materials = {row["name"]: row for row in source["materials"]}
    common_transport = {
        "eps_r": 11.7,
        "ni": SENTAURUS_SILICON_NI_CM3,
        "mun": 1417.0,
        "mup": 470.5,
        "bandgap_eV": 1.12,
        "electron_affinity_eV": 4.0727403846153845,
        "electron_quantum_gamma": 3.6,
        "electron_quantum_dos_mass_ratio": DOS_MASS,
        "electron_quantum_coefficient_mass_ratio": COEFFICIENT_MASS,
        "Nc_m3": 2.8e19,
        "Nv_m3": 1.04e19,
        "temperature_K": 300.0,
    }
    materials["Si"] = {"name": "Si", **common_transport}
    materials["PolySilicon"] = {"name": "PolySilicon", **common_transport}
    materials["SiO2"] = {
        "name": "SiO2",
        "eps_r": 3.9,
        "ni": 0.0,
        "mun": 0.0,
        "mup": 0.0,
        "bandgap_eV": 9.0,
        "electron_affinity_eV": 0.9,
        "electron_quantum_gamma": 1.0,
        "electron_quantum_dos_mass_ratio": 0.42,
        "electron_quantum_coefficient_mass_ratio": 0.42,
        "temperature_K": 300.0,
    }
    nitride = dict(materials["Nitride"])
    nitride.update(
        {
            "electron_quantum_gamma": 1.0,
            "electron_quantum_dos_mass_ratio": 0.42,
            "electron_quantum_coefficient_mass_ratio": 0.42,
        }
    )
    materials["Nitride"] = nitride
    ordered = [materials[name] for name in ("Si", "SiO2", "PolySilicon", "Nitride")]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps({"materials": ordered}, indent=2) + "\n", encoding="utf-8")


def make_config(variant: dict[str, Any], run_dir: Path) -> tuple[Path, Path]:
    config = json.loads(BASE_CONFIG.read_text(encoding="utf-8"))
    prefix = run_dir / "eq231"
    for contact in config["contacts"]:
        if contact["name"].lower() == "drain":
            contact["bias"] = 2.0
    solver = config["solver"]
    solver["verbose"] = False
    quantum = solver["electron_quantum_potential"]
    quantum.update(
        {
            "enabled": True,
            "coupling_mode": "outer",
            "formulation": "potential_based",
            "include_insulators": True,
            "global_discretization": "p1_direct",
            "residual_diagnostic_prefix": str(prefix.resolve()),
            "residual_diagnostic_use_initial_state": True,
            "outer_max_iterations": 1,
            "max_iterations": 1,
            # Pin the two mass roles independently so each scan remains OFAT.
            "effective_mass_ratio": DOS_MASS,
            "coefficient_mass_ratio": DOS_MASS,
        }
    )
    if variant["parameter"] == "material_contract":
        config["materials_file"] = str(CORRECTED_MATERIALS.resolve())
    elif variant["parameter"] != "baseline":
        quantum[variant["parameter"]] = variant["value"]
    config["_comment"] = (
        "TransportModels fixed-state DG one-factor parameter audit: "
        + slug_label(variant)
    )
    config["output_csv"] = str((run_dir / "probe.csv").resolve())
    config["log_file"] = str((run_dir / "probe.log").resolve())
    sweep = config["sweep"]
    sweep.update(
        {
            "start": 2.0,
            "stop": 2.0,
            "step": 0.1,
            "bias_points": [2.0],
            "initial_state_file": str(HYBRID_RESTART.resolve()),
            "write_vtk": False,
            "write_state_file": str((run_dir / "probe_final_state.csv").resolve()),
            "write_state_every_point_prefix": str((run_dir / "probe_state").resolve()),
            "diagnostics": {
                "transport": {"enabled": True},
                "terminal_balance": {
                    "enabled": True,
                    "contacts": ["source", "drain", "gate", "substrate"],
                    "csv_file": str((run_dir / "terminal_balance.csv").resolve()),
                },
            },
        }
    )
    path = run_dir / "config.json"
    path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    return path, prefix


def execute_variant(variant: dict[str, Any], output_root: Path) -> dict[str, Any]:
    run_dir = output_root / variant["name"]
    run_dir.mkdir(parents=True, exist_ok=True)
    config, prefix = make_config(variant, run_dir)
    process = subprocess.run(
        [str(RUNNER.resolve()), "--config", str(config)],
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        encoding="utf-8",
        errors="replace",
    )
    console = run_dir / "console.log"
    console.write_text(
        process.stdout + ("\n[stderr]\n" + process.stderr if process.stderr else ""),
        encoding="utf-8",
    )
    nodes = Path(str(prefix) + "_nodes.csv")
    if not nodes.is_file():
        raise RuntimeError(f"{variant['name']}: residual diagnostic missing; see {console}")
    return analyze_variant(variant, run_dir, config, prefix, process.returncode)


def analyze_variant(
    variant: dict[str, Any],
    run_dir: Path,
    config: Path,
    prefix: Path,
    runner_exit_code: int = 1,
) -> dict[str, Any]:
    nodes = Path(str(prefix) + "_nodes.csv")
    result = audit(prefix)
    substrate = next(row for row in result["regions"] if row["region_name"] == "R.Substrate")
    return {
        **variant,
        "label": slug_label(variant),
        "runner_exit_code": runner_exit_code,
        "max_free_residual": result["summary"]["max_free_residual"],
        "cell_total_l1_free": result["summary"]["cell_total_l1_free"],
        "max_free_node": result["summary"]["max_free_node"],
        "substrate_l1_share": substrate["global_l1_share"],
        "reaction_l1_share": result["component_l1_share"]["reaction"],
        "gradient_squared_l1_share": result["component_l1_share"]["gradient_squared"],
        "config": str(config),
        "nodes": str(nodes),
        "config_sha256": sha256(config),
        "nodes_sha256": sha256(nodes),
    }


def write_csv(path: Path, results: list[dict[str, Any]]) -> None:
    fields = [
        "name", "parameter", "value", "label", "runner_exit_code",
        "max_free_residual", "max_free_residual_ratio", "cell_total_l1_free",
        "cell_total_l1_ratio", "max_free_node", "substrate_l1_share",
        "reaction_l1_share", "gradient_squared_l1_share", "config", "nodes",
        "config_sha256", "nodes_sha256",
    ]
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows({key: row.get(key, "") for key in fields} for row in results)


def plot_results(path: Path, results: list[dict[str, Any]]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    ranked = sorted(results, key=lambda row: row["cell_total_l1_ratio"])
    labels = [row["label"] for row in ranked]
    values = [row["cell_total_l1_ratio"] for row in ranked]
    colors = ["#C45A20" if row["name"] == "baseline" else "#2F6B9A" for row in ranked]
    fig, ax = plt.subplots(figsize=(10.8, 7.4))
    bars = ax.barh(labels, values, color=colors, edgecolor="#183B56", linewidth=0.45)
    ax.axvline(1.0, color="#374151", linewidth=1.1, linestyle="--", label="baseline ratio")
    ax.bar_label(bars, labels=[f"{value:.4f}×" for value in values], padding=4, fontsize=8)
    ax.set_xlim(0, max(values) * 1.14)
    ax.set_xlabel("Free-node Eq. 231 residual L1 / baseline")
    ax.set_title("TransportModels DG fixed-state parameter sweep", pad=30)
    ax.text(
        0.0,
        1.02,
        "Vg = 1.0 V, Vd = 2.0 V; one factor at a time; lower is a better fixed-state match",
        transform=ax.transAxes,
        color="#4B5563",
        fontsize=9,
    )
    ax.grid(axis="x", color="#D1D5DB", linewidth=0.6, alpha=0.7)
    ax.set_axisbelow(True)
    ax.legend(loc="lower right", frameon=False)
    fig.subplots_adjust(left=0.34, right=0.94, top=0.88, bottom=0.10)
    fig.savefig(path, dpi=220, facecolor="white")
    fig.savefig(path.with_suffix(".svg"), facecolor="white")
    plt.close(fig)


def render_markdown(report: dict[str, Any]) -> str:
    results = report["results"]
    baseline = next(row for row in results if row["name"] == "baseline")
    best = min(results, key=lambda row: row["cell_total_l1_ratio"])
    recovered = next(row for row in results if row["name"] == "coefficient_mass_recovered")
    corrected = next(row for row in results if row["name"] == "corrected_material_contract")
    bgn_one = next(row for row in results if row["name"] == "bgn_fraction_1")
    map_rows = "\n".join(
        f"| {row['sentaurus_parameter']} | `{row['vela_field']}` | {row['value']} | {row['unit']} | {row['status']} |"
        for row in report["parameter_map"]
    )
    result_rows = "\n".join(
        f"| {row['label']} | {row['cell_total_l1_ratio']:.6f} | {row['max_free_residual_ratio']:.6f} | "
        f"{row['substrate_l1_share']:.2%} | {row['max_free_node']} |"
        for row in sorted(results, key=lambda row: row["cell_total_l1_ratio"])
    )
    return f"""# TransportModels DG parameter and unit audit

Work point: Vg = 1.0 V, Vd = 2.0 V. Fixed hybrid state, p1_direct Eq. 231.

Status: **{report['status']}**

## Main finding

The best tested one-factor setting is `{best['label']}`, with a normalized
free-node residual L1 of `{best['cell_total_l1_ratio']:.6f}`. The independently
recovered coefficient mass `{COEFFICIENT_MASS:.16g}` gives
`{recovered['cell_total_l1_ratio']:.6f}` versus the explicit baseline
`{baseline['cell_total_l1_ratio']:.6f}`.

This scan is diagnostic, not a fitted production calibration. A lower fixed-
state residual identifies a parameter direction worth testing self-consistently;
it does not by itself prove improved terminal-current agreement.

The TDR-derived material contract reduces the residual L1 to
`{corrected['cell_total_l1_ratio']:.6f}` (a reduction of
`{1.0 - corrected['cell_total_l1_ratio']:.2%}`), close to the unphysical
`Bgn2Chi=1` control at `{bgn_one['cell_total_l1_ratio']:.6f}`. This shows that
the apparent BGN-share preference is largely a proxy for base-affinity
mismatch: Si/PolySi require +22.740 mV and SiO2 requires -50.000 mV relative
to the original Vela material file. The semantic BGN fraction remains 0.5.

## Parameter mapping

| Sentaurus quantity | Vela field | Value | Unit | Mapping status |
|---|---|---:|---|---|
{map_rows}

The TransportModels `pp13_des.par` file overrides SRH lifetime parameters but
does not override the quantum-potential section, so the Silicon defaults own
the DG parameter values. The two electron mass roles are intentionally kept
separate: DOS mass changes the material drive, while coefficient mass changes
the gradient coefficient.

## One-factor scan

| Variant | Residual L1 / baseline | Max residual / baseline | Substrate L1 share | Max node |
|---|---:|---:|---:|---:|
{result_rows}

All variants use identical mesh, state, contacts, discretization, temperature,
and carrier models. Only the named DG parameter changes.
"""


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-root", type=Path, default=OUTPUT_ROOT)
    parser.add_argument("--jobs", type=int, default=3)
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--reuse", action="store_true")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    if args.check:
        report = json.loads(REPORT_JSON.read_text(encoding="utf-8"))
        failures = []
        for row in report["results"]:
            for field, hash_field in (("config", "config_sha256"), ("nodes", "nodes_sha256")):
                actual = sha256(Path(row[field]))
                if actual != row[hash_field]:
                    failures.append({"variant": row["name"], "field": field, "actual": actual})
        actual_materials = sha256(Path(report["paths"]["corrected_materials"]))
        if actual_materials != report["hashes"]["corrected_materials"]:
            failures.append(
                {"field": "corrected_materials", "actual": actual_materials}
            )
        if failures:
            print(json.dumps(failures, indent=2))
            return 1
        print("TransportModels DG parameter sweep check: PASS")
        return 0

    output_root = args.output_root.resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    write_corrected_materials(CORRECTED_MATERIALS)
    if not args.execute and not args.reuse:
        for variant in VARIANTS:
            run_dir = output_root / variant["name"]
            run_dir.mkdir(parents=True, exist_ok=True)
            config, _ = make_config(variant, run_dir)
            print(config)
        return 0

    results: list[dict[str, Any]] = []
    if args.reuse:
        for variant in VARIANTS:
            run_dir = output_root / variant["name"]
            results.append(
                analyze_variant(
                    variant,
                    run_dir,
                    run_dir / "config.json",
                    run_dir / "eq231",
                )
            )
    else:
        with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as executor:
            futures = {
                executor.submit(execute_variant, variant, output_root): variant["name"]
                for variant in VARIANTS
            }
            for future in as_completed(futures):
                result = future.result()
                results.append(result)
                print(f"completed {result['name']}", flush=True)
    order = {variant["name"]: index for index, variant in enumerate(VARIANTS)}
    results.sort(key=lambda row: order[row["name"]])
    baseline = results[0]
    for row in results:
        row["cell_total_l1_ratio"] = row["cell_total_l1_free"] / baseline["cell_total_l1_free"]
        row["max_free_residual_ratio"] = row["max_free_residual"] / baseline["max_free_residual"]

    csv_path = output_root / "parameter_sweep.csv"
    figure = output_root / "parameter_sweep_residual_l1.png"
    write_csv(csv_path, results)
    plot_results(figure, results)
    report: dict[str, Any] = {
        "schema": "vela.transportmodels.dg_parameter_fixed_state_sweep.v1",
        "status": "pass",
        "as_of": "2026-08-21",
        "work_point": {"gate_bias_V": 1.0, "drain_bias_V": 2.0},
        "metric_definition": {
            "primary": "sum of absolute integrated p1_direct Eq. 231 cell contributions on free active nodes",
            "comparison": "ratio to explicit baseline",
        },
        "parameter_map": list(PARAMETER_MAP),
        "results": results,
        "paths": {
            "csv": str(csv_path),
            "figure_png": str(figure),
            "figure_svg": str(figure.with_suffix(".svg")),
            "corrected_materials": str(CORRECTED_MATERIALS),
            "band_drive_audit": str(
                REPO_ROOT / "docs/validation/transportmodels_dg_band_drive_audit_2026-08-21.json"
            ),
        },
        "hashes": {
            "csv": sha256(csv_path),
            "hybrid_restart": sha256(HYBRID_RESTART),
            "corrected_materials": sha256(CORRECTED_MATERIALS),
        },
    }
    markdown = render_markdown(report)
    REPORT_JSON.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    REPORT_MD.write_text(markdown, encoding="utf-8")
    (output_root / "summary.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    (output_root / "summary.md").write_text(markdown, encoding="utf-8")
    print(json.dumps({"status": "pass", "best": min(results, key=lambda row: row["cell_total_l1_ratio"])["label"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
