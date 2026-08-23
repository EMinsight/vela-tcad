#!/usr/bin/env python3
"""Fixed-state SRH model A/Bs for the TransportModels DD deep-off point."""

from __future__ import annotations

import csv
import importlib.util
import json
import math
import os
import subprocess
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[1]
AUDIT_SPEC = importlib.util.spec_from_file_location(
    "transportmodels_dd_deep_off_audit",
    REPO / "scripts/run_transportmodels_dd_deep_off_fixed_state_audit.py",
)
assert AUDIT_SPEC is not None and AUDIT_SPEC.loader is not None
audit = importlib.util.module_from_spec(AUDIT_SPEC)
AUDIT_SPEC.loader.exec_module(audit)

OUTPUT = audit.OUTPUT / "srh_model_ab"
STATE = audit.OUTPUT / "sentaurus_state_for_vela.csv"
FEEDBACK = audit.OUTPUT / "sentaurus_feedback_fields"
REPORT = REPO / "docs/validation/transportmodels_dd_deep_off_srh_ab_2026-08-23.md"
VT = 0.025851999786435535
SILICON_NI_CM3 = 1.0e10
SENTAURUS_SILICON_NI_CM3 = 1.4638914958767616e10


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def run_config(label: str, config: dict[str, Any]) -> dict[str, Any]:
    config_path = OUTPUT / f"{label}.json"
    config_path.write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    env = os.environ.copy()
    env["Path"] = r"D:\msys64\ucrt64\bin" + os.pathsep + env.get("Path", "")
    completed = subprocess.run(
        [str(audit.RUNNER), "--config", str(config_path),
         "--log", str(OUTPUT / f"{label}.log")],
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


def write_sentaurus_ni_materials() -> Path:
    base_config = json.loads((audit.VELA_ROOT / "config.json").read_text(encoding="utf-8"))
    source = Path(base_config["materials_file"])
    data = json.loads(source.read_text(encoding="utf-8"))
    for material in data["materials"]:
        if material["name"] in {"Si", "PolySilicon"}:
            material["ni"] = SENTAURUS_SILICON_NI_CM3
    output = OUTPUT / "materials_sentaurus_ni.json"
    output.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    return output


def configure_variant(config: dict[str, Any], statistics: str,
                      bandgap_narrowing: str | dict[str, Any],
                      materials_file: Path | None = None) -> None:
    config["solver"]["carrier_statistics"] = {"model": statistics}
    config["solver"]["bandgap_narrowing"] = bandgap_narrowing
    if materials_file is not None:
        config["materials_file"] = str(materials_file.resolve())


def run_variant(label: str, statistics: str,
                bandgap_narrowing: str | dict[str, Any],
                silicon_nodes: set[int], current_scale: float,
                materials_file: Path | None = None) -> dict[str, Any]:
    feedback_csv = OUTPUT / f"{label}_feedback.csv"
    feedback_cfg = audit.probe_config(
        "newton_feedback_substitution_probe",
        STATE,
        feedback_csv,
        srh_reference_internal=1.0e16,
    )
    feedback_cfg["feedback_state_fields_dir"] = str(FEEDBACK.resolve())
    configure_variant(feedback_cfg, statistics, bandgap_narrowing, materials_file)
    feedback_status = run_config(f"{label}_feedback", feedback_cfg)

    density_rows = [
        row for row in read_csv(feedback_csv)
        if row["variant"] == "density_only" and int(row["node_id"]) in silicon_nodes
    ]
    scaled_source = sum(float(row["electron_recombination"]) for row in density_rows)

    terms_csv = OUTPUT / f"{label}_terms.csv"
    terms_cfg = audit.probe_config(
        "newton_carrier_term_probe",
        STATE,
        terms_csv,
        srh_reference_internal=1.0e16,
    )
    configure_variant(terms_cfg, statistics, bandgap_narrowing, materials_file)
    terms_status = run_config(f"{label}_terms", terms_cfg)

    return {
        "statistics": statistics,
        "bandgap_narrowing": bandgap_narrowing,
        "materials_file": str(materials_file.resolve()) if materials_file else None,
        "fixed_exact_density_srh_A_per_um": scaled_source * current_scale,
        "feedback_status": feedback_status,
        "terms_status": terms_status,
        "feedback_csv": str(feedback_csv.resolve()),
        "terms_csv": str(terms_csv.resolve()),
    }


def weighted_metrics(values: list[tuple[float, float]]) -> dict[str, float]:
    total_weight = sum(weight for _, weight in values)
    if total_weight <= 0.0:
        return {"weighted_mean": math.nan, "weighted_rms": math.nan, "maximum": math.nan}
    return {
        "weighted_mean": sum(value * weight for value, weight in values) / total_weight,
        "weighted_rms": math.sqrt(
            sum(value * value * weight for value, weight in values) / total_weight
        ),
        "maximum": max(value for value, _ in values),
    }


def bgn_comparison(variant_terms: dict[str, dict[str, Any]]) -> dict[str, Any]:
    sent_bgn = audit.scalar_field("BandgapNarrowing", 3)
    sent_srh = audit.scalar_field("srhRecombination", 3)
    nodes = {
        int(row["id"]): (float(row["x_um"]), float(row["y_um"]))
        for row in read_csv(audit.SENTAURUS_ROOT / "nodes.csv")
    }
    control_area = {node: 0.0 for node in sent_bgn}
    for row in read_csv(audit.SENTAURUS_ROOT / "elements.csv"):
        if row["region"] != "R.Substrate":
            continue
        ids = (int(row["node0"]), int(row["node1"]), int(row["node2"]))
        area = audit.triangle_area_um2(nodes, ids)
        for node in ids:
            control_area[node] += area / 3.0

    baseline_rows = {
        int(row["node_id"]): row
        for row in read_csv(Path(variant_terms["fermi_old_slotboom"]["terms_csv"]))
    }
    corrected_rows = {
        int(row["node_id"]): row
        for row in read_csv(Path(
            variant_terms["fermi_old_slotboom_fermi_correction"]["terms_csv"]
        ))
    }
    output_rows: list[dict[str, Any]] = []
    baseline_errors: list[tuple[float, float]] = []
    corrected_errors: list[tuple[float, float]] = []
    for node in sorted(sent_bgn):
        baseline_ni = max(float(baseline_rows[node]["ni_eff_m3"]), SILICON_NI_CM3)
        corrected_ni = max(float(corrected_rows[node]["ni_eff_m3"]), SILICON_NI_CM3)
        baseline_delta = 2.0 * VT * math.log(baseline_ni / SILICON_NI_CM3)
        corrected_delta = 2.0 * VT * math.log(corrected_ni / SILICON_NI_CM3)
        weight = abs(sent_srh[node]) * control_area[node]
        baseline_error = abs(baseline_delta - sent_bgn[node])
        corrected_error = abs(corrected_delta - sent_bgn[node])
        baseline_errors.append((baseline_error, weight))
        corrected_errors.append((corrected_error, weight))
        output_rows.append({
            "node_id": node,
            "sentaurus_bgn_eV": sent_bgn[node],
            "vela_old_slotboom_bgn_eV": baseline_delta,
            "vela_fermi_corrected_bgn_eV": corrected_delta,
            "sentaurus_srh_cm3_s": sent_srh[node],
            "barycentric_control_area_um2": control_area[node],
            "srh_absolute_weight": weight,
        })
    audit.write_csv(OUTPUT / "bgn_node_comparison.csv", output_rows)
    return {
        "old_slotboom": weighted_metrics(baseline_errors),
        "old_slotboom_fermi_correction": weighted_metrics(corrected_errors),
        "node_csv": str((OUTPUT / "bgn_node_comparison.csv").resolve()),
    }


def main() -> int:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    parent_summary = json.loads((audit.OUTPUT / "summary.json").read_text(encoding="utf-8"))
    current_scale = float(
        parent_summary["srh_integral_A_per_um"]["carrier_term_to_A_per_um_scale"]
    )
    silicon_nodes = set(audit.scalar_field("srhRecombination", 3))
    sentaurus_ni_materials = write_sentaurus_ni_materials()
    variants_spec: dict[str, tuple[str, str | dict[str, Any], Path | None]] = {
        "fermi_old_slotboom": ("fermi_dirac", "old_slotboom", None),
        "boltzmann_old_slotboom": ("boltzmann", "old_slotboom", None),
        "fermi_no_bgn": ("fermi_dirac", "none", None),
        "boltzmann_no_bgn": ("boltzmann", "none", None),
        "fermi_old_slotboom_fermi_correction": (
            "fermi_dirac",
            {"model": "old_slotboom", "fermi_statistics_correction": True},
            None,
        ),
        "fermi_corrected_bgn_sentaurus_ni": (
            "fermi_dirac",
            {"model": "old_slotboom", "fermi_statistics_correction": True},
            sentaurus_ni_materials,
        ),
    }
    variants = {
        label: run_variant(
            label, statistics, bgn, silicon_nodes, current_scale, materials_file
        )
        for label, (statistics, bgn, materials_file) in variants_spec.items()
    }
    sentaurus_source = abs(float(
        parent_summary["srh_integral_A_per_um"]["sentaurus_exported_field"]
    ))
    for result in variants.values():
        result["magnitude_ratio_to_sentaurus"] = (
            abs(result["fixed_exact_density_srh_A_per_um"]) / sentaurus_source
        )

    baseline = abs(variants["fermi_old_slotboom"]["fixed_exact_density_srh_A_per_um"])
    boltzmann = abs(variants["boltzmann_old_slotboom"]["fixed_exact_density_srh_A_per_um"])
    no_bgn = abs(variants["fermi_no_bgn"]["fixed_exact_density_srh_A_per_um"])
    fermi_bgn = abs(
        variants["fermi_old_slotboom_fermi_correction"]
        ["fixed_exact_density_srh_A_per_um"]
    )
    sentaurus_ni = abs(
        variants["fermi_corrected_bgn_sentaurus_ni"]
        ["fixed_exact_density_srh_A_per_um"]
    )
    sentaurus_terminal = audit.SENTAURUS_DD_CURRENT_A_PER_UM
    sentaurus_field = sentaurus_source
    summary = {
        "schema": "vela.transportmodels_dd_deep_off_srh_ab.v1",
        "status": "complete",
        "bias": {"gate_V": -1.0, "drain_V": 1.1},
        "variants": variants,
        "effects": {
            "fermi_factor_relative_change_vs_boltzmann": baseline / boltzmann - 1.0,
            "old_slotboom_relative_change_vs_no_bgn": baseline / no_bgn - 1.0,
            "fermi_bgn_correction_relative_change_vs_old_slotboom": fermi_bgn / baseline - 1.0,
            "sentaurus_ni_relative_change_vs_1e10": sentaurus_ni / fermi_bgn - 1.0,
        },
        "bgn_field_comparison": bgn_comparison(variants),
        "sentaurus_silicon_intrinsic_density_cm3": SENTAURUS_SILICON_NI_CM3,
        "quadrature_closure": {
            "sentaurus_exported_srh_integral_A_per_um": sentaurus_field,
            "sentaurus_terminal_current_A_per_um": sentaurus_terminal,
            "relative_difference": abs(sentaurus_terminal - sentaurus_field)
            / abs(sentaurus_terminal),
            "interpretation": (
                "The barycentric integral of the exported Sentaurus SRH field "
                "already closes to the terminal current, bounding quadrature "
                "and non-SRH terminal contributions at this bias."
            ),
        },
    }
    summary_path = OUTPUT / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    lines = [
        "# TransportModels DD deep-off SRH model A/B",
        "",
        "Bias: `Vg=-1 V`, `Vd=1.1 V`; exact Sentaurus psi, quasi-Fermi potentials, n, and p are held fixed.",
        "",
        "| Variant | SRH integral (A/um) | Sentaurus ratio |",
        "|---|---:|---:|",
    ]
    for label, result in variants.items():
        lines.append(
            f"| {label} | {result['fixed_exact_density_srh_A_per_um']:.9e} | "
            f"{100.0 * result['magnitude_ratio_to_sentaurus']:.4f}% |"
        )
    lines.extend([
        "",
        "## Isolated effects",
        "",
        f"- Generalized Fermi factors relative to Boltzmann: {100.0 * summary['effects']['fermi_factor_relative_change_vs_boltzmann']:.4f}%.",
        f"- OldSlotboom relative to no BGN: {100.0 * summary['effects']['old_slotboom_relative_change_vs_no_bgn']:.4f}%.",
        f"- Fermi BGN correction relative to OldSlotboom: {100.0 * summary['effects']['fermi_bgn_correction_relative_change_vs_old_slotboom']:.4f}%.",
        f"- Sentaurus silicon ni relative to 1e10 cm^-3: {100.0 * summary['effects']['sentaurus_ni_relative_change_vs_1e10']:.4f}%.",
        f"- Sentaurus SRH-field integral versus terminal-current closure: {100.0 * summary['quadrature_closure']['relative_difference']:.4f}%.",
        "",
        f"Raw summary: `{summary_path}`",
    ])
    REPORT.parent.mkdir(parents=True, exist_ok=True)
    REPORT.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
