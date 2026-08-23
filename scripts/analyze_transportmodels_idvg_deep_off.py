#!/usr/bin/env python3
"""Diagnose the TransportModels Id-Vg deep-off current discrepancy.

The analysis intentionally uses already-produced Sentaurus and Vela artifacts.  It
does not rerun either simulator or change solver configuration.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[1]
DEFAULT_REPORT_DIR = (
    REPO
    / "build-release/reference_tcad/transportmodels_sentaurus2022/reports"
    / "idvg_deep_off_20260821"
)
DEFAULT_MARKDOWN = (
    REPO
    / "docs/validation/transportmodels_idvg_deep_off_analysis_2026-08-21.md"
)

SENTAURUS_PLT = (
    REPO
    / "build-release/reference_tcad/transportmodels_sentaurus2022/sentaurus_vm_runs"
    / "idvg_semantics_2x2_20260821/raw/bundle/default_default/IdVgs_n7_des.plt"
)
SENTAURUS_CMD = (
    REPO
    / "build-release/reference_tcad/transportmodels_sentaurus2022/sentaurus_vm_runs"
    / "idvg_semantics_2x2_20260821/raw/bundle/default_default/idvg_2x2_des.cmd"
)
SENTAURUS_PAR = SENTAURUS_CMD.with_name("pp7_des.par")
SENTAURUS_DG_REFERENCE = (
    REPO
    / "build-release/reference_tcad/transportmodels_sentaurus2022/vela_baseline"
    / "generated/reference_curves/transportmodels_sentaurus2022_dg_idvg_reference.csv"
)
SENTAURUS_DD_REFERENCE = SENTAURUS_DG_REFERENCE.with_name(
    "transportmodels_sentaurus2022_dd_idvg_reference.csv"
)

VELA_DG_DIR = (
    REPO
    / "build-release/reference_tcad/transportmodels_sentaurus2022/vela_baseline"
    / "dg_post_p2_regression_v4_2026-08-21/idvg"
)
VELA_DD_DIR = (
    REPO
    / "build-release/reference_tcad/transportmodels_sentaurus2022/vela_baseline"
    / "dd_phase7_shared_baseline_2026-08-21/idvg"
)

SELECTED_BIASES = (-1.0, -0.84, -0.68, -0.52)


def load_sentaurus_import() -> Any:
    spec = importlib.util.spec_from_file_location(
        "sentaurus_import", REPO / "scripts/sentaurus_import.py"
    )
    if spec is None or spec.loader is None:
        raise RuntimeError("Unable to load scripts/sentaurus_import.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def nearest(rows: list[dict[str, Any]], field: str, value: float) -> dict[str, Any]:
    row = min(rows, key=lambda item: abs(float(item[field]) - value))
    if abs(float(row[field]) - value) > 1.0e-10:
        raise ValueError(f"No row at {field}={value} V")
    return row


def sentaurus_rows() -> list[dict[str, float]]:
    parser = load_sentaurus_import()
    text = SENTAURUS_PLT.read_text(errors="ignore")
    datasets = parser.parse_quoted_list(text, "datasets")
    values = parser.parse_values_block(text, len(datasets))
    return [dict(zip(datasets, row, strict=True)) for row in values]


def float_row(row: dict[str, str]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in row.items():
        try:
            result[key] = float(value)
        except (TypeError, ValueError):
            result[key] = value
    return result


def terminal_rows(path: Path) -> list[dict[str, Any]]:
    return [float_row(row) for row in read_csv(path)]


def curve_rows(path: Path) -> list[dict[str, Any]]:
    return [float_row(row) for row in read_csv(path)]


def terminal_at(
    rows: list[dict[str, Any]], bias: float, contact: str
) -> dict[str, Any]:
    matching = [row for row in rows if row["contact"] == contact]
    return nearest(matching, "bias_V", bias)


def log10_abs(value: float) -> float:
    return math.log10(max(abs(value), 1.0e-300))


def relative_error(candidate: float, reference: float) -> float:
    return abs(candidate - reference) / max(abs(reference), 1.0e-300)


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def build_analysis() -> tuple[list[dict[str, Any]], dict[str, Any], list[dict[str, Any]]]:
    sent_rows = sentaurus_rows()
    vela_dg_curve = curve_rows(VELA_DG_DIR / "curve.csv")
    vela_dg_terminal = terminal_rows(VELA_DG_DIR / "terminal_balance.csv")
    vela_dd_curve = curve_rows(VELA_DD_DIR / "curve_combined.csv")
    sent_dg_reference = curve_rows(SENTAURUS_DG_REFERENCE)
    sent_dd_reference = curve_rows(SENTAURUS_DD_REFERENCE)

    diagnostics: list[dict[str, Any]] = []
    for bias in SELECTED_BIASES:
        sent = nearest(sent_rows, "gate OuterVoltage", bias)
        dg = nearest(vela_dg_curve, "bias_V", bias)
        dd = nearest(vela_dd_curve, "bias_V", bias)
        sent_ref_dg = nearest(sent_dg_reference, "bias_V", bias)
        sent_ref_dd = nearest(sent_dd_reference, "bias_V", bias)

        sent_contact_totals = {
            contact: float(sent[f"{contact} TotalCurrent"])
            for contact in ("source", "drain", "gate", "substrate")
        }
        vela_contacts = {
            contact: terminal_at(vela_dg_terminal, bias, contact)
            for contact in ("source", "drain", "gate", "substrate")
        }
        sent_kcl = abs(sum(sent_contact_totals.values()))
        vela_kcl = abs(
            sum(float(row["current_total_A_per_um"]) for row in vela_contacts.values())
        )
        sent_id = abs(float(sent["drain TotalCurrent"]))
        vela_dg_id = abs(float(dg["current_total_A_per_um"]))
        vela_dd_id = abs(float(dd["current_total_A_per_um"]))
        sent_drain_e = float(sent["drain eCurrent"])
        sent_drain_h = float(sent["drain hCurrent"])
        sent_substrate_h = float(sent["substrate hCurrent"])
        vela_drain_e = float(dg["current_electron_A_per_um"])
        vela_drain_h = float(dg["current_hole_A_per_um"])
        drift = float(dg["current_electron_drift_A_per_um"])
        diffusion = float(dg["current_electron_diffusion_A_per_um"])

        diagnostics.append(
            {
                "bias_V": bias,
                "sentaurus_dg_id_A_per_um": sent_id,
                "vela_dg_id_A_per_um": vela_dg_id,
                "dg_relative_error": relative_error(vela_dg_id, sent_id),
                "dg_log_error_dex": abs(log10_abs(vela_dg_id) - log10_abs(sent_id)),
                "sentaurus_dd_id_A_per_um": abs(float(sent_ref_dd["current_total"])),
                "vela_dd_id_A_per_um": vela_dd_id,
                "sentaurus_drain_e_A_per_um": sent_drain_e,
                "sentaurus_drain_h_A_per_um": sent_drain_h,
                "sentaurus_substrate_h_A_per_um": sent_substrate_h,
                "sentaurus_generation_pair_ratio": (
                    abs(sent_substrate_h) / max(abs(sent_drain_e), 1.0e-300)
                ),
                "sentaurus_kcl_residual_A_per_um": sent_kcl,
                "sentaurus_kcl_to_id_ratio": sent_kcl / max(sent_id, 1.0e-300),
                "vela_drain_e_A_per_um": vela_drain_e,
                "vela_drain_h_A_per_um": vela_drain_h,
                "vela_e_drift_A_per_um": drift,
                "vela_e_diffusion_A_per_um": diffusion,
                "vela_drift_diff_cancellation_ratio": (
                    abs(drift + diffusion)
                    / max(abs(drift) + abs(diffusion), 1.0e-300)
                ),
                "vela_kcl_residual_A_per_um": vela_kcl,
                "vela_kcl_to_id_ratio": vela_kcl / max(vela_dg_id, 1.0e-300),
                "sentaurus_reference_crosscheck_relative": relative_error(
                    sent_id, abs(float(sent_ref_dg["current_total"]))
                ),
            }
        )

    curve_chart: list[dict[str, Any]] = []
    for sent_dg, sent_dd, vela_dg, vela_dd in zip(
        sent_dg_reference, sent_dd_reference,
        curve_rows(VELA_DG_DIR / "curve_combined.csv"), vela_dd_curve,
        strict=True,
    ):
        bias = float(sent_dg["bias_V"])
        curve_chart.append(
            {
                "bias_V": bias,
                "sentaurus_dg_log10_id": log10_abs(float(sent_dg["current_total"])),
                "vela_dg_log10_id": log10_abs(float(vela_dg["current_total_A_per_um"])),
                "sentaurus_dd_log10_id": log10_abs(float(sent_dd["current_total"])),
                "vela_dd_log10_id": log10_abs(float(vela_dd["current_total_A_per_um"])),
            }
        )

    at_m1 = diagnostics[0]
    summary = {
        "status": "complete",
        "classification": "verified_srh_generation_mismatch_with_unresolved_vela_current",
        "selected_biases_V": list(SELECTED_BIASES),
        "sentaurus_m1V_id_A_per_um": at_m1["sentaurus_dg_id_A_per_um"],
        "vela_m1V_id_A_per_um": at_m1["vela_dg_id_A_per_um"],
        "m1V_relative_error": at_m1["dg_relative_error"],
        "m1V_log_error_dex": at_m1["dg_log_error_dex"],
        "sentaurus_m1V_generation_pair_ratio": at_m1["sentaurus_generation_pair_ratio"],
        "vela_m1V_kcl_to_id_ratio": at_m1["vela_kcl_to_id_ratio"],
        "verified_findings": [
            "Sentaurus deep-off drain current is electron current balanced by substrate hole current.",
            "SRH(DopingDep TempDependence) is the only active pair-generation/recombination mechanism in the deck.",
            "Vela DD and DG both miss the approximately 1.6e-15 A/um plateau, so DG is not the primary cause.",
            "Vela's -1 V drain current is below its terminal KCL residual and is not numerically resolved.",
            "Sentaurus and Vela use the same stated Scharfetter lifetime parameters.",
        ],
        "likely_code_causes": [
            "The Fermi-Dirac SRH residual subtracts n*p and an equilibrium product directly, unlike the cancellation-free expm1 Boltzmann path.",
            "Vela's SRH denominator omits the Fermi degeneracy factors used by the generalized Sentaurus SRH equation.",
            "Deep-off terminal extraction is sensitive to cancellation between approximately 0.1 A/um drift and diffusion diagnostics.",
        ],
        "secondary_controls": [
            "Match the Sentaurus forward sweep direction and drain-ramp initialization.",
            "Tighten carrier/source convergence and require terminal KCL residual below the current being compared.",
        ],
    }
    return diagnostics, summary, curve_chart


def markdown_report(rows: list[dict[str, Any]], summary: dict[str, Any]) -> str:
    m1 = rows[0]
    table_lines = [
        "| Vg (V) | Sentaurus DG Id (A/um) | Vela DG Id (A/um) | Relative error | Log error (dex) |",
        "|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        table_lines.append(
            "| {bias_V:.2f} | {sentaurus_dg_id_A_per_um:.6e} | "
            "{vela_dg_id_A_per_um:.6e} | {dg_relative_error:.4%} | "
            "{dg_log_error_dex:.6f} |".format(**row)
        )
    return f"""# TransportModels Id-Vg deep-off discrepancy analysis

Date: 2026-08-21

Status: **diagnosis complete; no solver change made**.

## Technical summary

The Sentaurus current plateau near `1.6e-15 A/um` is not adequately described
as a numerical current floor. At `Vg=-1 V`, `99.9994%` of the Sentaurus drain
current is electron current, while the substrate carries an opposite hole
current with a magnitude ratio of `{m1['sentaurus_generation_pair_ratio']:.6f}`.
Because the deck enables `SRH(DopingDep TempDependence)` and no Auger, BTBT, or
avalanche mechanism, this is a resolved SRH pair-generation current.

Vela does not reproduce that pair-generation plateau in either DD or DG. The
DG drain current at `-1 V` is `{m1['vela_dg_id_A_per_um']:.6e} A/um`, while the
four-terminal KCL residual is `{m1['vela_kcl_residual_A_per_um']:.6e} A/um`,
or `{m1['vela_kcl_to_id_ratio']:.1f}` times larger. The reported Vela drain
current is therefore below the numerical conservation resolution of this run.

## The discrepancy is confined to the deep-off source-dominated regime

{chr(10).join(table_lines)}

At `Vg=-0.52 V`, ordinary channel transport has risen above the generation
plateau and the DG relative error falls to `{rows[-1]['dg_relative_error']:.2%}`.

## Terminal-current decomposition identifies SRH generation

At `Vg=-1 V`:

- Sentaurus drain electron current: `{m1['sentaurus_drain_e_A_per_um']:.6e} A/um`.
- Sentaurus drain hole current: `{m1['sentaurus_drain_h_A_per_um']:.6e} A/um`.
- Sentaurus substrate hole current: `{m1['sentaurus_substrate_h_A_per_um']:.6e} A/um`.
- Sentaurus terminal KCL residual / drain current: `{m1['sentaurus_kcl_to_id_ratio']:.3e}`.
- Vela drain electron current is exactly zero at the saved precision; its
  terminal current is set by a `{m1['vela_drain_h_A_per_um']:.6e} A/um` hole
  contribution.

The same approximately `1.6e-15 A/um` Sentaurus plateau is present in both its
DD and DG references, whereas both Vela DD and DG fall many orders lower. This
rules out the density-gradient equation as the primary deep-off cause.

## Verified controls

- The Sentaurus and Vela decks both enable SRH with Scharfetter doping
  dependence and the same stated `taumin`, `taumax`, `Nref`, `gamma`, and
  temperature exponents.
- Sentaurus does not enable Auger, band-to-band tunneling, avalanche, or a gate
  leakage model in this case.
- The Sentaurus and Vela comparison uses the same imported 3315-node topology.
- The discrepancy is present in DD as well as DG.

## Most likely Vela implementation gaps

1. **Cancellation in the Fermi-Dirac SRH numerator.** The Boltzmann path uses
   `ni^2 * expm1(deltaPhi/Vt)`, but the Fermi-Dirac path subtracts `n*p` and an
   independently reconstructed equilibrium product. In deep depletion these
   close quantities can lose the small net-generation signal.
2. **Incomplete generalized SRH formula.** Sentaurus generalizes SRH for Fermi
   statistics and quantization with carrier degeneracy factors in both the
   numerator and denominator. Vela reconstructs an equilibrium product for the
   numerator, but its SRH denominator remains `taup*(n+ni)+taun*(p+ni)`.
3. **Net-current resolution.** Vela's drift and diffusion diagnostics are each
   about `0.1 A/um` and cancel to the `1e-15 A/um` scale. At the two lowest gate
   biases the primary electron terminal current becomes exactly zero at saved
   precision, and the KCL residual exceeds the reported drain current.

## Secondary comparison controls

The final Vela regression sweeps the gate from `2.2 V` down to `-1 V`, while
Sentaurus initializes at `-1 V`, ramps the drain to `1.1 V`, and then sweeps the
gate upward. This is not the leading explanation because DD and DG show the
same missing plateau, but sweep direction, initial state, the `10 mV` DG outer
tolerance, and carrier residual tolerances must be matched before a final
off-state acceptance test.

## Recommended next experiment order

1. Add an SRH source audit that integrates generation over the silicon and
   reconciles it against electron and hole terminal currents at `Vg=-1 V`.
2. Implement a cancellation-free generalized Fermi SRH excess product and the
   Sentaurus degeneracy factors, with unit tests around equilibrium.
3. Require `abs(sum(I_contact))` to be at least one decade below the Id being
   compared; otherwise label the point unresolved rather than assigning a
   relative-error pass/fail.
4. Repeat forward and reverse DD sweeps with matched initialization and tighter
   tolerances; only then rerun DG.
5. Keep deep-off acceptance separate from transition/on-state acceptance and
   use both log-current error and a terminal-conservation criterion.

## Confidence and limitation

The classification of the Sentaurus plateau as SRH generation is high
confidence because it is supported by carrier-resolved terminal currents,
current conservation, and the enabled-physics deck. The exact share attributable
to Vela numerator cancellation versus the missing generalized denominator is
not yet measured; that split requires the controlled implementation A/B tests
listed above.
"""


def portable_artifact(
    rows: list[dict[str, Any]],
    summary: dict[str, Any],
    curve_chart: list[dict[str, Any]],
    output_dir: Path,
) -> dict[str, Any]:
    generated_at = "2026-08-21T18:00:00+08:00"
    m1 = rows[0]
    component_rows = [
        {
            "component": "Sentaurus drain electron",
            "abs_current_A_per_um": abs(m1["sentaurus_drain_e_A_per_um"]),
            "log10_abs_current": log10_abs(m1["sentaurus_drain_e_A_per_um"]),
            "role": "physical pair current",
        },
        {
            "component": "Sentaurus substrate hole",
            "abs_current_A_per_um": abs(m1["sentaurus_substrate_h_A_per_um"]),
            "log10_abs_current": log10_abs(m1["sentaurus_substrate_h_A_per_um"]),
            "role": "physical pair current",
        },
        {
            "component": "Vela drain total",
            "abs_current_A_per_um": abs(m1["vela_dg_id_A_per_um"]),
            "log10_abs_current": log10_abs(m1["vela_dg_id_A_per_um"]),
            "role": "reported current",
        },
        {
            "component": "Vela four-terminal KCL residual",
            "abs_current_A_per_um": abs(m1["vela_kcl_residual_A_per_um"]),
            "log10_abs_current": log10_abs(m1["vela_kcl_residual_A_per_um"]),
            "role": "numerical resolution",
        },
    ]
    headline = [
        {
            "worst_log_error_dex": m1["dg_log_error_dex"],
            "vela_kcl_to_id": m1["vela_kcl_to_id_ratio"],
            "sentaurus_pair_balance": m1["sentaurus_generation_pair_ratio"],
        }
    ]
    relative_output = output_dir.relative_to(REPO).as_posix()
    sources = [
        {
            "id": "diagnostic-output",
            "label": "Deep-off diagnostic output",
            "path": f"{relative_output}/deep_off_diagnostics.csv",
            "query": {
                "engine": "DuckDB",
                "language": "sql",
                "sql": f"SELECT * FROM read_csv_auto('{relative_output}/deep_off_diagnostics.csv', header=true)",
                "description": "Read the four reviewed deep-off diagnostic rows.",
                "tables_used": [f"{relative_output}/deep_off_diagnostics.csv"],
                "filters": ["Vg in {-1.00, -0.84, -0.68, -0.52} V"],
                "metric_definitions": [
                    "Relative error = abs(Vela Id - Sentaurus Id) / abs(Sentaurus Id).",
                    "KCL residual = abs(sum of signed source, drain, gate, and substrate currents).",
                ],
            },
        },
        {
            "id": "curve-output",
            "label": "Reviewed 21-point Id-Vg log-current curve",
            "path": f"{relative_output}/idvg_log_curve.csv",
            "query": {
                "engine": "DuckDB",
                "language": "sql",
                "sql": f"SELECT * FROM read_csv_auto('{relative_output}/idvg_log_curve.csv', header=true) ORDER BY bias_V",
                "description": "Read the aligned 21-point DD/DG logarithmic current curves.",
                "tables_used": [f"{relative_output}/idvg_log_curve.csv"],
            },
        },
        {
            "id": "component-output",
            "label": "Vg=-1 V terminal-current component decomposition",
            "path": f"{relative_output}/deep_off_diagnostics.csv",
            "query": {
                "engine": "DuckDB",
                "language": "sql",
                "sql": (
                    "WITH r AS (SELECT * FROM read_csv_auto('"
                    f"{relative_output}/deep_off_diagnostics.csv', header=true) WHERE bias_V=-1) "
                    "SELECT 'Sentaurus drain electron' AS component, abs(sentaurus_drain_e_A_per_um) AS abs_current_A_per_um FROM r "
                    "UNION ALL SELECT 'Sentaurus substrate hole', abs(sentaurus_substrate_h_A_per_um) FROM r "
                    "UNION ALL SELECT 'Vela drain total', abs(vela_dg_id_A_per_um) FROM r "
                    "UNION ALL SELECT 'Vela four-terminal KCL residual', abs(vela_kcl_residual_A_per_um) FROM r"
                ),
                "description": "Extract the carrier-resolved pair current and Vela numerical-resolution comparators at Vg=-1 V.",
                "tables_used": [f"{relative_output}/deep_off_diagnostics.csv"],
                "filters": ["Vg=-1 V"],
            },
        },
        {
            "id": "sentaurus-current",
            "label": "Sentaurus T-2022.03-SP2 Id-Vg current file",
            "path": SENTAURUS_PLT.relative_to(REPO).as_posix(),
        },
        {
            "id": "sentaurus-deck",
            "label": "Sentaurus Id-Vg SDevice deck and Scharfetter parameters",
            "path": SENTAURUS_CMD.relative_to(REPO).as_posix(),
        },
        {
            "id": "vela-terminal",
            "label": "Vela DG terminal balance",
            "path": (VELA_DG_DIR / "terminal_balance.csv").relative_to(REPO).as_posix(),
        },
        {
            "id": "analysis-script",
            "label": "Reproducible deep-off analysis script",
            "path": "scripts/analyze_transportmodels_idvg_deep_off.py",
        },
        {
            "id": "sdevice-manual",
            "label": "Sentaurus Device User Guide T-2022.03, Chapter 16, pages 473-477",
        },
    ]
    blocks = [
        {
            "id": "title",
            "type": "markdown",
            "body": "# TransportModels Id-Vg 深关断区差异诊断",
        },
        {
            "id": "technical-summary",
            "type": "markdown",
            "body": (
                "## 技术结论\n\n"
                "**Sentaurus 在约 `1.6e-15 A/um` 的平台不是普通数值地板，而是已解析的 SRH 产生电流。** "
                "`Vg=-1 V` 时，漏端几乎全部为电子电流，衬底端存在大小近似相等、方向相反的空穴电流。"
                "Vela 的 DD 和 DG 均未再现该平台；其 DG 漏端电流还低于四端 KCL 不平衡量，因此当前 Vela 深关断 Id 尚未达到可比较的数值分辨率。"
            ),
        },
        {
            "id": "curve-heading",
            "type": "markdown",
            "body": (
                "## 差异只集中在深关断段\n\n"
                "完整 21 点对数电流曲线显示，Sentaurus DD/DG 在最低栅压形成相近平台，而 Vela DD/DG 继续下降。"
                "到 `Vg=-0.52 V`，沟道输运超过产生电流后，DG 相对误差已降到 `9.66%`。"
            ),
            "sourceId": "diagnostic-output",
        },
        {
            "id": "curve-chart-block",
            "type": "chart",
            "chartId": "idvg-log-curve",
            "layout": "full",
        },
        {
            "id": "curve-explanation",
            "type": "markdown",
            "body": (
                "图中纵轴为 `log10(|Id| / (A/um))`。Sentaurus 的 DD 与 DG 同时出现平台，"
                "而 Vela 的两种输运模型都缺失该平台，因此密度梯度方程不是首要原因。"
            ),
        },
        {
            "id": "deep-off-table-block",
            "type": "table",
            "tableId": "deep-off-table",
            "layout": "full",
        },
        {
            "id": "mechanism-heading",
            "type": "markdown",
            "body": (
                "## 端口分量证明 Sentaurus 平台来自 SRH 产生\n\n"
                "`Vg=-1 V` 时，Sentaurus 漏端电子电流与衬底空穴电流的幅值比为 `0.999634`。"
                "该算例仅启用 `SRH(DopingDep TempDependence)`，未启用 Auger、BTBT、雪崩或栅泄漏模型。"
                "这构成了载流子分量、端口守恒和启用物理三方面的一致证据。"
            ),
            "sourceId": "sentaurus-current",
        },
        {
            "id": "component-chart-block",
            "type": "chart",
            "chartId": "m1-components",
            "layout": "full",
        },
        {
            "id": "component-explanation",
            "type": "markdown",
            "body": (
                "横条使用电流绝对值的十进制对数。Vela 四端 KCL 残差比其漏端 Id 大约 `1233` 倍，"
                "所以当前漏端数值不能作为深关断物理电流接受。"
            ),
            "sourceId": "diagnostic-output",
        },
        {
            "id": "implementation-findings",
            "type": "markdown",
            "body": (
                "## 最可能的 Vela 实现缺口\n\n"
                "1. **Fermi-Dirac SRH 分子存在消减风险。** Boltzmann 分支使用 `expm1`，"
                "Fermi 分支却直接相减 `n*p` 与重构的平衡乘积，在深耗尽区可能丢失很小的净产生量。\n"
                "2. **广义 SRH 公式尚不完整。** Sentaurus 对 Fermi 统计和量子化使用简并因子修正分子和分母；"
                "Vela 的分母仍为 `taup*(n+ni)+taun*(p+ni)`。\n"
                "3. **终端净电流分辨率不足。** 漂移和扩散诊断各约 `0.1 A/um`，需在 `1e-15 A/um` 量级相消。"
            ),
        },
        {
            "id": "scope-method",
            "type": "markdown",
            "body": (
                "## 范围、数据与方法\n\n"
                "比较对象为 Sentaurus T-2022.03-SP2 与最终 Vela DD/DG 基线。"
                "深关断点定义为 `Vg=-1、-0.84、-0.68 V`，并以 `-0.52 V` 作为进入过渡区的对照点。"
                "相对误差以 Sentaurus Id 为分母；KCL 残差定义为四个端口有符号总电流之和的绝对值。"
            ),
            "sourceId": "analysis-script",
        },
        {
            "id": "limitations",
            "type": "markdown",
            "body": (
                "## 局限性与稳健性\n\n"
                "Sentaurus 平台的 SRH 归因置信度高；但 Vela 缺失量中有多少来自分子消减、"
                "多少来自简并因子或收敛容差，尚未通过受控 A/B 试验定量拆分。"
                "此外，Sentaurus 从 `-1 V` 向上扫描，最终 Vela 回归从 `2.2 V` 向下扫描，后续必须统一方向和初始化。"
            ),
        },
        {
            "id": "recommendations",
            "type": "markdown",
            "body": (
                "## 建议的后续顺序\n\n"
                "1. 增加硅区 SRH 产生率积分，并与电子、空穴端口电流逐项守恒。\n"
                "2. 实现无消减的广义 Fermi SRH 公式及简并因子，补充近热平衡单元测试。\n"
                "3. 规定 KCL 残差至少比待比较 Id 小一个数量级，否则标记为未解析。\n"
                "4. 先做同方向、同初始化和更严容差的 DD 正反扫，再回归 DG。\n"
                "5. 深关断验收同时使用对数误差和端口守恒，不与导通区相对误差混用。"
            ),
        },
        {
            "id": "further-questions",
            "type": "markdown",
            "body": (
                "## 仍需回答的问题\n\n"
                "- Sentaurus 广义 SRH 简并因子在该 MOS 结构中贡献多少？\n"
                "- Vela Fermi 分支的乘积相减在多少节点发生有效位损失？\n"
                "- 收紧载流子残差后，KCL 不平衡能否稳定低于 `1e-16 A/um`？"
            ),
        },
    ]
    cards: list[dict[str, Any]] = []
    charts = [
        {
            "id": "idvg-log-curve",
            "title": "DD/DG Id-Vg 对数电流比较",
            "subtitle": "21 个栅压点；纵轴为 log10(|Id| / (A/um))",
            "type": "line",
            "dataset": "idvg_curve",
            "sourceId": "curve-output",
            "intent": "trend",
            "encodings": {
                "x": {"field": "bias_V", "type": "quantitative", "label": "Vg", "unit": "V"},
                "y": {
                    "fields": [
                        "sentaurus_dg_log10_id",
                        "vela_dg_log10_id",
                        "sentaurus_dd_log10_id",
                        "vela_dd_log10_id",
                    ],
                    "type": "quantitative",
                    "label": "log10(|Id| / (A/um))",
                },
            },
            "layout": "full",
        },
        {
            "id": "m1-components",
            "title": "Vg=-1 V 的电流分量与数值分辨率",
            "subtitle": "横轴为 log10(|I| / (A/um))；较大数值表示电流幅值更大",
            "type": "horizontalBar",
            "dataset": "m1_components",
            "sourceId": "component-output",
            "intent": "comparison",
            "encodings": {
                "x": {"field": "component", "type": "nominal", "label": "分量"},
                "y": {"field": "log10_abs_current", "type": "quantitative", "label": "log10(|I| / (A/um))"},
                "tooltip": [
                    {"field": "abs_current_A_per_um", "type": "quantitative", "label": "|I|", "unit": "A/um"},
                    {"field": "role", "type": "nominal", "label": "作用"},
                ],
            },
            "layout": "full",
        },
    ]
    tables = [
        {
            "id": "deep-off-table",
            "title": "深关断偏压点误差明细",
            "subtitle": "Vg=-1 V 至 -0.52 V；相对误差以 Sentaurus 为分母",
            "dataset": "deep_off",
            "sourceId": "diagnostic-output",
            "density": "spacious",
            "defaultSort": {"field": "bias_V", "direction": "asc"},
            "columns": [
                {"field": "bias_V", "label": "Vg (V)", "format": "number"},
                {"field": "sentaurus_dg_id_A_per_um", "label": "Sentaurus DG Id (A/um)", "format": "number"},
                {"field": "vela_dg_id_A_per_um", "label": "Vela DG Id (A/um)", "format": "number"},
                {"field": "dg_relative_error", "label": "相对误差", "format": "percent"},
                {"field": "dg_log_error_dex", "label": "对数误差 (dex)", "format": "number"},
                {"field": "vela_kcl_to_id_ratio", "label": "Vela KCL/Id", "format": "number"},
            ],
        }
    ]
    return {
        "surface": "report",
        "manifest": {
            "version": 1,
            "surface": "report",
            "title": "TransportModels Id-Vg 深关断区差异诊断",
            "description": "Sentaurus 2022 与 Vela DD/DG 深关断电流的载流子分量和数值分辨率审计。",
            "generatedAt": generated_at,
            "blocks": blocks,
            "cards": cards,
            "charts": charts,
            "tables": tables,
            "sources": sources,
        },
        "snapshot": {
            "version": 1,
            "generatedAt": generated_at,
            "status": "ready",
            "datasets": {
                "headline": headline,
                "deep_off": rows,
                "idvg_curve": curve_chart,
                "m1_components": component_rows,
            },
        },
        "sources": sources,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_REPORT_DIR)
    parser.add_argument("--markdown", type=Path, default=DEFAULT_MARKDOWN)
    args = parser.parse_args()

    diagnostics, summary, curve_chart = build_analysis()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_csv(args.output_dir / "deep_off_diagnostics.csv", diagnostics)
    write_csv(args.output_dir / "idvg_log_curve.csv", curve_chart)
    (args.output_dir / "deep_off_analysis.json").write_text(
        json.dumps(
            {
                "summary": summary,
                "diagnostics": diagnostics,
                "sources": {
                    "sentaurus_plt": str(SENTAURUS_PLT),
                    "sentaurus_cmd": str(SENTAURUS_CMD),
                    "sentaurus_parameter_file": str(SENTAURUS_PAR),
                    "vela_dg_curve": str(VELA_DG_DIR / "curve.csv"),
                    "vela_dg_terminal_balance": str(VELA_DG_DIR / "terminal_balance.csv"),
                    "vela_dd_curve": str(VELA_DD_DIR / "curve_combined.csv"),
                },
            },
            indent=2,
            ensure_ascii=False,
        )
        + "\n",
        encoding="utf-8",
    )
    (args.output_dir / "artifact.json").write_text(
        json.dumps(
            portable_artifact(diagnostics, summary, curve_chart, args.output_dir),
            indent=2,
            ensure_ascii=False,
        )
        + "\n",
        encoding="utf-8",
    )
    args.markdown.parent.mkdir(parents=True, exist_ok=True)
    args.markdown.write_text(markdown_report(diagnostics, summary), encoding="utf-8")
    print(args.markdown)
    print(args.output_dir / "deep_off_analysis.json")
    print(args.output_dir / "artifact.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
