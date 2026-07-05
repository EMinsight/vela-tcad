#!/usr/bin/env python3
"""Rebuild the PN2D coarse BV baseline after contact/handoff fixes."""

from __future__ import annotations

import csv
import json
import math
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any


REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "scripts"))

import run_pn2d_coarse7x3_previous_full20_compare as compare  # noqa: E402


OUT_DIR = (
    REPO
    / "build-release"
    / "reference_tcad"
    / "pn2d_sentaurus2018_coarse7x3"
    / "reports"
    / "rebaseline_post_contact_handoff_fix_20260705"
)
IMPORT_DIR = (
    REPO
    / "build-release"
    / "reference_tcad"
    / "pn2d_sentaurus2018_coarse7x3"
    / "imported_reference"
)
LEGACY_SENT_EXPORT = (
    REPO
    / "build-release"
    / "reference_tcad"
    / "pn2d_sentaurus2018_coarse7x3"
    / "reports"
    / "coarse_vm_vector_compare"
    / "sentaurus_multibias"
)
RUNNER = REPO / "build-release" / ("vela_example_runner.exe" if os.name == "nt" else "vela_example_runner")
ANCHOR_BIASES = [0.0, -1.0, -5.0, -10.0, -18.0, -20.0]
NODE_COMPARE_BIASES = [-1.0, -5.0, -10.0, -18.0, -20.0]
OLD_VELA_ALPHA_FLUX_OVER_SENT_FULL_18V = 3.23e-72
OLD_FLUX_EXPANSION_18V = 2.74e9
Q_C = 1.602176634e-19


def read_csv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return [
            {str(k).strip(): (v.strip() if isinstance(v, str) else "") for k, v in row.items()}
            for row in csv.DictReader(handle, skipinitialspace=True)
        ]


def write_csv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: "" if row.get(field) is None else row.get(field) for field in fields})


def write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def fnum(value: Any) -> float | None:
    if value in (None, ""):
        return None
    try:
        parsed = float(value)
    except (TypeError, ValueError):
        return None
    return parsed if math.isfinite(parsed) else None


def median(values: list[float]) -> float | None:
    clean = sorted(v for v in values if math.isfinite(v))
    if not clean:
        return None
    mid = len(clean) // 2
    return clean[mid] if len(clean) % 2 else 0.5 * (clean[mid - 1] + clean[mid])


def quantile(values: list[float], q: float) -> float | None:
    clean = sorted(v for v in values if math.isfinite(v))
    if not clean:
        return None
    index = round((len(clean) - 1) * q)
    return clean[max(0, min(len(clean) - 1, index))]


def bias_close(value: Any, target: float, tol: float = 1.0e-6) -> bool:
    parsed = fnum(value)
    return parsed is not None and abs(parsed - target) <= tol


def full_bias_ladder() -> list[float]:
    values = {round(float(v), 12) for v in compare.PREVIOUS_FULL20_BIAS_POINTS}
    values.update(round(v, 12) for v in ANCHOR_BIASES)
    return sorted(values, reverse=True)


def ensure_sentaurus_exports(out_dir: Path) -> tuple[dict[float, Path], dict[str, Any]]:
    export_root = out_dir / "sentaurus_multibias"
    exports: dict[float, Path] = {}
    source_note: dict[str, Any] = {
        "source": "legacy_raw_sentaurus_multibias_copy",
        "reason": "source fixture has no pn2d_bv_multibias_*.tdr files; derived CSV/PNG/JSON are recomputed in this directory",
        "legacy_source_dir": str(LEGACY_SENT_EXPORT),
    }
    for bias in NODE_COMPARE_BIASES:
        name = f"sentaurus_{compare.signed_bias_token(bias)}v"
        src = LEGACY_SENT_EXPORT / name
        dst = export_root / name
        if not src.exists():
            continue
        if not dst.exists():
            shutil.copytree(src, dst)
        exports[compare.bias_key(bias)] = dst
    return exports, source_note


def run_runner(config_path: Path, cwd: Path) -> dict[str, Any]:
    proc = subprocess.run(
        [str(RUNNER), "--config", str(config_path)],
        cwd=cwd,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    (cwd / "runner_stdout.txt").write_text(proc.stdout, encoding="utf-8")
    (cwd / "runner_stderr.txt").write_text(proc.stderr, encoding="utf-8")
    return {"returncode": proc.returncode, "stdout_tail": proc.stdout[-2000:], "stderr_tail": proc.stderr[-2000:]}


def node_state_rows(node_rows: list[dict[str, Any]], target_bias: float) -> dict[int, dict[str, float]]:
    by_node: dict[int, dict[str, float]] = {}
    for row in node_rows:
        if not bias_close(row.get("bias_V"), target_bias):
            continue
        node_raw = row.get("node_id")
        quantity = str(row.get("quantity", "")).strip()
        if node_raw in (None, ""):
            continue
        node = int(node_raw)
        item = by_node.setdefault(node, {
            "node_id": node,
            "x_um": fnum(row.get("x_um")) or 0.0,
            "y_um": fnum(row.get("y_um")) or 0.0,
        })
        sent = fnum(row.get("sentaurus_value"))
        vela = fnum(row.get("vela_value_scaled_to_sentaurus_units"))
        if sent is not None:
            item[f"sent_{quantity}"] = sent
        if vela is not None:
            item[f"vela_{quantity}"] = vela
    return by_node


def bias_summary_rows(node_rows: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    summary: list[dict[str, Any]] = []
    spatial: list[dict[str, Any]] = []
    for bias in NODE_COMPARE_BIASES:
        by_node = node_state_rows(node_rows, bias)
        deltas: list[float] = []
        logn: list[float] = []
        abs_deltas_high: list[float] = []
        abs_deltas_other: list[float] = []
        sent_fields = [
            item["sent_electric_field"]
            for item in by_node.values()
            if "sent_electric_field" in item
        ]
        high_cut = quantile(sent_fields, 0.75)
        for item in by_node.values():
            if all(k in item for k in ("sent_potential", "sent_electron_qf", "vela_potential", "vela_electron_qf")):
                delta = (
                    item["vela_potential"] - item["vela_electron_qf"]
                    - (item["sent_potential"] - item["sent_electron_qf"])
                )
                item["delta_psi_minus_phin_V"] = delta
                deltas.append(delta)
                if high_cut is not None and item.get("sent_electric_field", -math.inf) >= high_cut:
                    abs_deltas_high.append(abs(delta))
                else:
                    abs_deltas_other.append(abs(delta))
            if item.get("sent_electron_density", 0.0) > 0.0 and item.get("vela_electron_density", 0.0) > 0.0:
                value = math.log10(item["vela_electron_density"] / item["sent_electron_density"])
                item["log10_n_vela_over_sentaurus"] = value
                logn.append(value)
            if bias == -18.0:
                spatial.append({
                    "bias_V": bias,
                    "node_id": item["node_id"],
                    "x_um": item["x_um"],
                    "y_um": item["y_um"],
                    "sent_electric_field_V_cm": item.get("sent_electric_field"),
                    "delta_psi_minus_phin_V": item.get("delta_psi_minus_phin_V"),
                    "log10_n_vela_over_sentaurus": item.get("log10_n_vela_over_sentaurus"),
                })
        summary.append({
            "bias_V": bias,
            "node_count": len(by_node),
            "delta_node_count": len(deltas),
            "median_delta_psi_minus_phin_V": median(deltas),
            "p95_abs_delta_psi_minus_phin_V": quantile([abs(v) for v in deltas], 0.95),
            "median_log10_n_vela_over_sentaurus": median(logn),
            "p95_abs_log10_n_vela_over_sentaurus": quantile([abs(v) for v in logn], 0.95),
            "sent_electric_field_top_quartile_cut_V_cm": high_cut,
            "median_abs_delta_high_field_top_quartile_V": median(abs_deltas_high),
            "median_abs_delta_other_nodes_V": median(abs_deltas_other),
        })
    return summary, spatial


def sent_flux_from_current_a_cm2(value: float | None) -> float | None:
    return abs(value) * 1.0e4 / Q_C if value is not None else None


def replay_ratio_rows(support_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for bias in (-18.0, -20.0):
        sums = {
            "vela_alpha_vela_flux": 0.0,
            "sent_alpha_vela_flux": 0.0,
            "vela_alpha_sent_flux": 0.0,
            "sent_full": 0.0,
            "vela_flux": 0.0,
            "sent_flux": 0.0,
        }
        edge_count = 0
        for row in support_rows:
            if not bias_close(row.get("bias_V"), bias) or not bias_close(row.get("nearest_sentaurus_bias_V"), bias):
                continue
            e_flux = fnum(row.get("electron_flux_abs")) or 0.0
            h_flux = fnum(row.get("hole_flux_abs")) or 0.0
            e_alpha = fnum(row.get("electron_alpha_m_inv")) or 0.0
            h_alpha = fnum(row.get("hole_alpha_m_inv")) or 0.0
            sent_e_alpha = fnum(row.get("sent_e_alpha")) or 0.0
            sent_h_alpha = fnum(row.get("sent_h_alpha")) or 0.0
            sent_e_flux = sent_flux_from_current_a_cm2(fnum(row.get("sent_e_current"))) or 0.0
            sent_h_flux = sent_flux_from_current_a_cm2(fnum(row.get("sent_h_current"))) or 0.0
            sums["vela_alpha_vela_flux"] += e_alpha * e_flux + h_alpha * h_flux
            sums["sent_alpha_vela_flux"] += sent_e_alpha * e_flux + sent_h_alpha * h_flux
            sums["vela_alpha_sent_flux"] += e_alpha * sent_e_flux + h_alpha * sent_h_flux
            sums["sent_full"] += sent_e_alpha * sent_e_flux + sent_h_alpha * sent_h_flux
            sums["vela_flux"] += e_flux + h_flux
            sums["sent_flux"] += sent_e_flux + sent_h_flux
            edge_count += 1
        sent_full = sums["sent_full"]
        rows.append({
            "bias_V": bias,
            "edge_count": edge_count,
            "vela_alpha_vela_flux_over_sentaurus_full": sums["vela_alpha_vela_flux"] / sent_full if sent_full else None,
            "sentaurus_alpha_vela_flux_over_sentaurus_full": sums["sent_alpha_vela_flux"] / sent_full if sent_full else None,
            "vela_alpha_sentaurus_flux_over_sentaurus_full": sums["vela_alpha_sent_flux"] / sent_full if sent_full else None,
            "vela_flux_over_sentaurus_flux": sums["vela_flux"] / sums["sent_flux"] if sums["sent_flux"] else None,
            **sums,
            "old_18v_vela_alpha_vela_flux_over_sentaurus_full": OLD_VELA_ALPHA_FLUX_OVER_SENT_FULL_18V if bias == -18.0 else None,
            "old_18v_vela_flux_over_sentaurus_flux": OLD_FLUX_EXPANSION_18V if bias == -18.0 else None,
        })
    return rows


def load_contact_node_ids(export_dir: Path, contact_name: str = "Anode") -> list[int]:
    for row in read_csv(export_dir / "contacts.csv"):
        if row.get("name") == contact_name:
            return [int(token) for token in row.get("node_ids", "").split(";") if token.strip()]
    return []


def sentaurus_current_proxy_curve(sentaurus_exports: dict[float, Path]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for bias in NODE_COMPARE_BIASES:
        export = sentaurus_exports.get(compare.bias_key(bias))
        if export is None:
            continue
        loaded = compare.load_sentaurus_scalar(export, ["TotalCurrentDensity"])
        contact_nodes = load_contact_node_ids(export, "Anode")
        values = []
        if loaded:
            _, data = loaded
            values = [abs(data[node]) for node in contact_nodes if node in data]
        rows.append({
            "bias_V": bias,
            "anode_contact_node_count": len(contact_nodes),
            "sentaurus_anode_total_current_density_abs_median_A_cm2": median(values),
            "sentaurus_anode_total_current_density_abs_mean_A_cm2": (sum(values) / len(values) if values else None),
        })
    return rows


def merge_iv_rows(vela_csv: Path, sent_proxy: list[dict[str, Any]]) -> list[dict[str, Any]]:
    sent_by_bias = {compare.bias_key(float(row["bias_V"])): row for row in sent_proxy}
    rows: list[dict[str, Any]] = []
    for row in read_csv(vela_csv):
        bias = fnum(row.get("bias_V"))
        if bias is None:
            continue
        sent = sent_by_bias.get(compare.bias_key(bias), {})
        vela_i = fnum(row.get("current_total_A_per_um"))
        sent_median = fnum(sent.get("sentaurus_anode_total_current_density_abs_median_A_cm2"))
        rows.append({
            "bias_V": bias,
            "converged": row.get("converged"),
            "iterations": row.get("iterations"),
            "newton_iterations": row.get("newton_iterations"),
            "failure_reason": row.get("failure_reason"),
            "vela_current_total_A_per_um": vela_i,
            "abs_vela_current_total_A_per_um": abs(vela_i) if vela_i is not None else None,
            "sentaurus_anode_total_current_density_abs_median_A_cm2": sent_median,
            "proxy_ratio_abs_vela_A_per_um_over_sentaurus_A_cm2": (
                abs(vela_i) / sent_median if vela_i is not None and sent_median not in (None, 0.0) else None
            ),
        })
    return rows


def load_plotting() -> Any:
    mpl_config = OUT_DIR / ".matplotlib"
    mpl_config.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_config))
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker

    plt.rcParams.update(
        {
            "figure.facecolor": "#FCFCFD",
            "axes.facecolor": "#FFFFFF",
            "axes.edgecolor": "#D7DBE7",
            "axes.labelcolor": "#1F2430",
            "grid.color": "#E6E8F0",
            "axes.grid": True,
            "grid.alpha": 0.45,
            "font.family": "sans-serif",
            "font.sans-serif": ["Aptos", "Inter", "Segoe UI", "DejaVu Sans", "Arial", "sans-serif"],
            "font.size": 12,
            "axes.titlesize": 14,
            "axes.labelsize": 12,
            "legend.fontsize": 10,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
        }
    )
    return plt, mticker


def save_scatter_map(rows: list[dict[str, Any]], value_col: str, title: str, subtitle: str, path: Path) -> None:
    plt, _ = load_plotting()
    xs = [fnum(row.get("x_um")) for row in rows]
    ys = [fnum(row.get("y_um")) for row in rows]
    vals = [fnum(row.get(value_col)) for row in rows]
    clean = [(x, y, v) for x, y, v in zip(xs, ys, vals) if x is not None and y is not None and v is not None]
    fig, ax = plt.subplots(figsize=(8.4, 5.0), constrained_layout=True)
    scatter = ax.scatter(
        [x for x, _, _ in clean],
        [y for _, y, _ in clean],
        c=[v for _, _, v in clean],
        s=120,
        cmap="coolwarm",
        edgecolors="#464C55",
        linewidths=0.7,
    )
    fig.colorbar(scatter, ax=ax, label=value_col)
    ax.set_xlabel("x (um)")
    ax.set_ylabel("y (um)")
    ax.set_title(title)
    ax.text(0.0, 1.02, subtitle, transform=ax.transAxes, fontsize=9, color="#6F768A")
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=220)
    plt.close(fig)


def save_cutline_plot(spatial_by_bias: dict[float, dict[int, dict[str, float]]], path: Path) -> None:
    plt, _ = load_plotting()
    fig, ax = plt.subplots(figsize=(8.8, 5.2), constrained_layout=True)
    colors = ["#5477C4", "#B8A037", "#71B436", "#CC6F47", "#BD569B"]
    for color, bias in zip(colors, NODE_COMPARE_BIASES):
        nodes = list(spatial_by_bias[bias].values())
        if not nodes:
            continue
        y_counts: dict[float, int] = {}
        for node in nodes:
            y_counts[node["y_um"]] = y_counts.get(node["y_um"], 0) + 1
        cut_y = max(y_counts.items(), key=lambda item: item[1])[0]
        cut = sorted(
            (node for node in nodes if abs(node["y_um"] - cut_y) <= 1.0e-9 and "delta_psi_minus_phin_V" in node),
            key=lambda node: node["x_um"],
        )
        ax.plot(
            [node["x_um"] for node in cut],
            [node["delta_psi_minus_phin_V"] for node in cut],
            marker="o",
            linewidth=1.0,
            color=color,
            label=f"{bias:g} V, y={cut_y:g} um",
        )
    ax.axhline(0.0, color="#464C55", linestyle=":", linewidth=1.0)
    ax.set_xlabel("x (um)")
    ax.set_ylabel("Delta (psi - phin), Vela - Sentaurus (V)")
    ax.set_title("PN2D center cutline residual across five reverse biases")
    ax.text(0.0, 1.02, "Exact Vela bias points after contact and handoff fixes", transform=ax.transAxes, fontsize=9, color="#6F768A")
    ax.legend(loc="best", fontsize=8)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=220)
    plt.close(fig)


def save_replay_plot(rows: list[dict[str, Any]], path: Path) -> None:
    plt, mticker = load_plotting()
    labels = [
        ("vela_alpha_vela_flux_over_sentaurus_full", "Vela alpha x Vela flux"),
        ("sentaurus_alpha_vela_flux_over_sentaurus_full", "Sentaurus alpha x Vela flux"),
        ("vela_alpha_sentaurus_flux_over_sentaurus_full", "Vela alpha x Sentaurus flux"),
    ]
    fig, ax = plt.subplots(figsize=(9.0, 5.2), constrained_layout=True)
    x = list(range(len(rows)))
    width = 0.24
    colors = ["#A3BEFA", "#FFE15B", "#A3D576"]
    for offset, (key, label), color in zip([-width, 0.0, width], labels, colors):
        ax.bar(
            [i + offset for i in x],
            [max(fnum(row.get(key)) or 1.0e-300, 1.0e-300) for row in rows],
            width=width,
            label=label,
            color=color,
            edgecolor="#464C55",
            linewidth=0.8,
        )
    ax.set_yscale("log")
    ax.set_xticks(x, [f"{row['bias_V']:g} V" for row in rows])
    ax.yaxis.set_major_formatter(mticker.LogFormatterSciNotation())
    ax.set_ylabel("Ratio to Sentaurus alpha x Sentaurus flux")
    ax.set_title("Replay ratios after contact and handoff fixes")
    ax.text(0.0, 1.02, "Edge-summed ratios at exact -18 V and -20 V anchors", transform=ax.transAxes, fontsize=9, color="#6F768A")
    ax.legend(loc="best", fontsize=8)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=220)
    plt.close(fig)


def save_iv_plot(rows: list[dict[str, Any]], path: Path) -> None:
    plt, mticker = load_plotting()
    fig, ax = plt.subplots(figsize=(9.0, 5.2), constrained_layout=True)
    sorted_rows = sorted(rows, key=lambda row: abs(float(row["bias_V"])))
    ax.semilogy(
        [abs(float(row["bias_V"])) for row in sorted_rows],
        [max(fnum(row.get("abs_vela_current_total_A_per_um")) or 1.0e-300, 1.0e-300) for row in sorted_rows],
        marker="o",
        color="#5477C4",
        label="Vela terminal |I| (A/um)",
    )
    proxy = [row for row in sorted_rows if fnum(row.get("sentaurus_anode_total_current_density_abs_median_A_cm2")) is not None]
    if proxy:
        ax2 = ax.twinx()
        ax2.semilogy(
            [abs(float(row["bias_V"])) for row in proxy],
            [max(fnum(row.get("sentaurus_anode_total_current_density_abs_median_A_cm2")) or 1.0e-300, 1.0e-300) for row in proxy],
            marker="s",
            color="#CC6F47",
            label="Sentaurus anode |TotalCurrentDensity| proxy (A/cm2)",
        )
        ax2.set_ylabel("Sentaurus proxy |J| (A/cm2)")
        ax2.yaxis.set_major_formatter(mticker.LogFormatterSciNotation())
        lines, labels = ax.get_legend_handles_labels()
        lines2, labels2 = ax2.get_legend_handles_labels()
        ax.legend(lines + lines2, labels + labels2, loc="best", fontsize=8)
    else:
        ax.legend(loc="best", fontsize=8)
    ax.set_xlabel("Reverse bias |V| (V)")
    ax.set_ylabel("Vela terminal |I| (A/um)")
    ax.yaxis.set_major_formatter(mticker.LogFormatterSciNotation())
    ax.set_title("BV current curve with Sentaurus field proxy")
    ax.text(0.0, 1.02, "Sentaurus .plt terminal current is unavailable in the fixture; proxy uses anode node TotalCurrentDensity median", transform=ax.transAxes, fontsize=9, color="#6F768A")
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=220)
    plt.close(fig)


def save_alpha_flux_scatter(support_rows: list[dict[str, Any]], bias: float, path: Path) -> None:
    plt, _ = load_plotting()
    xs: list[float] = []
    ys: list[float] = []
    for row in support_rows:
        if not bias_close(row.get("bias_V"), bias) or not bias_close(row.get("nearest_sentaurus_bias_V"), bias):
            continue
        e_flux = fnum(row.get("electron_flux_abs")) or 0.0
        h_flux = fnum(row.get("hole_flux_abs")) or 0.0
        e_alpha = fnum(row.get("electron_alpha_m_inv")) or 0.0
        h_alpha = fnum(row.get("hole_alpha_m_inv")) or 0.0
        sent_e_flux = sent_flux_from_current_a_cm2(fnum(row.get("sent_e_current"))) or 0.0
        sent_h_flux = sent_flux_from_current_a_cm2(fnum(row.get("sent_h_current"))) or 0.0
        sent_e_alpha = fnum(row.get("sent_e_alpha")) or 0.0
        sent_h_alpha = fnum(row.get("sent_h_alpha")) or 0.0
        sent_full = sent_e_alpha * sent_e_flux + sent_h_alpha * sent_h_flux
        vela_full = e_alpha * e_flux + h_alpha * h_flux
        if sent_full > 0.0 and vela_full > 0.0:
            xs.append(sent_full)
            ys.append(vela_full)
    fig, ax = plt.subplots(figsize=(7.8, 5.2), constrained_layout=True)
    ax.scatter(xs, ys, s=72, color="#A3BEFA", edgecolors="#2E4780", linewidths=0.8, alpha=0.8)
    if xs and ys:
        lo = min(min(xs), min(ys))
        hi = max(max(xs), max(ys))
        ax.plot([lo, hi], [lo, hi], color="#464C55", linestyle=":", linewidth=1.0)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Sentaurus alpha x Sentaurus flux")
    ax.set_ylabel("Vela alpha x Vela flux")
    ax.set_title(f"Edge alpha-flux scatter at {bias:g} V")
    ax.text(0.0, 1.02, "Each point is one SG avalanche diagnostic edge; diagonal marks parity", transform=ax.transAxes, fontsize=9, color="#6F768A")
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=220)
    plt.close(fig)


def write_report(path: Path, summary: dict[str, Any]) -> None:
    invalid_dirs = summary["invalid_old_report_dirs"]
    lines = [
        "# PN2D coarse BV rebaseline after contact and handoff fixes",
        "",
        "## Technical summary",
        "",
        "- The repaired baseline reaches -20 V with exact anchor points and recomputed node/edge diagnostics in this directory.",
        "- Old derived report data are excluded from the new evidence set; raw Sentaurus field exports were copied into this directory because the fixture does not include the original multibias TDR files.",
        "- Sentaurus .plt terminal current is unavailable in the current fixture, so terminal-current comparison is limited to Vela terminal current plus a Sentaurus contact-field current-density proxy.",
        "",
        "## Invalid old report directories",
        "",
    ]
    lines.extend(f"- `{item}`" for item in invalid_dirs)
    lines.extend([
        "",
        "## Outputs",
        "",
        f"- Bias summary: `{summary['bias_summary_csv']}`",
        f"- Replay ratios: `{summary['replay_ratios_csv']}`",
        f"- BV current proxy: `{summary['bv_curve_proxy_csv']}`",
        f"- Spatial residual CSV: `{summary['spatial_residual_csv']}`",
        f"- Plot directory: `{summary['plot_dir']}`",
        "",
    ])
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    plot_dir = OUT_DIR / "plots"
    sentaurus_exports, sentaurus_source_note = ensure_sentaurus_exports(OUT_DIR)
    ladder = full_bias_ladder()

    config_path = compare.write_previous_full20_config(
        base_config=IMPORT_DIR / "vela" / "simulation_bv.json",
        out_dir=OUT_DIR,
        output_csv_name="coarse_rebaseline_post_contact_handoff_fix.csv",
        bias_points=ladder,
        config_name="simulation_rebaseline_post_contact_handoff_fix.json",
        vtk_subdir="vtk_rebaseline",
        diagnostics_suffix="_rebaseline",
    )
    run_status = run_runner(config_path, OUT_DIR)
    vela_csv = OUT_DIR / "coarse_rebaseline_post_contact_handoff_fix.csv"
    vtk_dir = OUT_DIR / "vtk_rebaseline"
    vela_vtks = compare.discover_vela_vtks(vtk_dir)

    node_rows = compare.node_field_compare(
        sentaurus_exports=sentaurus_exports,
        vela_vtks=vela_vtks,
        biases=NODE_COMPARE_BIASES,
        exact_vela_bias=True,
    )
    node_csv = OUT_DIR / "coarse_node_field_compare_rebaseline.csv"
    compare.write_csv_rows(node_csv, node_rows, [
        "bias_V", "vela_bias_V", "node_id", "x_um", "y_um", "quantity", "sentaurus_field",
        "sentaurus_value", "vela_field", "vela_value_scaled_to_sentaurus_units",
        "comparison_basis", "vela_over_sentaurus", "diff", "abs_diff",
    ])

    support_rows = compare.current_support_compare(
        sg_edges_csv=OUT_DIR / "sg_avalanche_edges_rebaseline.csv",
        sentaurus_exports=sentaurus_exports,
    )
    support_csv = OUT_DIR / "coarse_current_support_compare_rebaseline.csv"
    compare.write_csv_rows(support_csv, support_rows, [
        "bias_V", "nearest_sentaurus_bias_V", "edge_id", "node0", "node1", "edge_class",
        "source_integral_total", "electron_source_integral", "hole_source_integral",
        "current_comparison_basis",
        "electric_field_V_m", "electron_qf_field_V_m", "hole_qf_field_V_m",
        "electron_flux_abs", "hole_flux_abs", "electron_alpha_m_inv", "hole_alpha_m_inv",
        "sent_e_velocity", "sent_h_velocity", "sent_e_alpha", "sent_h_alpha",
        "sent_e_ion_integral", "sent_h_ion_integral", "sent_mean_ion_integral",
        "sent_e_current", "sent_h_current",
    ])

    bias_rows, spatial_rows = bias_summary_rows(node_rows)
    bias_summary_csv = OUT_DIR / "rebaseline_bias_summary.csv"
    write_csv(bias_summary_csv, bias_rows, [
        "bias_V", "node_count", "delta_node_count", "median_delta_psi_minus_phin_V",
        "p95_abs_delta_psi_minus_phin_V", "median_log10_n_vela_over_sentaurus",
        "p95_abs_log10_n_vela_over_sentaurus", "sent_electric_field_top_quartile_cut_V_cm",
        "median_abs_delta_high_field_top_quartile_V", "median_abs_delta_other_nodes_V",
    ])
    spatial_csv = OUT_DIR / "rebaseline_spatial_residual_minus18.csv"
    write_csv(spatial_csv, spatial_rows, [
        "bias_V", "node_id", "x_um", "y_um", "sent_electric_field_V_cm",
        "delta_psi_minus_phin_V", "log10_n_vela_over_sentaurus",
    ])

    replay_rows = replay_ratio_rows(support_rows)
    replay_csv = OUT_DIR / "rebaseline_replay_ratios.csv"
    write_csv(replay_csv, replay_rows, [
        "bias_V", "edge_count",
        "vela_alpha_vela_flux_over_sentaurus_full",
        "sentaurus_alpha_vela_flux_over_sentaurus_full",
        "vela_alpha_sentaurus_flux_over_sentaurus_full",
        "vela_flux_over_sentaurus_flux",
        "vela_alpha_vela_flux", "sent_alpha_vela_flux", "vela_alpha_sent_flux",
        "sent_full", "vela_flux", "sent_flux",
        "old_18v_vela_alpha_vela_flux_over_sentaurus_full",
        "old_18v_vela_flux_over_sentaurus_flux",
    ])

    sent_proxy = sentaurus_current_proxy_curve(sentaurus_exports)
    iv_rows = merge_iv_rows(vela_csv, sent_proxy)
    iv_csv = OUT_DIR / "rebaseline_bv_curve_proxy.csv"
    write_csv(iv_csv, iv_rows, [
        "bias_V", "converged", "iterations", "newton_iterations", "failure_reason",
        "vela_current_total_A_per_um", "abs_vela_current_total_A_per_um",
        "sentaurus_anode_total_current_density_abs_median_A_cm2",
        "proxy_ratio_abs_vela_A_per_um_over_sentaurus_A_cm2",
    ])

    spatial_by_bias = {bias: node_state_rows(node_rows, bias) for bias in NODE_COMPARE_BIASES}
    save_replay_plot(replay_rows, plot_dir / "replay_ratios_minus18_minus20.png")
    save_iv_plot(iv_rows, plot_dir / "bv_current_proxy.png")
    save_scatter_map(spatial_rows, "delta_psi_minus_phin_V", "-18 V delta(psi - phin) node map", "Vela minus Sentaurus at recomputed fixed baseline", plot_dir / "delta_psi_phin_map_minus18.png")
    save_scatter_map(spatial_rows, "log10_n_vela_over_sentaurus", "-18 V log10 electron-density ratio node map", "Positive means Vela higher; negative means Vela lower", plot_dir / "log10_n_ratio_map_minus18.png")
    save_cutline_plot(spatial_by_bias, plot_dir / "cutline_delta_psi_phin_five_biases.png")
    save_alpha_flux_scatter(support_rows, -18.0, plot_dir / "alpha_flux_scatter_minus18.png")
    save_alpha_flux_scatter(support_rows, -20.0, plot_dir / "alpha_flux_scatter_minus20.png")

    invalid_old_dirs = [
        "build-release/reference_tcad/pn2d_sentaurus2018_coarse7x3/reports/coarse_vm_vector_compare",
        "build-release/reference_tcad/pn2d_sentaurus2018_coarse7x3/reports/coarse_psi_gradient_proxy_20260705",
        "build-release/reference_tcad/pn2d_sentaurus2018_coarse7x3/reports/coarse_psi_gradient_proxy_1vgrid_20260705",
        "build-release/reference_tcad/pn2d_sentaurus2018_coarse7x3/reports/qf_cap_warmstart_branch_20260705",
        "build-release/reference_tcad/pn2d_sentaurus2018_coarse7x3/reports/active_region_branch_feedback_20260705",
    ]
    invalid_path = OUT_DIR / "invalid_old_report_dirs.md"
    invalid_path.write_text(
        "# Invalidated old PN2D coarse report directories\n\n" +
        "\n".join(f"- `{item}`" for item in invalid_old_dirs) +
        "\n",
        encoding="utf-8",
    )

    summary = {
        "out_dir": str(OUT_DIR),
        "config": str(config_path),
        "runner": str(RUNNER),
        "run_status": run_status,
        "bias_ladder_count": len(ladder),
        "bias_ladder": ladder,
        "vela_point_count": len(read_csv(vela_csv)),
        "all_vela_points_converged": all(row.get("converged") == "1" for row in read_csv(vela_csv)),
        "sentaurus_source_note": sentaurus_source_note,
        "sentaurus_exports": {str(k): str(v) for k, v in sentaurus_exports.items()},
        "node_compare_csv": str(node_csv),
        "support_compare_csv": str(support_csv),
        "bias_summary_csv": str(bias_summary_csv),
        "replay_ratios_csv": str(replay_csv),
        "bv_curve_proxy_csv": str(iv_csv),
        "spatial_residual_csv": str(spatial_csv),
        "plot_dir": str(plot_dir),
        "invalid_old_report_dirs": invalid_old_dirs,
        "invalid_old_report_dirs_md": str(invalid_path),
    }
    write_json(OUT_DIR / "rebaseline_summary.json", summary)
    write_report(OUT_DIR / "rebaseline_report.md", summary)
    print(json.dumps(summary, indent=2), flush=True)
    return 0 if run_status["returncode"] == 0 else run_status["returncode"]


if __name__ == "__main__":
    raise SystemExit(main())
