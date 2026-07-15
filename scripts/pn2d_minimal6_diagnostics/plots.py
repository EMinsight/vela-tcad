"""Static, contract-checked figures for PN2D minimal6 diagnostics."""
from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .schemas import DISCLAIMER

FIGURES = (
    ("gradient", "Recovered electrostatic-gradient components", "V/m"),
    ("current_alpha", "Signed carrier-current projections and alpha", "A/cm^2 and cm^-1"),
    ("source_waterfall", "Native and reconstructed avalanche-source anchors", "s^-1 per 1 cm depth"),
    ("interaction", "Counterfactual interaction matrix", "dex"),
    ("topology_symmetry", "Sketch/mirror source symmetry", "s^-1 per 1 cm depth"),
)


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _number(value: str | None) -> float | None:
    try:
        return float(value) if value not in (None, "") else None
    except ValueError:
        return None


def _finish(fig: plt.Figure, ax: plt.Axes, title: str, unit: str) -> None:
    ax.set_title(title)
    ax.set_ylabel(unit)
    ax.grid(True, alpha=0.25)
    fig.text(0.01, 0.01, DISCLAIMER, fontsize=8)
    fig.tight_layout(rect=(0, 0.04, 1, 1))


def _plot_gradient(ax: plt.Axes, ledger: list[dict[str, str]]) -> None:
    rows = [row for row in ledger if row["record_kind"] == "cell_replay" and row["quantity"] == "minus_grad_psi" and row["component"] in ("x", "y")]
    for component, marker in (("x", "o"), ("y", "s")):
        selected = [row for row in rows if row["component"] == component]
        ax.scatter(range(len(selected)), [_number(row["value"]) for row in selected], marker=marker, label=f"-grad(psi) {component}")
    ax.set_xlabel("cell replay row (signed components retained)")
    ax.legend(loc="best")


def _plot_current_alpha(ax: plt.Axes, ledger: list[dict[str, str]]) -> None:
    rows = [row for row in ledger if row["record_kind"] == "edge_raw" and row["component"] == "signed_projection"]
    for carrier, color in (("electron", "tab:blue"), ("hole", "tab:orange")):
        selected = [row for row in rows if row["quantity"] == f"{carrier}_current"]
        ax.plot(range(len(selected)), [_number(row["value"]) for row in selected], marker="o", color=color, label=f"{carrier} signed J projection")
    ax.set_xlabel("edge row (edge direction defines the sign)")
    ax.legend(loc="best")


def _source_rows(ledger: list[dict[str, str]]) -> list[dict[str, str]]:
    return [row for row in ledger if row["record_kind"] == "source_integral"]


def _plot_sources(ax: plt.Axes, ledger: list[dict[str, str]]) -> None:
    rows = _source_rows(ledger)
    labels: list[str] = []
    values: list[float] = []
    colors: list[str] = []
    palette = {
        "sentaurus_native_avalanche_generation": "tab:green",
        "sentaurus_alpha_current_reconstruction": "tab:purple",
        "vela_alpha_flux_partial_volume_reconstruction": "tab:red",
    }
    for row in rows:
        value = _number(row["value_s_inv_per_unit_depth"])
        labels.append(f"{row['topology']} {row['bias_V']} V\n{row['source']}")
        values.append(0.0 if value is None else value)
        colors.append(palette.get(row["source"], "tab:gray"))
    ax.bar(range(len(values)), values, color=colors)
    zero_positions = [index for index, value in enumerate(values) if value == 0.0]
    if zero_positions:
        ax.scatter(zero_positions, [0.0] * len(zero_positions), marker="x", color="black", label="geometric zero (not log floored)")
        ax.legend(loc="best")
    ax.set_xticks(range(len(labels)), labels, rotation=90, fontsize=7)
    ax.set_xlabel("native Sentaurus (green), Sentaurus reconstruction (purple), Vela reconstruction (red)")


def _plot_interactions(ax: plt.Axes, report: dict[str, Any]) -> None:
    interactions = report.get("interactions", [])
    if not interactions:
        ax.text(0.5, 0.5, "No interaction crossed the 0.3 dex trigger.\nUnavailable is reported explicitly.", ha="center", va="center")
        ax.set_xticks([])
        return
    labels = [f"{item['factor_a']} / {item['factor_b']}" for item in interactions]
    values = [float(item["interaction_dex"]) for item in interactions]
    ax.bar(range(len(values)), values, color="tab:cyan")
    ax.set_xticks(range(len(labels)), labels, rotation=30, ha="right")


def _plot_symmetry(ax: plt.Axes, ledger: list[dict[str, str]]) -> None:
    rows = [row for row in _source_rows(ledger) if row["source"] == "sentaurus_native_avalanche_generation"]
    by_key = {(row["topology"], row["bias_V"]): _number(row["value_s_inv_per_unit_depth"]) or 0.0 for row in rows}
    biases = sorted({bias for _, bias in by_key}, key=float, reverse=True)
    x = list(range(len(biases)))
    sketch = [by_key.get(("sketch", bias), 0.0) for bias in biases]
    mirror = [by_key.get(("mirror", bias), 0.0) for bias in biases]
    ax.plot(x, sketch, marker="o", label="sketch native source")
    ax.plot(x, mirror, marker="s", label="mirror native source")
    ax.set_xticks(x, [f"{bias} V" for bias in biases])
    ax.set_xlabel("exact fixed-state bias")
    ax.legend(loc="best")


def render_formula_difference_figures(
    *, ledger_path: Path, waterfall_path: Path, report_path: Path, out_dir: Path,
    reviewer: str = "unreviewed", reviewed_on: str = "", qa_status: str = "pending_visual_inspection",
) -> dict[str, Any]:
    """Render the fixed figure set and return its self-describing manifest."""
    ledger = _read_csv(ledger_path)
    # Reading the waterfall is intentional: it makes the figure build fail closed if that
    # required diagnostic artifact was not produced.
    _read_csv(waterfall_path)
    report = json.loads(report_path.read_text(encoding="utf-8"))
    plotters = (_plot_gradient, _plot_current_alpha, _plot_sources, _plot_interactions, _plot_symmetry)
    entries: list[dict[str, Any]] = []
    for (stem, title, unit), plotter in zip(FIGURES, plotters, strict=True):
        fig, ax = plt.subplots(figsize=(8, 4.5), dpi=120)
        if stem == "interaction":
            plotter(ax, report)
        else:
            plotter(ax, ledger)
        _finish(fig, ax, title, unit)
        artifacts = [f"{stem}.png", f"{stem}.pdf"]
        fig.savefig(out_dir / artifacts[0], dpi=120)
        fig.savefig(out_dir / artifacts[1])
        plt.close(fig)
        entries.append({"stem": stem, "title": title, "unit": unit, "artifacts": artifacts, "contains_disclaimer": True})
    return {
        "schema": "vela.pn2d_minimal6_figure_manifest.v1",
        "diagnostic_disclaimer": DISCLAIMER,
        "figures": entries,
        "manual_qa": {
            "reviewer": reviewer,
            "reviewed_on": reviewed_on,
            "checklist": [
                "node_edge_triangle_identities", "mirror_symmetry", "carrier_signs", "units",
                "geometric_zero_classification", "waterfall_closure", "native_reconstructed_labels",
            ],
            "status": qa_status,
        },
    }