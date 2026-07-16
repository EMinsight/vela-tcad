"""Static, contract-checked figures for PN2D minimal6 diagnostics."""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .schemas import DISCLAIMER

FIGURES = (
    ("gradient", "Recovered electrostatic-gradient components", "V/m"),
    ("current_alpha", "Signed carrier-current projections and alpha", "A/m^2 and m^-1"),
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
    if not ax.get_ylabel():
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
    alpha_rows = [
        row for row in ledger
        if row["record_kind"] == "edge_replay" and row["quantity"].endswith("_alpha_per_m")
    ]
    alpha_ax = ax.twinx()
    for carrier, color in (("electron", "tab:blue"), ("hole", "tab:orange")):
        selected = [row for row in rows if row["quantity"] == f"{carrier}_current"]
        ax.plot(range(len(selected)), [_number(row["value"]) for row in selected], marker="o", color=color, label=f"{carrier} signed J projection")
        selected_alpha = [row for row in alpha_rows if row["quantity"].endswith(f"{carrier}_alpha_per_m")]
        alpha_ax.plot(range(len(selected_alpha)), [_number(row["value"]) for row in selected_alpha], linestyle="--", color=color, label=f"{carrier} alpha")
    ax.set_xlabel("edge row (edge direction defines the sign)")
    ax.set_ylabel("A/m^2")
    alpha_ax.set_ylabel("m^-1")
    lines, labels = ax.get_legend_handles_labels()
    alpha_lines, alpha_labels = alpha_ax.get_legend_handles_labels()
    ax.legend(lines + alpha_lines, labels + alpha_labels, loc="best")


def _source_rows(ledger: list[dict[str, str]]) -> list[dict[str, str]]:
    return [row for row in ledger if row["record_kind"] == "source_integral"]


def _plot_sources(ax: plt.Axes, ledger: list[dict[str, str]]) -> None:
    rows = _source_rows(ledger)
    labels: list[str] = []
    palette = {
        "sentaurus_native_avalanche_generation": "tab:green",
        "sentaurus_alpha_current_reconstruction": "tab:purple",
        "vela_alpha_flux_partial_volume_reconstruction": "tab:red",
    }
    available_positions: list[int] = []
    available_values: list[float] = []
    available_colors: list[str] = []
    available_zero_positions: list[int] = []
    geometric_zero_positions: list[int] = []
    unavailable_positions: list[int] = []
    for position, row in enumerate(rows):
        labels.append(f"{row['topology']} {row['bias_V']} V\n{row['source']}")
        status = row.get("status")
        if status not in {"available", "geometric_zero", "unavailable"}:
            raise ValueError(f"invalid source-integral status: {status!r}")
        if status == "available":
            value = _number(row.get("value_s_inv_per_unit_depth"))
            if value is None or not math.isfinite(value):
                raise ValueError("available source integral requires a finite numeric value")
            available_positions.append(position)
            available_values.append(value)
            available_colors.append(palette.get(row["source"], "tab:gray"))
            if value == 0.0:
                available_zero_positions.append(position)
        elif status == "geometric_zero":
            geometric_zero_positions.append(position)
        else:
            unavailable_positions.append(position)

    ax.bar(available_positions, available_values, color=available_colors)
    if available_zero_positions:
        ax.scatter(
            available_zero_positions, [0.0] * len(available_zero_positions),
            marker="o", facecolors="none", edgecolors="black",
            label="available zero (not geometric)",
        )
    if geometric_zero_positions:
        ax.scatter(
            geometric_zero_positions, [0.0] * len(geometric_zero_positions),
            marker="x", color="black",
            label="explicit geometric zero (not log floored)",
        )
    if unavailable_positions:
        for position in unavailable_positions:
            ax.text(
                position, 0.02, "N/A", ha="center", va="bottom", rotation=90,
                transform=ax.get_xaxis_transform(),
            )
        ax.scatter(
            [], [], marker="s", facecolors="none", edgecolors="tab:gray",
            label="unavailable (not plotted)",
        )
    if available_zero_positions or geometric_zero_positions or unavailable_positions:
        ax.legend(loc="best")
    ax.set_xticks(range(len(labels)), labels, rotation=90, fontsize=7)
    ax.set_xlabel("native Sentaurus (green), Sentaurus reconstruction (purple), Vela reconstruction (red)")


def _plot_interactions(ax: plt.Axes, report: dict[str, Any]) -> None:
    interactions = report.get("interactions", [])
    if not interactions:
        ax.text(0.5, 0.5, "No interaction crossed the 0.3 dex trigger.\nUnavailable is reported explicitly.", ha="center", va="center")
        ax.set_xticks([])
        ax.set_yticks([])
        return
    states = list(dict.fromkeys(
        (str(item["topology"]), float(item["bias_V"])) for item in interactions
    ))
    columns = list(dict.fromkeys(
        (str(item["path_identity"]), str(item["first_factor"]), str(item["second_factor"]))
        for item in interactions
    ))
    state_index = {key: index for index, key in enumerate(states)}
    column_index = {key: index for index, key in enumerate(columns)}
    matrix = np.full((len(states), len(columns)), np.nan)
    for item in interactions:
        state = (str(item["topology"]), float(item["bias_V"]))
        column = (
            str(item["path_identity"]), str(item["first_factor"]),
            str(item["second_factor"]),
        )
        row_index = state_index[state]
        col_index = column_index[column]
        if not np.isnan(matrix[row_index, col_index]):
            raise ValueError(f"duplicate interaction figure identity: {state} {column}")
        matrix[row_index, col_index] = float(item["interaction_dex"])
    finite = np.abs(matrix[np.isfinite(matrix)])
    scale = max(float(finite.max()) if finite.size else 0.0, 1.0e-12)
    image = ax.imshow(matrix, cmap="coolwarm", vmin=-scale, vmax=scale, aspect="auto")
    for row_index, col_index in zip(*np.where(np.isfinite(matrix))):
        value = matrix[row_index, col_index]
        color = "white" if abs(value) >= 0.55 * scale else "black"
        ax.text(col_index, row_index, f"{value:.3g}", ha="center", va="center", color=color)
    path_labels = {"forward_adjacent": "F", "reverse_adjacent": "R"}
    labels = [
        f"{path_labels.get(path, path)}: {first}\n-> {second}"
        for path, first, second in columns
    ]
    ax.set_xticks(range(len(columns)), labels, rotation=0, ha="center", fontsize=8)
    ax.set_yticks(
        range(len(states)), [f"{topology} {bias:g} V" for topology, bias in states]
    )
    ax.set_xlabel("ordered adjacent factor pair and path direction")
    ax.set_ylabel("exact fixed state")
    ax.figure.colorbar(image, ax=ax, label="interaction (dex)", pad=0.02)


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
        figsize = (10, 6) if stem == "interaction" else (8, 4.5)
        fig, ax = plt.subplots(figsize=figsize, dpi=120)
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
