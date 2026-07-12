#!/usr/bin/env python3
"""Generate static scientific figures for the PN2D BV comparison matrix."""

from __future__ import annotations

import argparse
import base64
import csv
import json
import math
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
from PIL import Image


VARIANTS = (
    "legacy_density_gradient",
    "legacy_gss_midpoint",
    "legacy_triangle_gss_gradqf",
    "reported_density_gradient",
    "reported_gss_midpoint",
    "reported_triangle_gss_gradqf",
)

DISPLAY_NAMES = {
    "legacy_density_gradient": "Legacy / density gradient",
    "legacy_gss_midpoint": "Legacy / GSS midpoint",
    "legacy_triangle_gss_gradqf": "Legacy / triangle GSS",
    "reported_density_gradient": "Reported / density gradient",
    "reported_gss_midpoint": "Reported / GSS midpoint",
    "reported_triangle_gss_gradqf": "Reported / triangle GSS",
}

MODE_COLORS = {
    "density_gradient": "#0072B2",
    "gss_midpoint": "#009E73",
    "triangle_gss_gradqf": "#D55E00",
}


@dataclass
class MeshState:
    points: np.ndarray
    triangles: np.ndarray
    point_data: dict[str, np.ndarray]
    cell_data: dict[str, np.ndarray]


def parse_args() -> argparse.Namespace:
    repo = Path(__file__).resolve().parents[1]
    report_root = repo / (
        "build-release/reference_tcad/pn2d_sentaurus2018_coarse7x3/reports/"
        "pn2d_bv_compensated_gss_matrix_v3_full_status_20260712"
    )
    reference = repo / (
        "build-release/reference_tcad/pn2d_sentaurus2018_coarse7x3/"
        "imported_reference/reference_curves/"
        "pn2d_sentaurus2018_coarse7x3_bv_reference.csv"
    )
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=report_root)
    parser.add_argument("--sentaurus-reference", type=Path, default=reference)
    parser.add_argument("--out-dir", type=Path)
    parser.add_argument("--visual-html", type=Path)
    return parser.parse_args()


def configure_matplotlib() -> None:
    plt.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["DejaVu Serif"],
            "font.size": 8.5,
            "axes.titlesize": 9.5,
            "axes.labelsize": 9,
            "axes.linewidth": 0.7,
            "legend.fontsize": 7.5,
            "xtick.labelsize": 7.5,
            "ytick.labelsize": 7.5,
            "savefig.dpi": 300,
            "savefig.bbox": "tight",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def read_numeric_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as stream:
        return list(csv.DictReader(stream))


def truthy(value: str | None) -> bool:
    return str(value).strip().lower() in {"1", "true", "yes"}


def find_iv_csv(variant_root: Path) -> Path:
    for path in sorted((variant_root / "run").glob("pn2d_bv_*.csv")):
        with path.open(encoding="utf-8-sig") as stream:
            header = stream.readline()
        if "bias_V" in header and "current_total_A_per_um" in header:
            return path
    raise FileNotFoundError(f"No I-V CSV found under {variant_root / 'run'}")


def load_iv(path: Path, current_column: str) -> tuple[np.ndarray, np.ndarray]:
    biases: list[float] = []
    currents: list[float] = []
    for row in read_numeric_csv(path):
        if "converged" in row and not truthy(row["converged"]):
            continue
        try:
            bias = float(row["bias_V"])
            current = float(row[current_column])
        except (KeyError, TypeError, ValueError):
            continue
        if math.isfinite(bias) and math.isfinite(current):
            biases.append(-bias)
            currents.append(max(abs(current), 1e-300))
    order = np.argsort(biases)
    return np.asarray(biases)[order], np.asarray(currents)[order]


def variant_mode(variant: str) -> str:
    for mode in MODE_COLORS:
        if variant.endswith(mode):
            return mode
    raise ValueError(f"Unknown variant: {variant}")


def save_figure(fig: plt.Figure, out_dir: Path, stem: str) -> list[Path]:
    outputs = []
    for suffix in (".png", ".pdf"):
        path = out_dir / f"{stem}{suffix}"
        fig.savefig(path, facecolor="white")
        outputs.append(path)
    plt.close(fig)
    return outputs


def plot_iv(root: Path, reference: Path, out_dir: Path) -> tuple[list[Path], dict]:
    sent_x, sent_y = load_iv(reference, "current_total")
    series: dict[str, tuple[np.ndarray, np.ndarray]] = {}
    for variant in VARIANTS:
        series[variant] = load_iv(find_iv_csv(root / "variants" / variant), "current_total_A_per_um")

    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.2), constrained_layout=True)
    for ax in axes:
        ax.semilogy(sent_x, sent_y, color="#222222", linewidth=1.8, label="Sentaurus 2018")
        for variant, (x, y) in series.items():
            mode = variant_mode(variant)
            linestyle = "-" if variant.startswith("legacy_") else "--"
            marker = "o" if mode == "triangle_gss_gradqf" else None
            markevery = max(1, len(x) // 12) if marker else None
            ax.semilogy(
                x,
                y,
                color=MODE_COLORS[mode],
                linestyle=linestyle,
                linewidth=1.35,
                marker=marker,
                markersize=2.8,
                markevery=markevery,
                label=DISPLAY_NAMES[variant],
            )
        ax.set_xlabel("Reverse bias |V| (V)")
        ax.set_ylabel("|I| (A, 1 um device depth)")
        ax.grid(True, which="both", color="#d0d0d0", linewidth=0.45, alpha=0.7)

    axes[0].set_xlim(0.0, 20.05)
    axes[0].set_title("Full reverse-bias sweep")
    axes[1].set_xlim(18.0, 20.05)
    axes[1].set_title("Breakdown-region detail")
    axes[1].axvline(19.4, color="#777777", linestyle=":", linewidth=1.0)
    axes[1].annotate(
        "triangle GSS last converged point",
        xy=(19.4, 0.96),
        xycoords=("data", "axes fraction"),
        xytext=(-4, -2),
        textcoords="offset points",
        ha="right",
        va="top",
        fontsize=7.2,
        color="#555555",
    )
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="outside lower center", ncol=4, frameon=False)
    fig.suptitle("PN2D BV current-voltage comparison", fontweight="normal")
    outputs = save_figure(fig, out_dir, "pn2d-bv-iv-comparison")
    endpoints = {
        variant: {
            "last_reverse_bias_V": float(values[0][-1]),
            "last_abs_current_A_per_um": float(values[1][-1]),
            "points": int(len(values[0])),
        }
        for variant, values in series.items()
    }
    return outputs, endpoints


def _read_values(lines: list[str], start: int, count: int) -> tuple[np.ndarray, int]:
    values: list[float] = []
    index = start
    while len(values) < count:
        if index >= len(lines):
            raise ValueError("Unexpected end of VTK file")
        values.extend(float(value) for value in lines[index].split())
        index += 1
    if len(values) != count:
        raise ValueError(f"Expected {count} VTK values, found {len(values)}")
    return np.asarray(values), index


def load_legacy_vtk(path: Path) -> MeshState:
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    points: np.ndarray | None = None
    triangles: np.ndarray | None = None
    point_data: dict[str, np.ndarray] = {}
    cell_data: dict[str, np.ndarray] = {}
    data_kind: str | None = None
    data_count = 0
    index = 0
    while index < len(lines):
        parts = lines[index].split()
        token = parts[0]
        if token == "POINTS":
            count = int(parts[1])
            values, index = _read_values(lines, index + 1, count * 3)
            points = values.reshape(count, 3)[:, :2]
            continue
        if token == "CELLS":
            count = int(parts[1])
            rows = []
            index += 1
            for _ in range(count):
                row = [int(value) for value in lines[index].split()]
                if row[0] == 3:
                    rows.append(row[1:4])
                index += 1
            triangles = np.asarray(rows, dtype=int)
            continue
        if token == "POINT_DATA":
            data_kind = "point"
            data_count = int(parts[1])
            index += 1
            continue
        if token == "CELL_DATA":
            data_kind = "cell"
            data_count = int(parts[1])
            index += 1
            continue
        if token == "SCALARS" and data_kind:
            name = parts[1]
            components = int(parts[3]) if len(parts) > 3 else 1
            index += 1
            if lines[index].startswith("LOOKUP_TABLE"):
                index += 1
            values, index = _read_values(lines, index, data_count * components)
            values = values.reshape(data_count, components)
            target = point_data if data_kind == "point" else cell_data
            target[name] = values[:, 0] if components == 1 else values
            continue
        if token == "VECTORS" and data_kind:
            name = parts[1]
            values, index = _read_values(lines, index + 1, data_count * 3)
            target = point_data if data_kind == "point" else cell_data
            target[name] = values.reshape(data_count, 3)
            continue
        index += 1
    if points is None or triangles is None:
        raise ValueError(f"Missing mesh in {path}")
    return MeshState(points, triangles, point_data, cell_data)


def find_vtk(root: Path, variant: str, bias: float) -> Path:
    label = f"-{bias:g}V"
    matches = sorted((root / "variants" / variant / "run" / "vtk").rglob(f"dc_sweep_*_{label}.vtk"))
    if len(matches) != 1:
        raise FileNotFoundError(f"Expected one VTK for {variant} at {label}, found {len(matches)}")
    return matches[0]


def point_magnitude(state: MeshState, name: str) -> np.ndarray:
    values = state.point_data[name]
    return np.linalg.norm(values[:, :2], axis=1) if values.ndim == 2 else np.abs(values)


def positive_log(values: np.ndarray) -> np.ndarray:
    absolute = np.abs(np.asarray(values, dtype=float))
    positive = absolute[np.isfinite(absolute) & (absolute > 0)]
    floor = max(float(positive.min()) * 1e-3, 1e-300) if positive.size else 1e-300
    return np.log10(np.maximum(absolute, floor))


def shared_limits(grids: list[list[np.ndarray]]) -> list[tuple[float, float]]:
    limits = []
    for column in zip(*grids):
        merged = np.concatenate([values[np.isfinite(values)] for values in column])
        low, high = np.nanpercentile(merged, [0.0, 100.0])
        if math.isclose(low, high):
            high = low + 1.0
        limits.append((float(low), float(high)))
    return limits


def draw_heatmap_grid(
    states: list[MeshState],
    row_labels: list[str],
    data: list[list[np.ndarray]],
    column_labels: list[str],
    out_dir: Path,
    stem: str,
    title: str,
) -> list[Path]:
    rows, columns = len(states), len(column_labels)
    limits = shared_limits(data)
    fig, axes = plt.subplots(
        rows,
        columns,
        figsize=(3.0 * columns, 1.55 * rows + 0.8),
        constrained_layout=True,
        squeeze=False,
    )
    artists: list[list] = [[] for _ in range(columns)]
    for row, state in enumerate(states):
        triangulation = mtri.Triangulation(state.points[:, 0], state.points[:, 1], state.triangles)
        for column in range(columns):
            ax = axes[row, column]
            values = data[row][column]
            artist = ax.tripcolor(
                triangulation,
                values,
                shading="gouraud",
                cmap="viridis",
                vmin=limits[column][0],
                vmax=limits[column][1],
            )
            ax.triplot(triangulation, color="#333333", linewidth=0.18, alpha=0.45)
            artists[column].append(artist)
            ax.set_aspect("equal", adjustable="box")
            ax.set_xlabel("x (um)" if row == rows - 1 else "")
            ax.set_ylabel("y (um)" if column == 0 else "")
            if row == 0:
                ax.set_title(column_labels[column])
            if column == 0:
                ax.text(
                    -0.38,
                    0.5,
                    row_labels[row],
                    transform=ax.transAxes,
                    ha="right",
                    va="center",
                    rotation=90,
                    fontsize=7.5,
                )
    for column in range(columns):
        colorbar = fig.colorbar(artists[column][0], ax=axes[:, column], shrink=0.82, pad=0.02)
        colorbar.ax.set_ylabel("log10 magnitude", rotation=90)
    fig.suptitle(title, fontweight="normal")
    return save_figure(fig, out_dir, stem)


def plot_six_variant_fields(root: Path, out_dir: Path) -> list[Path]:
    states = [load_legacy_vtk(find_vtk(root, variant, 19.0)) for variant in VARIANTS]
    data = []
    for state in states:
        data.append(
            [
                positive_log(point_magnitude(state, "ElectricField") * 1e4),
                positive_log(state.point_data["Electrons"]),
                positive_log(state.point_data["ElectronAlphaAvalanche"] * 1e4),
                positive_log(state.point_data["AvalancheGeneration"]),
            ]
        )
    return draw_heatmap_grid(
        states,
        [DISPLAY_NAMES[variant] for variant in VARIANTS],
        data,
        ["|E| (V/cm)", "n (cm^-3)", "alpha_n (cm^-1)", "G_ava (cm^-3 s^-1)"],
        out_dir,
        "pn2d-bv-six-variant-fields-19v",
        "Six-variant physical-field comparison at -19 V",
    )


def plot_triangle_critical(root: Path, out_dir: Path) -> list[Path]:
    pairs = [
        ("legacy_triangle_gss_gradqf", 19.0),
        ("legacy_triangle_gss_gradqf", 19.4),
        ("reported_triangle_gss_gradqf", 19.0),
        ("reported_triangle_gss_gradqf", 19.4),
    ]
    states = [load_legacy_vtk(find_vtk(root, variant, bias)) for variant, bias in pairs]
    data = []
    for state in states:
        data.append(
            [
                positive_log(state.point_data["ElectronImpactIonizationDrive"] * 1e4),
                positive_log(state.point_data["ElectronAlphaAvalanche"] * 1e4),
                positive_log(state.point_data["AvalancheGeneration"]),
                positive_log(point_magnitude(state, "ElectronCurrentDensityVector")),
            ]
        )
    labels = [f"{variant.split('_')[0].title()} / -{bias:g} V" for variant, bias in pairs]
    return draw_heatmap_grid(
        states,
        labels,
        data,
        ["|grad phi_Fn| (V/cm)", "alpha_n (cm^-1)", "G_ava (cm^-3 s^-1)", "|J_n| (A/cm^2)"],
        out_dir,
        "pn2d-bv-triangle-gss-critical-fields",
        "Triangle-GSS field evolution near the convergence boundary",
    )


def load_sentaurus_mesh(export: Path) -> MeshState:
    node_rows = read_numeric_csv(export / "nodes.csv")
    node_ids = [int(row["id"]) for row in node_rows]
    index_by_id = {node_id: index for index, node_id in enumerate(node_ids)}
    points = np.asarray([[float(row["x_um"]), float(row["y_um"])] for row in node_rows])
    triangles = []
    for row in read_numeric_csv(export / "elements.csv"):
        triangles.append([index_by_id[int(row[f"node{i}"])] for i in range(3)])
    point_data: dict[str, np.ndarray] = {}
    for name in ("eDensity", "eCurrentDensity", "eAlphaAvalanche", "ElectricField"):
        rows = read_numeric_csv(export / "fields" / f"{name}_region0.csv")
        components = sorted(key for key in rows[0] if key.startswith("component"))
        values = np.full((len(node_ids), len(components)), np.nan)
        for row in rows:
            values[index_by_id[int(row["node_id"])] ] = [float(row[key]) for key in components]
        point_data[name] = values[:, 0] if len(components) == 1 else values
    return MeshState(points, np.asarray(triangles, dtype=int), point_data, {})


def plot_sentaurus_vela_fields(root: Path, out_dir: Path) -> list[Path]:
    sentaurus = load_sentaurus_mesh(root / "sentaurus_exports" / "sentaurus_-19v")
    variants = (
        "reported_density_gradient",
        "reported_gss_midpoint",
        "reported_triangle_gss_gradqf",
    )
    states = [sentaurus] + [load_legacy_vtk(find_vtk(root, variant, 19.0)) for variant in variants]
    data = [
        [
            positive_log(sentaurus.point_data["eDensity"]),
            positive_log(point_magnitude(sentaurus, "eCurrentDensity")),
            positive_log(sentaurus.point_data["eAlphaAvalanche"]),
            positive_log(point_magnitude(sentaurus, "ElectricField")),
        ]
    ]
    for state in states[1:]:
        data.append(
            [
                positive_log(state.point_data["Electrons"]),
                positive_log(point_magnitude(state, "ElectronCurrentDensityVector")),
                positive_log(state.point_data["ElectronAlphaAvalanche"] * 1e4),
                positive_log(point_magnitude(state, "ElectricField") * 1e4),
            ]
        )
    return draw_heatmap_grid(
        states,
        ["Sentaurus 2018"] + [DISPLAY_NAMES[variant] for variant in variants],
        data,
        ["n (cm^-3)", "|J_n| (A/cm^2)", "alpha_n (cm^-1)", "|E| (V/cm)"],
        out_dir,
        "pn2d-bv-sentaurus-vela-fields-19v",
        "Sentaurus-Vela physical-field comparison at -19 V",
    )


def write_visual_html(path: Path, png_paths: list[Path]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    captions = {
        "pn2d-bv-iv-comparison.png": "Current-voltage comparison and breakdown-region detail",
        "pn2d-bv-six-variant-fields-19v.png": "Six-variant physical fields at -19 V",
        "pn2d-bv-triangle-gss-critical-fields.png": "Triangle-GSS evolution from -19 V to -19.4 V",
        "pn2d-bv-sentaurus-vela-fields-19v.png": "Sentaurus-Vela field comparison at -19 V",
    }
    blocks = []
    for index, png_path in enumerate(png_paths):
        preview_path = png_path.with_name(f"{png_path.stem}-preview.jpg")
        with Image.open(png_path) as image:
            image = image.convert("RGB")
            image.thumbnail((1400, 1400), Image.Resampling.LANCZOS)
            image.save(preview_path, quality=82, optimize=True)
        encoded = base64.b64encode(preview_path.read_bytes()).decode("ascii")
        blocks.append(
            f'<figure><img src="data:image/jpeg;base64,{encoded}" '
            f'alt="{captions[png_path.name]}" loading="{"eager" if index == 0 else "lazy"}">'
            f'<figcaption class="text-small text-muted">{captions[png_path.name]}</figcaption></figure>'
        )
    markup = "\n".join(
        [
            '<div id="pn2d-bv-static-figures" class="pn2d-bv-static-figures">',
            "<style>",
            "#pn2d-bv-static-figures { display: grid; gap: 1rem; width: 100%; }",
            "#pn2d-bv-static-figures figure { margin: 0; }",
            "#pn2d-bv-static-figures img { display: block; width: 100%; height: auto; }",
            "#pn2d-bv-static-figures figcaption { margin-top: 0.35rem; color: var(--muted-foreground); }",
            "</style>",
            *blocks,
            "</div>",
        ]
    )
    path.write_text(markup + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    root = args.root.resolve()
    out_dir = (args.out_dir or (root / "figures_static")).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    configure_matplotlib()

    generated: list[Path] = []
    iv_outputs, endpoints = plot_iv(root, args.sentaurus_reference.resolve(), out_dir)
    generated.extend(iv_outputs)
    generated.extend(plot_six_variant_fields(root, out_dir))
    generated.extend(plot_triangle_critical(root, out_dir))
    generated.extend(plot_sentaurus_vela_fields(root, out_dir))

    manifest = {
        "report_root": str(root),
        "sentaurus_reference": str(args.sentaurus_reference.resolve()),
        "figures": [str(path) for path in generated],
        "iv_endpoints": endpoints,
        "notes": [
            "All heatmaps use base-10 logarithms of field magnitudes.",
            "Vela electric field and alpha are converted from per-um internal units to per-cm units.",
            "Triangle-GSS curves stop at the last converged bias and are not extrapolated.",
        ],
    }
    (out_dir / "figure_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=True) + "\n", encoding="utf-8"
    )
    if args.visual_html:
        write_visual_html(args.visual_html.resolve(), [path for path in generated if path.suffix == ".png"])
    print(json.dumps(manifest, indent=2, ensure_ascii=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
