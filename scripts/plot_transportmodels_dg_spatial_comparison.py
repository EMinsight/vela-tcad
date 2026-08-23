#!/usr/bin/env python3
"""Compare Sentaurus and Vela DG spatial fields at Vg=1 V, Vd=2 V."""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np
from matplotlib.colors import Normalize, TwoSlopeNorm


REPO_ROOT = Path(__file__).resolve().parents[1]
BASELINE_DIR = (
    REPO_ROOT
    / "build-release/reference_tcad/transportmodels_sentaurus2022/vela_baseline"
)
DEFAULT_STATE = (
    BASELINE_DIR
    / "workflow_dg_idvd_resume_1p9_direct2p0_fixed1p5_outer200_run01"
    / "dg_idvd_curve_state_bias_2p000000.csv"
)
DEFAULT_REFERENCE_FIELDS = BASELINE_DIR / "generated/sim_fields/dg_idvd/fields"
DEFAULT_OUTPUT_DIR = (
    REPO_ROOT
    / "docs/progress_report_2026Q3/2026-08-19_transportmodels_dg_daily_report"
    / "figures/spatial_fields"
)


@dataclass(frozen=True)
class FieldSpec:
    key: str
    title: str
    sentaurus_file: str
    vela_column: str
    absolute_label: str
    difference_label: str
    logarithmic: bool = False
    vela_scale: float = 1.0
    difference_scale: float = 1.0


FIELD_SPECS = (
    FieldSpec(
        "electrostatic_potential",
        "静电势 ψ",
        "ElectrostaticPotential_region3.csv",
        "psi",
        "ψ  (V)",
        "Vela − Sentaurus  (mV)",
        difference_scale=1.0e3,
    ),
    FieldSpec(
        "electron_quasi_fermi_potential",
        "电子准费米势 φn",
        "eQuasiFermiPotential_region3.csv",
        "phin",
        "φn  (V)",
        "Vela − Sentaurus  (mV)",
        difference_scale=1.0e3,
    ),
    FieldSpec(
        "hole_quasi_fermi_potential",
        "空穴准费米势 φp",
        "hQuasiFermiPotential_region3.csv",
        "phip",
        "φp  (V)",
        "Vela − Sentaurus  (mV)",
        difference_scale=1.0e3,
    ),
    FieldSpec(
        "electron_density",
        "电子浓度 n",
        "eDensity_region3.csv",
        "electrons_m3",
        "log10(n / cm⁻³)",
        "log10(Vela / Sentaurus)  (dex)",
        logarithmic=True,
        vela_scale=1.0e-6,
    ),
    FieldSpec(
        "hole_density",
        "空穴浓度 p",
        "hDensity_region3.csv",
        "holes_m3",
        "log10(p / cm⁻³)",
        "log10(Vela / Sentaurus)  (dex)",
        logarithmic=True,
        vela_scale=1.0e-6,
    ),
    FieldSpec(
        "electron_quantum_potential",
        "电子量子势 Qn",
        "eQuantumPotential_region3.csv",
        "electron_quantum_potential_V",
        "Qn  (V；Sentaurus eV/q)",
        "Vela − Sentaurus  (mV)",
        difference_scale=1.0e3,
    ),
)


def configure_style() -> None:
    plt.rcParams.update(
        {
            "font.family": ["Microsoft YaHei", "DejaVu Sans"],
            "axes.unicode_minus": False,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "axes.edgecolor": "#343A40",
            "axes.labelcolor": "#252A31",
            "axes.titlecolor": "#20252B",
            "xtick.color": "#4E5661",
            "ytick.color": "#4E5661",
            "savefig.facecolor": "white",
            "savefig.bbox": "tight",
        }
    )


def read_scalar_csv(path: Path, value_column: str) -> dict[int, float]:
    values: dict[int, float] = {}
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        for row in csv.DictReader(stream):
            node_id = int(row["node_id"])
            if node_id in values:
                raise RuntimeError(f"Duplicate node_id={node_id} in {path}")
            values[node_id] = float(row[value_column])
    return values


def load_geometry(mesh_path: Path) -> tuple[
    np.ndarray, np.ndarray, np.ndarray, np.ndarray, list[int]
]:
    with mesh_path.open("r", encoding="utf-8") as stream:
        mesh = json.load(stream)
    nodes = sorted(mesh["nodes"], key=lambda row: int(row["id"]))
    node_ids = [int(node["id"]) for node in nodes]
    if node_ids != list(range(len(node_ids))):
        raise RuntimeError("Spatial comparison requires contiguous mesh node IDs")

    # Sentaurus x is surface-normal depth and y is lateral position.
    lateral = np.asarray([float(node["y"]) for node in nodes])
    depth = np.asarray([float(node["x"]) for node in nodes])
    substrate_triangles = np.asarray(
        [
            [int(node_id) for node_id in triangle["node_ids"]]
            for triangle in mesh["triangles"]
            if int(triangle["region_id"]) == 3
        ],
        dtype=int,
    )
    substrate_node_ids = sorted(set(substrate_triangles.ravel().tolist()))
    return lateral, depth, substrate_triangles, np.asarray(substrate_node_ids), node_ids


def dense_values(
    values: dict[int, float], node_count: int, required_node_ids: np.ndarray
) -> np.ndarray:
    missing = [int(node_id) for node_id in required_node_ids if int(node_id) not in values]
    if missing:
        raise RuntimeError(
            f"Field is missing {len(missing)} substrate nodes; first={missing[:5]}"
        )
    result = np.zeros(node_count, dtype=float)
    for node_id, value in values.items():
        if not np.isfinite(value):
            raise RuntimeError(f"Non-finite value at node_id={node_id}")
        result[node_id] = value
    return result


def robust_symmetric_limit(values: np.ndarray) -> float:
    finite = np.abs(values[np.isfinite(values)])
    if not len(finite):
        return 1.0
    limit = float(np.percentile(finite, 99.0))
    if limit <= 0.0:
        limit = float(np.max(finite))
    return limit if limit > 0.0 else 1.0


def render_field(
    output_dir: Path,
    spec: FieldSpec,
    triangulation: mtri.Triangulation,
    lateral: np.ndarray,
    depth: np.ndarray,
    substrate_node_ids: np.ndarray,
    reference: np.ndarray,
    vela: np.ndarray,
    vela_label: str,
) -> dict[str, float | int | str]:
    visible = substrate_node_ids[depth[substrate_node_ids] <= 0.22]
    reference_visible = reference[visible]
    vela_visible = vela[visible]

    if spec.logarithmic:
        if np.any(reference_visible <= 0.0) or np.any(vela_visible <= 0.0):
            raise RuntimeError(f"{spec.key} contains nonpositive carrier density")
        reference_display = np.zeros_like(reference)
        vela_display = np.zeros_like(vela)
        reference_display[substrate_node_ids] = np.log10(
            reference[substrate_node_ids]
        )
        vela_display[substrate_node_ids] = np.log10(vela[substrate_node_ids])
        difference = vela_display - reference_display
    else:
        reference_display = reference
        vela_display = vela
        difference = spec.difference_scale * (vela - reference)

    shared_values = np.concatenate(
        (reference_display[visible], vela_display[visible])
    )
    shared_norm = Normalize(
        vmin=float(np.min(shared_values)), vmax=float(np.max(shared_values))
    )
    difference_limit = robust_symmetric_limit(difference[visible])
    difference_norm = TwoSlopeNorm(
        vmin=-difference_limit, vcenter=0.0, vmax=difference_limit
    )

    fig, axes = plt.subplots(1, 3, figsize=(14.2, 5.2), constrained_layout=True)
    absolute_mappable = None
    difference_mappable = None
    for index, (axis, title, values, cmap, norm) in enumerate(
        (
            (axes[0], "Sentaurus T-2022.03-SP2", reference_display, "viridis", shared_norm),
            (axes[1], vela_label, vela_display, "viridis", shared_norm),
            (axes[2], "差值", difference, "RdBu_r", difference_norm),
        )
    ):
        mappable = axis.tripcolor(
            triangulation,
            values,
            shading="gouraud",
            cmap=cmap,
            norm=norm,
            rasterized=True,
        )
        axis.triplot(
            triangulation,
            color="#5D6670",
            linewidth=0.12,
            alpha=0.28,
        )
        axis.axhline(0.0, color="#343A40", linewidth=0.9)
        axis.set_title(title, fontsize=12, pad=8)
        axis.set_xlabel("横向位置  (μm)", fontsize=10)
        if index == 0:
            axis.set_ylabel("表面法向深度  (μm)", fontsize=10)
        axis.set_xlim(-0.2, 0.2)
        axis.set_ylim(0.22, -0.005)
        axis.set_aspect("equal", adjustable="box")
        axis.grid(False)
        if index < 2:
            absolute_mappable = mappable
        else:
            difference_mappable = mappable

    assert absolute_mappable is not None and difference_mappable is not None
    absolute_colorbar = fig.colorbar(
        absolute_mappable, ax=axes[:2], location="bottom", shrink=0.82, pad=0.08
    )
    absolute_colorbar.set_label(spec.absolute_label, fontsize=10)
    difference_colorbar = fig.colorbar(
        difference_mappable, ax=axes[2], location="bottom", shrink=0.82, pad=0.08
    )
    difference_colorbar.set_label(spec.difference_label, fontsize=10)

    absolute_difference = np.abs(difference[visible])
    fig.suptitle(
        f"TransportModels DG 空间场对比：{spec.title}", fontsize=16, y=1.035
    )
    fig.text(
        0.01,
        -0.015,
        (
            "工作点：Vg = 1.0 V，Vd = 2.0 V；区域：R.Substrate 表面 0–0.22 μm。"
            f"节点数={len(visible)}；差值色标按 99% 分位截断，实际最大绝对差值="
            f"{float(np.max(absolute_difference)):.4g} {spec.difference_label.split('(')[-1].rstrip(')')}。"
        ),
        fontsize=9.0,
        color="#626A74",
    )

    stem = f"transportmodels_dg_{spec.key}_spatial_comparison"
    for suffix in ("png", "svg"):
        fig.savefig(
            output_dir / f"{stem}.{suffix}",
            dpi=240 if suffix == "png" else None,
        )
    plt.close(fig)

    return {
        "nodes_compared": int(len(visible)),
        "absolute_min": float(np.min(shared_values)),
        "absolute_max": float(np.max(shared_values)),
        "difference_median_absolute": float(np.median(absolute_difference)),
        "difference_rmse": float(np.sqrt(np.mean(np.square(difference[visible])))),
        "difference_max_absolute": float(np.max(absolute_difference)),
        "difference_color_limit_99pct": difference_limit,
        "absolute_unit": spec.absolute_label,
        "difference_unit": spec.difference_label,
        "png": f"{stem}.png",
        "svg": f"{stem}.svg",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--state", type=Path, default=DEFAULT_STATE)
    parser.add_argument(
        "--reference-fields", type=Path, default=DEFAULT_REFERENCE_FIELDS
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--summary-json", type=Path)
    parser.add_argument("--vela-label", default="Vela DG")
    args = parser.parse_args()

    state_path = args.state.resolve()
    reference_fields = args.reference_fields.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    configure_style()

    mesh_path = (BASELINE_DIR / "generated/vela/mesh.json").resolve()
    lateral, depth, triangles, substrate_node_ids, node_ids = load_geometry(mesh_path)
    triangulation = mtri.Triangulation(lateral, depth, triangles)
    state_rows = read_scalar_csv(state_path, "psi")
    if sorted(state_rows) != node_ids:
        raise RuntimeError("Vela state node IDs do not exactly cover the mesh")

    metrics: dict[str, dict[str, float | int | str]] = {}
    for spec in FIELD_SPECS:
        reference_map = read_scalar_csv(
            reference_fields / spec.sentaurus_file, "component0"
        )
        vela_map = read_scalar_csv(state_path, spec.vela_column)
        reference = dense_values(reference_map, len(node_ids), substrate_node_ids)
        vela = spec.vela_scale * dense_values(
            vela_map, len(node_ids), substrate_node_ids
        )
        metrics[spec.key] = render_field(
            output_dir,
            spec,
            triangulation,
            lateral,
            depth,
            substrate_node_ids,
            reference,
            vela,
            args.vela_label,
        )

    summary = {
        "scope": {
            "device": "TransportModels DG MOS",
            "gate_bias_V": 1.0,
            "drain_bias_V": 2.0,
            "region": "R.Substrate",
            "depth_window_um": [0.0, 0.22],
            "comparison_grain": "shared global mesh node_id",
            "vela_label": args.vela_label,
        },
        "sources": {
            "mesh": str(mesh_path),
            "vela_state": str(state_path),
            "sentaurus_fields": str(reference_fields),
        },
        "unit_policy": {
            "carrier_density": "Sentaurus cm^-3; Vela m^-3 converted to cm^-3",
            "electron_quantum_potential": (
                "Sentaurus eV divided by elementary charge is numerically compared "
                "with Vela volts"
            ),
        },
        "fields": metrics,
        "chart_map": [
            {
                "segment": spec.key,
                "question": f"Where does Vela reproduce or differ from Sentaurus for {spec.title}?",
                "family": "spatial matrix and benchmark",
                "type": "three-panel triangular field map with shared absolute scale and signed difference",
                "fields": [spec.sentaurus_file, spec.vela_column, "node_id"],
                "palette_policy": "single-root absolute maps plus one diverging difference root",
                "outputs": [metrics[spec.key]["png"], metrics[spec.key]["svg"]],
            }
            for spec in FIELD_SPECS
        ],
    }
    summary_path = (
        args.summary_json
        or output_dir / "transportmodels_dg_spatial_comparison_summary.json"
    ).resolve()
    summary_path.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
