#!/usr/bin/env python3
"""Compare matched Sentaurus/Vela DD/DG fields in three Id-Vg regimes."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from pathlib import Path
from typing import Any, Iterable


REPO_ROOT = Path(__file__).resolve().parents[1]
SENTAURUS_MANIFEST = (
    REPO_ROOT
    / "build-release/reference_tcad/transportmodels_sentaurus2022/sentaurus_vm_runs"
    / "three_regime_spatial_oracles_20260824/three_regime_spatial_oracles_manifest.json"
)
VELA_MANIFEST = (
    REPO_ROOT
    / "build-release/reference_tcad/transportmodels_sentaurus2022/vela_baseline"
    / "three_regime_spatial_vtk_2026-08-24/three_regime_spatial_manifest.json"
)
OUTPUT_ROOT = (
    REPO_ROOT
    / "build-release/reference_tcad/transportmodels_sentaurus2022/reports"
    / "transportmodels_three_regime_spatial_20260824"
)
REPORT_JSON = REPO_ROOT / "docs/validation/transportmodels_three_regime_spatial_2026-08-24.json"
REPORT_MD = REPO_ROOT / "docs/validation/transportmodels_three_regime_spatial_2026-08-24.md"

VECTOR_FIELDS = (
    ("P0", "current", "Electron current density", "eCurrentDensity", "SentaurusElectronCurrentDensityVector", "A/cm2"),
    ("P0", "current", "Hole current density", "hCurrentDensity", "SentaurusHoleCurrentDensityVector", "A/cm2"),
    ("P0", "current", "Total current density", "TotalCurrentDensity", "SentaurusTotalCurrentDensityVector", "A/cm2"),
    ("P1", "transport", "Electric field", "ElectricField", "ElectricFieldVector", "V/cm"),
    ("P1", "transport", "Electron GradQuasiFermi", "eGradQuasiFermi", "ElectronGradQuasiFermiVector", "V/cm"),
    ("P1", "transport", "Hole GradQuasiFermi", "hGradQuasiFermi", "HoleGradQuasiFermiVector", "V/cm"),
)
SCALAR_FIELDS = (
    ("P1", "transport", "Electron Eparallel", "eEparallel", "SentaurusElectronEparallel", "V/cm", "magnitude", 1.0),
    ("P1", "transport", "Hole Eparallel", "hEparallel", "SentaurusHoleEparallel", "V/cm", "magnitude", 1.0),
    ("P1", "transport", "Electron Enormal", "eEnormal", "ElectronEnormal", "V/cm", "magnitude", 1.0),
    ("P1", "transport", "Hole Enormal", "hEnormal", "HoleEnormal", "V/cm", "magnitude", 1.0),
    ("P1", "transport", "Electron mobility", "eMobility", "ElectronMobilityCm2PerVs", "cm2/(V s)", "magnitude", 1.0),
    ("P1", "transport", "Hole mobility", "hMobility", "HoleMobilityCm2PerVs", "cm2/(V s)", "magnitude", 1.0),
    ("P1", "source_band", "SRH recombination", "srhRecombination", "SRHRecombinationCm3PerS", "cm-3 s-1", "signed", 1.0),
    ("P1", "source_band", "Space charge", "SpaceCharge", "SpaceCharge", "q/cm3", "signed", 1.0),
    ("P1", "source_band", "Band gap", "BandGap", "BandGap", "eV", "magnitude", 1.0),
    ("P1", "source_band", "Bandgap narrowing", "BandgapNarrowing", "BandgapNarrowing", "eV", "magnitude", 1.0),
    ("P1", "source_band", "Electron affinity", "ElectronAffinity", "ElectronAffinity", "eV", "magnitude", 1.0),
    ("P1", "source_band", "Conduction band", "ConductionBandEnergy", "ConductionBandEnergy", "eV", "signed", 1.0),
    ("P1", "source_band", "Valence band", "ValenceBandEnergy", "ValenceBandEnergy", "eV", "signed", 1.0),
)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise RuntimeError(f"Refusing to write empty CSV: {path}")
    fieldnames = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def percentile(values: Iterable[float], fraction: float) -> float:
    ordered = sorted(float(value) for value in values if math.isfinite(float(value)))
    if not ordered:
        return math.nan
    position = (len(ordered) - 1) * fraction
    low, high = math.floor(position), math.ceil(position)
    if low == high:
        return ordered[low]
    return ordered[low] * (high - position) + ordered[high] * (position - low)


def read_scalar(path: Path) -> dict[int, float]:
    return {int(row["node_id"]): float(row["component0"]) for row in read_csv(path)}


def read_vector(path: Path) -> dict[int, tuple[float, float]]:
    return {
        int(row["node_id"]): (float(row["component0"]), float(row["component1"]))
        for row in read_csv(path)
    }


def read_nodes(path: Path) -> dict[int, tuple[float, float]]:
    return {
        int(row["id"]): (float(row["x_um"]), float(row["y_um"]))
        for row in read_csv(path)
    }


def silicon_interior_nodes(export_dir: Path) -> set[int]:
    regions: dict[int, set[str]] = {}
    for element in read_csv(export_dir / "elements.csv"):
        for local in range(3):
            regions.setdefault(int(element[f"node{local}"]), set()).add(element["region"])
    return {node for node, names in regions.items() if names == {"R.Substrate"}}


def parse_vtk(path: Path) -> tuple[list[tuple[float, float]], dict[str, list[Any]]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    points_index = next(i for i, line in enumerate(lines) if line.startswith("POINTS "))
    point_count = int(lines[points_index].split()[1])
    points = [
        tuple(float(value) for value in lines[index].split()[:2])
        for index in range(points_index + 1, points_index + 1 + point_count)
    ]
    point_data = next(i for i, line in enumerate(lines) if line.startswith("POINT_DATA "))
    count = int(lines[point_data].split()[1])
    fields: dict[str, list[Any]] = {}
    index = point_data + 1
    while index < len(lines):
        parts = lines[index].split()
        if not parts:
            index += 1
        elif parts[0] == "SCALARS":
            name = parts[1]
            index += 2
            fields[name] = [float(lines[index + offset].split()[0]) for offset in range(count)]
            index += count
        elif parts[0] == "VECTORS":
            name = parts[1]
            index += 1
            fields[name] = [
                tuple(float(value) for value in lines[index + offset].split()[:3])
                for offset in range(count)
            ]
            index += count
        else:
            index += 1
    return points, fields


def magnitude(vector: tuple[float, ...]) -> float:
    return math.hypot(vector[0], vector[1])


def symmetric_percent(reference: float, candidate: float, floor: float) -> float:
    return 200.0 * abs(candidate - reference) / max(abs(candidate) + abs(reference), floor)


def vector_angle_deg(reference: tuple[float, float], candidate: tuple[float, float]) -> float:
    ref_norm, cand_norm = magnitude(reference), magnitude(candidate)
    if ref_norm == 0.0 or cand_norm == 0.0:
        return math.nan
    cosine = (reference[0] * candidate[0] + reference[1] * candidate[1]) / (ref_norm * cand_norm)
    return math.degrees(math.acos(max(-1.0, min(1.0, cosine))))


def vector_summary(rows: list[dict[str, Any]]) -> dict[str, Any]:
    peak = max(abs(float(row["sentaurus_magnitude"])) for row in rows)
    active_floor = max(peak * 1.0e-3, 1.0e-300)
    active = [row for row in rows if abs(float(row["sentaurus_magnitude"])) >= active_floor]
    for row in rows:
        row["symmetric_magnitude_error_percent"] = symmetric_percent(
            float(row["sentaurus_magnitude"]), float(row["vela_magnitude"]), active_floor)
        row["direction_error_deg"] = vector_angle_deg(
            (float(row["sentaurus_x"]), float(row["sentaurus_y"])),
            (float(row["vela_x"]), float(row["vela_y"])))
        row["active"] = row in active
    percentages = [float(row["symmetric_magnitude_error_percent"]) for row in active]
    angles = [float(row["direction_error_deg"]) for row in active if math.isfinite(float(row["direction_error_deg"]))]
    return {
        "nodes": len(rows), "active_nodes": len(active), "active_floor": active_floor,
        "primary_metric": "active-node symmetric magnitude error percent",
        "p50_percent": percentile(percentages, 0.50),
        "p95_percent": percentile(percentages, 0.95),
        "max_percent": max(percentages) if percentages else math.nan,
        "direction_p95_deg": percentile(angles, 0.95),
    }


def scalar_summary(rows: list[dict[str, Any]], metric_kind: str) -> dict[str, Any]:
    reference_abs = [abs(float(row["sentaurus_value"])) for row in rows]
    peak = max(reference_abs)
    active_floor = max(peak * 1.0e-3, 1.0e-300)
    active = [row for row in rows if abs(float(row["sentaurus_value"])) >= active_floor]
    if metric_kind == "magnitude":
        for row in rows:
            row["sentaurus_value"] = abs(float(row["sentaurus_value"]))
            row["vela_value"] = abs(float(row["vela_value"]))
            row["error_percent"] = symmetric_percent(
                float(row["sentaurus_value"]), float(row["vela_value"]), active_floor)
        values = [float(row["error_percent"]) for row in active]
        label = "active-node symmetric error percent"
    else:
        scale = max(percentile(reference_abs, 0.95), active_floor)
        for row in rows:
            row["error_percent"] = 100.0 * abs(
                float(row["vela_value"]) - float(row["sentaurus_value"])) / scale
        values = [float(row["error_percent"]) for row in rows]
        label = "p95-reference-scale normalized absolute error percent"
    for row in rows:
        row["active"] = row in active
    return {
        "nodes": len(rows), "active_nodes": len(active), "active_floor": active_floor,
        "primary_metric": label,
        "p50_percent": percentile(values, 0.50),
        "p95_percent": percentile(values, 0.95),
        "max_percent": max(values) if values else math.nan,
    }


def compare_case(sent_state: dict[str, Any], vela_state: dict[str, Any]) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    export_dir = Path(sent_state["export_dir"])
    vtk_path = Path(vela_state["vtk"])
    nodes = read_nodes(export_dir / "nodes.csv")
    vtk_points, vtk = parse_vtk(vtk_path)
    if len(nodes) != len(vtk_points) or any(
        max(abs(nodes[node][axis] - vtk_points[node][axis]) for axis in (0, 1)) > 1.0e-12
        for node in nodes
    ):
        raise RuntimeError(f"Non-exact Sentaurus/Vela node mapping for {sent_state['mode']}/{sent_state['regime']}")
    interior = silicon_interior_nodes(export_dir)
    fields_dir = export_dir / "fields"
    long_rows: list[dict[str, Any]] = []
    summaries: list[dict[str, Any]] = []
    common = {
        "mode": sent_state["mode"], "regime": sent_state["regime"],
        "gate_bias_V": sent_state["gate_bias_V"], "drain_bias_V": sent_state["drain_bias_V"],
    }
    for priority, category, label, sent_name, vela_name, unit in VECTOR_FIELDS:
        sent = read_vector(fields_dir / f"{sent_name}_region3.csv")
        candidate = vtk[vela_name]
        rows = []
        for node in sorted(interior & sent.keys()):
            sx, sy = sent[node]
            vx, vy = float(candidate[node][0]), float(candidate[node][1])
            rows.append({
                **common, "priority": priority, "category": category, "field": label,
                "unit": unit, "node_id": node, "x_um": nodes[node][0], "y_um": nodes[node][1],
                "sentaurus_x": sx, "sentaurus_y": sy, "vela_x": vx, "vela_y": vy,
                "sentaurus_magnitude": math.hypot(sx, sy), "vela_magnitude": math.hypot(vx, vy),
            })
        summary = vector_summary(rows)
        long_rows.extend(rows)
        summaries.append({**common, "priority": priority, "category": category, "field": label, "unit": unit, **summary})

    scalar_cache: dict[str, dict[int, float]] = {}
    for _, _, _, sent_name, _, _, _, _ in SCALAR_FIELDS:
        scalar_cache.setdefault(sent_name, read_scalar(fields_dir / f"{sent_name}_region3.csv"))
    # Sentaurus Ec/Ev carry one arbitrary global origin; align both with the Ec offset.
    sent_ec = scalar_cache["ConductionBandEnergy"]
    vela_ec = vtk["ConductionBandEnergy"]
    band_origin = statistics.median(
        sent_ec[node] - float(vela_ec[node]) for node in sorted(interior & sent_ec.keys())
    )
    for priority, category, label, sent_name, vela_name, unit, metric_kind, scale in SCALAR_FIELDS:
        sent = scalar_cache[sent_name]
        candidate = vtk[vela_name]
        rows = []
        for node in sorted(interior & sent.keys()):
            sent_value = sent[node]
            if sent_name in {"ConductionBandEnergy", "ValenceBandEnergy"}:
                sent_value -= band_origin
            rows.append({
                **common, "priority": priority, "category": category, "field": label,
                "unit": unit, "node_id": node, "x_um": nodes[node][0], "y_um": nodes[node][1],
                "sentaurus_value": sent_value, "vela_value": float(candidate[node]) * scale,
                "band_origin_removed_eV": band_origin if sent_name in {"ConductionBandEnergy", "ValenceBandEnergy"} else "",
            })
        summary = scalar_summary(rows, metric_kind)
        long_rows.extend(rows)
        summaries.append({**common, "priority": priority, "category": category, "field": label, "unit": unit, **summary})
    return long_rows, summaries


def make_heatmap(summary: list[dict[str, Any]], category: str, filename: str, title: str) -> tuple[Path, Path]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    rows = [row for row in summary if row["category"] == category]
    case_order = [(mode, regime) for mode in ("dd", "dg") for regime in ("deep_off", "threshold", "on")]
    field_order = list(dict.fromkeys(row["field"] for row in rows))
    values = np.array([
        [next(float(row["p95_percent"]) for row in rows if row["mode"] == mode and row["regime"] == regime and row["field"] == field)
         for field in field_order]
        for mode, regime in case_order
    ])
    color_max = 100.0 if category == "source_band" else 200.0
    clipped = np.minimum(values, color_max)
    width = max(10.0, 1.45 * len(field_order))
    fig, axis = plt.subplots(figsize=(width, 5.8))
    image = axis.imshow(clipped, aspect="auto", cmap="YlOrRd", vmin=0.0, vmax=color_max)
    axis.set_xticks(range(len(field_order)), field_order, rotation=35, ha="right")
    axis.set_yticks(range(len(case_order)), [f"{mode.upper()} {regime}" for mode, regime in case_order])
    for i in range(values.shape[0]):
        for j in range(values.shape[1]):
            label = f">{color_max:.0f}" if values[i, j] > color_max else f"{values[i, j]:.1f}"
            axis.text(j, i, label, ha="center", va="center", fontsize=8,
                      color="white" if clipped[i, j] > 0.60 * color_max else "black")
    fig.colorbar(image, ax=axis, label="P95 error (%)")
    axis.set_title(title)
    fig.tight_layout()
    png, svg = OUTPUT_ROOT / f"{filename}.png", OUTPUT_ROOT / f"{filename}.svg"
    fig.savefig(png, dpi=190)
    fig.savefig(svg)
    plt.close(fig)
    return png, svg


def make_current_vector_figure(rows: list[dict[str, Any]]) -> tuple[Path, Path]:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    current = [row for row in rows if row["field"] == "Total current density"]
    fig, axes = plt.subplots(2, 3, figsize=(15.2, 7.8), sharex=True, sharey=True)
    for axis, mode, regime in zip(
        axes.flat,
        ("dd", "dd", "dd", "dg", "dg", "dg"),
        ("deep_off", "threshold", "on", "deep_off", "threshold", "on"),
    ):
        subset = [row for row in current if row["mode"] == mode and row["regime"] == regime]
        x = np.array([float(row["y_um"]) for row in subset])
        y = np.array([float(row["x_um"]) for row in subset])
        error = np.array([float(row["symmetric_magnitude_error_percent"]) for row in subset])
        scatter = axis.scatter(x, y, c=np.minimum(error, 200.0), s=6, cmap="YlOrRd", vmin=0, vmax=200)
        stride = max(1, len(subset) // 140)
        sample = subset[::stride]
        sx = np.array([float(row["sentaurus_y"]) for row in sample])
        sy = np.array([float(row["sentaurus_x"]) for row in sample])
        vx = np.array([float(row["vela_y"]) for row in sample])
        vy = np.array([float(row["vela_x"]) for row in sample])
        smag = np.hypot(sx, sy); vmag = np.hypot(vx, vy)
        sx, sy = sx / np.maximum(smag, 1e-300), sy / np.maximum(smag, 1e-300)
        vx, vy = vx / np.maximum(vmag, 1e-300), vy / np.maximum(vmag, 1e-300)
        axis.quiver([float(row["y_um"]) for row in sample], [float(row["x_um"]) for row in sample],
                    sx, sy, color="black", alpha=0.65, scale=32, width=0.0022, label="Sentaurus direction")
        axis.quiver([float(row["y_um"]) for row in sample], [float(row["x_um"]) for row in sample],
                    vx, vy, color="#0072B2", alpha=0.65, scale=32, width=0.0022, label="Vela direction")
        axis.invert_yaxis()
        axis.set_title(f"{mode.upper()} {regime}, Vg={subset[0]['gate_bias_V']:.2f} V")
        axis.grid(alpha=0.15)
    axes[0, 0].legend(loc="upper right", fontsize=7)
    for axis in axes[-1, :]: axis.set_xlabel("Lateral coordinate (um)")
    for axis in axes[:, 0]: axis.set_ylabel("Depth coordinate (um)")
    color_axis = fig.add_axes((0.905, 0.16, 0.014, 0.66))
    fig.colorbar(scatter, cax=color_axis, label="Total-current magnitude error (%)")
    fig.suptitle("Total current-density vectors: direction overlay and magnitude error")
    fig.subplots_adjust(left=0.07, right=0.88, bottom=0.09, top=0.90, wspace=0.18, hspace=0.22)
    png, svg = OUTPUT_ROOT / "total_current_density_vectors.png", OUTPUT_ROOT / "total_current_density_vectors.svg"
    fig.savefig(png, dpi=190)
    fig.savefig(svg)
    plt.close(fig)
    return png, svg


def markdown(report: dict[str, Any]) -> str:
    lines = [
        "# TransportModels DD/DG three-regime spatial-field comparison",
        "",
        "Status: **complete**. The comparison uses exact shared node IDs and only strict interior nodes of `R.Substrate`.",
        "",
        "Operating points: deep off `Vg=-1.00 V`, threshold/transition `Vg=0.12 V`, and on state `Vg=0.92 V`; all use `Vd=1.10 V`.",
        "",
        "Vector and positive-magnitude fields use active-node symmetric percentage error. Signed/cross-zero fields use absolute error normalized by the Sentaurus p95 absolute field scale. Ec/Ev are compared after removing the one arbitrary global Sentaurus energy origin.",
        "",
        "| Mode | Regime | Priority | Field | P95 error (%) | Max error (%) | Direction P95 (deg) |",
        "|---|---|---|---|---:|---:|---:|",
    ]
    for row in report["summary"]:
        direction = row.get("direction_p95_deg")
        direction_text = f"{direction:.3g}" if direction is not None and math.isfinite(float(direction)) else "-"
        lines.append(
            f"| {row['mode'].upper()} | {row['regime']} | {row['priority']} | {row['field']} | "
            f"{row['p95_percent']:.4g} | {row['max_percent']:.4g} | {direction_text} |")
    lines.extend([
        "", "## Interpretation limits", "",
        "- Vela current vectors are node-reconstructed diagnostics from the solved state; terminal-current/KCL acceptance remains the conservative current criterion.",
        "- Local magnitude percentages are evaluated only where the Sentaurus magnitude exceeds `1e-3` of its case maximum; values below that active-field floor are retained in CSV but excluded from the primary percentage statistic.",
        "- SRH and SpaceCharge cross zero, so their reported percentages are full-field scale-normalized errors rather than pointwise relative errors.",
        "", "## Artifacts", "",
        f"- Total-current vector figure: `{report['artifacts']['current_vector_png']}`",
        f"- Transport-chain heatmap: `{report['artifacts']['transport_heatmap_png']}`",
        f"- SRH/charge/band heatmap: `{report['artifacts']['source_band_heatmap_png']}`",
        f"- Summary CSV: `{report['artifacts']['summary_csv']}`",
        f"- Long-form comparison CSV: `{report['artifacts']['long_csv']}`",
        "",
    ])
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.check:
        report = json.loads(REPORT_JSON.read_text(encoding="utf-8"))
        for name in ("summary_csv", "long_csv", "current_vector_png", "transport_heatmap_png", "source_band_heatmap_png"):
            artifact = Path(report["artifacts"][name])
            if not artifact.is_file():
                raise RuntimeError(f"Missing artifact: {artifact}")
        print("TransportModels three-regime spatial-field comparison: PASS")
        return 0

    sent = json.loads(SENTAURUS_MANIFEST.read_text(encoding="utf-8"))["states"]
    vela = json.loads(VELA_MANIFEST.read_text(encoding="utf-8"))["states"]
    sent_by_key = {(row["mode"], row["regime"]): row for row in sent}
    vela_by_key = {(row["mode"], row["regime"]): row for row in vela}
    if sent_by_key.keys() != vela_by_key.keys():
        raise RuntimeError("Sentaurus/Vela case coverage differs")
    all_rows: list[dict[str, Any]] = []
    summary: list[dict[str, Any]] = []
    for key in sorted(sent_by_key):
        rows, summaries = compare_case(sent_by_key[key], vela_by_key[key])
        all_rows.extend(rows); summary.extend(summaries)
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    long_csv, summary_csv = OUTPUT_ROOT / "spatial_field_comparison.csv", OUTPUT_ROOT / "field_error_summary.csv"
    write_csv(long_csv, all_rows); write_csv(summary_csv, summary)
    current_png, current_svg = make_current_vector_figure(all_rows)
    transport_png, transport_svg = make_heatmap(
        summary, "transport", "transport_chain_error_heatmap",
        "ElectricField to GradQuasiFermi / Eparallel / Enormal / Mobility: P95 error")
    source_png, source_svg = make_heatmap(
        summary, "source_band", "source_charge_band_error_heatmap",
        "SRH, SpaceCharge and band quantities: P95 error")
    report = {
        "schema": "vela.transportmodels.three_regime_spatial_compare.v1",
        "as_of": "2026-08-24", "status": "complete",
        "node_mapping": "exact; strict R.Substrate interior nodes only",
        "sentaurus_manifest": str(SENTAURUS_MANIFEST.resolve()),
        "vela_manifest": str(VELA_MANIFEST.resolve()),
        "summary": summary,
        "artifacts": {
            "long_csv": str(long_csv.resolve()), "long_csv_sha256": sha256(long_csv),
            "summary_csv": str(summary_csv.resolve()), "summary_csv_sha256": sha256(summary_csv),
            "current_vector_png": str(current_png.resolve()), "current_vector_svg": str(current_svg.resolve()),
            "transport_heatmap_png": str(transport_png.resolve()), "transport_heatmap_svg": str(transport_svg.resolve()),
            "source_band_heatmap_png": str(source_png.resolve()), "source_band_heatmap_svg": str(source_svg.resolve()),
        },
    }
    # Keep generated validation artifacts stable across Windows and POSIX so
    # Git does not treat CR characters as trailing whitespace.
    with REPORT_JSON.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write(json.dumps(report, indent=2, allow_nan=False) + "\n")
    with REPORT_MD.open("w", encoding="utf-8", newline="\n") as stream:
        stream.write(markdown(report))
    print(json.dumps({"status": "complete", "cases": len(sent_by_key), "fields": len(summary), "output": str(OUTPUT_ROOT)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
